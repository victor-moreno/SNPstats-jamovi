'use strict';

// snpImport.js — the browser half of the import.
//
// Named after the analysis on purpose: the jamovi UI compiler emits
// `this.handlers = require('./snpImport')` only for that name (same convention
// as SNPstats' snpPGS.js).
//
// ── What this file is responsible for ───────────────────────────────────────
//
// The selection. jamovi never opens the .bed; this does, and sends only the
// variants the user asked for:
//
//   1. read the .fam in full          -> sample count N
//   2. read the .bim in full          -> variant ID -> row index
//   3. match the user's SNP list      -> indices i1..iK          <- the selection
//   4. slice bytes [3 + i*ceil(N/4), +ceil(N/4)) out of the .bed
//   5. base64 the concatenation into a hidden option
//
// For 1000 SNPs of a 10000-sample panel that is 2.5 MB out of a 1.2 GB file.
// Steps 1-2 are the only large reads and the .bim is ~1% of the .bed.
//
// This file does NO semantic parsing. It finds bytes; R decides what they mean
// and re-validates everything it receives. The one exception is
// splitting .bim lines on whitespace to read column 2, which is byte-finding,
// not interpretation.
//
// ── The limit, and why it is enforced here ──────────────────────────────────
//
// An oversized results message is dropped by nanomsg with no error: the tables
// still render and the genotype columns silently never arrive, so the analysis
// looks like it worked. R cannot detect this. Refusing here,
// before anything is sent, is the only real prevention.

// The chosen genotype files, between picking them and pressing Load.
//
// Keyed by the panel's own DOM node rather than stored on `ui`: assigning an
// unknown property to a jamovi view or control object throws under 'use strict'
// if the object is sealed, and a throw in the pick handler is invisible — the
// files silently never arrive and Load reports that none were chosen. A
// WeakMap keeps no reference alive once the panel is gone.
var PICKED = new WeakMap();

function _panelKey(ui) {
    var ctrl = ui.genoFilename;
    return (ctrl && ctrl.$input && ctrl.$input.length) ? ctrl.$input[0] : null;
}

function _setPicked(ui, g) {
    var k = _panelKey(ui);
    if (k) PICKED.set(k, g);
}

function _getPicked(ui) {
    var k = _panelKey(ui);
    return k ? PICKED.get(k) : null;
}

var CELL_LIMIT   = 3.0e6;             // measured transport cap, ~85% of 3.54M
var PAYLOAD_SAFE = 3.6 * 1024 * 1024; // base64 bytes; websocket ceiling is 4 MiB
// Decompressed bytes. Must not exceed MAX_PAYLOAD_BYTES in R/payload.R, which
// refuses to inflate past it: gzip makes it possible to send a payload that
// clears PAYLOAD_SAFE and is then thrown out by R, and a refusal here says so
// while the file name is still on screen.
var PAYLOAD_RAW_MAX = 8 * 1024 * 1024;
var BED_MAGIC    = [0x6c, 0x1b, 0x01];
// The longest line the scanner will hold before deciding the file is not made
// of lines. A VCF row is one field per sample, so ~40 kB at 10 000 samples and
// ~400 kB at 100 000; this leaves an order of magnitude over the widest cohort
// anyone would import. Without it a file containing no newline at all — a
// binary picked by mistake, or something built to be one long line — is
// accumulated whole in the tab, which is the one thing the chunked scan exists
// to avoid.
var MAX_LINE_CHARS = 8 * 1024 * 1024;

module.exports = {
    view_updated: function(ui) { _inject(ui); },
    view_loaded:  function(ui) { _inject(ui); }
};

// ── option plumbing ─────────────────────────────────────────────────────────

function _setOpt(ui, name, value) {
    if (ui[name] && typeof ui[name].setValue === 'function') {
        ui[name].setValue(value); return;
    }
    if (typeof ui.setOptionValue === 'function') {
        ui.setOptionValue(name, value); return;
    }
    if (typeof ui.getOption === 'function') {
        var o = ui.getOption(name);
        if (o && typeof o.setValue === 'function') o.setValue(value);
    }
}

function _getOpt(ui, name) {
    if (typeof ui.getOptionValue === 'function') return ui.getOptionValue(name);
    if (ui[name] && typeof ui[name].value === 'function') return ui[name].value();
    return undefined;
}

// Clearing the payload also clears the analysis: R sees empty carriers and
// shows the instructions again. The message goes to the status line so the
// file field keeps showing which files are selected.
//
// `kind` says whether the message is a refusal ('error', the default — nothing
// here clears the payload for fun), something the user should notice but did
// not do wrong ('warn'), or an ordinary state change ('ok'). R restates
// anything that is not 'ok' in the results panel, because the status line is a
// one-line text box: "50 000 SNPs × 2 000 samples is too much to transfer"
// scrolls out of sight at the first space, and the explanation of what to do
// instead is entirely invisible.
function _clearGeno(ui, message, kind) {
    _setOpt(ui, 'genoContent', '');
    _setOpt(ui, 'variantContent', '');
    _setOpt(ui, 'sampleContent', '');
    _setOpt(ui, 'sourceDims', '');
    _setOpt(ui, 'loadStatus', message);
    _setOpt(ui, 'loadProblem', kind === 'ok' ? '' : (kind || 'error'));
}

// ── UI ──────────────────────────────────────────────────────────────────────

// Each control is injected independently. An exception in one used to abort
// the whole function, and because Load is built last, a throw while building a
// browse button left the panel with no Load button at all — a far worse failure
// than the missing widget that caused it.
function _guard(what, fn) {
    try { fn(); }
    catch (e) { if (window.console) console.error('snpImport: ' + what, e); }
}

// Same, but for the click handlers, where a silent throw is indistinguishable
// from "the button does nothing". The message goes to the status line so it is
// visible without opening the browser console.
function _guardUI(ui, what, fn) {
    try { fn(); }
    catch (e) {
        if (window.console) console.error('snpImport: ' + what, e);
        try {
            _setOpt(ui, 'loadStatus', what + ' failed: ' + (e && e.message));
            _setOpt(ui, 'loadProblem', 'error');
        } catch (ignored) { /* nothing left to report with */ }
    }
}

function _inject(ui) {
    // Every readable extension. Leaving the older, shorter list here meant the
    // file dialog silently hid .vcf, .ped and .tped, so the formats looked
    // unimplemented even though the readers were in place.
    _guard('genotype browse', function() {
        _browseButton(ui, 'genoFilename', 'snpi-geno',
                      '.bed,.bim,.fam,.ped,.map,.tped,.tfam,.vcf,.gz', true,
                      function(files) { _pickGenotypes(ui, files); });
    });
    // No .gz here: _loadSnpList decodes the *file* as text before it gzips it
    // for transport, so a .gz picked here would be read as mojibake and match
    // nothing. An ID list is small; decompress it first.
    _guard('SNP list browse', function() {
        _browseButton(ui, 'snpListFilename', 'snpi-snps', '.txt,.csv,.tsv', false,
                      function(files) { _loadSnpList(ui, files[0]); },
                      ['snpListContent', 'snpListFilename']);
    });
    // Covariate file browsing is now the native FileSelector control on
    // covFile (jamovi 28.3's File option type) -- no custom JS needed for
    // reading or storing the file. _styleCovFileSelector below only makes
    // its picked-file box sit beside Browse, cosmetic parity with the two
    // custom controls above; see its own comment for why and its limits.
    _guard('covariate file styling', function() { _styleCovFileSelector(ui); });
    _guard('load button', function() { _loadButton(ui); });
}

var COV_STYLE_ID = 'snpimport-inline-fileselector-css';
var COV_ROW_CLASS = 'snpimport-inline-fs';

// Repositions covFile's picked-file box beside its Browse button instead of
// FileSelector's own default layout (button, then a plain list below it),
// to match genoFilename/snpListFilename's look above. A class-based CSS
// rule, not per-element inline styles, because FileSelector redraws its own
// list on every value change (client/analysisui/fileselector.ts:update())
// and a class rule keeps applying to whatever it redraws; only 'body' and
// 'list' being public fields on the live control instance (exposed here as
// ui.covFile.body) is relied on, which is not a documented/stable API — a
// future jamovi client update could silently drop this styling, but nothing
// here touches data or options, so it cannot break the analysis itself.
function _styleCovFileSelector(ui) {
    var ctrl = ui.covFile;
    if (!ctrl || !ctrl.body) return;

    if (!document.getElementById(COV_STYLE_ID)) {
        var style = document.createElement('style');
        style.id = COV_STYLE_ID;
        style.textContent =
            '.' + COV_ROW_CLASS + ' { display: flex; flex-direction: row; ' +
            'align-items: center; flex-wrap: wrap; gap: 6px; }' +
            '.' + COV_ROW_CLASS + ' .jmv-file-selector-list { flex: 1 1 auto; min-width: 0; }' +
            '.' + COV_ROW_CLASS + ' .jmv-file-selector-item { width: 100%; box-sizing: border-box; ' +
            'border: 1px solid #bbb; border-radius: 3px; padding: 3px 6px; ' +
            'background: #f0f0f0; min-height: 14px; }';
        document.head.appendChild(style);
    }

    ctrl.body.classList.add(COV_ROW_CLASS);
}

// The Load button. Choosing files only remembers them; this is what reads
// them, so the panel can ask for the files first and the SNP list second.
//
// The handles live on the `ui` object rather than in a module variable,
// because module scope is shared by every Import analysis on the page and two
// of them would otherwise overwrite each other's files.
function _loadButton(ui) {
    var ctrl = ui.loadStatus;
    if (!ctrl) return;
    var $input = ctrl.$input;
    if (!$input || $input.length === 0) return;

    $input.prop('readonly', true);
    $input.prop('disabled', true);

    if ($input.prev('.snpi-load').length !== 0) return;

    var jq = $input.constructor;
    var $btn = jq('<button type="button" class="snpi-load">Load genotypes</button>').css({
        flexShrink: '0', cursor: 'pointer', padding: '3px 12px', fontSize: '13px',
        border: '1px solid #999', borderRadius: '3px', background: '#e8e8e8',
        whiteSpace: 'nowrap', fontWeight: '600'
    });

    $input.wrap(jq('<div></div>').css({
        display: 'flex', alignItems: 'center', width: '100%', gap: '6px'
    }));
    $input.css({ flex: '1 1 auto', minWidth: 0 });
    $input.before($btn);

    // Nothing else may go between here and .on('click'): a throw after the
    // button is in the DOM but before it is wired leaves a button that looks
    // right and does nothing, which is exactly how a stray copy of
    // _browseButton's clear-button block (referring to `clears`, undefined
    // here) killed Load. tests/js/test-panel.js clicks it to keep it wired.
    $btn.on('click', function(e) {
        e.preventDefault();
        e.stopPropagation();
        _guardUI(ui, 'load', function() { _doLoad(ui); });
    });
}

// Remember the chosen files; read nothing yet.
function _pickGenotypes(ui, files) {
    var g = _groupFiles(files);
    if (g.error) {
        _setPicked(ui, null);
        _setOpt(ui, 'genoFilename', '');
        _clearGeno(ui, g.error);
        return;
    }
    _setPicked(ui, g);
    // Drop whatever the previous Load left behind before sourceFormat moves:
    // the carriers and the format are read together, so bytes from a .bed with
    // sourceFormat 'vcf' would be parsed as VCF text and fail on nonsense.
    _clearGeno(ui, 'ready — press Load genotypes', 'ok');
    _setOpt(ui, 'sourceFormat', g.format);
    _setOpt(ui, 'genoFilename', g.label);
}

// Three shapes are accepted, told apart by extension:
//   .bed + .bim + .fam   binary, read by seeking to per-variant offsets
//   .tped + .tfam        variant-major text
//   .vcf / .vcf.gz       variant-major text, samples named in its own header
//
// All three must be selected in one go where more than one file is involved: a
// web page cannot open a file the user did not pick.
function _groupFiles(files) {

    var vcf = files.filter(function(f) { return /\.vcf(\.gz)?$/i.test(f.name); });
    if (vcf.length === 1 && files.length === 1)
        return { format: 'vcf', files: { vcf: vcf[0] },
                 label: vcf[0].name, gz: /\.gz$/i.test(vcf[0].name) };
    if (vcf.length > 1)
        return { error: 'choose one VCF at a time' };
    if (vcf.length === 1)
        return { error: 'a VCF is self-contained — select it on its own' };

    var stems = {};
    files.forEach(function(file) {
        var m = /^(.*)\.(bed|bim|fam|tped|tfam|ped|map)$/i.exec(file.name);
        if (!m) return;
        var stem = m[1], ext = m[2].toLowerCase();
        if (!stems[stem]) stems[stem] = {};
        stems[stem][ext] = file;
    });

    var names = Object.keys(stems);
    if (names.length === 0)
        return { error: 'choose .bed/.bim/.fam, .tped/.tfam, or a .vcf' };
    if (names.length > 1)
        return { error: 'those files are from different datasets ('
                        + names.join(', ') + ') — select one set' };

    var f = stems[names[0]], stem = names[0];

    if (f.ped || f.map) {
        var pm = ['ped', 'map'].filter(function(e) { return !f[e]; });
        if (pm.length)
            return { error: 'also select ' + stem + '.' + pm.join(' and ' + stem + '.')
                            + ' — pick both at once' };
        return { format: 'ped', files: f, label: stem + '.ped / .map' };
    }

    if (f.tped || f.tfam) {
        var miss = ['tped', 'tfam'].filter(function(e) { return !f[e]; });
        if (miss.length)
            return { error: 'also select ' + stem + '.' + miss.join(' and ' + stem + '.')
                            + ' — pick both at once' };
        return { format: 'tped', files: f, label: stem + '.tped / .tfam' };
    }

    var missing = ['bed', 'bim', 'fam'].filter(function(e) { return !f[e]; });
    if (missing.length)
        return { error: 'also select ' + missing.map(function(e) {
            return stem + '.' + e;
        }).join(' and ') + ' — pick all three at once' };

    return { format: 'bed', files: f, label: stem + '.bed / .bim / .fam' };
}

// Stream a file line by line, keeping only what `keep` accepts.
//
// Chunked rather than read whole: a VCF can be gigabytes, and only the matched
// lines are ever held. Handles a line straddling a chunk boundary, and
// transparently gunzips a .gz through the browser's own DecompressionStream.
// `cap`, when given, is a function returning the largest number of lines worth
// keeping. It is consulted as the scan runs, so a selection that is already over
// budget stops the read instead of finishing it: with no SNP list every line of
// a VCF matches, and a multi-gigabyte file was accumulated whole in the tab
// before the cell count was so much as looked at. The .bed path checks before it
// reads anything, which is the shape this brings the text formats closer to —
// they cannot know the sample count until the header arrives, so the cap is a
// callback rather than a number.
function _scanLines(file, gz, keep, onHeader, cap) {
    var stream;
    if (gz) {
        if (typeof DecompressionStream === 'undefined')
            return Promise.reject(new Error(
                'this browser cannot read .gz — decompress the file first'));
        stream = _gunzip(file);
    } else {
        stream = file.stream();
    }

    var reader = stream.getReader();
    var dec = new TextDecoder('utf-8');
    var tail = '', out = [], scanned = 0;

    function over() {
        if (!cap) return false;
        var max = cap();
        return max !== null && out.length > max;
    }

    function fail() {
        out = []; tail = '';
        throw new Error('no line break in the first ' +
            Math.round(MAX_LINE_CHARS / 1024 / 1024) +
            ' MB — this does not look like a text file');
    }

    function pump() {
        return reader.read().then(function(r) {
            if (r.done) {
                if (tail.length) { scanned++; if (keep(tail)) out.push(tail); }
                return { lines: out, scanned: scanned, truncated: false };
            }
            var text = tail + dec.decode(r.value, { stream: true });
            var parts = text.split('\n');
            tail = parts.pop();               // may be a partial line
            if (tail.length > MAX_LINE_CHARS) {
                return reader.cancel().then(fail, fail);
            }
            for (var i = 0; i < parts.length; i++) {
                var line = parts[i];
                if (line.charCodeAt(line.length - 1) === 13)
                    line = line.slice(0, -1);
                scanned++;
                if (onHeader && onHeader(line)) { out.push(line); continue; }
                if (keep(line)) out.push(line);
            }
            if (over()) {
                // Stop reading and let go of what was collected: the point is
                // not to hold a refused selection in memory either.
                return reader.cancel().then(function() {
                    return { lines: [], scanned: scanned, truncated: true,
                             kept: out.length };
                }, function() {
                    return { lines: [], scanned: scanned, truncated: true,
                             kept: out.length };
                });
            }
            return pump();
        });
    }
    return pump();
}

// Decompress a .gz, one gzip member at a time.
//
// A .vcf.gz almost never has just one. bgzip — what htslib, bcftools, tabix
// and plink2 all write, and what every indexed VCF is — emits an independent
// gzip member per ~64 kB block (BGZF). DecompressionStream('gzip') decodes the
// first member and then fails the stream with "trailing junk found after the
// end of the compressed stream", so piping the file through one of them read
// the first block of a bgzipped VCF and reported a load failure for the rest.
//
// Each member is therefore sliced out and decompressed on its own. A plain
// single-member .gz has exactly one member spanning the file, which is the
// behaviour this replaced.
//
// A member's output is passed on chunk by chunk rather than collected. That is
// invisible for BGZF, whose members are ~64 kB, and it is the whole file for a
// plain .gz — `_memberEnd` reports one member spanning it, so buffering the
// member meant materialising an entire decompressed VCF in the tab before a
// single line had been looked at. Nothing here needs a member whole: the
// consumer is _scanLines, which splits on newlines and keeps almost nothing.
function _gunzip(file) {
    var off = 0, reader = null;
    function next(ctrl) {
        return reader.read().then(function(r) {
            if (!r.done) { ctrl.enqueue(r.value); return; }
            reader = null;
            return pull(ctrl);              // on to the next member, if any
        });
    }
    function pull(ctrl) {
        if (reader) return next(ctrl);
        if (off >= file.size) { ctrl.close(); return; }
        var start = off;
        return _memberEnd(file, start).then(function(end) {
            off = end;
            reader = file.slice(start, end).stream()
                .pipeThrough(new DecompressionStream('gzip')).getReader();
            return next(ctrl);
        });
    }
    return new ReadableStream({
        pull: pull,
        cancel: function(reason) {
            if (reader) { var r = reader; reader = null; return r.cancel(reason); }
        }
    });
}

// Where the gzip member starting at `off` ends.
//
// BGZF states its own block size in a 'BC' extra subfield (SAM spec §4.1), so
// the boundary is read rather than searched for. Anything else is treated as a
// single member running to the end of the file — which is what a plain gzip is.
function _memberEnd(file, off) {
    return file.slice(off, off + 18).arrayBuffer().then(function(b) {
        var h = new Uint8Array(b);
        if (h.length < 12 || h[0] !== 0x1f || h[1] !== 0x8b)
            throw new Error('not a gzip file');
        if ((h[3] & 0x04) === 0) return file.size;      // no FEXTRA: plain gzip
        var xlen = h[10] | (h[11] << 8);
        if (h.length < 18 || xlen < 6) return file.size;
        // The BC subfield is first in every bgzip block; anything else is not
        // BGZF and is left to the single-member path.
        if (h[12] !== 0x42 || h[13] !== 0x43) return file.size;
        var bsize = h[16] | (h[17] << 8);
        return Math.min(off + bsize + 1, file.size);
    });
}

// The ID field differs by format, and is read without splitting the whole line:
// a VCF line for 10000 samples is ~40 kB and only its third field matters.
function _idAt(line, col) {
    var start = 0;
    for (var c = 1; c < col; c++) {
        var nx = _nextSep(line, start);
        if (nx < 0) return null;
        start = nx + 1;
    }
    var end = _nextSep(line, start);
    return end < 0 ? line.slice(start) : line.slice(start, end);
}

function _nextSep(line, from) {
    for (var i = from; i < line.length; i++) {
        var c = line.charCodeAt(i);
        if (c === 9 || c === 32) return i;
    }
    return -1;
}

function _doLoad(ui) {
    var g = _getPicked(ui);
    if (!g) {
        _setOpt(ui, 'loadStatus', 'choose the genotype files first');
        return;
    }
    _setOpt(ui, 'loadStatus', 'reading …');
    // Every load starts here, so clearing the flag here is what retires the
    // previous run's notice — the success paths do not each have to remember.
    _setOpt(ui, 'loadProblem', '');
    if (g.format === 'bed')      _loadGenotypes(ui, g.files);
    else if (g.format === 'ped') _loadPed(ui, g);
    else if (g.format === 'vcf') _loadText(ui, g);
    else {
        // The .tfam is small and it settles the sample count, which is what
        // turns the cell budget into a line budget for the .tped scan. Read it
        // first so the scan can stop early instead of discovering afterwards
        // that it collected far too much.
        _countLines(g.files.tfam).then(function(n) {
            g.tfamCount = n;
            _loadText(ui, g);
        }).catch(function(err) {
            _clearGeno(ui, '(could not read ' + g.files.tfam.name + ': '
                       + err.message + ')');
        });
    }
}

// The text formats: no byte offsets without an index, so the file is scanned
// and only matching lines are kept. That is the whole difference from the .bed
// path — which seeks — and it is why .bed is the format to recommend: this
// reads every byte of the source, however few lines come back.
function _loadText(ui, g) {
    _wantedIds(ui)
        .then(function(wanted) { return _loadTextSel(ui, g, wanted); })
        .catch(function(err) { _clearGeno(ui, '(failed: ' + err.message + ')'); });
}

function _loadTextSel(ui, g, wanted) {

    var isVcf  = g.format === 'vcf';
    var idCol  = isVcf ? 3 : 2;           // VCF ID column, .tped variant id
    var want   = wanted ? new Set(wanted) : null;

    var header = [];
    var onHeader = isVcf ? function(line) {
        if (line.charAt(0) !== '#') return false;
        if (line.charAt(1) !== '#') header.push(line);   // keep only #CHROM
        return false;                                    // not a data line
    } : null;

    var t0 = performance.now();
    var noId = 0;             // lines whose ID column is '.', i.e. unmatchable

    var dosage = _getOpt(ui, 'dosage') === true;
    var cellLimit = dosage ? Math.floor(CELL_LIMIT / 8) : CELL_LIMIT;

    // How many variant lines are still worth keeping, given the samples this
    // file has. Null until the sample count is known: for a VCF that is the
    // moment the #CHROM header goes by, for a .tped it is known before the scan
    // starts because the .tfam is read first.
    var nSamples = isVcf ? null : g.tfamCount;
    function maxLines() {
        if (nSamples === null && isVcf) nSamples = _vcfSampleCount(header) || null;
        return nSamples ? Math.floor(cellLimit / nSamples) : null;
    }

    _scanLines(g.files[isVcf ? 'vcf' : 'tped'], !!g.gz, function(line) {
        if (!line.length || line.charAt(0) === '#') return false;
        if (!want) return true;
        var id = _idAt(line, idCol);
        if (id === '.') noId++;
        return id !== null && want.has(id);
    }, onHeader, maxLines).then(function(res) {

        // Stopped part-way because the selection was already over budget. Said
        // here rather than after the whole file has been read into the tab.
        if (res.truncated) {
            _clearGeno(ui, 'more than ' + maxLines() + ' SNPs match in '
                + g.files[isVcf ? 'vcf' : 'tped'].name + ' (' + res.scanned
                + ' lines read before stopping) — at ' + nSamples
                + ' samples that is past what jamovi can carry'
                + (dosage ? ', and dosage costs 8x more' : '')
                + '. Give a SNP list, or a shorter one.');
            return;
        }

        var lines = res.lines;
        if (lines.length === 0) {
            // A VCF is routinely written with '.' in its ID column, and then no
            // list of rsIDs can ever match it. Without saying so the message
            // reads as "your SNPs are not in this file", which is wrong.
            _clearGeno(ui, 'none of the requested SNPs are in ' +
                       g.files[isVcf ? 'vcf' : 'tped'].name +
                       ' (' + res.scanned + ' lines scanned)' +
                       (noId > 0 ? ' — ' + noId + ' of them carry no ID ("."), '
                                 + 'so they can only be selected after '
                                 + 'plink2 --set-missing-var-ids @:#' : ''));
            return;
        }

        if (isVcf && header.length === 0) {
            _clearGeno(ui, 'no #CHROM header found in the VCF');
            return;
        }

        var payload = (isVcf ? header.concat(lines) : lines).join('\n');

        // The sample count is already known — from the .tfam, read before the
        // scan, or from the #CHROM header the scan went past — so this is the
        // final check on a selection that was capped as it was collected, not
        // the first look at it.
        return Promise.resolve(isVcf ? _vcfSampleCount(header) : nSamples)
          .then(function(nSamples) {
            var cells = nSamples * lines.length;
            var limit = cellLimit;
            if (cells > limit) {
                _clearGeno(ui, lines.length + ' SNPs × ' + nSamples
                    + ' samples is too much to transfer — use at most '
                    + Math.floor(limit / nSamples) + ' SNPs at this sample size'
                    + (dosage ? ' (or untick dosage for 8× more)' : '') + '.');
                return;
            }

            // Gzip is not optional here: genotype text compresses 13-19x and is
            // far over every ceiling raw.
            return _gzipB64(payload).then(function(b64) {
                if (b64.length > PAYLOAD_SAFE) {
                    _clearGeno(ui, 'the selected genotypes are '
                        + (b64.length / 1048576).toFixed(1)
                        + ' MB even compressed — select fewer SNPs.');
                    return;
                }
                _setOpt(ui, 'genoContent', b64);
                _setOpt(ui, 'variantContent', '');
                if (isVcf) _setOpt(ui, 'sampleContent', '');
                else return _textOption(ui, 'sampleContent', g.files.tfam);
            }).then(function() {
                if (!_getOpt(ui, 'genoContent')) return;
                _setOpt(ui, 'loadStatus', 'loaded ' + lines.length
                    + (wanted ? ' of ' + wanted.length : '') + ' SNPs, '
                    + nSamples + ' samples — scanned ' + res.scanned
                    + ' lines in ' + ((performance.now() - t0) / 1000).toFixed(1) + ' s');
            });
        });
    }).catch(function(err) {
        _clearGeno(ui, '(failed: ' + err.message + ')');
    });
}

// .ped is the sample-major one, so selecting means dropping *columns* from
// every line rather than skipping lines. Nothing can make that cheap — every
// byte is read and every field split, however few variants come back — which
// is why the panel recommends converting instead. It is still the simplest
// format to load for a small file, so it is supported with an honest warning.
function _loadPed(ui, g) {

    var t0 = performance.now();

    Promise.all([_text(g.files.map), _wantedIds(ui)]).then(function(r) {

        var mapTxt = r[0], wanted = r[1];
        var mapLines = _lines(mapTxt).filter(function(l) { return l.trim().length > 0; });

        var pick = [];                       // indices into the .map
        if (wanted) {
            var at = new Map();
            for (var i = 0; i < mapLines.length; i++) {
                var id = mapLines[i].trim().split(/[ \t]+/)[1];
                if (id !== undefined && !at.has(id)) at.set(id, i);
            }
            for (var w = 0; w < wanted.length; w++) {
                var k = at.get(wanted[w]);
                if (k !== undefined) pick.push(k);
            }
        } else {
            for (var m = 0; m < mapLines.length; m++) pick.push(m);
        }

        if (pick.length === 0) {
            _clearGeno(ui, 'none of the requested SNPs are in ' + g.files.map.name);
            return;
        }

        // 1-based field positions of each selected variant's allele pair
        var cols = [];
        for (var c = 0; c < pick.length; c++) {
            cols.push(6 + 2 * pick[c]);      // 0-based after the 6 sample fields
            cols.push(6 + 2 * pick[c] + 1);
        }

        var big = g.files.ped.size > 200 * 1024 * 1024;
        if (big)
            _setOpt(ui, 'loadStatus', 'reading a '
                + (g.files.ped.size / 1048576).toFixed(0)
                + ' MB .ped — every line must be split; this may take a while');

        var nFields = 6 + 2 * mapLines.length;
        var out = [], mismatched = 0;
        return _scanLines(g.files.ped, false, function(line) {
            if (!line.trim().length) return false;
            var f = line.trim().split(/[ \t]+/);
            // A line with fewer fields than the .map accounts for would take
            // `undefined` for the missing calls and reach R as a third allele,
            // which is reported as a parse error about the wrong thing.
            if (f.length !== nFields) { mismatched++; return false; }
            var row = f.slice(0, 6);
            for (var q = 0; q < cols.length; q++) row.push(f[cols[q]]);
            out.push(row.join(' '));
            return false;                    // collected here, not by _scanLines
        }, null).then(function(res) {

            if (mismatched > 0) {
                _clearGeno(ui, mismatched + ' line'
                    + (mismatched === 1 ? '' : 's')
                    + ' of ' + g.files.ped.name + ' do not have '
                    + nFields + ' fields (6 + 2 × ' + mapLines.length
                    + ' variants) — the .ped and the .map do not match');
                return;
            }

            if (out.length === 0) { _clearGeno(ui, 'the .ped has no data lines'); return; }

            var cells = out.length * pick.length;
            var dosage = _getOpt(ui, 'dosage') === true;
            var limit = dosage ? Math.floor(CELL_LIMIT / 8) : CELL_LIMIT;
            if (cells > limit) {
                _clearGeno(ui, pick.length + ' SNPs × ' + out.length
                    + ' samples is too much to transfer — use at most '
                    + Math.floor(limit / out.length) + ' SNPs at this sample size.');
                return;
            }

            return _gzipB64(out.join('\n')).then(function(b64) {
                if (b64.length > PAYLOAD_SAFE) {
                    _clearGeno(ui, 'the selected genotypes are too large even '
                        + 'compressed — select fewer SNPs.');
                    return;
                }
                _setOpt(ui, 'genoContent', b64);
                _setOpt(ui, 'sampleContent', '');       // a .ped carries its own
                return _gzipB64(pick.map(function(i) { return mapLines[i]; }).join('\n'))
                    .then(function(mb64) { _setOpt(ui, 'variantContent', mb64); });
            }).then(function() {
                if (!_getOpt(ui, 'genoContent')) return;
                _setOpt(ui, 'loadStatus', 'loaded ' + pick.length
                    + (wanted ? ' of ' + wanted.length : '') + ' SNPs, '
                    + out.length + ' samples from .ped in '
                    + ((performance.now() - t0) / 1000).toFixed(1) + ' s'
                    + ' — .bed would be far faster for this');
            });
        });
    }).catch(function(err) {
        _clearGeno(ui, '(failed: ' + err.message + ')');
    });
}

function _vcfSampleCount(header) {
    if (header.length === 0) return 0;
    var f = header[header.length - 1].split(/[\t ]+/);
    return Math.max(0, f.length - 9);
}

function _countLines(file) {
    return _scanLines(file, false, function() { return true; }, null)
        .then(function(r) { return r.lines.filter(function(l) {
            return l.trim().length > 0; }).length; });
}

function _textOption(ui, name, file) {
    return _text(file).then(function(txt) {
        return _gzipB64(txt).then(function(b64) { _setOpt(ui, name, b64); });
    });
}

// gzip through the browser's own CompressionStream, then base64. Falls back to
// plain base64 where the API is missing, which only affects how much fits.
function _gzipB64(text) {
    var bytes = new TextEncoder().encode(text);
    if (typeof CompressionStream === 'undefined')
        return Promise.resolve(_b64(bytes));
    var stream = new Blob([bytes]).stream()
                    .pipeThrough(new CompressionStream('gzip'));
    return new Response(stream).arrayBuffer().then(function(buf) {
        return _b64(new Uint8Array(buf));
    });
}

// The requested IDs, from whichever source is filled; null means "everything".
//
// A promise because snpListContent is gzipped and DecompressionStream is async.
// Every caller awaits it before touching a genotype file, which is the same
// ordering as before -- nothing is read until the selection is known.
function _wantedIds(ui) {
    var listB64 = _getOpt(ui, 'snpListContent');
    var listTxt = _getOpt(ui, 'snpListText');
    if (listB64 && listB64.length)
        return _decodeB64Text(listB64).then(function(txt) { return _idsFrom(txt); });
    if (listTxt && String(listTxt).trim().length)
        return Promise.resolve(String(listTxt).split(/[\s,;]+/)
                   .filter(function(s) { return s.length > 0; }));
    return Promise.resolve(null);
}

// Shows $input and its $clr together when ctrlName has a value, hides both
// (Browse stays visible) when it does not -- matching covFile's native
// FileSelector, whose file box+remove-button are absent until a file is
// picked (fileselector.css: '.jmv-file-selector-list:empty { display: none }').
//
// Called both from _browseButton on every view_updated/view_loaded, and
// directly from wherever ctrlName's value is set (the $clr click handler,
// _loadSnpList's callbacks) -- not only the former. An earlier version of
// this control hid the box until a file was chosen and relied on
// view_updated alone to un-hide it, which never ran in time: the client
// fires view_updated from its `view.ready` event, only when the options are
// re-initialised from the server under a new id, not when this file's own
// code sets an option. The box stayed hidden display:none for the whole
// session and there was no way to see or remove a chosen file at all.
// Calling this at the point of change, synchronously, is what fixes that.
function _syncClearable(ui, ctrlName, $input, $clr) {
    var val = _getOpt(ui, ctrlName);
    var has = val !== undefined && val !== null && String(val).length > 0;
    $input.css('display', has ? '' : 'none');
    if ($clr) $clr.css('display', has ? '' : 'none');
}

// `clears` names the options a ✕ button should empty. Without it there is no
// way to unpick a file: the field is deliberately read-only, so a wrong choice
// can only be replaced, never removed.
function _browseButton(ui, ctrlName, cls, accept, multiple, onPick, clears) {

    var ctrl = ui[ctrlName];
    if (!ctrl) return;
    var $input = ctrl.$input;
    if (!$input || $input.length === 0) return;

    // The field is a passive readout of what the button loaded. A typed path
    // would be saved into the .omv and resolved on whoever opened it, so
    // picking a file is the only way in.
    $input.prop('readonly', true);
    $input.prop('disabled', true);
    $input.css('cursor', 'default');

    if ($input.prev('.' + cls).length !== 0) {
        // Already injected on an earlier cycle -- just keep visibility current.
        if (clears && clears.length)
            _syncClearable(ui, ctrlName, $input, $input.next('.' + cls + '-clr'));
        return;
    }

    // Plain text, matching jamovi's own native FileSelector button (which
    // reads 'Browse…' too, hardcoded in client/analysisui/fileselector.ts)
    // rather than a hand-drawn folder icon — this control can't actually be
    // a FileSelector (see the header comment on why), but it should still
    // look like jamovi's own file pickers, not a one-off.
    var jq = $input.constructor;
    var $btn = jq('<button type="button" class="' + cls + '">Browse…</button>').css({
        flexShrink: '0', cursor: 'pointer', padding: '3px 12px', fontSize: '14px',
        lineHeight: '1.4', border: '1px solid #bbb', borderRadius: '3px',
        background: '#f0f0f0', whiteSpace: 'nowrap'
    });

    $input.wrap(jq('<div></div>').css({
        display: 'flex', alignItems: 'center', width: '100%', gap: '4px'
    }));
    $input.css({ flex: '1 1 auto', minWidth: 0 });
    $input.before($btn);

    if (clears && clears.length) {
        var $clr = jq('<button type="button" class="' + cls + '-clr" title="Remove this file">✕</button>').css({
            flexShrink: '0', cursor: 'pointer', padding: '1px 7px', fontSize: '13px',
            lineHeight: '1.4', border: '1px solid #bbb', borderRadius: '3px',
            background: '#f0f0f0', whiteSpace: 'nowrap'
        });
        $input.after($clr);

        $clr.on('click', function(e) {
            e.preventDefault();
            e.stopPropagation();
            for (var i = 0; i < clears.length; i++) _setOpt(ui, clears[i], '');
            // The SNP list decided which variants are in the payload, so the
            // payload has to go with it — otherwise removing the list leaves
            // the analysis running on its selection with nothing on screen
            // saying so.
            if (ctrlName === 'snpListFilename')
                _clearGeno(ui, 'SNP list removed — press Load genotypes to re-read', 'ok');
            _syncClearable(ui, ctrlName, $input, $clr);
        });

        _syncClearable(ui, ctrlName, $input, $clr);
    }

    $btn.on('click', function(e) {
        e.preventDefault();
        e.stopPropagation();
        var fi = document.createElement('input');
        fi.type = 'file';
        fi.accept = accept;
        if (multiple) fi.multiple = true;
        fi.style.display = 'none';
        document.body.appendChild(fi);
        fi.addEventListener('change', function() {
            if (fi.files && fi.files.length)
                _guardUI(ui, 'choosing files', function() {
                    onPick(Array.prototype.slice.call(fi.files));
                });
            document.body.removeChild(fi);
        });
        fi.click();
    });
}

// ── helpers ─────────────────────────────────────────────────────────────────

function _b64(bytes) {
    // btoa takes a binary string, and String.fromCharCode blows the argument
    // limit past ~100k, so build it in chunks.
    var CHUNK = 0x8000, parts = [];
    for (var i = 0; i < bytes.length; i += CHUNK)
        parts.push(String.fromCharCode.apply(null, bytes.subarray(i, i + CHUNK)));
    return btoa(parts.join(''));
}

function _text(file) {
    return file.arrayBuffer().then(function(b) {
        return new TextDecoder('utf-8').decode(new Uint8Array(b));
    });
}

function _lines(txt) {
    return txt.replace(/\r\n?/g, '\n').split('\n');
}

// ── the SNP selection list ──────────────────────────────────────────────────

// Gzipped and size-checked like every other upload. A published PGS Catalog
// weights file is the case this exists for and the one most likely to be big:
// hundreds of thousands of rows, plain text, several MB. Sent raw and over the
// ceiling it is dropped by nanomsg in silence — no error in the browser, none
// in R, and an analysis that looks like it worked. Gzip buys ~5-10x on an
// ID/weights list; the two refusals below are what happens past that.
function _loadSnpList(ui, file) {
    // The box+\u2715 are hidden until snpListFilename has a value (_syncClearable);
    // _inject() re-runs it right away in every branch below rather than
    // waiting for the next view_updated, which is what makes the box appear
    // the moment a name (or an error) is actually there to show.
    var refuse = function(why) {
        _setOpt(ui, 'snpListContent', '');
        _setOpt(ui, 'snpListFilename', file.name + ' \u2014 ' + why);
        _inject(ui);
    };
    _text(file).then(function(txt) {
        var raw = new TextEncoder().encode(txt);
        if (raw.length > PAYLOAD_RAW_MAX) {
            refuse('too large (' + (raw.length / 1048576).toFixed(1)
                 + ' MB of text) \u2014 cut it to the SNPs you need');
            return;
        }
        return _gzipB64(txt).then(function(b64) {
            if (b64.length > PAYLOAD_SAFE) {
                refuse('too large (' + (b64.length / 1048576).toFixed(1)
                     + ' MB even compressed) \u2014 cut it to the SNPs you need');
                return;
            }
            _setOpt(ui, 'snpListContent', b64);
            _setOpt(ui, 'snpListFilename', file.name);
            _inject(ui);
        });
    }).catch(function(err) {
        _setOpt(ui, 'snpListContent', '');
        _setOpt(ui, 'snpListFilename', '(could not read ' + file.name + ': ' + err.message + ')');
        _inject(ui);
    });
}

// Header names that can hold the variant ID, best first. Must stay identical,
// in this order, to SNP_ID_NAMES in R/payload.R: both halves read the same file
// and have to reach the same column. Preference, not position — on a header
// like `id, rsid, effect_allele` R used to take the leftmost match while this
// took the preferred one, so the browser sliced the rsIDs and the report said
// the contents of `id` had all gone missing. test-groupfiles.js checks the two
// lists against each other.
var SNP_ID_NAMES = ['rsid', 'variant_id', 'snpid', 'snp', 'id', 'name'];

// Same recognition rules as parse_snp_selection() in R, kept deliberately
// crude: this only needs to know *which IDs* to slice. R re-parses the file
// itself for everything else, including the effect alleles.
function _idsFrom(txt) {
    var lines = _lines(txt).map(function(l) { return l.trim(); })
                           .filter(function(l) { return l.length > 0; });
    var body = lines.filter(function(l) { return l.charAt(0) !== '#'; });
    if (body.length === 0) return [];

    var hdr = body[0].toLowerCase().split(/[\t,;]/).map(function(h) { return h.trim(); });
    var idCol = -1;
    SNP_ID_NAMES.forEach(function(k) {
        if (idCol < 0) idCol = hdr.indexOf(k);
    });

    if (idCol < 0)
        return body.join(' ').split(/[\s,;]+/).filter(function(s) { return s.length > 0; });

    var out = [];
    for (var i = 1; i < body.length; i++) {
        var f = body[i].split(/[\t,;]/);
        if (f.length > idCol) {
            var v = f[idCol].trim();
            if (v.length > 0) out.push(v);
        }
    }
    return out;
}

// ── the genotype trio ───────────────────────────────────────────────────────

function _loadGenotypes(ui, f) {

    // The SNP list, from whichever source is filled. Needed *before* the .bed
    // is touched — that is the whole point of the design.
    var wanted = _wantedIds(ui);                   // resolves null = everything

    Promise.all([wanted, _text(f.fam), _text(f.bim), _magic(f.bed)])
        .then(function(r) {
            var ids = r[0], famTxt = r[1], bimTxt = r[2], magicErr = r[3];
            if (magicErr) { _clearGeno(ui, magicErr); return; }
            // No list means every variant in the .bim. Fine for a curated
            // panel; the capacity check below is what stops a whole GWAS.
            return _slice(ui, f, ids, famTxt, bimTxt);
        })
        .catch(function(err) {
            _clearGeno(ui, '(failed: ' + err.message + ')');
        });
}

// The inverse of _gzipB64: base64 out, gunzip if the bytes say gzip. Detection
// is by magic rather than by a flag, so a payload written by the fallback path
// (no CompressionStream, plain base64) reads back the same way.
function _decodeB64Text(b64) {
    var bin = atob(b64);
    var arr = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) arr[i] = bin.charCodeAt(i);
    var dec = new TextDecoder('utf-8');
    if (arr.length < 2 || arr[0] !== 0x1f || arr[1] !== 0x8b ||
        typeof DecompressionStream === 'undefined')
        return Promise.resolve(dec.decode(arr));
    var stream = new Blob([arr]).stream()
                    .pipeThrough(new DecompressionStream('gzip'));
    return new Response(stream).arrayBuffer().then(function(buf) {
        return dec.decode(new Uint8Array(buf));
    });
}

// Check the .bed header before reading anything else, so a wrong file is
// rejected in milliseconds rather than after parsing a 20 MB .bim.
function _magic(bed) {
    return bed.slice(0, 3).arrayBuffer().then(function(b) {
        var m = new Uint8Array(b);
        if (m.length < 3 || m[0] !== BED_MAGIC[0] || m[1] !== BED_MAGIC[1])
            return f_err(bed.name + ' is not a .bed file');
        if (m[2] === 0x00)
            return f_err('this .bed is sample-major; convert with '
                       + 'plink --bfile <stem> --make-bed --out <stem>');
        if (m[2] !== 0x01)
            return f_err('unrecognised .bed layout byte');
        return null;
    });
    function f_err(s) { return s; }
}

function _slice(ui, f, wantedIds, famTxt, bimTxt) {

    var famLines = _lines(famTxt).filter(function(l) { return l.trim().length > 0; });
    var n = famLines.length;
    if (n === 0) throw new Error('.fam is empty');

    var bpv = Math.ceil(n / 4);

    // .bim -> variant ID -> row index. This is the only large parse, and at
    // ~1% of the .bed it is what buys us random access into the genotypes.
    // Blank lines are dropped first: a line's position in this array is the
    // variant's block number in the .bed, and a blank one occupies no block.
    var bimLines = _lines(bimTxt).filter(function(l) { return l.trim().length > 0; });

    // The three files have to be the same three files. Nothing downstream can
    // find out that they are not: the offset of variant i is computed from the
    // .fam's sample count and the .bim's line number, and R re-derives the
    // expected payload length from those same two numbers, so a wrong .fam or a
    // stale .bim agrees with itself perfectly and decodes a different variant's
    // bytes into a confident-looking genotype. A .bed states its own size, and
    // that is the one number in the trio neither of the other two can fake.
    var expected = 3 + bimLines.length * bpv;
    if (f.bed.size !== expected) {
        _clearGeno(ui, f.bed.name + ' is ' + f.bed.size + ' bytes, but '
            + f.bim.name + ' (' + bimLines.length + ' variants) and '
            + f.fam.name + ' (' + n + ' samples) describe ' + expected
            + ' — these files are not from the same dataset. Re-export the set '
            + 'with plink --bfile <stem> --make-bed --out <stem>.');
        return;
    }

    var index = new Map();
    for (var i = 0; i < bimLines.length; i++) {
        var parts = bimLines[i].trim().split(/[ \t]+/);
        if (parts.length >= 6 && !index.has(parts[1])) index.set(parts[1], i);
    }

    // The selection. A null list means take every variant in the .bim.
    var requested = wantedIds ? wantedIds.length : bimLines.length;
    var idx = [], keptBim = [], missing = 0;

    if (wantedIds) {
        for (var w = 0; w < wantedIds.length; w++) {
            var at = index.get(wantedIds[w]);
            if (at === undefined) { missing++; continue; }
            idx.push(at);
            keptBim.push(bimLines[at].trim());
        }
    } else {
        // By line, not by distinct ID. Walking the index instead meant that a
        // .bim converted from a VCF without --set-missing-var-ids — every
        // variant named '.' — imported exactly one variant and reported it as
        // a success. R gives the nameless ones a chr:bp column name.
        for (var m = 0; m < bimLines.length; m++) {
            idx.push(m);
            keptBim.push(bimLines[m].trim());
        }
    }

    if (idx.length === 0) {
        _clearGeno(ui, 'none of the ' + requested + ' requested SNPs are in ' + f.bim.name);
        return;
    }

    // Refuse before sending, not after. Both walls, whichever binds first.
    var cells = n * idx.length;
    var dosage = _getOpt(ui, 'dosage') === true;
    var limit = dosage ? Math.floor(CELL_LIMIT / 8) : CELL_LIMIT;
    if (cells > limit) {
        var maxSnps = Math.floor(limit / n);
        _clearGeno(ui, idx.length + ' SNPs × ' + n + ' samples is too much to '
            + 'transfer — jamovi would drop the columns without an error. '
            + 'Use at most ' + maxSnps + ' SNPs at this sample size'
            + (dosage ? ' (or untick dosage for 8× more)' : '') + '.');
        return;
    }

    var payload = idx.length * bpv * 4 / 3;
    if (payload > PAYLOAD_SAFE) {
        _clearGeno(ui, 'the selected genotypes are '
            + (payload / 1048576).toFixed(1) + ' MB encoded, over what jamovi '
            + 'can carry. Select fewer SNPs.');
        return;
    }

    // Read only the selected blocks. File.slice is random access: the browser
    // pages in these ranges and never touches the rest of the file.
    var reads = idx.map(function(i) {
        var off = 3 + i * bpv;
        return f.bed.slice(off, off + bpv).arrayBuffer();
    });

    return Promise.all(reads).then(function(buffers) {
        var all = new Uint8Array(idx.length * bpv);
        for (var j = 0; j < buffers.length; j++) {
            var b = new Uint8Array(buffers[j]);
            if (b.length !== bpv)
                throw new Error('.bed is truncated at variant ' + (idx[j] + 1)
                              + ' (got ' + b.length + ' of ' + bpv + ' bytes)');
            all.set(b, j * bpv);
        }

        var enc = new TextEncoder();
        // What the three files said about each other, so R can re-check the
        // arithmetic instead of only re-checking the slice it was handed. The
        // slice always agrees with itself; this does not.
        _setOpt(ui, 'sourceDims', f.bed.size + ',' + bimLines.length + ',' + n);
        _setOpt(ui, 'genoContent', _b64(all));
        _setOpt(ui, 'variantContent', _b64(enc.encode(keptBim.join('\n'))));
        _setOpt(ui, 'sampleContent', _b64(enc.encode(famLines.join('\n'))));
        _setOpt(ui, 'loadStatus',
            'loaded ' + idx.length
            + (wantedIds ? ' of ' + requested : '') + ' SNPs, ' + n + ' samples'
            + (missing ? ', ' + missing + ' not found' : ''));
    });
}
