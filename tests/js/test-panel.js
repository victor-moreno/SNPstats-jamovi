'use strict';

// The panel, injected and clicked.
//
//   node tests/js/test-panel.js
//
// test-groupfiles.js reads the source; this one runs it against the real
// fixtures, so a handler that throws before it is attached shows up as a dead
// button rather than as a passing regex.

const path = require('path');
const assert = require('assert');
const H = require('./harness');

const FX = path.join(__dirname, '..', '..', 'data-raw', 'fixtures');
const small = (ext) => path.join(FX, 'small', 'small' + ext);

let passed = 0, failed = 0, skipped = 0;
const only = process.argv[2];

// Only the `small` tier is committed; `medium` and `large` run to gigabytes and
// are gitignored, so a fresh clone has neither. Tests that need them say so and
// are skipped, the way the R suite's skip_without_fixture() does -- otherwise
// the whole JS suite fails on any machine that has not run make_fixtures.sh,
// which includes anyone who has just cloned this to look at it.
class Skip extends Error {}
function needFixture(...files) {
    const missing = files.filter(f => !require('fs').existsSync(f));
    if (missing.length)
        throw new Skip(path.relative(FX, missing[0]) +
                       ' not generated — bash data-raw/make_fixtures.sh medium large');
}

const tests = [];
function test(what, fn) { tests.push([what, fn]); }

async function main() {
    for (const [what, fn] of tests) {
        if (only && !what.includes(only)) continue;
        try { await fn(); passed++; console.log('ok   ' + what); }
        catch (e) {
            if (e instanceof Skip) {
                skipped++; console.log('skip ' + what + ' (' + e.message + ')');
                continue;
            }
            failed++;
            console.error('FAIL ' + what + '\n     ' + (e && e.stack || e));
        }
    }
    console.log(`\n${passed} pass, ${failed} fail` + (skipped ? `, ${skipped} skip` : ""));
    process.exit(failed > 0 ? 1 : 0);
}

// A fresh module per test: PICKED is module state keyed by the panel's node,
// but the option values are not, and a stale require would share them.
function load(queue) {
    H.installDocument(queue);
    delete require.cache[require.resolve('../../jamovi/js/snpimport.js')];
    return require('../../jamovi/js/snpimport.js');
}

// Inject the panel, pick `files` through the genotype browse button, press
// Load, and hand back the option values.
async function loadFiles(files, opts) {
    const queue = [files];
    const mod = load(queue);
    const ui = H.makeUi(opts);
    mod.view_loaded(ui);
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);
    return ui;
}

// ── the button exists and is wired ──────────────────────────────────────────

test('every browse button and Load are injected', () => {
    const mod = load([]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    // Covariate browsing is now the native FileSelector control (jamovi
    // 28.3's File option type), not a custom-JS button, so it is not in
    // this list any more.
    for (const cls of ['snpi-geno', 'snpi-snps', 'snpi-load'])
        assert.strictEqual(H.find(ui.root, cls).length, 1, 'missing ' + cls);
});

test('the Load button has a click handler', () => {
    // The regression: a stray reference in _loadButton threw after the button
    // was inserted but before .on('click') ran, so the panel looked complete
    // and pressing Load did nothing at all.
    const mod = load([]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    const btn = H.find(ui.root, 'snpi-load')[0];
    assert.ok(btn.handlers.click && btn.handlers.click.length,
              'Load was injected without a click handler');
});

test('pressing Load with no files chosen says so', () => {
    const mod = load([]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-load');
    assert.match(ui.values.loadStatus, /choose the genotype files/);
});

test('view_updated does not inject a second set of buttons', () => {
    const mod = load([]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    mod.view_updated(ui);
    mod.view_updated(ui);
    assert.strictEqual(H.find(ui.root, 'snpi-load').length, 1);
    assert.strictEqual(H.find(ui.root, 'snpi-geno').length, 1);
});

test('choosing files fills the name and arms Load without reading', () => {
    const mod = load([[H.file(small('.bed')), H.file(small('.bim')),
                       H.file(small('.fam'))]]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-geno');
    assert.match(ui.values.genoFilename, /small\.bed/);
    assert.strictEqual(ui.values.sourceFormat, 'bed');
    assert.match(ui.values.loadStatus, /press Load/);
    assert.strictEqual(ui.values.genoContent, '',
                       'picking files must not read them');
});

// ── every format loads ──────────────────────────────────────────────────────

const FORMATS = [
    ['bed',  () => [H.file(small('.bed')), H.file(small('.bim')),
                    H.file(small('.fam'))]],
    ['tped', () => [H.file(small('.tped')), H.file(small('.tfam'))]],
    ['ped',  () => [H.file(small('.ped')), H.file(small('.map'))]],
    ['vcf',  () => [H.file(small('.vcf'))]],
    ['vcf',  () => [H.file(small('.vcf.gz'))]]
];

for (const [fmt, files] of FORMATS) {
    test('a ' + fmt + ' set loads every variant when no list is given', async () => {
        const ui = await loadFiles(files());
        assert.strictEqual(ui.values.sourceFormat, fmt);
        assert.ok(ui.values.genoContent.length > 0,
                  'nothing was loaded: ' + ui.values.loadStatus);
        assert.match(ui.values.loadStatus, /^loaded /);
    });
}

test('a SNP list restricts what is sent, in every format', async () => {
    const ids = ['snp2', 'snp5', 'snp9'];
    for (const [fmt, files] of FORMATS) {
        const ui = await loadFiles(files(), { snpListText: ids.join('\n') });
        assert.ok(ui.values.genoContent.length > 0,
                  fmt + ': nothing loaded — ' + ui.values.loadStatus);
        assert.match(ui.values.loadStatus, new RegExp('loaded 3 of 3 SNPs'),
                     fmt + ': ' + ui.values.loadStatus);
    }
});

test('a bgzipped VCF is read whole, not just its first block', async () => {
    // bgzip writes one gzip member per ~64 kB block, and a single
    // DecompressionStream stops at the first. Comparing the payloads is what
    // makes a partial read visible: a truncated one still looks like a
    // successful load.
    // medium is over the cell budget whole, so it is selected down to a few
    // SNPs — which also puts the matches well past the first bgzip block.
    needFixture(path.join(FX, 'medium', 'medium.vcf.gz'));
    for (const [tier, opts] of [['small', {}],
                                ['medium', { snpListText: 'snp900 snp4500' }]]) {
        const vcf = path.join(FX, tier, tier + '.vcf');
        const plain = await loadFiles([H.file(vcf)], opts);
        const gzip  = await loadFiles([H.file(vcf + '.gz')], opts);
        assert.ok(plain.values.genoContent.length > 0,
                  tier + ': nothing loaded — ' + plain.values.loadStatus);
        assert.strictEqual(gzip.values.genoContent, plain.values.genoContent,
            tier + ': the gzipped VCF did not decode to the same bytes');
        assert.strictEqual(gzip.values.loadStatus.replace(/ in .*/, ''),
                           plain.values.loadStatus.replace(/ in .*/, ''));
    }
});

// Both committed .gz fixtures are BGZF, because that is what bgzip and plink2
// write. A file someone gzipped themselves is one member spanning the whole
// file, and that path -- the one where a member is the entire dataset -- had no
// fixture of its own. These build theirs.
const TMP = path.join(__dirname, '..', '..', '.tmp', 'js');

function plainGz(name, text) {
    require('fs').mkdirSync(TMP, { recursive: true });
    const at = path.join(TMP, name);
    require('fs').writeFileSync(at, require('zlib').gzipSync(Buffer.from(text)));
    return at;
}

// The same, for a file too large to hold: `line(i)` is written `n` times.
//
// Backpressure is waited on rather than ignored, which for a fixture writer
// looks like fussing over nothing and is not: writing 64 MB into a gzip stream
// without it queues all 64 MB inside the stream, and the test that follows
// cannot tell that allocation apart from the one it is looking for.
async function bigGz(name, header, line, n) {
    const fs = require('fs');
    fs.mkdirSync(TMP, { recursive: true });
    const at = path.join(TMP, name);
    const gz = require('zlib').createGzip();
    const out = fs.createWriteStream(at);
    gz.pipe(out);

    const put = (s) => new Promise((ok, no) => {
        gz.write(s, (err) => err ? no(err) : ok());
    });
    await put(header.join('\n') + '\n');
    let batch = [];
    for (let i = 0; i < n; i++) {
        batch.push(line(i));
        if (batch.length === 1000) { await put(batch.join('\n') + '\n'); batch = []; }
    }
    if (batch.length) await put(batch.join('\n') + '\n');

    gz.end();
    await new Promise((ok, no) => { out.on('finish', ok); out.on('error', no); });
    return at;
}

test('a plain (non-bgzip) .vcf.gz decodes to the same bytes as the .vcf', async () => {
    const vcf = small('.vcf');
    const gz = plainGz('plain.vcf.gz',
                       require('fs').readFileSync(vcf, 'utf8'));
    const plain = await loadFiles([H.file(vcf)]);
    const zipped = await loadFiles([H.file(gz)]);
    assert.ok(plain.values.genoContent.length > 0,
              'nothing loaded — ' + plain.values.loadStatus);
    assert.strictEqual(zipped.values.genoContent, plain.values.genoContent);
});

test('a plain .vcf.gz is streamed, not held whole while it is scanned', async () => {
    // A peak is a property of the process, and by the time the suite reaches
    // here it is carrying tens of MB of fixtures it has already read and not
    // yet collected — enough to swamp the reading. So this one test measures in
    // a child of its own, which starts out holding nothing.
    if (!process.env.SNPI_MEM_CHILD) {
        try {
            require('child_process').execFileSync(
                process.execPath, [__filename, 'streamed'],
                { env: Object.assign({}, process.env, { SNPI_MEM_CHILD: '1' }) });
        } catch (e) {
            throw new Error(String(e.stdout || '') + String(e.stderr || ''));
        }
        return;
    }

    // The single-member path used to collect the member through
    // Response.arrayBuffer() before a line of it was looked at, so opening a
    // gzipped VCF cost the tab the whole *decompressed* file: 64 MB of VCF in
    // 0.5 MB of .gz measured at 198 MB peak, and a real one is far larger.
    // Bounded memory is the claim, so it is measured rather than assumed — and
    // the fixture is written through a stream rather than built in a string,
    // because a 64 MB Buffer left uncollected in this process reads exactly
    // like the leak being tested for.
    const src = require('fs').readFileSync(small('.vcf'), 'utf8').split('\n');
    const one = src.filter(l => l.length && l.charAt(0) !== '#')[0].split('\t');
    const rows = Math.ceil(64 * 1024 * 1024 / one.join('\t').length);
    const gz = await bigGz('big.vcf.gz', src.filter(l => l.charAt(0) === '#'),
                           function (i) { one[2] = 'bigsnp' + i; return one.join('\t'); },
                           rows);

    let peak = 0;
    const tick = setInterval(() => {
        const m = process.memoryUsage();
        peak = Math.max(peak, m.arrayBuffers + m.external);
    }, 2);
    let ui;
    try { ui = await loadFiles([H.file(gz)], { snpListText: 'bigsnp7' }); }
    finally { clearInterval(tick); }

    assert.match(ui.values.loadStatus, /loaded 1 of 1 SNPs/);
    assert.ok(peak < 32 * 1024 * 1024,
              'held ' + Math.round(peak / 1048576) + ' MB of a 64 MB file at once');
});

test('a file with no line breaks is refused instead of accumulated', async () => {
    // The scanner holds the unfinished last line between chunks, so a file that
    // never ends one -- a binary picked by mistake -- is read into that buffer
    // in its entirety, which is what the chunked scan exists to prevent.
    require('fs').mkdirSync(TMP, { recursive: true });
    const at = path.join(TMP, 'oneline.vcf');
    require('fs').writeFileSync(at, Buffer.alloc(9 * 1024 * 1024, 0x41));
    const ui = await loadFiles([H.file(at)]);
    assert.match(ui.values.loadStatus, /does not look like a text file/);
});

test('an unknown SNP is reported, not silently dropped', async () => {
    const ui = await loadFiles(
        [H.file(small('.bed')), H.file(small('.bim')), H.file(small('.fam'))],
        { snpListText: 'snp2\nrsNoSuchThing' });
    assert.match(ui.values.loadStatus, /loaded 1 of 2 SNPs.*1 not found/);
});

test('a list matching nothing clears the payload rather than sending all', async () => {
    for (const [fmt, files] of FORMATS) {
        const ui = await loadFiles(files(), { snpListText: 'rsNope1\nrsNope2' });
        assert.strictEqual(ui.values.genoContent, '',
                           fmt + ' sent a payload for an empty selection');
        assert.match(ui.values.loadStatus, /none of the .*requested SNPs|none of the requested SNPs/,
                     fmt + ': ' + ui.values.loadStatus);
    }
});

// ── the carriers must never describe the wrong file ─────────────────────────

test('picking a new set discards the previous payload', async () => {
    // sourceFormat follows the pick, so leaving the old bytes in place hands R
    // a .bed and tells it to read a VCF.
    const queue = [[H.file(small('.bed')), H.file(small('.bim')),
                    H.file(small('.fam'))],
                   [H.file(small('.vcf'))]];
    const mod = load(queue);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);
    assert.ok(ui.values.genoContent.length > 0, 'the .bed did not load');

    H.click(ui, 'snpi-geno');                   // now a VCF, not yet loaded
    assert.strictEqual(ui.values.sourceFormat, 'vcf');
    assert.strictEqual(ui.values.genoContent, '',
        'the .bed payload survived a re-pick — R would read it as a VCF');
    assert.strictEqual(ui.values.variantContent, '');
    assert.strictEqual(ui.values.sampleContent, '');
});

test('an incomplete set clears the payload and names what is missing', async () => {
    const mod = load([[H.file(small('.bed'))]]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-geno');
    assert.strictEqual(ui.values.genoContent, '');
    assert.match(ui.values.loadStatus, /small\.bim/);
});

test('a .fam that does not belong to the .bed is refused before slicing', async () => {
    // The one mismatch nothing downstream can catch. Variant i's offset is
    // 3 + i * ceil(N/4) with N from the .fam, and R re-derives the expected
    // payload length from the same N, so a stale .fam agrees with itself and
    // decodes a different variant's bytes into a confident genotype — measured
    // at ~30% of calls disagreeing with plink. Only the .bed's own size, which
    // neither of the other two files can fake, settles it.
    const fs = require('fs');
    const famLines = fs.readFileSync(small('.fam'), 'utf8')
                       .split('\n').filter(l => l.trim().length);
    // 196 rather than 199: ceil(196/4) is 49, so the stride actually moves and
    // the wrong bytes really are read. At 199 the stride is unchanged.
    const stale = new File([famLines.slice(0, 196).join('\n') + '\n'], 'small.fam');

    const mod = load([[H.file(small('.bed')), H.file(small('.bim')), stale]]);
    const ui = H.makeUi({ snpListText: 'snp0 snp1 snp2 snp3' });
    mod.view_loaded(ui);
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);

    assert.strictEqual(ui.values.genoContent, '');
    assert.match(ui.values.loadStatus, /not from the same dataset/);
    assert.strictEqual(ui.values.loadProblem, 'error');
});

test('a good trio carries its dimensions across for R to re-check', async () => {
    const fs = require('fs');
    const ui = await loadFiles([H.file(small('.bed')), H.file(small('.bim')),
                                H.file(small('.fam'))]);
    const nBim = fs.readFileSync(small('.bim'), 'utf8')
                   .split('\n').filter(l => l.trim().length).length;
    const nFam = fs.readFileSync(small('.fam'), 'utf8')
                   .split('\n').filter(l => l.trim().length).length;
    assert.strictEqual(ui.values.sourceDims,
        fs.statSync(small('.bed')).size + ',' + nBim + ',' + nFam);
});

// ── the other two slots ─────────────────────────────────────────────────────

test('the SNP list box and ✕ are hidden until a file is picked, then shown immediately', async () => {
    // Deliberately the opposite of what this test used to assert: the box+✕
    // now start hidden (matching covFile's native FileSelector, whose file
    // box is absent until something is picked) and must appear the instant
    // _loadSnpList sets a value, not on some later view_updated. A version of
    // this control that hides on empty but only re-checks from view_updated
    // would fail here exactly as it failed the old (inverted) version of
    // this test — see _syncClearable's own comment.
    const mod = load([[H.file(path.join(FX, 'small', 'select.txt'))]]);
    const ui = H.makeUi();
    mod.view_loaded(ui);

    let btn = H.find(ui.root, 'snpi-snps-clr')[0];
    assert.ok(btn, 'missing snpi-snps-clr');
    assert.strictEqual(btn.style.display, 'none',
        'snpi-snps-clr should start hidden, nothing is loaded yet');

    H.click(ui, 'snpi-snps');
    await H.until(() => ui.values.snpListFilename !== '', 'the SNP list');
    btn = H.find(ui.root, 'snpi-snps-clr')[0];
    assert.notStrictEqual(btn.style.display, 'none',
        'snpi-snps-clr should show as soon as a file is loaded');

    H.click(ui, 'snpi-snps-clr');
    assert.strictEqual(ui.values.snpListFilename, '');
    btn = H.find(ui.root, 'snpi-snps-clr')[0];
    assert.strictEqual(btn.style.display, 'none',
        'snpi-snps-clr should hide again once the file is cleared');
});

test('clearing the SNP list drops the payload it selected', async () => {
    const queue = [[H.file(path.join(FX, 'small', 'select.txt'))],
                   [H.file(small('.bed')), H.file(small('.bim')),
                    H.file(small('.fam'))]];
    const mod = load(queue);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-snps');
    await H.until(() => ui.values.snpListContent !== '', 'the SNP list');
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);
    assert.ok(ui.values.genoContent.length > 0);

    H.click(ui, 'snpi-snps-clr');
    assert.strictEqual(ui.values.genoContent, '',
        'the payload outlived the list that chose it');
    assert.match(ui.values.loadStatus, /press Load genotypes/);
});

test('the SNP list file loads and can be cleared', async () => {
    // Covariate browsing moved to the native FileSelector control and is no
    // longer exercised by this JS harness — see R/snpimport.b.R and
    // jamovi/snpimport.a.yaml for covFile.
    const mod = load([[H.file(path.join(FX, 'small', 'select.txt'))]]);
    const ui = H.makeUi();
    mod.view_loaded(ui);

    H.click(ui, 'snpi-snps');
    await H.until(() => ui.values.snpListFilename !== '', 'the SNP list');
    assert.strictEqual(ui.values.snpListFilename, 'select.txt');
    assert.ok(ui.values.snpListContent.length > 0);

    H.click(ui, 'snpi-snps-clr');
    assert.strictEqual(ui.values.snpListContent, '');
    assert.strictEqual(ui.values.snpListFilename, '');
});

test('a loaded SNP list selects the variants, without a pasted list', async () => {
    const queue = [[H.file(path.join(FX, 'small', 'select.txt'))],
                   [H.file(small('.bed')), H.file(small('.bim')),
                    H.file(small('.fam'))]];
    const mod = load(queue);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-snps');
    await H.until(() => ui.values.snpListContent !== '', 'the SNP list');
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);
    const n = require('fs').readFileSync(path.join(FX, 'small', 'select.txt'),
        'utf8').trim().split(/\s+/).length;
    assert.match(ui.values.loadStatus, new RegExp('loaded ' + n + ' of ' + n));
});

// ── the SNP list is size-checked like every other upload ────────────────────
//
// It was the one loader that sent unconditionally. Over the transport ceiling
// nanomsg drops the message with no error anywhere, so the failure to test for
// is silence: snpListContent set and the analysis looking fine.

// Write a SNP-list file of `n` ids into .tmp and hand back a File for it.
function snpListFile(name, n, pad) {
    const fs = require('fs');
    const tmp = path.join(__dirname, '..', '..', '.tmp');
    fs.mkdirSync(tmp, { recursive: true });
    const dst = path.join(tmp, name);
    const out = fs.createWriteStream(dst);
    for (let i = 0; i < n; i++)
        out.write('rs' + i + (pad ? '\t' + pad : '') + '\n');
    out.end();
    return new Promise((res, rej) => {
        out.on('finish', () => res(H.file(dst)));
        out.on('error', rej);
    });
}

async function pickSnpList(f) {
    const mod = load([[f]]);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-snps');
    await H.until(() => ui.values.snpListFilename !== '', 'the SNP list');
    return ui;
}

test('a SNP list is gzipped, not sent as plain base64', async () => {
    const ui = await pickSnpList(await snpListFile('list-small.txt', 20000));
    const bin = Buffer.from(ui.values.snpListContent, 'base64');
    assert.strictEqual(bin[0], 0x1f, 'not a gzip stream');
    assert.strictEqual(bin[1], 0x8b, 'not a gzip stream');
    assert.strictEqual(ui.values.snpListFilename, 'list-small.txt');
});

test('a SNP list past R\'s decompression cap is refused, not sent', async () => {
    // Compresses to almost nothing, so only the raw-size check can catch it.
    const ui = await pickSnpList(
        await snpListFile('list-huge.txt', 300000, 'A'.repeat(40)));
    assert.strictEqual(ui.values.snpListContent, '',
                       'an oversized list was sent anyway');
    assert.match(ui.values.snpListFilename, /too large \(1[0-9]\.[0-9] MB of text\)/);
});

test('a SNP list past the transport ceiling is refused, not sent', async () => {
    // Random ids, so gzip cannot get it under PAYLOAD_SAFE either.
    const fs = require('fs');
    const tmp = path.join(__dirname, '..', '..', '.tmp');
    fs.mkdirSync(tmp, { recursive: true });
    const dst = path.join(tmp, 'list-random.txt');
    const out = fs.createWriteStream(dst);
    for (let i = 0; i < 400000; i++)
        out.write(require('crypto').randomBytes(9).toString('hex') + '\n');
    out.end();
    await new Promise((r, j) => { out.on('finish', r); out.on('error', j); });

    const ui = await pickSnpList(H.file(dst));
    assert.strictEqual(ui.values.snpListContent, '');
    assert.match(ui.values.snpListFilename, /even compressed/);
});

test('a gzipped SNP list still selects the right variants', async () => {
    // The round trip: _loadSnpList gzips, _wantedIds gunzips, the .bed path
    // slices what came back. A gunzip that quietly returned mojibake would
    // match nothing and show up here rather than as an empty spreadsheet.
    const queue = [[H.file(path.join(FX, 'small', 'select.txt'))],
                   [H.file(small('.bed')), H.file(small('.bim')),
                    H.file(small('.fam'))]];
    const mod = load(queue);
    const ui = H.makeUi();
    mod.view_loaded(ui);
    H.click(ui, 'snpi-snps');
    await H.until(() => ui.values.snpListContent !== '', 'the SNP list');
    assert.strictEqual(Buffer.from(ui.values.snpListContent, 'base64')[0], 0x1f);
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);
    const n = require('fs').readFileSync(path.join(FX, 'small', 'select.txt'),
        'utf8').trim().split(/\s+/).length;
    assert.match(ui.values.loadStatus, new RegExp('loaded ' + n + ' of ' + n));
});

test('a .ped that does not match its .map is refused, not half-sent', async () => {
    const fs = require('fs');
    const tmp = path.join(__dirname, '..', '..', '.tmp');
    fs.mkdirSync(tmp, { recursive: true });
    const ped = fs.readFileSync(small('.ped'), 'utf8').trim().split('\n');
    ped[3] = ped[3].split(/[ \t]+/).slice(0, -4).join(' ');   // two calls short
    fs.writeFileSync(path.join(tmp, 'bad.ped'), ped.join('\n') + '\n');
    fs.copyFileSync(small('.map'), path.join(tmp, 'bad.map'));

    const ui = await loadFiles([H.file(path.join(tmp, 'bad.ped')),
                                H.file(path.join(tmp, 'bad.map'))]);
    assert.strictEqual(ui.values.genoContent, '');
    assert.match(ui.values.loadStatus, /do not have .* fields/);
});

test('a .bim with no variant IDs still sends every variant', async () => {
    // plink2 converting a VCF without --set-missing-var-ids writes '.' for
    // every ID. Selecting by distinct ID sent one variant out of 500 and
    // called it a successful load.
    const fs = require('fs');
    const tmp = path.join(__dirname, '..', '..', '.tmp');
    fs.mkdirSync(tmp, { recursive: true });
    const bim = fs.readFileSync(small('.bim'), 'utf8').trim().split('\n')
        .map(l => { const f = l.split('\t'); f[1] = '.'; return f.join('\t'); });
    fs.writeFileSync(path.join(tmp, 'dot.bim'), bim.join('\n') + '\n');
    fs.copyFileSync(small('.bed'), path.join(tmp, 'dot.bed'));
    fs.copyFileSync(small('.fam'), path.join(tmp, 'dot.fam'));

    const ui = await loadFiles(['bed', 'bim', 'fam'].map(
        e => H.file(path.join(tmp, 'dot.' + e))));
    assert.match(ui.values.loadStatus,
                 new RegExp('loaded ' + bim.length + ' SNPs'),
                 'only some variants were sent: ' + ui.values.loadStatus);
});

// ── the capacity wall ───────────────────────────────────────────────────────

test('an oversized selection is refused before anything is sent', async () => {
    // The whole point of the browser-side check: an oversized results message
    // is dropped with no error anywhere, so it has to be stopped here.
    needFixture(path.join(FX, 'large', 'large.bed'));
    const ui = await loadFiles(
        [H.file(path.join(FX, 'large', 'large.bed')),
         H.file(path.join(FX, 'large', 'large.bim')),
         H.file(path.join(FX, 'large', 'large.fam'))]);
    assert.strictEqual(ui.values.genoContent, '');
    assert.match(ui.values.loadStatus, /too much to transfer|too large/);
    // The refusal has to reach the results panel: the status box is one line
    // wide and R cannot tell a refusal from an import nobody started yet.
    assert.strictEqual(ui.values.loadProblem, 'error');
});

test('a .tped with no SNP list stops scanning instead of eating the file', async () => {
    // Every line matches when no list is given, so the whole 400 MB of
    // medium.tped used to be split, kept and joined in the tab before the cell
    // count was so much as looked at. 50 000 variants x 2 000 samples is 100M
    // cells against a 3M budget, so the answer is knowable after 1 500 lines.
    const tped = path.join(FX, 'medium', 'medium.tped');
    needFixture(tped);

    const t0 = Date.now();
    const ui = await loadFiles([H.file(tped),
                                H.file(path.join(FX, 'medium', 'medium.tfam'))]);
    assert.strictEqual(ui.values.genoContent, '');
    assert.match(ui.values.loadStatus, /before stopping/);
    assert.strictEqual(ui.values.loadProblem, 'error');
    // Reading all of it takes far longer than this; the assertion is that it
    // gave up rather than that it is fast.
    assert.ok(Date.now() - t0 < 60000,
              'took ' + (Date.now() - t0) + ' ms — did it read the whole file?');
});

test('a load that works leaves no problem behind for R to report', async () => {
    const ui = await loadFiles([H.file(small('.bed')), H.file(small('.bim')),
                                H.file(small('.fam'))]);
    assert.match(ui.values.loadStatus, /^loaded /);
    assert.strictEqual(ui.values.loadProblem, '');
});

main();
