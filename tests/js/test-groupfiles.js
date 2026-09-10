'use strict';

// Format dispatch, tested without a browser.
//
// Written after shipping the .ped/.tped/VCF readers with the file dialog still
// filtering for '.bed,.bim,.fam': every reader worked, every R test passed, and
// the formats were unreachable because the picker would not show the files.
// Nothing in the R suite could have caught that — this is the layer where it
// lives.
//
//   node tests/js/test-groupfiles.js

const fs = require('fs');
const path = require('path');
const assert = require('assert');

const src = fs.readFileSync(
  path.join(__dirname, '..', '..', 'jamovi', 'js', 'snpimport.js'), 'utf8');

// The module is written for the jamovi UI, so it is not requireable here. Pull
// out the pure functions and the accept list, which is all this file needs.
function extract(name) {
  const at = src.indexOf('function ' + name + '(');
  assert.ok(at >= 0, 'function not found: ' + name);
  let depth = 0, i = src.indexOf('{', at);
  const start = at;
  for (; i < src.length; i++) {
    if (src[i] === '{') depth++;
    else if (src[i] === '}') { depth--; if (depth === 0) break; }
  }
  return src.slice(start, i + 1);
}

const sandbox = {};
new Function('exports', extract('_groupFiles') + '\nexports._groupFiles = _groupFiles;')(sandbox);
const groupFiles = sandbox._groupFiles;

const f = (name) => ({ name, size: 1024 });
let passed = 0, failed = 0;

function check(what, fn) {
  try { fn(); passed++; }
  catch (e) { failed++; console.error('FAIL: ' + what + '\n      ' + e.message); }
}

// ── the accept list must cover every format the readers handle ───────────────

check('the file dialog offers every readable extension', () => {
  // third argument of the call: (ui, name, cssClass, accept, ...)
  const m = /_browseButton\(ui,\s*'genoFilename',[\s\S]{0,300}?'[^']*',\s*'([^']*)'/
    .exec(src);
  assert.ok(m, 'could not find the genotype browse button');
  const accept = m[1];
  for (const ext of ['.bed', '.bim', '.fam', '.ped', '.map',
                     '.tped', '.tfam', '.vcf', '.gz'])
    assert.ok(accept.includes(ext),
      'accept list is missing ' + ext + ' — the picker will hide those files');
});

check('the removable file slots get a clear button', () => {
  // The filename fields are read-only by design, so without a clear button a
  // wrongly chosen file can only be replaced, never removed.
  assert.ok(src.includes("['snpListContent', 'snpListFilename']"),
    'the SNP list slot has no clear list, so its file cannot be removed');
  assert.ok(src.includes("['covContent', 'covFilename']"),
    'the covariate slot has no clear list, so its file cannot be removed');

  // and _browseButton must actually honour that argument
  assert.ok(/if \(clears && clears\.length\)/.test(src),
    '_browseButton ignores its clears argument');
});

check('no injection can take the Load button down with it', () => {
  // The regression this pins: _browseButton wrote a property onto a jamovi
  // control object, which throws under 'use strict' when the object is sealed.
  // _inject built the Load button last, so that throw left the panel with no
  // way to load anything — the browse buttons injected earlier still worked,
  // which made it look like Load itself had broken.
  // Nothing may be stashed on a jamovi object at all — control *or* view. The
  // second occurrence of this bug was `ui.__snpImportFiles`, which made picking
  // files fail silently and left Load with nothing to do.
  assert.ok(!/\b(ctrl|ui)\.__snpi/i.test(src),
    'nothing may be assigned to a jamovi control or view object');
  assert.ok(/new WeakMap\(\)/.test(src),
    'the picked files should live in a WeakMap keyed by a DOM node');

  const inject = /function _inject\(ui\)[\s\S]*?\n}/.exec(src);
  assert.ok(inject, '_inject not found');
  const body = inject[0];

  for (const call of ['_browseButton', '_loadButton'])
    assert.ok(body.includes(call), '_inject no longer calls ' + call);

  // Every injection must sit inside its own _guard. Counting rather than
  // pattern-matching per line: the calls legitimately appear indented *within*
  // the guard closures, so a line-shape test cannot tell guarded from bare.
  const guards = (body.match(/_guard\(/g) || []).length;
  const calls  = (body.match(/_(browseButton|loadButton)\(/g) || []).length;
  assert.ok(guards >= calls,
    `${calls} injection calls but only ${guards} guards — one throw kills the rest`);
});

check('click handlers report failures instead of doing nothing', () => {
  // A throw inside a click handler is indistinguishable from a dead button, so
  // every handler routes through _guardUI, which writes to the status line.
  assert.ok(/function _guardUI\(/.test(src), '_guardUI is missing');
  assert.ok(/_guardUI\(ui, 'load'/.test(src), 'the Load click is unguarded');
  assert.ok(/_guardUI\(ui, 'choosing files'/.test(src),
    'the file-picker change handler is unguarded');
});

// ── dispatch ────────────────────────────────────────────────────────────────

check('a complete .bed trio is recognised', () => {
  const g = groupFiles([f('x.bed'), f('x.bim'), f('x.fam')]);
  assert.strictEqual(g.format, 'bed');
  assert.ok(!g.error);
});

check('a .ped pair is recognised', () => {
  const g = groupFiles([f('x.ped'), f('x.map')]);
  assert.strictEqual(g.format, 'ped');
});

check('a .tped pair is recognised', () => {
  const g = groupFiles([f('x.tped'), f('x.tfam')]);
  assert.strictEqual(g.format, 'tped');
});

check('a lone VCF is recognised, plain or gzipped', () => {
  assert.strictEqual(groupFiles([f('x.vcf')]).format, 'vcf');
  const gz = groupFiles([f('x.vcf.gz')]);
  assert.strictEqual(gz.format, 'vcf');
  assert.strictEqual(gz.gz, true);
});

check('extensions are matched case-insensitively', () => {
  assert.strictEqual(groupFiles([f('X.BED'), f('X.BIM'), f('X.FAM')]).format, 'bed');
  assert.strictEqual(groupFiles([f('X.VCF')]).format, 'vcf');
});

// ── the errors have to name what is missing ─────────────────────────────────

check('an incomplete set names the files still needed', () => {
  const g = groupFiles([f('x.bed')]);
  assert.ok(g.error, 'expected an error');
  assert.ok(g.error.includes('x.bim') && g.error.includes('x.fam'),
            'error should name the missing files: ' + g.error);

  const p = groupFiles([f('y.ped')]);
  assert.ok(p.error.includes('y.map'), p.error);
});

check('files from different datasets are refused, not mixed', () => {
  const g = groupFiles([f('a.bed'), f('a.bim'), f('b.fam')]);
  assert.ok(g.error && /different datasets/.test(g.error), g.error);
});

check('a VCF alongside other files is refused', () => {
  const g = groupFiles([f('x.vcf'), f('x.bed')]);
  assert.ok(g.error && /on its own/.test(g.error), g.error);
});

check('two VCFs are refused', () => {
  const g = groupFiles([f('a.vcf'), f('b.vcf')]);
  assert.ok(g.error && /one VCF/.test(g.error), g.error);
});

check('nothing recognisable gives a useful message', () => {
  const g = groupFiles([f('notes.txt')]);
  assert.ok(g.error && /\.vcf/.test(g.error), g.error);
});

check('the SNP-ID header names match R, in the same order', () => {
    // Both halves read the user's selection file: the browser to decide which
    // variants to slice, R to decide what was requested and what the effect
    // alleles are. They must reach the same column. When R matched by position
    // and this matched by preference, a header of `id, rsid, effect_allele`
    // had the browser slicing the rsIDs while the report listed the contents of
    // `id` as requested and every one of them as not found.
    const js = /var SNP_ID_NAMES = \[([^\]]*)\]/.exec(src);
    assert.ok(js, 'SNP_ID_NAMES not found in snpimport.js');

    const rsrc = fs.readFileSync(
        path.join(__dirname, '..', '..', 'R', 'payload.R'), 'utf8');
    const r = /SNP_ID_NAMES <- c\(([^)]*)\)/.exec(rsrc);
    assert.ok(r, 'SNP_ID_NAMES not found in payload.R');

    const norm = (s) => s.split(',').map(x => x.trim().replace(/^['"]|['"]$/g, ''))
                         .filter(x => x.length);
    assert.deepStrictEqual(norm(js[1]), norm(r[1]));
});

console.log(`\n${passed} pass, ${failed} fail`);
process.exit(failed > 0 ? 1 : 0);
