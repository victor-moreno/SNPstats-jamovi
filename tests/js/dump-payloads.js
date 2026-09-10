'use strict';

// Produce, for each supported format, exactly the options the browser would
// send — and write them where the R suite can read them.
//
//   node tests/js/dump-payloads.js <tier> <out.tsv> [id ...]
//
// The two halves are otherwise tested against each other only by hand: the R
// tests build their own payloads, so a browser that slices the wrong bytes
// still passes them. tests/testthat/test-browser-payloads.R reads this file and
// runs the analysis on it, which closes that gap.

const fs = require('fs');
const path = require('path');
const H = require('./harness');

const tier = process.argv[2] || 'small';
const out  = process.argv[3];
const ids  = process.argv.slice(4);

const FX = path.join(__dirname, '..', '..', 'data-raw', 'fixtures', tier);
const f  = (ext) => path.join(FX, tier + ext);

const SETS = {
    bed:  ['.bed', '.bim', '.fam'],
    tped: ['.tped', '.tfam'],
    ped:  ['.ped', '.map'],
    vcf:  ['.vcf'],
    vcfgz: ['.vcf.gz']
};

// Tab-separated rather than JSON: base64 contains neither tab nor newline, and
// the R side then needs no parser and no extra package.
const CARRIERS = ['sourceFormat', 'genoContent', 'variantContent',
                  'sampleContent', 'loadStatus'];

async function one(files, opts) {
    H.installDocument([files.map(e => H.file(f(e)))]);
    delete require.cache[require.resolve('../../jamovi/js/snpimport.js')];
    const mod = require('../../jamovi/js/snpimport.js');
    const ui = H.makeUi(opts);
    mod.view_loaded(ui);
    H.click(ui, 'snpi-geno');
    H.click(ui, 'snpi-load');
    await H.settled(ui);
    const o = {};
    CARRIERS.forEach(k => { o[k] = ui.values[k]; });
    return o;
}

(async () => {
    const opts = ids.length ? { snpListText: ids.join('\n') } : {};
    const rows = [['set'].concat(CARRIERS).join('\t')];
    for (const [name, files] of Object.entries(SETS)) {
        const o = await one(files, opts);
        rows.push([name].concat(CARRIERS.map(k => o[k])).join('\t'));
    }
    const tsv = rows.join('\n') + '\n';
    if (out) fs.writeFileSync(out, tsv);
    else process.stdout.write(tsv);
})();
