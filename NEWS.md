History:
- 260930 v1.2.1 jamovi 28.4 compatibility. jamovi 28.4 renamed the .u.yaml
                control of a File option from FileSelector to File
                (upstream ac3126b1). 1.2.0 was built with the old name, which
                jamovi 28.4's client no longer knows, so the options panels
                of snpPGS (weights file) and snpImport (covariate file) cannot
                be built there; jmvtools 28.4 also refuses to compile it
                ("Cannot read properties of undefined (reading
                'isOptionControl')"). Both controls are renamed; nothing else
                changes (the regenerated .h.R files differ only in the
                analysis version). minApp raised 28.3.0 -> 28.4.0: no single
                .jmo serves both 28.3 and 28.4, and jmc 28.4 will not install
                into jamovi 28.0-28.3 anyway. On jamovi 28.3, stay on 1.2.0.
                Panel JS: the file boxes and buttons drawn by snpPGS.js and
                snpimport.js use translucent greys instead of hex greys
                (identical on jamovi's white panel, legible on any darker
                background) and buttons inherit the text colour.
- 260920 v1.2.0 snpPGS's weights file and snpImport's covariate file now use
                jamovi 28.3's native File option (FileSelector control)
                instead of a custom browse button that embedded base64
                content into hidden String options. jamovi carries the picked
                file as a resource inside the saved .omv, so a reopened
                analysis still works without the original file present --
                verified end to end (pick, save, quit jamovi, reopen) before
                this change was made. snpImport's SNP-list file and genotype
                trio (.bed/.bim/.fam etc.) are UNCHANGED: the genotype-trio
                slicing JS (jamovi/js/snpimport.js) reads the SNP list's raw
                text client-side before Load is pressed, which a native
                FileSelector cannot support (it never exposes file bytes to
                analysis JS, only {path, filename}), so both stay on the
                existing mechanism. jamovi/js/snpPGS.js (183 lines) is
                deleted entirely; jamovi/js/snpimport.js loses only its
                covariate-browsing code.
                minApp raised 28.1.0 -> 28.3.0 (the File option type needs
                jamovi 28.3). Needs jmvcore 28.3 (on CRAN since 2026-09-24;
                CRAN's earlier 2.7.38 has no OptionFile, which is what turned
                CI red on this release until then) and builds with the
                released jmvtools 28.3 -- re-verified 2026-09-27: full test
                suite green in CI and a clean .jmo build.
                A .omv saved before this change carries the old
                weightsContent/weightsFilename or covContent/covFilename
                values, which are simply unknown options now; re-pick the
                weights/covariate file once after opening it in this version.
- 260420 v0.2.0 First public release
- 260501 v0.3.0 Added PGS submodule
- 260508 v0.4.0 Combined submenus for SNPstats and categorical response
- 260717 v0.5.0 Refactored to eliminate table refresh
- 260727 v1.0.0 Audit fixes
- 260910 v1.1.0 Third analysis: Import genotypes (menu SNPstats > Data), merged
                in from the standalone snpImport module after its jamovi review.
                Reads PLINK .bed/.bim/.fam, .ped/.map, .tped/.tfam and VCF
                .vcf/.vcf.gz for a chosen SNP list, merges covariates by sample
                ID, applies QC filters and opens the result as a new dataset.
                Also in this release (previously unreleased): PGS plots export
                correctly (were blank); all plots follow the jamovi theme;
                snpStats declares weightsSupport: none.
                minApp raised 1.0.8 -> 28.1.0. The old value was a compiler
                default, never a tested claim; the import analysis genuinely
                needs jamovi 28.1 server behaviour (the Action option and
                Output columns), so pre-28.1 jamovi is no longer offered this
                module.

Issues:
- NAMESPACE must keep `import(jmvcore, except = format)`. jmvcore exports its
  own format(str, ..., context) for {} interpolation, which is not
  base::format(x, big.mark=, scientific=); the import analysis uses the base
  one two dozen times. Under a bare import(jmvcore) every import fails with
  "missing value where TRUE/FALSE needed" and the panel just says "Import
  failed" — the message never names format, so this is expensive to rediscover.
- snpPGS gets argument caseLevel without default: a `type: Level` option cannot
  carry a yaml `default:` (the compiler rejects it), so every run of the jamovi
  UI compiler re-emits a bare `caseLevel,` formal on the public snpPGS()
  function and R-side calls that omit it fail. NOT resolved at the source —
  it is re-patched automatically by `tools/patch_h.sh`, which both
  `tools/install.sh` and `tests/run_tests.sh` call after any rebuild.
  A bare `jmvtools::prepare()` on its own still leaves the tree broken.

Plan:
- parallel speed-up (LD / association).
- Haplotype tables recompute `haplo.em` on every option click
  (~1.6 s for 4 SNPs) because their rows cannot be predicted in `.init()`. 
- String translation.
