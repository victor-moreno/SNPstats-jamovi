History:
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
  `tools/install_jamovi.sh` and `tests/run_tests.sh` call after any rebuild.
  A bare `jmvtools::prepare()` on its own still leaves the tree broken.

Plan:
- parallel speed-up (LD / association).
- Haplotype tables recompute `haplo.em` on every option click
  (~1.6 s for 4 SNPs) because their rows cannot be predicted in `.init()`. 
- String translation.
