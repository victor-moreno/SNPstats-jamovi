History:
- 260420 v0.2.0 First public release
- 260501 v0.3.0 Added PGS submodule
- 260508 v0.4.0 Combined submenus for SNPstats and categorical response
- 260717 v0.5.0 Refactored to eliminate table refresh
- 260727 v1.0.0 Audit fixes (see 2026-07-27 SNPstats-response.md)
  - SECURITY: the `weightsPath` option is gone. A weights file path was saved
    into the .omv and re-resolved when the file was opened, so a crafted
    analysis could read an arbitrary file from a collaborator's machine. Weights
    now travel as embedded content only — the file-browse button in the UI, or
    the new exported `pgs_weights()` helper from R.
  - `ggplot2` and `base64enc` declared in Imports; `genetics` and `haplo.stats`
    moved from Depends to Imports.
  - Status messages no longer hardcode colours (unreadable in dark theme).
  - `.init()` reads the dataset once per option click instead of six times.
  - Missingness plot renders from image state instead of a private field.
  - `menuSubgroup` added to snpStats so both analyses sit together in the menu.

Issues:
- snpPGS gets argument caseLevel without default: a `type: Level` option cannot
  carry a yaml `default:` (the compiler rejects it), so every run of the jamovi
  UI compiler re-emits a bare `caseLevel,` formal on the public snpPGS()
  function and R-side calls that omit it fail. NOT resolved at the source —
  it is re-patched automatically by `tools/patch_h.sh`, which both
  `tools/install_jamovi.sh` and `tests/run_tests.sh` call after any rebuild.
  A bare `jmvtools::prepare()` on its own still leaves the tree broken.

Plan:
- remove genetics dependency because package is marked obsolete. Approach
  agreed 2026-07-27: vendor the four functions actually used
  (`genotype`, `allele`, `HWE.exact`, `LD`) into an auxiliary R file so the
  functionality survives if the package leaves CRAN, then drop the Imports
  entry. Golden values must be re-verified against the current outputs.
- parallel speed-up (LD / association / haplotype).
- Themed status messages via `jmvcore::Notice` + `results$insert()`; `type:
  Notice` is not accepted in the .r.yaml by this compiler version.
- Haplotype interaction tables recompute `haplo.glm` on every option click
  because their rows cannot be predicted in `.init()`.
- String translation (`jmvcore::.()`) if Spanish/Catalan versions are wanted.
