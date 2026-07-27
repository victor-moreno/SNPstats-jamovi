# Response to the SNPstats Module Audit

**Module:** SNPstats — Basic Genetic Epidemiology Analyses of SNP data
**Audit reviewed:** `2026-07-27 SNPstats.md` (Claudia)
**Version before:** 0.5.0 → **after:** 1.0.0
**Date:** 2026-07-27

---

## Summary

All 17 findings were investigated:

- **11 fixed** outright.
- **1 partly fixed by agreement** — `genetics` moved out of `Depends:` now, its
  removal planned with an agreed approach.
- **2 deliberately not done** — the long-function refactor and translation, both
  optional and both large; recorded in `NEWS.md` rather than silently dropped.
- **3 turned out to be incorrect as stated.** In each case the underlying
  mechanism was verified against the jamovi compiler or jmvcore rather than
  assumed; the evidence is below. One of them matters: acting on it as written
  would have broken a working feature.

The full test suite (`bash tests/run_tests.sh`) is green after every phase:
**9 files, 474 assertions, 0 failures**. Five new tests were added — four
pinning the security properties of the weights-file change (two of which
exercise the new public helper) and one pinning the missingness-plot state fix.

Work landed in five commits:

| Commit | Scope |
|---|---|
| `29d25e9` | Dependencies, metadata, labels |
| `a28b960` | **CRITICAL** — weights-file path removal |
| `da72bce` | Status messages, `.init` read, plot state |
| `ebb822b` | Version 1.0.0, NEWS/CLAUDE.md |
| `c57b09d` | This document, tutorial update |

---

## Findings status at a glance

| # | Severity | Finding | Outcome |
|---|---|---|---|
| 1 | CRITICAL | `.omv` can read arbitrary files via `weightsPath` | **Fixed** — option removed |
| 2 | HIGH | Version below library minimum | **Fixed** — 1.0.0 |
| 3 | HIGH | `ggplot2` / `base64enc` undeclared | **Fixed** |
| 4 | HIGH | Browse button never wired up | **Not a defect** — it was wired; see below |
| 5 | MEDIUM | `menuSubgroup` missing from snpStats | **Fixed** |
| 6 | MEDIUM | Haplotype interaction columns built in `.run()` | **Not the stated defect** — see below |
| 7 | MEDIUM | `genetics` is a `Depends:` | **Partly fixed** — moved to `Imports:`; removal planned |
| 8 | MEDIUM | Hand-styled HTML with hardcoded colours | **Fixed** — but not via `type: Notice`, which is unavailable |
| 9 | LOW | `.init()` reads the dataset six times | **Fixed** |
| 10 | LOW | Checkbox labels start with "Show" | **Fixed** |
| 11 | LOW | Section-heading casing inconsistent | **Fixed** |
| 12 | LOW | Missingness plot renders from a private field | **Fixed** |
| 13 | LOW | Several functions are very long | **Deliberately not done** — see below |
| 14 | LOW | Two analyses declare different `version:` | **Fixed** |
| 15 | INFO | No strings marked for translation | **Not done** — out of scope, recorded in NEWS |
| 16 | INFO | `References:` not a `DESCRIPTION` field | **Fixed** — replaced by `inst/CITATION` |
| 17 | INFO | `caseLevel` in NEWS.md looks resolved | **Not resolved** — NEWS.md corrected; see below |

---

## [CRITICAL] The weights-file path — removed

### The threat model, confirmed

Before fixing this, one premise was worth testing: *does an attacker actually
need the UI to produce a malicious `.omv`?* No. An `.omv` is a zip archive and
the analysis options are data inside it — anyone can unzip, edit the stored
option, and rezip. Hiding a control does not help either: `hidden:` affects
rendering only, and the option remains part of the serialised analysis and is
replayed on open. The audit's reading was correct.

### What was considered

A middle route was investigated and rejected. jmvcore's `Analysis` tracks
whether the caller handed it a data frame in-process: `private$.dataProvided`
defaults to `TRUE` and is flipped to `FALSE` inside `init()`/`run()` whenever
`private$.data` is `NULL`, which is exactly the engine path. Gating the path
read on that flag would have closed the document-open hole while keeping
`snpPGS(data = df, weightsPath = "x.csv")` working from R. It works, but it
leans on a jmvcore private field, and the decision was to take the stronger,
simpler route instead.

### What was done

**`weightsPath` is gone entirely.** There is no path-reading code left in the
module: no `file.exists()`, no `readLines()`, no `file.info()`.

- `jamovi/snpPGS.a.yaml` — option deleted. `weightsFilename` is now visible (it
  backs the display box); `weightsContent` documents *why* there is no path
  option so the next person does not helpfully add one back.
- `jamovi/snpPGS.u.yaml` — the free-text TextBox becomes a **read-only display
  of the chosen file's name**; `jamovi/js/snpPGS.js` sets `readonly` on the
  input and anchors the 📁 button to it.
- `R/snpPGS.b.R` — `.hasWeights`, `.weightsRawLines`, `.weightsSig` and
  `.weightsLabel` lost their path branches.
- `jamovi/snpPGS.r.yaml` — the 13 now-dead `weightsPath` `clearWith` entries
  removed.
- `jamovi/js/snpPGS.js` — the `reader.onerror` handler no longer writes
  `file.path || file.name` into an option; it reports the read failure.

### Keeping R usable — the `pgs_weights()` helper

Removing the option would otherwise have made R scripting awkward, so the module
now exports a helper (`R/pgs_weights.R`):

```r
w <- pgs_weights("pgs_catalog.csv")
do.call(snpPGS, c(list(data = mydata, snpCols = c("rs1", "rs2")), w))
```

It reads the file **in the caller's own session** — the same act as `read.csv()`
— and returns the `weightsContent` / `weightsFilename` pair. Nothing about the
path enters the analysis, so a saved `.omv` carries the weights themselves.

### The other hardening items

- **Gzip bomb.** `memDecompress(bytes, "gzip")` had no size cap. Both
  `.weightsRawLines()` and `pgs_weights()` now enforce
  `PGS_MAX_WEIGHTS_BYTES` (64 MB); an oversized payload is refused rather than
  expanded. A test compresses 70 MB into a <1 MB payload and asserts the
  analysis survives.
- **No content echo on parse failure.** A weights file with no recognisable
  rsID column used to render its header line into the results. It now reports
  the *expected* column names only.
- **The SNP Weights grid was deliberately left showing all catalog rows.** With
  the path gone, the file can only be one the user picked in-session or bytes
  the sender already possessed — so whole-file rendering is no longer a
  disclosure channel, and capping it would remove a feature for no security
  gain. The gzip cap bounds the pathological case.

### Tests added

```
snpPGS exposes no file-path option
pgs_weights round-trips a file into the content/name pair
a parse failure reports expectations, not file content
an oversized gunzipped payload is refused rather than expanded
```

The existing ~30 PGS test call sites were kept intact by giving the `run_pgs()`
test helper a `weightsFile =` shim that expands through the public
`pgs_weights()`, so the suite exercises the new helper on every run.

---

## [HIGH] Findings 2, 3 — version and dependencies

**Version → 1.0.0** in all four places that carry it: `DESCRIPTION`,
`jamovi/0000.yaml`, `jamovi/snpStats.a.yaml` and `jamovi/snpPGS.a.yaml` (the
last was still on `0.1.0`, finding #14). Applied last, on a clean tree.

**Dependencies.** `ggplot2` and `base64enc` added to `Imports:`, with
`importFrom(base64enc, base64decode, base64encode)` in `NAMESPACE`. This
required installing `ggplot2` into the project-local test library;
`tests/setup_test_env.sh` was updated so a fresh setup gets it. CI reads
`DESCRIPTION` automatically and needed no change.

Note this reverses a decision `CLAUDE.md` had recorded as intentional
("supplied by the jamovi runtime … intentionally not declared"). The audit is
right that R guarantees only what is in `Imports:`, and jamovi's own `jmv`
module lists `ggplot2` under `Imports:`. `CLAUDE.md` was updated to match.

---

## [HIGH] Finding 4 — the browse button *was* wired up

**This finding is incorrect**, and the correction matters because acting on it
as written would have broken a working feature.

jamovi has two conventions for view handlers, and they are not interchangeable:

1. **`jamovi/js/<analysis>.js`** — auto-loaded. The UI compiler emits
   `this.handlers = require('./<analysis>')` into the generated `.src.js` purely
   because the file is named after the analysis. `view_loaded` / `view_updated`
   resolve against it. This is what jmv's own `anova.js` uses.
2. **`jamovi/js/<analysis>.events.js`** — loaded **only** when the `.u.yaml`
   declares an explicit `events:` block with `./<file>::<export>` references,
   as jmv's `ancova.u.yaml` does.

The audit assumed (2) was the auto-loaded one and that `snpPGS.js` was dead.
It is the other way round. Verified directly:

```
$ grep -o "this.handlers = [^;]*" build/js/snpPGS.src.js
this.handlers = require('./snpPGS')          # with jamovi/js/snpPGS.js present
this.handlers = { }                          # with it absent
```

and against the jamovi-shipped `jmv` build:

```
$ grep -o "this.handlers = [^;]*" .../modules/jmv/ui/anova.js
this.handlers = require('./anova')
```

So the 📁 browse button, `weightsContent` and the whole cloud-compatible path
were live all along. The empty `snpPGS.events.js` was the file that did nothing.

**Action taken:** the empty `snpPGS.events.js` was **deleted** — it was the
thing that made the live file look dead — and `snpPGS.js` gained a header
comment explaining the naming rule. The convention is now also recorded in the
`CLAUDE.md` file layout so the misdiagnosis is not repeated.

---

## [MEDIUM] Finding 6 — haplotype interaction tables: analysed, not changed

The proposed fix (pre-create the per-level columns in `.init()`) would not
change what users see, and would introduce a new failure mode.

The three tables get **both** their columns *and their rows* in `.run()`. Rows
must be there, because the haplotype set is unknowable until `haplo.em` has run.
Consequently:

- after restore, `rowCount == 0`;
- `.need_fill(tbl)` is `rowCount == 0 || isNotFilled()`, so it is **always
  `TRUE`**;
- `.compute_haplo_interaction()` therefore runs on every click and rebuilds the
  tables completely.

They do not blank — they **recompute**. Pre-creating the columns in `.init()`
would additionally make `.run()`'s six unguarded `addColumn()` calls collide
with the pre-created ones, which would need guarding for no user-visible gain.

The real cost here is a repeated `haplo.glm` fit on every unrelated option
click. That is a genuine performance issue and is recorded in `NEWS.md` under
Plan, alongside the parallel speed-up item.

---

## [MEDIUM] Finding 8 — colours fixed; `type: Notice` is not available

The substantive complaint is right: 17 sites across `snpStats.b.R`,
`snpPGS.b.R` and `snp_compute.R` hardcoded `red` / `orange` / `#c0392b` /
`#555` / `#666`, which are unreadable against jamovi's dark theme.

The suggested remedy is not available in this jamovi version. Changing
`validationMsg` to `type: Notice` fails the compiler outright:

```
Unable to compile 'snpStats.r.yaml':
	results.items[1].type is not one of enum values: Table,Group,Array,Image,
	Preformatted,Html,State,Property,Output,Notification,Action
```

`Notification` is in the enum but jmvcore ships no `Notification` class, so it
is not usable from R either. jamovi's own `jmv` reaches Notices the only
remaining way — constructing `jmvcore::Notice$new(...)` in R and calling
`self$results$insert(i, notice)`.

**What was done instead:** two helpers in `snp_compute.R`, `msg_warn()` and
`msg_info()`, that emit no colour at all and let the platform pick one. All 17
sites converted, plus the `#666` on both `helpBanner`s and the `#555` on the
PGS "Getting started" block. This fully resolves the dark-theme defect.

**Why not the programmatic route now:** it would replace a declared,
`clearWith`-managed, restore-aware results element with elements inserted
dynamically in `.run()` — precisely the pattern this module documents at length
as not surviving jamovi's restore. Recorded in `NEWS.md` as a follow-up.

---

## [MEDIUM] Finding 7 — `genetics`

Done now: `genetics` and `haplo.stats` moved from `Depends:` to `Imports:`.
Zero code impact — every call already goes through `::` or `importFrom()`. The
module no longer attaches anything to the user's search path.

**Flagged for the next round, as agreed.** The plan recorded in `NEWS.md` is to
**vendor** the four functions actually used — `genotype`, `allele`,
`HWE.exact`, `LD` — into an auxiliary R file, so the functionality survives if
the package leaves CRAN, and then drop the `Imports:` entry. This is preferable
to reimplementing from scratch: it preserves numeric behaviour, so the golden
values stay meaningful. They will still need re-verification against current
outputs.

---

## LOW findings

**9 — `.init()` dataset reads.** `.init_data()` is now memoised per analysis
object (`.init_dat` + an `.init_dat_read` flag so a genuine `NULL` read is not
retried). Every caller is inside `.init()`, and jamovi builds a fresh object per
click, so the cache cannot go stale. Six full dataset reads per option click
become one. Applied to `snpPGS.b.R` too, which called it twice.

**10 — checkbox labels.** Seven titles reworded in `jamovi/snpStats.a.yaml`,
including the two the audit called marginal: "Missing values in tables", "SNP
missingness plot", "Interaction table", "Adjustment covariate parameters",
"Cross-classification table", "Stratified by response groups", "SNPs with
missing values removed".

**11 — heading casing.** Eight group headings title-cased in
`snpStats.u.yaml`; "Show in tables" → "Table Options"; "SNP QC: exclude SNPs" →
"SNP QC: Exclude SNPs" in `snpPGS.u.yaml`.

**12 — missingness plot.** Now uses `image$setState()` / `image$state`, matching
`.render_ld_plot`; the `.miss_cache` private field is gone. Verified by test
that a fresh object which never ran renders correctly from state alone — the
case that previously produced a blank panel. The five PGS plots were left on
`private$.cache` as the audit advised, since their cached object includes the
dosage matrix.

**13 — long functions: deliberately not done.** Seven functions over 200 lines
is real, but the audit itself says "no rush", and splitting well-tested
numerical code produces a large diff whose only benefit is future
maintainability. Doing it in the same pass as a security fix would make the
security change much harder to review. It is not recorded as a defect because
nothing is wrong with the code — it is a refactor worth doing before the next
feature lands in one of them.

---

## INFO findings

**16 — `References:`** deleted from `DESCRIPTION` (not an R field; contained
HTML and a non-ASCII en dash). Replaced with a proper `inst/CITATION`, so
`citation("SNPstats")` now works. `jamovi/00refs.yaml` was already correct.

**17 — `caseLevel`: NEWS.md was right, the audit was wrong to call it resolved.**
The generated header does contain `caseLevel = NULL`, but it is **not** what the
compiler emits. Verified by running `jmvtools::prepare()` on a clean tree: the
compiler emits a bare `caseLevel,` every time, and the whole suite then fails
with `argument "caseLevel" is missing`. The `NULL` is re-applied by a scripted
patch. The audit's alternative — adding `default: ''` to the option — is exactly
what the compiler rejects for a `Level` type.

This bit twice during this work, because the patch lived only inside
`install_jamovi.sh` and a bare `prepare()` bypassed it. It is now extracted to
**`tools/patch_h.sh`** and called by both `install_jamovi.sh` and
`tests/run_tests.sh`, so the tree self-heals. `NEWS.md` was corrected to say the
issue is *worked around*, not fixed.

**15 — translation.** Not done; genuinely optional and a large retrofit.
Recorded in `NEWS.md` under Plan.

---

## Incidental fixes

Two small things found while working, neither in the audit:

- **`.Rbuildignore` was missing `^\.Rlib-arm$`** (it listed `.Rlib` and
  `.Rlib-x64` only), so the arm64 project library would have been rolled into
  `R CMD build` output.
- **`docs/TUTORIAL.md`** described typing a path into the Weights file field.
  Rewritten around the 📁 button and `pgs_weights()`, with a short note on why
  the module stores contents rather than a location.

## Observation, not fixed

When a weights file fails to parse, the informative message ("no recognisable
rsID column") is immediately overwritten by a downstream one ("No SNPs passed QC
filters") from `.applyDosageFilter`, because the fallback unit-weight table has
no allele information. This is pre-existing and unrelated to the security
change, so it was left alone — but it means a user with a malformed weights file
gets a misleading diagnosis. Worth a look.

---

## Verification

```
$ bash tools/install_jamovi.sh
wrote module: SNPstats_1.0.0.jmo
Module installed successfully

$ bash tests/run_tests.sh
association:      177 ✓
descriptive:       42 ✓
edgecases:         21 ✓
golden-external:   23 ✓
golden:            38 ✓
ldhaplo:           76 ✓
pgs:               97 ✓
refresh-pgs:        S   (RProtoBuf not installed on arm64)
refresh:            S   (RProtoBuf not installed on arm64)
DONE — 0 failures
```

The two refresh suites skip on this machine for want of `RProtoBuf`; they run on
the x86_64 machine and in CI. **They should be run there before release**, since
the `weightsContent`-only option set and the new `.init_data` memoisation both
touch the restore path they exercise.

Also worth doing by hand before release, as the audit suggested: install into a
clean jamovi and confirm the 📁 button appears next to the (now read-only)
Weights file box and that picking a file populates the analysis. The wiring is
verified at the compiler level here, but not clicked.
