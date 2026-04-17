# SNPstats Module — Refactoring Guide
## Monolith → Three-Analysis Architecture

---

## What changed and why

The original module had a single `snpAnalysis` analysis that handled everything:
descriptives, association, LD, and haplotypes. At ~2500 lines of R and four large
YAML files it became hard to isolate bugs and test changes without risk of breaking
unrelated features.

The refactored module splits this into **three independent analyses** that share a
single helper library:

| Menu item                 | Files                                    | Responsibility                                      |
|---------------------------|------------------------------------------|-----------------------------------------------------|
| SNP Descriptive           | `snpDesc.{a,r,u}.yaml` / `snpDesc.b.R`  | Frequencies, HWE, covariate table, SNP summary      |
| SNP Association Analysis  | `snpAssoc.{a,r,u}.yaml` / `snpAssoc.b.R`| Genetic-model OR/β tables, interaction table        |
| LD and Haplotype Analysis | `snpLDHaplo.{a,r,u}.yaml` / `snpLDHaplo.b.R` | LD tables, heatmap, haplotype EM, GLM       |
| **Shared library**        | `snp_helpers.R`                          | All pure functions; no `self` / `jmvcore` calls     |

---

## File layout inside the module

```
SNPstats/
├── jamovi/
│   ├── 0000.yaml               # module descriptor (updated version + 3 analyses)
│   ├── 00refs.yaml             # references (unchanged)
│   ├── snpDesc.a.yaml          # options
│   ├── snpDesc.r.yaml          # results schema
│   ├── snpDesc.u.yaml          # UI layout
│   ├── snpAssoc.a.yaml
│   ├── snpAssoc.r.yaml
│   ├── snpAssoc.u.yaml
│   ├── snpLDHaplo.a.yaml
│   ├── snpLDHaplo.r.yaml
│   └── snpLDHaplo.u.yaml
└── R/
    ├── snp_helpers.R           # ← ALL shared pure functions
    ├── snpDesc.b.R             # ← R6 class for Descriptive
    ├── snpAssoc.b.R            # ← R6 class for Association
    └── snpLDHaplo.b.R          # ← R6 class for LD/Haplotype
```

jamovi auto-generates the `*Base` classes from the `*.a.yaml` and `*.r.yaml` files.
Each `*.b.R` inherits from its generated base (`snpDescBase`, `snpAssocBase`,
`snpLDHaploBase`). The `source("snp_helpers.R")` call at the top of each `*.b.R`
file makes all shared functions available; jamovi resolves that path relative to `R/`.

---

## The shared library: `snp_helpers.R`

Everything in this file is a **pure R function** — no `self`, no `jmvcore`, no
side-effects. This makes it testable outside jamovi with plain `Rscript`.

### Function groups

**Genotype parsing & validation**
- `split_alleles(g)` — split "A/B" string
- `check_biallelic(vals)` — returns `list($ok, $reason, $alleles)`
- `detect_snp_sep(x)` — separator detection (`/`, `|`, `>`, or `""`)
- `snp_biallelic_check(x)` — full check with reason, used for error messages
- `is_snp_column(x)` — convenience TRUE/FALSE
- `get_snp_level_order(x)` — extract user-ordered factor levels from jamovi column
- `parse_genotype(x, user_levels)` — normalise + call `genetics::genotype()`
- `get_ref_genotype(geno, user_levels)` — reference genotype selection
- `reorder_geno(gf, ref, user_levels)` — sort genotype frequency table rows

**Genetic model encoding**
- `encode_model(geno_char, ref, model, user_levels)` — codominant / dominant /
  recessive / overdominant / logadditive

**Model fitting**
- `fit_model(...)` — `glm`/`lm` for association; returns per-comparison result lists
- `fit_interaction_model(...)` — full × main LRT, returns term-level result lists

**Shared validation helpers**
- `validate_snp_vars(snp_vars, data)` → `list($valid_snps, $bad_html)`
- `detect_response_type(response_raw, responseType_opt)` → `"binary"` / `"quantitative"` / `NULL`
- `prepare_response(response_raw, response_type)` → integer / numeric vector
- `prepare_covariates(data, covariate_vars)` → factor-encoded data frame / `NULL`

---

## The danger you identified: cross-cutting changes

When a change must be applied to more than one analysis, the risk is editing one
file and forgetting the others. The architecture addresses this at different levels:

### Changes that only ever touch `snp_helpers.R` (safest)
Any bug or improvement to these areas is in one file only:
- Genotype format detection or parsing
- Reference allele / reference genotype logic
- Genetic model encoding (dosage logic, factor levels)
- `fit_model` / `fit_interaction_model` internals
- Shared validation rules (`validate_snp_vars`, `detect_response_type`, etc.)

**Rule:** if you find yourself about to edit the same logic in two `*.b.R` files,
move it to `snp_helpers.R` first.

### Changes local to one analysis (independent, no risk)
- Result table column sets in `*.r.yaml`
- UI layout in `*.u.yaml`
- Options added to `*.a.yaml` and consumed only in that analysis
- Private `.fill_*` methods in one `*.b.R`

### Changes that must propagate across analyses
Some options are deliberately **duplicated** across analyses because all three
panels expose them independently (each has its own response/SNP/covariate boxes
and `responseType`/`subpop` selectors). When you change one of these, update all
three:

| Option                              | Touches                              |
|-------------------------------------|--------------------------------------|
| `responseType` list values          | all three `*.a.yaml` + `*.u.yaml`    |
| `subpop` behaviour                  | `snpDesc.b.R`, `snpLDHaplo.b.R`      |
| `ciWidth` range / default           | `snpAssoc.a.yaml`, `snpLDHaplo.a.yaml` |
| New shared validation rule          | `snp_helpers.R` only (propagates automatically) |

To manage propagation risk, keep a short **CHANGES.md** at the module root. When
you touch a shared option, record which files need the parallel update.

---

## Adding a new option to one analysis

Example: adding a Bonferroni correction toggle to SNP Association only.

1. Add the option to `snpAssoc.a.yaml`.
2. Add the UI control to `snpAssoc.u.yaml`.
3. Add the result column (if any) to `snpAssoc.r.yaml`.
4. Consume `opts$bonferroni` inside `snpAssoc.b.R` → `.fill_assoc()`.
5. No other files change.

---

## Adding a new option that all three analyses share

Example: adding a `minN` minimum-sample-size guard.

1. Add the option to all three `*.a.yaml` files.
2. Add the UI control to all three `*.u.yaml` files.
3. Add a helper `check_min_n(n, min_n)` to `snp_helpers.R`.
4. Call it from the `.run()` method of each `*.b.R`.

---

## Testing outside jamovi

Because `snp_helpers.R` has no jamovi dependency, you can unit-test it with:

```r
source("R/snp_helpers.R")

# test genotype parsing
x <- c("A/A","A/T","T/T","A/A","A/T")
geno <- parse_genotype(x)
stopifnot(!is.null(geno))

# test model encoding
ref  <- "A/A"
enc  <- encode_model(x, ref, "dominant")
stopifnot(all(enc %in% c(0L, 1L, NA_integer_)))

# test fit_model with a toy dataset
set.seed(1)
n   <- 100
resp <- sample(0:1, n, replace=TRUE)
enc2 <- sample(c(0L,1L,2L), n, replace=TRUE)
res  <- fit_model(enc2, resp, NULL, "logadditive", "binary", 95)
stopifnot(!is.null(res))
```

---

## Checklist when editing the module

- [ ] Is the change in shared parsing/model logic? → edit `snp_helpers.R` only.
- [ ] Is the change analysis-specific? → edit only the relevant `*.b.R` + YAML pair.
- [ ] Does the change affect a duplicated option (`responseType`, `subpop`, `ciWidth`)?
      → update all affected `*.a.yaml` + `*.u.yaml` and note in `CHANGES.md`.
- [ ] Did you run the unit tests against `snp_helpers.R`?
- [ ] Did you rebuild the module (`jmvtools::prepare()`) and test in jamovi?
