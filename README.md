# snpstats — jamovi plugin

A jamovi module for genetic epidemiology SNP analysis, replicating the
functionality of the [SNPStats web application](https://www.snpstats.net)

---

## Features

| Analysis | Status |
|---|---|
| Allele frequencies | ✅ Phase 1 |
| Genotype frequencies | ✅ Phase 1 |
| Hardy-Weinberg equilibrium (exact test) | ✅ Phase 1 |
| SNP association (5 genetic models) | ✅ Phase 1 |
| Covariate adjustment | ✅ Phase 1 |
| Covariate descriptives | ✅ Phase 1 |
| Stratification by response | ✅ Phase 1 |
| Linkage disequilibrium (D, D′, r) | ✅ Phase 1 |
| Haplotype frequency estimation | ✅ Phase 1 |
| Haplotype association | ✅ Phase 1 |
| SNP × covariate interaction | 🔲 Phase 2 |
| Haplotype × covariate interaction | 🔲 Phase 2 |
| Reference category customisation per SNP | 🔲 Phase 2 |

---

## Data format

SNP columns must contain diploid genotypes in **slash-separated format**:
`C/C`, `C/T`, `T/T`. Missing values should be `NA` or empty.

Example:

```
id,snp1,snp2,status,sex,age
1,C/C,G/G,Case,Male,55
2,C/T,G/A,Control,Female,48
3,T/T,NA,Case,Male,70
```

---

## Installation (development)

```r
# Install dependencies
install.packages(c("jmvcore", "genetics", "haplo.stats", "R6"))

# Install module from source
install.packages("path/to/snpstats", repos = NULL, type = "source")

# Or via jamovi: Modules > Install from file > select the .jmo bundle
```

---

## Usage from R

```r
library(snpstats)

result <- snpAnalysis(
  data       = mydata,
  response   = "status",
  snps       = c("snp1", "snp2", "snp3"),
  covariates = c("age", "sex"),
  hweTest    = TRUE,
  snpAssoc   = TRUE,
  ldAnalysis = TRUE,
  haploFreq  = TRUE,
  haploAssoc = TRUE
)
```

---

## Five genetic models

| Model | Encoding | Tests |
|---|---|---|
| Codominant | Factor (3 levels) | Het vs ref; Alt/Alt vs ref |
| Dominant | Binary (0/1): het + alt/alt vs ref/ref | 1 df |
| Recessive | Binary (0/1): alt/alt vs ref/ref + het | 1 df |
| Overdominant | Binary (0/1): het vs both homozygotes | 1 df |
| Log-additive | Numeric 0/1/2 (alt allele dosage) | 1 df, per-allele OR |

For **binary** response: logistic regression, results as OR (95% CI, p-value).
For **quantitative** response: linear regression, results as β (95% CI, p-value).

---

## R package dependencies

- `genetics` — genotype objects, HWE.exact, LD
- `haplo.stats` — setupGeno, haplo.em, haplo.glm
- `jmvcore` — jamovi framework
- `R6` — class system

---

## File structure

```
snpstats/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── R/
│   ├── snpanalysis.b.R    # Analysis backend
│   └── snpanalysis.h.R    # Public API header
├── jamovi/
│   ├── 0000.yaml          # Module metadata
│   ├── snpanalysis.a.yaml # UI / options definition
│   └── snpanalysis.r.yaml # Results definition
└── data-raw/
    └── make_example.R     # Script to generate example dataset
```

---

## Reference

Solé X, Guinó E, Valls J, Iniesta R, Moreno V. *SNPStats: a web tool for the
analysis of association studies.* Bioinformatics. 2006;22(15):1928–9.
