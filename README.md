# SNPstats — jamovi plugin

A jamovi module for genetic epidemiology SNP analysis, replicating the functionality of the [SNPStats web application (https://www.snpstats.net)](https://www.snpstats.net).

## Overview

The **SNPstats** module provides an interface for conducting single-SNP and multi-SNP (haplotype) association studies. It handles the complexities of genetic data, including automated format detection, HWE testing, linkage disequilimium calculation and the estimation of haplotype phases via the EM algorithm.

See the mini [tutorial](docs/TUTORIAL.md) for more detailed information.

## 
    Features

* **Descriptives:** Allele and genotype frequencies with subpopulation stratification.
* **Quality Control:** Hardy-Weinberg equilibrium (exact test) per SNP.
* **Association:** Analysis of SNP-response associations under multiple genetic models.
* **Response type:** Binary (logistic regression), Quantitative (linear regression).
* **Genetic models:** Codominant, Dominant, Recessive, Overdominant, Log-additive.
* **Covariates:** Adjustment for continuous/categorical variables and covariate descriptive summaries.
* **Multi-SNP analysis:** Linkage disequilibrium (D, D′, r²) statistics, matrices, and heatmaps.
* **Haplotypes:** Frequency estimation (EM algorithm) and association testing with phase uncertainty propagation, including haplotype x covariate interactions.
* **Interaction testing:** SNP × covariate and Haplotype × covariate interaction testing.

## Limitations

This module is not intended for GWAS analysis. It is aimed to the detailed analysis of a few SNPs as those involved in a gene candidate analysis.

Only biallelic SNPs are supported. Other polymorphism types cannot be analyzed.

---

## Data Format

SNP columns must use diploid notation. The module automatically detects:

* **Slash-separated:** `C/C`, `C/T`, `T/T`
* **Pipe-separated:** `C|C`, `C|T`, `T|T`
* **No separator:** `CC`, `CT`, `TT`

---

## R package dependencies

- `genetics` — genotype objects, HWE.exact, LD
- `haplo.stats` — setupGeno, haplo.em, haplo.glm
- `jmvcore` — jamovi framework
- `R6` — class system

---

## Reference

Solé X, Guinó E, Valls J, Iniesta R, Moreno V. *SNPStats: a web tool for the analysis of association studies.* Bioinformatics. 2006;22(15):1928–9.

## Acknowledments

This tool has been developed with support of the Instituto de Salud Carlos III (ISCIII), “Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health (CIBERESP), action Genrisk.

This tool was migrated from SNPstats.net using Claude AI.

## License

This code is released under MIT license

## Report issues

The code has not been extensively tested and could have bugs. Please use github Issues to report any problems and enhancement requests.
