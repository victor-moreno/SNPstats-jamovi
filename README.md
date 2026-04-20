# SNPstats — jamovi plugin

A jamovi module for genetic epidemiology SNP analysis, replicating the functionality of the [SNPStats web application (https://www.snpstats.net)](https://www.snpstats.net).

## Overview
The **SNPstats** module provides an interface for conducting single-SNP and multi-SNP (haplotype) association studies. It handles the complexities of genetic data, including automated format detection, HWE testing, linkage disequilimium calculation and the estimation of haplotype phases via the EM algorithm.

## Features

* **Descriptives:** Allele and genotype frequencies with subpopulation stratification.
* **Quality Control:** Hardy-Weinberg equilibrium (exact test) per SNP.
* **Association:** Analysis of SNP-response associations under multiple genetic models.
* **Response type:** Binary (logistic regression), Quantitative (linear regression).
* **Genetic models:** Codominant, Dominant, Recessive, Overdominant, Log-additive.
* **Covariates:** Adjustment for continuous/categorical variables and covariate descriptive summaries.
* **Multi-SNP analysis:** Linkage disequilibrium (D, D′, r²) statistics, matrices, and heatmaps.
* **Haplotypes:** Frequency estimation (EM algorithm) and association testing with phase uncertainty propagation, including haplotype x covariate interactions.
* **Interaction testing:** SNP × covariate and Haplotype × covariate interaction testing.


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

Solé X, Guinó E, Valls J, Iniesta R, Moreno V. *SNPStats: a web tool for the
analysis of association studies.* Bioinformatics. 2006;22(15):1928–9.
