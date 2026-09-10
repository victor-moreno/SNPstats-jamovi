# SNPstats — jamovi plugin

A jamovi module for genetic epidemiology SNP analysis, replicating and enhancing the functionality of the [SNPStats web application (https://www.snpstats.net)](https://www.snpstats.net).

## Overview

The **SNPstats** module provides an interface for conducting single-SNP and multi-SNP (haplotype) association studies. It handles the complexities of genetic data, including automated format detection, HWE testing, linkage disequilimium calculation and the estimation of haplotype phases via the EM algorithm. A submodule calculates polygenic risk scores, and a third one imports genotypes directly from PLINK and VCF files so the data never has to be reshaped by hand.

See the mini [tutorial](https://victor-moreno.github.io/SNPstats-jamovi/TUTORIAL.html) for more detailed information.

## Features

* **Import:** Genotypes read straight from PLINK (`.bed`/`.bim`/`.fam`, `.ped`/`.map`, `.tped`/`.tfam`) and VCF (`.vcf`, `.vcf.gz`) files and opened as a new jamovi dataset, for a chosen list of SNPs.
* **Descriptives:** Allele and genotype frequencies with subpopulation stratification.
* **Quality Control:** Hardy-Weinberg equilibrium (exact test) per SNP.
* **Association:** Analysis of SNP-response associations under multiple genetic models.
* **Response type:** Binary (logistic regression), Quantitative (linear regression).
* **Genetic models:** Codominant, Dominant, Recessive, Overdominant, Log-additive.
* **Covariates:** Adjustment for continuous/categorical variables and covariate descriptive summaries.
* **Multi-SNP analysis:** Linkage disequilibrium (D, D′, r²) statistics, matrices, and heatmaps.
* **Haplotypes:** Frequency estimation (EM algorithm) and association testing with phase uncertainty propagation, including haplotype x covariate interactions.
* **Interaction testing:** SNP × covariate and Haplotype × covariate interaction testing.
* **Polygenic Risk Score:** Unweighted and weighted PGS using an auxiliary file in PGS Catalog format, loaded with the file-browse button (works on jamovi desktop and cloud) or, from R, via the exported `pgs_weights()` helper.

## Limitations

This module is not intended for GWAS analysis. It is aimed to the detailed analysis of a few SNPs as those involved in a gene candidate analysis.

Only biallelic SNPs are supported. Other polymorphism types cannot be analyzed.

---

## Importing genotypes

The **Import genotypes** analysis (menu *SNPstats → Data*) reads the standard
genotype formats and writes the result out as a new jamovi dataset, laid out in
the notation the analyses below expect:

* PLINK binary `.bed` + `.bim` + `.fam`
* PLINK text `.ped` + `.map`, and transposed `.tped` + `.tfam`
* VCF `.vcf` and `.vcf.gz`

**Only the SNPs you ask for are read.** For a `.bed` the byte offset of each
requested variant is computed from the `.bim` and just those ranges are sliced
out, so the source file can be arbitrarily large — 1 000 variants out of a
1.19 GB `.bed` takes 0.09 s. The text formats have no index to seek with, so
they are scanned once and only matching lines are kept.

Covariates can be merged in by sample ID, and the usual QC filters (MAF, HWE,
call rate) are applied and reported before anything is written.

Files are chosen with a browse button and read in the browser — there is no
file-path option, which is what lets the same analysis work on jamovi desktop
and in jamovi cloud. One consequence worth knowing: the selected genotypes live
in the analysis's own options while they are loaded, so an `.omv` saved with an
import still in it contains that genotype data. Clear the file selection before
sharing such a file.

---

## Data Format

The importer above produces this layout automatically; these rules matter when
you bring a spreadsheet in by other means.

SNP columns must use diploid notation. The module automatically detects:

* **Slash-separated:** `C/C`, `C/T`, `T/T`
* **Pipe-separated:** `C|C`, `C|T`, `T|T`
* **No separator:** `CC`, `CT`, `TT`
* **Dosage:** (only for PGS). Values in `[0,2]`, possibly decimal from imputation dosage estimates. 

Missing values: `'', NA, 'NA', 'N/A', 'N|A', '0/0'`

---

## R package dependencies

- `jmvcore` — jamovi framework
- `R6` — class system
- `nnet` — multinomial regression (categorical responses)
- `haplo.stats` — `setupGeno`, `haplo.em`, `haplo.glm` (haplotype estimation
  and association)
- `ggplot2` — LD heatmap and all PGS plots
- `base64enc` — decoding the embedded PGS weights file and the genotype
  payloads the import panel sends from the browser

Genotype parsing, the Hardy-Weinberg exact test and pairwise LD are implemented
in the module itself (`R/snp_genetics.R`); the `genetics` package was dropped in
v1.0.0 as it is marked obsolete upstream. The PLINK and VCF readers
(`R/plink_read.R`, `R/text_formats.R`) are likewise implemented from the format
specifications rather than by calling PLINK or htslib, so importing needs no
external binary.

---

## Reference

Solé X, Guinó E, Valls J, Iniesta R, Moreno V. *SNPStats: a web tool for the analysis of association studies.* Bioinformatics. 2006;22(15):1928–9.

## Acknowledments

This tool has been developed with support of the Instituto de Salud Carlos III (ISCIII), “Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health (CIBERESP), action Genrisk.

This tool was migrated from SNPstats.net using Claude AI.

Earlier versions used the [`genetics`](https://CRAN.R-project.org/package=genetics)
R package for genotype handling, the Hardy-Weinberg exact test and linkage
disequilibrium. Those computations are now implemented directly from the
original publications (Wigginton *et al.* 2005; Excoffier & Slatkin 1995;
Lewontin 1964; Hill & Robertson 1968), with thanks to that package's authors for
the reference implementation.

## License

This project is licensed under the **GNU General Public License v3** — see
[LICENSE.md](LICENSE.md) for the summary and [COPYING](COPYING) for the full
text.


## Report issues

The code has not been extensively tested and could have bugs. Please use github Issues to report any problems and enhancement requests.

Often problems are related to data, mainly when sample size is small and zeros appear in cells. Then, models may fail.

© [Catalan Institute of Oncology](http://iconcologia.net/)

