# ![snpstats logo](snpstats_logo.svg)
# SNP Analysis Tutorial

## Data format

Each row in the dataset represents one individual. SNP genotypes must be stored as character columns using diploid notation — two alleles separated by a delimiter, or concatenated as two characters. The module accepts the following formats:

Any allele names are accepted (single nucleotides, insertion/deletion codes, etc.), provided that exactly two distinct alleles appear across the column. Columns with three or more alleles, or inconsistent genotype combinations (e.g. A/A, A/C, B/C) are flagged and excluded with an informative message.

The response variable should be a binary (case/control coded 0/1 or as two-level factor) or continuous column. Covariates can be numeric or categorical.

| Format | Separator | Example values |
| --- | --- | --- |
| A/B | slash: / | T/T, T/C, C/C |
| A\|B | pipe: \| | T\|T, T\|C, C\|C |
| A>B | greater-than > | T>T, T>C, C&gt;C |
| AB | none (2 chars) | TT, TC, CC |

## Covariates

Enable Covariate summary (requires at least one covariate) to obtain a descriptive table of all covariate variables.

When Subpopulation is active and the response is binary, the table is stratified by response group. Separate columns appear for each group, and p-values are added: χ² test for categorical variables, independent-samples t-test for continuous variables.


## Single-SNP analysis

A polymorphism occurs when different individuals carry different genetic variants at the same genomic location (locus). When only one nucleotide differs, the variant is called a SNP (Single Nucleotide Polymorphism). For a biallelic locus with alleles T and C, three genotypes are possible: T/T, T/C, and C/C.

Enable Genotype frequencies and Allele frequencies under SNP descriptive to obtain per-SNP tables. Counts and proportions are computed over typed individuals only (missing genotypes are excluded from denominators but reported separately).

Allele frequencies are computed by counting each allele separately across both chromosomes of every typed individual:

Enable Subpopulation (requires a binary response) to show frequencies stratified by case/control status.

| Genotype | All (n) | All (%) | Controls (n) | Cases (n) |
| --- | --- | --- | --- | --- |
| T/T | 255 | 40.2% | 123 | 132 |
| T/C | 264 | 41.6% | 126 | 138 |
| C/C | 116 | 18.3% | 52 | 64 |
| Missing | 71 | — | 28 | 43 |

$$\text{Allele freq}(T) = \frac{2 \times N_{TT} + N_{TC}}{2 \times N_{\text{typed}}}$$

### Descriptive statistics

Before an association analysis, it is good practice to verify that each SNP is in Hardy–Weinberg equilibrium (HWE). Under HWE, allele frequencies within a large, randomly mating population remain constant across generations. A significant deviation may indicate genotyping error, population stratification, or selection pressure.

Enable HWE test to obtain an exact test p-value for each SNP, overall and (when Subpopulation is active) per stratum.

| Group | NTT | NTC | NCC | P-value (exact) |
| --- | --- | --- | --- | --- |
| All | 255 | 264 | 116 | 0.0015 |
| Controls | 123 | 126 | 52 | 0.051 |
| Cases | 132 | 138 | 64 | 0.012 |

### Hardy–Weinberg equilibrium

Enable SNP association under SNP association to fit regression models for each SNP against the response. The module fits logistic regression for binary outcomes and linear regression for continuous outcomes.

For each selected genetic model, the table reports:

A note at the bottom of the table shows the likelihood ratio test p-value for the codominant model, which serves as the overall test for association.

When covariates are present, a note lists the adjustment variables, and a separate note indicates the number of observations excluded due to missing covariate values.

### Association with outcome

Five inheritance models are available. Let A be the reference allele (most frequent homozygote) and B the variant allele:

| Model | Comparison | OR (95% CI) | P-value | AIC |
| --- | --- | --- | --- | --- |
| Codominant | T/C vs T/T | 1.02 (0.72–1.44) | 0.82 | 884.2 |
|  | C/C vs T/T | 1.15 (0.74–1.78) |  |  |
| Dominant | T/C+C/C vs T/T | 1.06 (0.77–1.45) | 0.73 | 882.5 |
| Recessive | C/C vs T/T+T/C | 1.14 (0.76–1.70) | 0.54 | 882.2 |
| Overdominant | T/C vs T/T+C/C | 0.98 (0.71–1.34) | 0.89 | 882.6 |
| Log-additive | per B allele | 1.06 (0.86–1.31) | 0.58 | 882.3 |

### Genetic models

Enable SNP × covariate interaction (requires at least one covariate and SNP association to be enabled) to test whether the SNP effect differs across levels of the first covariate. The model fitted is:

The table shows regression coefficients (or ORs) for SNP main effects and SNP × covariate interaction terms. The P (interaction) column shows the likelihood ratio test p-value comparing the model with interaction to the model without; this value appears once per genetic model on the first interaction term row.

$$logit(p) = β₀ + β_{SNP}·G + β_{Z}·Z + β_{int}·(G × Z)$$

### SNP × covariate interaction


## Multi-SNP analysis

When multiple SNPs are placed in the SNPs box, additional analyses become available under Linkage disequilibrium and Haplotype analysis.

Linkage disequilibrium (LD) is the non-random association between alleles at different loci on the same chromosome. It arises because chromosomes are transmitted in blocks, with recombination between neighbouring loci being rare. SNPs that are close to a causal variant will tend to be in LD with it, making LD analysis central to fine-mapping.

Enable Pairwise LD to obtain a table of pairwise LD statistics for all SNP pairs. Enable LD matrix to display them in matrix form; the LD metric dropdown selects which statistic to use for the matrix. Enable LD heatmap for a colour-coded plot.

Three statistics are reported:

A high D′ with a lower r² is common when allele frequencies differ markedly between loci: the two SNPs always co-occur on the same ancestral haplotype (D′ ≈ 1) but because one allele is rare, the statistical correlation r² is attenuated.

| SNP pair | D | D′ | r² | P-value |
| --- | --- | --- | --- | --- |
| SNP1 – SNP2 | 0.216 | 0.981 | 0.921 | &lt;0.001 |

### Linkage disequilibrium

A haplotype is the combination of alleles across multiple loci on a single chromosome. Because standard genotyping does not determine which alleles are on the same chromosome (phase is unknown for doubly heterozygous individuals), haplotypes must be inferred statistically.

Enable Haplotype frequencies to run the EM algorithm, which iterates between:

The Minimum haplotype frequency option (default 0.05) pools all haplotypes below the threshold into a single Rare category.

| # | SNP1 | SNP2 | Frequency | Controls | Cases |
| --- | --- | --- | --- | --- | --- |
| 1 | C | T | 0.606 | 0.617 | 0.597 |
| 2 | G | C | 0.356 | 0.361 | 0.350 |
| 3 | C | C | 0.034 | 0.020 | 0.047 |
| Rare | * | * | 0.004 | — | — |

### Haplotype frequencies

Enable Haplotype association (requires a response variable) to test each haplotype for association with the outcome. The method uses haplo.glm from the haplo.stats package, which fits a regression model weighted by the posterior probability of each haplotype assignment for each individual — correctly propagating phase uncertainty into the estimates.

The most frequent haplotype serves as the reference category. The table reports OR (binary) or β (continuous) with confidence interval and p-value for each common haplotype, plus the pooled Rare term.

| Haplotype | Frequency | OR (95% CI) | P-value |
| --- | --- | --- | --- |
| C-T (Ref) | 0.606 | 1.00 | — |
| G-C | 0.356 | 0.99 (0.80–1.23) | 0.94 |
| C-C | 0.034 | 2.41 (1.22–4.75) | 0.011 |
| Rare | 0.004 | 3.78 (0.42–34.2) | 0.24 |

### Haplotype association and interaction with covariates

Enable Haplotype × covariate interaction (requires a response and at least one covariate) to test whether haplotype effects differ across levels of the first covariate. The interaction model is:
$$logit(p) = β₀ + \sum{ β_{h}H_{h}} + β_{Z}Z + \sum{ β_{h×Z}(H_{h} × Z)}$$

The table lists the main haplotype effects and their interaction terms. A table note reports the overall likelihood ratio test for interaction, computed from the deviance difference between the interaction and main-effects models.


## References

Sole X, Guino E, Valls J, Iniesta R, Moreno V (2006). SNPStats: a web tool for the analysis of association studies. Bioinformatics 22(15):1928–1929.

Schaid DJ, Rowland CM, Tines DE, Jacobson RM, Poland GA (2002). Score tests for association between traits and haplotypes when linkage phase is ambiguous. Am J Hum Genet 70:425–434.

Akaike H (1974). A new look at the statistical model identification. IEEE Trans Autom Control 19(6):716–723.

