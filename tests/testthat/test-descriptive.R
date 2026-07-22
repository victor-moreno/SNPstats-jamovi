# Tab 1: Descriptive results — verified against independent allele counting and
# genetics::HWE.exact.

suppressMessages(library(genetics))

# Oracle: count alleles / genotypes / HWE directly from the raw "A/B" column.
desc_oracle <- function(col) {
  s <- as.character(col); s[grepl("0", s)] <- NA
  typed <- s[!is.na(s)]
  go    <- genetics::genotype(typed, sep = "/")
  af    <- summary(go)$allele.freq
  props <- af[rownames(af) != "NA", "Proportion"]
  list(
    n       = length(typed),
    missing = sum(is.na(s)),
    maf     = unname(min(props)),
    hwe     = genetics::HWE.exact(go)$p.value)
}

# ══════════════════════════════════════════════════════════════════════════════
# Response level order — a character response must follow data order of
# appearance, not R's alphabetical as.factor() default (homogeneous with snpPGS).
# ══════════════════════════════════════════════════════════════════════════════

test_that("snp_prepare: character response honors data order of appearance", {
  sp <- getFromNamespace("snp_prepare", "SNPstats")
  # appearance order Control, Case (alphabetical would be Case, Control)
  d  <- data.frame(y = c("Control", "Control", "Case", "Case"),
                   g = c("A/A", "A/G", "G/G", "A/G"), stringsAsFactors = FALSE)
  p  <- sp(d, snps = "g", response = "y")
  expect_identical(levels(p$response_raw), c("Control", "Case"))
  expect_identical(p$response_enc, c(0L, 0L, 1L, 1L))   # Control (first seen) = 0
})

# ══════════════════════════════════════════════════════════════════════════════
# snpSummaryTable — N, missing, MAF, genotype counts, HWE p-value
# ══════════════════════════════════════════════════════════════════════════════

test_that("snpSummary: values match independent allele counting + HWE.exact", {
  result <- run_snp(data = .test_data, snps = .snps2, snpSummary = TRUE)
  tbl <- as_df(result$descGroup$snpSummaryTablesGroup$snpSummaryTable)

  expect_equal(nrow(tbl), 2L)
  expect_setequal(tbl$snp, .snps2)

  for (snp in .snps2) {
    o   <- desc_oracle(.test_data[[snp]])
    row <- tbl[tbl$snp == snp, ]
    expect_equal(as.integer(row$n),       o$n,       label = paste(snp, "N"))
    expect_equal(as.integer(row$missing), o$missing, label = paste(snp, "missing"))
    expect_close(num(row$maf),     o$maf, tol = 0.0005, label = paste(snp, "MAF"))
    expect_close(num(row$hwePval), o$hwe, tol = 0.01,   label = paste(snp, "HWE p"))

    # genotype counts AA / AB / BB must sum to N
    counts <- as.integer(strsplit(as.character(row$genoCounts), "\\s*/\\s*")[[1]])
    expect_equal(length(counts), 3L)
    expect_equal(sum(counts), o$n, label = paste(snp, "geno counts sum"))
    # minor-allele count from genotypes reproduces the MAF
    maf_from_counts <- (counts[2] + 2 * counts[3]) / (2 * o$n)
    expect_close(maf_from_counts, o$maf, tol = 0.0005, label = paste(snp, "MAF from counts"))
  }
})

test_that("snpSummary: stratified by response adds one row per group per SNP", {
  result <- run_snp(data = .test_data, snps = .snps2, response = .resp,
                    snpSummary = TRUE, subpop = TRUE)
  tbl <- as_df(result$descGroup$snpSummaryTablesGroup$snpSummaryTable)

  expect_true("group" %in% names(tbl))
  # overall + 2 phenotype groups, per SNP
  expect_equal(nrow(tbl), 2L * 3L)
  expect_setequal(unique(tbl$group[tbl$group != ""]),
                  c("All", levels(.test_data$phenotype)))
})

# ══════════════════════════════════════════════════════════════════════════════
# allFreqTable / genoFreqTable / hweTable (per-SNP arrays)
# ══════════════════════════════════════════════════════════════════════════════

test_that("allFreq: allele counts sum to 2N and match minor-allele frequency", {
  result <- run_snp(data = .test_data, snps = .snps2, allFreq = TRUE)

  for (snp in .snps2) {
    o   <- desc_oracle(.test_data[[snp]])
    tab <- as_df(result$descGroup$descSnpResults$get(key = snp)$allFreqTable)
    expect_equal(nrow(tab), 2L, label = paste(snp, "two alleles"))
    counts <- as.integer(sub("\\s*\\(.*$", "", tab$stat))
    expect_equal(sum(counts), 2L * o$n, label = paste(snp, "alleles sum to 2N"))
    expect_close(min(counts) / (2 * o$n), o$maf, tol = 0.0005,
                 label = paste(snp, "MAF from allele table"))
  }
})

test_that("genoFreq: genotype counts sum to N", {
  result <- run_snp(data = .test_data, snps = .snps2, genoFreq = TRUE)

  for (snp in .snps2) {
    o   <- desc_oracle(.test_data[[snp]])
    tab <- as_df(result$descGroup$descSnpResults$get(key = snp)$genoFreqTable)
    counts <- as.integer(sub("\\s*\\(.*$", "", tab$stat))
    expect_equal(sum(counts), o$n, label = paste(snp, "geno counts sum to N"))
  }
})

test_that("hweTest: per-SNP HWE p-value matches HWE.exact", {
  result <- run_snp(data = .test_data, snps = .snps2, hweTest = TRUE)

  for (snp in .snps2) {
    o   <- desc_oracle(.test_data[[snp]])
    tab <- as_df(result$descGroup$descSnpResults$get(key = snp)$hweTable)
    overall <- tab[tab$group %in% c("Overall", "All", "") | nrow(tab) == 1, ][1, ]
    expect_close(num(overall$pval), o$hwe, tol = 0.01, label = paste(snp, "HWE p"))
    expect_equal(overall$n11 + overall$n12 + overall$n22, o$n,
                 label = paste(snp, "HWE table counts sum to N"))
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# covDescTable
# ══════════════════════════════════════════════════════════════════════════════

test_that("covDesc: non-stratified lists every covariate", {
  result <- run_snp(data = .test_data, snps = .snps2, response = .resp,
                    covariates = .covars, covDesc = TRUE)
  tbl <- as_df(result$descGroup$covDescGroup$covDescTable)

  expect_gt(nrow(tbl), 0L)
  expect_true(all(c("variable", "level", "stat_overall") %in% names(tbl)))
  expect_true(all(.covars %in% tbl$variable))
})

test_that("covDesc: stratified group difference p-value matches a t-test for age", {
  result <- run_snp(data = .test_data, snps = .snps2, response = .resp,
                    covariates = .covars, covDesc = TRUE, subpop = TRUE)
  tbl <- as_df(result$descGroup$covDescGroup$covDescTable)

  expect_true("pval" %in% names(tbl))
  age_p <- num(tbl$pval[tbl$variable == "age"])
  age_p <- age_p[!is.na(age_p)][1]
  ref_p <- t.test(age ~ phenotype, data = .test_data, var.equal = TRUE)$p.value
  expect_close(age_p, ref_p, tol = 0.01)
})
