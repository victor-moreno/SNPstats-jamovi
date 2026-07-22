# Degenerate-input robustness tests. These assert the module degrades
# gracefully (no uncaught R error, sensible structure) on inputs a public user
# can hit: monomorphic SNPs, all-missing columns, single-SNP LD, triallelic
# markers, degenerate covariates/response. Several were crashes before the
# monomorphic / empty-genotype-vector fixes.

.edge_data <- local({
  set.seed(11); n <- 140
  g <- function(p) vapply(seq_len(n), function(i)
    paste(sample(c("A", "G"), 2, TRUE, c(1 - p, p)), collapse = "/"), "")
  y <- factor(sample(c("Case", "Control"), n, TRUE))
  data.frame(
    s1        = g(0.30),
    s2        = g(0.20),
    mono      = rep("A/A", n),
    allhet    = rep("A/G", n),
    tri       = sample(c("A/C", "A/G", "C/G"), n, TRUE),
    monoInOne = ifelse(y == "Case", "A/A", g(0.25)),
    dose      = sample(0:2, n, TRUE),
    y         = y,
    age       = round(rnorm(n, 50, 10), 1),
    allNAcov  = rep(NA_real_, n),
    constResp = rep(5, n),
    stringsAsFactors = FALSE)
})

# Degenerate fits legitimately warn (separation / non-convergence) — the module
# captures these as table notes; suppress them here so they don't drown the log.
expect_runs <- function(expr) testthat::expect_error(suppressWarnings(expr), NA)

# ══════════════════════════════════════════════════════════════════════════════
# snpStats — degenerate SNPs
# ══════════════════════════════════════════════════════════════════════════════

test_that("monomorphic SNP: descriptives run and both SNPs appear", {
  result <- run_snp(data = .edge_data, snps = c("s1", "mono"),
                    snpSummary = TRUE, allFreq = TRUE, genoFreq = TRUE, hweTest = TRUE)
  tbl <- as_df(result$descGroup$snpSummaryTablesGroup$snpSummaryTable)
  expect_setequal(tbl$snp, c("s1", "mono"))
})

test_that("monomorphic SNP: association runs (empty assoc, no crash)", {
  expect_runs(run_snp(data = .edge_data, snps = "mono", response = "y",
                      snpAssoc = TRUE, modelCodominant = TRUE, modelLogAdditive = TRUE))
})

test_that("monomorphic SNP: stratified descriptives and interaction run", {
  expect_runs(run_snp(data = .edge_data, snps = c("s1", "mono"), response = "y",
                      subpop = TRUE, snpSummary = TRUE, allFreq = TRUE,
                      genoFreq = TRUE, hweTest = TRUE))
  expect_runs(run_snp(data = .edge_data, snps = "mono", response = "y",
                      covariates = "age", snpInteraction = TRUE, modelCodominant = TRUE))
})

test_that("SNP monomorphic within one response stratum runs all models", {
  expect_runs(run_snp(data = .edge_data, snps = "monoInOne", response = "y",
                      subpop = TRUE, snpAssoc = TRUE,
                      modelCodominant = TRUE, modelDominant = TRUE, modelRecessive = TRUE))
})

test_that("all-heterozygote SNP: descriptives and association run", {
  expect_runs(run_snp(data = .edge_data, snps = "allhet",
                      snpSummary = TRUE, allFreq = TRUE, genoFreq = TRUE, hweTest = TRUE))
  expect_runs(run_snp(data = .edge_data, snps = "allhet", response = "y", snpAssoc = TRUE,
                      modelCodominant = TRUE, modelRecessive = TRUE, modelOverdominant = TRUE))
})

test_that("triallelic SNP is skipped: valid SNP still computed, no crash", {
  result <- suppressWarnings(
    run_snp(data = .edge_data, snps = c("s1", "tri"), snpSummary = TRUE))
  tbl <- as_df(result$descGroup$snpSummaryTablesGroup$snpSummaryTable)
  expect_true("s1" %in% as.character(tbl$snp))
  expect_false(is.na(as.integer(tbl$n[tbl$snp == "s1"])))
})

test_that("all-missing SNP is skipped: valid SNP still computed, no crash", {
  d <- .edge_data; d$allmiss <- NA_character_
  result <- suppressWarnings(
    run_snp(data = d, snps = c("s1", "allmiss"), snpSummary = TRUE))
  tbl <- as_df(result$descGroup$snpSummaryTablesGroup$snpSummaryTable)
  expect_true("s1" %in% as.character(tbl$snp))
  expect_false(is.na(as.integer(tbl$n[tbl$snp == "s1"])))
})

test_that("numeric 0/1/2 SNP column does not crash snpStats", {
  expect_runs(run_snp(data = .edge_data, snps = "dose", snpSummary = TRUE))
})

test_that("duplicate SNP selection does not crash", {
  expect_runs(run_snp(data = .edge_data, snps = c("s1", "s1"), snpSummary = TRUE))
})

# ══════════════════════════════════════════════════════════════════════════════
# snpStats — LD / haplotype with too few or degenerate SNPs
# ══════════════════════════════════════════════════════════════════════════════

test_that("single SNP with LD options on does not crash", {
  expect_runs(run_snp(data = .edge_data, snps = "s1",
                      ldAnalysis = TRUE, ldMatrix = TRUE, ldPlot = TRUE))
})

test_that("monomorphic + valid SNP pair: LD and haplotype run", {
  expect_runs(run_snp(data = .edge_data, snps = c("s1", "mono"), response = "y",
                      ldAnalysis = TRUE, ldMatrix = TRUE,
                      haploFreq = TRUE, haploAssoc = TRUE))
})

# ══════════════════════════════════════════════════════════════════════════════
# snpStats — degenerate covariate / response
# ══════════════════════════════════════════════════════════════════════════════

test_that("all-NA covariate and constant response do not crash association", {
  expect_runs(run_snp(data = .edge_data, snps = "s1", response = "y",
                      covariates = "allNAcov", snpAssoc = TRUE))
  expect_runs(run_snp(data = .edge_data, snps = "s1", response = "constResp",
                      snpAssoc = TRUE))
})

# ══════════════════════════════════════════════════════════════════════════════
# snpPGS — degenerate SNPs / response
# ══════════════════════════════════════════════════════════════════════════════

test_that("PGS handles all-missing and monomorphic SNPs (excluded, no crash)", {
  d <- .edge_data; d$allmiss <- NA_character_
  expect_runs(run_pgs(data = d, snpCols = c("s1", "allmiss", "mono"),
                      showCoverage = TRUE, showAssoc = TRUE, responseCol = "y"))
})

test_that("PGS handles out-of-range numeric dosage and triallelic markers", {
  d <- .edge_data; d$oor <- sample(c(0, 1, 2, 3), nrow(d), TRUE)
  expect_runs(run_pgs(data = d, snpCols = c("s1", "oor", "tri"), showCoverage = TRUE))
})

test_that("PGS constant response and response==covariate do not crash", {
  expect_runs(run_pgs(data = .edge_data, snpCols = c("s1", "s2"),
                      responseCol = "constResp", showAssoc = TRUE))
  expect_runs(run_pgs(data = .edge_data, snpCols = c("s1", "s2"),
                      responseCol = "y", covCols = "y", showInteraction = TRUE))
})
