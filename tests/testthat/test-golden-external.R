# External-oracle golden test: values captured from the reference SNPStats web
# tool (https://www.snpstats.net) run on data/sample3.tsv (response = group,
# adjusted by sex + age + family). Unlike the oracle tests — which recompute
# with base R — these pin the module against the *independent reference tool* it
# claims compatibility with, so a future change that drifts away from snpstats.net
# is caught. See docs: association p-values are the model-level LRT (binary),
# matching the web tool.
#
# If one fails after an intentional change, re-verify against the web tool before
# updating the value.

.ext_data <- read.delim(
  Filter(file.exists, c("data/sample3.tsv", file.path("..", "..", "data", "sample3.tsv")))[[1]],
  header = TRUE, stringsAsFactors = TRUE)

# Degenerate genotype cells (snp7/snp8 have a genotype absent in one group)
# legitimately warn about separation; the module reports the LRT correctly and
# flags it, so suppress the console noise here.
.ext_assoc <- function(snp) {
  res <- suppressWarnings(run_snp(
    data = .ext_data, snps = snp, response = "group",
    covariates = c("sex", "age", "family"), snpAssoc = TRUE,
    showAIC = TRUE, modelCodominant = TRUE, modelDominant = TRUE,
    modelRecessive = TRUE, modelOverdominant = TRUE, modelLogAdditive = TRUE))
  as_df(res$assocGroup$assocSnpResults$get(key = snp)$assocTable)
}

# Single row for a model's first (labelled) row; which() drops the NA-model rows.
.model_row <- function(tbl, model_label) tbl[which(as.character(tbl$model) == model_label)[1], ]
.model_p   <- function(tbl, model_label) num(.model_row(tbl, model_label)$pval)

# ── Descriptives vs web tool ──────────────────────────────────────────────────
test_that("EXTERNAL snp1 descriptives match snpstats.net", {
  st <- as_df(run_snp(data = .ext_data, snps = c("snp1", "snp7"),
                      snpSummary = TRUE)$descGroup$snpSummaryTablesGroup$snpSummaryTable)
  r1 <- st[st$snp == "snp1", ]
  expect_equal(as.integer(r1$n), 536L)
  expect_equal(as.integer(r1$missing), 30L)
  expect_close(num(r1$maf), 0.16, tol = 0.005)                 # web: G 0.16
  expect_equal(as.character(r1$genoCounts), "382 / 136 / 18")  # web: 382/136/18
  expect_close(num(r1$hwePval), 0.20, tol = 0.01)              # web: 0.2
})

# ── Association p-values (LRT) vs web tool, strong signal (snp7) ───────────────
test_that("EXTERNAL snp7 association reproduces snpstats.net (OR, CI, LRT p)", {
  t <- .ext_assoc("snp7")
  # Web log-additive: OR 0.37 (0.18-0.78), p 0.0046. CI bounds differ by up to
  # ~0.04 because the module uses profile-likelihood CIs and the web tool Wald.
  la <- .model_row(t, "Additive")
  expect_close(num(la$effect), 0.37, tol = 0.01)
  expect_close(num(la$ciLow),  0.18, tol = 0.04)
  expect_close(num(la$ciHigh), 0.78, tol = 0.04)
  expect_close(num(la$pval),   0.0046, tol = 0.001)
  # Web model-level LRT p per model
  expect_close(.model_p(t, "Codominant"),   0.012, tol = 0.001)
  expect_close(.model_p(t, "Dominant"),     0.008, tol = 0.001)
  expect_close(.model_p(t, "Recessive"),    0.056, tol = 0.002)  # separation: LRT (Wald was wrong)
  expect_close(.model_p(t, "Overdominant"), 0.025, tol = 0.001)
  # Web codominant AIC/BIC
  cod <- .model_row(t, "Codominant")
  expect_close(num(cod$AIC), 743.6, tol = 0.1)
  expect_close(num(cod$BIC), 773.8, tol = 0.1)
})

# ── Association p-values (LRT) vs web tool, second signal (snp8) ───────────────
test_that("EXTERNAL snp8 association reproduces snpstats.net", {
  t <- .ext_assoc("snp8")
  la <- .model_row(t, "Additive")
  expect_close(num(la$effect), 0.41, tol = 0.01)               # web 0.41 (0.19-0.90)
  expect_close(num(la$pval),   0.020, tol = 0.001)             # web 0.02
  expect_close(.model_p(t, "Codominant"),   0.058, tol = 0.002)
  expect_close(.model_p(t, "Dominant"),     0.024, tol = 0.001)
  expect_close(.model_p(t, "Overdominant"), 0.033, tol = 0.001)
})

# ── Null-signal SNP (snp5) matches the web tool's ~1 ORs and large p ───────────
test_that("EXTERNAL snp5 null association reproduces snpstats.net", {
  t <- .ext_assoc("snp5")
  la <- .model_row(t, "Additive")
  expect_close(num(la$effect), 0.98, tol = 0.01)               # web 0.98 (0.73-1.32)
  expect_close(.model_p(t, "Codominant"), 0.99, tol = 0.02)    # web 0.99
  expect_close(.model_p(t, "Additive"), 0.90, tol = 0.02)  # web 0.9
})
