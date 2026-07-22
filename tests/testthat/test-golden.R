# Golden / regression values — captured from a known-good run and verified
# against glm/lm/multinom/genetics/haplo.stats (see the other test files).
#
# Purpose: these pin the *exact* numbers the module currently produces, so a
# future change in R or in a dependency (e.g. a different confint() or EM
# implementation) is caught even when the module's own logic is unchanged.
# If one fails after upgrading R/packages, re-verify against the oracle tests,
# then update the value here.

# Helper: numeric coefficient rows of a one-model assocTable (accepts a Table or
# an already-converted data.frame).
.golden_coef <- function(x) {
  df <- if (is.data.frame(x)) x else as_df(x)
  df[!is.na(num(df$ciLow)), , drop = FALSE]
}

# ── Descriptive ──────────────────────────────────────────────────────────────
test_that("GOLDEN snpSummary", {
  tbl <- as_df(run_snp(data = .test_data, snps = .snps2,
                       snpSummary = TRUE)$descGroup$snpSummaryTablesGroup$snpSummaryTable)
  r1 <- tbl[tbl$snp == "rs12080929", ]; r2 <- tbl[tbl$snp == "rs10911251", ]
  expect_equal(as.integer(r1$n), 2827L);   expect_equal(as.integer(r1$missing), 11L)
  expect_close(num(r1$maf), 0.2706, tol = 0.0005)
  expect_close(num(r1$hwePval), 0.8863, tol = 0.001)
  expect_equal(as.character(r1$genoCounts), "1502 / 1120 / 205")
  expect_close(num(r2$maf), 0.4110, tol = 0.0005)
  expect_close(num(r2$hwePval), 0.5865, tol = 0.001)
  expect_equal(as.character(r2$genoCounts), "973 / 1384 / 470")
})

# ── Association, binary — associated SNP rs10936599, all five models ──────────
test_that("GOLDEN association rs10936599 (binary)", {
  res <- run_snp(data = .test_data, snps = "rs10936599", response = .resp,
                 snpAssoc = TRUE, modelCodominant = TRUE, modelDominant = TRUE,
                 modelRecessive = TRUE, modelOverdominant = TRUE,
                 modelLogAdditive = TRUE, showAIC = TRUE)
  tbl <- as_df(res$assocGroup$assocSnpResults$get(key = "rs10936599")$assocTable)
  co <- .golden_coef(tbl)
  expect_close(num(co$effect[1]), 1.169, tol = 0.002)   # C/T
  expect_close(num(co$ciLow[1]),  0.996, tol = 0.002)
  expect_close(num(co$ciHigh[1]), 1.372, tol = 0.002)
  expect_close(num(co$effect[2]), 1.346, tol = 0.002)   # T/T
  # dominant carrier, recessive, overdominant het, log-additive
  expect_close(num(co$effect[3]), 1.192, tol = 0.002)
  expect_close(num(co$effect[4]), 1.273, tol = 0.002)
  expect_close(num(co$effect[5]), 1.141, tol = 0.002)
  expect_close(num(co$effect[6]), 1.165, tol = 0.002)   # per-allele
  expect_close(num(co$pval[6]),   0.017, tol = 0.001)   # log-additive LRT p (was Wald 0.0174)
  # log-additive AIC/BIC
  la <- tbl[which(tbl$model == "Additive"), ]
  expect_close(num(la$AIC), 3865.49, tol = 0.05)
  expect_close(num(la$BIC), 3877.39, tol = 0.05)
})

# ── Association, quantitative (age) and categorical (bmiOMS) ──────────────────
test_that("GOLDEN association rs10936599 (quantitative + categorical)", {
  q <- .golden_coef(run_snp(data = .test_data, snps = "rs10936599", response = "age",
        responseType = "quantitative", snpAssoc = TRUE, modelLogAdditive = TRUE,
        modelCodominant = FALSE)$assocGroup$assocSnpResults$get(key = "rs10936599")$assocTable)
  expect_close(num(q$effect[1]), -1.294, tol = 0.002)
  expect_close(num(q$ciLow[1]),  -2.386, tol = 0.002)
  expect_close(num(q$pval[1]),    0.0202, tol = 0.001)

  cc <- as_df(run_snp(data = .test_data, snps = "rs10936599", response = "bmiOMS",
        responseType = "categorical", snpAssoc = TRUE, modelLogAdditive = TRUE,
        modelCodominant = FALSE)$assocGroup$assocSnpResults$get(key = "rs10936599")$assocTable)
  cc <- cc[cc$genotype == "Per allele", ]
  expect_close(num(cc$effect[1]), 0.759, tol = 0.005)   # Obese
  expect_close(num(cc$effect[2]), 1.079, tol = 0.005)   # Overweight
  expect_close(num(cc$effect[3]), 1.189, tol = 0.005)   # Underweight
})

# ── Interaction (multiplicative) ─────────────────────────────────────────────
test_that("GOLDEN multiplicative interaction rs10936599 x sex", {
  res <- run_snp(data = .test_data, snps = "rs10936599", response = .resp,
                 covariates = "sex", responseType = "binary", snpInteraction = TRUE,
                 interactionType = "multiplicative", modelLogAdditive = TRUE,
                 modelCodominant = FALSE, showInteractionTable = TRUE)
  it <- as_df(res$assocGroup$assocSnpResults$get(key = "rs10936599")$interactionTable)
  expect_close(num(it$effect[it$term == "snp"]),         1.114, tol = 0.002)
  expect_close(num(it$effect[grepl(":", it$term)]),      1.087, tol = 0.002)
  expect_close(num(it$effect[it$term == "sexMale"]),     0.533, tol = 0.002)
})

# ── LD ───────────────────────────────────────────────────────────────────────
test_that("GOLDEN LD pair rs72647484 / rs6691170", {
  tbl <- as_df(run_snp(data = .test_data, snps = .snps4,
               ldAnalysis = TRUE)$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldTable)
  row <- tbl[tbl$snp1 == "rs72647484" & tbl$snp2 == "rs6691170", ]
  expect_close(num(row$r2),     0.001, tol = 0.0015)
  expect_close(num(row$Dprime), 0.070, tol = 0.0015)
  expect_close(num(row$pval),   0.0585, tol = 0.002)
})

# ── Haplotype frequencies (deterministic) ────────────────────────────────────
test_that("GOLDEN haploFreq (min 0.05)", {
  tbl <- as_df(run_snp(data = .test_data, snps = .snps4, response = .resp,
               haploFreq = TRUE, haploFreqMin = 0.05)$ldHaploGroup$haploGroup$haploFreqTable)
  f <- setNames(num(tbl$freq), tbl$haplotype)
  expect_close(f[["T-A-T-G"]], 0.258, tol = 0.002)
  expect_close(f[["T-C-T-G"]], 0.178, tol = 0.002)
  expect_close(f[["T-A-T-T"]], 0.140, tol = 0.002)
  expect_close(f[["Rare (<0.05)"]], 0.110, tol = 0.002)
})

# ── Haplotype association ─────────────────────────────────────────────────────
# The module seeds haplo.glm's EM internally, so this is deterministic without
# any external set.seed (tight tolerance is therefore appropriate).
test_that("GOLDEN haploAssoc snps2", {
  tbl <- as_df(run_snp(data = .test_data, snps = .snps2, response = .resp,
               haploAssoc = TRUE, haploFreqMin = 0.05)$ldHaploGroup$haploGroup$haploAssocTable)
  expect_gt(nrow(tbl), 1L)                                   # regression: rows without covariates
  base <- tbl[grepl("Ref", tbl$haplotype), ]
  expect_close(num(base$effect), 1.0, tol = 1e-6)           # reference OR = 1
  tc <- tbl[tbl$haplotype == "T-C", ]
  expect_close(num(tc$effect), 1.014, tol = 0.002)   # full-sample EM (haplo.em handles partial missingness); complete-case gave 1.008
})
