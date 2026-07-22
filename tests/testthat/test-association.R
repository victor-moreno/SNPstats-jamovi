# Tab 2: Association — every genetic model, response type and covariate setting
# is cross-checked against an independent glm / lm / multinom fit.

.models  <- c("codominant", "dominant", "recessive", "overdominant", "logadditive")
.modelOpt <- c(codominant = "modelCodominant", dominant = "modelDominant",
               recessive = "modelRecessive", overdominant = "modelOverdominant",
               logadditive = "modelLogAdditive")

# Run the module with exactly one genetic model enabled and return its assocTable.
run_assoc <- function(snp, response, type, covariates = NULL,
                      model = "logadditive", ci = 95, showAIC = FALSE) {
  flags <- as.list(setNames(rep(FALSE, length(.modelOpt)), unname(.modelOpt)))
  flags[[.modelOpt[[model]]]] <- TRUE
  args <- c(list(data = .test_data, snps = snp, response = response,
                 responseType = type, snpAssoc = TRUE, ciWidth = ci,
                 showAIC = showAIC), flags)
  if (!is.null(covariates)) args$covariates <- covariates
  res <- do.call(run_snp, args)
  res$assocGroup$assocSnpResults$get(key = snp)$assocTable
}

# Compare a module assocTable (one model) against the oracle coefficient rows.
# The module reports the effect / CI on the coefficient rows and one model-level
# p-value (LRT for binary/categorical, F for quantitative) on the model's first
# (labelled) row — matching the SNPstats reference tool. Per-genotype rows carry
# no p-value.
compare_assoc <- function(tbl, oracle, ptol = 0.01) {
  df   <- as_df(tbl)
  coef <- assoc_coef_rows(tbl)
  expect_equal(nrow(coef), nrow(oracle),
               label = "number of coefficient rows")
  for (i in seq_len(nrow(oracle))) {
    expect_close(num(coef$effect[i]), oracle$effect[i], label = paste0("effect[", i, "]"))
    expect_close(num(coef$ciLow[i]),  oracle$ciLow[i],  label = paste0("ciLow[", i, "]"))
    expect_close(num(coef$ciHigh[i]), oracle$ciHigh[i], label = paste0("ciHigh[", i, "]"))
  }
  model_p <- num(df$pval[nzchar(as.character(df$model))][1])
  gp <- attr(oracle, "global_p")
  if (is.na(model_p)) {          # module printed "< 0.001"
    expect_lt(gp, 0.001)
  } else {
    expect_close(model_p, gp, tol = ptol, label = "model LRT/F p")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Binary response (logistic) — all 5 models, with and without covariates
# ══════════════════════════════════════════════════════════════════════════════

test_that("binary assoc matches glm for every model (unadjusted)", {
  for (snp in .snps2) for (m in .models) {
    tbl <- run_assoc(snp, .resp, "binary", model = m)
    orc <- assoc_oracle(.test_data[[snp]], .test_data[[.resp]], "binary", model = m)
    compare_assoc(tbl, orc)
  }
})

test_that("binary assoc matches glm for every model (adjusted for sex+age+bmiOMS)", {
  covdf <- .test_data[, .covars]
  for (snp in .snps2) for (m in .models) {
    tbl <- run_assoc(snp, .resp, "binary", covariates = .covars, model = m)
    orc <- assoc_oracle(.test_data[[snp]], .test_data[[.resp]], "binary",
                        covdf = covdf, model = m)
    compare_assoc(tbl, orc)
  }
})

test_that("binary codominant reference-row p-value equals the global LRT", {
  for (snp in .snps2) {
    tbl <- as_df(run_assoc(snp, .resp, "binary", model = "codominant"))
    ref <- tbl[which(tbl$model == "Codominant"), ]
    orc <- assoc_oracle(.test_data[[snp]], .test_data[[.resp]], "binary",
                        model = "codominant")
    expect_close(num(ref$pval), attr(orc, "global_p"), tol = 0.01,
                 label = paste(snp, "global LRT"))
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# Quantitative response (linear) — effect = beta, CI = confint, p = t-test
# ══════════════════════════════════════════════════════════════════════════════

test_that("quantitative assoc matches lm for every model", {
  for (m in .models) {
    tbl <- run_assoc(.snps2[1], "age", "quantitative", model = m)
    orc <- assoc_oracle(.test_data[[.snps2[1]]], .test_data[["age"]],
                        "quantitative", model = m)
    compare_assoc(tbl, orc)
  }
})

test_that("quantitative assoc matches lm when adjusted", {
  covdf <- .test_data[, c("sex", "bmiOMS")]
  tbl <- run_assoc(.snps2[1], "age", "quantitative",
                   covariates = c("sex", "bmiOMS"), model = "logadditive")
  orc <- assoc_oracle(.test_data[[.snps2[1]]], .test_data[["age"]],
                      "quantitative", covdf = covdf, model = "logadditive")
  compare_assoc(tbl, orc)
})

# ══════════════════════════════════════════════════════════════════════════════
# Categorical response (multinomial) — one OR per non-reference category
# ══════════════════════════════════════════════════════════════════════════════

test_that("categorical assoc matches nnet::multinom (log-additive)", {
  tbl <- as_df(run_assoc(.snps2[1], "bmiOMS", "categorical", model = "logadditive"))
  orc <- assoc_oracle(.test_data[[.snps2[1]]], .test_data[["bmiOMS"]],
                      "categorical", model = "logadditive")
  coef <- tbl[!is.na(num(tbl$effect)) & tbl$genotype == "Per allele", ]
  expect_equal(nrow(coef), nrow(orc))
  for (i in seq_len(nrow(orc))) {
    expect_close(num(coef$effect[i]), orc$effect[i], tol = 0.005,
                 label = paste0("cat effect[", i, "]"))
    mp <- num(coef$pval[i])
    if (!is.na(mp)) expect_close(mp, orc$pval[i], tol = 0.02,
                                 label = paste0("cat pval[", i, "]"))
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# CI width and AIC / BIC
# ══════════════════════════════════════════════════════════════════════════════

test_that("ciWidth changes the interval to match confint at that level", {
  tbl <- run_assoc(.snps2[2], .resp, "binary", model = "logadditive", ci = 90)
  orc <- assoc_oracle(.test_data[[.snps2[2]]], .test_data[[.resp]], "binary",
                      model = "logadditive", ci = 90)
  compare_assoc(tbl, orc)
})

test_that("AIC and BIC match the fitted glm", {
  tbl <- as_df(run_assoc(.snps2[1], .resp, "binary", model = "logadditive",
                         showAIC = TRUE))
  orc <- assoc_oracle(.test_data[[.snps2[1]]], .test_data[[.resp]], "binary",
                      model = "logadditive")
  expect_true(all(c("AIC", "BIC") %in% names(tbl)))
  expect_close(num(tbl$AIC[1]), attr(orc, "aic"), tol = 0.05)
  expect_close(num(tbl$BIC[1]), attr(orc, "bic"), tol = 0.05)
})

# ══════════════════════════════════════════════════════════════════════════════
# SNP x covariate interaction (multiplicative) — verify the interaction term
# ══════════════════════════════════════════════════════════════════════════════

test_that("multiplicative interaction term matches glm(resp ~ snp * sex)", {
  snp <- .snps2[1]
  res <- run_snp(data = .test_data, snps = snp, response = .resp,
                 covariates = "sex", responseType = "binary",
                 snpInteraction = TRUE, interactionType = "multiplicative",
                 modelLogAdditive = TRUE, modelCodominant = FALSE,
                 showInteractionTable = TRUE)
  item <- res$assocGroup$assocSnpResults$get(key = snp)
  itbl <- as_df(item$interactionTable)

  g   <- dose_minor(.test_data[[snp]])
  df  <- data.frame(resp = bin01(.test_data[[.resp]]), snp = g, sex = .test_data$sex)
  df  <- df[complete.cases(df), ]
  fit <- glm(resp ~ snp * sex, data = df, family = binomial())
  co  <- summary(fit)$coefficients
  irow <- grep(":", rownames(co))
  or_int <- exp(co[irow, "Estimate"])

  mod_int <- num(itbl$effect[grepl(":", itbl$term)])
  expect_equal(length(mod_int), 1L)
  expect_close(mod_int, or_int, tol = 0.005)
})

test_that("conditional-on-covariate interaction matches glm(resp ~ sex/snp)", {
  snp <- .snps2[1]
  res <- run_snp(data = .test_data, snps = snp, response = .resp,
                 covariates = "sex", responseType = "binary",
                 snpInteraction = TRUE, interactionType = "conditional_on_covar",
                 modelLogAdditive = TRUE, modelCodominant = FALSE,
                 showInteractionTable = TRUE)
  itbl <- as_df(res$assocGroup$assocSnpResults$get(key = snp)$interactionTable)

  df  <- data.frame(resp = bin01(.test_data[[.resp]]),
                    snp = dose_minor(.test_data[[snp]]), sex = .test_data$sex)
  df  <- df[complete.cases(df), ]
  co  <- summary(glm(resp ~ sex/snp, data = df, family = binomial()))$coefficients
  or  <- exp(co[grep("snp", rownames(co)), "Estimate"])   # snp within each sex level

  mod <- num(itbl$effect[grepl(":snp$|snp$", itbl$term) & grepl(":", itbl$term)])
  expect_equal(length(mod), length(or))
  for (i in seq_along(or)) expect_close(sort(mod)[i], sort(or)[i], tol = 0.005)
})

# ══════════════════════════════════════════════════════════════════════════════
# Response-type auto-detection
# ══════════════════════════════════════════════════════════════════════════════

test_that("auto-detect gives the same result as explicit response types", {
  auto <- .golden <- NULL
  a <- as_df(run_assoc(.snps2[1], .resp, "auto", model = "logadditive"))
  b <- as_df(run_assoc(.snps2[1], .resp, "binary", model = "logadditive"))
  expect_equal(num(a$effect), num(b$effect))     # phenotype -> binary

  a2 <- as_df(run_assoc(.snps2[1], "age", "auto", model = "logadditive"))
  b2 <- as_df(run_assoc(.snps2[1], "age", "quantitative", model = "logadditive"))
  expect_equal(num(a2$effect), num(b2$effect))   # age -> quantitative
})

# ══════════════════════════════════════════════════════════════════════════════
# Stratified interaction tables — structural (populate without error)
# ══════════════════════════════════════════════════════════════════════════════

test_that("stratified-by-covariate / by-genotype / cross-class tables populate", {
  snp <- .snps2[1]
  res <- run_snp(data = .test_data, snps = snp, response = .resp,
                 covariates = "sex", responseType = "binary",
                 snpInteraction = TRUE, modelCodominant = TRUE,
                 showStratByCovariate = TRUE, showStratByGenotype = TRUE,
                 showCrossClassTable = TRUE)
  item <- res$assocGroup$assocSnpResults$get(key = snp)
  expect_gt(nrow(as_df(item$stratByCovariate)), 0L)
  expect_gt(nrow(as_df(item$stratByGenotype)), 0L)
  expect_gt(nrow(as_df(item$crossClassTable)), 0L)
})

# A zero-observation continuous cell must render the same em dash in every
# stratified table (stratByGenotype used to call fmt_cont directly and print a
# bare "NA"). Force it: the rare C/C genotype is blanked for one sex, so its cell
# in the other sex has no observations.
test_that("empty continuous cells render '—' consistently across stratified tables", {
  d <- .test_data
  g <- as.character(d[["rs72647484"]])
  d[["rs72647484"]][g == "C/C" & d$sex == "Female"] <- NA
  res <- run_snp(data = d, snps = "rs72647484", response = "age",
                 responseType = "quantitative", covariates = "sex",
                 snpInteraction = TRUE, interactionType = "conditional_on_covar",
                 modelCodominant = TRUE, showStratByCovariate = TRUE,
                 showStratByGenotype = TRUE, showCrossClassTable = TRUE)
  item <- res$assocGroup$assocSnpResults$get(key = "rs72647484")
  dash <- intToUtf8(0x2014)
  # No stratified table may show a bare "NA" for an empty continuous cell.
  for (nm in c("stratByCovariate", "stratByGenotype", "crossClassTable")) {
    s0 <- as_df(item[[nm]])$stat0
    expect_false(any(s0 == "NA", na.rm = TRUE), label = paste0(nm, " has a bare NA cell"))
  }
  # The previously-broken table (and the cross-classification) render the em dash.
  expect_true(any(as_df(item$stratByGenotype)$stat0 == dash, na.rm = TRUE))
  expect_true(any(as_df(item$crossClassTable)$stat0 == dash, na.rm = TRUE))
})
