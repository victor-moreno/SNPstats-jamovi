# snpPGS — polygenic score. Scores, scaling, missing handling, association and
# percentiles are cross-checked against an independent base-R reimplementation.

# ── Fixture: a small PGS-Catalog weights file for four dataset SNPs ───────────
# effect_allele = minor allele, other_allele = major (kept stable/hard-coded so
# the fixture does not depend on frequency recomputation).
.pgs_snps    <- c("rs12080929", "rs10911251", "rs10936599", "rs6691170")
.pgs_effect  <- c(rs12080929 = "C", rs10911251 = "C", rs10936599 = "T", rs6691170 = "T")
.pgs_other   <- c(rs12080929 = "T", rs10911251 = "A", rs10936599 = "C", rs6691170 = "G")
.pgs_weights <- c(rs12080929 = 0.5, rs10911251 = -0.3, rs10936599 = 0.8, rs6691170 = 0.2)

.pgs_weightsfile <- local({
  f <- tempfile(fileext = ".tsv")
  writeLines("# test PGS weights file", f)
  df <- data.frame(rsID = .pgs_snps,
                   effect_allele = .pgs_effect[.pgs_snps],
                   other_allele  = .pgs_other[.pgs_snps],
                   effect_weight = .pgs_weights[.pgs_snps],
                   chr_name = seq_along(.pgs_snps),
                   chr_position = seq_along(.pgs_snps) * 100L)
  suppressWarnings(write.table(df, f, sep = "\t", row.names = FALSE,
                               quote = FALSE, append = TRUE))
  f
})

# Effect-allele dosage matrix (rows = individuals, NA where genotype missing).
.pgs_dosage <- function() {
  sapply(.pgs_snps, function(s) {
    ea <- .pgs_effect[[s]]
    g  <- as.character(.test_data[[s]]); g[grepl("0", g)] <- NA
    p  <- strsplit(g, "/", fixed = TRUE)
    vapply(p, function(x)
      if (length(x) != 2 || any(is.na(x))) NA_real_ else sum(x == ea), numeric(1))
  })
}

# Independent score oracle replicating .computeScores for the options we test.
pgs_oracle <- function(unweighted = FALSE, scale = "proportion", factor = 10,
                       corrected = TRUE, standardize = FALSE,
                       missing = "SNP-wise") {
  D   <- .pgs_dosage()
  w   <- if (unweighted) rep(1, ncol(D)) else .pgs_weights[.pgs_snps]
  obs <- !is.na(D)
  Dimp <- D
  if (missing == "mean") {
    for (j in seq_len(ncol(D))) Dimp[!obs[, j], j] <- mean(D[, j], na.rm = TRUE)
  } else {
    Dimp[!obs] <- 0                       # SNP-wise / zero
  }
  keep_row <- if (missing == "exclude") rowSums(!obs) == 0 else rep(TRUE, nrow(D))

  num <- as.numeric(Dimp %*% w)
  den <- if (unweighted) 2 * rowSums(obs) else 2 * as.numeric(obs %*% pmax(w, 0))
  den[den == 0] <- NA
  score <- if (corrected) num / den else num
  score <- switch(scale,
                  percent  = score * 100,
                  multiply = score * factor,
                  score)
  if (standardize) score <- score / sd(score[keep_row], na.rm = TRUE)
  score[keep_row]
}

smry <- function(res, type) {
  st <- as_df(res$summaryTable)
  st[st$score_type == type & st$group == "Overall", ]
}

# ══════════════════════════════════════════════════════════════════════════════
# Scoring — summary statistics vs the oracle
# ══════════════════════════════════════════════════════════════════════════════

test_that("weighted proportion score matches the oracle", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                 weightsPath = .pgs_weightsfile, weightingMode = "weighted")
  o <- pgs_oracle()
  r <- smry(res, "Weighted")
  expect_equal(as.integer(r$n), sum(!is.na(o)))
  expect_close(num(r$mean), mean(o, na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$sd),   sd(o,   na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$min),  min(o,  na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$max),  max(o,  na.rm = TRUE), tol = 5e-4)
})

test_that("unweighted proportion score matches the oracle", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                 weightsPath = .pgs_weightsfile, weightingMode = "unweighted")
  o <- pgs_oracle(unweighted = TRUE)
  r <- smry(res, "Unweighted")
  expect_close(num(r$mean), mean(o, na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$sd),   sd(o,   na.rm = TRUE), tol = 5e-4)
})

test_that("scale methods (none/percent/multiply) match the oracle", {
  base <- list(data = .test_data, snpCols = .pgs_snps,
               weightsPath = .pgs_weightsfile, weightingMode = "weighted")
  r_none <- smry(do.call(run_pgs, c(base, scaleMethod = "none", missingCorrection = FALSE)), "Weighted")
  expect_close(num(r_none$mean), mean(pgs_oracle(scale = "none", corrected = FALSE), na.rm = TRUE), tol = 5e-4)

  r_pct <- smry(do.call(run_pgs, c(base, scaleMethod = "percent")), "Weighted")
  expect_close(num(r_pct$mean), mean(pgs_oracle(scale = "percent"), na.rm = TRUE), tol = 5e-3)

  r_mul <- smry(do.call(run_pgs, c(base, scaleMethod = "multiply", scaleFactor = 5)), "Weighted")
  expect_close(num(r_mul$mean), mean(pgs_oracle(scale = "multiply", factor = 5), na.rm = TRUE), tol = 5e-4)
})

test_that("standardize gives SD = 1 and the oracle mean", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", scaleMethod = "proportion", standardize = TRUE)
  r <- smry(res, "Weighted")
  expect_close(num(r$sd), 1, tol = 1e-3)
  expect_close(num(r$mean), mean(pgs_oracle(standardize = TRUE), na.rm = TRUE), tol = 5e-4)
})

test_that("missing-genotype strategies behave correctly", {
  base <- list(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
               weightingMode = "weighted", scaleMethod = "none", missingCorrection = FALSE)
  n_full <- as.integer(smry(do.call(run_pgs, c(base, missingStrategy = "SNP-wise")), "Weighted")$n)
  n_excl <- as.integer(smry(do.call(run_pgs, c(base, missingStrategy = "exclude")),  "Weighted")$n)
  expect_equal(n_full, nrow(.test_data))
  expect_lt(n_excl, n_full)               # individuals with any missing SNP dropped
  expect_equal(n_excl, sum(rowSums(is.na(.pgs_dosage())) == 0))
})

# ══════════════════════════════════════════════════════════════════════════════
# Association — vs glm / lm
# ══════════════════════════════════════════════════════════════════════════════

test_that("logistic PGS-response association matches glm (Wald CI)", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", responseCol = "phenotype", showAssoc = TRUE)
  at <- as_df(res$assocTable)
  lr <- at[at$test == "Logistic regression" & at$score_type == "Weighted", ]

  o  <- pgs_oracle()
  y  <- bin01(.test_data$phenotype); cc <- !is.na(o) & !is.na(y)
  fit <- glm(y[cc] ~ o[cc], family = binomial())
  co  <- summary(fit)$coefficients
  ci  <- suppressWarnings(confint.default(fit))
  expect_close(num(lr$estimate), exp(co[2, 1]), tol = 5e-4)
  expect_close(num(lr$ci_low),   exp(ci[2, 1]), tol = 5e-4)
  expect_close(num(lr$ci_high),  exp(ci[2, 2]), tol = 5e-4)
  expect_close(num(lr$p),        co[2, 4],      tol = 1e-3)
})

test_that("character response honors data order of appearance, not alphabetical", {
  # A character response carries no stored levels; the reference must be the
  # first OBSERVED level (data order), not R's alphabetical factor() default.
  # Reorder so Control appears first (alphabetical would make Case the reference).
  d <- .test_data[order(.test_data$phenotype != "Control"), ]
  d$rc <- as.character(d$phenotype)
  res  <- run_pgs(data = d, snpCols = .pgs_snps, responseCol = "rc", showAssoc = TRUE)
  note <- res$assocTable$notes[["respNote"]]$note
  expect_match(note, "Case vs Control", fixed = TRUE)   # Control (first seen) = reference
})

test_that("linear PGS-response association matches lm", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", responseCol = "age", showAssoc = TRUE)
  at <- as_df(res$assocTable)
  lr <- at[at$test == "Linear regression" & at$score_type == "Weighted", ]

  o  <- pgs_oracle(); cc <- !is.na(o) & !is.na(.test_data$age)
  fit <- lm(.test_data$age[cc] ~ o[cc]); co <- summary(fit)$coefficients; ci <- confint(fit)
  expect_close(num(lr$estimate), co[2, 1], tol = 5e-4)
  expect_close(num(lr$ci_low),   ci[2, 1], tol = 5e-3)
  expect_close(num(lr$p),        co[2, 4], tol = 1e-3)
})

test_that("PGS x covariate interaction matches glm and populates (regression)", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", responseCol = "phenotype",
                 covCols = "sex", showInteraction = TRUE)
  it <- as_df(res$interactionTable)
  expect_gt(nrow(it), 0L)                       # regression: was empty before the gate fix

  o  <- pgs_oracle()
  df <- data.frame(y = factor(.test_data$phenotype), pgs = o, sex = .test_data$sex)
  df <- df[complete.cases(df), ]
  co <- summary(glm(y ~ pgs * sex, data = df, family = binomial()))$coefficients

  expect_close(num(it$estimate[it$term == "PGS (main)"]),       exp(co["pgs", 1]),          tol = 2e-3)
  expect_close(num(it$estimate[grepl("×", it$term)]),      exp(co["pgs:sexMale", 1]),  tol = 2e-3)
  expect_close(num(it$p[grepl("×", it$term)]),             co["pgs:sexMale", 4],       tol = 2e-3)
})

# ══════════════════════════════════════════════════════════════════════════════
# Other tables and the documented calling convention
# ══════════════════════════════════════════════════════════════════════════════

test_that("documented no-response call works (regression for missing formal defaults)", {
  expect_error(run_pgs(data = .test_data, snpCols = .pgs_snps,
                       weightsPath = .pgs_weightsfile), NA)
  expect_error(run_pgs(data = .test_data, snpCols = .pgs_snps), NA)   # unweighted fallback
})

test_that("percentile category counts sum to N with monotonic ranges", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", responseCol = "phenotype",
                 showPercentiles = TRUE, percentileBreaks = "20,40,60,80")
  tt <- as_df(res$percentileThreshTable)
  counts <- as.integer(sub("\\s*\\(.*$", "", tt$n_overall))
  expect_equal(length(counts), 5L)                      # 4 breaks -> 5 bands
  expect_equal(sum(counts), sum(!is.na(pgs_oracle())))
})

test_that("coverage and SNP-grid tables report correct AF / matching", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", showCoverage = TRUE, showSnpGrid = TRUE)
  grid <- as_df(res$snpGridTable)
  grid <- grid[grid$rsid %in% .pgs_snps, ]
  expect_equal(nrow(grid), length(.pgs_snps))
  # effect-allele frequency from the grid matches the dosage matrix
  D <- .pgs_dosage()
  for (s in .pgs_snps) {
    af_oracle <- mean(D[, s], na.rm = TRUE) / 2
    expect_close(num(grid$effect_af[grid$rsid == s]), af_oracle, tol = 2e-3,
                 label = paste(s, "effect AF"))
  }
})

test_that("plots render without error (incl. calibration with tied predictions)", {
  skip_if_not_installed("ggplot2")
  grDevices::png(tempfile()); on.exit(grDevices::dev.off())
  # weak score + binary covariate -> tied predicted probabilities, which used to
  # crash the calibration plot with "'breaks' are not unique".
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "both", responseCol = "phenotype", covCols = "sex",
                 showDistPlot = TRUE, showRocPlot = TRUE, showCalibPlot = TRUE,
                 showForestPlot = TRUE)
  for (p in c("distPlot", "rocPlot", "calibPlot", "forestPlot"))
    expect_error(res[[p]]$.render(), NA, label = p)
})

test_that("plot image size is set from plotWidth/plotHeight in .init", {
  # Regression: size was only set inside the render function (after jamovi had
  # already created the device at the default 400x300), so the first render was
  # upscaled/blurry and ignored the options. It is now set in .init.
  opt <- snpPGSOptions$new(snpCols = as.list(.pgs_snps), responseCol = "phenotype",
           weightingMode = "unweighted", showDistPlot = TRUE, showRocPlot = TRUE,
           plotWidth = 800, plotHeight = 500)
  a <- snpPGSClass$new(options = opt, data = .test_data, analysisId = 1, revision = 1)
  a$init()                                  # no run/render yet
  expect_equal(c(a$results$distPlot$width, a$results$distPlot$height), c(800, 500))
  expect_equal(c(a$results$rocPlot$width,  a$results$rocPlot$height),  c(800, 500))
})

test_that("plot visibility is set in .init so .run does not re-touch (no re-render)", {
  # Regression: forest/roc/calib visibility was set in .run every run, and an
  # image touched in .run is re-rendered by the engine — so toggling any plot
  # re-rendered the others (slow, they refit models). Visibility is now predicted
  # in .init; .run only corrects a genuine mismatch. Here the prediction is exact,
  # so .run's guard is a no-op (the plots are never touched in .run).
  chk <- function(respCol, exp) {
    opt <- snpPGSOptions$new(snpCols = as.list(.pgs_snps), responseCol = respCol,
             weightingMode = "unweighted", showDistPlot = TRUE, showForestPlot = TRUE,
             showRocPlot = TRUE, showCalibPlot = TRUE)
    a <- snpPGSClass$new(options = opt, data = .test_data, analysisId = 1, revision = 1)
    a$init()
    vis0 <- a$.__enclos_env__$private$.plotVis
    a$run()
    vis1 <- a$.__enclos_env__$private$.plotVis
    expect_identical(vis0, vis1)                 # .run guard was a no-op
    expect_equal(vis1[names(exp)], exp)
  }
  # binary phenotype: forest/roc/calib shown, scatter (strat) hidden
  chk("phenotype", list(stratPlot = FALSE, forestPlot = TRUE,
                        rocPlot = TRUE, calibPlot = TRUE))
  # continuous age: scatter + forest shown, roc/calib hidden (not categorical)
  chk("age", list(stratPlot = TRUE, forestPlot = TRUE,
                  rocPlot = FALSE, calibPlot = FALSE))
})

test_that("plot clearWith lists only the options each plot depends on", {
  # A plot re-renders when any clearWith option changes, so listing options it
  # does not use makes it refresh on unrelated clicks. distPlot/stratPlot draw
  # scores by response and use neither covariates nor percentile options;
  # stratPlot (continuous) is also independent of caseLevel. forest/roc/calib do
  # use covariates and the case level.
  opt <- snpPGSOptions$new(snpCols = list("rs12080929"))
  a <- snpPGSClass$new(options = opt,
                       data = data.frame(rs12080929 = factor("A/A")),
                       analysisId = 1, revision = 1)
  cw <- function(nm) unlist(a$results$get(nm)$.__enclos_env__$private$.clearWith)

  expect_false("covCols"          %in% cw("distPlot"))
  expect_true ("caseLevel"        %in% cw("distPlot"))
  expect_false(any(c("covCols", "percentileBreaks", "pgsRefCategory", "caseLevel")
                   %in% cw("stratPlot")))
  for (nm in c("forestPlot", "rocPlot", "calibPlot")) {
    expect_true("covCols"   %in% cw(nm), label = paste(nm, "covCols"))
    expect_true("caseLevel" %in% cw(nm), label = paste(nm, "caseLevel"))
  }
  expect_true(all(c("percentileBreaks", "pgsRefCategory") %in% cw("forestPlot")))
})

test_that("reference-level tables clear on caseLevel", {
  # caseLevel relevels the response reference, which flips the binary OR/t-test
  # direction and sets the polytomous baseline. Tables whose numbers depend on it
  # must list caseLevel in clearWith or a reference change silently shows stale
  # estimates (the gated fill is skipped when isNotFilled() stays FALSE).
  opt <- snpPGSOptions$new(snpCols = list("rs12080929"))
  a <- snpPGSClass$new(options = opt,
                       data = data.frame(rs12080929 = factor("A/A")),
                       analysisId = 1, revision = 1)
  cw <- function(nm) unlist(a$results$get(nm)$.__enclos_env__$private$.clearWith)
  for (nm in c("assocTable", "interactionTable", "percentileTable", "summaryTable"))
    expect_true("caseLevel" %in% cw(nm), label = paste(nm, "caseLevel"))
})

test_that("distPlotType switches the distribution plot geometry", {
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  geoms <- function(ptype) {
    r <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", responseCol = "phenotype",
                 showDistPlot = TRUE, distPlotType = ptype)
    r$distPlot$.render()
    lp <- ggplot2::last_plot()
    vapply(lp$layers, function(L) class(L$geom)[1], "")
  }
  expect_true("GeomDensity" %in% geoms("density") && !("GeomBar" %in% geoms("density")))
  expect_true("GeomBar"     %in% geoms("histogram") && !("GeomDensity" %in% geoms("histogram")))
  expect_true(all(c("GeomBar", "GeomDensity") %in% geoms("both")))
})

test_that("skewness returns NA (not NaN) for a constant score", {
  expect_true(is.na(SNPstats:::skewness(rep(0.5, 100))))   # zero variance
  expect_true(is.na(SNPstats:::skewness(c(1, 2))))         # n < 3
  expect_gt(SNPstats:::skewness(c(0, 0, 0, 0, 10)), 0)     # right-skewed -> positive
})

test_that("QC missingness filter excludes SNPs above threshold", {
  D <- .pgs_dosage()
  pct_miss <- colMeans(is.na(D)) * 100
  thr <- min(pct_miss) + (max(pct_miss) - min(pct_miss)) / 2   # between min and max
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "weighted", qcFilterMissing = TRUE, qcMaxMissingPct = thr)
  grid <- as_df(res$snpGridTable)
  excluded <- grid$rsid[grepl("excl \\(missing", grid$allele_status)]
  expect_setequal(excluded, names(pct_miss)[pct_miss > thr])
})

# ══════════════════════════════════════════════════════════════════════════════
# Golden values (verified above; detect future R/package changes)
# ══════════════════════════════════════════════════════════════════════════════

test_that("GOLDEN pgs scores and association", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsPath = .pgs_weightsfile,
                 weightingMode = "both", responseCol = "phenotype", showAssoc = TRUE)
  w <- smry(res, "Weighted"); u <- smry(res, "Unweighted")
  expect_close(num(w$mean), 0.17651, tol = 5e-4)
  expect_close(num(w$sd),   0.20799, tol = 5e-4)
  expect_close(num(u$mean), 0.31730, tol = 5e-4)
  at <- as_df(res$assocTable)
  lr <- at[at$test == "Logistic regression" & at$score_type == "Weighted", ]
  expect_close(num(lr$estimate), 1.3327, tol = 2e-3)
  expect_close(num(lr$ci_low),   0.9304, tol = 2e-3)
  expect_close(num(lr$ci_high),  1.9089, tol = 2e-3)
})
