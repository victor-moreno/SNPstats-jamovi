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
    Dimp[!obs] <- 0                       # SNP-wise / zero both impute 0
  }
  keep_row <- if (missing == "exclude") rowSums(!obs) == 0 else rep(TRUE, nrow(D))

  # Which SNPs count toward the missingness-correction denominator. Only
  # 'SNP-wise' uses an observed-only denominator — that is what defines it.
  # 'zero' asserts the missing genotype is a real dosage 0 and 'mean' substitutes
  # a dosage the numerator then uses, so under both the SNP counts. (Getting this
  # wrong is what made 'zero' identical to 'SNP-wise' and let 'mean' produce
  # "proportion of maximum" scores above 1.)
  cnt <- if (missing %in% c("zero", "mean")) matrix(TRUE, nrow(obs), ncol(obs)) else obs

  num <- as.numeric(Dimp %*% w)
  # Attainable bounds over the SNPs that count. A negative weight is the positive
  # weight |w| on the other allele, so it widens the range just as much as a
  # positive one; the corrected score is the position within [lo, hi]. Using
  # 2*sum(pmax(w,0)) instead (the old formula) drops negative-weight SNPs from
  # the denominator, which puts "proportion" scores outside [0,1] and can leave
  # an individual with a denominator of 0. For all-positive weights lo is 0 and
  # this reduces to the old expression; unweighted (w = 1) likewise.
  lo  <- 2 * as.numeric(cnt %*% pmin(w, 0))
  rng <- 2 * as.numeric(cnt %*% abs(w))
  rng[rng == 0] <- NA
  score <- if (corrected) (num - lo) / rng else num
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
                 weightsFile = .pgs_weightsfile, weightingMode = "weighted")
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
                 weightsFile = .pgs_weightsfile, weightingMode = "unweighted")
  o <- pgs_oracle(unweighted = TRUE)
  r <- smry(res, "Unweighted")
  expect_close(num(r$mean), mean(o, na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$sd),   sd(o,   na.rm = TRUE), tol = 5e-4)
})

test_that("hand-rolled base64 weights match the pgs_weights() helper", {
  b64 <- base64enc::base64encode(readBin(.pgs_weightsfile, "raw",
                                         file.info(.pgs_weightsfile)$size))
  res_helper <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                       weightsFile = .pgs_weightsfile, weightingMode = "weighted")
  res_embed <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                       weightsContent = b64, weightsFilename = "weights.tsv",
                       weightingMode = "weighted")
  expect_equal(as_df(res_embed$summaryTable), as_df(res_helper$summaryTable))
})

test_that("gzipped embedded weights match the uncompressed route", {
  raw <- readBin(.pgs_weightsfile, "raw", file.info(.pgs_weightsfile)$size)
  gz  <- memCompress(raw, "gzip")
  b64 <- base64enc::base64encode(gz)
  res_plain <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                       weightsFile = .pgs_weightsfile, weightingMode = "weighted")
  res_gz    <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                       weightsContent = b64, weightsFilename = "weights.tsv.gz",
                       weightingMode = "weighted")
  expect_equal(as_df(res_gz$summaryTable), as_df(res_plain$summaryTable))
})

test_that("orientation note is set only when no weights file is loaded", {
  note_txt <- function(res, tbl) {
    n <- res[[tbl]]$notes[["orientNote"]]
    if (is.null(n)) NULL else n$note
  }
  res_file <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                      weightsFile = .pgs_weightsfile, weightingMode = "weighted")
  res_nofile <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                        weightingMode = "unweighted")
  expect_null(note_txt(res_file, "summaryTable"))
  expect_match(note_txt(res_nofile, "summaryTable"), "risk", ignore.case = TRUE)
  expect_match(note_txt(res_nofile, "assocTable"), "frequency", ignore.case = TRUE)
})

test_that("no-weights effect allele defaults to the minor allele (snpStats-consistent)", {
  res  <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                  weightingMode = "unweighted", showSnpGrid = TRUE)
  grid <- as_df(res$snpGridTable)
  for (s in .pgs_snps) {
    g <- as.character(.test_data[[s]]); g[grepl("0", g)] <- NA
    b <- unlist(strsplit(g[!is.na(g)], "/", fixed = TRUE))
    b <- b[b %in% c("A", "C", "G", "T")]
    minor <- names(sort(table(b)))[1]                 # least frequent allele
    expect_equal(grid$effect_allele[grid$rsid == s], minor,
                 label = paste(s, "effect allele = data minor allele"))
  }
})

test_that("no-weights orientation respects a reordered genotype factor", {
  gts <- c(rep("G/G", 30), rep("A/G", 12), rep("A/A", 2))   # G major, A minor
  d   <- data.frame(snp1 = gts, stringsAsFactors = FALSE)
  # Alphabetical (auto) levels -> no user intent -> frequency default: effect = A
  res_def <- run_pgs(data = d, snpCols = "snp1",
                     weightingMode = "unweighted", showSnpGrid = TRUE)
  g_def   <- as_df(res_def$snpGridTable)
  expect_equal(g_def$effect_allele[g_def$rsid == "snp1"], "A")
  # User reorders levels non-alphabetically so A/A is the reference homozygote;
  # the effect allele then flips to the major allele G.
  d$snp1  <- factor(d$snp1, levels = c("A/A", "G/G", "A/G"))
  res_usr <- run_pgs(data = d, snpCols = "snp1",
                     weightingMode = "unweighted", showSnpGrid = TRUE)
  g_usr   <- as_df(res_usr$snpGridTable)
  expect_equal(g_usr$effect_allele[g_usr$rsid == "snp1"], "G")
})

test_that("scale methods (none/percent/multiply) match the oracle", {
  base <- list(data = .test_data, snpCols = .pgs_snps,
               weightsFile = .pgs_weightsfile, weightingMode = "weighted")
  r_none <- smry(do.call(run_pgs, c(base, scaleMethod = "none", missingCorrection = FALSE)), "Weighted")
  expect_close(num(r_none$mean), mean(pgs_oracle(scale = "none", corrected = FALSE), na.rm = TRUE), tol = 5e-4)

  r_pct <- smry(do.call(run_pgs, c(base, scaleMethod = "percent")), "Weighted")
  expect_close(num(r_pct$mean), mean(pgs_oracle(scale = "percent"), na.rm = TRUE), tol = 5e-3)

  r_mul <- smry(do.call(run_pgs, c(base, scaleMethod = "multiply", scaleFactor = 5)), "Weighted")
  expect_close(num(r_mul$mean), mean(pgs_oracle(scale = "multiply", factor = 5), na.rm = TRUE), tol = 5e-4)
})

test_that("standardize gives SD = 1 and the oracle mean", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
                 weightingMode = "weighted", scaleMethod = "proportion", standardize = TRUE)
  r <- smry(res, "Weighted")
  expect_close(num(r$sd), 1, tol = 1e-3)
  expect_close(num(r$mean), mean(pgs_oracle(standardize = TRUE), na.rm = TRUE), tol = 5e-4)
})

test_that("missing-genotype strategies behave correctly", {
  base <- list(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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
                       weightsFile = .pgs_weightsfile), NA)
  expect_error(run_pgs(data = .test_data, snpCols = .pgs_snps), NA)   # unweighted fallback
})

test_that("percentile category counts sum to N with monotonic ranges", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
                 weightingMode = "weighted", responseCol = "phenotype",
                 showPercentiles = TRUE, percentileBreaks = "20,40,60,80")
  tt <- as_df(res$percentileThreshTable)
  counts <- as.integer(sub("\\s*\\(.*$", "", tt$n_overall))
  expect_equal(length(counts), 5L)                      # 4 breaks -> 5 bands
  expect_equal(sum(counts), sum(!is.na(pgs_oracle())))
})

test_that("coverage and SNP-grid tables report correct AF / matching", {
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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

test_that("coverage reports count of SNPs whose risk allele is the major allele", {
  cov_val <- function(res) {
    cv  <- as_df(res$coverageTable)
    row <- cv[grepl("Risk allele", cv$field), ]
    as.integer(row$value[1])
  }
  # Fixture effect_allele = minor allele -> effect_af <= 0.5 -> count 0.
  res0 <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
                  weightingMode = "weighted", showCoverage = TRUE)
  expect_equal(cov_val(res0), 0L)

  # Swapped file: effect_allele = major allele -> effect_af > 0.5 for every SNP.
  fswap <- tempfile(fileext = ".tsv")
  writeLines("# swapped effect/other alleles", fswap)
  dsw <- data.frame(rsID = .pgs_snps,
                    effect_allele = .pgs_other[.pgs_snps],   # major
                    other_allele  = .pgs_effect[.pgs_snps],  # minor
                    effect_weight = .pgs_weights[.pgs_snps])
  suppressWarnings(write.table(dsw, fswap, sep = "\t", row.names = FALSE,
                               quote = FALSE, append = TRUE))
  res1 <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = fswap,
                  weightingMode = "weighted", showCoverage = TRUE)
  # Independent expectation: SNPs whose major-allele frequency exceeds 0.5.
  D <- .pgs_dosage()
  expected <- sum(vapply(.pgs_snps, function(s) (1 - mean(D[, s], na.rm = TRUE) / 2) > 0.5,
                         logical(1)))
  expect_equal(cov_val(res1), as.integer(expected))
})

test_that("plots render without error (incl. calibration with tied predictions)", {
  skip_if_not_installed("ggplot2")
  grDevices::png(tempfile()); on.exit(grDevices::dev.off())
  # weak score + binary covariate -> tied predicted probabilities, which used to
  # crash the calibration plot with "'breaks' are not unique".
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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
    r <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
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
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
                 weightingMode = "weighted", qcFilterMissing = TRUE, qcMaxMissingPct = thr)
  grid <- as_df(res$snpGridTable)
  excluded <- grid$rsid[grepl("excl \\(missing", grid$allele_status)]
  expect_setequal(excluded, names(pct_miss)[pct_miss > thr])
})

# ══════════════════════════════════════════════════════════════════════════════
# Golden values (verified above; detect future R/package changes)
# ══════════════════════════════════════════════════════════════════════════════

test_that("GOLDEN pgs scores and association", {
  # The WEIGHTED values here were re-pinned when the missingness correction began
  # treating a negative weight as the positive weight on the other allele
  # (was mean 0.17651 / sd 0.20799, OR 1.3327 [0.9304, 1.9089]). The fixture
  # carries one negative weight, rs10911251 = -0.3, which the old denominator
  # dropped: the corrected score was raw / 2*sum(pmax(w,0)), so it was divided by
  # a maximum that ignored that SNP. The new values were verified against the
  # independent oracle in this file and against a hand-oriented score (flip the
  # negative-weight SNP to its other allele, use |w|), which they reproduce
  # exactly. The UNWEIGHTED golden is unchanged, as it must be — with all weights
  # equal to 1 the two formulas coincide.
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = .pgs_weightsfile,
                 weightingMode = "both", responseCol = "phenotype", showAssoc = TRUE)
  w <- smry(res, "Weighted"); u <- smry(res, "Unweighted")
  expect_close(num(w$mean), 0.31355, tol = 5e-4)
  expect_close(num(w$sd),   0.17333, tol = 5e-4)
  expect_close(num(u$mean), 0.31730, tol = 5e-4)
  at <- as_df(res$assocTable)
  lr <- at[at$test == "Logistic regression" & at$score_type == "Weighted", ]
  expect_close(num(lr$estimate), 1.3997, tol = 2e-3)
  expect_close(num(lr$ci_low),   0.9095, tol = 2e-3)
  expect_close(num(lr$ci_high),  2.1541, tol = 2e-3)
})

# ══════════════════════════════════════════════════════════════════════════════
# Weights source — content only, never a filesystem path
#
# Analysis options are serialised into the .omv and re-run when it is opened, so
# a path option would be resolved against the *opener's* filesystem. These pin
# the property that there is no path to resolve.
# ══════════════════════════════════════════════════════════════════════════════

test_that("snpPGS exposes no file-path option", {
  expect_false("weightsPath" %in% names(formals(SNPstats::snpPGS)))
  expect_false("weightsPath" %in% SNPstats::snpPGSOptions$new()$names)
})

test_that("pgs_weights round-trips a file into the content/name pair", {
  w <- SNPstats::pgs_weights(.pgs_weightsfile)
  expect_named(w, c("weightsContent", "weightsFilename"))
  expect_equal(w$weightsFilename, basename(.pgs_weightsfile))
  expect_equal(rawToChar(base64enc::base64decode(w$weightsContent)),
               paste0(paste(readLines(.pgs_weightsfile), collapse = "\n"), "\n"))
  expect_error(SNPstats::pgs_weights(file.path(tempdir(), "no-such-file.csv")),
               "file not found")
})

test_that("a parse failure reports expectations, not file content", {
  secret <- tempfile(fileext = ".csv")
  writeLines(c("TOP_SECRET_HEADER,another_private_column", "1,2"), secret)
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = secret,
                 weightingMode = "weighted")
  msg <- res$validationMsg$content

  # The actionable diagnosis names what to fix and what was expected...
  expect_match(msg, "rsID")
  expect_match(msg, "no recognisable", ignore.case = TRUE)
  # ...and survives the downstream, more generic message rather than being
  # overwritten by it (validation messages accumulate within a run).
  expect_match(msg, "No SNPs passed QC filters")
  expect_lt(regexpr("rsID", msg, fixed = TRUE),
            regexpr("No SNPs passed QC", msg, fixed = TRUE))

  # And nothing of the file's own content is echoed back, anywhere.
  all_out <- paste(capture.output(print(res)), collapse = "\n")
  expect_false(grepl("TOP_SECRET_HEADER", all_out, fixed = TRUE))
  expect_false(grepl("another_private_column", all_out, fixed = TRUE))
  expect_false(grepl("TOP_SECRET_HEADER", msg, fixed = TRUE))
})

test_that("validation messages do not leak between runs", {
  # .run() resets the accumulator; a clean run must not inherit the previous
  # run's complaint.
  secret <- tempfile(fileext = ".csv")
  writeLines(c("nothing_useful,at_all", "1,2"), secret)
  bad <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = secret,
                 weightingMode = "weighted")
  expect_match(bad$validationMsg$content, "rsID")

  ok <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                weightsFile = .pgs_weightsfile, weightingMode = "weighted")
  txt <- ok$validationMsg$content
  expect_false(grepl("rsID", if (is.null(txt)) "" else txt, fixed = TRUE))
})

test_that("an oversized gunzipped payload is refused rather than expanded", {
  # A '.gz' name makes .weightsRawLines gunzip the embedded bytes; the result is
  # rejected above PGS_MAX_WEIGHTS_BYTES so a crafted .omv cannot exhaust the
  # engine. Scored output falls back to unit weights.
  bomb <- memCompress(charToRaw(strrep("A", 70 * 1024^2)), "gzip")
  expect_true(length(bomb) < 1e6)                     # small on the wire
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                 weightsContent  = base64enc::base64encode(bomb),
                 weightsFilename = "bomb.tsv.gz",
                 weightingMode   = "weighted")
  expect_gt(nrow(as_df(res$summaryTable)), 0L)
})

# ══════════════════════════════════════════════════════════════════════════════
# Missing-genotype strategies
#
# Each of the four must do something distinct, and a "proportion of maximum"
# score must stay inside [0, 1]. 'zero' used to be byte-identical to 'SNP-wise'
# (both imputed 0 AND both kept the observed-only denominator, so the option did
# nothing), and 'mean' added an imputed dosage to the numerator whose SNP the
# denominator left out, pushing the score past 1.
# ══════════════════════════════════════════════════════════════════════════════

.pgs_gappy <- local({
  d <- .test_data
  set.seed(1)
  for (s in .pgs_snps) d[[s]][sample(nrow(d), 800)] <- NA
  d
})

test_that("each missing-genotype strategy matches the oracle", {
  for (ms in c("SNP-wise", "zero", "mean", "exclude")) {
    res <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                   weightsFile = .pgs_weightsfile, weightingMode = "weighted",
                   missingStrategy = ms)
    o <- pgs_oracle(missing = ms)
    r <- smry(res, "Weighted")
    expect_close(num(r$mean), mean(o, na.rm = TRUE), tol = 5e-4,
                 label = paste(ms, "mean"))
    expect_close(num(r$sd), sd(o, na.rm = TRUE), tol = 5e-4,
                 label = paste(ms, "sd"))
  }
})

test_that("'zero' and 'SNP-wise' are not the same strategy", {
  g <- function(ms) num(smry(run_pgs(data = .pgs_gappy, snpCols = .pgs_snps,
                                     weightsFile = .pgs_weightsfile,
                                     weightingMode = "weighted",
                                     missingStrategy = ms), "Weighted")$mean)
  snpwise <- g("SNP-wise"); zero <- g("zero")
  expect_false(isTRUE(all.equal(snpwise, zero)))
  # zero counts the un-typed SNPs in the denominator, so it can only be lower
  expect_lt(zero, snpwise)
})

test_that("proportion-scaled scores never exceed 1, for every strategy", {
  # The upper bound is the invariant: the correction divides by the maximum the
  # individual could have scored, so the ratio cannot exceed 1. ('mean' used to
  # reach ~3 because the imputed dosage was in the numerator but its SNP was not
  # in the denominator.)
  #
  # There is no matching lower bound for the WEIGHTED score: the denominator uses
  # positive weights only, so a negative effect weight — the fixture has one,
  # rs10911251 = -0.3 — legitimately drives the score below 0. Unweighted scores
  # are sums of non-negative dosages and cannot be.
  for (ms in c("SNP-wise", "zero", "mean", "exclude")) {
    for (wm in c("weighted", "unweighted")) {
      res <- run_pgs(data = .pgs_gappy, snpCols = .pgs_snps,
                     weightsFile = .pgs_weightsfile, weightingMode = wm,
                     missingStrategy = ms, scaleMethod = "proportion")
      st <- as_df(res$summaryTable)
      st <- st[st$group == "Overall", ]
      lab <- paste(ms, wm)
      expect_lte(max(num(st$max)), 1 + 1e-9, label = paste(lab, "max"))
      if (wm == "unweighted")
        expect_gte(min(num(st$min)), -1e-9, label = paste(lab, "min"))
    }
  }
})

test_that("HWE QC filter tests the group named by caseLevel, not the raw first level", {
  # Standard case-control QC tests HWE in CONTROLS: departure from equilibrium in
  # cases can be a real association signal, so filtering on cases discards true
  # positives. The filter used to read the response column's raw first factor
  # level and ignore caseLevel entirely — and with the usual "Case"/"Control"
  # labels the alphabetical first level is Case, so it ran on cases by default.
  G  <- asNamespace("SNPstats")
  pg <- get("parse_genotype", envir = G)
  hx <- get("snp_hwe_exact",  envir = G)
  cl <- get("clean_null_alleles", envir = G)

  hand <- function(level) vapply(.pgs_snps, function(s) {
    sub <- .test_data[!is.na(.test_data$phenotype) & .test_data$phenotype == level, ]
    gg  <- tryCatch(pg(cl(as.character(sub[[s]])), NULL), error = function(e) NULL)
    if (is.null(gg)) return(NA_real_)
    tryCatch(hx(gg)$p.value, error = function(e) NA_real_)
  }, 0)
  grid_hwe <- function(...) {
    r <- run_pgs(data = .test_data, snpCols = .pgs_snps, responseCol = "phenotype",
                 weightsFile = .pgs_weightsfile, qcFilterHwe = TRUE,
                 showSnpGrid = TRUE, ...)
    g <- as_df(r$snpGridTable)
    num(g$hwe_p[match(.pgs_snps, g$rsid)])
  }

  expect_equal(levels(.test_data$phenotype)[1], "Case")   # the trap this guards
  expect_equal(grid_hwe(caseLevel = "Control"), unname(hand("Control")), tolerance = 1e-8)
  expect_equal(grid_hwe(caseLevel = "Case"),    unname(hand("Case")),    tolerance = 1e-8)
  # caseLevel must actually change the answer — it used to be ignored
  expect_false(isTRUE(all.equal(grid_hwe(caseLevel = "Control"),
                                grid_hwe(caseLevel = "Case"))))
})

test_that("the group HWE was tested in is always stated, not only on an exclusion", {
  on_res <- run_pgs(data = .test_data, snpCols = .pgs_snps, responseCol = "phenotype",
                    weightsFile = .pgs_weightsfile, qcFilterHwe = TRUE,
                    showSnpGrid = TRUE)
  note <- on_res$snpGridTable$notes[["hweGroupNote"]]
  expect_false(is.null(note))
  expect_match(note$note, "phenotype")
  expect_match(note$note, "Case")          # the level actually used

  off_res <- run_pgs(data = .test_data, snpCols = .pgs_snps, responseCol = "phenotype",
                     weightsFile = .pgs_weightsfile, showSnpGrid = TRUE)
  n_off <- off_res$snpGridTable$notes[["hweGroupNote"]]
  expect_true(is.null(n_off) || is.null(n_off$note))
})

# ══════════════════════════════════════════════════════════════════════════════
# Negative effect weights
#
# A negative weight on the effect allele is the positive weight |w| on the OTHER
# allele, because dosage_other = 2 - dosage_effect. Such a SNP therefore widens
# the range of attainable scores exactly as much as a positive one. The
# correction denominator used to be 2*sum(pmax(w, 0)), which dropped those SNPs.
# ══════════════════════════════════════════════════════════════════════════════

.pgs_mixed_wfile <- local({
  f <- tempfile(fileext = ".tsv")
  writeLines("# mixed-sign weights", f)
  w  <- c(rs12080929 = 0.5, rs10911251 = -0.3, rs10936599 = 0.8, rs6691170 = -0.2)
  df <- data.frame(rsID = .pgs_snps,
                   effect_allele = .pgs_effect[.pgs_snps],
                   other_allele  = .pgs_other[.pgs_snps],
                   effect_weight = w[.pgs_snps],
                   chr_name = seq_along(.pgs_snps),
                   chr_position = seq_along(.pgs_snps) * 100L)
  suppressWarnings(write.table(df, f, sep = "\t", row.names = FALSE,
                               quote = FALSE, append = TRUE))
  f
})
.pgs_mixed_w <- c(rs12080929 = 0.5, rs10911251 = -0.3, rs10936599 = 0.8, rs6691170 = -0.2)

test_that("a negative weight scores as the positive weight on the other allele", {
  # Hand-orient: flip the negative-weight SNPs to their other allele (2 - d) and
  # score with |w|. The module must reproduce this exactly.
  D    <- .pgs_dosage()
  w    <- .pgs_mixed_w[.pgs_snps]
  flip <- w < 0
  Dor  <- D; Dor[, flip] <- 2 - D[, flip]
  obs  <- !is.na(D)
  Dz   <- Dor; Dz[!obs] <- 0
  oriented <- as.numeric(Dz %*% abs(w)) / (2 * as.numeric(obs %*% abs(w)))

  res <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                 weightsFile = .pgs_mixed_wfile, weightingMode = "weighted")
  r <- smry(res, "Weighted")
  expect_close(num(r$mean), mean(oriented, na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$sd),   sd(oriented,   na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$min),  min(oriented,  na.rm = TRUE), tol = 5e-4)
  expect_close(num(r$max),  max(oriented,  na.rm = TRUE), tol = 5e-4)
})

test_that("mixed-sign weights keep proportion scores inside [0, 1]", {
  for (ms in c("SNP-wise", "zero", "mean", "exclude")) {
    res <- run_pgs(data = .pgs_gappy, snpCols = .pgs_snps,
                   weightsFile = .pgs_mixed_wfile, weightingMode = "weighted",
                   missingStrategy = ms, scaleMethod = "proportion")
    st <- as_df(res$summaryTable); st <- st[st$group == "Overall", ]
    expect_gte(min(num(st$min)), -1e-9, label = paste(ms, "min"))
    expect_lte(max(num(st$max)), 1 + 1e-9, label = paste(ms, "max"))
  }
})

test_that("nobody is dropped for having only negative-weight SNPs observed", {
  # The old denominator was 2*sum(pmax(w,0)) over observed SNPs, so an individual
  # whose every typed SNP carried a negative weight got 0 -> NA and vanished.
  d <- .test_data
  # leave only the two negative-weight SNPs typed for the first 300 people
  pos <- names(.pgs_mixed_w)[.pgs_mixed_w > 0]
  for (s in pos) d[[s]][seq_len(300)] <- NA
  res <- run_pgs(data = d, snpCols = .pgs_snps, weightsFile = .pgs_mixed_wfile,
                 weightingMode = "weighted")
  n_scored <- as.integer(smry(res, "Weighted")$n)

  D <- sapply(.pgs_snps, function(s) {
    g <- as.character(d[[s]]); g[grepl("0", g)] <- NA
    p <- strsplit(g, "/", fixed = TRUE)
    vapply(p, function(x) if (length(x) != 2 || any(is.na(x))) NA_real_ else 1, 0)
  })
  expect_equal(n_scored, sum(rowSums(!is.na(D)) > 0))
})

test_that("all-positive weights are unaffected by the negative-weight handling", {
  # lo = 0 and range = 2*sum(w) when no weight is negative, so the corrected
  # score is the old raw/max_possible exactly.
  D   <- .pgs_dosage()
  w   <- abs(.pgs_weights[.pgs_snps])          # same magnitudes, all positive
  f   <- tempfile(fileext = ".tsv"); writeLines("# all positive", f)
  df  <- data.frame(rsID = .pgs_snps, effect_allele = .pgs_effect[.pgs_snps],
                    other_allele = .pgs_other[.pgs_snps], effect_weight = w[.pgs_snps],
                    chr_name = seq_along(.pgs_snps),
                    chr_position = seq_along(.pgs_snps) * 100L)
  suppressWarnings(write.table(df, f, sep = "\t", row.names = FALSE,
                               quote = FALSE, append = TRUE))
  obs <- !is.na(D); Dz <- D; Dz[!obs] <- 0
  expected <- as.numeric(Dz %*% w) / (2 * as.numeric(obs %*% w))
  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = f,
                 weightingMode = "weighted")
  expect_close(num(smry(res, "Weighted")$mean), mean(expected, na.rm = TRUE), tol = 5e-4)
})

test_that("a SNP with no weight is reported, not silently dropped", {
  # It passes QC and still counts toward the UNWEIGHTED score, but .run drops it
  # from the weighted one. Coverage used to report the full count regardless, so
  # the table said 4 SNPs were used while the weighted score had summed 3.
  f  <- tempfile(fileext = ".tsv"); writeLines("# one blank weight", f)
  w  <- .pgs_weights; w["rs10911251"] <- NA
  df <- data.frame(rsID = .pgs_snps, effect_allele = .pgs_effect[.pgs_snps],
                   other_allele = .pgs_other[.pgs_snps], effect_weight = w[.pgs_snps],
                   chr_name = seq_along(.pgs_snps),
                   chr_position = seq_along(.pgs_snps) * 100L)
  suppressWarnings(write.table(df, f, sep = "\t", row.names = FALSE,
                               quote = FALSE, append = TRUE))

  res <- run_pgs(data = .test_data, snpCols = .pgs_snps, weightsFile = f,
                 weightingMode = "both", showCoverage = TRUE, showSnpGrid = TRUE)
  cov <- as_df(res$coverageTable)
  used <- cov$value[cov$field == "SNPs used in score"]
  expect_match(used, "3 weighted")
  expect_match(used, "4 unweighted")
  grid <- as_df(res$snpGridTable)
  expect_match(grid$allele_status[grid$rsid == "rs10911251"], "no weight")

  # the weighted score must equal one computed with that SNP simply removed
  D <- .pgs_dosage(); keep <- setdiff(.pgs_snps, "rs10911251")
  wk  <- .pgs_weights[keep]
  obs <- !is.na(D[, keep]); Dz <- D[, keep]; Dz[!obs] <- 0
  lo  <- 2 * as.numeric(obs %*% pmin(wk, 0))
  rg  <- 2 * as.numeric(obs %*% abs(wk))
  expected <- (as.numeric(Dz %*% wk) - lo) / rg
  expect_close(num(smry(res, "Weighted")$mean), mean(expected, na.rm = TRUE), tol = 5e-4)

  # and with every weight present the row stays a bare count
  ok <- run_pgs(data = .test_data, snpCols = .pgs_snps,
                weightsFile = .pgs_weightsfile, weightingMode = "both",
                showCoverage = TRUE)
  cov_ok <- as_df(ok$coverageTable)
  expect_equal(cov_ok$value[cov_ok$field == "SNPs used in score"], "4")
})
