# Tab 3: LD and haplotype — LD verified against genetics::LD, haplotype
# frequencies against haplo.stats::haplo.em.

suppressMessages({ library(genetics); library(haplo.stats) })

# genetics::LD on the pairwise-complete genotype objects (matches the backend).
ld_oracle <- function(snp1, snp2) {
  s1 <- as.character(.test_data[[snp1]]); s1[grepl("0", s1)] <- NA
  s2 <- as.character(.test_data[[snp2]]); s2[grepl("0", s2)] <- NA
  m  <- !is.na(s1) & !is.na(s2)
  g1 <- genetics::genotype(s1[m], sep = "/")
  g2 <- genetics::genotype(s2[m], sep = "/")
  genetics::LD(g1, g2)
}

# ══════════════════════════════════════════════════════════════════════════════
# ldTable — r², D', D and p-value for every SNP pair
# ══════════════════════════════════════════════════════════════════════════════

test_that("ldAnalysis: every pair matches genetics::LD", {
  result <- run_snp(data = .test_data, snps = .snps4, ldAnalysis = TRUE, ldMetric = "r2")
  tbl <- as_df(result$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldTable)

  expect_equal(nrow(tbl), choose(length(.snps4), 2))

  for (i in seq_len(nrow(tbl))) {
    o <- ld_oracle(tbl$snp1[i], tbl$snp2[i])
    lab <- paste(tbl$snp1[i], tbl$snp2[i])
    expect_close(num(tbl$r2[i]),     as.numeric(o$`r`)^2, tol = 0.0015, label = paste(lab, "r2"))
    expect_close(num(tbl$Dprime[i]), as.numeric(o$`D'`),  tol = 0.0015, label = paste(lab, "D'"))
    expect_close(num(tbl$D[i]),      as.numeric(o$`D`),   tol = 0.0015, label = paste(lab, "D"))
    mp <- num(tbl$pval[i])
    if (is.na(mp)) expect_lt(o$`P-value`, 0.001)
    else expect_close(mp, o$`P-value`, tol = 0.01, label = paste(lab, "p"))
  }
})

test_that("ldMatrix: square matrix with one row/column per SNP", {
  result <- run_snp(data = .test_data, snps = .snps4, ldMatrix = TRUE, ldMetric = "r2")
  tbl <- as_df(result$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldMatrixTable)

  expect_equal(nrow(tbl), length(.snps4))
  expect_equal(ncol(tbl), length(.snps4) + 1L)   # label column + one per SNP
})

# ══════════════════════════════════════════════════════════════════════════════
# haploFreqTable — frequencies match haplo.em; rare haplotypes pooled
# ══════════════════════════════════════════════════════════════════════════════

test_that("haploFreq: frequencies match haplo.em and sum to ~1", {
  result <- run_snp(data = .test_data, snps = .snps4, response = .resp,
                    haploFreq = TRUE, haploFreqMin = 0.01)
  tbl <- as_df(result$ldHaploGroup$haploGroup$haploFreqTable)
  freq <- num(tbl$freq)

  expect_close(sum(freq), 1, tol = 0.01)
  expect_true(any(grepl("Rare", tbl$haplotype)))   # rare haplotypes pooled

  # independent EM estimate
  mat <- do.call(cbind, lapply(.snps4, function(s) {
    p  <- strsplit(as.character(.test_data[[s]]), "/", fixed = TRUE)
    a1 <- sapply(p, `[`, 1); a2 <- sapply(p, `[`, 2)
    a1[grepl("0", a1)] <- NA; a2[grepl("0", a2)] <- NA
    cbind(a1, a2)
  }))
  em <- haplo.stats::haplo.em(geno = mat, locus.label = .snps4)
  ref <- setNames(em$hap.prob, apply(em$haplotype, 1, paste, collapse = "-"))

  common <- tbl[!grepl("Rare", tbl$haplotype), ]
  for (i in seq_len(nrow(common))) {
    h <- common$haplotype[i]
    expect_true(h %in% names(ref), label = paste("haplotype", h, "present in EM"))
    expect_close(num(common$freq[i]), ref[[h]], tol = 0.005,
                 label = paste("freq", h))
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# haploAssocTable / haploInteraction — structure (haplo.glm based)
# ══════════════════════════════════════════════════════════════════════════════

test_that("haploAssoc: reference haplotype has effect 1 and others are estimated", {
  result <- run_snp(data = .test_data, snps = .snps4, response = .resp,
                    covariates = .covars, haploAssoc = TRUE)
  tbl <- as_df(result$ldHaploGroup$haploGroup$haploAssocTable)

  expect_gt(nrow(tbl), 1L)
  expect_true(all(c("haplotype", "effect", "pval") %in% names(tbl)))
  effects <- num(tbl$effect)
  expect_true(any(abs(effects - 1) < 1e-6, na.rm = TRUE))   # reference OR = 1
  expect_true(all(effects > 0, na.rm = TRUE))               # ORs are positive
})

test_that("haploInteraction: conditional tables populate without error", {
  result <- run_snp(data = .test_data, snps = .snps4, response = .resp,
                    covariates = .covars, haploInteraction = TRUE)
  expect_gt(nrow(as_df(result$ldHaploGroup$haploGroup$haploInteractionTable)), 0L)
  expect_gt(nrow(as_df(result$ldHaploGroup$haploGroup$haploCondCovarTable)), 0L)
  expect_gt(nrow(as_df(result$ldHaploGroup$haploGroup$haploCondHaploTable)), 0L)
})

# The interaction LRT df is the number of geno x covariate interaction
# coefficients — NOT the haplo.glm df.residual difference (the EM expands each
# subject into weighted haplotype-pair rows, so that difference is meaningless).
test_that("haploInteraction LRT p-value uses the interaction-term df", {
  res <- run_snp(data = .test_data, snps = .snps2, response = .resp,
                 covariates = .covars, haploInteraction = TRUE, haploFreqMin = 0.05)
  notes <- res$ldHaploGroup$haploGroup$haploInteractionTable$notes
  note  <- Filter(function(n) n$key == "lrt_inter", notes)[[1]]$note
  mod_p <- num(sub(".*: *", "", note))

  na.geno.keep <- getFromNamespace("na.geno.keep", "SNPstats")
  mat <- do.call(cbind, lapply(.snps2, function(s) {
    p <- strsplit(as.character(.test_data[[s]]), "/", fixed = TRUE)
    a1 <- sapply(p, `[`, 1); a2 <- sapply(p, `[`, 2)
    a1[grepl("0", a1)] <- NA; a2[grepl("0", a2)] <- NA; cbind(a1, a2)
  }))
  geno <- haplo.stats::setupGeno(mat, locus.label = .snps2)
  md <- data.frame(y = as.numeric(as.factor(.test_data[[.resp]])) - 1L,
                   sex = .test_data$sex, age = .test_data$age, bmiOMS = .test_data$bmiOMS)
  md$geno <- geno
  ctl <- haplo.stats::haplo.glm.control(haplo.effect = "additive", haplo.freq.min = 0.05)
  set.seed(20240920L)
  fm <- haplo.stats::haplo.glm(y ~ geno * sex + age + bmiOMS, family = binomial,
                               data = md, na.action = na.geno.keep, control = ctl)
  set.seed(20240920L)
  fa <- haplo.stats::haplo.glm(y ~ geno + sex + age + bmiOMS, family = binomial,
                               data = md, na.action = na.geno.keep, control = ctl)
  inter_df <- sum(grepl(":", rownames(summary(fm)$coefficients)) &
                  grepl("geno", rownames(summary(fm)$coefficients)))
  ora_p <- pchisq(fa$deviance - fm$deviance, df = inter_df, lower.tail = FALSE)

  if (ora_p < 0.001) expect_match(note, "<")
  else               expect_close(mod_p, ora_p, tol = 0.005, label = "interaction LRT p")
})

# haplo.em / haplo.glm handle partial missingness via EM, so the default keeps
# every subject typed at >= 1 SNP; completeCases restricts to all-SNP-typed rows.
# With informative missingness (SNP2 blanked for every SNP1 C-carrier), the
# default EM still recovers the C-* haplotypes from SNP1, while completeCases
# drops those subjects and loses them — so the two must differ markedly.
test_that("completeCases restricts the haplotype EM to complete cases", {
  dat <- .test_data
  dat$rs10911251[grepl("C", as.character(dat$rs12080929))] <- NA
  freq_of <- function(cc) {
    tbl <- as_df(run_snp(data = dat, snps = .snps2, response = .resp,
                         haploFreq = TRUE, haploFreqMin = 0.05, completeCases = cc)$
                 ldHaploGroup$haploGroup$haploFreqTable)
    setNames(num(tbl$freq), tbl$haplotype)
  }
  f_default <- freq_of(FALSE)
  f_cc      <- freq_of(TRUE)
  expect_true(any(grepl("^C", names(f_default))))   # default recovers C-* haplotypes
  expect_false(any(grepl("^C", names(f_cc))))        # complete-case drops all C-carriers
})

# ══════════════════════════════════════════════════════════════════════════════
# LD metric option (r2 / D' / D) drives the matrix
# ══════════════════════════════════════════════════════════════════════════════

test_that("ldMatrix diagonal reflects the chosen metric", {
  for (metric in c("r2", "Dprime")) {
    res <- run_snp(data = .test_data, snps = .snps4, ldMatrix = TRUE, ldMetric = metric)
    tbl <- as_df(res$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldMatrixTable)
    expect_equal(nrow(tbl), length(.snps4))
    expect_equal(ncol(tbl), length(.snps4) + 1L)
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# Regression: haploAssoc must work WITHOUT covariates (bug fix) and is
# reproducible run-to-run (module seeds haplo.glm's EM internally).
# ══════════════════════════════════════════════════════════════════════════════

test_that("haploAssoc returns rows without covariates (regression)", {
  res <- run_snp(data = .test_data, snps = .snps2, response = .resp,
                 haploAssoc = TRUE, haploFreqMin = 0.05)
  tbl <- as_df(res$ldHaploGroup$haploGroup$haploAssocTable)
  expect_gt(nrow(tbl), 1L)
  expect_true(any(abs(num(tbl$effect) - 1) < 1e-6, na.rm = TRUE))   # reference OR = 1
})

test_that("haploAssoc is reproducible without external seeding", {
  f <- function()
    as_df(run_snp(data = .test_data, snps = .snps4, response = .resp,
          haploAssoc = TRUE, haploFreqMin = 0.05)$ldHaploGroup$haploGroup$haploAssocTable)$effect
  expect_equal(num(f()), num(f()))
})

test_that("haploAssoc does not disturb the caller's RNG stream", {
  set.seed(7); a1 <- runif(1)
  invisible(run_snp(data = .test_data, snps = .snps2, response = .resp,
                    haploAssoc = TRUE, haploFreqMin = 0.05))
  a2 <- runif(1)                                   # stream continues if undisturbed
  set.seed(7); b <- runif(2)
  expect_equal(c(a1, a2), b)
})

# Regression: a quantitative response used to leave the overall-association LRT
# footnote empty ("P = ") because the lm null model has no $df.null. The df now
# comes from the geno coefficient count and the gaussian deviance is divided by
# the model dispersion (chi-square), mirroring haplo.glm's own deviance analysis.
test_that("haploAssoc quantitative reports a scaled-deviance LRT p (not empty)", {
  res  <- run_snp(data = .test_data, snps = .snps4, response = "age",
                  responseType = "quantitative", haploAssoc = TRUE, haploFreqMin = 0.05)
  htbl <- res$ldHaploGroup$haploGroup$haploAssocTable
  tbl  <- as_df(htbl)

  expect_equal(htbl$getColumn("effect")$title, intToUtf8(0x3B2))  # beta, not OR
  base <- tbl[grepl("Ref", tbl$haplotype), ]
  expect_close(num(base$effect), 0, tol = 1e-6)                    # reference beta = 0

  note  <- Filter(function(n) n$key == "lrt_assoc", htbl$notes)[[1]]$note
  mod_p <- num(sub(".*P = *", "", note))
  expect_false(is.na(mod_p))                                       # was empty "P = "

  na.geno.keep <- getFromNamespace("na.geno.keep", "SNPstats")
  mat <- do.call(cbind, lapply(.snps4, function(s) {
    p <- strsplit(as.character(.test_data[[s]]), "/", fixed = TRUE)
    a1 <- sapply(p, `[`, 1); a2 <- sapply(p, `[`, 2)
    a1[grepl("0", a1)] <- NA; a2[grepl("0", a2)] <- NA; cbind(a1, a2)
  }))
  geno <- haplo.stats::setupGeno(mat, locus.label = .snps4)
  md <- data.frame(y = .test_data$age); md$geno <- geno
  set.seed(20240920L)
  hf <- haplo.stats::haplo.glm(y ~ geno, family = gaussian, data = md,
          na.action = na.geno.keep,
          control = haplo.stats::haplo.glm.control(haplo.freq.min = 0.05))
  df_geno <- sum(grepl("^geno", names(coef(hf))))
  ora_p   <- pchisq((hf$null.deviance - hf$deviance) / summary(hf)$dispersion,
                    df = df_geno, lower.tail = FALSE)
  expect_close(mod_p, ora_p, tol = 0.005, label = "quant haplo LRT p")
})

# Categorical responses are intentionally not supported for haplotype
# association/interaction: the module hides the tables and shows a notice.
test_that("haploAssoc is blocked (not implemented) for a categorical response", {
  res <- run_snp(data = .test_data, snps = .snps2, response = "bmiOMS",
                 responseType = "categorical", haploAssoc = TRUE)
  hg  <- res$ldHaploGroup$haploGroup
  expect_true(hg$haploNotImplMsg$visible)
  expect_match(hg$haploNotImplMsg$content, "only implemented for binary and quantitative")
  expect_equal(nrow(as_df(hg$haploAssocTable)), 0L)
})
