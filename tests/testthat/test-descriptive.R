# Tab 1: Descriptive results — verified against independent allele counting and
# genetics::HWE.exact.

# genetics is a Suggests-only oracle now: the module no longer uses it at
# runtime (see R/snp_genetics.R), but cross-checking against it is still the
# strongest available check that the reimplementation stayed faithful. Only the
# tests that actually call it skip when it is absent — the rest of this file
# tests module internals and must keep running either way.

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
  skip_if_not_installed("genetics")
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
  skip_if_not_installed("genetics")
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
  skip_if_not_installed("genetics")
  result <- run_snp(data = .test_data, snps = .snps2, genoFreq = TRUE)

  for (snp in .snps2) {
    o   <- desc_oracle(.test_data[[snp]])
    tab <- as_df(result$descGroup$descSnpResults$get(key = snp)$genoFreqTable)
    counts <- as.integer(sub("\\s*\\(.*$", "", tab$stat))
    expect_equal(sum(counts), o$n, label = paste(snp, "geno counts sum to N"))
  }
})

test_that("hweTest: per-SNP HWE p-value matches HWE.exact", {
  skip_if_not_installed("genetics")
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

test_that("missingness plot renders from image state, not a private field", {
  # jamovi rebuilds the analysis object on every option click, so a private R6
  # field does not survive; if the image is redrawn without a preceding .run()
  # the plot used to return FALSE and the user got a blank panel.
  opts <- SNPstats::snpStatsOptions$new(snps = .snps4, response = .resp,
                                        showMissingnessPlot = TRUE,
                                        missingnessThreshold = 0)
  a <- SNPstats::snpStatsClass$new(options = opts, data = .test_data)
  a$init(); a$run()
  img <- a$results$descGroup$missingnessPlot
  expect_false(is.null(img$state))
  expect_equal(length(img$state), length(.snps4))

  render <- function(analysis, image) {
    png(tempfile(fileext = ".png"), width = 560, height = 500)
    on.exit(dev.off(), add = TRUE)
    analysis$.__enclos_env__$private$.plotMissingness(image)
  }
  expect_true(render(a, img))

  # A fresh object that never ran: state carries the plot, a field could not.
  b <- SNPstats::snpStatsClass$new(options = opts, data = .test_data)
  b$init()
  b$results$descGroup$missingnessPlot$setState(img$state)
  expect_true(render(b, b$results$descGroup$missingnessPlot))
})

# ══════════════════════════════════════════════════════════════════════════════
# snp_genetics.R — the in-house replacements for the `genetics` package
#
# The module no longer depends on `genetics` (obsolete upstream). These pin the
# reimplementation against it directly, over every SNP in the shipped dataset,
# so a divergence shows up here rather than as a shifted result three tables
# away. genetics is Suggests-only; the whole block skips without it.
# ══════════════════════════════════════════════════════════════════════════════

test_that("snp_genotype reproduces genetics::genotype tables on every SNP", {
  skip_if_not_installed("genetics")
  G     <- asNamespace("SNPstats")
  pg    <- get("parse_genotype", envir = G)
  clean <- get("clean_null_alleles", envir = G)

  snps <- names(.test_data)[grepl("^rs", names(.test_data))]
  expect_gt(length(snps), 50)

  for (s in snps) {
    x  <- clean(as.character(.test_data[[s]]))
    x  <- x[!is.na(x)]
    gr <- genetics::genotype(x, sep = "/")
    gn <- pg(x, NULL)
    sr <- summary(gr); sn <- summary(gn)

    expect_identical(rownames(sn$allele.freq), rownames(sr$allele.freq), info = s)
    expect_equal(unname(sn$allele.freq[, "Count"]),
                 unname(sr$allele.freq[, "Count"]), info = s)
    expect_equal(sn$n.typed, sr$n.typed, info = s)

    if (is.matrix(sr$genotype.freq)) {
      expect_identical(rownames(sn$genotype.freq), rownames(sr$genotype.freq), info = s)
      expect_equal(unname(sn$genotype.freq[, "Count"]),
                   unname(sr$genotype.freq[, "Count"]), info = s)
    }
    # The allele matrix is what feeds haplo.stats::setupGeno. Compare the
    # values only: genetics::allele() tags the result with `which` /
    # `allele.names` attributes that cbind() drops in the backend anyway.
    am_r <- genetics::allele(gr); am_n <- get("snp_allele", envir = G)(gn)
    expect_identical(as.vector(am_n[, 1]), as.vector(am_r[, 1]), info = s)
    expect_identical(as.vector(am_n[, 2]), as.vector(am_r[, 2]), info = s)
  }
})

test_that("snp_hwe_exact matches genetics::HWE.exact", {
  skip_if_not_installed("genetics")
  G  <- asNamespace("SNPstats")
  hx <- get("snp_hwe_exact", envir = G)
  hp <- get("hwe_exact_p",   envir = G)
  gt <- function(a, b, c) genetics::genotype(rep(c("A/A","A/B","B/B"), c(a,b,c)), sep = "/")

  set.seed(9)
  cases <- list(c(95,92,13), c(10,2,8), c(1,1,1), c(0,10,0), c(1,20,1), c(100,0,100))
  for (i in 1:25) cases[[length(cases)+1]] <- as.integer(rmultinom(1, sample(5:500,1), runif(3)))

  for (cc in cases) {
    if (2*cc[1]+cc[2] == 0 || 2*cc[3]+cc[2] == 0) next     # monomorphic: both error
    expect_equal(hp(cc[1], cc[2], cc[3]),
                 genetics::HWE.exact(gt(cc[1],cc[2],cc[3]))$p.value,
                 tolerance = 1e-10, info = paste(cc, collapse = "/"))
  }

  # A monomorphic locus must ERROR, not return NA: callers tryCatch it and treat
  # the failure as "no HWE result", which is the correct output there.
  expect_error(hx(get("parse_genotype", envir = G)(rep("A/A", 20), NULL)),
               "2 alleles")
})

test_that("snp_ld reaches the haplotype-frequency MLE", {
  skip_if_not_installed("genetics")
  G     <- asNamespace("SNPstats")
  pg    <- get("parse_genotype", envir = G)
  clean <- get("clean_null_alleles", envir = G)
  sld   <- get("snp_ld", envir = G)

  snps <- names(.test_data)[grepl("^rs", names(.test_data))][1:8]
  npair <- 0
  for (i in 1:(length(snps)-1)) for (j in (i+1):length(snps)) {
    v1 <- clean(as.character(.test_data[[snps[i]]]))
    v2 <- clean(as.character(.test_data[[snps[j]]]))
    m  <- !is.na(v1) & !is.na(v2)
    gr1 <- tryCatch(genetics::genotype(v1[m], sep="/"), error = function(e) NULL)
    gr2 <- tryCatch(genetics::genotype(v2[m], sep="/"), error = function(e) NULL)
    if (is.null(gr1) || is.null(gr2)) next
    ref <- tryCatch(genetics::LD(gr1, gr2), error = function(e) NULL)
    if (is.null(ref) || is.na(ref$D)) next
    mine <- sld(pg(v1[m], NULL), pg(v2[m], NULL))
    npair <- npair + 1

    # Same D sign and the same value to well within the printed 3 decimals.
    expect_equal(sign(mine$D), sign(ref$D), info = paste(snps[i], snps[j]))
    expect_equal(mine$D,       ref$D,       tolerance = 1e-3)
    expect_equal(mine[["R^2"]], ref$r^2,    tolerance = 1e-3)
    expect_equal(mine[["D'"]], ref[["D'"]], tolerance = 5e-3)

    # snp_ld solves the EM to convergence, so its D must be at least as likely
    # as genetics' under the two-locus multinomial. (genetics stops early; on
    # this dataset it lands up to ~1.4e-5 short of the MLE.)
    a <- strsplit(v1[m], "/"); b <- strsplit(v2[m], "/")
    A <- names(sort(table(unlist(a)), decreasing = TRUE))[1]
    B <- names(sort(table(unlist(b)), decreasing = TRUE))[1]
    nA <- vapply(a, function(z) sum(z == A), 0L); nB <- vapply(b, function(z) sum(z == B), 0L)
    tb <- table(factor(nA, 0:2), factor(nB, 0:2))
    pA <- sum(nA)/(2*length(nA)); pB <- sum(nB)/(2*length(nB))
    ll <- function(D) {
      p <- c(pA*pB + D, pA*(1-pB) - D, (1-pA)*pB - D, (1-pA)*(1-pB) + D)
      if (min(p) <= 0) return(-Inf)
      P <- matrix(0, 3, 3)
      P[3,3] <- p[1]^2;      P[3,2] <- 2*p[1]*p[2];        P[3,1] <- p[2]^2
      P[2,3] <- 2*p[1]*p[3]; P[2,2] <- 2*(p[1]*p[4]+p[2]*p[3]); P[2,1] <- 2*p[2]*p[4]
      P[1,3] <- p[3]^2;      P[1,2] <- 2*p[3]*p[4];        P[1,1] <- p[4]^2
      sum(tb * log(P))
    }
    expect_gte(ll(mine$D), ll(ref$D) - 1e-9)
  }
  expect_gt(npair, 20)
})
