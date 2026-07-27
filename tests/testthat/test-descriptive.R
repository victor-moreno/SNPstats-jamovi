# Tab 1: Descriptive results — verified against independent allele counting and
# an independently-derived exact HWE test (see helper-data.R).

# Oracle: count alleles / genotypes / HWE directly from the raw "A/B" column.
# No package involved — `geno_counts_oracle` counts off the strings and
# `hwe_closed_oracle` writes the exact-test probability a different way from the
# implementation.
desc_oracle <- function(col) {
  o <- geno_counts_oracle(col)
  props <- o$allele_cnt / sum(o$allele_cnt)
  cnt3  <- {
    a1 <- o$alleles[1]; a2 <- if (length(o$alleles) > 1) o$alleles[2] else a1
    g  <- o$geno_cnt
    gv <- function(nm) if (nm %in% names(g)) as.integer(g[[nm]]) else 0L
    c(gv(paste0(a1, "/", a1)), gv(paste0(a1, "/", a2)), gv(paste0(a2, "/", a2)))
  }
  list(
    n       = o$n_typed,
    missing = o$n_missing,
    maf     = unname(min(props)),
    counts  = cnt3,
    hwe     = hwe_closed_oracle(cnt3[1], cnt3[2], cnt3[3]))
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
# snp_genetics.R — the in-house genotype / HWE / LD code
#
# The module used to get these from the `genetics` package, which is obsolete
# upstream and is no longer installed at all. The oracles here are independent
# of the implementation by construction (helper-data.R): direct string counting,
# a brute-force enumeration of the allele pairings, an alternative closed form,
# and a numerical likelihood maximisation.
# ══════════════════════════════════════════════════════════════════════════════

test_that("snp_genotype tables match direct counting on every SNP", {
  G     <- asNamespace("SNPstats")
  pg    <- get("parse_genotype", envir = G)
  clean <- get("clean_null_alleles", envir = G)

  snps <- names(.test_data)[grepl("^rs", names(.test_data))]
  expect_gt(length(snps), 50)

  for (s in snps) {
    x  <- clean(as.character(.test_data[[s]]))
    x  <- x[!is.na(x)]
    gn <- pg(x, NULL)
    sn <- summary(gn)
    o  <- geno_counts_oracle(.test_data[[s]])

    # Allele order is a module convention (descending count, ties by name
    # descending); the oracle reproduces it from the raw strings.
    expect_identical(rownames(sn$allele.freq), o$alleles, info = s)
    expect_equal(unname(sn$allele.freq[, "Count"]), o$allele_cnt, info = s)
    expect_equal(sn$n.typed, o$n_typed, info = s)
    expect_equal(sum(sn$allele.freq[, "Count"]), 2L * o$n_typed, info = s)

    # Genotype rows: same labels, same counts, summing to N.
    gf <- sn$genotype.freq
    expect_setequal(rownames(gf), names(o$geno_cnt))
    expect_equal(unname(gf[, "Count"]),
                 as.integer(o$geno_cnt[rownames(gf)]), info = s)
    expect_equal(sum(gf[, "Count"]), o$n_typed, info = s)

    # The allele matrix that feeds haplo.stats::setupGeno must re-split into the
    # same genotypes it came from.
    am <- get("snp_allele", envir = G)(gn)
    expect_equal(paste0(am[, 1], "/", am[, 2]), as.character(gn), info = s)
  }
})

test_that("snp_hwe_exact matches an independent enumeration of the allele pairings", {
  G  <- asNamespace("SNPstats")
  hp <- get("hwe_exact_p", envir = G)

  # Brute force: enumerate every pairing of the 2n alleles into n genotypes.
  # This is the combinatorial definition of the test with no formula at all, so
  # it shares nothing with the implementation. Feasible only for tiny n.
  for (cc in list(c(2,2,1), c(1,2,2), c(3,1,1), c(0,4,0), c(1,3,1), c(2,1,3),
                  c(4,1,1), c(1,1,4), c(2,2,2))) {
    expect_equal(hp(cc[1], cc[2], cc[3]),
                 hwe_bruteforce(cc[1], cc[2], cc[3]),
                 tolerance = 1e-9, info = paste(cc, collapse = "/"))
  }

  # Larger samples against the alternative closed form.
  set.seed(9)
  cases <- list(c(95,92,13), c(10,2,8), c(1,20,1), c(100,0,100), c(1502,1120,205))
  for (i in 1:25) cases[[length(cases)+1]] <- as.integer(rmultinom(1, sample(5:500,1), runif(3)))
  for (cc in cases) {
    if (2*cc[1]+cc[2] == 0 || 2*cc[3]+cc[2] == 0) next     # monomorphic
    expect_equal(hp(cc[1], cc[2], cc[3]),
                 hwe_closed_oracle(cc[1], cc[2], cc[3]),
                 tolerance = 1e-10, info = paste(cc, collapse = "/"))
  }

  # Large, near-equilibrium samples: the exact test must agree with the
  # asymptotic chi-square it converges to. Independent of both formulations.
  set.seed(3)
  for (i in 1:10) {
    n <- 4000; p <- runif(1, 0.2, 0.8)
    g <- rmultinom(1, n, c(p^2, 2*p*(1-p), (1-p)^2))
    n11 <- g[1]; n12 <- g[2]; n22 <- g[3]
    nn  <- n11+n12+n22; pa <- (2*n11+n12)/(2*nn)
    e   <- c(pa^2, 2*pa*(1-pa), (1-pa)^2) * nn
    x2  <- sum((c(n11,n12,n22) - e)^2 / e)
    expect_equal(hp(n11, n12, n22), pchisq(x2, 1, lower.tail = FALSE),
                 tolerance = 0.05)
  }

  # A monomorphic locus must ERROR, not return NA: callers tryCatch it and treat
  # the failure as "no HWE result", which is the correct output there.
  expect_error(get("snp_hwe_exact", envir = G)(
                 get("parse_genotype", envir = G)(rep("A/A", 20), NULL)),
               "2 alleles")
})

test_that("snp_ld reaches the haplotype-frequency MLE", {
  G   <- asNamespace("SNPstats")
  pg  <- get("parse_genotype", envir = G)
  cl  <- get("clean_null_alleles", envir = G)
  sld <- get("snp_ld", envir = G)

  snps <- names(.test_data)[grepl("^rs", names(.test_data))][1:8]
  npair <- 0
  for (i in 1:(length(snps)-1)) for (j in (i+1):length(snps)) {
    v1 <- cl(as.character(.test_data[[snps[i]]]))
    v2 <- cl(as.character(.test_data[[snps[j]]]))
    m  <- !is.na(v1) & !is.na(v2)
    mine <- sld(pg(v1[m], NULL), pg(v2[m], NULL))
    if (is.null(mine) || is.na(mine$D)) next
    ref <- ld_oracle_mle(.test_data[[snps[i]]], .test_data[[snps[j]]])
    npair <- npair + 1
    lab <- paste(snps[i], snps[j])

    # The oracle maximises the same likelihood by 1-D numerical search instead
    # of by EM, so agreement validates the EM's answer, not just its agreement
    # with another EM.
    expect_equal(mine$n, ref$n, info = lab)
    expect_equal(mine$D, ref$D, tolerance = 1e-6, info = lab)
    expect_equal(mine[["D'"]], ref$Dprime, tolerance = 1e-5, info = lab)
    expect_equal(mine[["R^2"]], ref$R2, tolerance = 1e-6, info = lab)
    # Looser on P: optimize() pins D to ~1e-8, and the chi-square amplifies that
    # residual. Still four orders tighter than the 3 decimals ever displayed.
    expect_equal(mine[["P-value"]], ref$p, tolerance = 1e-4, info = lab)
    # and it must actually be the maximum, not merely close to it
    expect_gte(ref$loglik(mine$D), ref$loglik(ref$D) - 1e-9)
  }
  expect_gt(npair, 20)
})
