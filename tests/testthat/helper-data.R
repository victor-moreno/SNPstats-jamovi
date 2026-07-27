# Shared fixtures and verification oracles — loaded automatically by testthat.
# Run the suite with tests/run_tests.R (see that file for the library setup).
#
# The oracles below re-derive every quantity the module reports using only base R
# (glm / lm / nnet::multinom) and haplo.stats, independently of the module's
# internal functions, so the comparison is a genuine cross-check.

.test_data <- local({
  candidates <- c(
    file.path("data", "CRCgenet-SNPs.tsv"),
    file.path("..", "..", "data", "CRCgenet-SNPs.tsv")
  )
  path <- Filter(file.exists, candidates)
  if (length(path) == 0L)
    stop("Test data not found — run tests from the package root.")
  read.delim(path[[1L]], header = TRUE, stringsAsFactors = TRUE)
})

# ── Variable aliases used across test files ──────────────────────────────────
.snps2  <- c("rs12080929", "rs10911251")
.snps4  <- c("rs12080929", "rs10911251", "rs72647484", "rs6691170")
.covars <- c("sex", "age", "bmiOMS")
.resp   <- "phenotype"

# ── Calling convention ───────────────────────────────────────────────────────
# snpStats() uses tidy-eval (jmvcore::resolveQuo) on `snps`, `response` and
# `covariates`: passing a *variable* deparses to its name and fails. do.call()
# inlines the values, so it accepts variables holding character vectors.
run_snp <- function(...) do.call(SNPstats::snpStats, list(...))

# snpPGS has no weights-file *path* option — the weights travel as base64
# content so that a saved .omv never re-reads a file from whoever opens it.
# `weightsFile = <path>` is a test-only convenience that expands to the real
# weightsContent / weightsFilename pair via the public pgs_weights() helper,
# which it therefore also exercises.
run_pgs <- function(...) {
  args <- list(...)
  # [[ ]], not $: `$` partial-matches, so args$weightsFile would pick up a
  # weightsFilename passed by the embedded-content tests.
  if (!is.null(args[["weightsFile"]])) {
    args <- c(args[names(args) != "weightsFile"],
              SNPstats::pgs_weights(args[["weightsFile"]]))
  }
  do.call(SNPstats::snpPGS, args)
}

# jmvcore >= 2 exposes Table$asDF as an active binding (a property, not a
# method). Tolerate both forms.
as_df <- function(tbl) {
  df <- try(tbl$asDF, silent = TRUE)
  if (inherits(df, "data.frame")) return(df)
  tbl$asDF()
}

# Module result cells are pre-formatted strings (3 decimals; "< 0.001" for tiny
# p-values). Parse to numeric; "< 0.001" becomes NA and is handled per-test.
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ── Independent genotype encoding ────────────────────────────────────────────
# Minor-allele dosage (0/1/2) parsed directly from the raw "A/B" strings.
# Null-allele genotypes (containing "0") and malformed cells become NA.
dose_minor <- function(col) {
  s <- as.character(col); s[grepl("0", s)] <- NA
  parts <- strsplit(s, "/", fixed = TRUE)
  minor <- names(sort(table(unlist(parts))))[1]
  vapply(parts, function(p)
    if (length(p) != 2 || any(is.na(p))) NA_integer_ else sum(p == minor),
    integer(1))
}

# Binary 0/1 coding matching the module: as.integer(as.factor(x)) - 1
# (first factor level -> 0). For phenotype: Case = 0, Control = 1.
bin01 <- function(x) { r <- as.integer(as.factor(x)) - 1L; r[is.na(x)] <- NA_integer_; r }

# Encode the genetic model from minor-allele dosage, matching encode_model().
encode_dose <- function(g, model) switch(model,
  codominant   = factor(g, levels = c(0L, 1L, 2L)),
  dominant     = as.integer(g >= 1),
  recessive    = as.integer(g == 2),
  overdominant = as.integer(g == 1),
  logadditive  = as.numeric(g))

# ── Association oracle ───────────────────────────────────────────────────────
# Returns a data.frame of the non-reference coefficient rows (in the same order
# the module reports them) plus global_p / aic / bic, replicating fit_model().
assoc_oracle <- function(snpcol, resp_raw, type, covdf = NULL,
                         model = "logadditive", ci = 95) {
  g   <- dose_minor(snpcol)
  enc <- encode_dose(g, model)
  resp <- switch(type,
    binary       = bin01(resp_raw),
    categorical  = as.factor(resp_raw),
    as.numeric(resp_raw))
  df <- data.frame(resp = resp, snp = enc)
  if (!is.null(covdf)) df <- cbind(df, covdf)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  z  <- stats::qnorm(1 - (1 - ci / 100) / 2)

  if (type == "categorical") {
    full <- nnet::multinom(resp ~ ., data = df, trace = FALSE)
    null <- nnet::multinom(resp ~ . - snp, data = df, trace = FALSE)
    co <- summary(full)$coefficients; se <- summary(full)$standard.errors
    cols <- grep("^snp", colnames(co), value = TRUE)
    out <- list()
    for (cat in rownames(co)) for (cc in cols) {
      b <- co[cat, cc]; s <- se[cat, cc]
      out[[length(out) + 1L]] <- data.frame(
        category = cat, term = cc,
        effect = exp(b), ciLow = exp(b - z * s), ciHigh = exp(b + z * s),
        pval = 2 * (1 - stats::pnorm(abs(b / s))), stringsAsFactors = FALSE)
    }
    res <- do.call(rbind, out)
    attr(res, "global_p") <- anova(null, full)[2, "Pr(Chi)"]
    attr(res, "aic") <- AIC(full)
    return(res)
  }

  if (type == "binary") {
    full <- glm(resp ~ ., data = df, family = binomial())
    null <- glm(resp ~ . - snp, data = df, family = binomial())
    pcol <- "Pr(>|z|)"; gp <- anova(null, full, test = "Chisq")[2, "Pr(>Chi)"]
    expo <- exp
  } else {
    full <- lm(resp ~ ., data = df)
    null <- lm(resp ~ . - snp, data = df)
    pcol <- "Pr(>|t|)"; gp <- anova(null, full, test = "F")[2, "Pr(>F)"]
    expo <- identity
  }
  co  <- summary(full)$coefficients
  rows <- grep("^snp", rownames(co))
  cis <- suppressMessages(suppressWarnings(
    stats::confint(full, level = ci / 100)[rows, , drop = FALSE]))
  res <- data.frame(
    term   = rownames(co)[rows],
    effect = expo(co[rows, "Estimate"]),
    ciLow  = expo(cis[, 1]), ciHigh = expo(cis[, 2]),
    pval   = co[rows, pcol], stringsAsFactors = FALSE)
  attr(res, "global_p") <- gp
  attr(res, "aic") <- AIC(full)
  attr(res, "bic") <- BIC(full)
  res
}

# Pull the non-reference coefficient rows from a module assocTable for one model
# (rows that carry a confidence interval, i.e. ciLow is non-empty).
assoc_coef_rows <- function(tbl) {
  df <- as_df(tbl)
  df[!is.na(num(df$ciLow)), , drop = FALSE]
}

# Absolute-tolerance comparison robust to the module's 3-decimal string output.
expect_close <- function(actual, expected, tol = 0.0015, label = NULL) {
  testthat::expect_lt(abs(actual - expected), tol,
                      label = label %||% deparse(substitute(actual)))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ══════════════════════════════════════════════════════════════════════════════
# Independent oracles for the in-house genotype / HWE / LD code
#
# The module used to lean on the `genetics` package for these, and the tests
# cross-checked against it. `genetics` is gone (obsolete upstream), so these
# oracles have to stand on their own. They deliberately do NOT mirror the
# implementation in R/snp_genetics.R:
#
#   geno_counts_oracle  counts straight off the raw "A/B" strings
#   hwe_bruteforce      enumerates every pairing of the allele pool - the
#                       combinatorial definition of the test, no formula at all
#   hwe_closed_oracle   the same probability written a different way
#                       (multinomial / C(2n,n1)) from the implementation's
#   ld_oracle_mle       maximises the two-locus likelihood by 1-D numerical
#                       search instead of by EM
#
# The last one is the strongest of the four: it answers the same estimation
# question by a different numerical method, so it validates the EM's *answer*
# rather than its agreement with some other implementation.
# ══════════════════════════════════════════════════════════════════════════════

# Allele and genotype counts read directly off the raw column.
geno_counts_oracle <- function(col) {
  s     <- as.character(col); s[grepl("0", s)] <- NA
  typed <- s[!is.na(s)]
  parts <- strsplit(typed, "/", fixed = TRUE)
  al    <- unlist(parts)
  acnt  <- table(al)
  # module convention: descending count, ties broken by allele name descending
  nms   <- sort(names(acnt), decreasing = TRUE)
  order_al <- nms[order(-as.integer(acnt[nms]))]
  canon <- vapply(parts, function(p) {
    i <- match(p, order_al)
    paste0(order_al[min(i)], "/", order_al[max(i)])
  }, character(1))
  list(alleles     = order_al,
       allele_cnt  = as.integer(acnt[order_al]),
       geno_cnt    = table(canon),
       n_typed     = length(typed),
       n_missing   = sum(is.na(s)))
}

# Exact HWE p-value by brute force: enumerate every way of pairing the 2n
# alleles into n genotypes and count the resulting heterozygote counts. This is
# the definition the exact test formalises, with no probability formula used.
# Only feasible for tiny n - the number of pairings grows as (2n-1)!!.
hwe_bruteforce <- function(n11, n12, n22) {
  pool <- c(rep("A", 2*n11 + n12), rep("B", 2*n22 + n12))
  n    <- (n11 + n12 + n22)
  stopifnot(n <= 6)
  counts <- new.env(parent = emptyenv())
  pair_up <- function(v, nhet) {
    if (length(v) == 0) {
      k <- as.character(nhet)
      assign(k, (if (exists(k, counts)) get(k, counts) else 0) + 1, counts)
      return(invisible(NULL))
    }
    first <- v[1]
    for (j in 2:length(v)) {
      rest <- v[-c(1, j)]
      pair_up(rest, nhet + (first != v[j]))
    }
  }
  pair_up(pool, 0)
  ks <- as.integer(ls(counts))
  w  <- vapply(ls(counts), function(k) get(k, counts), 0)
  p  <- w / sum(w)
  obs <- p[as.character(n12)]
  sum(p[p <= obs * (1 + 1e-9)])
}

# Same exact test, written as 2^k * multinomial(n; n11,k,n22) / C(2n, n1) -
# algebraically identical to the implementation but a different code path, so a
# transcription slip in either shows up.
hwe_closed_oracle <- function(n11, n12, n22) {
  n  <- n11 + n12 + n22
  n1 <- 2*n11 + n12; n2 <- 2*n22 + n12
  if (n1 == 0 || n2 == 0) return(NA_real_)
  ks <- seq(n1 %% 2, min(n1, n2), by = 2)
  lp <- vapply(ks, function(k) {
    a <- (n1 - k)/2; b <- (n2 - k)/2
    lfactorial(n) - lfactorial(a) - lfactorial(k) - lfactorial(b) +
      k * log(2) - lchoose(2*n, n1)
  }, 0)
  p   <- exp(lp)
  obs <- p[ks == n12]
  sum(p[p <= obs * (1 + 1e-9)])
}

# Two-locus LD by direct numerical maximisation of the multinomial likelihood
# over D, with the allele frequencies fixed at their observed values. No EM.
ld_oracle_mle <- function(col1, col2) {
  s1 <- as.character(col1); s1[grepl("0", s1)] <- NA
  s2 <- as.character(col2); s2[grepl("0", s2)] <- NA
  m  <- !is.na(s1) & !is.na(s2)
  a  <- strsplit(s1[m], "/", fixed = TRUE)
  b  <- strsplit(s2[m], "/", fixed = TRUE)
  n  <- length(a)
  A  <- names(sort(table(unlist(a)), decreasing = TRUE))[1]
  B  <- names(sort(table(unlist(b)), decreasing = TRUE))[1]
  nA <- vapply(a, function(z) sum(z == A), 0L)
  nB <- vapply(b, function(z) sum(z == B), 0L)
  tb <- table(factor(nA, 0:2), factor(nB, 0:2))
  pA <- sum(nA) / (2*n); pB <- sum(nB) / (2*n)

  ll <- function(D) {
    p <- c(pA*pB + D, pA*(1-pB) - D, (1-pA)*pB - D, (1-pA)*(1-pB) + D)
    if (min(p) <= 0) return(-Inf)
    P <- matrix(0, 3, 3)
    P[3,3] <- p[1]^2;      P[3,2] <- 2*p[1]*p[2];              P[3,1] <- p[2]^2
    P[2,3] <- 2*p[1]*p[3]; P[2,2] <- 2*(p[1]*p[4] + p[2]*p[3]); P[2,1] <- 2*p[2]*p[4]
    P[1,3] <- p[3]^2;      P[1,2] <- 2*p[3]*p[4];              P[1,1] <- p[4]^2
    sum(tb * log(P))
  }
  lo <- max(-pA*pB, -(1-pA)*(1-pB)) + 1e-12
  hi <- min(pA*(1-pB), (1-pA)*pB) - 1e-12
  D  <- stats::optimize(ll, c(lo, hi), maximum = TRUE, tol = 1e-14)$maximum

  Dmax <- if (D > 0) min(pA*(1-pB), (1-pA)*pB) else max(-pA*pB, -(1-pA)*(1-pB))
  r    <- D / sqrt(pA*(1-pA)*pB*(1-pB))
  X2   <- 2 * n * r^2
  list(D = D, Dprime = D/Dmax, r = r, R2 = r^2, n = n,
       p = stats::pchisq(X2, 1, lower.tail = FALSE), loglik = ll)
}
