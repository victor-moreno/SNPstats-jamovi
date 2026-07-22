# Shared fixtures and verification oracles — loaded automatically by testthat.
# Run the suite with tests/run_tests.R (see that file for the library setup).
#
# The oracles below re-derive every quantity the module reports using only base R
# (glm / lm / nnet::multinom) and genetics/haplo.stats, independently of the
# module's internal functions, so the comparison is a genuine cross-check.

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
run_pgs <- function(...) do.call(SNPstats::snpPGS,   list(...))

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
