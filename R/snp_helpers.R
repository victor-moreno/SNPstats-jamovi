# ── snp_helpers.R ─────────────────────────────────────────────────────────────
#
# Shared utility functions for all three SNPstats analyses:
#   snpDescriptive  (snpDesc)
#   snpAssociation  (snpAssoc)
#   snpLDHaplotype  (snpLDHaplo)
#
# These functions handle:
#   - Genotype parsing and validation
#   - SNP column format detection
#   - Reference allele / genotype resolution
#   - Genetic model encoding
#   - Association and interaction model fitting
#
# Source this file at the top of each *Class.R:
#   source("snp_helpers.R")   (jamovi resolves relative to R/)
# ──────────────────────────────────────────────────────────────────────────────

#' @importFrom genetics genotype allele HWE.exact LD
#' @importFrom R6 R6Class
#' @import jmvcore


# ── Allele / genotype string utilities ────────────────────────────────────────

#' Split a normalised "A/B" genotype string into its two alleles.
split_alleles <- function(g) strsplit(g, "/", fixed = TRUE)[[1]]

#' Check that a vector of unique "A/B" genotype strings is biallelic and
#' well-formed.  Returns list($ok, $reason, $alleles).
check_biallelic <- function(vals) {
  pairs <- lapply(vals, split_alleles)
  bad   <- which(sapply(pairs, length) != 2)
  if (length(bad) > 0)
    return(list(ok = FALSE,
                reason = paste0("cannot split into two alleles: ",
                                paste(vals[bad], collapse = ", "))))

  alleles <- unique(unlist(pairs))
  if (length(alleles) > 2)
    return(list(ok = FALSE,
                reason = paste0("more than 2 alleles found (",
                                paste(sort(alleles), collapse = ", "),
                                "); only biallelic SNPs are supported")))

  a <- alleles[1]; b <- if (length(alleles) == 2) alleles[2] else alleles[1]
  valid    <- c(paste0(a,"/",a), paste0(a,"/",b),
                paste0(b,"/",a), paste0(b,"/",b))
  bad_geno <- vals[!vals %in% valid]
  if (length(bad_geno) > 0)
    return(list(ok = FALSE,
                reason = paste0("unexpected genotype(s): ",
                                paste(bad_geno, collapse = ", "))))

  list(ok = TRUE, reason = NULL, alleles = alleles)
}

#' Extract the user-defined genotype level order from a jamovi factor column.
get_snp_level_order <- function(x) {
  if (!is.factor(x)) return(NULL)
  lvls <- levels(x)
  if (length(lvls) == 0) return(NULL)

  norm <- lvls
  for (sep in c("|", ">")) {
    pat <- paste0("^.+", if (sep == "|") "\\|" else sep, ".+$")
    if (all(grepl(pat, norm))) {
      norm <- sub(sep, "/", norm, fixed = TRUE); break
    }
  }
  if (all(nchar(norm) == 2) && all(grepl("^[A-Za-z0-9]{2}$", norm)))
    norm <- paste0(substr(norm, 1, 1), "/", substr(norm, 2, 2))
  if (!all(grepl("^.+/.+$", norm))) return(NULL)

  ref_allele <- NULL
  for (g in norm) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[1] == parts[2]) { ref_allele <- parts[1]; break }
  }

  if (!is.null(ref_allele)) {
    norm <- sapply(norm, function(g) {
      parts <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2 && parts[1] != parts[2] && parts[2] == ref_allele)
        paste0(parts[2], "/", parts[1])
      else g
    }, USE.NAMES = FALSE)
  }

  is_het <- function(g) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    length(parts) == 2 && parts[1] != parts[2]
  }

  if (length(norm) == 3) {
    is_het_vec <- sapply(norm, is_het)
    hom_levels <- norm[!is_het_vec]; het_levels <- norm[is_het_vec]
    if (length(het_levels) == 1 && length(hom_levels) == 2) {
      ref_hom <- if (!is.null(ref_allele))
        hom_levels[sapply(hom_levels, function(g) {
          parts <- strsplit(g, "/", fixed = TRUE)[[1]]
          length(parts) == 2 && parts[1] == ref_allele })]
      else hom_levels[1]
      alt_hom <- hom_levels[hom_levels != ref_hom[1]]
      if (length(ref_hom) == 1 && length(alt_hom) == 1)
        norm <- c(ref_hom, het_levels, alt_hom)
    }
  } else {
    if (length(norm) >= 2 && is_het(norm[1])) norm[c(1, 2)] <- norm[c(2, 1)]
  }
  norm
}

#' Detect genotype separator; returns NULL if column is not valid biallelic SNP.
detect_snp_sep <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  if (length(vals) == 0 || length(vals) > 10) return(NULL)

  for (sep in c("/", "|", ">")) {
    pat  <- paste0("^.+", if (sep == "|") "\\|" else sep, ".+$")
    if (all(grepl(pat, vals))) {
      norm <- if (sep == "/") vals else sub(sep, "/", vals, fixed = TRUE)
      if (!isTRUE(check_biallelic(norm)$ok)) return(NULL)
      return(sep)
    }
  }
  if (all(nchar(vals) == 2) && all(grepl("^[A-Za-z0-9]{2}$", vals))) {
    norm <- paste0(substr(vals,1,1), "/", substr(vals,2,2))
    if (!isTRUE(check_biallelic(norm)$ok)) return(NULL)
    return("")
  }
  NULL
}

#' Full biallelic check (returns reason); used for validation messages.
snp_biallelic_check <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  sep  <- NULL
  for (s in c("/", "|", ">")) {
    pat <- paste0("^.+", if (s == "|") "\\|" else s, ".+$")
    if (all(grepl(pat, vals))) { sep <- s; break }
  }
  if (is.null(sep) && all(nchar(vals) == 2) &&
      all(grepl("^[A-Za-z0-9]{2}$", vals))) sep <- ""
  if (is.null(sep)) return(list(ok = FALSE, reason = "unrecognised format"))
  norm <- if (sep == "") paste0(substr(vals,1,1),"/",substr(vals,2,2)) else
          if (sep == "/") vals else sub(sep, "/", vals, fixed = TRUE)
  check_biallelic(norm)
}

#' TRUE if column looks like valid biallelic genotype data.
is_snp_column <- function(x) !is.null(detect_snp_sep(x))

#' Normalise and parse a raw genotype vector via genetics::genotype().
parse_genotype <- function(x, user_levels = NULL) {
  sep <- detect_snp_sep(x)
  if (is.null(sep)) return(NULL)

  x_chr <- as.character(x)
  if (sep == "") {
    x_norm <- ifelse(is.na(x_chr), NA_character_,
                     paste0(substr(x_chr,1,1), "/", substr(x_chr,2,2)))
  } else if (sep == "/") {
    x_norm <- x_chr
  } else {
    x_norm <- ifelse(is.na(x_chr), NA_character_,
                     sub(sep, "/", x_chr, fixed = TRUE))
  }

  ref_allele <- NULL
  if (!is.null(user_levels)) {
    for (g in user_levels) {
      parts <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2 && parts[1] == parts[2]) { ref_allele <- parts[1]; break }
    }
  }
  if (is.null(ref_allele)) {
    all_pairs   <- strsplit(as.character(na.omit(x_norm)), "/", fixed = TRUE)
    all_alleles <- unique(unlist(all_pairs))
    for (pr in all_pairs) {
      if (length(pr) == 2 && pr[1] == pr[2]) { ref_allele <- pr[1]; break }
    }
    if (is.null(ref_allele) && length(all_alleles) >= 1) ref_allele <- all_alleles[1]
  }

  if (!is.null(ref_allele)) {
    x_norm <- ifelse(is.na(x_norm), NA_character_,
      sapply(x_norm, function(g) {
        parts <- strsplit(g, "/", fixed = TRUE)[[1]]
        if (length(parts) == 2 && parts[1] != parts[2] && parts[2] == ref_allele)
          paste0(parts[2], "/", parts[1])
        else g
      }, USE.NAMES = FALSE))
  }

  tryCatch(genetics::genotype(x_norm, sep = "/"), error = function(e) NULL)
}

#' Determine reference genotype (user-specified first, then most-frequent homozygote).
get_ref_genotype <- function(geno, user_levels = NULL) {
  if (is.null(geno)) return(NULL)
  if (!is.null(user_levels) && length(user_levels) > 0) {
    sm  <- summary(geno)
    obs <- rownames(sm$genotype.freq)[rownames(sm$genotype.freq) != "NA"]
    for (lvl in user_levels) if (lvl %in% obs) return(lvl)
  }
  sm      <- summary(geno)
  gf      <- sm$genotype.freq
  alleles <- rownames(gf)
  is_hom  <- sapply(alleles, function(g) {
    p <- strsplit(g,"/")[[1]]; length(p)==2 && p[1]==p[2] })
  homz_gf <- gf[is_hom, , drop = FALSE]
  if (nrow(homz_gf) == 0) return(alleles[1])
  rownames(homz_gf)[which.max(homz_gf[,"Count"])]
}

#' Reorder genotype frequency table: user order first, then ref/het/alt fallback.
reorder_geno <- function(gf, ref, user_levels = NULL) {
  alleles <- rownames(gf)
  na_row  <- alleles == "NA"
  other   <- alleles[!na_row]

  canon_key <- function(g) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    if (length(parts) == 2) paste(sort(parts), collapse = "/") else g
  }

  if (!is.null(user_levels) && length(user_levels) > 0) {
    other_keys     <- setNames(sapply(other, canon_key), other)
    canon_to_actual <- setNames(names(other_keys), other_keys)
    ordered <- character(0)
    for (ul in user_levels) {
      ck <- canon_key(ul)
      if (ck %in% names(canon_to_actual)) {
        actual_nm <- canon_to_actual[[ck]]
        if (!actual_nm %in% ordered) ordered <- c(ordered, actual_nm)
      }
    }
    ordered <- c(ordered, other[!other %in% ordered])
  } else {
    is_hom  <- sapply(other, function(g) { p <- strsplit(g,"/")[[1]]; length(p)==2 && p[1]==p[2] })
    ordered <- c(ref,
                 other[!is_hom & other != ref],
                 other[ is_hom & other != ref])
    ordered <- unique(ordered[ordered %in% other])
  }

  final <- c(ordered, alleles[na_row])
  gf[final[final %in% alleles], , drop = FALSE]
}



# ── Shared validation helper ───────────────────────────────────────────────────

#' Validate SNP variables and return a list:
#'   $valid_snps  – character vector of SNP names passing validation
#'   $bad_html    – HTML string for the validation message (empty string if none)
validate_snp_vars <- function(snp_vars, data) {
  bad_snps <- character(0); bad_msgs <- character(0)
  for (v in snp_vars) {
    chk <- snp_biallelic_check(data[[v]])
    if (!isTRUE(chk$ok)) {
      bad_snps <- c(bad_snps, v)
      bad_msgs <- c(bad_msgs, paste0("<b>", v, "</b>: ", chk$reason))
    }
  }
  html <- if (length(bad_snps) > 0)
    paste0("<p style='color:red;'>The following SNP columns were skipped ",
           "(accepted formats: A/B, A|B, A>B, or AB; exactly 2 alleles required):</p>",
           "<ul>", paste0("<li>", bad_msgs, "</li>", collapse = ""), "</ul>")
  else ""

  list(valid_snps = setdiff(snp_vars, bad_snps), bad_html = html)
}

#' Detect response type ("binary" / "quantitative" / "none").
detect_response_type <- function(response_raw, responseType_opt) {
  if (responseType_opt != "auto") return(responseType_opt)
  if (is.null(response_raw)) return("none")
  n_unique <- length(unique(na.omit(response_raw)))
  if (n_unique == 2) "binary"
  else if (is.numeric(response_raw)) "quantitative"
  else if (n_unique > 2 & n_unique <= 6) "categorical"
  else "none"
}

#' Prepare response as integer (binary) or numeric (quantitative).
prepare_response <- function(response_raw, response_type) {
  if (is.null(response_raw) || is.null(response_type)) return(NULL)
  if (response_type == "binary") {
    r <- as.integer(as.factor(response_raw)) - 1L
    r[is.na(response_raw)] <- NA_integer_
    r
  } else {
    as.numeric(response_raw)
  }
}

#' Prepare covariate data frame (factor-encode character columns).
prepare_covariates <- function(data, covariate_vars) {
  if (length(covariate_vars) == 0) return(NULL)
  cov_df <- data[, covariate_vars, drop = FALSE]
  for (v in covariate_vars)
    if (!is.numeric(cov_df[[v]])) cov_df[[v]] <- as.factor(cov_df[[v]])
  cov_df
}



# ── Genetic model encoding ─────────────────────────────────────────────────────

#' Encode SNP under a given genetic model as a numeric/factor vector.
encode_model <- function(geno_char, ref, model, user_levels = NULL) {
  ref_allele <- strsplit(ref, "/")[[1]][1]
  dosage <- sapply(geno_char, function(g) {
    if (is.na(g) || g == "NA") return(NA_integer_)
    sum(strsplit(g, "/")[[1]] == ref_allele)
  })

  switch(model,
    codominant = {
      if (!is.null(user_levels) && length(user_levels) > 0) {
        obs_genos <- unique(geno_char[!is.na(geno_char) & geno_char != "NA"])
        lvls <- user_levels[user_levels %in% obs_genos]
        lvls <- c(lvls, obs_genos[!obs_genos %in% lvls])
      } else {
        lvls <- c(ref, unique(geno_char[geno_char != ref & !is.na(geno_char)]))
      }
      factor(geno_char, levels = lvls)
    },
    dominant    = ifelse(is.na(dosage), NA_integer_, as.integer(dosage < 2)),
    recessive   = ifelse(is.na(dosage), NA_integer_, as.integer(dosage == 0)),
    overdominant = {
      is_het <- sapply(geno_char, function(g) {
        if (is.na(g)) return(NA)
        length(unique(strsplit(g, "/")[[1]])) > 1
      })
      as.integer(is_het)
    },
    logadditive = 2L - as.integer(dosage)
  )
}

# ── Model fitting ──────────────────────────────────────────────────────────────

#' Fit association model for one SNP under one genetic model.
#' Returns list of per-comparison result lists.
fit_model <- function(snp_enc, response, covariates_df, model_name,
                      response_type, ci_width) {
  df <- data.frame(resp = response, snp = snp_enc)
  cov_formula <- ""
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df          <- cbind(df, covariates_df)
    cov_formula <- paste("+", paste(names(covariates_df), collapse = "+"))
  }
  df <- df[complete.cases(df), , drop = FALSE]

  formula_full <- as.formula(paste("resp ~ snp", cov_formula))
  formula_null <- as.formula(paste("resp ~ 1",   cov_formula))

  tryCatch({
    if (response_type == "binary") {
      fit_full   <- glm(formula_full, data = df, family = binomial())
      fit_null   <- glm(formula_null, data = df, family = binomial())
      lrtest     <- "Chisq"; lrtest_label <- "Pr(>Chi)"; pval_col <- "Pr(>|z|)"
    } else {
      fit_full   <- lm(formula_full, data = df)
      fit_null   <- lm(formula_null, data = df)
      lrtest     <- "F";     lrtest_label <- "Pr(>F)";   pval_col <- "Pr(>|t|)"
    }

    lrt      <- tryCatch(anova(fit_null, fit_full, test = lrtest), error = function(e) NULL)
    global_p <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
    aic_val  <- AIC(fit_full)
    coefs    <- summary(fit_full)$coefficients
    snp_rows <- grep("^snp", rownames(coefs))
    if (length(snp_rows) == 0) return(NULL)

    ci <- tryCatch(
      confint(fit_full, level = ci_width / 100)[snp_rows, , drop = FALSE],
      error = function(e) matrix(NA, nrow = length(snp_rows), ncol = 2))

    lapply(seq_along(snp_rows), function(i) {
      row  <- snp_rows[i]
      beta <- coefs[row, "Estimate"]
      pval <- coefs[row, pval_col]
      ci_lo <- ci[i, 1]; ci_hi <- ci[i, 2]
      if (response_type == "binary")
        list(effect = exp(beta), ci_low = exp(ci_lo), ci_high = exp(ci_hi),
             pval = pval, global_p = global_p, aic = aic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
      else
        list(effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, global_p = global_p, aic = aic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
    })
  }, error = function(e) NULL)
}

#' Fit SNP × covariate interaction model under one genetic model.
fit_interaction_model <- function(snp_enc, response, covariates_df,
                                  interaction_var, model_name,
                                  response_type, ci_width) {
  df <- data.frame(resp = response, snp = snp_enc)
  adj_covs <- character(0)
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df       <- cbind(df, covariates_df)
    adj_covs <- setdiff(names(covariates_df), interaction_var)
  }
  if (!(interaction_var %in% names(df))) return(NULL)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5) return(NULL)

  adj_part     <- if (length(adj_covs) > 0) paste("+", paste(adj_covs, collapse = "+")) else ""
  formula_int  <- as.formula(paste("resp ~ snp *", interaction_var, adj_part))
  formula_main <- as.formula(paste("resp ~ snp +", interaction_var, adj_part))

  tryCatch({
    if (response_type == "binary") {
      fit_int  <- glm(formula_int,  data = df, family = binomial())
      fit_main <- glm(formula_main, data = df, family = binomial())
      pval_col <- "Pr(>|z|)"; lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"
    } else {
      fit_int  <- lm(formula_int,  data = df)
      fit_main <- lm(formula_main, data = df)
      pval_col <- "Pr(>|t|)"; lrtest <- "F"; lrtest_label <- "Pr(>F)"
    }

    lrt     <- tryCatch(anova(fit_main, fit_int, test = lrtest), error = function(e) NULL)
    p_inter <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
    aic_val <- AIC(fit_int)
    coefs   <- summary(fit_int)$coefficients
    ci      <- tryCatch(confint(fit_int, level = ci_width / 100),
                        error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2,
                                                   dimnames = list(rownames(coefs), c("lo","hi"))))

    all_rows   <- rownames(coefs)
    snp_rows   <- grep("^snp", all_rows)
    inter_rows <- grep(paste0("^snp.*:", interaction_var, "|^", interaction_var, ":.*snp"), all_rows)
    keep_rows  <- unique(c(snp_rows, inter_rows))
    if (length(keep_rows) == 0) return(NULL)

    lapply(keep_rows, function(r) {
      beta  <- coefs[r, "Estimate"]
      pval  <- coefs[r, pval_col]
      ci_lo <- ci[r, 1]; ci_hi <- ci[r, 2]
      is_inter <- r %in% inter_rows
      if (response_type == "binary")
        list(term = all_rows[r], effect = exp(beta), ci_low = exp(ci_lo), ci_high = exp(ci_hi),
             pval = pval, pval_interaction = if (is_inter) p_inter else NA_real_,
             aic = aic_val, is_first = (r == keep_rows[1]))
      else
        list(term = all_rows[r], effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, pval_interaction = if (is_inter) p_inter else NA_real_,
             aic = aic_val, is_first = (r == keep_rows[1]))
    })
  }, error = function(e) NULL)
}

# formatting helper for categorical descriptives (N and %)
fmt_cat <- function(n, total) sprintf("%d (%.1f%%)", n, if (total > 0) 100*n/total else 0)
fmt_catpct <- function(n, pct) sprintf("%d (%.1f%%)", n, pct)
fmt_catn <- function(n) sprintf("%d", n)
fmt_cont <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%.2f \u00B1 %.2f", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE))
}
