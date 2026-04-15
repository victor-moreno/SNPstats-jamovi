#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype allele HWE.exact LD
#' @importFrom haplo.stats setupGeno hapl.em, haplo.glm
#' @import ggplot2


# ── Helpers ────────────────────────────────────────────────────────────────────

#' Split a normalised "A/B" genotype string into its two alleles.
split_alleles <- function(g) strsplit(g, "/", fixed = TRUE)[[1]]

#' Check that a vector of unique "A/B" genotype strings is biallelic and
#' well-formed (only AA, AB/BA, BB combinations for exactly two distinct
#' alleles).  Returns a list:
#'   $ok      TRUE/FALSE
#'   $reason  human-readable reason string when ok == FALSE
#'   $alleles character(2) canonical allele pair when ok == TRUE
check_biallelic <- function(vals) {
  # vals: unique non-NA genotype strings already normalised to "A/B" form
  pairs  <- lapply(vals, split_alleles)
  bad    <- which(sapply(pairs, length) != 2)
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

  # Every observed genotype must be one of AA, AB, BA, BB
  # (BA will be canonicalised to AB downstream)
  a <- alleles[1]; b <- if (length(alleles) == 2) alleles[2] else alleles[1]
  valid <- c(paste0(a,"/",a), paste0(a,"/",b),
             paste0(b,"/",a), paste0(b,"/",b))
  bad_geno <- vals[!vals %in% valid]
  if (length(bad_geno) > 0)
    return(list(ok = FALSE,
                reason = paste0("unexpected genotype(s): ",
                                paste(bad_geno, collapse = ", "))))

  list(ok = TRUE, reason = NULL, alleles = sort(alleles))
}

#' Detect if a character vector looks like diploid genotypes and return the
#' separator used, or NULL if the column does not look like genotype data or
#' fails the biallelic consistency check.
#'
#' Supported formats:
#'   A/B  A|B  A>B   (any allele names, separator is /, | or >)
#'   AB                (exactly 2 characters, no separator — single-char alleles)
detect_snp_sep <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  if (length(vals) == 0 || length(vals) > 10) return(NULL)

  # Explicit separators: any allele names allowed on either side
  for (sep in c("/", "|", ">")) {
    pat <- paste0("^.+", if (sep == "|") "\\|" else sep, ".+$")
    if (all(grepl(pat, vals))) {
      # Normalise to "/" for the biallelic check
      norm <- if (sep == "/") vals else sub(sep, "/", vals, fixed = TRUE)
      if (!isTRUE(check_biallelic(norm)$ok)) return(NULL)
      return(sep)
    }
  }

  # No-separator two-character format: each value is exactly 2 non-space chars
  if (all(nchar(vals) == 2) && all(grepl("^[A-Za-z0-9]{2}$", vals))) {
    norm <- paste0(substr(vals, 1, 1), "/", substr(vals, 2, 2))
    if (!isTRUE(check_biallelic(norm)$ok)) return(NULL)
    return("")
  }

  NULL
}

#' Return the biallelic check result (including reason) for a raw column.
#' Used by the validation step to produce a specific error message.
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

#' Convenience wrapper: TRUE if column looks like valid biallelic genotype data.
is_snp_column <- function(x) !is.null(detect_snp_sep(x))

#' Normalise a raw genotype vector to canonical "A/B" format (A <= B
#' alphabetically, so B/A becomes A/B), then parse via genetics::genotype().
#' Returns NULL if the format cannot be determined or the column is not
#' biallelic.
parse_genotype <- function(x) {
  sep <- detect_snp_sep(x)
  if (is.null(sep)) return(NULL)

  x_chr <- as.character(x)

  # Step 1: convert to slash-separated
  if (sep == "") {
    x_norm <- ifelse(is.na(x_chr), NA_character_,
                     paste0(substr(x_chr, 1, 1), "/", substr(x_chr, 2, 2)))
  } else if (sep == "/") {
    x_norm <- x_chr
  } else {
    x_norm <- ifelse(is.na(x_chr), NA_character_,
                     sub(sep, "/", x_chr, fixed = TRUE))
  }

  # Step 2: canonicalise allele order so A/B and B/A become the same genotype.
  # Use alphabetical order of the two alleles as the canonical form.
  alleles <- sort(unique(unlist(strsplit(
    as.character(na.omit(x_norm)), "/", fixed = TRUE))))
  if (length(alleles) == 2) {
    a1 <- alleles[1]; a2 <- alleles[2]   # a1 <= a2 alphabetically
    x_norm <- ifelse(
      is.na(x_norm), NA_character_,
      ifelse(x_norm == paste0(a2, "/", a1),
             paste0(a1, "/", a2),
             x_norm)
    )
  }

  tryCatch(
    genetics::genotype(x_norm, sep = "/"),
    error = function(e) NULL
  )
}

#' Determine reference genotype (most frequent homozygote)
get_ref_genotype <- function(geno) {
  if (is.null(geno)) return(NULL)
  sm <- summary(geno)
  gf <- sm$genotype.freq
  # Homozygotes: allele1 == allele2
  alleles <- rownames(gf)
  is_homoz <- sapply(alleles, function(g) {
    parts <- strsplit(g, "/")[[1]]
    length(parts) == 2 && parts[1] == parts[2]
  })
  homoz_gf <- gf[is_homoz, , drop = FALSE]
  if (nrow(homoz_gf) == 0) return(alleles[1])
  rownames(homoz_gf)[which.max(homoz_gf[, "Count"])]
}

#' Reorder genotype frequency table: ref homozygote first, then het, then alt
reorder_geno <- function(gf, ref) {
  alleles <- rownames(gf)
  na_row  <- alleles == "NA"
  other   <- alleles[!na_row]
  # ref first, then hets, then other homozygotes
  is_homoz <- sapply(other, function(g) {
    parts <- strsplit(g, "/")[[1]]
    length(parts) == 2 && parts[1] == parts[2]
  })
  ordered <- c(
    ref,
    other[!is_homoz & other != ref],
    other[is_homoz  & other != ref]
  )
  ordered <- unique(ordered[ordered %in% other])
  final <- c(ordered, alleles[na_row])
  gf[final[final %in% alleles], , drop = FALSE]
}

#' Encode SNP under a given genetic model as a numeric/factor vector
encode_model <- function(geno_char, ref, model) {
  # alleles of reference homozygote
  ref_allele <- strsplit(ref, "/")[[1]][1]
  # count ref alleles per genotype
  dosage <- sapply(geno_char, function(g) {
    if (is.na(g) || g == "NA") return(NA_integer_)
    parts <- strsplit(g, "/")[[1]]
    sum(parts == ref_allele)
  })

  switch(model,
    codominant = {
      # factor with ref as first level
      lvls <- c(ref,
                unique(geno_char[geno_char != ref & !is.na(geno_char)]))
      factor(geno_char, levels = lvls)
    },
    dominant = {
      # 0 = ref homozygote, 1 = het or alt homozygote
      ifelse(is.na(dosage), NA_integer_, as.integer(dosage < 2))
    },
    recessive = {
      # 0 = ref + het, 1 = alt homozygote
      ifelse(is.na(dosage), NA_integer_, as.integer(dosage == 0))
    },
    overdominant = {
      # 0 = homozygote (either), 1 = heterozygote
      is_het <- sapply(geno_char, function(g) {
        if (is.na(g)) return(NA)
        parts <- strsplit(g, "/")[[1]]
        length(unique(parts)) > 1
      })
      as.integer(is_het)
    },
    logadditive = {
      # 0/1/2 copies of alt allele
      2L - as.integer(dosage)
    }
  )
}

#' Fit association model for one SNP under one genetic model
#' Returns list(effect, ci_low, ci_high, pval, global_p, comparison, aic)
fit_model <- function(snp_enc, response, covariates_df, model_name,
                      response_type, ci_width) {
  alpha <- 1 - ci_width / 100

  df <- data.frame(resp = response, snp = snp_enc)
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df <- cbind(df, covariates_df)

    cov_formula <- paste("+", paste(names(covariates_df), collapse = "+"))
  } else {
    cov_formula <- ""
  }

  # remove missings
  df <- df[complete.cases(df), , drop = FALSE]

  formula_full <- as.formula(paste("resp ~ snp", cov_formula))
  formula_null <- as.formula(paste("resp ~ 1", cov_formula))

  tryCatch({
    if (response_type == "binary") {
      fit_full <- glm(formula_full, data = df, family = binomial())
      fit_null <- glm(formula_null, data = df, family = binomial())
      lrtest <-'Chisq'
      lrtest_label <- 'Pr(>Chi)'
    } else {
      fit_full <- lm(formula_full, data = df)
      fit_null <- lm(formula_null, data = df)
      lrtest <-'F'
      lrtest_label <- 'Pr(>F)'
    }

    # Global LRT p-value
    lrt <- tryCatch(
      anova(fit_null, fit_full, test = lrtest),
      error = function(e) NULL
    )
    global_p <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
    
    aic_val <- AIC(fit_full)

    coefs <- summary(fit_full)$coefficients
    snp_rows <- grep("^snp", rownames(coefs))

    if (length(snp_rows) == 0) {
      return(NULL)
    }

    ci <- tryCatch(
      confint(fit_full, level = ci_width / 100)[snp_rows, , drop = FALSE],
      error = function(e) matrix(NA, nrow = length(snp_rows), ncol = 2)
    )

    results <- lapply(seq_along(snp_rows), function(i) {
      row <- snp_rows[i]
      beta <- coefs[row, "Estimate"]
      pval <- coefs[row, ifelse(response_type == "binary",
                                 "Pr(>|z|)", "Pr(>|t|)")]
      ci_lo <- ci[i, 1]
      ci_hi <- ci[i, 2]

      if (response_type == "binary") {
        list(
          effect   = exp(beta),
          ci_low   = exp(ci_lo),
          ci_high  = exp(ci_hi),
          pval     = pval,
          global_p = global_p,
          aic      = aic_val,
          comparison = sub("^snp", "", rownames(coefs)[row])
        )
      } else {
        list(
          effect   = beta,
          ci_low   = ci_lo,
          ci_high  = ci_hi,
          pval     = pval,
          global_p = global_p,
          aic      = aic_val,
          comparison = sub("^snp", "", rownames(coefs)[row])
        )
      }
    })
    results
  }, error = function(e) NULL)
}


#' Fit SNP × covariate interaction model under one genetic model.
#' Returns a list of rows: main SNP terms, interaction term(s), plus
#' a p_interaction (LRT of model-with-interaction vs model-without).
fit_interaction_model <- function(snp_enc, response, covariates_df,
                                  interaction_var, model_name,
                                  response_type, ci_width) {
  alpha <- 1 - ci_width / 100

  df <- data.frame(resp = response, snp = snp_enc)
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df <- cbind(df, covariates_df)
    adj_covs <- setdiff(names(covariates_df), interaction_var)
  } else {
    adj_covs <- character(0)
  }

  if (!(interaction_var %in% names(df))) return(NULL)

  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5) return(NULL)

  adj_part <- if (length(adj_covs) > 0)
    paste("+", paste(adj_covs, collapse = "+"))
  else ""

  # model with interaction
  formula_int  <- as.formula(
    paste("resp ~ snp *", interaction_var, adj_part))
  # model without interaction (for LRT)
  formula_main <- as.formula(
    paste("resp ~ snp +", interaction_var, adj_part))

  tryCatch({
    if (response_type == "binary") {
      fit_int  <- glm(formula_int,  data = df, family = binomial())
      fit_main <- glm(formula_main, data = df, family = binomial())
      pval_col <- "Pr(>|z|)"
      lrtest   <- "Chisq"; lrtest_label <- "Pr(>Chi)"
    } else {
      fit_int  <- lm(formula_int,  data = df)
      fit_main <- lm(formula_main, data = df)
      pval_col <- "Pr(>|t|)"
      lrtest   <- "F"; lrtest_label <- "Pr(>F)"
    }

    lrt      <- tryCatch(anova(fit_main, fit_int, test = lrtest),
                         error = function(e) NULL)
    p_inter  <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
    aic_val  <- AIC(fit_int)

    coefs <- summary(fit_int)$coefficients
    ci    <- tryCatch(
      confint(fit_int, level = ci_width / 100),
      error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2,
                                  dimnames = list(rownames(coefs), c("lo","hi")))
    )

    # rows of interest: SNP main effect(s) and interaction term(s)
    all_rows    <- rownames(coefs)
    snp_rows    <- grep("^snp",                  all_rows)
    inter_rows  <- grep(paste0("^snp.*:", interaction_var,
                               "|^", interaction_var, ":.*snp"),
                        all_rows)
    keep_rows   <- unique(c(snp_rows, inter_rows))
    if (length(keep_rows) == 0) return(NULL)

    results <- lapply(keep_rows, function(r) {
      beta  <- coefs[r, "Estimate"]
      pval  <- coefs[r, pval_col]
      ci_lo <- ci[r, 1]
      ci_hi <- ci[r, 2]
      term  <- all_rows[r]
      is_inter_term <- r %in% inter_rows

      if (response_type == "binary") {
        list(term            = term,
             effect          = exp(beta),
             ci_low          = exp(ci_lo),
             ci_high         = exp(ci_hi),
             pval            = pval,
             pval_interaction = if (is_inter_term) p_inter else NA_real_,
             aic             = aic_val,
             is_first        = (r == keep_rows[1]))
      } else {
        list(term            = term,
             effect          = beta,
             ci_low          = ci_lo,
             ci_high         = ci_hi,
             pval            = pval,
             pval_interaction = if (is_inter_term) p_inter else NA_real_,
             aic             = aic_val,
             is_first        = (r == keep_rows[1]))
      }
    })
    results
  }, error = function(e) NULL)
}


# ── Main analysis class ────────────────────────────────────────────────────────

snpAnalysisClass <- if (requireNamespace("jmvcore", quietly=TRUE)) R6::R6Class(
  "snpAnalysisClass",
  inherit = snpAnalysisBase,
  private = list(

    .init = function() {
      # Hide all optional groups immediately — shown only when options are on
      self$results$covDescGroup$setVisible(FALSE)
      self$results$snpSummaryTablesGroup$setVisible(isTRUE(self$options$snpSummary))
      self$results$ldGroup$setVisible(FALSE)
      self$results$haploGroup$setVisible(FALSE)
      
      # Initialise the per-SNP array when SNPs are assigned
      snp_names <- self$options$snps
      if (length(snp_names) == 0) return()

      arr <- self$results$snpResults
      for (nm in snp_names) {
        arr$addItem(key = nm)
      }
    },

    .run = function() {
      data      <- self$data
      opts      <- self$options
      response_var   <- opts$response
      snp_vars       <- opts$snps
      covariate_vars <- opts$covariates

      # Initialize flags for what will actually run
      run_snpSummary <- isTRUE(opts$snpSummary)
      run_allFreq <- isTRUE(opts$allFreq)
      run_genoFreq <- isTRUE(opts$genoFreq)
      run_hweTest <- isTRUE(opts$hweTest)
      run_ldAnalysis <- isTRUE(opts$ldAnalysis)
      run_ldMatrix <- isTRUE(opts$ldMatrix)
      run_ldPlot <- isTRUE(opts$ldPlot)
      run_snpAssoc <- isTRUE(opts$snpAssoc)
      run_snpInteraction <- isTRUE(opts$snpInteraction)
      run_haploFreq <- isTRUE(opts$haploFreq)
      run_haploAssoc <- isTRUE(opts$haploAssoc)
      run_haploInteraction <- isTRUE(opts$haploInteraction)
      run_subpop <- isTRUE(opts$subpop)
      run_covDesc <- isTRUE(opts$covDesc)
      run_showMissing <- isTRUE(opts$showMissing)

      # ── Validation 1: SNPs required ──────────────────────────────
      if (length(snp_vars) == 0) {
        self$results$validationMsgSNP$setContent(paste0(
            "<p style='color:red;'> Please add at least one SNP variable.</p>"))
        self$results$validationMsgSNP$setVisible(TRUE)
        run_snpSummary <- FALSE
        run_ldAnalysis <- FALSE
        run_ldMatrix <- FALSE
        run_ldPlot <- FALSE
        run_snpAssoc <- FALSE
        run_snpInteraction <- FALSE
        run_haploFreq <- FALSE
        run_haploAssoc <- FALSE
        run_haploInteraction <- FALSE
      } else {
        self$results$validationMsgSNP$setVisible(FALSE)
      }

      # ── Validation 2: Minimum SNPs for LD/haplotype ──────────────
      if (length(snp_vars) < 2) {
        if (run_ldAnalysis || run_ldMatrix || run_ldPlot || 
            run_haploFreq || run_haploAssoc) {
          self$results$validationMsg$setContent(
            "<p style='color:red;'> LD and haplotype analyses require at least 2 SNPs.</p>")
          self$results$validationMsg$setVisible(TRUE)
          # Disable by setting flags to FALSE
          run_ldAnalysis <- run_ldMatrix <- run_ldPlot <- 
            run_haploFreq <- run_haploAssoc <- FALSE
        }
      } else {
        # Hide validation message if conditions are met
        if (!run_ldAnalysis && !run_ldMatrix && !run_ldPlot &&
            !run_haploFreq && !run_haploAssoc) {
          self$results$validationMsg$setVisible(FALSE)
        }
      }

      # ── Validation 3: Response required for certain tests ────────
      if ((run_snpAssoc || run_subpop) && 
          (is.null(response_var) || response_var == "")) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'> Association tests and stratification require a response variable.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpAssoc <- run_subpop <- FALSE
      }

      # ── Determine response type ──────────────────────────────────
      response_raw <- if (!is.null(response_var) && response_var != "") 
                        data[[response_var]] else NULL
        
      response_type <- opts$responseType
      if (!is.null(response_raw) && response_type == "auto") {
        n_unique <- length(unique(na.omit(response_raw)))
        if (n_unique == 2) {
          response_type <- "binary"
        } else if (is.numeric(response_raw)) {
          response_type <- "quantitative"
        } else {
          response_type <- NULL
        }
      }

      # ── Validation 4: Subpop compatibility ───────────────────────
      if (run_subpop && response_type == "quantitative") {
        self$results$validationMsg$setContent(
          "<p style='color:red;'> Subpopulation analysis is only available for binary responses.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_subpop <- FALSE
      }

      # ── Validation 5: Covariates required for interactions ───────
      if (run_snpInteraction && length(covariate_vars) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'> SNP × covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }
      
      if (run_haploInteraction && length(covariate_vars) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'> Haplotype × covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploInteraction <- FALSE
      }

      # ── Validation 6: Haplotype association requires response ────
      if (run_haploAssoc && (is.null(response_var) || response_var == "")) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'> Haplotype association requires a response variable.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploAssoc <- FALSE
      }

      # ── Validation 7: Check SNP columns format and biallelic consistency ──
      bad_snps  <- character(0)
      bad_msgs  <- character(0)
      for (v in snp_vars) {
        chk <- snp_biallelic_check(data[[v]])
        if (!isTRUE(chk$ok)) {
          bad_snps <- c(bad_snps, v)
          bad_msgs <- c(bad_msgs,
                        paste0("<b>", v, "</b>: ", chk$reason))
        }
      }
      if (length(bad_snps) > 0) {
        self$results$validationMsgGeno$setContent(paste0(
          "<p style='color:red;'>The following SNP columns were skipped ",
          "(accepted formats: A/B, A|B, A>B, or AB; exactly 2 alleles required):</p>",
          "<ul>", paste0("<li>", bad_msgs, "</li>", collapse = ""), "</ul>"))
        snp_vars <- setdiff(snp_vars, bad_snps)
        self$results$validationMsgGeno$setVisible(TRUE)

        if (length(snp_vars) == 0) {
          self$results$validationMsgSNP$setContent(
            "<p style='color:red;'> No valid SNP columns found. Please check your data format.</p>")
          self$results$validationMsgSNP$setVisible(TRUE)
          return()
        }
      } else {
        self$results$validationMsgGeno$setVisible(FALSE)
      }

      # ── Prepare response ─────────────────────────────────────────
      if (!is.null(response_raw) && response_var != "") {
        if (response_type == "binary") {
          response <- as.integer(as.factor(response_raw)) - 1L
          response[is.na(response_raw)] <- NA_integer_
        } else {
          response <- as.numeric(response_raw)
        }
      } else {
        response <- NULL
      }

      # ── Covariates ───────────────────────────────────────────────
      if (length(covariate_vars) > 0) {
        cov_df <- data[, covariate_vars, drop = FALSE]
        for (v in covariate_vars) {
          if (!is.numeric(cov_df[[v]])) {
            cov_df[[v]] <- as.factor(cov_df[[v]])
          }
        }
      } else {
        cov_df <- NULL
      }

      # ── Show/hide optional result groups based on flags ──────────
      self$results$covDescGroup$setVisible(run_covDesc && length(covariate_vars) > 0)
      self$results$snpSummaryTablesGroup$setVisible(run_snpSummary)
      self$results$ldGroup$setVisible(run_ldAnalysis || run_ldMatrix || run_ldPlot)
      self$results$haploGroup$setVisible(run_haploFreq || run_haploAssoc || run_haploInteraction)

      # ── Covariate descriptives ───────────────────────────────────
      if (run_covDesc && !is.null(cov_df) && length(covariate_vars) > 0) {
        private$.run_cov_desc(cov_df, response_raw, response_type, run_subpop)
      }

      # ── SNP summary table ─────────────────────────────────────────
      if (run_snpSummary) {
        private$.fill_snp_summary(data, snp_vars, response_raw, response_type, run_subpop, cov_df)
      }

      # ── Complete-case mask across response + covariates ─────────────
      # SNP-level NAs are handled per-SNP below; this mask covers everything
      # except SNP missingness so all descriptive tables use the same sample
      # as the models.
      n_rows <- nrow(data)
      complete_mask <- rep(TRUE, n_rows)
      if (!is.null(response))                    complete_mask <- complete_mask & !is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0) complete_mask <- complete_mask & complete.cases(cov_df)

      # ── Per-SNP analyses ─────────────────────────────────────────
      arr <- self$results$snpResults
      geno_list <- list()

      for (snp_nm in snp_vars) {
        snp_raw  <- data[[snp_nm]]
        geno_obj <- parse_genotype(snp_raw)
        if (is.null(geno_obj)) next

        geno_list[[snp_nm]] <- geno_obj

        # Per-SNP complete-case mask: response + covariates + this SNP
        snp_complete_mask <- complete_mask & !is.na(snp_raw)

        # Number of obs that have valid response+covariates but missing SNP
        # VM calculate a vector for each strata of response or total if no stratification
        if (run_subpop && response_type == "binary" && !is.null(response_raw)) {
          n_snp_missing <- tapply(snp_raw, response_raw, function(x) sum(is.na(x)))
          names(n_snp_missing) <- levels(as.factor(response_raw))
        } else {
          n_snp_missing <- sum(complete_mask & is.na(snp_raw))
        }
        

        # Restrict all descriptive objects to the analysis sample
        snp_raw_cc  <- snp_raw[snp_complete_mask]
        geno_obj_cc <- parse_genotype(snp_raw_cc)
        response_cc <- if (!is.null(response))     response[snp_complete_mask]     else NULL
        resp_raw_cc <- if (!is.null(response_raw)) response_raw[snp_complete_mask] else NULL
        cov_df_cc   <- if (!is.null(cov_df))       cov_df[snp_complete_mask, , drop = FALSE] else NULL

        if (is.null(geno_obj_cc)) next

        item <- arr$get(key = snp_nm)

        # Typed samples reflects the analysis sample (all variables non-missing)
        n_typed <- sum(snp_complete_mask)
        n_total <- n_rows
        pct     <- round(n_typed / n_total * 100, 1)
        item$typingRate$setContent(sprintf(
          "<b>Typed samples:</b> %d / %d (%.1f%%)", n_typed, n_total, pct))

        snp_summary_cc <- summary(geno_obj_cc)
        ref <- get_ref_genotype(geno_obj_cc)

        # Allele frequencies
        if (run_allFreq) {
          private$.fill_allele_freq(item$allFreqTable, snp_summary_cc,
                                    snp_nm, resp_raw_cc, run_subpop,
                                    response_type, snp_raw_cc,
                                    run_showMissing, n_snp_missing)
        }

        # Genotype frequencies
        if (run_genoFreq) {
          private$.fill_geno_freq(item$genoFreqTable, snp_summary_cc, ref,
                                  snp_raw_cc, response_cc, response_type,
                                  run_subpop, resp_raw_cc,
                                  run_showMissing, n_snp_missing)
        }

        # HWE
        if (run_hweTest) {
          private$.fill_hwe(item$hweTable, geno_obj_cc, snp_nm,
                            resp_raw_cc, run_subpop,
                            run_showMissing, n_snp_missing)
        }

        # Association (fit_model does its own complete.cases internally,
        # but inputs are already restricted so the sample is consistent)
        if (run_snpAssoc && !is.null(response_cc)) {
          private$.fill_assoc(item$assocTable, snp_raw_cc, ref, response_cc,
                              cov_df_cc, response_type, opts, run_snpAssoc)
        }

        # SNP x covariate interaction
        if (run_snpInteraction && !is.null(cov_df_cc) && ncol(cov_df_cc) >= 1 &&
            !is.null(response_cc)) {
          private$.fill_interaction(item$interactionTable, snp_raw_cc, ref,
                                    response_cc, cov_df_cc, names(cov_df_cc)[1],
                                    response_type, opts)
        }
      }

      # ── LD analysis ──────────────────────────────────────────────
      if ((run_ldAnalysis || run_ldMatrix || run_ldPlot) && length(geno_list) >= 2) {
        private$.run_ld(geno_list, opts, run_ldAnalysis, run_ldMatrix, run_ldPlot)
      }

      # ── Haplotype analysis ───────────────────────────────────────
      if ((run_haploFreq || run_haploAssoc || run_haploInteraction) && length(geno_list) >= 2) {
        private$.run_haplo(geno_list, data, response, response_raw, response_type, cov_df,
                          opts, run_haploFreq, run_haploAssoc, run_haploInteraction,
                          run_subpop, complete_mask)
      }
    },

    # ── Covariate descriptives ────────────────────────────────────
    .run_cov_desc = function(cov_df, response_raw, response_type, subpop = FALSE) {
      tbl      <- self$results$covDescGroup$covDescTable
      do_strat <- isTRUE(subpop) && response_type == "binary" && !is.null(response_raw)

      if (do_strat) {
        grp_lvls <- levels(as.factor(response_raw))
        if (length(grp_lvls) == 2) {
          tbl$getColumn("stat_g0")$setTitle(as.character(grp_lvls[1]))
          tbl$getColumn("stat_g1")$setTitle(as.character(grp_lvls[2]))
          tbl$getColumn("stat_g0")$setVisible(TRUE)
          tbl$getColumn("stat_g1")$setVisible(TRUE)
          tbl$getColumn("pval")$setVisible(TRUE)
        } else {
          do_strat <- FALSE
        }
      }

      fmt_cat  <- function(n, total)
        sprintf("%d (%.1f%%)", n, if (total > 0) n / total * 100 else 0)
      fmt_cont <- function(x)
        sprintf("%.2f \u00B1 %.2f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

      has_cont <- FALSE   # track whether any continuous variable is present

      for (v in names(cov_df)) {
        col    <- cov_df[[v]]
        is_cat <- is.factor(col) || is.character(col)
        if (is_cat) col <- as.factor(col)
        n_miss <- sum(is.na(col))

        if (is_cat) {
          lvls     <- levels(col)
          grp_fac  <- if (do_strat) as.factor(response_raw) else NULL
          # p-value: chi-square excluding missing as a category
          pval_cat <- if (do_strat) tryCatch({
            ct <- table(col[!is.na(response_raw)],
                        grp_fac[!is.na(response_raw)],
                        useNA = "no")
            suppressWarnings(chisq.test(ct)$p.value)
          }, error = function(e) NA_real_) else NA_real_

          first_row <- TRUE
          # ── non-missing levels ────────────────────────────────────
          for (lvl in lvls) {
            mask    <- !is.na(col) & col == lvl
            tot_all <- length(col)            # denominator = all rows incl. missing

            row_vals <- list(
              variable     = if (first_row) v else "",
              level        = lvl,
              stat_overall = fmt_cat(sum(mask), tot_all)
            )

            if (do_strat) {
              mask0 <- mask & !is.na(response_raw) & grp_fac == grp_lvls[1]
              mask1 <- mask & !is.na(response_raw) & grp_fac == grp_lvls[2]
              tot0  <- sum(!is.na(response_raw) & grp_fac == grp_lvls[1])
              tot1  <- sum(!is.na(response_raw) & grp_fac == grp_lvls[2])
              row_vals$stat_g0 <- fmt_cat(sum(mask0), tot0)
              row_vals$stat_g1 <- fmt_cat(sum(mask1), tot1)
              row_vals$pval    <- if (first_row) pval_cat else ''
            }

            tbl$addRow(rowKey = paste0(v, "_", lvl), values = row_vals)
            first_row <- FALSE
          }

          # ── Missing category (if any) ─────────────────────────────
          if (n_miss > 0) {
            miss_vals <- list(
              variable     = "",
              level        = "Missing",
              stat_overall = fmt_cat(n_miss, length(col))
            )
            if (do_strat) {
              miss_mask0 <- is.na(col) & !is.na(response_raw) & grp_fac == grp_lvls[1]
              miss_mask1 <- is.na(col) & !is.na(response_raw) & grp_fac == grp_lvls[2]
              tot0 <- sum(!is.na(response_raw) & grp_fac == grp_lvls[1])
              tot1 <- sum(!is.na(response_raw) & grp_fac == grp_lvls[2])
              miss_vals$stat_g0 <- fmt_cat(sum(miss_mask0), tot0)
              miss_vals$stat_g1 <- fmt_cat(sum(miss_mask1), tot1)
              miss_vals$pval    <- ''
            }
            tbl$addRow(rowKey = paste0(v, "_missing"), values = miss_vals)
          }

        } else {
          has_cont <- TRUE
          row_vals <- list(
            variable     = v,
            level        = "Mean \u00B1 SD",
            stat_overall = fmt_cont(col)
          )
          if (do_strat) {
            grp_fac  <- as.factor(response_raw)
            g0       <- col[!is.na(response_raw) & grp_fac == grp_lvls[1]]
            g1       <- col[!is.na(response_raw) & grp_fac == grp_lvls[2]]
            row_vals$stat_g0 <- fmt_cont(g0)
            row_vals$stat_g1 <- fmt_cont(g1)
            row_vals$pval    <- tryCatch(t.test(g0, g1)$p.value,
                                         error = function(e) NA_real_)
          }
          tbl$addRow(rowKey = v, values = row_vals)

          # ── Missing row for continuous (if any) ───────────────────
          if (n_miss > 0) {
            miss_vals <- list(
              variable     = "",
              level        = "Missing",
              stat_overall = fmt_cat(n_miss, length(col))
            )
            if (do_strat) {
              grp_fac   <- as.factor(response_raw)
              miss_mask <- is.na(col)
              miss_mask0 <- miss_mask & !is.na(response_raw) & grp_fac == grp_lvls[1]
              miss_mask1 <- miss_mask & !is.na(response_raw) & grp_fac == grp_lvls[2]
              tot0 <- sum(!is.na(response_raw) & grp_fac == grp_lvls[1])
              tot1 <- sum(!is.na(response_raw) & grp_fac == grp_lvls[2])
              miss_vals$stat_g0 <- fmt_cat(sum(miss_mask0), tot0)
              miss_vals$stat_g1 <- fmt_cat(sum(miss_mask1), tot1)
              miss_vals$pval    <- ''
            }
            tbl$addRow(rowKey = paste0(v, "_missing"), values = miss_vals)
          }
        }
      }

      # ── Table note for continuous variables ───────────────────────
      if (has_cont)
        tbl$setNote(note = "Continuous variables: mean \u00B1 SD.",
                    key  = "cont_fmt")
      else
        tbl$setNote(note = NULL, key = "cont_fmt")
    },

    # ── SNP summary table ─────────────────────────────────────────
    .fill_snp_summary = function(data, snp_vars, response_raw, response_type, subpop,
                                  cov_df = NULL) {
      tbl <- self$results$snpSummaryTablesGroup$snpSummaryTable

      do_strat <- isTRUE(subpop) &&
                  !is.null(response_raw) &&
                  identical(response_type, "binary")

      grp_levels <- character(0)
      if (do_strat) {
        grp_levels <- sort(unique(na.omit(as.character(response_raw))))
        if (length(grp_levels) != 2L) do_strat <- FALSE
      }

      tbl$getColumn("group")$setVisible(do_strat)

      has_response <- !is.null(response_raw)
      has_cov      <- !is.null(cov_df) && ncol(cov_df) > 0
      n_total_rows <- nrow(data)

      # Base complete-case mask: response + covariates (SNP added per-SNP below)
      base_cc <- rep(TRUE, n_total_rows)
      if (has_response) base_cc <- base_cc & !is.na(response_raw)
      if (has_cov)      base_cc <- base_cc & complete.cases(cov_df)

      row_key <- 0L

      for (snp_nm in snp_vars) {
        snp_raw  <- data[[snp_nm]]
        geno_obj <- parse_genotype(snp_raw)
        if (is.null(geno_obj)) next

        # Full complete-case mask for this SNP
        cc_mask    <- base_cc & !is.na(snp_raw)
        n_cc       <- sum(cc_mask)
        n_excluded <- n_total_rows - n_cc   # total excluded (SNP + resp + cov)
        snp_cc  <- snp_raw[cc_mask]
        geno_cc <- parse_genotype(snp_cc)
        if (is.null(geno_cc)) next

        resp_raw_cc <- if (has_response) response_raw[cc_mask] else NULL

        # for strata missingness counts below; also used for stratified stats if do_strat
        stratum_totals <- if (do_strat) {
          table(factor(response_raw[base_cc], levels = grp_levels))
        } else NULL
        
        sm_cc <- summary(geno_cc)
        ref   <- get_ref_genotype(geno_cc)

        # Allele labels from the complete-case summary
        af_all     <- sm_cc$allele.freq
        allele_nms <- rownames(af_all)
        ref_allele <- strsplit(ref, "/")[[1]][1]
        alt_allele <- allele_nms[allele_nms != ref_allele]
        alt_allele <- if (length(alt_allele) > 0) alt_allele[1] else "?"
        alleles_label <- paste0(ref_allele, "/", alt_allele)

        # Helper: stats for a subset of the already-cc vectors
        compute_row <- function(sub_mask) {
          snp_sub  <- if (is.null(sub_mask)) snp_cc else snp_cc[sub_mask]
          geno_sub <- if (is.null(sub_mask)) geno_cc else parse_genotype(snp_sub)
          if (is.null(geno_sub)) return(NULL)

          sm_sub  <- summary(geno_sub)
          n_typed <- sm_sub$n.typed

          af    <- sm_sub$allele.freq
          props <- af[, "Proportion"]
          maf   <- if (alt_allele %in% rownames(af)) {
            af[alt_allele, "Proportion"]
          } else if (length(props) >= 2) {
            min(props, na.rm = TRUE)
          } else NA_real_

          gf <- sm_sub$genotype.freq
          gf <- tryCatch(reorder_geno(gf, ref), error = function(e) gf)
          gf <- gf[rownames(gf) != "NA", , drop = FALSE]
          counts <- as.integer(gf[, "Count"])
          geno_str <- if (length(counts) == 3) {
            paste(counts, collapse = " / ")
          } else if (length(counts) == 2) {
            paste(c(counts, 0L), collapse = " / ")
          } else {
            paste(counts, collapse = " / ")
          }

          hwe_p <- tryCatch({
            hw <- genetics::HWE.exact(geno_sub)
            hw$p.value
          }, error = function(e) NA_real_)
          
          list(n = n_typed, maf = maf, genoCounts = geno_str, hwePval = hwe_p)
        }

        # Overall row (complete cases only)
        res_all <- compute_row(NULL)
        if (!is.null(res_all)) {
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            snp        = snp_nm,
            alleles    = alleles_label,
            group      = if (do_strat) "All" else "",
            n          = res_all$n,
            missing    = if (n_excluded > 0L) n_excluded else '',
            maf        = round(res_all$maf, 4),
            genoCounts = res_all$genoCounts,
            hwePval    = res_all$hwePval
          ))
        }

        # Stratified rows — subset the already-cc vectors by response level
        if (do_strat) {
          for (lvl in grp_levels) {
            sub_mask <- !is.na(resp_raw_cc) & as.character(resp_raw_cc) == lvl
            res_lvl  <- compute_row(sub_mask)
            if (is.null(res_lvl)) next

            # Calculate & attach missing cases for this stratum
            res_lvl$n_excluded_strata <- stratum_totals[lvl] - res_lvl$n
            res_lvl$n_excluded_strata <- max(0L, as.integer(res_lvl$n_excluded_strata))

            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              snp        = "",
              alleles    = "",
              group      = as.character(lvl),
              n          = res_lvl$n,
              missing    = if (res_lvl$n_excluded_strata > 0L) res_lvl$n_excluded_strata else '',
              maf        = round(res_lvl$maf, 4),
              genoCounts = res_lvl$genoCounts,
              hwePval    = res_lvl$hwePval
            ))
          }
        }
      }

      # ── Table-level note ─────────────────────────────────────────
      if (has_response || has_cov) {
        parts <- c()
        if (has_cov)      parts <- c(parts, "covariates")
        if (has_response) parts <- c(parts, "response")
        tbl$setNote(
          note = paste0("Complete cases used: rows missing any ",
                        paste(parts, collapse = " or "),
                        " or SNP value are excluded from all statistics."),
          key  = "missing_resp_cov")
      } else {
        tbl$setNote(note = NULL, key = "missing_resp_cov")
      }
    },


    # ── Allele frequencies ────────────────────────────────────────
    .fill_allele_freq = function(tbl, sm, snp_nm, response_raw, subpop,
                                  response_type, snp_raw,
                                  show_missing = FALSE, n_missing = 0L) {
      af <- sm$allele.freq


      # Stratify by binary response only
      do_strat <- isTRUE(subpop) && !is.null(response_raw) &&
                  identical(response_type, "binary")
      grp_levels <- character(0)
      if (do_strat) {
        grp_levels <- sort(unique(na.omit(as.character(response_raw))))
        if (length(grp_levels) < 2L || length(grp_levels) > 5L)
          do_strat <- FALSE
      }

      if (do_strat) {
        for (g in grp_levels) {
          safe <- gsub("[^A-Za-z0-9]", "_", g)
          tbl$addColumn(name   = paste0("count_", safe),
                        title  = paste0("N (", g, ")"),
                        type   = "integer")
          tbl$addColumn(name   = paste0("prop_", safe),
                        title  = paste0("% (", g, ")"),
                        type   = "number",
                        format = "dp=1")
        }
      }

      # n_nonmiss is used as denominator for % (already excludes NA since
      # snp_raw_cc contains no NAs, but stated explicitly for clarity)
      n_nonmiss <- sum(!is.na(snp_raw)) * 2L   # total allele count (diploid)

      allele_names <- rownames(af)
      for (i in seq_len(sm$nallele)) {
        al <- allele_names[i]
        row_vals <- list(
          allele = al,
          count  = as.integer(af[i, "Count"]),
          prop   = round(af[i, "Proportion"]*100, 1)
        )
        if (do_strat) {
          for (g in grp_levels) {
            safe   <- gsub("[^A-Za-z0-9]", "_", g)
            mask_g <- as.character(response_raw) == g
            # split every diploid genotype into its two alleles
            # (inputs are already complete-case so no NA guards needed)
            all_alleles <- unlist(strsplit(as.character(snp_raw)[mask_g], "/"))
            n_al  <- sum(all_alleles == al, na.rm = TRUE)
            n_tot <- length(all_alleles)
            row_vals[[paste0("count_", safe)]] <- as.integer(n_al)
            row_vals[[paste0("prop_",  safe)]] <-
              if (n_tot > 0L) round(n_al / n_tot * 100, 1) else NA_real_
          }
        }
        tbl$addRow(rowKey = al, values = row_vals)
      }

      # ── Missing row ──────────────────────────────────────────────
      if (isTRUE(show_missing) && sum(n_missing) > 0L) {
        miss_vals <- list(
          allele   = "Missing",
          count    = as.integer(sum(n_missing)),
          prop     = ''   # % not applicable for missing row
        )
        if (do_strat) {
          for (g in grp_levels) {
            safe <- gsub("[^A-Za-z0-9]", "_", g)
            miss_vals[[paste0("count_", safe)]] <- n_missing[g]
            miss_vals[[paste0("prop_",  safe)]] <- ''
          }
        }
        tbl$addRow(rowKey = "missing", values = miss_vals)
      }
    },

# ── Genotype frequencies ──────────────────────────────────────
    .fill_geno_freq = function(tbl, sm, ref, snp_raw, response,
                                response_type, subpop, response_raw,
                                show_missing = FALSE, n_missing = 0L) {
      

      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
        
        # Fix: Ensure response is numeric
        if (!is.numeric(response)) {
          response <- as.numeric(as.character(response))
        }
      }
      
      gf <- sm$genotype.freq
      gf <- tryCatch(reorder_geno(gf, ref), error = function(e) gf)

      # Stratify by binary response only
      do_strat <- isTRUE(subpop) && !is.null(response_raw) &&
                  identical(response_type, "binary")
      grp_levels <- character(0)
      if (do_strat) {
        grp_levels <- sort(unique(na.omit(as.character(response_raw))))
        if (length(grp_levels) < 2L || length(grp_levels) > 5L)
          do_strat <- FALSE
      }

      if (do_strat) {
        for (g in grp_levels) {
          safe <- gsub("[^A-Za-z0-9]", "_", g)
          tbl$addColumn(name   = paste0("count_", safe),
                        title  = paste0("N (", g, ")"),
                        type   = "integer")
          tbl$addColumn(name   = paste0("prop_", safe),
                        title  = paste0("% (", g, ")"),
                        type   = "number",
                        format = "dp=1")
        }
      }

      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
      }

      # n_nonmiss is the denominator for % (snp_raw_cc has no NAs)
      n_nonmiss <- sum(!is.na(snp_raw))

      snp_char <- as.character(snp_raw)
      for (i in seq_len(nrow(gf))) {
        geno_name <- rownames(gf)[i]
        if (geno_name == "NA") next

        cnt  <- as.integer(gf[i, "Count"])
        # Proportion from genetics::summary is already computed on non-missing
        prop <- if (is.na(gf[i, "Proportion"])) NA_real_
                else round(gf[i, "Proportion"]*100, 1)

        resp_stat <- ""
        if (response_type == "quantitative") {
          mask <- (snp_char == geno_name) & !is.na(response)
          # FIX: Added na.rm = TRUE to prevent crash on missing SNP values
          if (sum(mask, na.rm = TRUE) > 0) {
            mn <- mean(response[mask], na.rm = TRUE)
            se <- sd(response[mask],   na.rm = TRUE) / sqrt(sum(mask, na.rm = TRUE))
            resp_stat <- sprintf("%.2f (%.2f)", mn, se)
          }
        }

        row_vals <- list(
          genotype     = geno_name,
          count        = cnt,
          prop         = prop,
          responseStat = resp_stat
        )

        if (do_strat) {
          for (g in grp_levels) {
            safe   <- gsub("[^A-Za-z0-9]", "_", g)
            mask_g <- as.character(response_raw) == g
            # inputs are already complete-case; snp_char has no NAs
            n_g    <- sum(mask_g & snp_char == geno_name)
            n_tot  <- sum(mask_g)
            row_vals[[paste0("count_", safe)]] <- as.integer(n_g)
            row_vals[[paste0("prop_",  safe)]] <-
              if (n_tot > 0L) round(n_g / n_tot * 100, 1) else NA_real_
          }
        }
        tbl$addRow(rowKey = geno_name, values = row_vals)
      }

      # ── Missing row ──────────────────────────────────────────────
      if (isTRUE(show_missing) && sum(n_missing) > 0L) {
        miss_vals <- list(
          genotype     = "Missing",
          count        = as.integer(sum(n_missing)),
          prop         = '',
          responseStat = ""
        )
        if (do_strat) {
          for (g in grp_levels) {
            safe <- gsub("[^A-Za-z0-9]", "_", g)
            miss_vals[[paste0("count_", safe)]] <-  n_missing[g]
            miss_vals[[paste0("prop_",  safe)]] <- ''
          }
        }
        tbl$addRow(rowKey = "missing", values = miss_vals)
      }
    },
    
    # ── Hardy-Weinberg test ───────────────────────────────────────
    .fill_hwe = function(tbl, geno_obj, snp_nm, response_raw, subpop,
                         show_missing = FALSE, n_missing = 0L) {
      # Show/hide missing column
      tbl$getColumn("missing")$setVisible(isTRUE(show_missing))

      # Overall
      hw <- tryCatch(genetics::HWE.exact(geno_obj), error = function(e) NULL)
      if (is.null(hw)) return()

      st <- hw$statistic
      pr <- hw$parameter
      tbl$addRow(rowKey = "All", values = list(
        group   = "All subjects",
        n11     = as.integer(st["N11"]),
        n12     = as.integer(st["N12"]),
        n22     = as.integer(st["N22"]),
        missing = if (isTRUE(show_missing)) as.integer( sum(n_missing)) else '',
        pval    = hw$p.value
      ))

      # Stratified by response
      if (subpop && !is.null(response_raw)) {
        lvls <- unique(as.character(response_raw))   # no NAs: inputs are cc
        if (length(lvls) <= 5) {
          for (lvl in lvls) {
            mask <- as.character(response_raw) == lvl
            hw_sub <- tryCatch(
              genetics::HWE.exact(geno_obj[mask]),
              error = function(e) NULL)
            if (is.null(hw_sub)) next
            st2 <- hw_sub$statistic
            pr2 <- hw_sub$parameter
            tbl$addRow(rowKey = lvl, values = list(
              group   = paste0("Response = ", lvl),
              n11     = as.integer(st2["N11"]),
              n12     = as.integer(st2["N12"]),
              n22     = as.integer(st2["N22"]),
              missing = if (isTRUE(show_missing)) as.integer(n_missing[lvl]) else '',
              pval    = hw_sub$p.value
            ))
          }
        }
      }
    },

    # ── Linkage disequilibrium ────────────────────────────────────
    .run_ld = function(geno_list, opts, run_ldAnalysis, run_ldMatrix, run_ldPlot) {
      nms   <- names(geno_list)
      n     <- length(nms)
      pairs <- combn(nms, 2, simplify = FALSE)

      # Compute all pairwise LD results once
      ld_store <- list()
      for (pair in pairs) {
        key    <- paste(pair, collapse = "___")
        ld_res <- tryCatch(genetics::LD(geno_list[[pair[1]]], geno_list[[pair[2]]]),
                          error = function(e) NULL)
        if (!is.null(ld_res)) ld_store[[key]] <- ld_res
      }

      # Pairwise table
      if (run_ldAnalysis) {
        tbl <- self$results$ldGroup$ldTable
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          tbl$addRow(rowKey = paste(pair, collapse = "_"), values = list(
            snp1   = pair[1],
            snp2   = pair[2],
            r2     = round(ld_res$`r`^2,  3),
            Dprime = round(ld_res$`D'`,   3),
            D      = round(ld_res$`D`,    3),
            pval   = ld_res$`P-value`
          ))
        }
      }

      # ── LD matrix ─────────────────────────────────────────────
      if (run_ldMatrix) {
        mtbl   <- self$results$ldGroup$ldMatrixTable
        metric <- opts$ldMetric   # "r2", "Dprime", or "D"

        # Add one column per SNP (beyond the row-label column already defined)
        for (nm in nms) {
          safe_nm <- gsub("[^A-Za-z0-9_]", "_", nm)
          mtbl$addColumn(name = safe_nm, title = nm, type = "text")
        }

        # Build n×n value matrices
        upper_mat <- matrix("", n, n, dimnames = list(nms, nms))
        lower_mat <- matrix("", n, n, dimnames = list(nms, nms))
        diag(upper_mat) <- nms
        diag(lower_mat) <- nms

        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next

          p_val <- ld_res$`P-value`
          p_str <- if (!is.na(p_val)) {
            if (p_val < 0.001) "< .001" else sprintf("%.3f", p_val)
          } else ""

          up_val <- switch(metric,
            Dprime = sprintf("%.3f", round(ld_res$`D'`, 3)),
            r2      = sprintf("%.3f", round(ld_res$`r`^2,   3)),
            D      = sprintf("%.3f", round(ld_res$`D`,   3))
          )

          upper_mat[pair[1], pair[2]] <- up_val
          upper_mat[pair[2], pair[1]] <- ""          # will be lower
          lower_mat[pair[1], pair[2]] <- ""
          lower_mat[pair[2], pair[1]] <- p_str
        }

        for (i in seq_len(n)) {
          row_vals <- list(snp = nms[i])
          for (j in seq_len(n)) {
            safe_nm <- gsub("[^A-Za-z0-9_]", "_", nms[j])
            if (i == j) {
              row_vals[[safe_nm]] <- nms[i]
            } else if (j > i) {
              row_vals[[safe_nm]] <- upper_mat[i, j]
            } else {
              row_vals[[safe_nm]] <- lower_mat[i, j]
            }
          }
          mtbl$addRow(rowKey = paste0("row_", i), values = row_vals)
        }
        # Add footnote explaining the layout
        metric_label <- switch(metric,
          Dprime = "D\'", r2 = "r²", D = "D")
        mtbl$setNote(
          key  = "layout",
          note = paste0("Upper triangle: ", metric_label,
                        ". Lower triangle: P-value. Diagonal: SNP name."))
      }

      # ── Store LD data for heatmap plot ──────────────────────────
      if (run_ldPlot) {
        private$.ld_store  <- ld_store
        private$.ld_nms    <- nms
        private$.ld_metric <- opts$ldMetric
        self$results$ldGroup$ldPlotImage$setState(list(
          ld_store = ld_store, nms = nms, metric = opts$ldMetric))
      }
    },

    # ── LD heatmap render ─────────────────────────────────────────
    .render_ld_plot = function(image, ggtheme, theme, ...) {
      state <- image$state
      if (is.null(state)) return(FALSE)

      ld_store <- state$ld_store
      nms      <- state$nms
      metric   <- state$metric
      n        <- length(nms)

      # Build a data frame for ggplot (full symmetric matrix for heatmap)
      metric_label <- switch(metric, Dprime = "D'", r2 = "r²", D = "D")
      df_rows <- list()
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i == j) {
            val <- 1.0
          } else {
            key <- if (i < j) paste(c(nms[i], nms[j]), collapse = "___")
                   else        paste(c(nms[j], nms[i]), collapse = "___")
            ld_res <- ld_store[[key]]
            val <- if (!is.null(ld_res)) {
              switch(metric,
                Dprime = abs(as.numeric(ld_res$`D'`)),
                r2     = as.numeric(ld_res$`r`)^2,
                D      = abs(as.numeric(ld_res$`D`))
              )
            } else NA_real_
          }
          df_rows[[length(df_rows) + 1L]] <- data.frame(
            SNP1  = factor(nms[i], levels = rev(nms)),
            SNP2  = factor(nms[j], levels = nms),
            value = val,
            stringsAsFactors = FALSE
          )
        }
      }
      df <- do.call(rbind, df_rows)

      # Value string for lower-triangle annotation (p-values)
      p_mat <- matrix(NA_real_, n, n, dimnames = list(nms, nms))
      for (pair_key in names(ld_store)) {
        parts <- strsplit(pair_key, "___")[[1]]
        pv    <- ld_store[[pair_key]]$`P-value`
        p_mat[parts[1], parts[2]] <- pv
        p_mat[parts[2], parts[1]] <- pv
      }

      df$label <- ""
      for (k in seq_len(nrow(df))) {
        i_nm <- as.character(df$SNP1[k])
        j_nm <- as.character(df$SNP2[k])
        i_idx <- which(nms == i_nm)
        j_idx <- which(nms == j_nm)
        if (i_idx > j_idx) {   # lower triangle → p-value
          pv <- p_mat[i_nm, j_nm]
          df$label[k] <- if (!is.na(pv)) {
            if (pv < 0.001) "<.001" else sprintf("%.3f", pv)
          } else ""
        } else if (i_idx < j_idx) {  # upper triangle → metric value
          key <- paste(c(nms[min(i_idx, j_idx)], nms[max(i_idx, j_idx)]),
                       collapse = "___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) {
            raw <- switch(metric,
              r2      = ld_res$`r`^2,
              Dprime = ld_res$`D'`,
              D      = ld_res$`D`
            )
            df$label[k] <- sprintf("%.3f", round(as.numeric(raw), 3))
          }
        } else {
          df$label[k] <- i_nm   # diagonal
        }
      }

      colour_label <- switch(metric,
        Dprime = "|D'|", r2 = "r²", D = "|D|")

      # Use scale_fill_gradientn with explicit continuous colours to avoid
      # jamovi theme interfering with the fill scale via ggPalette()
      p <- ggplot2::ggplot(df, ggplot2::aes(x = SNP2, y = SNP1, fill = value)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = label),
                           size = 3, colour = "grey10") +
        ggplot2::scale_fill_gradientn(
          colours  = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"),
          limits   = c(0, 1),
          na.value = "grey85",
          name     = colour_label
        ) +
        ggplot2::scale_x_discrete(position = "bottom") +
        ggplot2::scale_y_discrete() +
        ggplot2::labs(
          title = paste0("LD Heatmap  •  upper: ", metric_label,
                         "  |  lower: p-value"),
          x = NULL, y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y     = ggplot2::element_text(hjust = 1),
          panel.grid      = ggplot2::element_blank(),
          legend.position = "right",
          plot.title      = ggplot2::element_text(size = 11, face = "bold",
                                                  margin = ggplot2::margin(b = 8))
        )

      print(p)
      TRUE
    },

    # ── SNP association ───────────────────────────────────────────
    .fill_assoc = function(tbl, snp_raw, ref, response, cov_df,
                        response_type, opts, run_snpAssoc) {

      if (!run_snpAssoc) return()
      
      # Validate response for binary models
      if (response_type == "binary") {
        # Check if response has exactly 2 levels after removing NAs
        resp_clean <- response[!is.na(response)]
        if (length(unique(resp_clean)) != 2) {
          tbl$setNote(
            key = "response_error",
            note = "Binary response requires exactly 2 categories. ",
            "Consider using quantitative response type or check your data.")
          return()
        }
      }

      # ── Dynamic column header ─────────────────────────
      effect_col <- tbl$getColumn("effect")
      if (response_type == "binary") {
        effect_col$setTitle("OR")
      } else {
        effect_col$setTitle("\u03B2")
      }

      # ── Add covariate note ─────────────────────────────
      note_key <- "covariates"
      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        cov_names <- names(cov_df)
        cov_names <- sapply(cov_names, function(x) {
          if (!is.null(self$data[[x]])) {
            attr(self$data[[x]], "label") %||% x
          } else x
        })
        note_txt <- paste0("Model adjusted for: ", paste(cov_names, collapse = ", "))
        tbl$setNote(note = note_txt, key = note_key)
      } else {
        tbl$setNote(note = NULL, key = note_key)
      }

      # ── Missing covariate note ─────────────────────────
      {
        snp_enc_tmp <- encode_model(as.character(snp_raw), ref, "logadditive")
        complete_full <- !is.na(response) & !is.na(snp_enc_tmp)
        if (!is.null(cov_df) && ncol(cov_df) > 0)
          complete_full <- complete_full & complete.cases(cov_df)
        n_miss <- length(response) - sum(complete_full)
        if (n_miss > 0)
          tbl$setNote(
            note = paste0("Note: ", n_miss,
                          " observation(s) excluded due to missing values."),
            key  = "missing_cov")
        else
          tbl$setNote(note = NULL, key = "missing_cov")
      }
          models <- c()
      if (opts$modelCodominant)   models <- c(models, "codominant")
      if (opts$modelDominant)     models <- c(models, "dominant")
      if (opts$modelRecessive)    models <- c(models, "recessive")
      if (opts$modelOverdominant) models <- c(models, "overdominant")
      if (opts$modelLogAdditive)  models <- c(models, "logadditive")

      model_labels <- c(
        codominant   = "Codominant",
        dominant     = "Dominant",
        recessive    = "Recessive",
        overdominant = "Overdominant",
        logadditive  = "Log-additive"
      )

      note_key <- "lrt"
      lrt_note <- NULL
      
      row_key <- 0L
      for (mdl in models) {
        snp_enc <- encode_model(as.character(snp_raw), ref, mdl)
        res_list <- fit_model(snp_enc, response, cov_df, mdl,
                               response_type, opts$ciWidth)
        if (is.null(res_list)) next

        if (mdl == "codominant" && !is.null(res_list) && length(res_list) > 0) {
          gp <- res_list[[1]]$global_p
          if (!is.na(gp)) {
            lrt_note <- paste0(
              "Codominant model: likelihood ratio test P = ",
              format.pval(gp, digits = 3)
            )
          }
        }

        # ── Set / remove LRT note ──────────────────────────
        if (!is.null(lrt_note)) {
          tbl$setNote(note = lrt_note, key = note_key)
        } else {
          tbl$setNote(note = NULL, key = note_key)
        }

        first_row <- TRUE
        for (res in res_list) {
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model      = if (first_row) model_labels[mdl] else "",
            comparison = res$comparison,
            effect     = res$effect,
            ciLow      = res$ci_low,
            ciHigh     = res$ci_high,
            pval       = res$pval,
            AIC = if (first_row && !is.nan(res$aic)) res$aic else ""
          ))
          first_row <- FALSE
        }
      }
    },

    # ── SNP × covariate interaction ───────────────────────────────
    .fill_interaction = function(tbl, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts) {

      effect_col <- tbl$getColumn("effect")
      if (response_type == "binary") {
        effect_col$setTitle("OR")
      } else {
        effect_col$setTitle("\u03B2")
      }

      # Note: interaction variable + remaining adjusters
      adj_vars <- setdiff(names(cov_df), interaction_var)
      note_parts <- paste0("Interaction covariate: ", interaction_var)
      if (length(adj_vars) > 0)
        note_parts <- paste0(note_parts, ". Adjusted for: ",
                             paste(adj_vars, collapse = ", "))
      tbl$setNote(note = note_parts, key = "intcov")

      # ── Missing note (SNP + response + covariates) ────────────────
      {
        snp_enc_tmp <- encode_model(as.character(snp_raw), ref, "logadditive")
        complete_full <- !is.na(response) & !is.na(snp_enc_tmp) & complete.cases(cov_df)
        n_miss <- length(response) - sum(complete_full)
        if (n_miss > 0)
          tbl$setNote(
            note = paste0("Note: ", n_miss,
                          " observation(s) excluded due to missing values."),
            key  = "missing_cov")
        else
          tbl$setNote(note = NULL, key = "missing_cov")
      }

      models <- c()
      if (opts$modelCodominant)   models <- c(models, "codominant")
      if (opts$modelDominant)     models <- c(models, "dominant")
      if (opts$modelRecessive)    models <- c(models, "recessive")
      if (opts$modelOverdominant) models <- c(models, "overdominant")
      if (opts$modelLogAdditive)  models <- c(models, "logadditive")

      model_labels <- c(
        codominant   = "Codominant",
        dominant     = "Dominant",
        recessive    = "Recessive",
        overdominant = "Overdominant",
        logadditive  = "Log-additive"
      )

      row_key <- 0L
      for (mdl in models) {
        snp_enc  <- encode_model(as.character(snp_raw), ref, mdl)
        res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                          interaction_var, mdl,
                                          response_type, opts$ciWidth)
        if (is.null(res_list)) next

        first_row      <- TRUE
        first_inter    <- TRUE   # pvalInteraction shown only on first interaction row
        for (res in res_list) {
          row_key <- row_key + 1L
          is_inter <- !is.na(res$pval_interaction)
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model           = if (first_row) model_labels[mdl] else "",
            term            = res$term,
            effect          = res$effect,
            ciLow           = res$ci_low,
            ciHigh          = res$ci_high,
            pval            = res$pval,
            pvalInteraction = if (is_inter && first_inter) res$pval_interaction else "",
            AIC             = if (first_row && !is.nan(res$aic)) res$aic else ""
          ))
          first_row <- FALSE
          if (is_inter) first_inter <- FALSE
        }
      }
    },

    .run_haplo = function(geno_list, data, response, response_raw, response_type, cov_df, 
                       opts, run_haploFreq, run_haploAssoc, run_haploInteraction,
                       run_subpop = FALSE, complete_mask = NULL) {
            
      # ── Inline helper (hoisted so haploAssoc and haploInteraction can share) ──
      na.geno.keep <- function(m) {
        mf.gindx <- function(m) {
          nvars    <- length(m)
          typevars <- rep(0, nvars)
          for (i in seq_len(nvars)) typevars[i] <- data.class(m[[i]])
          gindx <- seq_len(nvars)[typevars == "model.matrix" | typevars == "matrix"]
          if (length(gindx) == 0) stop("No geno matrix in data frame")
          if (length(gindx) >  1) stop("More than 1 geno matrix in data frame")
          gindx
        }
        gindx    <- mf.gindx(m)
        yxmiss   <- apply(is.na(m[, -gindx, drop = FALSE]), 1, any)
        gmiss    <- apply(is.na(m[,  gindx, drop = FALSE]), 1, all)
        genoAttr <- attributes(m[, gindx])
        allmiss  <- yxmiss | gmiss
        m        <- m[!allmiss, ]
        genoAttr$dim[1] <- genoAttr$dim[1] - sum(allmiss)
        nloc <- ncol(m[, gindx]) / 2
        for (k in seq_len(nloc)) {
          ualleles <- unique(c(m[, gindx][, 2*k-1], m[, gindx][, 2*k]))
          nalleles <- length(genoAttr$unique.alleles[[k]])
          if (length(ualleles) < nalleles)
            genoAttr$unique.alleles[[k]] <-
              genoAttr$unique.alleles[[k]][!is.na(match(seq_len(nalleles), ualleles))]
        }
        for (att in names(genoAttr)) attr(m[, gindx], att) <- genoAttr[[att]]
        attr(m, "yxmiss") <- yxmiss
        attr(m, "gmiss")  <- gmiss
        m
      }
      # ── Inline helper end ───────────────────────────────────────────────

      snp_names <- names(geno_list)
      allele_list <- lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]]))
      allele_mat  <- do.call(cbind, allele_list)
      
      geno_setup <- tryCatch(
        haplo.stats::setupGeno(allele_mat, locus.label = snp_names),
        error = function(e) NULL
      )
      if (is.null(geno_setup)) return()
      
      u_alleles <- attr(geno_setup, "unique.alleles")

      decode_haplo_row <- function(codes, label_list) {
        # Check if this is the "rare" indicator from haplo.df
        if (all(codes == "*") || any(codes == "*")) return("Rare (combined)")
        
        parts <- character(length(codes))
        for (i in seq_along(codes)) {
          idx <- as.numeric(codes[i])
          parts[i] <- if (!is.na(idx)) label_list[[i]][idx] else "?"
        }
        paste(parts, collapse = "-")
      }

      # ── Missing counts and complete-case keep mask ──────────────
      n_total_rows <- nrow(allele_mat)

      # A row is SNP-missing if ANY SNP is NA (not just all)
      snp_miss_mask <- apply(is.na(allele_mat), 1, any)
      n_snp_miss    <- sum(snp_miss_mask)

      # Response/covariate missings among rows that have at least one valid SNP
      resp_cov_miss_mask <- rep(FALSE, n_total_rows)
      if (!is.null(response))
        resp_cov_miss_mask <- resp_cov_miss_mask | is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0)
        resp_cov_miss_mask <- resp_cov_miss_mask | !complete.cases(cov_df)
      n_resp_cov_miss <- sum(resp_cov_miss_mask & !snp_miss_mask)

      # Use complete_mask passed from .run() when available; otherwise rebuild.
      # In both cases also exclude rows with any SNP missing.
      if (!is.null(complete_mask)) {
        keep <- complete_mask & !snp_miss_mask
      } else {
        keep <- if (!is.null(response)) !is.na(response) else rep(TRUE, n_total_rows)
        if (!is.null(cov_df)) keep <- keep & complete.cases(cov_df)
        keep <- keep & !snp_miss_mask
      }

      # subset_geno: subset a setupGeno matrix by row while preserving all
      # attributes (unique.alleles, locus.label, etc.) that [.matrix strips.
      subset_geno <- function(gs, idx) {
        saved_attr <- attributes(gs)
        gs2        <- gs[idx, , drop = FALSE]
        # Restore every attribute except dim/dimnames which [.matrix sets correctly
        for (att in setdiff(names(saved_attr), c("dim", "dimnames"))) {
          attr(gs2, att) <- saved_attr[[att]]
        }
        gs2
      }

      # Helper to fill haplotype freq table for a given row-mask
      fill_haplo_freq_tbl <- function(tbl, row_mask, row_key_prefix = "") {
        em_sub <- tryCatch(
          haplo.stats::haplo.em(subset_geno(geno_setup, row_mask), locus.label = snp_names),
          error = function(e) NULL
        )
        if (is.null(em_sub)) return(invisible(NULL))
        freqs    <- em_sub$hap.prob
        rare_sum <- 0
        for (i in seq_along(freqs)) {
          if (freqs[i] < opts$haploFreqMin) {
            rare_sum <- rare_sum + freqs[i]
            next
          }
          label <- decode_haplo_row(as.numeric(em_sub$haplotype[i, ]), u_alleles)
          tbl$addRow(rowKey = paste0(row_key_prefix, "f", i), values = list(
            haplotype = label,
            freq      = round(freqs[i], 4)
          ))
        }
        if (rare_sum > 0) {
          tbl$addRow(rowKey = paste0(row_key_prefix, "rare"), values = list(
            haplotype = paste0("Rare (<", opts$haploFreqMin, ")"),
            freq      = round(rare_sum, 4)
          ))
        }
        invisible(em_sub)
      }

      # 1. Haplotype Frequencies (EM)
      if (run_haploFreq) {
        tbl <- self$results$haploGroup$haploFreqTable

        # ── Determine if we should stratify ──────────────────────────
        do_strat_haplo <- isTRUE(run_subpop) &&
                          !is.null(response_raw) &&
                          identical(response_type, "binary")
        grp_levels_haplo <- character(0)
        if (do_strat_haplo) {
          grp_levels_haplo <- sort(unique(na.omit(as.character(response_raw[keep]))))
          if (length(grp_levels_haplo) != 2L) do_strat_haplo <- FALSE
        }

        if (do_strat_haplo) {
          # Add one freq column per group (2-column layout, not extra rows)
          tbl$addColumn(name  = "freq_g0",
                        title = as.character(grp_levels_haplo[1]),
                        type  = "number")
          tbl$addColumn(name  = "freq_g1",
                        title = as.character(grp_levels_haplo[2]),
                        type  = "number")
        }

        # Overall EM on the full complete-case sample (keep mask)
        em_all <- tryCatch(
          haplo.stats::haplo.em(subset_geno(geno_setup, keep), locus.label = snp_names),
          error = function(e) NULL
        )

        if (!is.null(em_all)) {
          freqs    <- em_all$hap.prob
          rare_sum <- 0

          # Build per-group EM results up-front so we can fill same rows
          em_grp <- list()
          if (do_strat_haplo) {
            for (lvl in grp_levels_haplo) {
              # keep_lvl derives from keep, so response + covariates + SNPs are all complete
              keep_lvl <- keep & as.character(response_raw) == lvl
              if (sum(keep_lvl) < 5) next
              em_grp[[lvl]] <- tryCatch(
                haplo.stats::haplo.em(subset_geno(geno_setup, keep_lvl),
                                      locus.label = snp_names),
                error = function(e) NULL
              )
            }
            # Build label -> freq lookup per group
            grp_freq <- lapply(em_grp, function(em_g) {
              if (is.null(em_g)) return(list())
              setNames(
                as.list(round(em_g$hap.prob, 4)),
                sapply(seq_len(nrow(em_g$haplotype)), function(j)
                  decode_haplo_row(as.numeric(em_g$haplotype[j, ]), u_alleles))
              )
            })
          }

          for (i in seq_along(freqs)) {
            if (freqs[i] < opts$haploFreqMin) {
              rare_sum <- rare_sum + freqs[i]
              next
            }
            label    <- decode_haplo_row(as.numeric(em_all$haplotype[i, ]), u_alleles)
            row_vals <- list(haplotype = label, freq = round(freqs[i], 4))
            if (do_strat_haplo) {
              row_vals$freq_g0 <- grp_freq[[grp_levels_haplo[1]]][[label]] %||% NA_real_
              row_vals$freq_g1 <- grp_freq[[grp_levels_haplo[2]]][[label]] %||% NA_real_
            }
            tbl$addRow(rowKey = paste0("f", i), values = row_vals)
          }
          if (rare_sum > 0) {
            row_vals <- list(
              haplotype = paste0("Rare (<", opts$haploFreqMin, ")"),
              freq      = round(rare_sum, 4)
            )
            if (do_strat_haplo) {
              em0     <- em_grp[[grp_levels_haplo[1]]]
              em1     <- em_grp[[grp_levels_haplo[2]]]
              rare_g0 <- if (!is.null(em0))
                           round(sum(em0$hap.prob[em0$hap.prob < opts$haploFreqMin]), 4)
                         else NA_real_
              rare_g1 <- if (!is.null(em1))
                           round(sum(em1$hap.prob[em1$hap.prob < opts$haploFreqMin]), 4)
                         else NA_real_
              row_vals$freq_g0 <- if (!is.na(rare_g0) && rare_g0 > 0) rare_g0 else NA_real_
              row_vals$freq_g1 <- if (!is.na(rare_g1) && rare_g1 > 0) rare_g1 else NA_real_
            }
            tbl$addRow(rowKey = "rare_freq", values = row_vals)
          }
        }

        # ── Missing notes ──────────────────────────────────────────
        if (n_snp_miss > 0)
          tbl$setNote(
            note = paste0(n_snp_miss, " observation(s) with missing SNP genotype(s) excluded."),
            key  = "missing_snp")
        else
          tbl$setNote(note = NULL, key = "missing_snp")

        if (n_resp_cov_miss > 0) {
          parts <- c()
          if (!is.null(cov_df) && ncol(cov_df) > 0) parts <- c(parts, "covariate")
          if (!is.null(response))                    parts <- c(parts, "response")
          tbl$setNote(
            note = paste0(n_resp_cov_miss, " observation(s) with missing ",
                          paste(parts, collapse = "/"), " values excluded."),
            key  = "missing_resp_cov")
        } else {
          tbl$setNote(note = NULL, key = "missing_resp_cov")
        }

      }  # end if (run_haploFreq)

    # ── Haplotype association ─────────────────────────────────────
    if (run_haploAssoc && !is.null(response)) {

      family     <- if (response_type == "binary") "binomial" else "gaussian"
      y_sub      <- if (response_type == "binary") {
        as.numeric(as.factor(response[keep])) - 1L
      } else {
        response[keep]
      }

      m_model      <- data.frame(y = y_sub)
      m_model$geno <- subset_geno(geno_setup, keep)
      if (!is.null(cov_df)) {
        m_model    <- cbind(m_model, cov_df[keep, , drop = FALSE])
        formula_str <- paste("y ~ geno +", paste(names(cov_df), collapse = " + "))
      } else {
        formula_str <- "y ~ geno"
      }

      haplo_fit <- tryCatch(
        haplo.stats::haplo.glm(
          as.formula(formula_str),
          family    = family,
          data      = m_model,
          na.action = na.geno.keep,
          control   = haplo.stats::haplo.glm.control(
                        haplo.freq.min = opts$haploFreqMin)
        ),
        error = function(e) {
          self$results$validationMsg$setContent(
            paste0("<b>Haplotype GLM error:</b> ", e$message))
          NULL
        }
      )

      if (!is.null(haplo_fit)) {
        tbl <- self$results$haploGroup$haploAssocTable

        # Set effect column title to match response type
        tbl$getColumn("effect")$setTitle(
          if (response_type == "binary") "OR" else "β")

        # ── Decode haplotype label from haplo.unique row ────────────
        # haplo.unique stores allele characters directly (e.g. "C", "T", "A"),
        # one per locus — confirmed from diagnostic output.
        label_from_unique_row <- function(row_vec) {
          paste(as.character(row_vec), collapse = "-")
        }

        # ── Pull coefficients and CIs ───────────────────────────────
        # haplo.glm model matrix column names for haplotype terms are stored
        # in haplo_fit$haplo.names.  The actual rownames of coef() follow the
        # convention "geno" + haplo.names (the model frame column is "geno").
        # We match positionally rather than by name to be robust to separator
        # differences across haplo.stats versions.
        coef_sum <- tryCatch(summary(haplo_fit)$coefficients, error = function(e) NULL)
        ci_mat   <- tryCatch(confint(haplo_fit, level = opts$ciWidth / 100),
                             error = function(e) NULL)

        # summary(haplo_fit)$coefficients columns are: coef | se | t.stat | pval
        # (haplo.glm uses its own summary method, not summary.glm)
        # haplo_rows: positions of geno.* rows in coef_sum (intercept is row 1)
        haplo_rows <- if (!is.null(coef_sum)) {
          grep("^geno", rownames(coef_sum))
        } else integer(0)

        # Helper: get beta, se, pval, ci for haplo-term at position pos
        # (1-based within haplo_rows, matching order of haplo.common).
        get_stats <- function(pos) {
          row_idx <- if (!is.na(pos) && pos >= 1L && pos <= length(haplo_rows))
                       haplo_rows[pos] else NA_integer_
          if (is.na(row_idx) || is.null(coef_sum) ||
              row_idx < 1L || row_idx > nrow(coef_sum)) {
            return(list(beta = NA_real_, se = NA_real_, pval = NA_real_,
                        ci_lo = NA_real_, ci_hi = NA_real_))
          }
          rn   <- rownames(coef_sum)[row_idx]
          beta <- coef_sum[row_idx, "coef"]
          se   <- coef_sum[row_idx, "se"]
          pval <- coef_sum[row_idx, "pval"]
          # CI from ci_mat if available, otherwise Wald ± z * se
          if (!is.null(ci_mat) && rn %in% rownames(ci_mat)) {
            ci_lo <- ci_mat[rn, 1]
            ci_hi <- ci_mat[rn, 2]
          } else {
            z     <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
            ci_lo <- beta - z * se
            ci_hi <- beta + z * se
          }
          list(beta = beta, se = se, pval = pval, ci_lo = ci_lo, ci_hi = ci_hi)
        }

        make_row <- function(label, freq, stats) {
          b  <- stats$beta
          lo <- stats$ci_lo
          hi <- stats$ci_hi
          list(
            haplotype = label,
            freq      = round(freq, 4),
            effect    = if (response_type == "binary") exp(b)  else b,
            ciLow     = if (response_type == "binary") exp(lo) else lo,
            ciHigh    = if (response_type == "binary") exp(hi) else hi,
            pval      = stats$pval
          )
        }

        # ── Reference haplotype (OR = 1 by definition) ────────────
        base_idx   <- haplo_fit$haplo.base
        base_label <- label_from_unique_row(haplo_fit$haplo.unique[base_idx, ])
        base_freq  <- haplo_fit$haplo.freq[base_idx]
        tbl$addRow(rowKey = "base", values = list(
          haplotype = paste0(base_label, " (Ref)"),
          freq      = round(base_freq, 4),
          effect    = if (response_type == "binary") 1.0 else 0.0,
          ciLow     = '',
          ciHigh    = '',
          pval      = ''
        ))

        # ── Common haplotypes ──────────────────────────────────────
        # haplo.common: integer index vector into haplo.unique rows,
        # in the same order as haplo.names / the GLM coefficients.
        common_idx <- haplo_fit$haplo.common
        for (j in seq_along(common_idx)) {
          h_idx   <- common_idx[j]
          h_label <- label_from_unique_row(haplo_fit$haplo.unique[h_idx, ])
          h_freq  <- haplo_fit$haplo.freq[h_idx]
          stats   <- get_stats(j)
          tbl$addRow(rowKey = paste0("h", j),
                     values = make_row(h_label, h_freq, stats))
        }

        # ── Rare combined term ─────────────────────────────────────
        has_rare <- isTRUE(haplo_fit$haplo.rare.term) ||
                    (length(haplo_fit$haplo.rare) > 0)
        if (has_rare) {
          rare_freq <- sum(haplo_fit$haplo.freq[haplo_fit$haplo.rare])
          # rare term is the last haplotype coefficient
          stats     <- get_stats(length(common_idx) + 1L)
          tbl$addRow(rowKey = "rare",
                     values = make_row(
                       paste0("Rare (<", opts$haploFreqMin, ")"),
                       rare_freq, stats))
        }

              # ── Add covariate note ─────────────────────────────
        note_key <- "covariates"
        if (!is.null(cov_df) && ncol(cov_df) > 0) {
          cov_names <- names(cov_df)
          cov_names <- sapply(cov_names, function(x) {
            if (!is.null(self$data[[x]])) {
              attr(self$data[[x]], "label") %||% x
            } else x
          })
          note_txt <- paste0("Model adjusted for: ", paste(cov_names, collapse = ", "))
          tbl$setNote(note = note_txt, key = note_key)
        } else {
          tbl$setNote(note = NULL, key = note_key)
        }

        # ── Missing notes (SNPs separate; response/cov combined) ───
        if (n_snp_miss > 0)
          tbl$setNote(
            note = paste0(n_snp_miss,
                          " observation(s) with missing SNP genotype(s) excluded."),
            key  = "missing_snp")
        else
          tbl$setNote(note = NULL, key = "missing_snp")

        if (n_resp_cov_miss > 0) {
          ha_parts <- c()
          if (!is.null(cov_df) && ncol(cov_df) > 0) ha_parts <- c(ha_parts, "covariate")
          ha_parts <- c(ha_parts, "response")
          tbl$setNote(
            note = paste0(n_resp_cov_miss, " observation(s) with missing ",
                          paste(ha_parts, collapse = "/"), " values excluded."),
            key  = "missing_resp_cov")
        } else {
          tbl$setNote(note = NULL, key = "missing_resp_cov")
        }

      }
    }

# ── Haplotype × covariate interaction ─────────────────────────
if (run_haploInteraction && !is.null(cov_df) && ncol(cov_df) >= 1 && !is.null(response)) {

  int_var  <- names(cov_df)[1]   # first covariate is always the interaction term
  adj_vars <- setdiff(names(cov_df), int_var)

  tbl_int <- self$results$haploGroup$haploInteractionTable
  tbl_int$getColumn("effect")$setTitle(
    if (response_type == "binary") "OR" else "\u03B2")

  note_parts <- paste0("Interaction covariate: ", int_var)
  if (length(adj_vars) > 0)
    note_parts <- paste0(note_parts, ". Adjusted for: ",
                         paste(adj_vars, collapse = ", "))
  tbl_int$setNote(note = note_parts, key = "intcov")

  # ── Missing notes (SNPs separate; response/cov combined) ──────
  if (n_snp_miss > 0)
    tbl_int$setNote(
      note = paste0(n_snp_miss,
                    " observation(s) with missing SNP genotype(s) excluded."),
      key  = "missing_snp")
  else
    tbl_int$setNote(note = NULL, key = "missing_snp")

  if (n_resp_cov_miss > 0) {
    hi_parts <- c("covariate", "response")
    tbl_int$setNote(
      note = paste0(n_resp_cov_miss, " observation(s) with missing ",
                    paste(hi_parts, collapse = "/"), " values excluded."),
      key  = "missing_resp_cov")
  } else {
    tbl_int$setNote(note = NULL, key = "missing_resp_cov")
  }

  family_int <- if (response_type == "binary") "binomial" else "gaussian"
  y_int      <- if (response_type == "binary") {
    as.numeric(as.factor(response[keep])) - 1L
  } else {
    response[keep]
  }

  m_int      <- data.frame(y = y_int)
  m_int$geno <- subset_geno(geno_setup, keep)
  if (!is.null(cov_df))
    m_int <- cbind(m_int, cov_df[keep, , drop = FALSE])

  adj_part <- if (length(adj_vars) > 0)
    paste("+", paste(adj_vars, collapse = "+")) else ""

  formula_int_str  <- paste("y ~ geno *", int_var, adj_part)
  formula_main_str <- paste("y ~ geno +", int_var, adj_part)

  haplo_int_fit <- tryCatch(
    haplo.stats::haplo.glm(
      as.formula(formula_int_str),
      family    = family_int,
      data      = m_int,
      na.action = na.geno.keep,
      control   = haplo.stats::haplo.glm.control(
                    haplo.freq.min = opts$haploFreqMin)
    ),
    error = function(e) {
      self$results$validationMsg$setContent(
        paste0("<b>Haplotype interaction GLM error:</b> ", e$message))
      NULL
    }
  )
  haplo_main_fit <- tryCatch(
    haplo.stats::haplo.glm(
      as.formula(formula_main_str),
      family    = family_int,
      data      = m_int,
      na.action = na.geno.keep,
      control   = haplo.stats::haplo.glm.control(
                    haplo.freq.min = opts$haploFreqMin)
    ),
    error = function(e) NULL
  )

  p_inter_haplo <- NA_real_   # initialised here so note is always settable

  if (!is.null(haplo_int_fit) && !is.null(haplo_main_fit)) {

    # ── LRT for the overall interaction ───────────────────────
    # anova.haplo.glm fails when the two models have different effective sample
    # sizes (different EM convergence sets). Compute the LRT directly from the
    # deviances and residual df stored in each fit object instead.
    dev_diff <- haplo_main_fit$deviance - haplo_int_fit$deviance
    df_diff  <- haplo_main_fit$df.residual - haplo_int_fit$df.residual
    p_inter_haplo <- if (!is.na(dev_diff) && !is.na(df_diff) && df_diff > 0)
      pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
    else NA_real_

    coef_sum_int <- tryCatch(summary(haplo_int_fit)$coefficients,
                              error = function(e) NULL)
    ci_int       <- tryCatch(confint(haplo_int_fit, level = opts$ciWidth / 100),
                              error = function(e) NULL)

    if (!is.null(coef_sum_int)) {
      all_rows_int <- rownames(coef_sum_int)

      # ── Build coef-name → display-label map ───────────────────
      # haplo.glm names interaction-model coefficients as geno.1, geno.2, ...
      # (and geno.1:int_var, geno.2:int_var for interaction terms).
      # The numeric suffix is a 1-based position into haplo.common (excluding
      # the base haplotype), with an optional trailing "rare" term.
      # We reconstruct allele labels the same way the association table does:
      # positionally via haplo.common → haplo.unique.
      decode_haplo_label <- function(row_vec) {
        paste(as.character(row_vec), collapse = "-")
      }

      rare_label <- paste0("Rare (<", opts$haploFreqMin, ")")

      # haplo.names holds the full coefficient names, e.g. "geno.5", "geno.6",
      # "geno.rare".  The numeric suffix IS the row index into haplo.unique —
      # confirmed from debug output.  Build the map from every geno main-effect
      # row to its allele label, then add the matching interaction row.
      geno_main_rows <- grep("^geno[^:]+$", all_rows_int, value = TRUE)
      raw_to_label   <- character(0)

      for (rn in geno_main_rows) {
        suffix <- sub("^geno\\.", "", rn)   # "5", "6", "rare", …

        display_label <- if (grepl("^[0-9]+$", suffix)) {
          idx <- as.integer(suffix)
          if (!is.na(idx) && idx >= 1L && idx <= nrow(haplo_int_fit$haplo.unique)) {
            decode_haplo_label(haplo_int_fit$haplo.unique[idx, ])
          } else paste0("Haplotype ", suffix)
        } else if (grepl("rare", suffix, ignore.case = TRUE)) {
          rare_label
        } else {
          suffix   # already an allele string in older haplo.stats builds
        }

        raw_to_label[rn] <- display_label
        # Interaction rows are named "geno.5:SEXMale" — the covariate level is
        # appended, so we can't construct the exact name.  Match by prefix instead.
        inter_rns <- grep(paste0("^", rn, ":"), all_rows_int, value = TRUE)
        for (irn in inter_rns)
          raw_to_label[irn] <- paste0(display_label, " \u00D7 ", sub(paste0("^", rn, ":"), "", irn))
      }

      # ── Identify main-effect and interaction rows ──────────
      main_rows  <- grep("^geno[^:]+$",               all_rows_int)
      inter_rows <- grep(paste0("^geno.*:", int_var,
                                "|^", int_var, ":.*geno"),
                         all_rows_int)
      show_rows  <- c(main_rows, inter_rows)

      # reference haplotype row (OR = 1) 
      tbl_int$addRow(
        rowKey = "base", values = list(
          term      = paste0(base_label, " (Ref)"),
          effect    = if (response_type == "binary") 1.0 else 0.0,
          ciLow     = '',
          ciHigh    = '',
          pval      = ''
        )
      )

      for (r in show_rows) {
        raw_nm  <- all_rows_int[r]
        
        # Get the display label
        label <- raw_to_label[raw_nm]
        if (is.na(label) || length(label) == 0) {
          # Fallback: try to create a sensible label
          suffix <- sub("^geno", "", raw_nm)
          suffix <- sub(paste0(":", int_var, "$"), "", suffix)
          if (grepl("-", suffix)) {
            label <- suffix
          } else if (grepl("rare", suffix, ignore.case = TRUE)) {
            label <- rare_label
          } else {
            label <- suffix
          }
          if (grepl(paste0(":", int_var, "$"), raw_nm)) {
            label <- paste0(label, " \u00D7 ", int_var)
          }
        }

        beta <- coef_sum_int[r, "coef"]
        pval <- coef_sum_int[r, "pval"]
        if (!is.null(ci_int) && raw_nm %in% rownames(ci_int)) {
          ci_lo <- ci_int[raw_nm, 1]; ci_hi <- ci_int[raw_nm, 2]
        } else {
          z     <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
          se    <- coef_sum_int[r, "se"]
          ci_lo <- beta - z * se; ci_hi <- beta + z * se
        }

        is_inter_term <- r %in% inter_rows
        tbl_int$addRow(
          rowKey = paste0("hi", r),
          values = list(
            term   = label,
            effect = if (response_type == "binary") exp(beta)  else beta,
            ciLow  = if (response_type == "binary") exp(ci_lo) else ci_lo,
            ciHigh = if (response_type == "binary") exp(ci_hi) else ci_hi,
            pval   = pval
          )
        )
      }
    }
  }

  # ── LRT interaction note (set regardless of fit outcome) ──────
  if (!is.na(p_inter_haplo))
    tbl_int$setNote(
      note = paste0("Likelihood ratio test for interaction: P = ",
                    format.pval(p_inter_haplo, digits = 3)),
      key  = "lrt_inter")
  else
    tbl_int$setNote(note = NULL, key = "lrt_inter")
}
    },

    # ── Private LD storage for plot render ────────────────────────
    .ld_store  = NULL,
    .ld_nms    = NULL,
    .ld_metric = NULL
  )
)  # end R6Class
