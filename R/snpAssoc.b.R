#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype
source("R/snp_helpers.R")

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
#' 
fit_interaction_model <- function(snp_enc, response, covariates_df,
                                  interaction_var, model_name,
                                  response_type, ci_width,
                                  conditional = FALSE, cond_var = interaction_var) {
  # VM fixes table but breaks lables
  df <- data.frame(resp = response, snp = as.factor(snp_enc))
  adj_covs <- character(0)
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df       <- cbind(df, covariates_df)
    adj_covs <- setdiff(names(covariates_df), interaction_var)
  }
  if (!(interaction_var %in% names(df))) return(NULL)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5) return(NULL)
  adj_part <- if (length(adj_covs) > 0) paste("+", paste(adj_covs, collapse = "+")) else ""

  if (conditional) {
    # ── Build nested formula based on conditioning direction ──────────────────
    if (cond_var == "snp") {
      # covariate effect within each SNP level
      formula_fit  <- as.formula(paste("resp ~ snp /", interaction_var, adj_part))
    } else {
      # SNP effect within each covariate level
      formula_fit <- as.formula(paste("resp ~", interaction_var, "/ snp", adj_part))
    }
    # Additive (no-interaction) formula for LRT
    formula_add <- as.formula(paste("resp ~ snp +", interaction_var, adj_part))

    if (response_type == "binary") {
      fit     <- glm(formula_fit, data = df, family = binomial())
      fit_add <- glm(formula_add, data = df, family = binomial())
      pval_col <- "Pr(>|z|)"; lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"
    } else {
      fit     <- lm(formula_fit, data = df)
      fit_add <- lm(formula_add, data = df)
      pval_col <- "Pr(>|t|)"; lrtest <- "F"; lrtest_label <- "Pr(>F)"
    }

    # ── LRT for interaction (nested vs additive) ──────────────────────────────
    lrt_cond <- tryCatch(anova(fit_add, fit, test = lrtest), error = function(e) NULL)
    p_inter  <- if (!is.null(lrt_cond)) lrt_cond[2, lrtest_label] else NA_real_

    coefs   <- summary(fit)$coefficients
    ci_mat  <- tryCatch(confint(fit, level = ci_width / 100),
                        error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2))
    aic_val <- AIC(fit)
    all_rows <- rownames(coefs)

    # ── Classify rows ─────────────────────────────────────────────────────────
    # Nested/interaction terms always contain ":"
    inter_rows_idx <- grep(":", all_rows)
    # For snp*x (numeric snp): standalone covariate main-effect terms also belong
    # to the conditional set (they represent the covariate effect in the ref group)
    if (cond_var == "snp" && !is.factor(snp_enc)) {
      cond_extra_idx <- grep(paste0("^", interaction_var), all_rows)
    } else {
      cond_extra_idx <- integer(0)
    }
    # SNP main-effect rows
    snp_rows_idx   <- grep("^snp", all_rows)
    snp_rows_idx   <- setdiff(snp_rows_idx, inter_rows_idx)
    # Adjustment covariate rows: not intercept, not snp, not interaction_var, not ":"
    adj_rows_idx   <- setdiff(
      seq_along(all_rows),
      c(grep("^\\(Intercept\\)", all_rows),
        snp_rows_idx, inter_rows_idx, cond_extra_idx,
        grep(paste0("^", interaction_var), all_rows)))

    if (length(inter_rows_idx) == 0 && length(cond_extra_idx) == 0) return(NULL)

    # All rows we might return (filtering by show_* happens in .fill_interaction)
    all_keep <- unique(c(snp_rows_idx, cond_extra_idx, inter_rows_idx, adj_rows_idx))

    first_inter_done <- FALSE
    lapply(all_keep, function(r) {
      idx  <- r
      term <- all_rows[r]
      beta <- coefs[idx, "Estimate"]
      pval <- coefs[idx, pval_col]
      ci_lo <- ci_mat[idx, 1]; ci_hi <- ci_mat[idx, 2]

      is_inter_term <- r %in% inter_rows_idx
      row_type <- if (r %in% snp_rows_idx)     "snp"
                  else if (is_inter_term)        "interaction"
                  else if (r %in% cond_extra_idx) "covariate"
                  else                            "adjustment"

      # For snp*x with numeric snp: interaction terms represent the *additional*
      # covariate effect vs ref group — combine with main covariate term (delta method).
      if (cond_var == "snp" && !is.factor(snp_enc) && grepl(":", term, fixed = TRUE)) {
        main_term <- sub(".*:", "", term)   # "snp:xMale" -> "xMale"
        main_idx  <- match(main_term, all_rows)
        if (!is.na(main_idx)) {
          beta     <- beta + coefs[main_idx, "Estimate"]
          vcov_mat <- vcov(fit)
          se       <- sqrt(vcov_mat[idx, idx] + vcov_mat[main_idx, main_idx] +
                           2 * vcov_mat[idx, main_idx])
          z        <- qnorm(1 - (1 - ci_width / 100) / 2)
          ci_lo    <- beta - z * se
          ci_hi    <- beta + z * se
          pval     <- 2 * pnorm(-abs(beta / se))
        }
      }

      # Attach p_inter to the first interaction (":") term only
      attach_p <- is_inter_term && !first_inter_done
      if (attach_p) first_inter_done <<- TRUE

      if (response_type == "binary")
        list(term = term, effect = exp(beta), ci_low = exp(ci_lo), ci_high = exp(ci_hi),
             pval = pval,
             pval_interaction = if (attach_p) p_inter else NA_real_,
             aic = aic_val, row_type = row_type)
      else
        list(term = term, effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval,
             pval_interaction = if (attach_p) p_inter else NA_real_,
             aic = aic_val, row_type = row_type)
    })
  } else {
    # ── Original * interaction logic ─────────────────────────────────────
    formula_int  <- as.formula(paste("resp ~ snp *", interaction_var, adj_part))
    formula_main <- as.formula(paste("resp ~ snp +", interaction_var, adj_part))
    tryCatch({
      if (response_type == "binary") {
        fit_int   <- glm(formula_int,  data = df, family = binomial())
        fit_main  <- glm(formula_main, data = df, family = binomial())
        pval_col  <- "Pr(>|z|)"; lrtest  <- "Chisq"; lrtest_label  <- "Pr(>Chi)"
      } else {
        fit_int   <- lm(formula_int,  data = df)
        fit_main  <- lm(formula_main, data = df)
        pval_col  <- "Pr(>|t|)"; lrtest  <- "F"; lrtest_label  <- "Pr(>F)"
      }
      lrt      <- tryCatch(anova(fit_main, fit_int, test = lrtest), error = function(e) NULL)
      p_inter  <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
      aic_val  <- AIC(fit_int)
      coefs    <- summary(fit_int)$coefficients
      ci       <- tryCatch(confint(fit_int, level = ci_width / 100),
                          error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2,
                                                      dimnames = list(rownames(coefs), c("lo", "hi"))))
      all_rows    <- rownames(coefs)
      snp_rows    <- grep("^snp", all_rows)
      inter_rows  <- grep(paste0("^snp.*:", interaction_var, "|^", interaction_var, ":.*snp"), all_rows)
      # covariate main-effect rows (interaction_var but not interaction terms)
      covar_rows  <- grep(paste0("^", interaction_var), all_rows)
      covar_rows  <- setdiff(covar_rows, inter_rows)
      # additional adjustment covariate rows: everything else except intercept, snp, covar, inter
      adj_rows    <- setdiff(seq_along(all_rows),
                             c(grep("^\\(Intercept\\)", all_rows),
                               snp_rows, inter_rows, covar_rows))
      # always return SNP + interaction rows; flag covar/adj rows for optional display
      keep_rows   <- unique(c(snp_rows, inter_rows, covar_rows, adj_rows))
      if (length(keep_rows) == 0) return(NULL)

      lapply(keep_rows, function(r) {
        beta     <- coefs[r, "Estimate"]
        pval     <- coefs[r, pval_col]
        ci_lo    <- ci[r, 1]; ci_hi  <- ci[r, 2]
        is_inter <- r %in% inter_rows
        row_type <- if (r %in% snp_rows)   "snp"
                    else if (r %in% inter_rows)  "interaction"
                    else if (r %in% covar_rows)  "covariate"
                    else                         "adjustment"
        if (response_type == "binary")
          list(term = all_rows[r], effect = exp(beta), ci_low = exp(ci_lo), ci_high = exp(ci_hi),
               pval = pval, pval_interaction = if (is_inter) p_inter else NA_real_,
               aic = aic_val, is_first = (r == keep_rows[1]), row_type = row_type)
        else
          list(term = all_rows[r], effect = beta, ci_low = ci_lo, ci_high = ci_hi,
               pval = pval, pval_interaction = if (is_inter) p_inter else NA_real_,
               aic = aic_val, is_first = (r == keep_rows[1]), row_type = row_type)
      })
    }, error = function(e) NULL)
  }
}

snpAssocClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "snpAssocClass",
  inherit = snpAssocBase,
  private = list(

    .init = function() {
      snp_names <- self$options$snps
      if (length(snp_names) == 0) return()
      arr <- self$results$snpResults
      for (nm in snp_names) arr$addItem(key = nm)
    },

    .run = function() {
      data           <- self$data
      opts           <- self$options
      response_var   <- opts$response
      snp_vars       <- opts$snps
      covariate_vars <- opts$covariates

      run_snpAssoc       <- isTRUE(opts$snpAssoc)
      run_snpInteraction <- isTRUE(opts$snpInteraction)

      # ── Validate SNPs ─────────────────────────────────────────────────────
      if (length(snp_vars) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'>Please add at least one SNP variable.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      } else {
        self$results$validationMsg$setVisible(FALSE)
      }


      val      <- validate_snp_vars(snp_vars, data)
      snp_vars <- val$valid_snps
      if (nchar(val$bad_html) > 0) {
        self$results$validationMsg$setContent(val$bad_html)
        self$results$validationMsg$setVisible(TRUE)
      } else {
        self$results$validationMsg$setVisible(FALSE)
      }
      if (length(snp_vars) == 0) return()

      # ── Response required ─────────────────────────────────────────────────
      if (is.null(response_var) || response_var == "") {
        self$results$validationMsg$setContent(
          "<p style='color:red;'>A response variable is required for association analysis.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      # ── Prepare response / covariates ─────────────────────────────────────
      response_raw  <- data[[response_var]]
      response_type <- detect_response_type(response_raw, opts$responseType)
      response      <- prepare_response(response_raw, response_type)
      cov_df        <- prepare_covariates(data, covariate_vars)

      if (run_snpInteraction && length(covariate_vars) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>SNP \u00D7 covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }

      # ── Complete-case mask ────────────────────────────────────────────────
      n_rows        <- nrow(data)
      complete_mask <- !is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0)
        complete_mask <- complete_mask & complete.cases(cov_df)

      # ── Per-SNP loop ──────────────────────────────────────────────────────
      arr <- self$results$snpResults
      for (snp_nm in snp_vars) {
        snp_raw     <- data[[snp_nm]]

        user_levels <- get_snp_level_order(snp_raw)
        geno_obj    <- parse_genotype(snp_raw, user_levels)
        if (is.null(geno_obj)) next

        # For each SNP, we need complete cases: no missing in SNP, response, or covariates
        snp_complete_mask <- complete_mask & !is.na(snp_raw)
        n_miss_assoc      <- n_rows - sum(snp_complete_mask)

        # _cc variables are the complete-case versions
        snp_raw_cc  <- snp_raw[snp_complete_mask]
        response_cc <- response[snp_complete_mask]
        response_raw_cc <- response_raw[snp_complete_mask]
        cov_df_cc   <- if (!is.null(cov_df)) cov_df[snp_complete_mask, , drop = FALSE] else NULL

        user_levels <- get_snp_level_order(snp_raw) # from data ordering 
        geno_obj_cc <- parse_genotype(snp_raw_cc, user_levels) # get genotype object
        if (is.null(geno_obj_cc)) next

        item <- arr$get(key = snp_nm)
        ref  <- get_ref_genotype(geno_obj_cc, user_levels)

        item$typingRate$setContent(sprintf(
          "<b>Typed samples:</b> %d / %d (%.1f%%)",
          sum(snp_complete_mask), n_rows, sum(snp_complete_mask) / n_rows * 100))

        if (run_snpAssoc)
          private$.fill_assoc(item$assocTable, snp_raw_cc, ref, response_cc,
                              cov_df_cc, response_type, opts,
                              n_miss = n_miss_assoc, user_levels, response_raw, snp_nm)

        if (run_snpInteraction && !is.null(cov_df_cc) && ncol(cov_df_cc) >= 1) {
          interaction_var <- names(cov_df_cc)[1]

          # Vector of models for interaction — currently single element from
          # the dropdown.  When interactionModel becomes multi-select, replace
          # this one line with: private$.get_interaction_models(opts)
          int_models <- opts$interactionModel

          private$.fill_interaction(
            item$interactionTable, snp_raw_cc, ref,
            response_cc, cov_df_cc, interaction_var,
            response_type, opts, int_models, user_levels, response_raw, snp_nm)

          if (isTRUE(opts$showStratByCovariate))
            private$.fill_strat_by_covariate(
              item$stratByCovariate, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels, response_raw, snp_nm)

          if (isTRUE(opts$showStratByGenotype))
            private$.fill_strat_by_genotype(
              item$stratByGenotype, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels, response_raw_cc, snp_nm)

          if (isTRUE(opts$showCrossClassTable))
            private$.fill_cross_class(
              item$crossClassTable, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels, response_raw_cc, snp_nm)
        }
      }
    },

    # ── Helper: build interaction model vector from options ───────────────────
    # Currently wraps the single List option.  To support multi-select in the
    # future: add Bool options interactionModelCodominant etc., then replace
    # the body of this function with the c(if...) pattern used in .fill_assoc.
    .get_interaction_models = function(opts) {
      opts$interactionModel   # returns e.g. "logadditive"
    },

    # ── Shared helpers ────────────────────────────────────────────────────────

    # Returns ordered genotype labels for each model row
    .geno_labels_for_model = function(mdl, all_genos, ref) {
      if (mdl %in% c("codominant", "logadditive")) return(all_genos)
      het  <- all_genos[all_genos != ref & all_genos != all_genos[length(all_genos)]]
      hom2 <- all_genos[length(all_genos)]
      if (length(het) == 0) het <- hom2
      if (mdl == "dominant")     return(c(ref, paste(c(het, hom2), collapse = "-")))
      if (mdl == "recessive")    return(c(paste(c(ref, het), collapse = "-"), hom2))
      if (mdl == "overdominant") return(c(paste(c(ref, hom2), collapse = "-"), het))
      all_genos
    },

    # Split combined label "A/B-C/D" into constituent genotypes
    .split_genos = function(gl)
      unlist(strsplit(gl, "(?<=[A-Za-z0-9*])-(?=[A-Za-z0-9*])", perl = TRUE)),

    # Compute N(%) per group (binary) or mean(SD) (quantitative)
    .compute_stats = function(geno_labels, snp_char, response, response_type) {
      split_genos <- private$.split_genos
    if (response_type == "binary") {
      lv       <- levels(as.factor(response))
      # Calculate column totals (total N for each response level)
      n_col0   <- sum(response == lv[1] & !is.na(response))
      n_col1   <- sum(response == lv[2] & !is.na(response))
      
      stats0   <- character(length(geno_labels))
      stats1   <- character(length(geno_labels))
      for (i in seq_along(geno_labels)) {
        mask  <- snp_char %in% split_genos(geno_labels[i]) & !is.na(response)
        n0    <- sum(mask & response == lv[1])
        n1    <- sum(mask & response == lv[2])
        
        if ((n0 + n1) == 0) { stats0[i] <- "---"; stats1[i] <- "---"; next }
        
        # Compute column-wise percentages
        pct0 <- if (n_col0 > 0) n0 / n_col0 * 100 else 0
        pct1 <- if (n_col1 > 0) n1 / n_col1 * 100 else 0
        
        stats0[i] <- sprintf("%d (%.1f%%)", n0, pct0)
        stats1[i] <- sprintf("%d (%.1f%%)", n1, pct1)
      }
      list(s0 = stats0, s1 = stats1)
    } else {
        stats0 <- character(length(geno_labels))
        for (i in seq_along(geno_labels)) {
          vals <- response[snp_char %in% split_genos(geno_labels[i]) & !is.na(response)]
          stats0[i] <- if (length(vals) == 0) "---"
                       else sprintf("%.2f (%.2f)", mean(vals), sd(vals))
        }
        list(s0 = stats0, s1 = rep("", length(geno_labels)))
      }
    },

    # BIC from AIC: BIC = AIC + df*(log(n) - 2)
    .bic_from_aic = function(aic_val, mdl, n_fit, n_cov) {
      snp_df <- c(codominant = 2L, dominant = 1L, recessive = 1L,
                  overdominant = 1L, logadditive = 1L)
      if (is.null(aic_val) || is.na(aic_val) || is.nan(aic_val)) return(NA_real_)
      round(aic_val + (1L + n_cov + snp_df[[mdl]]) * (log(n_fit) - 2), 2)
    },

    # ── Association table ─────────────────────────────────────────────────────
    .fill_assoc = function(tbl, snp_raw, ref, response, cov_df,
                           response_type, opts, n_miss = 0L, user_levels = NULL, response_raw, snp_lbl) {

      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
      tbl$getColumn("genotype")$setTitle(snp_lbl)

      resp_lbl <- attr(self$data[[self$options$response]], "label") %||% self$options$response
      tbl$setTitle(paste0("Association with ", resp_lbl))

      if (response_type == "binary") {
        lv       <- levels(as.factor(response_raw))
        tbl$getColumn("stat0")$setTitle(lv[1])
        tbl$getColumn("stat1")$setTitle(lv[2])
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(TRUE)
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)")
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }

      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))
      tbl$getColumn("BIC")$setVisible(isTRUE(opts$showAIC))

      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        cov_names <- sapply(names(cov_df), function(x) attr(self$data[[x]], "label") %||% x)
        note_txt  <- paste0("Model adjusted for: ", paste(cov_names, collapse = ", "))
        if (!is.na(n_miss) && n_miss > 0)
          note_txt <- paste0(note_txt, ".  ", n_miss, " observation(s) excluded.")
        tbl$setNote(note = note_txt, key = "covariates")
      } else if (!is.na(n_miss) && n_miss > 0) {
        tbl$setNote(note = paste0(n_miss, " observation(s) excluded."), key = "covariates")
      }

      models <- c(
        if (opts$modelCodominant)   "codominant",
        if (opts$modelDominant)     "dominant",
        if (opts$modelRecessive)    "recessive",
        if (opts$modelOverdominant) "overdominant",
        if (opts$modelLogAdditive)  "logadditive"
      )
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")

      snp_char  <- as.character(snp_raw)
      all_genos <- c(ref, setdiff(
        if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])),
        ref))
      n_fit <- sum(!is.na(snp_char) & !is.na(response) &
                     (if (!is.null(cov_df) && ncol(cov_df) > 0) complete.cases(cov_df) else TRUE))
      n_cov <- if (!is.null(cov_df)) ncol(cov_df) else 0L

      row_key <- 0L
      for (mdl in models) {
        snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
        res_list <- fit_model(snp_enc, response, cov_df, mdl, response_type, opts$ciWidth)
        if (is.null(res_list)) next

        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
        st          <- private$.compute_stats(geno_labels, snp_char, response, response_type)
        aic_val     <- { a <- res_list[[1]]$aic; if (!is.null(a) && !is.nan(a)) round(a, 2) else NA_real_ }
        bic_val     <- private$.bic_from_aic(aic_val, mdl, n_fit, n_cov)

        if (mdl == "logadditive") {
          res <- res_list[[1]]
          row_key <- row_key + 1L
          if (response_type == "binary") {
            lv <- levels(as.factor(response))
            stat0_val <- sprintf("%d", sum(response == lv[1], na.rm = TRUE))
            stat1_val <- sprintf("%d", sum(response == lv[2], na.rm = TRUE))
          } else {
            stat0_val <- sprintf("%.2f (%.2f)", mean(response, na.rm = TRUE), sd(response, na.rm = TRUE))
            stat1_val <- " "
          }
          
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model = model_labels[mdl], genotype = "Per allele", stat0 = stat0_val, stat1 = stat1_val,
            effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval = res$pval, AIC = aic_val, BIC = bic_val))

          next
        }

        # Reference row
        pval_row1 <- if (mdl == "codominant") res_list[[1]]$global_p else ''
        row_key <- row_key + 1L
        tbl$addRow(rowKey = as.character(row_key), values = list(
          model    = model_labels[mdl],
          genotype = geno_labels[1],
          stat0    = st$s0[1], stat1 = st$s1[1],
          effect   = if (response_type == "binary") 1. else 0.,
          ciLow = '', ciHigh = '', pval = pval_row1,
          AIC = aic_val, BIC = bic_val))
        
        if (mdl == "codominant") tbl$setNote(key = "lrt", note = "First p-value in Codominant is LRT for overall association")

        # Non-reference rows
        for (i in seq_along(res_list)) {
          res <- res_list[[i]]
          gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else res$comparison
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model    = "",
            genotype = gl,
            stat0    = if ((i + 1) <= length(st$s0)) st$s0[i + 1] else "-",
            stat1    = if ((i + 1) <= length(st$s1)) st$s1[i + 1] else "",
            effect   = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval     = res$pval, AIC = '', BIC = ''))
        }
      }
    },

    # ── Interaction omnibus table ─────────────────────────────────────────────
    .fill_interaction = function(tbl, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  int_models, user_levels = NULL, response_raw, snp_lbl) {
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")

      adj_vars   <- setdiff(names(cov_df), interaction_var)
      int_lbl    <- attr(self$data[[interaction_var]], "label") %||% interaction_var

      # ── Interaction type: formula parameterisation ──────────────────────────
      # "multiplicative"       snp * covar   (default)
      # "conditional_on_snp"   covar / snp   (covar effect within each SNP stratum)
      # "conditional_on_covar" snp / covar   (SNP effect within each covariate stratum)
      int_type <- if (is.null(opts$interactionType)) "multiplicative" else opts$interactionType

      # Human-readable formula token for table title
      formula_token <- switch(int_type,
        multiplicative       = paste0(snp_lbl, " \u00D7 ", int_lbl),
        conditional_on_snp   = paste0(int_lbl,  " | ", snp_lbl),
        conditional_on_covar = paste0(snp_lbl, " | ", int_lbl))
      tbl$setTitle(paste0(formula_token, " interaction"))

      if (length(adj_vars) > 0) {
        note_parts <- paste0("Adjusted for: ", paste(sapply(adj_vars, function(x)
          attr(self$data[[x]], "label") %||% x), collapse = ", "))
        tbl$setNote(note = note_parts, key = "intcov")
      }


      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))
      tbl$getColumn("BIC")$setVisible(isTRUE(opts$showAIC))

      # display filters
      show_adj     <- isTRUE(opts$showInteractionAdjVars)

      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")

      # helper: pretty-print a model coefficient term for display only.
      # res$term (the raw R coefficient name) is NEVER modified here — this
      # function is purely cosmetic and is called only in .fill_interaction.
      #
      # Term shapes by formula and model type:
      #
      #  snp * covar  (multiplicative, any encoding)
      #    codominant:  "snpA/G"          "snpA/G:SEXmale"   "SEXmale" (covar main)
      #    collapsed:   "snp"             "snp:SEXmale"       "SEXmale"
      #
      #  snp / covar  (conditional_on_snp, factor snp)
      #    codominant:  "snpA/G:SEXmale"  (nested covar within genotype)
      #    collapsed:   "snp:SEXmale"
      #
      #  covar / snp  (conditional_on_covar)
      #    codominant:  "SEXmale:snpA/G"  "SEXmale:snpG/G"   (snp after colon)
      #    collapsed:   "SEXmale:snp"                         (bare snp after colon)
      #
      # In every case the SNP token to replace is either:
      #   (a) leading "snp" possibly followed by a genotype suffix then ":" or EOL
      #   (b) ":snp" anywhere in the string, possibly followed by a genotype suffix
      #
      # We detect whether a genotype suffix is already present (codominant) vs absent
      # (collapsed/logadditive), and in the latter case inject geno_labels[2].
      .label_term <- function(term, mdl, geno_labels) {
        geno_suffix <- if (mdl == "logadditive") "per allele"
                       else if (length(geno_labels) >= 2) geno_labels[2]
                       else ""
        tag_collapsed <- if (nchar(geno_suffix) > 0) paste0("(", geno_suffix, ")") else ""

        # ── Case A: leading "snp" (multiplicative or conditional_on_snp) ─────
        # Codominant: "snpA/G:SEX" or standalone "snpA/G"  → wrap suffix in parens
        # Collapsed:  "snp:SEX"    or standalone "snp"      → inject geno tag
        lbl <- gsub("^snp([^:]+)",               # codominant leading: has suffix
                    paste0(snp_lbl, "(\\1)"), term)
        if (lbl == term) {
          # No match above → collapsed leading: bare "snp" (possibly "snp:...")
          lbl <- gsub("^snp(?=:|$)",             # bare snp at start, before ":" or EOL
                      paste0(snp_lbl, tag_collapsed), term, perl = TRUE)
        }

        # ── Case B: ":snp" anywhere (conditional_on_covar, covar / snp) ──────
        # Codominant: "SEXmale:snpA/G"  → "SEXmale:SNP1(A/G)"
        # Collapsed:  "SEXmale:snp"     → "SEXmale:SNP1(A/G-G/G)"
        # (Also handles multiplicative covar-leading cross terms if they arise)
        lbl <- gsub(":snp([^:]+)",               # codominant trailing: has suffix
                    paste0(":", snp_lbl, "(\\1)"), lbl)
        lbl <- gsub(":snp(?=:|$)",               # collapsed trailing: bare snp
                    paste0(":", snp_lbl, tag_collapsed), lbl, perl = TRUE)

        lbl
      }

      row_key <- 0L
      for (mdl in int_models) { # for now only  model can be selected, but keep loop structure for future multi-select
        snp_enc  <- encode_model(as.character(snp_raw), ref, mdl, user_levels)

        # ── Choose conditional flag and cond_var based on interactionType ──
        if (int_type == "multiplicative") {
          res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                            interaction_var, mdl, response_type, opts$ciWidth,
                                            conditional = FALSE)
        } else if (int_type == "conditional_on_snp") {
          res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                            interaction_var, mdl, response_type, opts$ciWidth,
                                            conditional = TRUE, cond_var = "snp")
        } else {   # conditional_on_covar
          res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                            interaction_var, mdl, response_type, opts$ciWidth,
                                            conditional = TRUE, cond_var = interaction_var)
        }
        if (is.null(res_list)) next

        # geno_labels needed for term labelling in collapsed models
        snp_char_l  <- as.character(snp_raw)
        all_genos_l <- c(ref, setdiff(
          if (!is.null(user_levels)) user_levels
          else sort(unique(snp_char_l[!is.na(snp_char_l)])), ref))
        geno_labels_l <- private$.geno_labels_for_model(mdl, all_genos_l, ref)

        # BIC denominator
        n_fit_bic <- sum(!is.na(snp_enc) & !is.na(response) & complete.cases(cov_df))
        n_cov_bic <- ncol(cov_df)

        first_row <- TRUE; first_inter <- TRUE

        for (res in res_list) {
          # Determine whether to show this row based on its type
          rtype <- if (is.null(res$row_type)) "snp" else res$row_type
          # - "adjustment" rows: optional in all model types
          if (rtype == "adjustment" && !show_adj)   next

          row_key  <- row_key + 1L
          is_inter <- !is.na(res$pval_interaction)

          term_label <- .label_term(res$term, mdl, geno_labels_l)

          # p_interaction: show on the first term that carries it (any model type)
          pval_int_val <- if (!is.na(res$pval_interaction) && first_inter)
                            res$pval_interaction else ""

          vals <- list(
            model  = if (first_row) model_labels[mdl] else "",
            term   = term_label,
            effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval   = res$pval,
            pvalInteraction = pval_int_val)

          if (isTRUE(opts$showAIC)) {
            aic_val <- if (first_row && !is.nan(res$aic)) round(res$aic, 2) else ""
            bic_val <- if (first_row && !is.nan(res$aic))
              private$.bic_from_aic(res$aic, mdl, n_fit_bic, n_cov_bic) else ""
            vals[["AIC"]] <- aic_val
            vals[["BIC"]] <- bic_val
          }

          tbl$addRow(rowKey = as.character(row_key), values = vals)
          first_row <- FALSE
          if (!is.na(res$pval_interaction)) first_inter <- FALSE
        }
      }
    },

    .fill_strat_by_covariate = function(arr, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  int_models, user_levels = NULL, response_raw, snp_lbl) {
      int_var_data <- cov_df[[interaction_var]]
      if (length(table(int_var_data)) > 6) { # numerical covariate: skip stratified results
        return()
      }
      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      snp_char   <- as.character(snp_raw)
      all_genos  <- c(ref, setdiff(if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])), ref))
      model_labels <- c(codominant="Codominant", dominant="Dominant", recessive="Recessive", overdominant="Overdominant", logadditive="Log-additive")
      adj_vars     <- setdiff(names(cov_df), interaction_var)
      adj_cov_df   <- if (length(adj_vars) > 0) cov_df[, adj_vars, drop=FALSE] else NULL
      cov_levels   <- if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(int_var_data[!is.na(int_var_data)]))

      # VM
      # Fit to compute omnibus interaction p-value. Could be extracted from the conditional model fit below    
      mdl_for_p <- if ("codominant" %in% int_models) "codominant" else int_models[1]
      int_res_p <- fit_interaction_model(encode_model(snp_char, ref, mdl_for_p, user_levels), response, cov_df,
                                        interaction_var, mdl_for_p, response_type, opts$ciWidth, conditional = FALSE)
      p_inter <- if (!is.null(int_res_p)) int_res_p[[1]]$pval_interaction else NA_real_


      # VM
      # table is OK. bu why the interaction models is refitted for each cov_level?

      for (cl in cov_levels) {
        cl_label <- as.character(cl)
        key_k    <- paste0(int_lbl, ": ", cl_label)
        if (is.null(tryCatch(arr$get(key = key_k), error = function(e) NULL))) arr$addItem(key = key_k)
        tbl <- arr$get(key = key_k)
        tbl$getColumn("genotype")$setTitle(snp_lbl)
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
        if (response_type == "binary") {
          resp_lv <- levels(as.factor(response_raw))
          tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
          tbl$getColumn("stat0")$setVisible(TRUE); tbl$getColumn("stat1")$setVisible(TRUE)
        } else {
          tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
          tbl$getColumn("stat1")$setVisible(FALSE)
        }
        if (!is.na(p_inter)) tbl$setNote(key = "pinter", note = paste0("Interaction p-value: ", format.pval(p_inter, digits=3, eps=0.001)))

        # Descriptive stats for this stratum
        mask_k <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
        if (!is.null(adj_cov_df) && ncol(adj_cov_df) > 0) mask_k <- mask_k & complete.cases(adj_cov_df)
        snp_c_k    <- snp_char[mask_k]
        response_k <- response[mask_k]

        row_key <- 0L
        for (mdl in int_models) {
          # Fit conditional model on FULL data to keep adjustment, then filter
          res_list <- fit_interaction_model(
            snp_enc = encode_model(snp_char, ref, mdl, user_levels), response = response,
            covariates_df = cov_df, interaction_var = interaction_var, model_name = mdl,
            response_type = response_type, ci_width = opts$ciWidth, conditional = TRUE)
          if (is.null(res_list)) next

          # Extract only nested/interaction terms (row_type "interaction"), then
          # match by covariate level label.  "snp" main-effect and "covariate" rows
          # must be excluded — the same contamination fix applied to fill_strat_by_genotype.
          # For collapsed models with snp*covar the reference covariate level produces
          # no ":" term at all (its SNP effect is in the bare "snp" main-effect row);
          # level_res will legitimately be empty for that stratum, but we still want
          # to render the table (ref genotype OR=1, non-ref genotype rows get no estimate).
          inter_only <- res_list[sapply(res_list, function(r)
            is.null(r$row_type) || r$row_type == "interaction")]
          level_res  <- inter_only[grepl(cl_label, sapply(inter_only, `[[`, "term"), fixed = TRUE)]
          is_ref_cov <- cl == cov_levels[1]
          if (length(level_res) == 0 && !is_ref_cov) next

          geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
          st          <- private$.compute_stats(geno_labels, snp_c_k, response_k, response_type)

          if (mdl == "logadditive") {
            res <- level_res[[1]]
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              genotype = paste0(model_labels[mdl], " (per allele)"), stat0 = "---", stat1 = " ",
              effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high, pval = res$pval))
            next
          }

          # Reference row
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            genotype = geno_labels[1], stat0 = st$s0[1], stat1 = st$s1[1],
            effect = if (response_type == "binary") 1.0 else 0.0, ciLow = "", ciHigh = "", pval = ""))

          # Comparison rows
          # For the reference covariate level of collapsed models, level_res is empty
          # because the SNP effect at the ref covariate level is the bare "snp" main-effect
          # term (row_type "snp"), not a ":" term.  Fall back to those rows.
          cmp_res <- if (length(level_res) > 0) level_res else
            res_list[sapply(res_list, function(r) !is.null(r$row_type) && r$row_type == "snp")]
          for (i in seq_along(cmp_res)) {
            res <- cmp_res[[i]]
            gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else sub("snp", "", res$term)
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              genotype = gl,
              stat0    = if ((i + 1) <= length(st$s0)) st$s0[i + 1] else "-",
              stat1    = if ((i + 1) <= length(st$s1)) st$s1[i + 1] else " ",
              effect   = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high, pval = res$pval))
          }
        }
      }
    }, 
    .fill_strat_by_genotype = function(arr, snp_raw, ref, response, cov_df,
                                   interaction_var, response_type, opts,
                                   int_models, user_levels = NULL, response_raw, snp_lbl) {
      snp_char     <- as.character(snp_raw)
      all_genos    <- c(ref, setdiff(if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])), ref))
      int_var_data <- cov_df[[interaction_var]]
      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      resp_lv      <- levels(as.factor(response_raw))

      is_numerical <- length(unique(int_var_data)) > 6 && sum(is.na(as.numeric(int_var_data))) == 0
      if (!is_numerical) {
        cov_levels <- if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(as.character(int_var_data[!is.na(int_var_data)])))
      } else {
        cov_levels <- interaction_var
      }

      for (mdl in int_models) {
        if (mdl == "logadditive") next
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)

        snp_enc_m <- encode_model(snp_char, ref, mdl, user_levels)
        res_list  <- fit_interaction_model(snp_enc_m, response, cov_df,
                                           interaction_var, mdl, response_type, opts$ciWidth,
                                           conditional = TRUE, cond_var = "snp")
        if (is.null(res_list)) next

        # Number of res_list entries per genotype group:
        # categorical: one per non-reference covariate level; numerical: one
        n_cov_contrasts <- if (is_numerical) 1L else max(1L, length(cov_levels) - 1L)

        for (gl in geno_labels) {
          gl_idx <- match(gl, geno_labels)
          key_g  <- paste0(snp_lbl, ": ", gl)
          if (is.null(tryCatch(arr$get(key = key_g), error = function(e) NULL))) arr$addItem(key = key_g)
          tbl <- arr$get(key = key_g)

          tbl$getColumn("level")$setTitle(int_lbl)
          tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")

          if (response_type == "binary") {
            tbl$getColumn("stat0")$setTitle(resp_lv[1])
            tbl$getColumn("stat1")$setTitle(resp_lv[2])
            tbl$getColumn("stat0")$setVisible(TRUE)
            tbl$getColumn("stat1")$setVisible(TRUE)
          } else {
            tbl$getColumn("stat0")$setTitle("Mean (SD)")
            tbl$getColumn("stat0")$setVisible(TRUE)
            tbl$getColumn("stat1")$setVisible(FALSE)
          }

          # split_genos handles combined labels like "A/G-G/G" for aggregated models
          mask_g     <- snp_char %in% private$.split_genos(gl)
          int_g      <- int_var_data[mask_g]
          resp_g     <- response[mask_g]
          resp_raw_g <- response_raw[mask_g]

          if (response_type == "binary") {
            counts <- table(factor(int_g, levels = cov_levels),
                            factor(resp_raw_g, levels = resp_lv))
            totals <- colSums(counts)
          }

          # Filter res_list to this genotype group.
          # We want only the nested covariate-within-genotype terms (row_type "interaction"),
          # which contain ":". The "snp" main-effect rows (e.g. "snpA/G") must be excluded —
          # they are the SNP effect at the reference covariate level, not what this table shows.
          inter_only <- res_list[sapply(res_list, function(r)
            is.null(r$row_type) || r$row_type == "interaction")]
          gl_res <- inter_only[grepl(gl, sapply(inter_only, `[[`, "term"), fixed = TRUE)]
          if (length(gl_res) == 0) {
            # Aggregated models: no gl label in term names → positional slice over inter_only
            has_ref_terms <- length(inter_only) >= length(geno_labels) * n_cov_contrasts
            gl_offset <- if (has_ref_terms) gl_idx - 1L else gl_idx - 2L
            start  <- gl_offset * n_cov_contrasts + 1L
            end    <- min((gl_offset + 1L) * n_cov_contrasts, length(inter_only))
            gl_res <- if (start >= 1L && start <= length(inter_only)) inter_only[start:end] else list()
          }

          row_key <- 0L
          if (!is_numerical) {
            # Reference covariate level always OR=1 — no model term emitted for it
            cl_ref  <- cov_levels[1]
            stat0   <- if (response_type == "binary") fmt_cat(counts[cl_ref, 1], totals[1]) else { vals <- resp_g[as.character(int_g) == cl_ref]; fmt_cont(vals) }
            stat1   <- if (response_type == "binary") fmt_cat(counts[cl_ref, 2], totals[2]) else ""
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              level  = cl_ref, stat0 = stat0, stat1 = stat1,
              effect = if (response_type == "binary") 1.0 else 0.0,
              ciLow = "", ciHigh = "", pval = ""))

            # Non-reference covariate levels: positional match into gl_res
            for (i in seq_along(cov_levels[-1])) {
              cl    <- cov_levels[-1][i]
              res   <- if (i <= length(gl_res)) gl_res[[i]] else NULL
              stat0 <- if (response_type == "binary") fmt_cat(counts[cl, 1], totals[1]) else { vals <- resp_g[as.character(int_g) == cl]; fmt_cont(vals) }
              stat1 <- if (response_type == "binary") fmt_cat(counts[cl, 2], totals[2]) else ""
              row_key <- row_key + 1L
              tbl$addRow(rowKey = as.character(row_key), values = list(
                level  = cl, stat0 = stat0, stat1 = stat1,
                effect = if (!is.null(res)) res$effect else if (response_type == "binary") 1.0 else 0.0,
                ciLow  = if (!is.null(res)) res$ci_low  else "",
                ciHigh = if (!is.null(res)) res$ci_high else "",
                pval   = if (!is.null(res)) res$pval    else ""))
            }
          } else {
            # Numerical covariate: single summary row, one term per genotype group
            stat0 <- if (response_type == "binary") fmt_cont(int_g[resp_raw_g == resp_lv[1]]) else fmt_cont(resp_g)
            stat1 <- if (response_type == "binary") fmt_cont(int_g[resp_raw_g == resp_lv[2]]) else ""
            res   <- if (length(gl_res) > 0) gl_res[[1]] else NULL
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              level  = "Overall", stat0 = stat0, stat1 = stat1,
              effect = if (!is.null(res)) res$effect else if (response_type == "binary") 1.0 else 0.0,
              ciLow  = if (!is.null(res)) res$ci_low  else "",
              ciHigh = if (!is.null(res)) res$ci_high else "",
              pval   = if (!is.null(res)) res$pval    else ""))
          }
        }
      }
    }
  )
)
