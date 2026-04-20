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
  df <- data.frame(resp = response, snp = snp_enc)
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
    # encode_model always returns a factor (except logadditive), so snp / covar
    # always produces proper nested terms for every model.
    if (cond_var == "snp") {
      formula_fit <- as.formula(paste("resp ~ snp /", interaction_var, adj_part))
    } else {
      formula_fit <- as.formula(paste("resp ~", interaction_var, "/ snp", adj_part))
    }
    # Additive formula for LRT
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
    p_inter  <- lrt_cond[2, lrtest_label]

    coefs    <- summary(fit)$coefficients
    ci_mat   <- tryCatch(confint(fit, level = ci_width / 100),
                         error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2))
    aic_val  <- AIC(fit)
    all_rows <- rownames(coefs)

    # ── Classify rows ─────────────────────────────────────────────────────────
    # Nested terms always contain ":"
    inter_rows_idx <- grep(":", all_rows)
    # SNP main-effect rows (non-nested; represent SNP effect at ref covariate level)
    snp_rows_idx   <- setdiff(grep("^snp", all_rows), inter_rows_idx)
    # Adjustment covariate rows
    adj_rows_idx   <- setdiff(
      seq_along(all_rows),
      c(grep("^\\(Intercept\\)", all_rows),
        snp_rows_idx, inter_rows_idx,
        grep(paste0("^", interaction_var), all_rows)))

    if (length(inter_rows_idx) == 0) return(NULL)

    all_keep <- unique(c(snp_rows_idx, inter_rows_idx, adj_rows_idx))

    first_inter_done <- FALSE
    result <- lapply(all_keep, function(r) {
      term  <- all_rows[r]
      beta  <- coefs[r, "Estimate"]
      pval  <- coefs[r, pval_col]
      ci_lo <- ci_mat[r, 1]; ci_hi <- ci_mat[r, 2]

      is_inter_term <- r %in% inter_rows_idx
      row_type <- if (r %in% snp_rows_idx) "snp"
                  else if (is_inter_term)   "interaction"
                  else                      "adjustment"

      attach_p <- is_inter_term && !first_inter_done
      if (attach_p) first_inter_done <<- TRUE

      if (response_type == "binary")
        list(term = term, effect = exp(beta), ci_low = exp(ci_lo), ci_high = exp(ci_hi),
             pval = pval, pval_interaction = if (attach_p) p_inter else NA_real_,
             aic = aic_val, row_type = row_type)
      else
        list(term = term, effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, pval_interaction = if (attach_p) p_inter else NA_real_,
             aic = aic_val, row_type = row_type)
    })
    attr(result, "pval_interaction") <- p_inter
    result
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

      result <- lapply(keep_rows, function(r) {
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
      attr(result, "pval_interaction") <- p_inter
      result
    }, error = function(e) NULL)
    # VM
  } # capture return object and add pval_interaction as attribute to be used in notes for stratified tables
}

snpAssocClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "snpAssocClass",
  inherit = snpAssocBase,
  private = list(

    .init = function() {
      snp_names    <- self$options$snps
      if (length(snp_names) == 0) return()
      arr          <- self$results$snpResults
      int_models   <- private$.get_interaction_models(self$options)
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      for (nm in snp_names) {
        arr$addItem(key = nm)
        if (isTRUE(self$options$snpInteraction) && length(int_models) > 0) {
          int_arr <- arr$get(key = nm)$interactionResults
          for (mdl in int_models)
            int_arr$addItem(key = model_labels[[mdl]])
        }
      }
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
        # VM as.factor
        snp_raw_cc  <- as.factor(snp_raw[snp_complete_mask])
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
          int_models      <- private$.get_interaction_models(opts)
          model_labels    <- c(codominant = "Codominant", dominant = "Dominant",
                               recessive  = "Recessive",  overdominant = "Overdominant",
                               logadditive = "Log-additive")

          int_lbl <- attr(self$data[[interaction_var]], "label") %||% interaction_var

          for (mdl in int_models) {
            mdl_label  <- model_labels[[mdl]]
            mdl_item   <- item$interactionResults$get(key = mdl_label)

            if (isTRUE(opts$showInteractionTable))
              private$.fill_interaction(
                mdl_item$interactionTable, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw, snp_nm)

            if (isTRUE(opts$showStratByCovariate)) {
              mdl_item$stratByCovariateHeading$setContent(
                paste0("<h3>Stratified by Covariate: ", int_lbl, "</h3>"))
              private$.fill_strat_by_covariate(
                mdl_item$stratByCovariate, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw, snp_nm)
            }

            if (isTRUE(opts$showStratByGenotype)) {
              mdl_item$stratByGenotypeHeading$setContent(
                paste0("<h3>Stratified by Genotype: ", snp_nm, "</h3>"))
              private$.fill_strat_by_genotype(
                mdl_item$stratByGenotype, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw_cc, snp_nm)
            }

            if (isTRUE(opts$showCrossClassTable)) {
              mdl_item$crossClassHeading$setContent(
                paste0("<h3>Cross-Classification: ", snp_nm, " × ", int_lbl, "</h3>"))
              private$.fill_cross_class(
                mdl_item$crossClassTable, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw_cc, snp_nm)
            }
          }
        }
      }
    },

    # ── Helper: build interaction model vector from options ───────────────────
    # Reuses the same genetic model checkboxes as the association table.
    .get_interaction_models = function(opts) {
      c(
        if (isTRUE(opts$modelCodominant))   "codominant",
        if (isTRUE(opts$modelDominant))     "dominant",
        if (isTRUE(opts$modelRecessive))    "recessive",
        if (isTRUE(opts$modelOverdominant)) "overdominant",
        if (isTRUE(opts$modelLogAdditive))  "logadditive"
      )
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
      tbl$setTitle(paste0("<b>", formula_token, " interaction</b>"))

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
      # res$term (the raw R coefficient name) is NEVER modified here.
      #
      # Since encode_model now returns a named factor for all diploid models,
      # every term already embeds the full genotype label, e.g.:
      #   codominant:   "snpA/G:SEXmale"   "SEXmale:snpA/G"
      #   dominant:     "snpA/G-G/G:SEXmale"
      #   recessive:    "snpG/G:SEXmale"
      # We just replace the "snp" token with snp_lbl and wrap the suffix in parens.
      # logadditive is the only model that keeps a numeric snp (no suffix) —
      # it never appears in the interaction table so no special case is needed.
      .label_term <- function(term, mdl, geno_labels) {
        # Case A: leading "snp<suffix>" — multiplicative or conditional_on_snp
        lbl <- gsub("^snp([^:]+)", paste0(snp_lbl, "(\\1)"), term)
        # Case B: ":snp<suffix>" — conditional_on_covar  (covar / snp)
        lbl <- gsub(":snp([^:]+)", paste0(":", snp_lbl, "(\\1)"), lbl)
        lbl
      }

      row_key <- 0L
      last_pval_interaction <- NA_real_
      for (mdl in int_models) { # int_models may be a single string (called per-model from .run) or a vector
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

        p_inter <- attr(res_list, "pval_interaction")
        if (!is.null(p_inter) && !is.na(p_inter)) last_pval_interaction <- p_inter

        # geno_labels needed for term labelling in collapsed models
        snp_char_l  <- as.character(snp_raw)
        all_genos_l <- c(ref, setdiff(
          if (!is.null(user_levels)) user_levels
          else sort(unique(snp_char_l[!is.na(snp_char_l)])), ref))
        geno_labels_l <- private$.geno_labels_for_model(mdl, all_genos_l, ref)

        # BIC denominator
        n_fit_bic <- sum(!is.na(snp_enc) & !is.na(response) & complete.cases(cov_df))
        n_cov_bic <- ncol(cov_df)

        first_row <- TRUE

        for (res in res_list) {
          # Determine whether to show this row based on its type
          rtype <- if (is.null(res$row_type)) "snp" else res$row_type
          # - "adjustment" rows: optional in all model types
          if (rtype == "adjustment" && !show_adj)   next

          row_key  <- row_key + 1L

          term_label <- .label_term(res$term, mdl, geno_labels_l)

          vals <- list(
            model  = if (first_row) model_labels[mdl] else "",
            term   = term_label,
            effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval   = res$pval)

          if (isTRUE(opts$showAIC)) {
            aic_val <- if (first_row && !is.nan(res$aic)) round(res$aic, 2) else ""
            bic_val <- if (first_row && !is.nan(res$aic))
              private$.bic_from_aic(res$aic, mdl, n_fit_bic, n_cov_bic) else ""
            vals[["AIC"]] <- aic_val
            vals[["BIC"]] <- bic_val
          }

          tbl$addRow(rowKey = as.character(row_key), values = vals)
          first_row <- FALSE
        }
      }

      # ── Interaction p-value note ──────────────────────────────────────────────
      if (!is.na(last_pval_interaction)) {
        tbl$setNote(
          note = paste0("Interaction p-value (LRT): ", format.pval(last_pval_interaction, digits = 3, eps = 0.001)),
          key  = "interactionPval")
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

      for (cl in cov_levels) {
        cl_label <- as.character(cl)
        key_k    <- paste0(int_lbl, ": ", cl_label)
        if (is.null(tryCatch(arr$get(key = key_k), error = function(e) NULL))) arr$addItem(key = key_k)
        tbl <- arr$get(key = key_k)
        tbl$setTitle(paste0("<b>", int_lbl, ": ", cl_label, "</b>"))
        tbl$getColumn("genotype")$setTitle(paste0("<b>", snp_lbl, "</b>"))
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
        if (response_type == "binary") {
          resp_lv <- levels(as.factor(response_raw))
          tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
          tbl$getColumn("stat0")$setVisible(TRUE); tbl$getColumn("stat1")$setVisible(TRUE)
        } else {
          tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
          tbl$getColumn("stat1")$setVisible(FALSE)
        }
      }

      for (mdl in int_models) {
        res_list <- fit_interaction_model(
          snp_enc = encode_model(snp_char, ref, mdl, user_levels), response = response,
          covariates_df = cov_df, interaction_var = interaction_var, model_name = mdl,
          response_type = response_type, ci_width = opts$ciWidth, conditional = TRUE)
        if (is.null(res_list)) next

        inter_only  <- res_list[sapply(res_list, function(r)
          is.null(r$row_type) || r$row_type == "interaction")]
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)

        for (cl in cov_levels) {
          cl_label  <- as.character(cl)
          level_res <- inter_only[grepl(cl_label, sapply(inter_only, `[[`, "term"), fixed = TRUE)]
          if (length(level_res) == 0) next

          tbl <- arr$get(key = paste0(int_lbl, ": ", cl_label))

          mask_k     <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
          if (!is.null(adj_cov_df) && ncol(adj_cov_df) > 0) mask_k <- mask_k & complete.cases(adj_cov_df)
          st <- private$.compute_stats(geno_labels, snp_char[mask_k], response[mask_k], response_type)

          row_key <- 0L

          if (mdl == "logadditive") {
            tbl$addRow(rowKey = "1", values = list(
              genotype = paste0(model_labels[mdl], " (per allele)"), stat0 = "", stat1 = "",   # perhaps here I could add just N for the level?
              effect = level_res[[1]]$effect, ciLow = level_res[[1]]$ci_low,
              ciHigh = level_res[[1]]$ci_high, pval = level_res[[1]]$pval))
            tbl$getColumn("stat1")$setVisible(FALSE); tbl$getColumn("stat0")$setVisible(FALSE) # hide columns 
            next
          }

          # Reference genotype row (OR = 1 by definition)
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            genotype = geno_labels[1], stat0 = st$s0[1], stat1 = st$s1[1],
            effect = if (response_type == "binary") 1.0 else 0.0, ciLow = "", ciHigh = "", pval = ""))

          for (i in seq_along(level_res)) {
            res <- level_res[[i]]
            gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else sub("snp", "", res$term)
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              genotype = gl,
              stat0    = if ((i + 1) <= length(st$s0)) st$s0[i + 1] else "-",
              stat1    = if ((i + 1) <= length(st$s1)) st$s1[i + 1] else " ",
              effect   = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high, pval = res$pval))
          }
        }
        note_txt <- paste0("The reference group is <b>",  snp_lbl, ": ",geno_labels[1], "</b> across all strata.")
        tbl$setNote(note = note_txt, key = "interStratCov")
        pval_interaction <- attr(res_list, "pval_interaction")
        note_pval <- paste0("Interaction p-value: ", format.pval(pval_interaction, digits = 3, eps = 0.001))
        tbl$setNote(note = note_pval, key = "interStratCovPval")

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
          tbl$setTitle(paste0("<b>", snp_lbl, ": ", gl, "</b>"))
          tbl$getColumn("level")$setTitle(paste0("<b>", int_lbl, "</b>"))
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
        note_txt <- paste0("The reference group is <b>", interaction_var, ": ", cov_levels[1], "</b> across all strata.")
        tbl$setNote(note = note_txt, key = "interStatGeno")
        pval_interaction <- attr(res_list, "pval_interaction")
        note_pval <- paste0("Interaction p-value: ", format.pval(pval_interaction, digits = 3, eps = 0.001))
        tbl$setNote(note = note_pval, key = "interStratGenoPval")
      }
    },
   .fill_cross_class = function(arr, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  int_models, user_levels = NULL, response_raw, snp_lbl) {
      if (int_models == "logadditive") return() # not meaningful 
   
      int_var_data <- cov_df[[interaction_var]]
      if (length(table(int_var_data)) > 6) return() # Skip numerical

      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      snp_char     <- as.character(snp_raw)
      all_genos    <- c(ref, setdiff(if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])), ref))
      cov_levels   <- if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(int_var_data[!is.na(int_var_data)]))
      adj_vars     <- setdiff(names(cov_df), interaction_var)
      
      # 1. Setup Array Items (Tables)
      for (cl in cov_levels) {
        key_k <- paste0(int_lbl, ": ", as.character(cl))
        if (is.null(tryCatch(arr$get(key = key_k), error = function(e) NULL))) arr$addItem(key = key_k)
        tbl <- arr$get(key = key_k)
        tbl$setTitle(paste0("<b>", int_lbl, ": ", as.character(cl), "</b>"))
        tbl$getColumn("genotype")$setTitle(paste0("<b>", snp_lbl, "</b>"))
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
        if (response_type == "binary") {
          resp_lv <- levels(as.factor(response_raw))
          tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
        } else {
          tbl$getColumn("stat0")$setTitle("Mean (SD)")
        }
      }

      for (mdl in int_models) {
        # 2. Fit Multiplicative Model: resp ~ snp * interaction_var + adj_vars
        snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
        df_fit   <- data.frame(resp = response, snp = snp_enc, interaction_var = int_var_data)
        if (length(adj_vars) > 0) df_fit <- cbind(df_fit, cov_df[, adj_vars, drop=FALSE])
        
        adj_part <- if (length(adj_vars) > 0) paste("+", paste(adj_vars, collapse = "+")) else ""
        formula_str     <- paste("resp ~ snp * interaction_var", adj_part)
        formula_main_str <- paste("resp ~ snp + interaction_var", adj_part)
        
        fit <- if (response_type == "binary") {
          glm(as.formula(formula_str), data = df_fit, family = binomial())
        } else {
          lm(as.formula(formula_str), data = df_fit)
        }
        fit_main_cc <- if (response_type == "binary") {
          glm(as.formula(formula_main_str), data = df_fit, family = binomial())
        } else {
          lm(as.formula(formula_main_str), data = df_fit)
        }

        # Compute interaction p-value via LRT
        lrtest_str   <- if (response_type == "binary") "Chisq" else "F"
        lrtest_label <- if (response_type == "binary") "Pr(>Chi)" else "Pr(>F)"
        lrt_cc       <- tryCatch(anova(fit_main_cc, fit, test = lrtest_str), error = function(e) NULL)
        p_inter_cc   <- if (!is.null(lrt_cc)) lrt_cc[2, lrtest_label] else NA_real_

        betas <- coef(fit)
        v_cov <- vcov(fit)
        ci_z  <- qnorm(1 - (1 - opts$ciWidth/100)/2)
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)

        for (j in seq_along(cov_levels)) {
          cl <- cov_levels[j]
          tbl <- arr$get(key = paste0(int_lbl, ": ", as.character(cl)))
          
          # Stats for the specific stratum
          mask_k <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
          st <- private$.compute_stats(geno_labels, snp_char[mask_k], response[mask_k], response_type)

          for (i in seq_along(geno_labels)) {
            gl <- geno_labels[i]
            
            # Identify active terms in the model for this specific cell
            # Reference Cell (SNP ref & Covar ref): Effect is 0 (or OR 1)
            if (i == 1 && j == 1) {
              tbl$addRow(rowKey = paste0(mdl, i), values = list(
                genotype = gl, stat0 = st$s0[i], stat1 = st$s1[i],
                effect = if (response_type == "binary") 1.0 else 0.0, ciLow = "", ciHigh = "", pval = ""
              ))
              next
            }

            # Map cell to model terms (R factor naming: interaction_varLevelName)
            term_snp   <- paste0("snp", gl)
            term_cov   <- paste0("interaction_var", cl)
            term_inter <- paste0("snp", gl, ":interaction_var", cl)
            
            active_terms <- c(
              if (i > 1) term_snp,
              if (j > 1) term_cov,
              if (i > 1 && j > 1) term_inter
            )
            active_terms <- active_terms[active_terms %in% names(betas)]

            # Calculate combined effect: Beta_total = sum(betas_active)
            combined_beta <- sum(betas[active_terms])
            combined_se   <- sqrt(sum(v_cov[active_terms, active_terms]))
            
            z_val   <- combined_beta / combined_se
            p_val   <- 2 * (1 - pnorm(abs(z_val)))
            lo_beta <- combined_beta - ci_z * combined_se
            hi_beta <- combined_beta + ci_z * combined_se

            tbl$addRow(rowKey = paste0(mdl, i), values = list(
              genotype = gl, stat0 = st$s0[i], stat1 = st$s1[i],
              effect = if (response_type == "binary") exp(combined_beta) else combined_beta,
              ciLow  = if (response_type == "binary") exp(lo_beta) else lo_beta,
              ciHigh = if (response_type == "binary") exp(hi_beta) else hi_beta,
              pval   = p_val
            ))
          }
          tbl$getColumn("stat1")$setVisible(response_type == "binary")
        }
        note_txt <- paste0("The reference group is <b>", interaction_var, ": ", cov_levels[1], " and ", snp_lbl, ": ", geno_labels[1], "</b>")
        tbl$setNote(note = note_txt, key = "interCrossClass")
        note_pval <- paste0("Interaction p-value: ", format.pval(p_inter_cc, digits = 3, eps = 0.001))
        tbl$setNote(note = note_pval, key = "interCrossClassPval")
      }
    }
  )
)
