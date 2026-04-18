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
        self$results$validationMsg$setContent('')
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

        snp_complete_mask <- complete_mask & !is.na(snp_raw)
        n_miss_assoc      <- n_rows - sum(snp_complete_mask)

        snp_raw_cc  <- snp_raw[snp_complete_mask]
        geno_obj_cc <- parse_genotype(snp_raw_cc, user_levels)
        response_cc <- response[snp_complete_mask]
        cov_df_cc   <- if (!is.null(cov_df)) cov_df[snp_complete_mask, , drop = FALSE] else NULL
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

          if (isTRUE(opts$showStratByResponse) && response_type == "binary")
            private$.fill_strat_by_response(
              item$stratByResponse, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels, response_raw, snp_nm)

          if (isTRUE(opts$showStratByGenotype))
            private$.fill_strat_by_genotype(
              item$stratByGenotype, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels, response_raw, snp_nm)
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
        lv      <- levels(as.factor(response))
        n_total <- sum(!is.na(snp_char) & !is.na(response))
        stats0  <- character(length(geno_labels))
        stats1  <- character(length(geno_labels))
        for (i in seq_along(geno_labels)) {
          mask <- snp_char %in% split_genos(geno_labels[i]) & !is.na(response)
          n0   <- sum(mask & response == lv[1])
          n1   <- sum(mask & response == lv[2])
          if ((n0 + n1) == 0) { stats0[i] <- "---"; stats1[i] <- "---"; next }
          stats0[i] <- sprintf("%d (%.1f%%)", n0, n0 / n_total * 100)
          stats1[i] <- sprintf("%d (%.1f%%)", n1, n1 / n_total * 100)
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
          res <- res_list[[1]]; row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model = model_labels[mdl], genotype = "---", stat0 = "---", stat1 = "",
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

      tbl$setTitle(paste0(snp_lbl, " \u00D7 ", int_lbl, " interaction"))

      if (length(adj_vars) > 0){
        note_parts <- paste0(". Adjusted for: ", paste(sapply(adj_vars, function(x)
                               attr(self$data[[x]], "label") %||% x), collapse = ", "))
        tbl$setNote(note = note_parts, key = "intcov")
      }

      # report missings (not working previous filter)
      snp_enc_tmp   <- encode_model(as.character(snp_raw), ref, "logadditive", user_levels)
      complete_full <- !is.na(response) & !is.na(snp_enc_tmp) & complete.cases(cov_df)
      n_miss        <- length(response) - sum(complete_full)
      if (n_miss > 0)
        tbl$setNote(note = paste0(n_miss, " observation(s) excluded."), key = "missing_cov")
      else
        tbl$setNote(note = NULL, key = "missing_cov")

      
      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))

      # now only 1 model in int_models, but keep loop structure for future multi-select
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      row_key <- 0L
      for (mdl in int_models) {
        snp_enc  <- encode_model(as.character(snp_raw), ref, mdl, user_levels)
        res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                          interaction_var, mdl, response_type, opts$ciWidth)
        if (is.null(res_list)) next
        first_row <- TRUE; first_inter <- TRUE
        for (res in res_list) {
          row_key  <- row_key + 1L
          is_inter <- !is.na(res$pval_interaction)
          vals <- list(
            model  = if (first_row) model_labels[mdl] else "",
            term   = res$term,
            effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval   = res$pval,
            pvalInteraction = if (is_inter && first_inter) res$pval_interaction else "")
          if (isTRUE(opts$showAIC))
            vals[["AIC"]] <- if (first_row && !is.nan(res$aic)) round(res$aic, 2) else ""
          tbl$addRow(rowKey = as.character(row_key), values = vals)
          first_row <- FALSE
          if (is_inter) first_inter <- FALSE
        }
      }
    },

    # ── Stratified by response (binary: one table per response level) ─────────
    .fill_strat_by_response = function(arr, snp_raw, ref, response, cov_df,
                                        interaction_var, response_type, opts,
                                        int_models, user_levels = NULL, response_raw, snp_lbl) {
      
      # Factorize to get internal levels and user-facing labels
      resp_factor <- as.factor(response)
      lv          <- levels(resp_factor)
      
      # Get the actual labels (e.g., "Control", "Case")
      resp_raw_factor <- as.factor(response_raw)
      lv_names        <- levels(resp_raw_factor)
      
      resp_lbl <- attr(self$data[[self$options$response]], "label") %||% self$options$response
      snp_char <- as.character(snp_raw)

      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      all_genos <- c(ref, setdiff(
        if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])),
        ref))

      adj_cov_df <- if (!is.null(cov_df) && ncol(cov_df) > 0) cov_df else NULL

      # Loop using index so we can match lv with lv_names
      for (i_lv in seq_along(lv)) {
        current_lv    <- lv[i_lv]
        current_label <- lv_names[i_lv]
        
        # CHANGE: New key format "Variable: Label"
        key_k <- paste0(resp_lbl, ": ", current_label)
        
        if (is.null(tryCatch(arr$get(key = key_k), error = function(e) NULL)))
          arr$addItem(key = key_k)
        
        tbl <- arr$get(key = key_k)
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
        tbl$getColumn("genotype")$setTitle(snp_lbl)

        # Subset to this response stratum
        mask_k    <- !is.na(response) & response == current_lv & !is.na(snp_raw)
        if (!is.null(adj_cov_df) && ncol(adj_cov_df) > 0)
          mask_k <- mask_k & complete.cases(adj_cov_df)

        snp_k   <- snp_raw[mask_k]
        snp_c_k <- snp_char[mask_k]
        
        # Determine the interaction variable's effect within this stratum
        int_response_k <- cov_df[[interaction_var]][mask_k]
        int_resp_type  <- detect_response_type(int_response_k, opts$responseType)
        adj_cov_k      <- if (!is.null(adj_cov_df)) adj_cov_df[mask_k, setdiff(names(adj_cov_df), interaction_var), drop = FALSE] else NULL
        if (!is.null(adj_cov_k) && ncol(adj_cov_k) == 0) adj_cov_k <- NULL

        row_key <- 0L
        for (mdl in int_models) {
          snp_enc  <- encode_model(snp_c_k, ref, mdl, user_levels)
          res_list <- fit_model(snp_enc, int_response_k, adj_cov_k,
                                mdl, int_resp_type, opts$ciWidth)
          if (is.null(res_list)) next

          geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)

          if (mdl == "logadditive") {
            res <- res_list[[1]]; row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              genotype = paste0(model_labels[mdl], " (per allele)" ),
              n        = sum(!is.na(snp_enc) & !is.na(int_response_k)),
              effect   = res$effect, ciLow = res$ci_low,
              ciHigh   = res$ci_high, pval  = res$pval))
            next
          }

          # Reference row
          ref_mask <- snp_c_k %in% private$.split_genos(geno_labels[1]) & !is.na(int_response_k)
          row_key  <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            genotype = geno_labels[1],
            n        = sum(ref_mask),
            effect   = if (int_resp_type == "binary") 1.0 else 0.0,
            ciLow    = '', ciHigh = '', pval = ''))

          for (i in seq_along(res_list)) {
            res <- res_list[[i]]
            gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else res$comparison
            gl_mask  <- snp_c_k %in% private$.split_genos(gl) & !is.na(int_response_k)
            row_key  <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              genotype = gl, 
              n        = sum(gl_mask),
              effect   = res$effect, ciLow = res$ci_low,
              ciHigh   = res$ci_high, pval  = res$pval))
          }
        }
      }
    },

    # ── Stratified by genotype (one table per genotype group) ─────────────────
    # For each genotype group g: fit response ~ covariate(s) in the subset
    # where SNP == g, giving the effect of the covariate within that genotype.

    # ── Stratified by genotype (one table per genotype group) ─────────────────
    .fill_strat_by_genotype = function(arr, snp_raw, ref, response, cov_df,
                                       interaction_var, response_type, opts,
                                       int_models, user_levels = NULL, response_raw, snp_lbl) {
      snp_char <- as.character(snp_raw)
      all_genos <- c(ref, setdiff(
        if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])),
        ref))

      # Use the interaction variable name or its label as the column header for levels
      int_lbl <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      adj_vars <- setdiff(names(cov_df), interaction_var)
      int_response <- cov_df[[interaction_var]]
      int_resp_type <- detect_response_type(int_response, opts$responseType)

      for (mdl in int_models) {
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
        if (mdl == "logadditive") next

        for (gl in geno_labels) {
          key_g <- paste0(snp_lbl, ": ", gl)
#          key_g <- gl 
          
          if (is.null(tryCatch(arr$get(key = key_g), error = function(e) NULL)))
            arr$addItem(key = key_g)
          
          tbl <- arr$get(key = key_g)
          
          tbl$getColumn("level")$setTitle(int_lbl)
          
          tbl$getColumn("effect")$setTitle(
            if (response_type == "binary") "OR" else "\u03B2")

          # Subset to this genotype group
          genos_in_group <- private$.split_genos(gl)
          mask_g  <- snp_char %in% genos_in_group & !is.na(response) & !is.na(int_response)
          if (!is.null(cov_df) && ncol(cov_df) > 0)
            mask_g <- mask_g & complete.cases(cov_df)

          n_g          <- sum(mask_g)
          response_g   <- response[mask_g]
          int_resp_g   <- int_response[mask_g]
          adj_cov_g    <- if (length(adj_vars) > 0)
                            cov_df[mask_g, adj_vars, drop = FALSE] else NULL

          if (int_resp_type == "binary") {
            lv_int  <- levels(as.factor(int_resp_g))
            row_key <- 0L
            
            # 1. Reference Row
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              covariate = "", # Removed "Covariate" text
              level     = as.character(lv_int[1]),
              n         = sum(int_resp_g == lv_int[1]),
              effect    = if (response_type == "binary") 1.0 else 0.0,
              ciLow     = '', ciHigh = '', pval = ''
            ))

            # 2. Comparison Rows
            for (i in 2:length(lv_int)) {
              sub_mask    <- int_resp_g %in% c(lv_int[1], lv_int[i])
              snp_enc_g   <- as.integer(int_resp_g[sub_mask] == lv_int[i])
              
              res_list <- fit_model(snp_enc_g, response_g[sub_mask], 
                                   if(!is.null(adj_cov_g)) adj_cov_g[sub_mask, , drop=FALSE] else NULL,
                                   "logadditive", response_type, opts$ciWidth)
              
              if (!is.null(res_list)) {
                res <- res_list[[1]]
                row_key <- row_key + 1L
                tbl$addRow(rowKey = as.character(row_key), values = list(
                  covariate = "", 
                  level     = as.character(lv_int[i]),
                  n         = sum(int_resp_g == lv_int[i]),
                  effect    = res$effect, 
                  ciLow     = res$ci_low,
                  ciHigh    = res$ci_high, 
                  pval      = res$pval))
              }
            }
          } else {
            # Quantitative interaction: single row
            row_key   <- 1L
            snp_enc_g <- int_resp_g  
            res_list  <- fit_model(snp_enc_g, response_g, adj_cov_g,
                                   "logadditive", response_type, opts$ciWidth)
            if (!is.null(res_list)) {
              res <- res_list[[1]]
              tbl$addRow(rowKey = as.character(row_key), values = list(
                covariate = "", 
                level     = "continuous",
                n         = n_g,
                effect    = res$effect, ciLow = res$ci_low,
                ciHigh    = res$ci_high, pval  = res$pval))
            }
          }
        }
      }
    }
  )
)
