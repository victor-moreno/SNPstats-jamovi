#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype
source("R/snp_helpers.R")

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
                              n_miss = n_miss_assoc, user_levels = user_levels)

        if (run_snpInteraction && !is.null(cov_df_cc) && ncol(cov_df_cc) >= 1) {
          interaction_var <- names(cov_df_cc)[1]

          # Vector of models for interaction — currently single element from
          # the dropdown.  When interactionModel becomes multi-select, replace
          # this one line with: private$.get_interaction_models(opts)
          int_models <- opts$interactionModel

          private$.fill_interaction(
            item$interactionTable, snp_raw_cc, ref,
            response_cc, cov_df_cc, interaction_var,
            response_type, opts, int_models, user_levels = user_levels)

          if (isTRUE(opts$showStratByResponse) && response_type == "binary")
            private$.fill_strat_by_response(
              item$stratByResponse, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels = user_levels)

          if (isTRUE(opts$showStratByGenotype))
            private$.fill_strat_by_genotype(
              item$stratByGenotype, snp_raw_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              response_type, opts, int_models, user_levels = user_levels)
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
                           response_type, opts, n_miss = 0L, user_levels = NULL) {
      if (response_type == "binary") {
        if (length(unique(response[!is.na(response)])) != 2) {
          tbl$setNote(key = "response_error",
                      note = "Binary response requires exactly 2 categories.")
          return()
        }
      }

      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")

      if (response_type == "binary") {
        lv       <- levels(as.factor(response))
        resp_lbl <- attr(self$data[[self$options$response]], "label") %||% self$options$response
        tbl$getColumn("stat0")$setTitle(paste0(resp_lbl, "=", lv[1]))
        tbl$getColumn("stat1")$setTitle(paste0(resp_lbl, "=", lv[2]))
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

        if (mdl == "codominant") {
          gp <- res_list[[1]]$global_p
          if (!is.null(gp) && !is.na(gp))
            tbl$setNote(
              note = paste0("Codominant model: LRT P = ", format.pval(gp, digits = 3)),
              key  = "lrt")
        }

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
        row_key <- row_key + 1L
        tbl$addRow(rowKey = as.character(row_key), values = list(
          model    = model_labels[mdl],
          genotype = geno_labels[1],
          stat0    = st$s0[1], stat1 = st$s1[1],
          effect   = if (response_type == "binary") 1. else 0.,
          ciLow = '', ciHigh = '', pval = '',
          AIC = aic_val, BIC = bic_val))

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
                                  int_models, user_levels = NULL) {
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")

      adj_vars   <- setdiff(names(cov_df), interaction_var)
      int_lbl    <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      note_parts <- paste0("Interaction covariate: ", int_lbl)
      if (length(adj_vars) > 0)
        note_parts <- paste0(note_parts, ". Adjusted for: ",
                             paste(sapply(adj_vars, function(x)
                               attr(self$data[[x]], "label") %||% x), collapse = ", "))
      tbl$setNote(note = note_parts, key = "intcov")

      snp_enc_tmp   <- encode_model(as.character(snp_raw), ref, "logadditive", user_levels)
      complete_full <- !is.na(response) & !is.na(snp_enc_tmp) & complete.cases(cov_df)
      n_miss        <- length(response) - sum(complete_full)
      if (n_miss > 0)
        tbl$setNote(note = paste0(n_miss, " observation(s) excluded."), key = "missing_cov")
      else
        tbl$setNote(note = NULL, key = "missing_cov")

      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))

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
    # For each response level k: fit SNP ~ covariates in the subset where
    # response == k, giving OR/β of the SNP within that response stratum.
    .fill_strat_by_response = function(arr, snp_raw, ref, response, cov_df,
                                        interaction_var, response_type, opts,
                                        int_models, user_levels = NULL) {
      lv       <- levels(as.factor(response))
      resp_lbl <- attr(self$data[[self$options$response]], "label") %||% self$options$response
      snp_char <- as.character(snp_raw)

      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      all_genos <- c(ref, setdiff(
        if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])),
        ref))

      # Covariates to adjust for within stratum (exclude the interaction covariate
      # because we are now stratifying by response, and interaction_var is a covariate)
      adj_cov_df <- if (!is.null(cov_df) && ncol(cov_df) > 0) cov_df else NULL

      for (lv_k in lv) {
        key_k <- paste0(resp_lbl, "=", lv_k)
        # Add sub-table for this level if not yet present
        if (is.null(tryCatch(arr$get(key = key_k), error = function(e) NULL)))
          arr$addItem(key = key_k)
        tbl <- arr$get(key = key_k)
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")

        # Subset to this response stratum
        mask_k    <- !is.na(response) & response == lv_k & !is.na(snp_raw)
        if (!is.null(adj_cov_df) && ncol(adj_cov_df) > 0)
          mask_k <- mask_k & complete.cases(adj_cov_df)

        snp_k   <- snp_raw[mask_k]
        snp_c_k <- snp_char[mask_k]
        cov_k   <- if (!is.null(adj_cov_df)) adj_cov_df[mask_k, , drop = FALSE] else NULL

        # For binary response stratified by level we have a constant response —
        # use quantitative encoding (linear) or skip if degenerate.
        # More usefully: fit covariate association with SNP within stratum.
        # Standard interpretation: within cases / within controls, how does
        # the covariate (interaction_var) associate with the SNP?
        # We use the interaction_var as the new "response" within each stratum.
        int_response_k <- cov_df[[interaction_var]][mask_k]
        int_resp_type  <- detect_response_type(int_response_k, opts$responseType)
        adj_cov_k      <- if (!is.null(cov_k)) cov_k[, setdiff(names(cov_k), interaction_var), drop = FALSE] else NULL
        if (!is.null(adj_cov_k) && ncol(adj_cov_k) == 0) adj_cov_k <- NULL

        int_lbl <- attr(self$data[[interaction_var]], "label") %||% interaction_var
        tbl$setNote(note = paste0("Covariate: ", int_lbl, ".  n = ", sum(mask_k)),
                    key = "stratum_note")

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
              genotype = paste0(model_labels[mdl], " (---)" ),
              n        = sum(!is.na(snp_enc) & !is.na(int_response_k)),
              effect   = res$effect, ciLow = res$ci_low,
              ciHigh   = res$ci_high, pval  = res$pval))
            next
          }

          # Reference row
          ref_mask <- snp_c_k %in% private$.split_genos(geno_labels[1]) & !is.na(int_response_k)
          row_key  <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            genotype = paste0(model_labels[mdl], ": ", geno_labels[1], " (ref)"),
            n        = sum(ref_mask),
            effect   = if (int_resp_type == "binary") 1. else 0.,
            ciLow    = '', ciHigh = '', pval = ''))

          for (i in seq_along(res_list)) {
            res <- res_list[[i]]
            gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else res$comparison
            gl_mask  <- snp_c_k %in% private$.split_genos(gl) & !is.na(int_response_k)
            row_key  <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              genotype = if (i == 1) paste0(model_labels[mdl], ": ", gl) else gl,
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
    .fill_strat_by_genotype = function(arr, snp_raw, ref, response, cov_df,
                                        interaction_var, response_type, opts,
                                        int_models, user_levels = NULL) {
      snp_char <- as.character(snp_raw)
      all_genos <- c(ref, setdiff(
        if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)])),
        ref))

      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")

      int_lbl  <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      adj_vars <- setdiff(names(cov_df), interaction_var)
      int_response <- cov_df[[interaction_var]]
      int_resp_type <- detect_response_type(int_response, opts$responseType)

      for (mdl in int_models) {
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
        # For log-additive there are no discrete groups to stratify by
        if (mdl == "logadditive") next

        for (gl in geno_labels) {
          key_g <- paste0(model_labels[mdl], ": ", gl)
          if (is.null(tryCatch(arr$get(key = key_g), error = function(e) NULL)))
            arr$addItem(key = key_g)
          tbl <- arr$get(key = key_g)
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

          tbl$setNote(note = paste0("Covariate: ", int_lbl, ".  n = ", n_g),
                      key  = "stratum_note")

          if (n_g < 5) {
            tbl$setNote(note = "Too few observations in this genotype group.",
                        key  = "sparse")
            next
          }

          # Fit response ~ interaction_var (± other covariates) in this genotype stratum
          int_enc  <- if (int_resp_type == "binary")
                        as.integer(as.factor(int_resp_g)) - 1L
                      else
                        int_resp_g

          # If interaction_var is categorical, loop over its levels
          if (int_resp_type == "binary") {
            lv_int    <- levels(as.factor(int_resp_g))
            row_key   <- 0L

            # Encode cov as SNP-like predictor for fit_model
            snp_enc_g <- encode_model(as.character(as.factor(int_resp_g)),
                                      lv_int[1], "logadditive", NULL)
            res_list  <- fit_model(snp_enc_g, response_g, adj_cov_g,
                                   "logadditive", response_type, opts$ciWidth)
            if (!is.null(res_list)) {
              res <- res_list[[1]]; row_key <- row_key + 1L
              tbl$addRow(rowKey = as.character(row_key), values = list(
                covariate = int_lbl, level = paste(lv_int, collapse = " vs "),
                n         = n_g,
                effect    = res$effect, ciLow = res$ci_low,
                ciHigh    = res$ci_high, pval  = res$pval))
            }
          } else {
            # Quantitative interaction covariate: single linear row
            row_key  <- 0L
            snp_enc_g <- int_resp_g   # use raw numeric as predictor
            res_list  <- fit_model(snp_enc_g, response_g, adj_cov_g,
                                   "logadditive", response_type, opts$ciWidth)
            if (!is.null(res_list)) {
              res <- res_list[[1]]; row_key <- row_key + 1L
              tbl$addRow(rowKey = as.character(row_key), values = list(
                covariate = int_lbl, level = "continuous",
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
