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
      run_subpop         <- isTRUE(opts$subpop)

      # ── Validate SNPs ────────────────────────────────────────────
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

      # ── Response required ────────────────────────────────────────
      if (is.null(response_var) || response_var == "") {
        self$results$validationMsg$setContent(
          "<p style='color:red;'>A response variable is required for association analysis.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      # ── Prepare response / covariates ────────────────────────────
      response_raw  <- data[[response_var]]
      response_type <- detect_response_type(response_raw, opts$responseType)
      response      <- prepare_response(response_raw, response_type)
      cov_df        <- prepare_covariates(data, covariate_vars)

      if (run_subpop && response_type == "quantitative") run_subpop <- FALSE

      if (run_snpInteraction && length(covariate_vars) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>SNP \u00D7 covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }

      # ── Complete-case mask ───────────────────────────────────────
      n_rows        <- nrow(data)
      complete_mask <- rep(TRUE, n_rows)
      complete_mask <- complete_mask & !is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0)
        complete_mask <- complete_mask & complete.cases(cov_df)

      # ── Per-SNP association ──────────────────────────────────────
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
        cov_df_cc   <- if (!is.null(cov_df)) cov_df[snp_complete_mask, , drop=FALSE] else NULL
        if (is.null(geno_obj_cc)) next

        item <- arr$get(key = snp_nm)
        ref  <- get_ref_genotype(geno_obj_cc, user_levels)

        item$typingRate$setContent(sprintf(
          "<b>Typed samples:</b> %d / %d (%.1f%%)",
          sum(snp_complete_mask), n_rows, sum(snp_complete_mask)/n_rows*100))

        if (run_snpAssoc)
          private$.fill_assoc(item$assocTable, snp_raw_cc, ref, response_cc,
                              cov_df_cc, response_type, opts,
                              n_miss = n_miss_assoc, user_levels = user_levels)

        if (run_snpInteraction && !is.null(cov_df_cc) && ncol(cov_df_cc) >= 1)
          private$.fill_interaction(item$interactionTable, snp_raw_cc, ref,
                                    response_cc, cov_df_cc, names(cov_df_cc)[1],
                                    response_type, opts, user_levels = user_levels)
      }
    },

    # ── Association table ────────────────────────────────────────────────────
    .fill_assoc = function(tbl, snp_raw, ref, response, cov_df,
                           response_type, opts, n_miss = 0L,
                           user_levels = NULL) {
      if (response_type == "binary") {
        resp_clean <- response[!is.na(response)]
        if (length(unique(resp_clean)) != 2) {
          tbl$setNote(key="response_error",
                      note="Binary response requires exactly 2 categories.")
          return()
        }
      }

      tbl$getColumn("effect")$setTitle(if (response_type=="binary") "OR" else "\u03B2")

      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        cov_names <- sapply(names(cov_df), function(x) attr(self$data[[x]],"label") %||% x)
        note_txt  <- paste0("Model adjusted for: ", paste(cov_names, collapse=", "))
        if (!is.na(n_miss) && n_miss > 0)
          note_txt <- paste0(note_txt, ".  ", n_miss, " observation(s) excluded.")
        tbl$setNote(note=note_txt, key="covariates")
      } else if (!is.na(n_miss) && n_miss > 0) {
        tbl$setNote(note=paste0(n_miss," observation(s) excluded."), key="covariates")
      } else {
        tbl$setNote(note=NULL, key="covariates")
      }

      models <- c(
        if (opts$modelCodominant)   "codominant",
        if (opts$modelDominant)     "dominant",
        if (opts$modelRecessive)    "recessive",
        if (opts$modelOverdominant) "overdominant",
        if (opts$modelLogAdditive)  "logadditive"
      )
      model_labels <- c(codominant="Codominant", dominant="Dominant",
                        recessive="Recessive",   overdominant="Overdominant",
                        logadditive="Log-additive")

      if (isTRUE(opts$showAIC))
        tbl$addColumn(name="AIC", title="AIC", type="number", format="zto,dp=2")

      row_key <- 0L
      for (mdl in models) {
        snp_enc  <- encode_model(as.character(snp_raw), ref, mdl, user_levels)
        res_list <- fit_model(snp_enc, response, cov_df, mdl, response_type, opts$ciWidth)
        if (is.null(res_list)) next

        if (mdl == "codominant" && length(res_list) > 0) {
          gp <- res_list[[1]]$global_p
          if (!is.na(gp))
            tbl$setNote(note=paste0("Codominant model: LRT P = ", format.pval(gp, digits=3)),
                        key="lrt")
        }
        first_row <- TRUE
        for (res in res_list) {
          row_key <- row_key + 1L
          vals <- list(model=if(first_row) model_labels[mdl] else "",
                       comparison=res$comparison, effect=res$effect,
                       ciLow=res$ci_low, ciHigh=res$ci_high, pval=res$pval)
          if (isTRUE(opts$showAIC))
            vals[["AIC"]] <- if (first_row && !is.nan(res$aic)) round(res$aic,2) else ""
          tbl$addRow(rowKey=as.character(row_key), values=vals)
          first_row <- FALSE
        }
      }
    },

    # ── Interaction table ────────────────────────────────────────────────────
    .fill_interaction = function(tbl, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  user_levels = NULL) {
      tbl$getColumn("effect")$setTitle(if (response_type=="binary") "OR" else "\u03B2")
      adj_vars   <- setdiff(names(cov_df), interaction_var)
      note_parts <- paste0("Interaction covariate: ", interaction_var)
      if (length(adj_vars) > 0)
        note_parts <- paste0(note_parts, ". Adjusted for: ", paste(adj_vars, collapse=", "))
      tbl$setNote(note=note_parts, key="intcov")

      snp_enc_tmp   <- encode_model(as.character(snp_raw), ref, "logadditive", user_levels)
      complete_full <- !is.na(response) & !is.na(snp_enc_tmp) & complete.cases(cov_df)
      n_miss        <- length(response) - sum(complete_full)
      if (n_miss > 0)
        tbl$setNote(note=paste0(n_miss," observation(s) excluded."), key="missing_cov")
      else
        tbl$setNote(note=NULL, key="missing_cov")

      models <- c(
        if (opts$modelCodominant)   "codominant",
        if (opts$modelDominant)     "dominant",
        if (opts$modelRecessive)    "recessive",
        if (opts$modelOverdominant) "overdominant",
        if (opts$modelLogAdditive)  "logadditive"
      )
      model_labels <- c(codominant="Codominant", dominant="Dominant",
                        recessive="Recessive",   overdominant="Overdominant",
                        logadditive="Log-additive")

      if (isTRUE(opts$showAIC))
        tbl$addColumn(name="AIC", title="AIC", type="number", format="zto,dp=2")

      row_key <- 0L
      for (mdl in models) {
        snp_enc  <- encode_model(as.character(snp_raw), ref, mdl, user_levels)
        res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                          interaction_var, mdl, response_type, opts$ciWidth)
        if (is.null(res_list)) next
        first_row <- TRUE; first_inter <- TRUE
        for (res in res_list) {
          row_key  <- row_key + 1L
          is_inter <- !is.na(res$pval_interaction)
          vals <- list(model=if(first_row) model_labels[mdl] else "",
                       term=res$term, effect=res$effect,
                       ciLow=res$ci_low, ciHigh=res$ci_high, pval=res$pval,
                       pvalInteraction=if(is_inter && first_inter) res$pval_interaction else "")
          if (isTRUE(opts$showAIC))
            vals[["AIC"]] <- if (first_row && !is.nan(res$aic)) round(res$aic,2) else ""
          tbl$addRow(rowKey=as.character(row_key), values=vals)
          first_row <- FALSE
          if (is_inter) first_inter <- FALSE
        }
      }
    }
  )
)
