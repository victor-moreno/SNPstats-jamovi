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

      # ── Column titles for stat columns ───────────────────────────────────
      if (response_type == "binary") {
        resp_fac    <- as.factor(response)
        resp_levels <- levels(resp_fac)
        resp_lbl    <- attr(self$data[[self$options$response]], "label") %||% self$options$response
        lbl0 <- paste0(resp_lbl, "=", resp_levels[1])
        lbl1 <- paste0(resp_lbl, "=", resp_levels[2])
        tbl$getColumn("stat0")$setTitle(lbl0)
        tbl$getColumn("stat1")$setTitle(lbl1)
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(TRUE)
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)")
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }

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

      # ── AIC / BIC column visibility ──────────────────────────────────────
      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))
      tbl$getColumn("BIC")$setVisible(isTRUE(opts$showAIC))

      # ── BIC from AIC: BIC = AIC + df*(log(n) - 2) ───────────────────────
      # df = 1 (intercept) + n_covariates + SNP parameters for model
      n_fit  <- sum(!is.na(as.character(snp_raw)) & !is.na(response) &
                      (if (!is.null(cov_df) && ncol(cov_df) > 0) complete.cases(cov_df) else TRUE))
      n_cov  <- if (!is.null(cov_df)) ncol(cov_df) else 0L
      snp_df <- c(codominant=2L, dominant=1L, recessive=1L, overdominant=1L, logadditive=1L)
      bic_from_aic <- function(aic_val, mdl) {
        if (is.null(aic_val) || is.na(aic_val) || is.nan(aic_val)) return(NA_real_)
        df <- 1L + n_cov + snp_df[[mdl]]
        round(aic_val + df * (log(n_fit) - 2), 2)
      }

      # ── Genotype levels (ordered as ref, het, hom_alt) ──────────────────
      snp_char  <- as.character(snp_raw)
      all_genos <- if (!is.null(user_levels)) user_levels else sort(unique(snp_char[!is.na(snp_char)]))
      all_genos <- c(ref, setdiff(all_genos, ref))

      # ── Helper: genotype labels per model ───────────────────────────────
      geno_labels_for_model <- function(mdl, all_genos, ref) {
        if (mdl == "codominant" || mdl == "logadditive") return(all_genos)
        het  <- all_genos[all_genos != ref & all_genos != all_genos[length(all_genos)]]
        hom2 <- all_genos[length(all_genos)]
        if (length(het) == 0) het <- hom2
        if (mdl == "dominant")     return(c(ref, paste(c(het, hom2), collapse="-")))
        if (mdl == "recessive")    return(c(paste(c(ref, het), collapse="-"), hom2))
        if (mdl == "overdominant") return(c(paste(c(ref, hom2), collapse="-"), het))
        all_genos
      }

      # ── Helper: compute N(%) or mean(SD) per genotype group ─────────────
      compute_stats <- function(geno_labels, snp_char, response, response_type) {
        split_genos <- function(gl)
          unlist(strsplit(gl, "(?<=[A-Za-z0-9*])-(?=[A-Za-z0-9*])", perl=TRUE))

        if (response_type == "binary") {
          resp_fac <- as.factor(response)
          lv       <- levels(resp_fac)
          n_total  <- sum(!is.na(snp_char) & !is.na(response))
          stats0   <- character(length(geno_labels))
          stats1   <- character(length(geno_labels))
          for (i in seq_along(geno_labels)) {
            mask <- snp_char %in% split_genos(geno_labels[i]) & !is.na(response)
            n0   <- sum(mask & response == lv[1])
            n1   <- sum(mask & response == lv[2])
            if ((n0 + n1) == 0) { stats0[i] <- "---"; stats1[i] <- "---"; next }
            stats0[i] <- sprintf("%d (%.1f%%)", n0, n0/n_total*100)
            stats1[i] <- sprintf("%d (%.1f%%)", n1, n1/n_total*100)
          }
          list(s0=stats0, s1=stats1)
        } else {
          stats0 <- character(length(geno_labels))
          for (i in seq_along(geno_labels)) {
            mask <- snp_char %in% split_genos(geno_labels[i]) & !is.na(response)
            vals <- response[mask]
            if (length(vals) == 0) { stats0[i] <- "---"; next }
            stats0[i] <- sprintf("%.2f (%.2f)", mean(vals), sd(vals))
          }
          list(s0=stats0, s1=rep("", length(geno_labels)))
        }
      }

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

        geno_labels <- geno_labels_for_model(mdl, all_genos, ref)
        st          <- compute_stats(geno_labels, snp_char, response, response_type)

        # AIC and BIC are model-level — shown once on first (reference) row
        aic_val <- {
          a <- res_list[[1]]$aic
          if (!is.null(a) && !is.na(a) && !is.nan(a)) round(a, 2) else NA_real_
        }
        bic_val <- bic_from_aic(aic_val, mdl)

        if (mdl == "logadditive") {
          res     <- res_list[[1]]
          row_key <- row_key + 1L
          tbl$addRow(rowKey=as.character(row_key), values=list(
            model    = model_labels[mdl],
            genotype = "---",
            stat0    = "---",
            stat1    = "",
            effect   = res$effect,
            ciLow    = res$ci_low,
            ciHigh   = res$ci_high,
            pval     = res$pval,
            AIC      = aic_val,
            BIC      = bic_val
          ))
          next
        }

        # ── Reference row (OR = 1 / β = 0) — AIC & BIC shown here ──────
        row_key <- row_key + 1L
        tbl$addRow(rowKey=as.character(row_key), values=list(
          model    = model_labels[mdl],
          genotype = geno_labels[1],
          stat0    = st$s0[1],
          stat1    = st$s1[1],
          effect   = if (response_type == "binary") '1.' else '0.',
          ciLow    = '',
          ciHigh   = '',
          pval     = '',
          AIC      = aic_val,
          BIC      = bic_val
        ))

        # ── Non-reference rows ───────────────────────────────────────────
        for (i in seq_along(res_list)) {
          res     <- res_list[[i]]
          gl      <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else res$comparison
          row_key <- row_key + 1L
          tbl$addRow(rowKey=as.character(row_key), values=list(
            model    = "",
            genotype = gl,
            stat0    = if ((i+1) <= length(st$s0)) st$s0[i+1] else "-",
            stat1    = if ((i+1) <= length(st$s1)) st$s1[i+1] else "",
            effect   = res$effect,
            ciLow    = res$ci_low,
            ciHigh   = res$ci_high,
            pval     = res$pval,
            AIC      = '',
            BIC      = ''
          ))
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

      if (isTRUE(opts$showAIC)) {
        tbl$addColumn(name="AIC", title="AIC", type="number", format="zto,dp=2")
#        tbl$addColumn(name="BIC", title="BIC", type="number", format="zto,dp=2")
      }

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
          if (isTRUE(opts$showAIC)) {
            vals[["AIC"]] <- if (first_row && !is.nan(res$aic)) round(res$aic,2) else ""
#            vals[["BIC"]] <- if (first_row && !is.nan(res$bic)) round(res$bic,2) else ""
          }
          tbl$addRow(rowKey=as.character(row_key), values=vals)
          first_row <- FALSE
          if (is_inter) first_inter <- FALSE
        }
      }
    }
  )
)
