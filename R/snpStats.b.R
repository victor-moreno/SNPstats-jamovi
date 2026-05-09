#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype allele HWE.exact LD
#' @importFrom haplo.stats setupGeno hapl.em haplo.glm haplo.glm.control
#' @import ggplot2

# ── Formula safety helper ──────────────────────────────────────────────────
# Wraps a column name in backticks so it is safe to interpolate into a
# formula string even when the name contains spaces, operators, or other
# special characters.  This is the same escaping that jmvcore::composeTerm()
# applies and matches the recommendation in the jamovi developer docs.
safe_term <- function(x) paste0("`", gsub("`", "\\`", x, fixed = TRUE), "`")

# Escape a vector of names and collapse to a single "+"-joined string,
# suitable for the RHS of a formula.
safe_rhs <- function(nms) paste(sapply(nms, safe_term), collapse = " + ")



# ══════════════════════════════════════════════════════════════════════════════
# snpStatsClass
# ══════════════════════════════════════════════════════════════════════════════

snpStatsClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "snpStatsClass",
  inherit = snpStatsBase,
  private = list(

    # ── Private state ───────────────────────────────────────────────────────
    .miss_cache = NULL,   # populated during descriptive run; read by .plotMissingness
    .ld_store   = NULL,   # LD heatmap data
    .ld_nms     = NULL,
    .ld_metric  = NULL,

    # ════════════════════════════════════════════════════════════════════════
    # .init
    # ════════════════════════════════════════════════════════════════════════
    .init = function() {
      private$.miss_cache <- list()

      # Always reset group visibility
      self$results$descGroup$covDescGroup$setVisible(FALSE)
      self$results$descGroup$snpSummaryTablesGroup$setVisible(FALSE)
      self$results$ldHaploGroup$ldGroup$setVisible(FALSE)
      self$results$ldHaploGroup$haploGroup$setVisible(FALSE)

      snp_names <- self$options$snps
      if (length(snp_names) == 0) return()

      # Descriptive per-SNP array
      desc_arr <- self$results$descGroup$descSnpResults
      for (nm in snp_names) desc_arr$addItem(key = nm)

      # Association per-SNP array
      opts         <- self$options
      assoc_arr    <- self$results$assocGroup$assocSnpResults
      int_models   <- private$.get_interaction_models(opts)
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      for (nm in snp_names) {
        assoc_arr$addItem(key = nm)
        if (isTRUE(opts$snpInteraction) && length(int_models) > 0) {
          int_arr <- assoc_arr$get(key = nm)$interactionResults
          for (mdl in int_models) int_arr$addItem(key = model_labels[[mdl]])
        }
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # .run  — runs all enabled analyses based on option flags
    # ════════════════════════════════════════════════════════════════════════
    .run = function() {
      data           <- self$data
      opts           <- self$options
      response_var   <- opts$response
      snp_vars       <- opts$snps
      covariate_vars <- opts$covariates

      # ── Shared validation ──────────────────────────────────────────────
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

      # ── Shared data preparation ─────────────────────────────────────────
      response_raw  <- if (!is.null(response_var) && response_var != "")
                         data[[response_var]] else NULL
      response_type <- detect_response_type(response_raw, opts$responseType)
      response      <- prepare_response(response_raw, response_type)
      cov_df        <- prepare_covariates(data, covariate_vars)
      if (is.null(cov_df) && !is.null(response_raw))
        cov_df <- data.frame(row.names = seq_len(nrow(data)))

      n_rows        <- nrow(data)
      complete_mask <- rep(TRUE, n_rows)
      if (!is.null(response))                        complete_mask <- complete_mask & !is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0)     complete_mask <- complete_mask & complete.cases(cov_df)

      # ── Always run descriptive (gated internally by its option flags) ───
      private$.run_descriptive(data, snp_vars, response_var, response_raw,
                               response_type, response, cov_df,
                               complete_mask, n_rows, opts)

      # ── Association: run only when at least one assoc option is ticked ──
      any_assoc <- isTRUE(opts$snpAssoc) || isTRUE(opts$snpInteraction)
      if (any_assoc) {
        if (is.null(response_var) || response_var == "") {
          self$results$validationMsg$setContent(
            "<p style='color:red;'>A response variable is required for association analysis.</p>")
          self$results$validationMsg$setVisible(TRUE)
        } else {
          private$.run_association(data, snp_vars, response_var, response_raw,
                                   response_type, response, cov_df,
                                   complete_mask, n_rows, opts)
        }
      }

      # ── LD/Haplotype: run only when at least one LD/haplo option is ticked
      any_ld <- isTRUE(opts$ldAnalysis) || isTRUE(opts$ldMatrix) || isTRUE(opts$ldPlot) ||
                isTRUE(opts$haploFreq)  || isTRUE(opts$haploAssoc) || isTRUE(opts$haploInteraction)
      if (any_ld) {
        if (length(snp_vars) < 2) {
          self$results$validationMsg$setContent(
            "<p style='color:red;'>LD and haplotype analyses require at least 2 SNPs.</p>")
          self$results$validationMsg$setVisible(TRUE)
        } else {
          private$.run_ldhaplo(data, snp_vars, response_var, response_raw,
                               response_type, response, cov_df,
                               complete_mask, n_rows, opts)
        }
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # .run_descriptive
    # ════════════════════════════════════════════════════════════════════════
    .run_descriptive = function(data, snp_vars, response_var, response_raw,
                                response_type, response, cov_df,
                                complete_mask, n_rows, opts) {
      run_snpSummary   <- isTRUE(opts$snpSummary)
      run_allFreq      <- isTRUE(opts$allFreq)
      run_genoFreq     <- isTRUE(opts$genoFreq)
      run_hweTest      <- isTRUE(opts$hweTest)
      run_subpop       <- isTRUE(opts$subpop)
      run_covDesc      <- isTRUE(opts$covDesc)
      run_showMissing  <- isTRUE(opts$showMissing)
      run_rmSnpMissing <- isTRUE(opts$rmSnpMissing)

      private$.miss_cache <- list()

      res <- self$results$descGroup

      # ── Validation messages ───────────────────────────────────────────
      res$validationMsgSNP$setVisible(FALSE)
      val2 <- validate_snp_vars(snp_vars, data)
      if (nchar(val2$bad_html) > 0) {
        res$validationMsgGeno$setContent(val2$bad_html)
        res$validationMsgGeno$setVisible(TRUE)
      } else {
        res$validationMsgGeno$setVisible(FALSE)
      }

      if (run_subpop && (is.null(response_raw) || response_type == "quantitative"))
        run_subpop <- FALSE

      # ── Covariate descriptives ────────────────────────────────────────
      if (run_covDesc && (!is.null(cov_df) || !is.null(response_raw))){
        res$covDescGroup$setVisible(TRUE)
        private$.run_cov_desc(cov_df, response_raw, response_type, run_subpop,
                              opts$response, data, snp_vars, run_rmSnpMissing)
      }
      # ── SNP summary table ─────────────────────────────────────────────
      if (run_snpSummary) {
        res$snpSummaryTablesGroup$setVisible(run_snpSummary && length(snp_vars) > 0)
        private$.fill_snp_summary(data, snp_vars, response_raw, response_type,
                                  run_subpop, cov_df)
      }
      
      # ── Per-SNP descriptives ──────────────────────────────────────────
      null_pat               <- "^0[/|>]0$|^00$"
      total_null_across_snps <- 0L
      arr <- res$descSnpResults

      for (snp_nm in snp_vars) {
        snp_raw_chr     <- as.character(data[[snp_nm]])
        n_null_replaced <- sum(!is.na(snp_raw_chr) &
                                 grepl(null_pat, snp_raw_chr, ignore.case = TRUE))
        total_null_across_snps <- total_null_across_snps + n_null_replaced
        # Extract factor level order BEFORE clean_null_alleles drops the factor class
        user_levels <- get_snp_level_order(data[[snp_nm]])
        snp_raw     <- clean_null_alleles(snp_raw_chr)
        geno_obj    <- parse_genotype(snp_raw, user_levels)

        if (is.null(geno_obj)) {
          if (all(is.na(snp_raw))) {
            item <- arr$get(key = snp_nm)
            n_total_eligible <- sum(complete_mask)
            item$typingRate$setContent(sprintf(
              "<b>Typed samples:</b> 0 / %d (0.0%%) &nbsp;&mdash;&nbsp; <span style='color:orange;'>all genotypes missing</span>",
              n_total_eligible))
          }
          next
        }

        snp_complete_mask <- complete_mask & !is.na(snp_raw)
        total_missing     <- sum(is.na(snp_raw) & complete_mask)
        n_total_eligible  <- sum(complete_mask)

        if (run_subpop && (response_type == "binary" || response_type == "categorical")) {
          resp_levels <- levels(as.factor(response_raw))
          n_miss_by_level <- sapply(resp_levels, function(lvl)
            sum(is.na(snp_raw) & complete_mask & !is.na(response_raw) & response_raw == lvl))
          names(n_miss_by_level) <- resp_levels
        } else {
          n_miss_by_level <- NULL
        }

        snp_raw_cc  <- snp_raw[snp_complete_mask]
        geno_obj_cc <- parse_genotype(snp_raw_cc, user_levels)
        response_cc <- if (!is.null(response))     response[snp_complete_mask]     else NULL
        resp_raw_cc <- if (!is.null(response_raw)) response_raw[snp_complete_mask] else NULL
        if (is.null(geno_obj_cc)) next

        item    <- arr$get(key = snp_nm)
        n_typed <- sum(snp_complete_mask)

        typing_html <- sprintf("<b>Typed samples:</b> %d / %d (%.1f%%)",
          n_typed, n_total_eligible,
          if (n_total_eligible > 0) n_typed / n_total_eligible * 100 else 0)
        if (total_missing > 0)
          typing_html <- paste0(typing_html, sprintf(
            " &nbsp;&mdash;&nbsp; <b>Missing SNP:</b> %d (%.1f%%)",
            total_missing,
            if (n_total_eligible > 0) total_missing / n_total_eligible * 100 else 0))
        item$typingRate$setContent(typing_html)

        private$.miss_cache[[snp_nm]] <- list(
          n_total_eligible = n_total_eligible,
          total_missing    = total_missing,
          n_miss_by_level  = n_miss_by_level)

        snp_summary_cc <- summary(geno_obj_cc)
        ref            <- get_ref_genotype(geno_obj_cc, user_levels)

        if (run_allFreq)
          private$.fill_allele_freq(item$allFreqTable, snp_summary_cc, snp_nm,
                                    resp_raw_cc, run_subpop, response_type, snp_raw_cc,
                                    run_showMissing, n_miss_by_level,
                                    n_total_eligible, total_missing,
                                    user_levels = user_levels)
        if (run_genoFreq)
          private$.fill_geno_freq(item$genoFreqTable, snp_summary_cc, ref,
                                  geno_obj_cc, response_cc, response_type,
                                  run_subpop, resp_raw_cc,
                                  run_showMissing, n_miss_by_level,
                                  n_total_eligible, total_missing,
                                  user_levels = user_levels)
        if (run_hweTest)
          private$.fill_hwe(item$hweTable, geno_obj_cc, snp_nm,
                            resp_raw_cc, run_subpop,
                            run_showMissing, n_miss_by_level, n_total_eligible, total_missing,
                            ref = ref, user_levels = user_levels)
      }

      # Missingness plot visibility
      res$missingnessPlot$setVisible(
        isTRUE(opts$showMissingnessPlot) && length(private$.miss_cache) > 0)

      # Null-allele note on summary table
      if (run_snpSummary) {
        tbl <- res$snpSummaryTablesGroup$snpSummaryTable
        tbl$setNote(
          key  = "null_allele",
          note = if (total_null_across_snps > 0)
            paste0(total_null_across_snps,
                   " genotype(s) coded as 0/0 were treated as missing (NA).")
          else NULL)
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # .run_association
    # ════════════════════════════════════════════════════════════════════════
    .run_cov_desc = function(cov_df, response_raw, response_type,
                             subpop = FALSE, response_var = NULL,
                             data = NULL, snp_vars = NULL, rm_snp_missing = FALSE) {
      tbl <- self$results$descGroup$covDescGroup$covDescTable
      has_covariates <- !is.null(cov_df) && ncol(cov_df) > 0
      n_removed_snp <- 0L
      if (isTRUE(rm_snp_missing) && !is.null(data) && !is.null(snp_vars) && length(snp_vars) > 0) {
        snp_mat <- as.data.frame(
          lapply(data[, snp_vars, drop = FALSE], function(col) clean_null_alleles(as.character(col))),
          stringsAsFactors = FALSE)
        snp_cc <- complete.cases(snp_mat)
        n_removed_snp <- sum(!snp_cc)
        if (n_removed_snp > 0) {
          cov_df <- cov_df[snp_cc, , drop=FALSE]
          if (!is.null(response_raw)) response_raw <- response_raw[snp_cc]
        }
      }
      has_response <- !is.null(response_raw) && !is.null(response_type)
      is_binary    <- has_response && response_type == "binary"
      is_cat_resp  <- has_response && response_type == "categorical"
      is_cont_resp <- has_response && response_type == "quantitative"
      do_strat     <- isTRUE(subpop) && (is_binary || is_cat_resp)
      valid_resp   <- if (has_response) !is.na(response_raw) else rep(TRUE, nrow(cov_df))
      if (do_strat) {
        grp_fac  <- as.factor(response_raw)
        grp_lvls <- levels(grp_fac)
        mask_list <- lapply(grp_lvls, function(l) valid_resp & grp_fac == l)
        names(mask_list) <- grp_lvls
        totals <- sapply(mask_list, sum)
        for (i in seq_along(grp_lvls)) {
          nm <- paste0("stat_g", i-1)
          tbl$getColumn(nm)$setTitle(grp_lvls[i])
          tbl$getColumn(nm)$setVisible(TRUE)
        }
        tbl$getColumn("pval")$setVisible(TRUE)
        get_counts <- function(mask) sapply(mask_list, function(m) sum(mask & m))
      }
      has_cont <- FALSE
      if (has_response) {
        if (is_cont_resp) {
          tbl$addRow(rowKey = paste0(response_var, "_mean"),
            values = list(variable = response_var, level = "Mean \u00B1 SD",
                          stat_overall = fmt_cont(response_raw)))
          mask <- !is.na(response_raw)
          tbl$addRow(rowKey = paste0(response_var, "_valid"),
            values = list(variable = "", level = "Valid", stat_overall = sum(mask)))
          n_miss <- sum(!mask)
          if (n_miss > 0)
            tbl$addRow(rowKey = paste0(response_var, "_missing"),
              values = list(variable = "", level = "Missing",
                            stat_overall = fmt_cat(n_miss, length(response_raw))))
        } else {
          mask <- valid_resp
          row_vals <- list(variable = response_var, level = "Valid", stat_overall = sum(mask))
          if (do_strat) {
            cnt <- get_counts(mask)
            for (i in seq_along(cnt)) row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], sum(mask))
            row_vals$pval <- ''
          }
          tbl$addRow(rowKey = response_var, values = row_vals)
          n_miss <- sum(!valid_resp)
          if (n_miss > 0) {
            miss_vals <- list(variable = "", level = "Missing",
                              stat_overall = fmt_cat(n_miss, length(response_raw)))
            if (do_strat) {
              for (i in seq_along(grp_lvls)) miss_vals[[paste0("stat_g", i-1)]] <- ''
              miss_vals$pval <- ''
            }
            tbl$addRow(rowKey = paste0(response_var, "_missing"), values = miss_vals)
          }
        }
      }
      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        for (v in names(cov_df)) {
          col    <- cov_df[[v]]
          is_cat <- is.factor(col) || is.character(col)
          if (is_cat && !is.factor(col)) col <- factor(col)
          n      <- length(col)
          n_miss <- sum(is.na(col))
          if (is_cat) {
            lvls <- levels(col)
            pval <- if (do_strat) tryCatch({
              ct <- table(col[valid_resp], grp_fac[valid_resp])
              suppressWarnings(chisq.test(ct)$p.value)
            }, error = function(e) '') else ''
            first <- TRUE
            for (lvl in lvls) {
              mask <- !is.na(col) & col == lvl
              row_vals <- list(variable = if (first) v else "", level = lvl,
                               stat_overall = fmt_cat(sum(mask), n))
              if (do_strat) {
                cnt <- get_counts(mask)
                for (i in seq_along(cnt)) row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], totals[i])
                row_vals$pval <- if (first) pval else ''
              }
              tbl$addRow(rowKey = paste0(v, "_", lvl), values = row_vals)
              first <- FALSE
            }
            if (n_miss > 0) {
              mask <- is.na(col)
              row_vals <- list(variable = "", level = "Missing", stat_overall = fmt_cat(n_miss, n))
              if (do_strat) {
                cnt <- get_counts(mask)
                for (i in seq_along(cnt)) row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], totals[i])
                row_vals$pval <- ''
              }
              tbl$addRow(rowKey = paste0(v, "_missing"), values = row_vals)
            }
          } else {
            has_cont <- TRUE
            row_vals <- list(variable = v, level = "Mean \u00B1 SD", stat_overall = fmt_cont(col))
            if (do_strat) {
              for (i in seq_along(mask_list)) row_vals[[paste0("stat_g", i-1)]] <- fmt_cont(col[mask_list[[i]]])
              row_vals$pval <- tryCatch({
                groups <- split(col[valid_resp], grp_fac[valid_resp])
                if (length(groups) == 2) t.test(groups[[1]], groups[[2]])$p.value
                else summary(aov(col ~ grp_fac))[[1]][["Pr(>F)"]][1]
              }, error = function(e) '')
            }
            tbl$addRow(rowKey = v, values = row_vals)
            if (n_miss > 0) {
              mask <- is.na(col)
              row_vals <- list(variable = "", level = "Missing", stat_overall = fmt_cat(n_miss, n))
              if (do_strat) {
                cnt <- get_counts(mask)
                for (i in seq_along(cnt)) row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], totals[i])
                row_vals$pval <- ''
              }
              tbl$addRow(rowKey = paste0(v, "_missing"), values = row_vals)
            }
          }
        }
      }
      notes <- c()
      if (n_removed_snp > 0) notes <- c(notes, paste0("Removed ", n_removed_snp, " cases with missing SNP values."))
      tbl$setNote(note = if (length(notes)) paste(notes, collapse=" ") else NULL, key = "cov_desc_note")
    },

    .fill_snp_summary = function(data, snp_vars, response_raw, response_type,
                                 subpop, cov_df = NULL) {
      tbl <- self$results$descGroup$snpSummaryTablesGroup$snpSummaryTable
      do_strat     <- isTRUE(subpop) && !is.null(response_raw) &&
                       (response_type == "binary" || response_type == "categorical")
      has_response <- !is.null(response_raw)
      has_cov      <- !is.null(cov_df) && ncol(cov_df) > 0
      n_total      <- nrow(data)
      grp_levels   <- if (do_strat) levels(response_raw) else NULL
      tbl$getColumn("group")$setVisible(do_strat)
      base_cc <- rep(TRUE, n_total)
      if (has_response) base_cc <- base_cc & !is.na(response_raw)
      if (has_cov) base_cc <- base_cc & complete.cases(cov_df)
      if (do_strat) {
        resp_base      <- response_raw[base_cc]
        stratum_totals <- table(factor(resp_base, levels = grp_levels))
      }
      row_key <- 0L
      for (snp_nm in snp_vars) {
        user_levels_sum <- get_snp_level_order(data[[snp_nm]])
        snp_raw         <- clean_null_alleles(as.character(data[[snp_nm]]))
        geno_obj        <- parse_genotype(snp_raw, user_levels_sum)
        if (is.null(geno_obj)) next
        cc_mask <- base_cc & !is.na(snp_raw)
        n_cc    <- sum(cc_mask)
        n_excluded <- n_total - n_cc
        snp_cc  <- snp_raw[cc_mask]
        geno_cc <- parse_genotype(snp_cc, user_levels_sum)
        if (is.null(geno_cc)) next
        resp_cc <- if (has_response) response_raw[cc_mask] else NULL
        sm_cc   <- summary(geno_cc)
        ref     <- get_ref_genotype(geno_cc, user_levels_sum)
        af_all  <- sm_cc$allele.freq
        allele_nms  <- rownames(af_all)
        ref_allele  <- strsplit(ref, "/")[[1]][1]
        alt_allele  <- setdiff(allele_nms, ref_allele)
        alt_allele  <- if (length(alt_allele)) alt_allele[1] else "?"
        alleles_label <- paste0(ref_allele, "/", alt_allele)
        compute_row <- function(mask = NULL) {
          snp_sub  <- if (is.null(mask)) snp_cc else snp_cc[mask]
          geno_sub <- if (is.null(mask)) geno_cc else parse_genotype(snp_sub, user_levels_sum)
          if (is.null(geno_sub)) return(NULL)
          sm   <- summary(geno_sub)
          af   <- sm$allele.freq
          props <- af[, "Proportion"]
          maf  <- if (alt_allele %in% rownames(af)) af[alt_allele, "Proportion"]
                  else if (length(props) >= 2) min(props, na.rm = TRUE) else NA_real_
          gf   <- sm$genotype.freq
          gf   <- tryCatch(reorder_geno(gf, ref, user_levels_sum), error = function(e) gf)
          gf   <- gf[rownames(gf) != "NA", , drop = FALSE]
          counts   <- as.integer(gf[, "Count"])
          len      <- length(counts)
          geno_str <- if (len == 3) paste(counts, collapse = " / ")
                      else if (len == 2) paste(c(counts, 0L), collapse = " / ")
                      else paste(counts, collapse = " / ")
          hwe <- tryCatch(genetics::HWE.exact(geno_sub)$p.value, error = function(e) NA_real_)
          list(n = sm$n.typed, maf = maf, genoCounts = geno_str, hwePval = hwe)
        }
        res_all <- compute_row()
        if (!is.null(res_all)) {
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            snp = snp_nm, alleles = alleles_label,
            group = if (do_strat) "All" else "",
            n = res_all$n, missing = if (n_excluded > 0L) n_excluded else '',
            maf = round(res_all$maf, 4), genoCounts = res_all$genoCounts, hwePval = res_all$hwePval))
        }
        if (do_strat) {
          resp_cc_chr <- as.character(resp_cc)
          for (lvl in grp_levels) {
            mask <- !is.na(resp_cc_chr) & resp_cc_chr == lvl
            res  <- compute_row(mask)
            if (is.null(res)) next
            n_excl <- max(0L, as.integer(stratum_totals[lvl] - res$n))
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              snp = "", alleles = "", group = lvl, n = res$n,
              missing = if (n_excl > 0L) n_excl else '',
              maf = round(res$maf, 4), genoCounts = res$genoCounts, hwePval = res$hwePval))
          }
        }
      }
      if (has_response || has_cov) {
        parts <- c(if (has_cov) "covariates", if (has_response) "response")
        tbl$setNote(
          note = paste0("Complete cases used: rows missing any ",
                        paste(parts, collapse = " or "), " or SNP value are excluded."),
          key = "missing_resp_cov")
      } else {
        tbl$setNote(note = NULL, key = "missing_resp_cov")
      }
    },

    .fill_allele_freq = function(tbl, sm, snp_nm, response_raw, subpop,
                                  response_type, snp_raw, show_missing = FALSE,
                                  n_miss_by_level = NULL, n_total_eligible = 0L,
                                  total_missing = 0L, user_levels = NULL) {
      af           <- sm$allele.freq
      allele_names <- rownames(af)
      if (!is.null(user_levels) && length(user_levels)) {
        ref <- NULL
        for (g in user_levels) {
          p <- strsplit(g, "/", fixed=TRUE)[[1]]
          if (length(p)==2 && p[1]==p[2]) { ref <- p[1]; break }
        }
        if (!is.null(ref) && ref %in% allele_names)
          allele_names <- c(ref, setdiff(allele_names, ref))
      }
      do_strat <- isTRUE(subpop) && !is.null(response_raw) &&
                   (response_type == "binary" || response_type == "categorical")
      if (do_strat) {
        grp_levels    <- levels(response_raw)
        resp_chr      <- as.character(response_raw)
        alleles_split <- strsplit(as.character(snp_raw), "/")
        for (i in seq_along(grp_levels))
          tbl$addColumn(name=paste0("stat_g", i-1), title=grp_levels[i], type="string")
      }
      for (al in allele_names) {
        if (!al %in% rownames(af)) next
        count    <- as.integer(af[al, "Count"])
        prop     <- round(af[al, "Proportion"] * 100, 1)
        row_vals <- list(allele = al, stat = fmt_catpct(count, prop))
        if (do_strat) {
          for (i in seq_along(grp_levels)) {
            lvl    <- grp_levels[i]
            idx    <- resp_chr == lvl
            all_al <- unlist(alleles_split[idx])
            n_al   <- sum(all_al == al, na.rm=TRUE)
            row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(n_al, length(all_al))
          }
        }
        tbl$addRow(rowKey=al, values=row_vals)
      }
      if (isTRUE(show_missing) && total_missing > 0L) {
        miss_vals <- list(allele = "Missing", stat = fmt_cat(total_missing, n_total_eligible))
        if (do_strat && !is.null(n_miss_by_level)) {
          for (j in seq_along(grp_levels)) {
            lvl        <- grp_levels[j]
            miss_count <- if (lvl %in% names(n_miss_by_level)) n_miss_by_level[lvl] else 0
            miss_vals[[paste0("stat_g", j-1)]] <- fmt_catn(miss_count)
          }
        }
        tbl$addRow(rowKey="missing", values=miss_vals)
      }
    },

    .fill_geno_freq = function(tbl, sm, ref, geno_obj, response,
                                response_type, subpop, response_raw,
                                show_missing = FALSE, n_miss_by_level = NULL,
                                n_total_eligible = 0L, total_missing = 0L,
                                user_levels = NULL) {
      # as.character() on a genetics::genotype object returns the same normalised
      # strings that appear in sm$genotype.freq rownames, so comparisons are exact.
      snp_chr <- as.character(geno_obj)
      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
        if (!is.numeric(response)) response <- as.numeric(as.character(response))
      }
      gf <- tryCatch(reorder_geno(sm$genotype.freq, ref, user_levels), error=function(e) sm$genotype.freq)
      gf <- gf[rownames(gf) != "NA", , drop = FALSE]
      do_strat <- isTRUE(subpop) && !is.null(response_raw) &&
                   (response_type == "binary" || response_type == "categorical")
      if (do_strat) {
        grp_levels <- levels(response_raw)
        resp_chr   <- as.character(response_raw)
        # Per-stratum totals = typed (non-NA SNP) observations in each group
        strat_totals <- sapply(grp_levels, function(lvl)
          sum(resp_chr == lvl & !is.na(snp_chr), na.rm = TRUE))
        for (i in seq_along(grp_levels))
          tbl$addColumn(name=paste0("stat_g", i-1), title=grp_levels[i], type="string")
      }
      for (i in seq_len(nrow(gf))) {
        geno <- rownames(gf)[i]
        if (geno == "NA") next
        count    <- as.integer(gf[i, "Count"])
        prop     <- gf[i, "Proportion"] * 100
        row_vals <- list(genotype = geno, stat = fmt_catpct(count, prop), responseStat = "")
        if (response_type == "quantitative" && !is.null(response)) {
          mask  <- snp_chr == geno & !is.na(snp_chr) & !is.na(response)
          n_mask <- sum(mask)
          if (n_mask > 0) {
            mn <- mean(response[mask], na.rm=TRUE)
            se <- sd(response[mask], na.rm=TRUE) / sqrt(n_mask)
            row_vals$responseStat <- sprintf("%.2f (%.2f)", mn, se)
          }
        }
        if (do_strat) {
          for (j in seq_along(grp_levels)) {
            lvl   <- grp_levels[j]
            idx   <- resp_chr == lvl & !is.na(resp_chr)
            n_g   <- sum(idx & snp_chr == geno, na.rm = TRUE)
            n_tot <- strat_totals[[lvl]]
            row_vals[[paste0("stat_g", j-1)]] <- fmt_cat(n_g, n_tot)
          }
        }
        tbl$addRow(rowKey=geno, values=row_vals)
      }
      if (isTRUE(show_missing) && total_missing > 0) {
        miss_vals <- list(genotype = "Missing", stat = fmt_cat(total_missing, n_total_eligible), responseStat = "")
        if (do_strat && !is.null(n_miss_by_level)) {
          for (j in seq_along(grp_levels)) {
            lvl        <- grp_levels[j]
            miss_count <- if (lvl %in% names(n_miss_by_level)) n_miss_by_level[lvl] else 0
            miss_vals[[paste0("stat_g", j-1)]] <- fmt_catn(miss_count)
          }
        }
        tbl$addRow(rowKey="missing", values=miss_vals)
      }
    },

    .fill_hwe = function(tbl, geno_obj, snp_nm, response_raw, subpop,
                          show_missing = FALSE, n_miss_by_level = NULL,
                          n_total_eligible = 0L, total_missing = 0L,
                          ref = NULL, user_levels = NULL) {
      tbl$getColumn("missing")$setVisible(isTRUE(show_missing))
      hw <- tryCatch(genetics::HWE.exact(geno_obj), error=function(e) NULL)
      if (is.null(hw)) return()

      # Derive display order from reorder_geno so labels and counts are consistent.
      # HWE.exact assigns N11/N12/N22 by its own internal frequency order, which
      # may differ from the user's factor order. We read counts from genotype.freq
      # directly to guarantee alignment with the column labels.
      get_ordered_counts <- function(go) {
        gf <- tryCatch(
          reorder_geno(summary(go)$genotype.freq, ref, user_levels),
          error = function(e) summary(go)$genotype.freq)
        gf <- gf[rownames(gf) != "NA", , drop = FALSE]
        list(labels  = rownames(gf),
             counts  = as.integer(gf[, "Count"]))
      }
      info <- get_ordered_counts(geno_obj)

      if (length(info$labels) == 3) {
        tbl$getColumn("n11")$setTitle(info$labels[1])
        tbl$getColumn("n12")$setTitle(info$labels[2])
        tbl$getColumn("n22")$setTitle(info$labels[3])
      }

      add_row <- function(key, label, counts, miss, p) {
        tbl$addRow(rowKey = key, values = list(
          group   = label,
          n11     = counts[1L], n12 = counts[2L], n22 = counts[3L],
          missing = miss, pval = p))
      }

      add_row("All", "All subjects", info$counts,
              if (isTRUE(show_missing)) total_missing else '', hw$p.value)

      if (isTRUE(subpop) && !is.null(response_raw)) {
        lvls <- levels(response_raw)
        if (length(lvls) <= 5) {
          for (lvl in lvls) {
            mask <- response_raw == lvl & !is.na(response_raw)
            if (sum(mask) == 0) next
            hw_sub <- tryCatch(genetics::HWE.exact(geno_obj[mask]), error=function(e) NULL)
            if (is.null(hw_sub)) next
            sub_info <- get_ordered_counts(geno_obj[mask])
            miss_count <- if (isTRUE(show_missing) && !is.null(n_miss_by_level) &&
                               lvl %in% names(n_miss_by_level)) n_miss_by_level[lvl] else ''
            add_row(lvl, lvl, sub_info$counts, miss_count, hw_sub$p.value)
          }
        }
      }
    },

    .plotMissingness = function(image, ...) {
      cache <- private$.miss_cache
      if (is.null(cache) || length(cache) == 0) return(FALSE)
      run_subpop <- isTRUE(self$options$subpop)
      threshold  <- self$options$missingnessThreshold
      if (!is.numeric(threshold) || is.na(threshold)) threshold <- 0.1
      all_nms <- names(cache)
      pct_all <- vapply(all_nms, function(nm) {
        d <- cache[[nm]]
        if (d$n_total_eligible == 0) return(0)
        d$total_missing / d$n_total_eligible * 100
      }, numeric(1))
      keep        <- pct_all > threshold
      n_hidden    <- sum(!keep)
      snp_nms     <- all_nms[keep]
      pct_overall <- pct_all[keep]
      n_snps      <- length(snp_nms)
      if (n_snps == 0) {
        opar <- par(no.readonly = TRUE); on.exit(par(opar))
        par(bg = "white", mar = c(1, 1, 3, 1))
        plot.new(); title(main = "SNP Missingness")
        text(0.5, 0.5, sprintf("No SNPs have missingness > %.1f%%\n(threshold: %.1f%%)", threshold, threshold),
             cex = 1.1, col = "#555555")
        return(TRUE)
      }
      grp_data <- NULL
      if (run_subpop) {
        all_grps <- unique(unlist(lapply(cache[snp_nms], function(d) names(d$n_miss_by_level))))
        if (length(all_grps) > 0) {
          grp_data <- lapply(all_grps, function(g) {
            vapply(snp_nms, function(nm) {
              d <- cache[[nm]]
              if (is.null(d$n_miss_by_level) || !g %in% names(d$n_miss_by_level)) return(NA_real_)
              if (d$n_total_eligible == 0) return(0)
              d$n_miss_by_level[[g]] / d$n_total_eligible * 100
            }, numeric(1))
          })
          names(grp_data) <- all_grps
        }
      }
      n_grps <- if (!is.null(grp_data)) length(grp_data) else 0
      opar <- par(no.readonly = TRUE); on.exit(par(opar))
      note_lines  <- if (n_hidden > 0) 1 else 0
      left_margin <- max(nchar(snp_nms)) * 0.55 + 1
      par(bg = "white", mar = c(4.5 + note_lines, left_margin, 3, 1.5))
      grp_pal <- c("#2980B9", "#C0392B", "#27AE60", "#8E44AD", "#E67E22")
      x_max <- max(pct_overall, unlist(grp_data), na.rm = TRUE)
      x_max <- max(x_max * 1.15, threshold * 1.5, 0.5)
      y_pos <- seq_len(n_snps)
      plot(NULL, xlim = c(0, x_max), ylim = c(0.5, n_snps + 0.5),
           xlab = "Missing genotypes (%)", ylab = "", main = "SNP Missingness",
           yaxt = "n", las = 1, bty = "l")
      axis(2, at = y_pos, labels = snp_nms, las = 1, tick = FALSE,
           cex.axis = min(1, 14 / n_snps + 0.3))
      abline(v = threshold, lty = 3, col = "#AAAAAA", lwd = 1.2)
      bar_col <- adjustcolor("#2C3E50", 0.20)
      bar_brd <- adjustcolor("#2C3E50", 0.45)
      rect(rep(0, n_snps), y_pos - 0.35, pct_overall, y_pos + 0.35,
           col = bar_col, border = bar_brd, lwd = 0.8)
      points(pct_overall, y_pos, pch = 19, col = "#2C3E50", cex = 1.1)
      if (!is.null(grp_data) && n_grps > 0) {
        offsets <- seq(-0.18, 0.18, length.out = n_grps)
        for (gi in seq_len(n_grps)) {
          col_i <- grp_pal[(gi - 1) %% length(grp_pal) + 1]
          points(grp_data[[gi]], y_pos + offsets[gi],
                 pch = 21, bg = adjustcolor(col_i, 0.70), col = col_i, cex = 0.90, lwd = 0.8)
        }
        legend("topright",
               legend = c("Overall", names(grp_data)),
               pch    = c(19, rep(21, n_grps)),
               col    = c("#2C3E50", grp_pal[seq_len(n_grps)]),
               pt.bg  = c(NA, adjustcolor(grp_pal[seq_len(n_grps)], 0.70)),
               pt.cex = c(1.1, rep(0.90, n_grps)), bty = "n", cex = 0.80)
      }
      text(pct_overall + x_max * 0.015, y_pos,
           labels = sprintf("%.1f%%", pct_overall), adj = 0, cex = 0.72, col = "#333333")
      note_txt <- if (n_hidden > 0)
        sprintf("%d SNP(s) with missingness \u2264 %.1f%% not shown  |  dashed line = threshold", n_hidden, threshold)
      else sprintf("Dashed line = %.1f%% threshold", threshold)
      mtext(note_txt, side = 1, line = 3.6, cex = 0.72, col = "#666666", adj = 0)
      TRUE
    },

    # ════════════════════════════════════════════════════════════════════════
    # LD / Haplotype private methods (verbatim from snpLDHaplo_b.R)
    # ════════════════════════════════════════════════════════════════════════
    .run_association = function(data, snp_vars, response_var, response_raw,
                                response_type, response, cov_df,
                                complete_mask, n_rows, opts) {
      run_snpAssoc       <- isTRUE(opts$snpAssoc)
      run_snpInteraction <- isTRUE(opts$snpInteraction)

      if (run_snpInteraction && (is.null(cov_df) || ncol(cov_df) == 0)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>SNP \u00D7 covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }

      if (run_snpInteraction && response_type == "categorical") {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>SNP \u00D7 covariate interaction is not available for categorical responses. Use binary or quantitative.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }

      arr <- self$results$assocGroup$assocSnpResults

      for (snp_nm in snp_vars) {
        snp_raw <- data[[snp_nm]]
        user_levels     <- get_snp_level_order(snp_raw)
        geno_obj        <- parse_genotype(snp_raw, user_levels)
        if (is.null(geno_obj)) next

        snp_raw_clean     <- clean_null_alleles(as.character(snp_raw))
        snp_complete_mask <- complete_mask & !is.na(snp_raw_clean)
        n_miss_assoc      <- n_rows - sum(snp_complete_mask)
        snp_raw_cc        <- as.factor(snp_raw_clean[snp_complete_mask])
        response_cc       <- response[snp_complete_mask]
        response_raw_cc   <- response_raw[snp_complete_mask]
        cov_df_cc         <- if (!is.null(cov_df)) cov_df[snp_complete_mask, , drop = FALSE] else NULL
        user_levels       <- get_snp_level_order(snp_raw)
        geno_obj_cc       <- parse_genotype(snp_raw_cc, user_levels)
        if (is.null(geno_obj_cc)) next

        item <- arr$get(key = snp_nm)
        ref  <- get_ref_genotype(geno_obj_cc, user_levels)

        item$typingRate$setContent(sprintf(
          "<b>Typed samples:</b> %d / %d (%.1f%%)",
          sum(snp_complete_mask), n_rows, sum(snp_complete_mask) / n_rows * 100))

        if (run_snpAssoc)
          private$.fill_assoc(item$assocTable, snp_raw_cc, ref, response_cc,
                              cov_df_cc, response_type, opts,
                              n_miss = n_miss_assoc, user_levels, response_raw_cc, snp_nm)

        if (run_snpInteraction && !is.null(cov_df_cc) && ncol(cov_df_cc) >= 1) {
          interaction_var <- names(cov_df_cc)[1]
          int_models      <- private$.get_interaction_models(opts)
          model_labels    <- c(codominant = "Codominant", dominant = "Dominant",
                               recessive  = "Recessive",  overdominant = "Overdominant",
                               logadditive = "Log-additive")
          int_lbl <- attr(self$data[[interaction_var]], "label") %||% interaction_var

          for (mdl in int_models) {
            mdl_label <- model_labels[[mdl]]
            mdl_item  <- item$interactionResults$get(key = mdl_label)

            if (isTRUE(opts$showInteractionTable))
              private$.fill_interaction(
                mdl_item$interactionTable, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw_cc, snp_nm)

            if (isTRUE(opts$showStratByCovariate)) {
              mdl_item$stratByCovariateHeading$setContent(
                paste0("<h3>Stratified by Covariate: ", int_lbl, "</h3>"))
              private$.fill_strat_by_covariate(
                mdl_item$stratByCovariate, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw_cc, snp_nm)
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
                paste0("<h3>Cross-Classification: ", snp_nm, " \u00D7 ", int_lbl, "</h3>"))
              private$.fill_cross_class(
                mdl_item$crossClassTable, snp_raw_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                response_type, opts, mdl, user_levels, response_raw_cc, snp_nm)
            }
          }
        }
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # .run_ldhaplo
    # ════════════════════════════════════════════════════════════════════════
    .get_interaction_models = function(opts) {
      c(
        if (isTRUE(opts$modelCodominant))   "codominant",
        if (isTRUE(opts$modelDominant))     "dominant",
        if (isTRUE(opts$modelRecessive))    "recessive",
        if (isTRUE(opts$modelOverdominant)) "overdominant",
        if (isTRUE(opts$modelLogAdditive))  "logadditive"
      )
    },

    .fill_assoc = function(tbl, snp_raw, ref, response, cov_df,
                           response_type, opts, n_miss = 0L, user_levels = NULL,
                           response_raw, snp_lbl) {

      is_categorical <- (response_type == "categorical")

      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR"
                                       else if (is_categorical) "OR" else "\u03B2")
      tbl$getColumn("genotype")$setTitle(snp_lbl)
      resp_lbl <- attr(self$data[[self$options$response]], "label") %||% self$options$response
      tbl$setTitle(paste0("Association with ", resp_lbl))

      if (response_type == "binary") {
        lv <- levels(as.factor(response_raw))
        tbl$getColumn("stat0")$setTitle(lv[1]); tbl$getColumn("stat1")$setTitle(lv[2])
        tbl$getColumn("stat0")$setVisible(TRUE); tbl$getColumn("stat1")$setVisible(TRUE)
      } else if (is_categorical) {
        # For categorical: stat0 shows N (%) per genotype within each category block.
        # The reference category is the first factor level of response_raw.
        ref_cat <- levels(as.factor(response_raw))[1]
        tbl$getColumn("stat0")$setTitle("N (%)")
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
        tbl$setNote(key = "multinom_ref",
                    note = paste0("Multinomial logistic regression. Reference category: \u2018",
                                  ref_cat, "\u2019. OR and CI vs. reference category."))
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }
      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))
      tbl$getColumn("BIC")$setVisible(isTRUE(opts$showAIC))
      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        cov_names <- sapply(names(cov_df), function(x) attr(self$data[[x]], "label") %||% x)
        note_txt  <- paste0("Model adjusted for: ", paste(cov_names, collapse = ", "))
        if (!is.na(n_miss) && n_miss > 0) note_txt <- paste0(note_txt, ".  ", n_miss, " observation(s) excluded.")
        tbl$setNote(note = note_txt, key = "covariates")
      } else if (!is.na(n_miss) && n_miss > 0) {
        tbl$setNote(note = paste0(n_miss, " observation(s) excluded."), key = "covariates")
      }
      models <- c(
        if (opts$modelCodominant)   "codominant",
        if (opts$modelDominant)     "dominant",
        if (opts$modelRecessive)    "recessive",
        if (opts$modelOverdominant) "overdominant",
        if (opts$modelLogAdditive)  "logadditive")
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      snp_char  <- as.character(snp_raw)
      # Respect the user's factor level order. When user_levels exist use them
      # directly (they are already ref-normalised by get_snp_level_order).
      # Only fall back to ref-anchored sort when there are no user levels.
      all_genos <- if (!is.null(user_levels) && length(user_levels) > 0) {
        user_levels
      } else {
        c(ref, setdiff(sort(unique(snp_char[!is.na(snp_char)])), ref))
      }
      n_fit <- sum(!is.na(snp_char) & !is.na(response) &
                     (if (!is.null(cov_df) && ncol(cov_df) > 0) complete.cases(cov_df) else TRUE))
      n_cov <- if (!is.null(cov_df)) ncol(cov_df) else 0L
      row_key     <- 0L
      any_clamped <- FALSE   # set TRUE if any OR was suppressed due to separation

      for (mdl in models) {
        snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
        res_list <- fit_model(snp_enc, response, cov_df, mdl, response_type, opts$ciWidth)
        if (is.null(res_list)) next

        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
        aic_val     <- { a <- res_list[[1]]$aic; if (!is.null(a) && !is.nan(a)) round(a,2) else NA_real_ }
        bic_val     <- private$.bic_from_aic(aic_val, mdl, n_fit, n_cov)

        # ── Categorical: one block per response category ──────────────────
        if (is_categorical) {
          cats <- unique(sapply(res_list, `[[`, "category"))
          st   <- private$.compute_stats(geno_labels, snp_char, response,
                                         response_type, response_raw)
          for (cat in cats) {
            cat_res  <- res_list[sapply(res_list, function(r) r$category == cat)]
            cat_sts  <- if (!is.null(st$by_cat)) st$by_cat[[cat]] else rep("", length(geno_labels))
            n_cat    <- sum(as.character(response_raw) == cat, na.rm = TRUE)
            if (mdl == "logadditive") {
              # Log-additive: single per-allele row per category (total N, not per-genotype)
              res <- cat_res[[1]]
              row_key <- row_key + 1L
              tbl$addRow(rowKey = as.character(row_key), values = list(
                model    = paste0(model_labels[mdl], " \u2014 ", cat, " (n=", n_cat, ")"),
                genotype = "Per allele",
                stat0    = sprintf("%d", n_cat),
                stat1    = "",
                effect   = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
                pval     = res$pval, AIC = aic_val, BIC = bic_val))
              next
            }
            # Header / reference row for this category
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              model    = paste0(model_labels[mdl], " \u2014 ", cat, " (n=", n_cat, ")"),
              genotype = geno_labels[1],
              stat0    = if (length(cat_sts) >= 1) cat_sts[1] else "",
              stat1    = "",
              effect   = 1., ciLow = '', ciHigh = '',
              pval     = cat_res[[1]]$global_p,
              AIC = aic_val, BIC = bic_val))
            # One row per non-ref genotype
            for (i in seq_along(cat_res)) {
              res <- cat_res[[i]]
              gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else res$comparison
              row_key <- row_key + 1L
              tbl$addRow(rowKey = as.character(row_key), values = list(
                model = "", genotype = gl,
                stat0 = if ((i + 1) <= length(cat_sts)) cat_sts[i + 1] else "",
                stat1 = "",
                effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
                pval = res$pval, AIC = '', BIC = ''))
            }
          }
          if (mdl == "codominant")
            tbl$setNote(key = "lrt",
                        note = "P-value in first row of each category is LRT for overall association")
          next
        }

        # ── Binary / quantitative ────────────────────────────────────────────
        st <- private$.compute_stats(geno_labels, snp_char, response, response_type,
                                     response_raw)
        if (mdl == "logadditive") {
          res <- res_list[[1]]; row_key <- row_key + 1L
          if (response_type == "binary") {
            lv <- levels(as.factor(response_raw))
            stat0_val <- sprintf("%d", sum(response_raw == lv[1], na.rm = TRUE))
            stat1_val <- sprintf("%d", sum(response_raw == lv[2], na.rm = TRUE))
          } else {
            stat0_val <- sprintf("%.2f (%.2f)", mean(response, na.rm=TRUE), sd(response, na.rm=TRUE))
            stat1_val <- " "
          }
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model = model_labels[mdl], genotype = "Per allele",
            stat0 = stat0_val, stat1 = stat1_val,
            effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval = res$pval, AIC = aic_val, BIC = bic_val))
          next
        }
        pval_row1 <- if (mdl == "codominant") res_list[[1]]$global_p else ''
        row_key <- row_key + 1L
        tbl$addRow(rowKey = as.character(row_key), values = list(
          model = model_labels[mdl], genotype = geno_labels[1],
          stat0 = st$s0[1], stat1 = st$s1[1],
          effect = if (response_type == "binary") 1. else 0.,
          ciLow = '', ciHigh = '', pval = pval_row1, AIC = aic_val, BIC = bic_val))
        if (mdl == "codominant")
          tbl$setNote(key = "lrt", note = "First p-value in Codominant is LRT for overall association")
        for (i in seq_along(res_list)) {
          res <- res_list[[i]]
          gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else res$comparison
          if (is.na(res$effect)) any_clamped <- TRUE
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model = "", genotype = gl,
            stat0 = if ((i+1) <= length(st$s0)) st$s0[i+1] else "-",
            stat1 = if ((i+1) <= length(st$s1)) st$s1[i+1] else "",
            effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
            pval = res$pval, AIC = '', BIC = ''))
        }
      }
      if (any_clamped)
        tbl$setNote(key = "separation",
                    note = "One or more OR/CI suppressed (shown as blank) due to complete or quasi-complete separation.")
    },

    # .fill_interaction, .fill_strat_by_covariate, .fill_strat_by_genotype,
    # .fill_cross_class are copied verbatim from snpAssoc_b.R — paste below.
    # (Omitted here for brevity in this skeleton; the actual merge pastes them in full.)
    .fill_interaction = function(tbl, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  int_models, user_levels = NULL, response_raw, snp_lbl) {
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
      adj_vars <- setdiff(names(cov_df), interaction_var)
      int_lbl  <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      int_type <- if (is.null(opts$interactionType)) "multiplicative" else opts$interactionType
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
      show_adj     <- isTRUE(opts$showInteractionAdjVars)
      any_clamped  <- FALSE
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      .label_term <- function(term, mdl, geno_labels) {
        lbl <- gsub("^snp([^:]+)", paste0(snp_lbl, "(\\1)"), term)
        lbl <- gsub(":snp([^:]+)", paste0(":", snp_lbl, "(\\1)"), lbl)
        lbl
      }
      row_key <- 0L
      last_pval_interaction <- NA_real_
      for (mdl in int_models) {
        snp_enc <- encode_model(as.character(snp_raw), ref, mdl, user_levels)
        if (int_type == "multiplicative") {
          res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                            interaction_var, mdl, response_type, opts$ciWidth,
                                            conditional = FALSE)
        } else if (int_type == "conditional_on_snp") {
          res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                            interaction_var, mdl, response_type, opts$ciWidth,
                                            conditional = TRUE, cond_var = "snp")
        } else {
          res_list <- fit_interaction_model(snp_enc, response, cov_df,
                                            interaction_var, mdl, response_type, opts$ciWidth,
                                            conditional = TRUE, cond_var = interaction_var)
        }
        if (is.null(res_list)) next
        p_inter <- attr(res_list, "pval_interaction")
        if (!is.null(p_inter) && !is.na(p_inter)) last_pval_interaction <- p_inter
        snp_char_l  <- as.character(snp_raw)
        all_genos_l <- if (!is.null(user_levels) && length(user_levels) > 0) {
          user_levels
        } else {
          c(ref, setdiff(sort(unique(snp_char_l[!is.na(snp_char_l)])), ref))
        }
        geno_labels_l <- private$.geno_labels_for_model(mdl, all_genos_l, ref)
        n_fit_bic <- sum(!is.na(snp_enc) & !is.na(response) & complete.cases(cov_df))
        n_cov_bic <- ncol(cov_df)
        first_row   <- TRUE
        for (res in res_list) {
          rtype <- if (is.null(res$row_type)) "snp" else res$row_type
          if (rtype == "adjustment" && !show_adj) next
          if (is.na(res$effect)) any_clamped <- TRUE
          row_key    <- row_key + 1L
          term_label <- .label_term(res$term, mdl, geno_labels_l)
          vals <- list(model = if (first_row) model_labels[mdl] else "",
                       term  = term_label,
                       effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high,
                       pval   = res$pval)
          if (isTRUE(opts$showAIC)) {
            aic_val <- if (first_row && !is.nan(res$aic)) round(res$aic, 2) else ""
            bic_val <- if (first_row && !is.nan(res$aic))
              private$.bic_from_aic(res$aic, mdl, n_fit_bic, n_cov_bic) else ""
            vals[["AIC"]] <- aic_val; vals[["BIC"]] <- bic_val
          }
          tbl$addRow(rowKey = as.character(row_key), values = vals)
          first_row <- FALSE
        }
      }
      if (!is.na(last_pval_interaction))
        tbl$setNote(
          note = paste0("Interaction p-value (LRT): ",
                        format.pval(last_pval_interaction, digits = 3, eps = 0.001)),
          key  = "interactionPval")
      if (any_clamped)
        tbl$setNote(key = "separation",
                    note = "One or more OR/CI suppressed (shown as blank) due to complete or quasi-complete separation.")
    },
    .fill_strat_by_covariate = function(arr, snp_raw, ref, response, cov_df,
                                        interaction_var, response_type, opts,
                                        int_models, user_levels = NULL, response_raw, snp_lbl) {
      int_var_data <- cov_df[[interaction_var]]
      if (length(table(int_var_data)) > 6) return()
      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      snp_char     <- as.character(snp_raw)
      all_genos <- if (!is.null(user_levels) && length(user_levels) > 0) user_levels else c(ref, setdiff(sort(unique(snp_char[!is.na(snp_char)])), ref))
      model_labels <- c(codominant="Codominant", dominant="Dominant", recessive="Recessive",
                        overdominant="Overdominant", logadditive="Log-additive")
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
          tbl$getColumn("stat0")$setVisible(TRUE);     tbl$getColumn("stat1")$setVisible(TRUE)
        } else {
          tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
          tbl$getColumn("stat1")$setVisible(FALSE)
        }
      }
      for (mdl in int_models) {
        snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
        res_list <- fit_interaction_model(snp_enc, response, cov_df, interaction_var, mdl,
                                          response_type, opts$ciWidth, conditional = TRUE)
        if (is.null(res_list)) next
        inter_only  <- res_list[sapply(res_list, function(r)
          is.null(r$row_type) || r$row_type == "interaction")]
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
        for (cl in cov_levels) {
          cl_label  <- as.character(cl)
          level_res <- inter_only[grepl(cl_label, sapply(inter_only, `[[`, "term"), fixed = TRUE)]
          if (length(level_res) == 0) next
          tbl <- arr$get(key = paste0(int_lbl, ": ", cl_label))
          mask_k <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
          if (!is.null(adj_cov_df) && ncol(adj_cov_df) > 0) mask_k <- mask_k & complete.cases(adj_cov_df)
          st <- private$.compute_stats(geno_labels, snp_char[mask_k], response[mask_k], response_type)
          row_key <- 0L
          if (mdl == "logadditive") {
            tbl$addRow(rowKey = "1", values = list(
              genotype = paste0(model_labels[mdl], " (per allele)"), stat0 = "", stat1 = "",
              effect = level_res[[1]]$effect, ciLow = level_res[[1]]$ci_low,
              ciHigh = level_res[[1]]$ci_high, pval = level_res[[1]]$pval))
            tbl$getColumn("stat1")$setVisible(FALSE); tbl$getColumn("stat0")$setVisible(FALSE)
            next
          }
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
              stat0 = if ((i+1) <= length(st$s0)) st$s0[i+1] else "-",
              stat1 = if ((i+1) <= length(st$s1)) st$s1[i+1] else " ",
              effect = res$effect, ciLow = res$ci_low, ciHigh = res$ci_high, pval = res$pval))
          }
        }
        note_txt  <- paste0("The reference group is <b>", snp_lbl, ": ", geno_labels[1], "</b> across all strata.")
        pval_interaction <- attr(res_list, "pval_interaction")
        note_pval <- paste0("Interaction p-value: ", format.pval(pval_interaction, digits = 3, eps = 0.001))
        tbl$setNote(note = note_txt,  key = "interStratCov")
        tbl$setNote(note = note_pval, key = "interStratCovPval")
      }
    },
    .fill_strat_by_genotype = function(arr, snp_raw, ref, response, cov_df,
                                       interaction_var, response_type, opts,
                                       int_models, user_levels = NULL, response_raw, snp_lbl) {
      snp_char     <- as.character(snp_raw)
      all_genos <- if (!is.null(user_levels) && length(user_levels) > 0) user_levels else c(ref, setdiff(sort(unique(snp_char[!is.na(snp_char)])), ref))
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
        geno_labels     <- private$.geno_labels_for_model(mdl, all_genos, ref)
        snp_enc_m       <- encode_model(snp_char, ref, mdl, user_levels)
        res_list        <- fit_interaction_model(snp_enc_m, response, cov_df, interaction_var, mdl,
                                                 response_type, opts$ciWidth,
                                                 conditional = TRUE, cond_var = "snp")
        if (is.null(res_list)) next
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
            tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
            tbl$getColumn("stat0")$setVisible(TRUE);     tbl$getColumn("stat1")$setVisible(TRUE)
          } else {
            tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
            tbl$getColumn("stat1")$setVisible(FALSE)
          }
          mask_g     <- snp_char %in% private$.split_genos(gl)
          int_g      <- int_var_data[mask_g]
          resp_g     <- response[mask_g]
          resp_raw_g <- response_raw[mask_g]
          if (response_type == "binary") {
            counts <- table(factor(int_g, levels = cov_levels), factor(resp_raw_g, levels = resp_lv))
            totals <- colSums(counts)
          }
          inter_only <- res_list[sapply(res_list, function(r) is.null(r$row_type) || r$row_type == "interaction")]
          gl_res <- inter_only[grepl(gl, sapply(inter_only, `[[`, "term"), fixed = TRUE)]
          if (length(gl_res) == 0) {
            has_ref_terms <- length(inter_only) >= length(geno_labels) * n_cov_contrasts
            gl_offset <- if (has_ref_terms) gl_idx - 1L else gl_idx - 2L
            start <- gl_offset * n_cov_contrasts + 1L
            end   <- min((gl_offset + 1L) * n_cov_contrasts, length(inter_only))
            gl_res <- if (start >= 1L && start <= length(inter_only)) inter_only[start:end] else list()
          }
          row_key <- 0L
          if (!is_numerical) {
            cl_ref  <- cov_levels[1]
            stat0   <- if (response_type == "binary") fmt_cat(counts[cl_ref, 1], totals[1]) else { vals <- resp_g[as.character(int_g) == cl_ref]; fmt_cont(vals) }
            stat1   <- if (response_type == "binary") fmt_cat(counts[cl_ref, 2], totals[2]) else ""
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              level = cl_ref, stat0 = stat0, stat1 = stat1,
              effect = if (response_type == "binary") 1.0 else 0.0, ciLow = "", ciHigh = "", pval = ""))
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
        note_txt  <- paste0("The reference group is <b>", interaction_var, ": ", cov_levels[1], "</b> across all strata.")
        pval_interaction <- attr(res_list, "pval_interaction")
        note_pval <- paste0("Interaction p-value: ", format.pval(pval_interaction, digits = 3, eps = 0.001))
        tbl$setNote(note = note_txt,  key = "interStatGeno")
        tbl$setNote(note = note_pval, key = "interStratGenoPval")
      }
    },

    .fill_cross_class = function(arr, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  int_models, user_levels = NULL, response_raw, snp_lbl) {
      if (int_models == "logadditive") return()
      int_var_data <- cov_df[[interaction_var]]
      if (length(table(int_var_data)) > 6) return()
      int_lbl    <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      snp_char   <- as.character(snp_raw)
      all_genos <- if (!is.null(user_levels) && length(user_levels) > 0) user_levels else c(ref, setdiff(sort(unique(snp_char[!is.na(snp_char)])), ref))
      cov_levels <- if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(int_var_data[!is.na(int_var_data)]))
      adj_vars   <- setdiff(names(cov_df), interaction_var)
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
        snp_enc          <- encode_model(snp_char, ref, mdl, user_levels)
        df_fit           <- data.frame(resp = response, snp = snp_enc, interaction_var = int_var_data)
        if (length(adj_vars) > 0) df_fit <- cbind(df_fit, cov_df[, adj_vars, drop=FALSE])
        adj_part         <- if (length(adj_vars) > 0) paste("+", safe_rhs(adj_vars)) else ""
        formula_str      <- paste("resp ~ snp * interaction_var", adj_part)
        formula_main_str <- paste("resp ~ snp + interaction_var", adj_part)
        fit <- if (response_type == "binary")
          glm(as.formula(formula_str), data = df_fit, family = binomial())
        else lm(as.formula(formula_str), data = df_fit)
        fit_main_cc <- if (response_type == "binary")
          glm(as.formula(formula_main_str), data = df_fit, family = binomial())
        else lm(as.formula(formula_main_str), data = df_fit)
        lrtest_str   <- if (response_type == "binary") "Chisq" else "F"
        lrtest_label <- if (response_type == "binary") "Pr(>Chi)" else "Pr(>F)"
        lrt_cc       <- tryCatch(anova(fit_main_cc, fit, test = lrtest_str), error = function(e) NULL)
        p_inter_cc   <- if (!is.null(lrt_cc)) lrt_cc[2, lrtest_label] else NA_real_
        betas       <- coef(fit)
        v_cov       <- vcov(fit)
        ci_z        <- qnorm(1 - (1 - opts$ciWidth/100)/2)
        geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
        for (j in seq_along(cov_levels)) {
          cl  <- cov_levels[j]
          tbl <- arr$get(key = paste0(int_lbl, ": ", as.character(cl)))
          mask_k <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
          st <- private$.compute_stats(geno_labels, snp_char[mask_k], response[mask_k], response_type)
          for (i in seq_along(geno_labels)) {
            gl <- geno_labels[i]
            if (i == 1 && j == 1) {
              tbl$addRow(rowKey = paste0(mdl, i), values = list(
                genotype = gl, stat0 = st$s0[i], stat1 = st$s1[i],
                effect = if (response_type == "binary") 1.0 else 0.0, ciLow = "", ciHigh = "", pval = ""))
              next
            }
            term_snp   <- paste0("snp", gl)
            term_cov   <- paste0("interaction_var", cl)
            term_inter <- paste0("snp", gl, ":interaction_var", cl)
            active_terms <- c(if (i > 1) term_snp, if (j > 1) term_cov, if (i > 1 && j > 1) term_inter)
            active_terms <- active_terms[active_terms %in% names(betas)]
            combined_beta <- sum(betas[active_terms])
            combined_se   <- sqrt(sum(v_cov[active_terms, active_terms]))
            z_val   <- combined_beta / combined_se
            p_val   <- 2 * (1 - pnorm(abs(z_val)))
            lo_beta <- combined_beta - ci_z * combined_se
            hi_beta <- combined_beta + ci_z * combined_se
            tbl$addRow(rowKey = paste0(mdl, i), values = list(
              genotype = gl, stat0 = st$s0[i], stat1 = st$s1[i],
              effect = if (response_type == "binary") .exp_or(combined_beta) else combined_beta,
              ciLow  = if (response_type == "binary") .exp_or(lo_beta) else lo_beta,
              ciHigh = if (response_type == "binary") .exp_or(hi_beta) else hi_beta,
              pval   = p_val))
          }
          tbl$getColumn("stat1")$setVisible(response_type == "binary")
        }
        note_txt  <- paste0("The reference group is <b>", interaction_var, ": ", cov_levels[1],
                            " and ", snp_lbl, ": ", geno_labels[1], "</b>")
        note_pval <- paste0("Interaction p-value: ", format.pval(p_inter_cc, digits = 3, eps = 0.001))
        tbl$setNote(note = note_txt,  key = "interCrossClass")
        tbl$setNote(note = note_pval, key = "interCrossClassPval")
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # Descriptive private methods (verbatim from snpDesc_b.R)
    # ════════════════════════════════════════════════════════════════════════
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

    .split_genos = function(gl)
      unlist(strsplit(gl, "(?<=[A-Za-z0-9*])-(?=[A-Za-z0-9*])", perl = TRUE)),

    .compute_stats = function(geno_labels, snp_char, response, response_type,
                              response_raw = NULL) {
      split_genos <- private$.split_genos

      # Normalise genotype strings so ref allele is always first (e.g. "G/A" →
      # "A/G" when A is the ref).  This must match the orientation produced by
      # parse_genotype(), otherwise %in% lookups silently miss heterozygotes
      # and every cell falls through to the "---" / zero-count guard.
      norm_snp_char <- function(sc) {
        # Determine ref allele from the first homozygote among ALL geno_labels.
        # geno_labels[1] may be a compound label (e.g. "A/A-B/B" for overdominant),
        # so scan all single-genotype labels for a homozygote X/X pattern.
        ref_al <- NULL
        for (lbl in geno_labels) {
          parts <- strsplit(lbl, "/", fixed = TRUE)[[1]]
          if (length(parts) == 2 && parts[1] == parts[2]) { ref_al <- parts[1]; break }
        }
        if (is.null(ref_al)) return(sc)
        sapply(sc, function(g) {
          if (is.na(g)) return(NA_character_)
          p <- strsplit(g, "/", fixed = TRUE)[[1]]
          if (length(p) == 2 && p[1] != p[2] && p[2] == ref_al)
            paste0(p[2], "/", p[1])
          else g
        }, USE.NAMES = FALSE)
      }
      sc <- norm_snp_char(snp_char)

      if (response_type == "binary") {
        # Use response_raw (original labels) for column titles / grouping so
        # that the counts match the header labels set in .fill_assoc.
        resp_grp <- if (!is.null(response_raw)) response_raw else response
        lv       <- levels(as.factor(resp_grp))
        if (length(lv) < 2) lv <- c(lv, "")
        n_col0 <- sum(resp_grp == lv[1] & !is.na(resp_grp))
        n_col1 <- sum(resp_grp == lv[2] & !is.na(resp_grp))
        stats0 <- character(length(geno_labels))
        stats1 <- character(length(geno_labels))
        for (i in seq_along(geno_labels)) {
          mask <- sc %in% split_genos(geno_labels[i]) & !is.na(resp_grp)
          n0   <- sum(mask & resp_grp == lv[1])
          n1   <- sum(mask & resp_grp == lv[2])
          stats0[i] <- sprintf("%d (%.1f%%)", n0, if (n_col0 > 0) n0 / n_col0 * 100 else 0)
          stats1[i] <- sprintf("%d (%.1f%%)", n1, if (n_col1 > 0) n1 / n_col1 * 100 else 0)
        }
        list(s0 = stats0, s1 = stats1)
      } else if (response_type == "categorical") {
        # For categorical: compute N (%) per genotype × category.
        # Returns a list of per-category vectors; callers index by category.
        resp_grp <- if (!is.null(response_raw)) response_raw else response
        cats     <- levels(as.factor(resp_grp))
        n_cats   <- sapply(cats, function(c) sum(resp_grp == c & !is.na(resp_grp)))
        cat_stats <- lapply(cats, function(cat) {
          n_cat <- n_cats[[cat]]
          sapply(seq_along(geno_labels), function(i) {
            mask <- sc %in% split_genos(geno_labels[i]) & !is.na(resp_grp)
            n    <- sum(mask & resp_grp == cat)
            sprintf("%d (%.1f%%)", n, if (n_cat > 0) n / n_cat * 100 else 0)
          })
        })
        names(cat_stats) <- cats
        # Also provide overall N per genotype for the header rows
        stats_total <- sapply(seq_along(geno_labels), function(i) {
          mask <- sc %in% split_genos(geno_labels[i]) & !is.na(resp_grp)
          sprintf("%d", sum(mask))
        })
        list(s0 = stats_total, s1 = rep("", length(geno_labels)),
             by_cat = cat_stats, cats = cats)
      } else {
        stats0 <- character(length(geno_labels))
        for (i in seq_along(geno_labels)) {
          vals <- as.numeric(response[sc %in% split_genos(geno_labels[i]) & !is.na(response)])
          stats0[i] <- if (length(vals) == 0) "---"
                       else sprintf("%.2f (%.2f)", mean(vals), sd(vals))
        }
        list(s0 = stats0, s1 = rep("", length(geno_labels)))
      }
    },

    .bic_from_aic = function(aic_val, mdl, n_fit, n_cov) {
      snp_df <- c(codominant = 2L, dominant = 1L, recessive = 1L,
                  overdominant = 1L, logadditive = 1L)
      if (is.null(aic_val) || is.na(aic_val) || is.nan(aic_val)) return(NA_real_)
      round(aic_val + (1L + n_cov + snp_df[[mdl]]) * (log(n_fit) - 2), 2)
    },

    # ── Association, interaction, strat, cross-class fill methods ───────────
    # These are copied verbatim from snpAssoc_b.R private list.
    # Self-references to self$results$snpResults become self$results$assocGroup$assocSnpResults
    # where they appear (only in .fill_cross_class's null check for arr items).
    .run_ldhaplo = function(data, snp_vars, response_var, response_raw,
                            response_type, response, cov_df,
                            complete_mask, n_rows, opts) {
      run_ldAnalysis       <- isTRUE(opts$ldAnalysis)
      run_ldMatrix         <- isTRUE(opts$ldMatrix)
      run_ldPlot           <- isTRUE(opts$ldPlot)
      run_haploFreq        <- isTRUE(opts$haploFreq)
      run_haploAssoc       <- isTRUE(opts$haploAssoc)
      run_haploInteraction <- isTRUE(opts$haploInteraction)
      run_subpop           <- isTRUE(opts$ldSubpop)

      if (run_haploInteraction && (is.null(cov_df) || ncol(cov_df) == 0)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>Haplotype \u00D7 covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploInteraction <- FALSE
      }
      if (run_haploAssoc && is.null(response_raw)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>Haplotype association requires a response variable.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploAssoc <- FALSE
      }
      if (response_type == "categorical" && (run_haploAssoc || run_haploInteraction)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>Haplotype association and interaction are not available for categorical responses. Use binary or quantitative.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploAssoc       <- FALSE
        run_haploInteraction <- FALSE
      }

      ld_res_grp    <- self$results$ldHaploGroup
      ld_res_grp$ldGroup$setVisible(run_ldAnalysis || run_ldMatrix || run_ldPlot)
      ld_res_grp$haploGroup$setVisible(run_haploFreq || run_haploAssoc || run_haploInteraction)

      # ── Parse genotypes ────────────────────────────────────────────────
      geno_list <- list()
      for (snp_nm in snp_vars) {
        snp_raw     <- data[[snp_nm]]
        user_levels <- get_snp_level_order(snp_raw)
        geno_obj    <- parse_genotype(snp_raw, user_levels)
        if (!is.null(geno_obj)) geno_list[[snp_nm]] <- geno_obj
      }
      if (length(geno_list) < 2) return()

      if (run_ldAnalysis || run_ldMatrix || run_ldPlot)
        private$.run_ld(geno_list, opts, run_ldAnalysis, run_ldMatrix, run_ldPlot)

      if (run_haploFreq || run_haploAssoc || run_haploInteraction)
        private$.run_haplo(geno_list, data, response, response_raw, response_type,
                           cov_df, opts, run_haploFreq, run_haploAssoc,
                           run_haploInteraction, run_subpop, complete_mask)
    },

    # ════════════════════════════════════════════════════════════════════════
    # Association private methods (verbatim from snpAssoc_b.R)
    # ════════════════════════════════════════════════════════════════════════
    .run_ld = function(geno_list, opts, run_ldAnalysis, run_ldMatrix, run_ldPlot) {
      nms   <- names(geno_list)
      n     <- length(nms)
      pairs <- combn(nms, 2, simplify = FALSE)
      ld_store <- list()
      for (pair in pairs) {
        key    <- paste(pair, collapse = "___")
        ld_res <- tryCatch(genetics::LD(geno_list[[pair[1]]], geno_list[[pair[2]]]),
                           error = function(e) NULL)
        if (!is.null(ld_res)) ld_store[[key]] <- ld_res
      }
      if (run_ldAnalysis) {
        tbl <- self$results$ldHaploGroup$ldGroup$ldTable
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          tbl$addRow(rowKey = paste(pair, collapse="_"), values = list(
            snp1   = pair[1], snp2 = pair[2],
            r2     = round(ld_res$`r`^2,  3),
            Dprime = round(ld_res$`D'`,   3),
            D      = round(ld_res$`D`,    3),
            pval   = ld_res$`P-value`))
        }
      }
      if (run_ldMatrix) {
        mtbl   <- self$results$ldHaploGroup$ldGroup$ldMatrixTable
        metric <- opts$ldMetric
        for (nm in nms) {
          safe_nm <- gsub("[^A-Za-z0-9_]","_",nm)
          mtbl$addColumn(name = safe_nm, title = nm, type = "text")
        }
        upper_mat <- matrix("", n, n, dimnames = list(nms, nms))
        lower_mat <- matrix("", n, n, dimnames = list(nms, nms))
        diag(upper_mat) <- nms; diag(lower_mat) <- nms
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          p_val  <- ld_res$`P-value`
          p_str  <- if (!is.na(p_val)) { if (p_val < 0.001) "< .001" else sprintf("%.3f",p_val) } else ""
          up_val <- switch(metric,
            Dprime = sprintf("%.3f", round(ld_res$`D'`,  3)),
            r2     = sprintf("%.3f", round(ld_res$`r`^2, 3)),
            D      = sprintf("%.3f", round(ld_res$`D`,   3)))
          upper_mat[pair[1], pair[2]] <- up_val
          lower_mat[pair[2], pair[1]] <- p_str
        }
        for (i in seq_len(n)) {
          row_vals <- list(snp = nms[i])
          for (j in seq_len(n)) {
            safe_nm <- gsub("[^A-Za-z0-9_]","_",nms[j])
            row_vals[[safe_nm]] <- if(i==j) nms[i] else if(j>i) upper_mat[i,j] else lower_mat[i,j]
          }
          mtbl$addRow(rowKey = paste0("row_",i), values = row_vals)
        }
        metric_label <- switch(metric, Dprime="D'", r2="r²", D="D")
        mtbl$setNote(key="layout",
                     note=paste0("Upper triangle: ", metric_label,
                                 ". Lower triangle: P-value. Diagonal: SNP name."))
      }
      if (run_ldPlot) {
        private$.ld_store  <- ld_store
        private$.ld_nms    <- nms
        private$.ld_metric <- opts$ldMetric
        self$results$ldHaploGroup$ldGroup$ldPlotImage$setState(
          list(ld_store = ld_store, nms = nms, metric = opts$ldMetric))
      }
    },

    .render_ld_plot = function(image, ggtheme, theme, ...) {
      state <- image$state
      if (is.null(state)) return(FALSE)
      ld_store <- state$ld_store; nms <- state$nms; metric <- state$metric; n <- length(nms)
      metric_label <- switch(metric, Dprime="D'", r2="r²", D="D")
      df_rows <- list()
      for (i in seq_len(n)) for (j in seq_len(n)) {
        val <- if (i==j) 1.0 else {
          key <- paste(c(nms[min(i,j)], nms[max(i,j)]), collapse="___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) switch(metric,
            Dprime = abs(as.numeric(ld_res$`D'`)),
            r2     = as.numeric(ld_res$`r`)^2,
            D      = abs(as.numeric(ld_res$`D`))) else NA_real_
        }
        df_rows[[length(df_rows)+1L]] <- data.frame(
          SNP1  = factor(nms[i], levels=rev(nms)),
          SNP2  = factor(nms[j], levels=nms),
          value = val, stringsAsFactors=FALSE)
      }
      df <- do.call(rbind, df_rows)
      p_mat <- matrix(NA_real_, n, n, dimnames=list(nms,nms))
      for (pk in names(ld_store)) {
        parts <- strsplit(pk,"___")[[1]]
        pv    <- ld_store[[pk]]$`P-value`
        p_mat[parts[1],parts[2]] <- pv; p_mat[parts[2],parts[1]] <- pv
      }
      df$label <- ""
      for (k in seq_len(nrow(df))) {
        i_nm <- as.character(df$SNP1[k]); j_nm <- as.character(df$SNP2[k])
        i_idx <- which(nms==i_nm); j_idx <- which(nms==j_nm)
        if (i_idx > j_idx) {
          pv <- p_mat[i_nm, j_nm]
          df$label[k] <- if (!is.na(pv)) { if(pv<0.001) "<.001" else sprintf("%.3f",pv) } else ""
        } else if (i_idx < j_idx) {
          key <- paste(c(nms[min(i_idx,j_idx)], nms[max(i_idx,j_idx)]), collapse="___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) {
            raw <- switch(metric, r2=ld_res$`r`^2, Dprime=ld_res$`D'`, D=ld_res$`D`)
            df$label[k] <- sprintf("%.3f", round(as.numeric(raw),3))
          }
        } else { df$label[k] <- i_nm }
      }
      colour_label <- switch(metric, Dprime="|D'|", r2="r²", D="|D|")
      p <- ggplot2::ggplot(df, ggplot2::aes(x=SNP2, y=SNP1, fill=value)) +
        ggplot2::geom_tile(colour="white", linewidth=0.5) +
        ggplot2::geom_text(ggplot2::aes(label=label), size=3, colour="grey10") +
        ggplot2::scale_fill_gradientn(
          colours  = c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b"),
          limits   = c(0,1), na.value="grey85", name=colour_label) +
        ggplot2::scale_x_discrete(position="bottom") +
        ggplot2::labs(title=paste0("LD Heatmap  \u2022  upper: ",metric_label," | lower: p-value"),
                      x=NULL, y=NULL) +
        ggplot2::theme_minimal(base_size=11) +
        ggplot2::theme(
          axis.text.x=ggplot2::element_text(angle=45,hjust=1,vjust=1),
          axis.text.y=ggplot2::element_text(hjust=1),
          panel.grid=ggplot2::element_blank(),
          legend.position="right",
          plot.title=ggplot2::element_text(size=11,face="bold",
                                           margin=ggplot2::margin(b=8)))
      print(p); TRUE
    },

    .run_haplo = function(geno_list, data, response, response_raw, response_type,
                           cov_df, opts, run_haploFreq, run_haploAssoc,
                           run_haploInteraction, run_subpop, complete_mask) {
      snp_names   <- names(geno_list)
      allele_mat  <- do.call(cbind, lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]])))
      geno_setup  <- tryCatch(haplo.stats::setupGeno(allele_mat, locus.label = snp_names),
                              error = function(e) NULL)
      if (is.null(geno_setup)) return()
      u_alleles      <- attr(geno_setup, "unique.alleles")
      snp_miss_mask  <- apply(is.na(allele_mat), 1, all)
      keep           <- complete_mask & !snp_miss_mask
      n_miss         <- sum(snp_miss_mask & complete_mask)
      if (run_haploFreq)
        private$.compute_haplo_freqs(geno_setup, response_raw, response_type, keep,
                                     n_miss, opts, run_subpop, snp_names, u_alleles)
      if (run_haploAssoc && !is.null(response))
        private$.compute_haplo_assoc(geno_setup, response, response_type, cov_df, keep,
                                     n_miss, opts, snp_names, u_alleles)
      if (run_haploInteraction && !is.null(cov_df) && !is.null(response))
        private$.compute_haplo_interaction(geno_setup, response, response_type, cov_df, keep,
                                           n_miss, opts, snp_names, u_alleles)
    },

    .compute_haplo_freqs = function(geno_setup, response_raw, response_type, keep,
                                    n_miss, opts, run_subpop, snp_names, u_alleles) {
      tbl <- self$results$ldHaploGroup$haploGroup$haploFreqTable
      tbl$setTitle("<b>Haplotype Frequencies</b>")
      do_strat_haplo   <- isTRUE(run_subpop) && !is.null(response_raw) &&
                          identical(response_type, "binary")
      grp_levels_haplo <- levels(response_raw[keep])
      if (do_strat_haplo) {
        tbl$addColumn(name="freq_g0", title=as.character(grp_levels_haplo[1]), type="number", format="zto")
        tbl$addColumn(name="freq_g1", title=as.character(grp_levels_haplo[2]), type="number", format="zto")
      }
      em_all <- tryCatch(haplo.stats::haplo.em(subset_geno(geno_setup, keep), locus.label=snp_names),
                         error=function(e) NULL)
      if (!is.null(em_all)) {
        freqs    <- em_all$hap.prob
        rare_sum <- 0
        em_grp   <- list()
        grp_freq <- list()
        if (do_strat_haplo) {
          for (lvl in grp_levels_haplo) {
            keep_lvl <- keep & as.character(response_raw)==lvl
            if (sum(keep_lvl) < 5) next
            em_grp[[lvl]] <- tryCatch(
              haplo.stats::haplo.em(subset_geno(geno_setup, keep_lvl), locus.label=snp_names),
              error=function(e) NULL)
          }
          grp_freq <- lapply(em_grp, function(em_g) {
            if (is.null(em_g)) return(list())
            setNames(as.list(round(em_g$hap.prob,3)),
                     sapply(seq_len(nrow(em_g$haplotype)), function(j)
                       decode_haplo_row(as.numeric(em_g$haplotype[j,]), u_alleles)))
          })
          grp_levels <- levels(as.factor(response_raw[keep]))
          tbl$getColumn('freq_g0')$setVisible(TRUE)
          tbl$getColumn('freq_g0')$setTitle(as.character(grp_levels[1]))
          tbl$getColumn('freq_g1')$setVisible(TRUE)
          tbl$getColumn('freq_g1')$setTitle(as.character(grp_levels[2]))
        } else {
          tbl$getColumn('freq_g0')$setVisible(FALSE)
          tbl$getColumn('freq_g1')$setVisible(FALSE)
        }
        sorted_idx <- order(freqs, decreasing = TRUE)
        for (i in sorted_idx) {
          if (freqs[i] < opts$haploFreqMin) { rare_sum <- rare_sum + freqs[i]; next }
          label    <- decode_haplo_row(as.numeric(em_all$haplotype[i,]), u_alleles)
          row_vals <- list(haplotype=label, freq=round(freqs[i],3))
          if (do_strat_haplo) {
            row_vals$freq_g0 <- grp_freq[[grp_levels_haplo[1]]][[label]] %||% NA_real_
            row_vals$freq_g1 <- grp_freq[[grp_levels_haplo[2]]][[label]] %||% NA_real_
          }
          tbl$addRow(rowKey=paste0("f",i), values=row_vals)
        }
        if (rare_sum > 0) {
          row_vals <- list(haplotype=paste0("Rare (<",opts$haploFreqMin,")"), freq=round(rare_sum,3))
          if (do_strat_haplo) {
            em0 <- em_grp[[grp_levels_haplo[1]]]; em1 <- em_grp[[grp_levels_haplo[2]]]
            rare_g0 <- if (!is.null(em0)) round(sum(em0$hap.prob[em0$hap.prob < opts$haploFreqMin]),3) else NA_real_
            rare_g1 <- if (!is.null(em1)) round(sum(em1$hap.prob[em1$hap.prob < opts$haploFreqMin]),3) else NA_real_
            row_vals$freq_g0 <- if(!is.na(rare_g0) && rare_g0>0) rare_g0 else NA_real_
            row_vals$freq_g1 <- if(!is.na(rare_g1) && rare_g1>0) rare_g1 else NA_real_
          }
          tbl$addRow(rowKey="rare_freq", values=row_vals)
        }
      }
      if (n_miss > 0) tbl$setNote(note=paste0(n_miss," observation(s) with missing data excluded."), key="missing_snp")
      else tbl$setNote(note=NULL, key="missing_snp")
    },

    .compute_haplo_assoc = function(geno_setup, response, response_type, cov_df, keep,
                                    n_miss, opts, snp_names, u_alleles) {
      family <- if (response_type == "binary") "binomial" else "gaussian"
      y_sub  <- if (response_type == "binary") as.numeric(as.factor(response[keep])) - 1L
                else response[keep]
      m_model      <- data.frame(y = y_sub)
      m_model$geno <- subset_geno(geno_setup, keep)
      if (!is.null(cov_df)) {
        m_model    <- cbind(m_model, cov_df[keep, , drop = FALSE])
        formula_str <- paste("y ~ geno +", safe_rhs(names(cov_df)))
      } else {
        formula_str <- "y ~ geno"
      }
      haplo_fit <- tryCatch(
        haplo.stats::haplo.glm(as.formula(formula_str), family = family, data = m_model,
                               na.action = na.geno.keep,
                               control = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)),
        error = function(e) {
          self$results$validationMsg$setContent(paste0("<b>Haplotype GLM error:</b> ", e$message)); NULL
        })
      if (!is.null(haplo_fit)) {
        tbl <- self$results$ldHaploGroup$haploGroup$haploAssocTable
        tbl$setTitle("<b>Haplotype Association</b>")
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
        null_formula_str <- if (!is.null(cov_df) && ncol(cov_df) > 0)
          paste("y ~", safe_rhs(names(cov_df))) else "y ~ 1"
        haplo_null_fit <- tryCatch(
          if (family == "binomial") glm(as.formula(null_formula_str), family = binomial(), data = m_model)
          else lm(as.formula(null_formula_str), data = m_model),
          error = function(e) NULL)
        p_lrt_assoc <- NA_real_
        if (!is.null(haplo_null_fit)) {
          dev_diff <- deviance(haplo_null_fit) - haplo_fit$deviance
          df_diff  <- (haplo_fit$df.null - haplo_fit$df.residual) - (haplo_null_fit$df.null - haplo_null_fit$df.residual)
          p_lrt_assoc <- if (!is.na(dev_diff) && !is.na(df_diff) && df_diff > 0)
            pchisq(dev_diff, df = df_diff, lower.tail = FALSE) else NA_real_
        }
        tbl$setNote(note = paste0("Likelihood ratio test for overall haplotype association: P = ",
                                  format.pval(p_lrt_assoc, digits = 3)), key = "lrt_assoc")
        label_from_unique_row <- function(row_vec) paste(as.character(row_vec), collapse = "-")
        coef_sum  <- tryCatch(summary(haplo_fit)$coefficients, error = function(e) NULL)
        ci_mat    <- tryCatch(confint(haplo_fit, level = opts$ciWidth / 100), error = function(e) NULL)
        haplo_rows <- if (!is.null(coef_sum)) grep("^geno", rownames(coef_sum)) else integer(0)
        se_col    <- if (!is.null(coef_sum) && "SE" %in% colnames(coef_sum)) "SE" else "se"
        get_stats <- function(pos) {
          row_idx <- if (!is.na(pos) && pos >= 1L && pos <= length(haplo_rows)) haplo_rows[pos] else NA_integer_
          if (is.na(row_idx) || is.null(coef_sum) || row_idx < 1L || row_idx > nrow(coef_sum))
            return(list(beta = NA_real_, se = NA_real_, pval = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
          rn   <- rownames(coef_sum)[row_idx]
          beta <- coef_sum[row_idx, "coef"]
          se   <- coef_sum[row_idx, se_col]
          pval <- coef_sum[row_idx, "pval"]
          if (!is.null(ci_mat) && rn %in% rownames(ci_mat)) {
            ci_lo <- ci_mat[rn, 1]; ci_hi <- ci_mat[rn, 2]
          } else {
            z <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
            ci_lo <- beta - z * se; ci_hi <- beta + z * se
          }
          list(beta = beta, se = se, pval = pval, ci_lo = ci_lo, ci_hi = ci_hi)
        }
        make_row <- function(label, freq, stats) {
          b <- stats$beta; lo <- stats$ci_lo; hi <- stats$ci_hi
          list(haplotype = label, freq = round(freq, 4),
               effect = if (response_type == "binary") .exp_or(b)  else b,
               ciLow  = if (response_type == "binary") .exp_or(lo) else lo,
               ciHigh = if (response_type == "binary") .exp_or(hi) else hi,
               pval   = stats$pval)
        }
        base_idx   <- haplo_fit$haplo.base
        base_label <- label_from_unique_row(haplo_fit$haplo.unique[base_idx, ])
        base_freq  <- haplo_fit$haplo.freq[base_idx]
        tbl$addRow(rowKey = "base", values = list(
          haplotype = paste0(base_label, " (Ref)"), freq = round(base_freq, 4),
          effect = if (response_type == "binary") 1.0 else 0.0, ciLow = '', ciHigh = '', pval = ''))
        common_idx   <- haplo_fit$haplo.common
        common_freqs <- haplo_fit$haplo.freq[common_idx]
        sorted_j     <- order(common_freqs, decreasing = TRUE)
        for (j in sorted_j) {
          h_idx   <- common_idx[j]
          h_label <- label_from_unique_row(haplo_fit$haplo.unique[h_idx, ])
          h_freq  <- haplo_fit$haplo.freq[h_idx]
          tbl$addRow(rowKey = paste0("h", j), values = make_row(h_label, h_freq, get_stats(j)))
        }
        has_rare <- isTRUE(haplo_fit$haplo.rare.term) || (length(haplo_fit$haplo.rare) > 0)
        if (has_rare) {
          rare_freq <- sum(haplo_fit$haplo.freq[haplo_fit$haplo.rare])
          tbl$addRow(rowKey = "rare",
            values = make_row(paste0("Rare (<", opts$haploFreqMin, ")"), rare_freq,
                              get_stats(length(common_idx) + 1L)))
        }
        if (!is.null(cov_df) && ncol(cov_df) > 0) {
          cov_names <- sapply(names(cov_df), function(x) attr(self$data[[x]], "label") %||% x)
          tbl$setNote(note = paste0("Model adjusted for: ", paste(cov_names, collapse = ", ")), key = "covariates")
        } else tbl$setNote(note = NULL, key = "covariates")
        if (n_miss > 0) tbl$setNote(note = paste0(n_miss, " observation(s) with missing data excluded."), key = "missing_snp")
        else tbl$setNote(note = NULL, key = "missing_snp")
      }
    },

    .compute_haplo_interaction = function(geno_setup, response, response_type, cov_df, keep,
                                          n_miss, opts, snp_names, u_alleles) {
      int_var      <- names(cov_df)[1]
      int_var_vals <- cov_df[[int_var]]
      if (!is.factor(int_var_vals) && !is.character(int_var_vals)) {
        self$results$validationMsg$setContent(
          paste0("<p style='color:orange;'>Haplotype interaction tables require a categorical covariate. '",
                 int_var, "' is numeric \u2014 please convert it to a factor.</p>"))
        self$results$validationMsg$setVisible(TRUE)
        return()
      }
      if (length(unique(na.omit(as.character(int_var_vals)))) > 6) {
        self$results$validationMsg$setContent(
          paste0("<p style='color:orange;'>Haplotype interaction tables require a covariate with at most 6 categories. '",
                 int_var, "' has more.</p>"))
        self$results$validationMsg$setVisible(TRUE)
        return()
      }
      adj_vars   <- setdiff(names(cov_df), int_var)
      family_int <- if (response_type == "binary") "binomial" else "gaussian"
      is_binary  <- (response_type == "binary")
      y_int  <- if (is_binary) as.numeric(as.factor(response[keep])) - 1L else response[keep]
      m_int  <- data.frame(y = y_int)
      m_int$geno <- subset_geno(geno_setup, keep)
      if (!is.null(cov_df)) m_int <- cbind(m_int, cov_df[keep, , drop = FALSE])
      adj_part         <- if (length(adj_vars) > 0) paste("+", safe_rhs(adj_vars)) else ""
      formula_mult_str <- paste("y ~ geno *", safe_term(int_var), adj_part)
      formula_add_str  <- paste("y ~ geno +", safe_term(int_var), adj_part)
      haplo_fit_mult <- tryCatch(
        haplo.stats::haplo.glm(as.formula(formula_mult_str), family = family_int, data = m_int,
                               na.action = na.geno.keep,
                               control = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)),
        error = function(e) {
          self$results$validationMsg$setContent(paste0("<b>Haplotype interaction GLM error:</b> ", e$message)); NULL
        })
      haplo_fit_add <- tryCatch(
        haplo.stats::haplo.glm(as.formula(formula_add_str), family = family_int, data = m_int,
                               na.action = na.geno.keep,
                               control = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)),
        error = function(e) NULL)
      if (is.null(haplo_fit_mult)) return()
      coef_sum <- tryCatch(summary(haplo_fit_mult)$coefficients, error = function(e) NULL)
      if (is.null(coef_sum)) return()
      se_col   <- if ("SE" %in% colnames(coef_sum)) "SE" else "se"
      vcov_mat <- tryCatch(vcov(haplo_fit_mult), error = function(e) NULL)
      decode_haplo_label <- function(row_vec) paste(as.character(row_vec), collapse = "-")
      rare_label  <- paste0("Rare (<", opts$haploFreqMin, ")")
      z_crit      <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
      base_idx    <- haplo_fit_mult$haplo.base
      base_label  <- decode_haplo_label(haplo_fit_mult$haplo.unique[base_idx, ])
      int_var_factor <- as.factor(m_int[[int_var]])
      covar_levels   <- levels(int_var_factor)
      ref_covar_lvl  <- covar_levels[1]
      common_idx   <- haplo_fit_mult$haplo.common
      common_freqs <- haplo_fit_mult$haplo.freq[common_idx]
      has_rare     <- isTRUE(haplo_fit_mult$haplo.rare.term) || (length(haplo_fit_mult$haplo.rare) > 0)
      rare_freq    <- if (has_rare) sum(haplo_fit_mult$haplo.freq[haplo_fit_mult$haplo.rare]) else 0
      all_haplo_entries <- list()
      all_haplo_entries[["base"]] <- list(label = base_label, freq = haplo_fit_mult$haplo.freq[base_idx], coef_pos = NA_integer_)
      sorted_common_j <- order(common_freqs, decreasing = TRUE)
      for (j in sorted_common_j) {
        h_idx   <- common_idx[j]
        h_label <- decode_haplo_label(haplo_fit_mult$haplo.unique[h_idx, ])
        all_haplo_entries[[h_label]] <- list(label = h_label, freq = haplo_fit_mult$haplo.freq[h_idx], coef_pos = j)
      }
      if (has_rare) all_haplo_entries[["rare"]] <- list(label = rare_label, freq = rare_freq, coef_pos = length(common_idx) + 1L)
      all_coef_names <- rownames(coef_sum)
      haplo_rows     <- grep("^geno", all_coef_names)
      get_beta_se    <- function(nm) {
        idx <- match(nm, all_coef_names)
        if (is.na(idx)) return(list(beta = NA_real_, se = NA_real_))
        list(beta = coef_sum[idx, "coef"], se = coef_sum[idx, se_col])
      }
      combine_coefs <- function(nm1, nm2 = NULL, sign2 = 1) {
        i1 <- match(nm1, all_coef_names)
        if (is.na(i1)) return(list(beta = NA_real_, se = NA_real_))
        b1 <- coef_sum[i1, "coef"]
        if (is.null(nm2)) return(list(beta = b1, se = coef_sum[i1, se_col]))
        i2 <- match(nm2, all_coef_names)
        if (is.na(i2)) return(list(beta = b1, se = coef_sum[i1, se_col]))
        b2 <- coef_sum[i2, "coef"]
        beta_comb <- b1 + sign2 * b2
        if (!is.null(vcov_mat) && i1 <= nrow(vcov_mat) && i2 <= nrow(vcov_mat)) {
          var_comb <- vcov_mat[i1,i1] + vcov_mat[i2,i2] + 2*sign2*vcov_mat[i1,i2]
          se_comb  <- sqrt(max(0, var_comb))
        } else se_comb <- sqrt(coef_sum[i1, se_col]^2 + coef_sum[i2, se_col]^2)
        list(beta = beta_comb, se = se_comb)
      }
      fmt_or_ci <- function(beta, se) {
        if (is.na(beta) || is.na(se)) return(NA_character_)
        sprintf("%.2f (%.2f\u2013%.2f)", .exp_or(beta), .exp_or(beta - z_crit*se), .exp_or(beta + z_crit*se))
      }
      fmt_b_ci <- function(beta, se) {
        if (is.na(beta) || is.na(se)) return(NA_character_)
        sprintf("%.3f (%.3f\u2013%.3f)", beta, beta - z_crit*se, beta + z_crit*se)
      }
      fmt_effect_ci <- if (is_binary) fmt_or_ci else fmt_b_ci
      build_notes <- function(tbl) {
        if (length(adj_vars) > 0) {
          adj_lbl <- sapply(adj_vars, function(x) attr(self$data[[x]], "label") %||% x)
          tbl$setNote(note = paste0("Adjusted for: ", paste(adj_lbl, collapse = ", ")), key = "intcov")
        }
        if (n_miss > 0) tbl$setNote(note = paste0(n_miss, " observation(s) with missing data excluded."), key = "missing_snp")
        else tbl$setNote(note = NULL, key = "missing_snp")
      }
      p_inter <- NA_real_
      if (!is.null(haplo_fit_add)) {
        dev_diff <- haplo_fit_add$deviance - haplo_fit_mult$deviance
        df_diff  <- haplo_fit_add$df.residual - haplo_fit_mult$df.residual
        p_inter  <- if (!is.na(dev_diff) && !is.na(df_diff) && df_diff > 0)
          pchisq(dev_diff, df = df_diff, lower.tail = FALSE) else NA_real_
      }
      covar_main_nms <- grep(paste0("^", int_var, ".+$"), all_coef_names, value = TRUE)
      covar_main_nms <- covar_main_nms[!grepl(":", covar_main_nms)]
      covar_lvl_to_main <- setNames(covar_main_nms, sub(paste0("^", int_var), "", covar_main_nms))
      find_inter_nm <- function(geno_nm, covar_lvl_suffix) {
        cand1 <- paste0(int_var, covar_lvl_suffix, ":", geno_nm)
        cand2 <- paste0(geno_nm, ":", int_var, covar_lvl_suffix)
        if (cand1 %in% all_coef_names) return(cand1)
        if (cand2 %in% all_coef_names) return(cand2)
        NA_character_
      }
      get_geno_coef_nm <- function(coef_pos) {
        if (is.na(coef_pos) || coef_pos < 1L || coef_pos > length(haplo_rows)) return(NA_character_)
        all_coef_names[haplo_rows[coef_pos]]
      }
      # Table 1: Cross-classification
      tbl_cross <- self$results$ldHaploGroup$haploGroup$haploInteractionTable
      tbl_cross$setTitle(paste0("<b>Haplotype \u00D7 ", int_var, " (cross-classification)</b>"))
      build_notes(tbl_cross)
      if (!is.na(p_inter)) tbl_cross$setNote(note = paste0("Interaction p-value (LRT): ", format.pval(p_inter, digits = 3)), key = "lrt_inter")
      for (lvl in covar_levels) {
        col_nm <- paste0("cross_", make.names(lvl))
        tbl_cross$addColumn(name = col_nm, title = paste0(lvl, if (is_binary) " OR (95%CI)" else " \u03B2 (95%CI)"), type = "text")
      }
      for (entry_nm in names(all_haplo_entries)) {
        entry   <- all_haplo_entries[[entry_nm]]
        h_label <- entry$label; h_freq <- entry$freq; cp <- entry$coef_pos
        geno_nm <- get_geno_coef_nm(cp)
        row_vals <- list(term = h_label, freq = round(h_freq, 4))
        for (lvl in covar_levels) {
          col_nm       <- paste0("cross_", make.names(lvl))
          lvl_suf      <- sub(paste0("^", int_var), "", covar_lvl_to_main[lvl] %||% "")
          is_ref_covar <- (lvl == ref_covar_lvl)
          is_ref_haplo <- is.na(cp)
          val <- if (is_ref_haplo && is_ref_covar) { if (is_binary) "1.00 (Ref)" else "0 (Ref)"
          } else if (is_ref_haplo) {
            covar_main_nm <- covar_lvl_to_main[lvl_suf]
            if (!is.null(covar_main_nm) && !is.na(covar_main_nm)) { bs <- get_beta_se(covar_main_nm); fmt_effect_ci(bs$beta, bs$se) } else NA_character_
          } else if (is_ref_covar) { bs <- get_beta_se(geno_nm); fmt_effect_ci(bs$beta, bs$se)
          } else {
            covar_main_nm <- covar_lvl_to_main[lvl_suf]
            inter_nm      <- if (!is.na(geno_nm)) find_inter_nm(geno_nm, lvl_suf) else NA_character_
            i1 <- match(geno_nm, all_coef_names); i2 <- match(covar_main_nm, all_coef_names)
            i3 <- if (!is.na(inter_nm)) match(inter_nm, all_coef_names) else NA_integer_
            if (is.na(i1) || is.na(i2)) NA_character_
            else if (is.na(i3) || is.null(vcov_mat)) {
              beta_sum <- coef_sum[i1,"coef"] + coef_sum[i2,"coef"] + if (!is.na(i3)) coef_sum[i3,"coef"] else 0
              se_sum   <- sqrt(coef_sum[i1,se_col]^2 + coef_sum[i2,se_col]^2 + if (!is.na(i3)) coef_sum[i3,se_col]^2 else 0)
              fmt_effect_ci(beta_sum, se_sum)
            } else {
              idx_used <- c(i1, i2, i3); idx_used <- idx_used[!is.na(idx_used)]
              fmt_effect_ci(sum(coef_sum[idx_used, "coef"]), sqrt(max(0, sum(vcov_mat[idx_used, idx_used]))))
            }
          }
          row_vals[[col_nm]] <- val %||% ""
        }
        tbl_cross$addRow(rowKey = paste0("cross_", entry_nm), values = row_vals)
      }
      # Table 2: Haplotype effect conditional on covariate
      tbl_cond_covar <- self$results$ldHaploGroup$haploGroup$haploCondCovarTable
      tbl_cond_covar$setTitle(paste0("<b>Haplotype effect within ", int_var, " levels</b>"))
      build_notes(tbl_cond_covar)
      for (lvl in covar_levels)
        tbl_cond_covar$addColumn(name = paste0("condcovar_", make.names(lvl)),
          title = paste0(lvl, if (is_binary) " OR (95%CI)" else " \u03B2 (95%CI)"), type = "text")
      for (entry_nm in names(all_haplo_entries)) {
        entry   <- all_haplo_entries[[entry_nm]]
        cp      <- entry$coef_pos; geno_nm <- get_geno_coef_nm(cp); is_ref_haplo <- is.na(cp)
        row_vals <- list(term = entry$label, freq = round(entry$freq, 4))
        for (lvl in covar_levels) {
          col_nm       <- paste0("condcovar_", make.names(lvl))
          lvl_suf      <- sub(paste0("^", int_var), "", covar_lvl_to_main[lvl] %||% "")
          is_ref_covar <- (lvl == ref_covar_lvl)
          val <- if (is_ref_haplo) { if (is_binary) "1.00 (Ref)" else "0 (Ref)"
          } else if (is_ref_covar) { bs <- get_beta_se(geno_nm); fmt_effect_ci(bs$beta, bs$se)
          } else {
            inter_nm <- if (!is.na(geno_nm)) find_inter_nm(geno_nm, lvl_suf) else NA_character_
            bs <- if (!is.na(inter_nm)) combine_coefs(geno_nm, inter_nm) else get_beta_se(geno_nm)
            fmt_effect_ci(bs$beta, bs$se)
          }
          row_vals[[col_nm]] <- val %||% ""
        }
        tbl_cond_covar$addRow(rowKey = paste0("cc_", entry_nm), values = row_vals)
      }
      if (!is.na(p_inter)) tbl_cond_covar$setNote(note = paste0("Interaction p-value (LRT): ", format.pval(p_inter, digits = 3)), key = "lrt_inter2")
      # Table 3: Covariate effect conditional on haplotype
      tbl_cond_haplo <- self$results$ldHaploGroup$haploGroup$haploCondHaploTable
      tbl_cond_haplo$setTitle(paste0("<b>", int_var, " effect within haplotypes</b>"))
      build_notes(tbl_cond_haplo)
      non_ref_covar_levels <- covar_levels[-1]
      ref_col_nm <- paste0("condhaplo_", make.names(ref_covar_lvl))
      tbl_cond_haplo$addColumn(name = ref_col_nm, title = paste0(ref_covar_lvl, if (is_binary) " OR (95%CI)" else " \u03B2 (95%CI)"), type = "text")
      for (lvl in non_ref_covar_levels)
        tbl_cond_haplo$addColumn(name = paste0("condhaplo_", make.names(lvl)),
          title = paste0(lvl, if (is_binary) " OR (95%CI)" else " \u03B2 (95%CI)"), type = "text")
      for (entry_nm in names(all_haplo_entries)) {
        entry   <- all_haplo_entries[[entry_nm]]
        cp      <- entry$coef_pos; geno_nm <- get_geno_coef_nm(cp); is_ref_haplo <- is.na(cp)
        row_vals <- list(term = entry$label, freq = round(entry$freq, 4))
        row_vals[[ref_col_nm]] <- if (is_binary) "1.00 (Ref)" else "0 (Ref)"
        for (lvl in non_ref_covar_levels) {
          col_nm        <- paste0("condhaplo_", make.names(lvl))
          lvl_suf       <- sub(paste0("^", int_var), "", covar_lvl_to_main[lvl] %||% "")
          covar_main_nm <- covar_lvl_to_main[lvl_suf]
          val <- if (is_ref_haplo) {
            if (!is.null(covar_main_nm) && !is.na(covar_main_nm)) { bs <- get_beta_se(covar_main_nm); fmt_effect_ci(bs$beta, bs$se) } else NA_character_
          } else {
            inter_nm <- if (!is.na(geno_nm)) find_inter_nm(geno_nm, lvl_suf) else NA_character_
            if (!is.null(covar_main_nm) && !is.na(covar_main_nm) && !is.na(inter_nm)) {
              bs <- combine_coefs(covar_main_nm, inter_nm); fmt_effect_ci(bs$beta, bs$se)
            } else if (!is.null(covar_main_nm) && !is.na(covar_main_nm)) {
              bs <- get_beta_se(covar_main_nm); fmt_effect_ci(bs$beta, bs$se)
            } else NA_character_
          }
          row_vals[[col_nm]] <- val %||% ""
        }
        tbl_cond_haplo$addRow(rowKey = paste0("ch_", entry_nm), values = row_vals)
      }
      if (!is.na(p_inter)) tbl_cond_haplo$setNote(note = paste0("Interaction p-value (LRT): ", format.pval(p_inter, digits = 3)), key = "lrt_inter2")
    }

  )  # end private
)

#' Fit association model for one SNP under one genetic model.
#' For categorical response, fits a multinomial logistic model (nnet::multinom)
#' and returns one result list per response category (vs. reference).
fit_model <- function(snp_enc, response, covariates_df, model_name,
                      response_type, ci_width) {
  df <- data.frame(resp = response, snp = snp_enc)
  cov_formula <- ""
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df          <- cbind(df, covariates_df)
    cov_formula <- paste("+", safe_rhs(names(covariates_df)))
  }
  df <- df[complete.cases(df), , drop = FALSE]
  formula_full <- as.formula(paste("resp ~ snp", cov_formula))
  formula_null <- as.formula(paste("resp ~ 1",   cov_formula))

  # ── Categorical: multinomial logistic via nnet::multinom ──────────────────
  if (response_type == "categorical") {
    if (!requireNamespace("nnet", quietly = TRUE)) return(NULL)
    df$resp <- as.factor(df$resp)
    fit_full <- tryCatch(
      nnet::multinom(formula_full, data = df, trace = FALSE),
      error = function(e) NULL)
    fit_null <- tryCatch(
      nnet::multinom(formula_null, data = df, trace = FALSE),
      error = function(e) NULL)
    if (is.null(fit_full)) return(NULL)
    lrt      <- tryCatch(anova(fit_null, fit_full), error = function(e) NULL)
    global_p <- if (!is.null(lrt) && nrow(lrt) >= 2) lrt[2, "Pr(Chi)"] else NA_real_
    aic_val  <- AIC(fit_full)
    coefs    <- summary(fit_full)$coefficients   # matrix: categories × terms
    ses      <- summary(fit_full)$standard.errors
    cats     <- rownames(coefs)                  # response categories (not ref)
    snp_cols <- grep("^snp", colnames(coefs))
    if (length(snp_cols) == 0) return(NULL)
    z_crit   <- qnorm(1 - (1 - ci_width / 100) / 2)
    # Return a list-of-lists: one entry per (category × snp_col)
    result <- list()
    for (cat in cats) {
      for (j in snp_cols) {
        beta  <- coefs[cat, j]
        se    <- ses[cat, j]
        ci_lo <- beta - z_crit * se
        ci_hi <- beta + z_crit * se
        pval  <- 2 * (1 - pnorm(abs(beta / se)))
        result[[length(result) + 1L]] <- list(
          category   = cat,
          comparison = sub("^snp", "", colnames(coefs)[j]),
          effect     = .exp_or(beta),
          ci_low     = .exp_or(ci_lo),
          ci_high    = .exp_or(ci_hi),
          pval       = pval,
          global_p   = global_p,
          aic        = aic_val,
          is_categorical = TRUE
        )
      }
    }
    return(result)
  }

  # ── Binary / quantitative ─────────────────────────────────────────────────
  tryCatch({
    if (response_type == "binary") {
      fit_full  <- glm(formula_full, data = df, family = binomial())
      fit_null  <- glm(formula_null, data = df, family = binomial())
      lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"; pval_col <- "Pr(>|z|)"
    } else {
      fit_full  <- lm(formula_full, data = df)
      fit_null  <- lm(formula_null, data = df)
      lrtest <- "F"; lrtest_label <- "Pr(>F)"; pval_col <- "Pr(>|t|)"
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
        list(effect = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
             pval = pval, global_p = global_p, aic = aic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
      else
        list(effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, global_p = global_p, aic = aic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
    })
  }, error = function(e) NULL)
}

#' Fit SNP × covariate interaction model.
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
  # Escape all user-supplied names before interpolating into formula strings
  iv_safe  <- safe_term(interaction_var)
  adj_part <- if (length(adj_covs) > 0) paste("+", safe_rhs(adj_covs)) else ""

  if (conditional) {
    if (cond_var == "snp") {
      formula_fit <- as.formula(paste("resp ~ snp /", iv_safe, adj_part))
    } else {
      formula_fit <- as.formula(paste("resp ~", iv_safe, "/ snp", adj_part))
    }
    formula_add <- as.formula(paste("resp ~ snp +", iv_safe, adj_part))
    if (response_type == "binary") {
      fit     <- glm(formula_fit, data = df, family = binomial())
      fit_add <- glm(formula_add, data = df, family = binomial())
      pval_col <- "Pr(>|z|)"; lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"
    } else {
      fit     <- lm(formula_fit, data = df)
      fit_add <- lm(formula_add, data = df)
      pval_col <- "Pr(>|t|)"; lrtest <- "F"; lrtest_label <- "Pr(>F)"
    }
    lrt_cond <- tryCatch(anova(fit_add, fit, test = lrtest), error = function(e) NULL)
    p_inter  <- lrt_cond[2, lrtest_label]
    coefs    <- summary(fit)$coefficients
    ci_mat   <- tryCatch(confint(fit, level = ci_width / 100),
                         error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2))
    aic_val  <- AIC(fit)
    all_rows <- rownames(coefs)
    inter_rows_idx <- grep(":", all_rows)
    snp_rows_idx   <- setdiff(grep("^snp", all_rows), inter_rows_idx)
    adj_rows_idx   <- setdiff(
      seq_along(all_rows),
      c(grep("^\\(Intercept\\)", all_rows),
        snp_rows_idx, inter_rows_idx,
        which(startsWith(all_rows, interaction_var))))
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
        list(term = term, effect = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
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
    formula_int  <- as.formula(paste("resp ~ snp *", iv_safe, adj_part))
    formula_main <- as.formula(paste("resp ~ snp +", iv_safe, adj_part))
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
      lrt      <- tryCatch(anova(fit_main, fit_int, test = lrtest), error = function(e) NULL)
      p_inter  <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
      aic_val  <- AIC(fit_int)
      coefs    <- summary(fit_int)$coefficients
      ci       <- tryCatch(confint(fit_int, level = ci_width / 100),
                           error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2,
                                                      dimnames = list(rownames(coefs), c("lo","hi"))))
      all_rows   <- rownames(coefs)
      snp_rows   <- grep("^snp", all_rows)
      inter_rows <- which(
        (startsWith(all_rows, "snp") & grepl(":", all_rows, fixed = TRUE) &
           sapply(all_rows, function(x) any(grepl(interaction_var, strsplit(x, ":")[[1]], fixed = TRUE)))) |
        (startsWith(all_rows, interaction_var) & grepl(":snp", all_rows, fixed = TRUE)))
      covar_rows <- setdiff(which(startsWith(all_rows, interaction_var)), inter_rows)
      adj_rows   <- setdiff(seq_along(all_rows),
                            c(grep("^\\(Intercept\\)", all_rows),
                              snp_rows, inter_rows, covar_rows))
      keep_rows  <- unique(c(snp_rows, inter_rows, covar_rows, adj_rows))
      if (length(keep_rows) == 0) return(NULL)
      result <- lapply(keep_rows, function(r) {
        beta  <- coefs[r, "Estimate"]
        pval  <- coefs[r, pval_col]
        ci_lo <- ci[r, 1]; ci_hi <- ci[r, 2]
        is_inter <- r %in% inter_rows
        row_type <- if (r %in% snp_rows)  "snp"
                    else if (r %in% inter_rows) "interaction"
                    else if (r %in% covar_rows) "covariate"
                    else                        "adjustment"
        if (response_type == "binary")
          list(term = all_rows[r], effect = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
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
  }
}
