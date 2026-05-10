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

    # ── Private state ────────────────────────────────────────────────────────
    .miss_cache = NULL,   # populated during descriptive run; read by .plotMissingness
    .ld_store   = NULL,   # LD heatmap data
    .ld_nms     = NULL,
    .ld_metric  = NULL,

    # ══════════════════════════════════════════════════════════════════════════
    # .init
    # ══════════════════════════════════════════════════════════════════════════
    .init = function() {
      private$.miss_cache <- list()

      # ── Hide every output object at startup ───────────────────────────────
      self$results$descGroup$setVisible(FALSE)
      self$results$assocGroup$setVisible(FALSE)

      self$results$descGroup$covDescGroup$setVisible(FALSE)
      self$results$descGroup$snpSummaryTablesGroup$setVisible(FALSE)
      self$results$descGroup$snpSummaryTablesGroup$snpSummaryTable$setVisible(FALSE)
      self$results$descGroup$missingnessPlot$setVisible(FALSE)

      snp_names <- self$options$snps
      if (length(snp_names) == 0) return()

      # Descriptive per-SNP array: items added here, tables hidden per-item
      desc_arr <- self$results$descGroup$descSnpResults
      for (nm in snp_names) {
        desc_arr$addItem(key = nm)
        it <- desc_arr$get(key = nm)
        it$allFreqTable$setVisible(FALSE)
        it$genoFreqTable$setVisible(FALSE)
        it$hweTable$setVisible(FALSE)
      }

      # Association per-SNP array
      opts         <- self$options
      assoc_arr    <- self$results$assocGroup$assocSnpResults
      int_models   <- private$.get_interaction_models(opts)
      model_labels <- c(codominant = "Codominant", dominant = "Dominant",
                        recessive  = "Recessive",  overdominant = "Overdominant",
                        logadditive = "Log-additive")
      for (nm in snp_names) {
        assoc_arr$addItem(key = nm)
        assoc_arr$get(key = nm)$assocTable$setVisible(FALSE)
        if (isTRUE(opts$snpInteraction) && length(int_models) > 0) {
          int_arr <- assoc_arr$get(key = nm)$interactionResults
          for (mdl in int_models) {
            int_arr$addItem(key = model_labels[[mdl]])
            mdl_it <- int_arr$get(key = model_labels[[mdl]])
            mdl_it$interactionTable$setVisible(FALSE)
            mdl_it$stratByCovariate$setVisible(FALSE)
            mdl_it$stratByGenotype$setVisible(FALSE)
            mdl_it$crossClassTable$setVisible(FALSE)
          }
        }
      }

      # haplotype / LD group
      # self$results$ldHaploGroup$ldGroup$setVisible(FALSE)
      # self$results$ldHaploGroup$haploGroup$setVisible(FALSE)
      # self$results$ldHaploGroup$ldGroup$ldTable$setVisible(FALSE)
      # self$results$ldHaploGroup$ldGroup$ldMatrixTable$setVisible(FALSE)
      # self$results$ldHaploGroup$ldGroup$ldPlotImage$setVisible(FALSE)
      # self$results$ldHaploGroup$haploGroup$haploFreqTable$setVisible(FALSE)
      # self$results$ldHaploGroup$haploGroup$haploAssocTable$setVisible(FALSE)
      # self$results$ldHaploGroup$haploGroup$haploInteractionTable$setVisible(FALSE)
      # self$results$ldHaploGroup$haploGroup$haploCondCovarTable$setVisible(FALSE)
      # self$results$ldHaploGroup$haploGroup$haploCondHaploTable$setVisible(FALSE)

    },

    # ══════════════════════════════════════════════════════════════════════════
    # .run  — single entry point; delegates ALL preprocessing to snp_prepare()
    # ══════════════════════════════════════════════════════════════════════════
    .run = function() {

      opts <- self$options

      # ── All preprocessing in one place ────────────────────────────────────
      #' @return Named list with elements:
      #'   $data, $snp_vars, $snp_data (per-SNP parsed objects),
      #'   $response_var, $response_raw, $response_type, $response_enc,
      #'   $cov_df, $complete_mask, $n_rows, $warnings
      prep <- snp_prepare(
        data           = self$data,
        snps           = opts$snps,
        response       = opts$response,
        covariates     = opts$covariates,
        response_type  = opts$responseType %||% "auto",
        rm_snp_missing = isTRUE(opts$rmSnpMissing)
      )

      # Validation messages
      if (nchar(prep$warnings) > 0) {
        self$results$validationMsg$setContent(prep$warnings)
        self$results$validationMsg$setVisible(TRUE)
      } else {
        self$results$validationMsg$setVisible(FALSE)
      }

      # activate visualization
      #
      # ── Descriptive ───────────────────────────────────────────────────────
      private$.run_descriptive(prep, opts)

      if(isTRUE(opts$covDesc) || isTRUE(opts$snpSummary)) {
        self$results$descGroup$setVisible(TRUE)       
      }
      # Covariate descriptives subgroup
      if(isTRUE(opts$covDesc) && !is.null(prep$cov_df) && ncol(prep$cov_df) > 0) {
        self$results$descGroup$covDescGroup$setVisible(TRUE)
      } # else{
      #   self$results$descGroup$covDescGroup$setVisible(FALSE)
      # } 
      
      # ── Association ───────────────────────────────────────────────────────
      if (isTRUE(opts$snpAssoc) || isTRUE(opts$snpInteraction)) {
        self$results$assocGroup$setVisible(TRUE)

        if (is.null(prep$response_var) || prep$response_var == "") {
          self$results$validationMsg$setContent(
            "<p style='color:red;'>A response variable is required for association analysis.</p>")
          self$results$validationMsg$setVisible(TRUE)
        } else {
          private$.run_association(prep, opts)
        }
      }

      # ── LD / Haplotype ────────────────────────────────────────────────────
      any_ld <- isTRUE(opts$ldAnalysis) || isTRUE(opts$ldMatrix) || isTRUE(opts$ldPlot) ||
                isTRUE(opts$haploFreq)  || isTRUE(opts$haploAssoc) || isTRUE(opts$haploInteraction)
      if (any_ld) {
        if (length(prep$snp_vars) < 2) {
          self$results$validationMsg$setContent(
            "<p style='color:red;'>LD and haplotype analyses require at least 2 SNPs.</p>")
          self$results$validationMsg$setVisible(TRUE)
        } else {
          private$.run_ldhaplo(prep, opts)
        }
      }
    },

    # ══════════════════════════════════════════════════════════════════════════
    # .run_descriptive  — calls compute_* then writes to jamovi tables
    # ══════════════════════════════════════════════════════════════════════════
    .run_descriptive = function(prep, opts) {
      run_snpSummary  <- isTRUE(opts$snpSummary)
      run_allFreq     <- isTRUE(opts$allFreq)
      run_genoFreq    <- isTRUE(opts$genoFreq)
      run_hweTest     <- isTRUE(opts$hweTest)
      run_subpop      <- isTRUE(opts$subpop) &&
                         !is.null(prep$response_raw) &&
                         prep$response_type %in% c("binary", "categorical")
      run_covDesc     <- isTRUE(opts$covDesc)
      run_showMissing <- isTRUE(opts$showMissing)

      private$.miss_cache <- list()
      res <- self$results$descGroup

      # Validation
      res$validationMsgSNP$setVisible(FALSE)
      res$validationMsgGeno$setVisible(FALSE)

      # ── Covariate descriptives ─────────────────────────────────────────────
      if (run_covDesc && (!is.null(prep$cov_df) || !is.null(prep$response_raw))) {
        result <- compute_cov_desc(prep, subpop = run_subpop)
        if (!is.null(result)) {
          res$setVisible(TRUE)
          res$covDescGroup$setVisible(TRUE)
          res$covDescGroup$covDescTable$setVisible(TRUE)
          private$.write_cov_desc(res$covDescGroup$covDescTable, result, prep)
        }
      }

      # ── SNP summary table ──────────────────────────────────────────────────
      if (! length(prep$snp_vars) > 0) return() # no SNPs, nothing to do here
      
      if (run_snpSummary) {
        result <- compute_snp_summary(prep, subpop = run_subpop)
        if (!is.null(result)) {
          res$setVisible(TRUE)
          res$snpSummaryTablesGroup$setVisible(TRUE)
          res$snpSummaryTablesGroup$snpSummaryTable$setVisible(TRUE)
          private$.write_snp_summary(res$snpSummaryTablesGroup$snpSummaryTable,
                                     result, prep)
        }
      }

      # ── Per-SNP descriptives ───────────────────────────────────────────────
      if (run_allFreq || run_genoFreq || run_hweTest) {
        res$setVisible(TRUE)
      }

      null_pat               <- "^0[/|>]0$|^00$"
      total_null_across_snps <- sum(sapply(prep$snp_vars, function(nm) {
        raw <- as.character(prep$data[[nm]])
        sum(!is.na(raw) & grepl(null_pat, raw, ignore.case = TRUE))
      }))
      arr <- res$descSnpResults

      for (nm in prep$snp_vars) {
        sd   <- prep$snp_data[[nm]]
        item <- arr$get(key = nm)

        n_typed          <- sd$n_typed
        n_total_eligible <- sum(prep$complete_mask)
        total_missing    <- n_total_eligible - n_typed

        # Typing-rate HTML
        typing_html <- sprintf("<b>Typed samples:</b> %d / %d (%.1f%%)",
          n_typed, n_total_eligible,
          if (n_total_eligible > 0) n_typed / n_total_eligible * 100 else 0)
        if (total_missing > 0)
          typing_html <- paste0(typing_html, sprintf(
            " ── <b>Missing SNP:</b> %d (%.1f%%)",
            total_missing,
            if (n_total_eligible > 0) total_missing / n_total_eligible * 100 else 0))
        item$typingRate$setContent(typing_html)

        # Missingness cache (used by .plotMissingness)
        n_miss_by_level <- if (run_subpop && !is.null(prep$response_raw)) {
          resp_levels <- levels(as.factor(prep$response_raw))
          v <- sapply(resp_levels, function(lvl)
            sum(is.na(sd$clean) & prep$complete_mask &
                !is.na(prep$response_raw) & prep$response_raw == lvl))
          names(v) <- resp_levels
          v
        } else NULL

        private$.miss_cache[[nm]] <- list(
          n_total_eligible = n_total_eligible,
          total_missing    = total_missing,
          n_miss_by_level  = n_miss_by_level)

        # Allele frequency
        if (run_allFreq) {
          result <- compute_allele_freq(nm, prep, subpop = run_subpop,
                                        show_missing = run_showMissing)
          item$allFreqTable$setVisible(TRUE)
          private$.write_allele_freq(item$allFreqTable, result, run_subpop, prep)
        }

        # Genotype frequency
        if (run_genoFreq) {
          result <- compute_geno_freq(nm, prep, subpop = run_subpop,
                                      show_missing = run_showMissing)
          item$genoFreqTable$setVisible(TRUE)
          private$.write_geno_freq(item$genoFreqTable, result, run_subpop,
                                   prep, prep$response_type)
        }

        # HWE
        if (run_hweTest) {
          result <- compute_hwe(nm, prep, subpop = run_subpop,
                                show_missing = run_showMissing)
          if (!is.null(result)) {
            item$hweTable$setVisible(TRUE)
            private$.write_hwe(item$hweTable, result)
          }
        }
      }

      # Missingness plot
      res$missingnessPlot$setVisible(
        isTRUE(opts$showMissingnessPlot) && length(private$.miss_cache) > 0)

      # Null-allele note
      if (run_snpSummary) {
        tbl <- res$snpSummaryTablesGroup$snpSummaryTable
        tbl$setNote(key = "null_allele",
          note = if (total_null_across_snps > 0)
            paste0(total_null_across_snps,
                   " genotype(s) coded as 0/0 were treated as missing (NA).")
          else NULL)
      }
    },

    # ══════════════════════════════════════════════════════════════════════════
    # Table writers — thin jamovi-specific adapters over compute_* results
    # ══════════════════════════════════════════════════════════════════════════

    .write_cov_desc = function(tbl, result, prep) {
      tbl_df     <- result$table
      grp_levels <- result$grp_levels   # NULL when not stratified
      do_strat   <- !is.null(grp_levels) && length(grp_levels) > 0

      # Set group column titles and make them visible
      if (do_strat) {
        for (j in seq_along(grp_levels)) {
          nm <- paste0("stat_g", j - 1L)
          tbl$getColumn(nm)$setTitle(grp_levels[j])
          tbl$getColumn(nm)$setVisible(TRUE)
        }
        tbl$getColumn("pval")$setVisible(TRUE)
      }

      for (i in seq_len(nrow(tbl_df))) {
        # Column names in tbl_df already match jamovi column names
        # (stat_g0, stat_g1, … for groups; overall → stat_overall).
        row <- list(
          variable     = as.character(tbl_df[["variable"]][i]),
          level        = as.character(tbl_df[["level"]][i]),
          stat_overall = as.character(tbl_df[["overall"]][i])
        )
        if (do_strat) {
          for (j in seq_along(grp_levels)) {
            nm        <- paste0("stat_g", j - 1L)
            row[[nm]] <- as.character(tbl_df[[nm]][i])
          }
          row$pval <- fmt_pval(tbl_df[["pval"]][i])
        }
        tbl$addRow(rowKey = as.character(i), values = row)
      }

      if (length(result$notes) > 0)
        tbl$setNote(key = "cov_desc_note",
                    note = paste(result$notes, collapse = " "))
    },

    .write_snp_summary = function(tbl, result, prep) {
      do_strat <- prep$response_type %in% c("binary","categorical") &&
                  !is.null(prep$response_raw)
      tbl$getColumn("group")$setVisible(do_strat)
      for (i in seq_len(nrow(result))) {
        tbl$addRow(rowKey = as.character(i), values = list(
          snp        = as.character(result[["snp"]][i]),
          alleles    = as.character(result[["alleles"]][i]),
          group      = as.character(result[["group"]][i]),
          n          = as.integer(result[["n"]][i]),
          missing    = if (is.na(result[["missing"]][i])) '' else as.integer(result[["missing"]][i]),
          maf        = fmt3(result[["maf"]][i]),
          genoCounts = as.character(result[["geno_counts"]][i]),
          hwePval    = fmt_pval(result[["hwe_pval"]][i])))
      }
      parts <- c(if (!is.null(prep$cov_df) && ncol(prep$cov_df)>0) "covariates",
                 if (!is.null(prep$response_raw)) "response")
      if (length(parts) > 0)
        tbl$setNote(key = "missing_resp_cov",
          note = paste0("Complete cases used: rows missing any ",
                        paste(parts, collapse = " or "), " or SNP value are excluded."))
    },

    .write_allele_freq = function(tbl, result, do_strat, prep) {
      grp_levels <- if (do_strat && !is.null(prep$response_raw))
                      levels(as.factor(prep$response_raw)) else NULL
      if (do_strat && !is.null(grp_levels))
        for (j in seq_along(grp_levels))
          tbl$addColumn(name = paste0("stat_g", j-1L), title = grp_levels[j],
                        type = "text")
      for (i in seq_len(nrow(result))) {
        row <- list(
          allele = as.character(result[["allele"]][i]),
          stat   = as.character(result[["overall"]][i]))
        if (do_strat && !is.null(grp_levels))
          for (j in seq_along(grp_levels))
            row[[paste0("stat_g", j-1L)]] <-
              as.character(result[[paste0("stat_g", j-1L)]][i] %||% "")
        tbl$addRow(rowKey = as.character(i), values = row)
      }
    },

    .write_geno_freq = function(tbl, result, do_strat, prep, response_type) {
      grp_levels <- if (do_strat && !is.null(prep$response_raw))
                      levels(as.factor(prep$response_raw)) else NULL
      if (response_type == "quantitative")
        tbl$getColumn("responseStat")$setVisible(TRUE)
      if (do_strat && !is.null(grp_levels))
        for (j in seq_along(grp_levels))
          tbl$addColumn(name = paste0("stat_g", j-1L), title = grp_levels[j],
                        type = "text")
      for (i in seq_len(nrow(result))) {
        row <- list(
          genotype     = as.character(result[["genotype"]][i]),
          stat         = as.character(result[["overall"]][i]),
          responseStat = as.character(result[["response_stat"]][i] %||% ""))
        if (do_strat && !is.null(grp_levels))
          for (j in seq_along(grp_levels))
            row[[paste0("stat_g", j-1L)]] <-
              as.character(result[[paste0("stat_g", j-1L)]][i] %||% "")
        tbl$addRow(rowKey = as.character(i), values = row)
      }
    },

    .write_hwe = function(tbl, result) {
      labels <- result$col_labels
      if (length(labels) == 3) {
        tbl$getColumn("n11")$setTitle(labels[1])
        tbl$getColumn("n12")$setTitle(labels[2])
        tbl$getColumn("n22")$setTitle(labels[3])
      }
      rows <- result$rows
      for (i in seq_len(nrow(rows))) {
        tbl$addRow(rowKey = as.character(i), values = list(
          group   = as.character(rows[["group"]][i]),
          n11     = as.integer(rows[["n11"]][i]),
          n12     = as.integer(rows[["n12"]][i]),
          n22     = as.integer(rows[["n22"]][i]),
          missing = if (is.na(rows[["missing"]][i])) '' else as.integer(rows[["missing"]][i]),
          pval    = fmt_pval(rows[["pval"]][i])))
      }
    },

    # ══════════════════════════════════════════════════════════════════════════
    # .run_association  — calls compute_assoc then writes to jamovi tables
    # ══════════════════════════════════════════════════════════════════════════
    .run_association = function(prep, opts) {
      run_snpAssoc       <- isTRUE(opts$snpAssoc)
      run_snpInteraction <- isTRUE(opts$snpInteraction)

      if (run_snpInteraction && (is.null(prep$cov_df) || ncol(prep$cov_df) == 0)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>SNP \u00D7 covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }
      if (run_snpInteraction && prep$response_type == "categorical") {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>SNP \u00D7 covariate interaction is not available for categorical responses.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }

      int_models   <- private$.get_interaction_models(opts)
      model_labels <- c(codominant="Codominant", dominant="Dominant",
                        recessive="Recessive", overdominant="Overdominant",
                        logadditive="Log-additive")
      assoc_arr    <- self$results$assocGroup$assocSnpResults

      for (nm in prep$snp_vars) {
        sd      <- prep$snp_data[[nm]]
        item    <- assoc_arr$get(key = nm)
        n_miss  <- prep$n_rows - sd$n_typed

        # Typing-rate HTML
        typing_html <- sprintf("<b>Typed samples:</b> %d / %d (%.1f%%)",
          sd$n_typed, prep$n_rows,
          if (prep$n_rows > 0) sd$n_typed / prep$n_rows * 100 else 0)
        if (n_miss > 0)
          typing_html <- paste0(typing_html, sprintf(
            " ── <b>Missing:</b> %d", n_miss))
        item$typingRate$setContent(typing_html)

        # Main association table
        if (run_snpAssoc) {
          result <- compute_assoc(nm, prep,
                                  models   = int_models,
                                  ci_width = opts$ciWidth %||% 95)
          if (!is.null(result$rows) && nrow(result$rows) > 0) {
            item$assocTable$setVisible(TRUE)
            private$.write_assoc(item$assocTable, result, opts, nm)
          }
        }

        # Interaction tables
        if (run_snpInteraction && !is.null(prep$cov_df) && ncol(prep$cov_df) > 0) {
          for (mdl in int_models) {
            mdl_item    <- item$interactionResults$get(key = model_labels[[mdl]])
            interaction_var <- names(prep$cov_df)[1]  # first covariate as interaction var
            int_lbl     <- attr(self$data[[interaction_var]], "label") %||% interaction_var
            snp_char    <- as.character(sd$clean_cc)
            response_cc <- prep$response_enc[sd$snp_mask]
            cov_df_cc   <- prep$cov_df[sd$snp_mask, , drop = FALSE]
            ref         <- sd$ref
            user_levels <- sd$user_levels

            if (isTRUE(opts$showInteractionTable)) {
              mdl_item$interactionTable$setVisible(TRUE)
              private$.fill_interaction(
                mdl_item$interactionTable, sd$clean_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                prep$response_type, opts, mdl, user_levels,
                prep$response_raw[sd$snp_mask], nm)
            }
            if (isTRUE(opts$showStratByCovariate)) {
              mdl_item$stratByCovariate$setVisible(TRUE)
              mdl_item$stratByCovariateHeading$setContent(
                paste0("<h3>Stratified by Covariate: ", int_lbl, "</h3>"))
              private$.fill_strat_by_covariate(
                mdl_item$stratByCovariate, sd$clean_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                prep$response_type, opts, mdl, user_levels,
                prep$response_raw[sd$snp_mask], nm)
            }
            if (isTRUE(opts$showStratByGenotype)) {
              mdl_item$stratByGenotype$setVisible(TRUE)
              mdl_item$stratByGenotypeHeading$setContent(
                paste0("<h3>Stratified by Genotype: ", nm, "</h3>"))
              private$.fill_strat_by_genotype(
                mdl_item$stratByGenotype, sd$clean_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                prep$response_type, opts, mdl, user_levels,
                prep$response_raw[sd$snp_mask], nm)
            }
            if (isTRUE(opts$showCrossClassTable)) {
              mdl_item$crossClassTable$setVisible(TRUE)
              mdl_item$crossClassHeading$setContent(
                paste0("<h3>Cross-Classification: ", nm, " \u00D7 ", int_lbl, "</h3>"))
              private$.fill_cross_class(
                mdl_item$crossClassTable, sd$clean_cc, ref,
                response_cc, cov_df_cc, interaction_var,
                prep$response_type, opts, mdl, user_levels,
                prep$response_raw[sd$snp_mask], nm)
            }
          }
        }
      }
    },

    # Write compute_assoc result to a jamovi assocTable
    .write_assoc = function(tbl, result, opts, snp_nm) {
      rtype <- result$response_type
      tbl$getColumn("effect")$setTitle(result$col_titles$effect)
      tbl$getColumn("genotype")$setTitle(snp_nm)
      resp_lbl <- attr(self$data[[self$options$response]], "label") %||%
                  self$options$response
      tbl$setTitle(paste0("Association with ", resp_lbl))

      if (rtype == "binary") {
        tbl$getColumn("stat0")$setTitle(result$col_titles$stat0)
        tbl$getColumn("stat1")$setTitle(result$col_titles$stat1)
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(TRUE)
      } else if (rtype == "categorical") {
        tbl$getColumn("stat0")$setTitle("N (%)")
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)")
        tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }
      tbl$getColumn("AIC")$setVisible(isTRUE(opts$showAIC))
      tbl$getColumn("BIC")$setVisible(isTRUE(opts$showAIC))

      for (note in result$notes)
        tbl$setNote(key = note, note = note)

      rows <- result$rows
      for (i in seq_len(nrow(rows))) {
        tbl$addRow(rowKey = as.character(i), values = list(
          model    = as.character(rows[["model"]][i]),
          genotype = as.character(rows[["genotype"]][i]),
          stat0    = as.character(rows[["stat0"]][i]),
          stat1    = as.character(rows[["stat1"]][i]),
          effect   = fmt3(rows[["effect"]][i]),
          ciLow    = fmt3(rows[["ci_low"]][i]),
          ciHigh   = fmt3(rows[["ci_high"]][i]),
          pval     = fmt_pval(rows[["pval"]][i]),
          AIC      = fmt3(rows[["aic"]][i]),
          BIC      = fmt3(rows[["bic"]][i])))
      }
    },

    # ══════════════════════════════════════════════════════════════════════════
    # .run_ldhaplo  — calls compute_ld / compute_haplo_* then writes tables
    # ══════════════════════════════════════════════════════════════════════════
    .run_ldhaplo = function(prep, opts) {
      run_ldAnalysis      <- isTRUE(opts$ldAnalysis)
      run_ldMatrix        <- isTRUE(opts$ldMatrix)
      run_ldPlot          <- isTRUE(opts$ldPlot)
      run_haploFreq       <- isTRUE(opts$haploFreq)
      run_haploAssoc      <- isTRUE(opts$haploAssoc)
      run_haploInteraction <- isTRUE(opts$haploInteraction)
      run_subpop          <- isTRUE(opts$subpop) &&
                             !is.null(prep$response_raw) &&
                             prep$response_type %in% c("binary","categorical")

      ld_res_grp <- self$results$ldHaploGroup

      # ── LD ───────────────────────────────────────────────────────────────
      if (run_ldAnalysis || run_ldMatrix || run_ldPlot) {
        ld_res_grp$ldGroup$setVisible(TRUE)
        private$.run_ld(prep, opts,
                        run_ldAnalysis, run_ldMatrix, run_ldPlot)
      }

      # ── Haplotype ─────────────────────────────────────────────────────────
      if (run_haploFreq || run_haploAssoc || run_haploInteraction) {
        ld_res_grp$haploGroup$setVisible(TRUE)
        private$.run_haplo(prep, opts,
                           run_haploFreq, run_haploAssoc, run_haploInteraction,
                           run_subpop)
      }
    },

    # ══════════════════════════════════════════════════════════════════════════
    # Shared association utilities (unchanged from original)
    # ══════════════════════════════════════════════════════════════════════════
    .get_interaction_models = function(opts) {
      c(if (isTRUE(opts$modelCodominant))   "codominant",
        if (isTRUE(opts$modelDominant))     "dominant",
        if (isTRUE(opts$modelRecessive))    "recessive",
        if (isTRUE(opts$modelOverdominant)) "overdominant",
        if (isTRUE(opts$modelLogAdditive))  "logadditive")
    },

    .geno_labels_for_model = function(model, all_genos, ref) {
      .geno_labels_for_model(model, all_genos, ref)
    },

    .split_genos = function(gl) .split_genos(gl),

    .compute_stats = function(geno_labels, snp_char, response, response_type,
                              response_raw = NULL) {
      .compute_stats(geno_labels, snp_char, response, response_type, response_raw)
    },

    .bic_from_aic = function(aic_val, mdl, n_fit, n_cov) {
      .bic_from_aic(aic_val, mdl, n_fit, n_cov)
    },

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
            aic_val <- if (first_row && !is.nan(res$aic)) fmt3(res$aic) else ""
            bic_val <- if (first_row && !is.null(res$bic) && !is.nan(res$bic))
              fmt3(res$bic) else ""
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
              effect = fmt3(level_res[[1]]$effect), ciLow = fmt3(level_res[[1]]$ci_low),
              ciHigh = fmt3(level_res[[1]]$ci_high), pval = fmt_pval(level_res[[1]]$pval)))
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
              effect = fmt3(res$effect), ciLow = fmt3(res$ci_low), ciHigh = fmt3(res$ci_high), pval = fmt_pval(res$pval)))
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
                effect = if (!is.null(res)) fmt3(res$effect) else if (response_type == "binary") 1.0 else 0.0,
                ciLow  = if (!is.null(res)) fmt3(res$ci_low)  else "",
                ciHigh = if (!is.null(res)) fmt3(res$ci_high) else "",
                pval   = if (!is.null(res)) fmt_pval(res$pval) else ""))
            }
          } else {
            stat0 <- if (response_type == "binary") fmt_cont(int_g[resp_raw_g == resp_lv[1]]) else fmt_cont(resp_g)
            stat1 <- if (response_type == "binary") fmt_cont(int_g[resp_raw_g == resp_lv[2]]) else ""
            res   <- if (length(gl_res) > 0) gl_res[[1]] else NULL
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              level  = "Overall", stat0 = stat0, stat1 = stat1,
              effect = if (!is.null(res)) fmt3(res$effect) else if (response_type == "binary") 1.0 else 0.0,
              ciLow  = if (!is.null(res)) fmt3(res$ci_low)  else "",
              ciHigh = if (!is.null(res)) fmt3(res$ci_high) else "",
              pval   = if (!is.null(res)) fmt_pval(res$pval) else ""))
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
              effect = fmt3(if (response_type == "binary") .exp_or(combined_beta) else combined_beta),
              ciLow  = fmt3(if (response_type == "binary") .exp_or(lo_beta) else lo_beta),
              ciHigh = fmt3(if (response_type == "binary") .exp_or(hi_beta) else hi_beta),
              pval   = fmt_pval(p_val)))
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
            r2     = fmt3(ld_res$`r`^2),
            Dprime = fmt3(ld_res$`D'`),
            D      = fmt3(ld_res$`D`),
            pval   = fmt_pval(ld_res$`P-value`)))
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
          p_str  <- fmt_pval(p_val)
          up_val <- switch(metric,
            Dprime = fmt3(ld_res$`D'`),
            r2     = fmt3(ld_res$`r`^2),
            D      = fmt3(ld_res$`D`))
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
          df$label[k] <- fmt_pval(pv)
        } else if (i_idx < j_idx) {
          key <- paste(c(nms[min(i_idx,j_idx)], nms[max(i_idx,j_idx)]), collapse="___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) {
            raw <- switch(metric, r2=ld_res$`r`^2, Dprime=ld_res$`D'`, D=ld_res$`D`)
            df$label[k] <- fmt3(as.numeric(raw))
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
            setNames(as.list(lapply(em_g$hap.prob, fmt3)),
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
          row_vals <- list(haplotype=label, freq=fmt3(freqs[i]))
          if (do_strat_haplo) {
            row_vals$freq_g0 <- grp_freq[[grp_levels_haplo[1]]][[label]] %||% NA_real_
            row_vals$freq_g1 <- grp_freq[[grp_levels_haplo[2]]][[label]] %||% NA_real_
          }
          tbl$addRow(rowKey=paste0("f",i), values=row_vals)
        }
        if (rare_sum > 0) {
          row_vals <- list(haplotype=paste0("Rare (<",opts$haploFreqMin,")"), freq=fmt3(rare_sum))
          if (do_strat_haplo) {
            em0 <- em_grp[[grp_levels_haplo[1]]]; em1 <- em_grp[[grp_levels_haplo[2]]]
            rare_g0 <- if (!is.null(em0)) fmt3(sum(em0$hap.prob[em0$hap.prob < opts$haploFreqMin])) else NA_real_
            rare_g1 <- if (!is.null(em1)) fmt3(sum(em1$hap.prob[em1$hap.prob < opts$haploFreqMin])) else NA_real_
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
          list(haplotype = label, freq = fmt3(freq),
               effect = if (response_type == "binary") .exp_or(b)  else b,
               ciLow  = if (response_type == "binary") .exp_or(lo) else lo,
               ciHigh = if (response_type == "binary") .exp_or(hi) else hi,
               pval   = fmt_pval(stats$pval))
        }
        base_idx   <- haplo_fit$haplo.base
        base_label <- label_from_unique_row(haplo_fit$haplo.unique[base_idx, ])
        base_freq  <- haplo_fit$haplo.freq[base_idx]
        tbl$addRow(rowKey = "base", values = list(
          haplotype = paste0(base_label, " (Ref)"), freq = fmt3(base_freq),
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
        row_vals <- list(term = h_label, freq = fmt3(h_freq))
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
        row_vals <- list(term = entry$label, freq = fmt3(entry$freq))
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
        row_vals <- list(term = entry$label, freq = fmt3(entry$freq))
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
    bic_val  <- BIC(fit_full)
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
          bic        = bic_val,
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
    bic_val  <- BIC(fit_full)
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
             pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
      else
        list(effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
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
    bic_val  <- BIC(fit_full)
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
          bic        = bic_val,
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
    bic_val  <- BIC(fit_full)
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
             pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
      else
        list(effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
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
