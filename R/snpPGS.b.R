# snpPGS_b.R  — Analysis class for SNPstats PGS module
#
# Weight / allele information flow:
#   (a) weightsPath set  — parse the catalog file; unknown columns are
#       concatenated into the extra_cols field of the results table.
#   (b) No file          — unit weights (weight = 1) for all selected SNPs.
#
# The weights file is re-read on every .run() call, so toggling reloadWeights
# or editing weightsPath always picks up the latest file contents.
# ─────────────────────────────────────────────────────────────────────────────

source("R/snp_helpers.R")

snpPGSClass <- R6::R6Class(
  "snpPGSClass",
  inherit = snpPGSBase,

  private = list(

    .keepMask  = NULL,
    .cache     = new.env(parent = emptyenv()),

    # ════════════════════════════════════════════════════════════════════════
    # .run
    # ════════════════════════════════════════════════════════════════════════
    .run = function() {

      snpCols     <- self$options$snpCols
      covCols     <- self$options$covCols
      respCol     <- self$options$responseCol
      missing_st  <- self$options$missingStrategy
      wmode       <- self$options$weightingMode

      # ── Nothing selected yet: show guidance, hide everything ────────────
      if (is.null(snpCols) || length(snpCols) == 0) {
        self$results$validationMsg$setContent(
          "<div style='color:#555; padding:6px 0;'>
             <b>Getting started:</b><br>
             \u2022 Drag one or more SNP columns into <i>SNP columns</i>.<br>
             \u2022 Optionally load a PGS Catalog weights file (.csv / .tsv) for weighted scoring.<br>
             \u2022 Optionally select a response variable to test association.
           </div>")
        self$results$validationMsg$setVisible(TRUE)
        self$results$snpGridTable$setVisible(FALSE)
        self$results$coverageTable$setVisible(FALSE)
        self$results$summaryTable$setVisible(FALSE)
        self$results$percentileTable$setVisible(FALSE)
        self$results$assocTable$setVisible(FALSE)
        self$results$distPlot$setVisible(FALSE)
        self$results$stratPlot$setVisible(FALSE)
        return()
      }

      # Reset keepMask at the start of each run so a previous 'exclude'
      # strategy mask is never silently applied if buildDosageMatrix returns early.
      private$.keepMask <- rep(TRUE, nrow(self$data))

      # ── Build weight table from file (or unit weights) ───────────────────
      wtable <- private$.buildWeightTable(snpCols)
      if (is.null(wtable)) return()

      # ── Dosage matrix + allele QC ────────────────────────────────────────
      qc     <- private$.buildDosageMatrix(snpCols, wtable, missing_st)
      if (is.null(qc)) return()
      dosage <- qc$mat
      wtable <- qc$wtable
      valid_snps <- qc$valid_snps

      # Note: qc$valid_counts is already row-filtered inside buildDosageMatrix
      # for the 'exclude' strategy, so it always matches qc$mat row count on
      # return. No additional filtering is needed here.

      private$.fillSnpGridTable(wtable, valid_snps)

      # ── Weights aligned to dosage columns ───────────────────────────────
      matched_rows  <- wtable[wtable$rsid %in% colnames(dosage), , drop = FALSE]
      catalog_wvec  <- setNames(matched_rows$effect_weight, matched_rows$rsid)
      catalog_wvec  <- catalog_wvec[colnames(dosage)]
      unit_wvec     <- setNames(rep(1, ncol(dosage)), colnames(dosage))

      # Determine which modes to run
      has_file    <- !is.null(self$options$weightsPath) &&
                     nchar(trimws(self$options$weightsPath)) > 0 &&
                     file.exists(self$options$weightsPath)
      run_weighted   <- wmode %in% c("weighted", "both") && has_file
      run_unweighted <- wmode %in% c("unweighted", "both") || !has_file

      modes <- list()
      if (run_weighted)   modes[["Weighted"]]   <- catalog_wvec
      if (run_unweighted) modes[["Unweighted"]] <- unit_wvec

      # ── Row indices after keepMask ───────────────────────────────────────
      keep <- private$.keepMask  # NULL or logical vector length nrow(data)

      # ── Resolve response vector ──────────────────────────────────────────
      resp <- NULL
      if (!is.null(respCol) && nchar(respCol) > 0 && respCol %in% names(self$data)) {
        resp_full <- self$data[[respCol]]
        resp <- if (!is.null(keep) && length(keep) == length(resp_full))
                  resp_full[keep] else resp_full
        if (is.numeric(resp) && length(unique(resp[!is.na(resp)])) == 2 &&
            all(resp[!is.na(resp)] %in% c(0, 1)))
          resp <- factor(resp)
      }

      # ── Resolve covariate data frame ─────────────────────────────────────
      covs <- NULL
      if (!is.null(covCols) && length(covCols) > 0) {
        valid_covs <- covCols[covCols %in% names(self$data)]
        if (length(valid_covs) > 0) {
          cov_full <- self$data[, valid_covs, drop = FALSE]
          covs <- if (!is.null(keep) && length(keep) == nrow(cov_full))
                    cov_full[keep, , drop = FALSE]
                  else cov_full
        }
      }

      private$.fillCoverageTable(snpCols, wtable, missing_st, valid_snps,
                                  valid_counts  = qc$valid_counts,
                                  n_total       = nrow(self$data),
                                  has_file      = has_file,
                                  n_qc_excl     = qc$n_qc_excl %||% 0L,
                                  n_flt_excl    = qc$n_flt_excl %||% 0L)

      # ── Clear tables before multi-mode fill ─────────────────────────────
      self$results$summaryTable$deleteRows()
      self$results$assocTable$deleteRows()
      self$results$interactionTable$deleteRows()

      all_scores <- list()   # named by mode_label, for saveScores output

      for (mode_label in names(modes)) {
        wvec <- modes[[mode_label]]

        na_snps <- names(wvec)[is.na(wvec)]
        if (length(na_snps) > 0) {
          keep_cols <- names(wvec)[!is.na(wvec)]
          wvec  <- wvec[keep_cols]
          dos_m <- dosage[, keep_cols, drop = FALSE]
        } else {
          dos_m <- dosage
        }

        if (length(wvec) == 0) next

        # Slice valid_counts to exactly the SNPs in dos_m (weighted mode may
        # have dropped NA-weight SNPs). Pass as a lightweight list rather than
        # copying the full qc object.
        vc_sliced <- list(
          valid_counts = qc$valid_counts[,
            intersect(colnames(qc$valid_counts), colnames(dos_m)), drop = FALSE]
        )
        scores <- private$.computeScores(dos_m, wvec, vc_sliced,
                                          is_unweighted = (mode_label == "Unweighted"))
        all_scores[[mode_label]] <- scores


        private$.fillSummaryTable(scores, resp, mode_label)

        if (self$options$showAssoc && !is.null(resp))
          private$.fillAssocTable(scores, resp, respCol, covs, mode_label)

        if (self$options$showInteraction && !is.null(resp) &&
            !is.null(covs) && ncol(covs) > 0)
          private$.fillInteractionTable(scores, resp, respCol, covs, mode_label)
      }

      if (length(all_scores) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>No SNPs with valid weights — cannot compute PGS.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      # ── Cache everything the plot render functions need ─────────────────
      if (self$options$showDistPlot) {
        private$.cache$all_scores <- all_scores
        private$.cache$resp       <- resp
        private$.cache$respCol    <- respCol
      }
      # stratPlot (continuous scatter) visible only when response is continuous
      is_binary_resp <- !is.null(resp) && (is.factor(resp) ||
                        length(unique(resp[!is.na(resp)])) == 2)
      show_strat <- self$options$showDistPlot && !is.null(resp) && !is_binary_resp
      self$results$stratPlot$setVisible(show_strat)

      if (self$options$showPercentiles)
        private$.fillPercentileTable(all_scores, resp, respCol, covs)

      if (self$results$saveScores$isNotFilled())
        private$.saveScoresToData(all_scores)
    },


    # ════════════════════════════════════════════════════════════════════════
    # .scaleLabel
    # Returns a short human-readable description of the active scoring pipeline
    # used in table notes and plot axis labels.
    # ════════════════════════════════════════════════════════════════════════
    .scaleLabel = function(mode_label) {
      mc  <- isTRUE(self$options$missingCorrection)
      sm  <- self$options$scaleMethod
      sf  <- self$options$scaleFactor %||% 10
      sw  <- isTRUE(self$options$scaleWeights)
      std <- isTRUE(self$options$standardize)

      scale_str <- switch(sm,
        proportion  = if (mc) "proportion of max [0–1]"
                      else    "raw sum / theoretical max",
        percent     = if (mc) "percent of max [0–100]"
                      else    "raw sum / theoretical max × 100",
        multiply    = paste0(if (mc) "corrected" else "raw", " × ", sf),
        perNAlleles = paste0("per ", sf, " risk alleles"),
        "raw weighted sum"
      )
      wt_str  <- if (sw && mode_label == "Weighted") " [weights L1-scaled]" else ""
      std_str <- if (std) ", SD-standardized" else ""
      mc_str  <- if (!mc && sm %in% c("proportion", "percent"))
                   " [WARNING: missingness correction off]" else ""
      paste0(mode_label, " — ", scale_str, wt_str, std_str, mc_str)
    },

    # ════════════════════════════════════════════════════════════════════════
    # .buildWeightTable
    #
    # Returns a data.frame with columns:
    #   rsid, effect_allele, other_allele, effect_weight, chr, pos,
    #   matched (logical), extra_cols (character — concatenated extra fields)
    #
    # If no file is configured, returns unit-weight rows for snpCols.
    # ════════════════════════════════════════════════════════════════════════
    .buildWeightTable = function(snpCols) {
      path  <- self$options$weightsPath

      has_file <- !is.null(path) && nchar(trimws(path)) > 0 && file.exists(path)

      if (has_file)
        return(private$.parseCatalogFile(path, snpCols))
      
      private$.unitWeightTableFromData(snpCols)
    },

    # ════════════════════════════════════════════════════════════════════════
    # .computeScores
    #
    # Pipeline (in order):
    #   1. Optional weight scaling  — L1-normalise wvec to unit mean so the
    #      weighted score is on the same numeric scale as the unweighted score.
    #   2. Raw score                — dosage %*% wvec  (per-individual)
    #   3. Missing correction       — divide by per-individual theoretical max
    #      (positive weights only; falls back to global max if valid_counts
    #      unavailable). This is the prerequisite for proportional rescaling.
    #   4. Rescaling                — one of:
    #        none        : keep the missingness-corrected (or raw) score
    #        proportion  : score is already in [0,1] after correction; no-op
    #                      (correction IS the proportion step)
    #        percent     : × 100  → [0, 100]
    #        multiply    : × scaleFactor  (user-supplied multiplier)
    #        perNAlleles : raw_score / (scaleFactor × denominator_per_snp)
    #                      for unweighted: denominator_per_snp = 1
    #                      for weighted:   denominator_per_snp = mean(|pos weights|)
    #   5. Standardize              — divide by SD across individuals (SD = 1)
    #
    # is_unweighted: signals unit weights so we use SNP-count denominator.
    # ════════════════════════════════════════════════════════════════════════
    .computeScores = function(dosage, wvec, qc, is_unweighted = FALSE) {

      missing_corr <- isTRUE(self$options$missingCorrection)
      scale_method <- self$options$scaleMethod    # none/proportion/percent/multiply/perNAlleles
      scale_factor <- self$options$scaleFactor
      scale_wts    <- isTRUE(self$options$scaleWeights) && !is_unweighted
      standardize  <- isTRUE(self$options$standardize)

      if (is.null(scale_factor) || is.na(scale_factor) || scale_factor == 0)
        scale_factor <- 10

      # ── Step 1: optional weight scaling ──────────────────────────────────
      # L1-normalise to mean(|wvec|) = 1 so weighted and unweighted scores
      # share a common scale ("effective risk-allele count").
      if (scale_wts) {
        mean_abs <- mean(abs(wvec), na.rm = TRUE)
        if (!is.na(mean_abs) && mean_abs > 0)
          wvec <- wvec / mean_abs
      }

      # ── Step 2: raw weighted sum ──────────────────────────────────────────
      scores <- as.numeric(dosage %*% wvec)

      # ── Step 3: per-individual maximum (used by correction and proportion) ─
      # Always compute this; the rescaling branch needs it even when
      # missingCorrection = FALSE (for 'proportion'/'percent' to be valid).
      vc_cols <- intersect(names(wvec), colnames(qc$valid_counts))

      if (length(vc_cols) > 0) {
        vc <- qc$valid_counts[, vc_cols, drop = FALSE]   # logical TRUE/FALSE matrix
        if (is_unweighted || scale_wts) {
          # After weight-scaling, all effective weights are ≈ 1, so use SNP count.
          # For truly unweighted: max = 2 × n_observed_SNPs.
          max_possible <- 2 * rowSums(vc)
        } else {
          # Max weighted = 2 × sum of positive weights for observed SNPs.
          # Negative weights reduce the score, so they cannot be part of the
          # maximum — including them would underestimate the denominator.
          w_sub        <- wvec[vc_cols]
          w_pos        <- pmax(w_sub, 0)
          max_possible <- 2 * as.numeric(vc %*% w_pos)
        }
      } else {
        # valid_counts unavailable — fall back to global (same for all individuals)
        if (is_unweighted || scale_wts) {
          max_possible <- rep(2 * length(wvec), nrow(dosage))
        } else {
          max_possible <- rep(2 * sum(pmax(wvec, 0), na.rm = TRUE), nrow(dosage))
        }
      }
      max_possible[max_possible == 0] <- NA_real_

      # ── Step 4: missingness correction ───────────────────────────────────
      # Divide by per-individual maximum to remove the effect of missingness.
      # Required for 'proportion' and 'percent' to be valid.
      if (missing_corr) {
        scores <- scores / max_possible
        # After correction, scores ≈ proportion in [0, 1] (may exceed 1 slightly
        # when some weights are negative, which is mathematically correct).
      }

      # ── Step 5: rescaling ─────────────────────────────────────────────────
      scores <- switch(scale_method,

        # 'proportion': correction is the proportion step — already done.
        # No further multiplication needed; this is the identity branch.
        proportion = scores,

        # 'percent': convert proportion to 0–100.
        percent = scores * 100,

        # 'multiply': user-supplied multiplier applied to the corrected score.
        # If missingCorrection = FALSE, this multiplies the raw sum.
        multiply = scores * scale_factor,

        # 'perNAlleles': express the score per N risk alleles.
        #   Unweighted: raw_sum / N  →  "dosage per N SNPs"
        #   Weighted: raw_sum / (N × mean_pos_weight)  →  comparable unit
        # If missingCorrection was applied, undo it first so we work on the
        # raw scale, then apply the per-N denominator.
        perNAlleles = {
          # Reconstruct the raw (uncorrected) score so the per-N denominator
          # is applied to the actual weighted sum, not the proportion.
          # When missing_corr=TRUE, scores == raw/max_possible, so
          # raw = scores * max_possible.  NAs propagate correctly because
          # max_possible is NA wherever the score was already NA.
          raw <- if (missing_corr) scores * max_possible else scores
          if (is_unweighted || scale_wts) {
            raw / scale_factor
          } else {
            mean_pos_w <- mean(pmax(wvec, 0), na.rm = TRUE)
            if (is.na(mean_pos_w) || mean_pos_w == 0) mean_pos_w <- 1
            raw / (scale_factor * mean_pos_w)
          }
        },

        # 'none': return the (optionally corrected) raw score as-is.
        scores
      )

      # ── Step 6: standardize to SD = 1 ────────────────────────────────────
      if (standardize) {
        s <- sd(scores, na.rm = TRUE)
        if (!is.na(s) && s > 0) scores <- scores / s
      }

      scores
    },

    # ════════════════════════════════════════════════════════════════════════
    # .parseCatalogFile
    # ════════════════════════════════════════════════════════════════════════
    .parseCatalogFile = function(path, snpCols) {

      raw <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
      if (is.null(raw)) {
        self$results$validationMsg$setContent(
          paste0("<p style='color:#c0392b;'>Cannot read weights file: ", path, "</p>"))
        self$results$validationMsg$setVisible(TRUE)
        return(private$.unitWeightTable(snpCols))
      }

      meta      <- private$.parseFileMetadata(raw)
      dataLines <- raw[!grepl("^#", raw) & nchar(trimws(raw)) > 0]
      if (length(dataLines) < 2) return(private$.unitWeightTable(snpCols))

      # Auto-detect delimiter from the header line: tab > semicolon > comma
      first    <- dataLines[1]
      n_tabs   <- lengths(regmatches(first, gregexpr("\t", first)))
      n_semis  <- lengths(regmatches(first, gregexpr(";",  first)))
      n_commas <- lengths(regmatches(first, gregexpr(",",  first)))
      sep <- if (n_tabs >= n_semis && n_tabs >= n_commas) "\t"
             else if (n_semis >= n_commas)                ";"
             else                                         ","

      parse_err <- NULL
      df <- tryCatch(
        read.table(text = paste(dataLines, collapse = "\n"),
                   header = TRUE, sep = sep,
                   stringsAsFactors = FALSE, quote = "\"",
                   fill = TRUE, comment.char = "", check.names = FALSE),
        error = function(e) { parse_err <<- conditionMessage(e); NULL }
      )
      if (is.null(df) || nrow(df) == 0) {
        msg <- if (!is.null(parse_err))
          paste0("<p style='color:#c0392b;'>Failed to parse weights file (sep='",
                 sep, "'): ", parse_err, "</p>")
        else
          "<p style='color:#c0392b;'>Weights file parsed to an empty table.</p>"
        self$results$validationMsg$setContent(msg)
        self$results$validationMsg$setVisible(TRUE)
        return(private$.unitWeightTable(snpCols))
      }

      orig_names  <- names(df)
      lower_names <- tolower(orig_names)

      # Helper: find first matching column name (case-insensitive)
      find_col <- function(...) {
        for (cn in tolower(c(...)))
          if (cn %in% lower_names) return(orig_names[match(cn, lower_names)])
        NULL
      }

      c_rsid   <- find_col("rsid", "snp", "snp_id","marker")
      c_ea     <- find_col("effect_allele", "risk_allele", "effect","risk", "alt_allele")
      c_oa     <- find_col("other_allele", "ref_allele", "non_effect_allele", "reference_allele", "common_allele","reference")
      c_weight <- find_col("effect_weight", "beta", "weight", "or")
      c_chr    <- find_col("chr_name", "chromosome", "chr", "chrom")
      c_pos    <- find_col("chr_position", "position", "pos", "bp")

      if (is.null(c_rsid)) {
        self$results$validationMsg$setContent(paste0(
          "<p style='color:#c0392b;'>Weights file has no recognisable rsID column.</p>",
          "<p>Columns found: <code>", paste(orig_names, collapse = ", "), "</code></p>",
          "<p>Expected one of: rsID, variant_id, snp, snp_id, marker_name</p>"))
        self$results$validationMsg$setVisible(TRUE)
        return(private$.unitWeightTable(snpCols))
      }

      # Known columns that we handle explicitly
      known_cols <- c(c_rsid, c_ea, c_oa, c_weight, c_chr, c_pos)
      known_cols <- known_cols[!sapply(known_cols, is.null)]
      extra_cols <- setdiff(orig_names, known_cols)

      # Build extra_cols string per row: "col1=val1; col2=val2"
      extra_str <- if (length(extra_cols) > 0) {
        apply(df[, extra_cols, drop = FALSE], 1, function(row) {
          paste(paste0(extra_cols, "=", row), collapse = "; ")
        })
      } else rep("", nrow(df))

      ea_vec <- if (!is.null(c_ea)) toupper(trimws(as.character(df[[c_ea]]))) else rep("", nrow(df))
      oa_vec <- if (!is.null(c_oa)) toupper(trimws(as.character(df[[c_oa]]))) else rep("", nrow(df))

      catalog <- data.frame(
        rsid          = as.character(df[[c_rsid]]),
        effect_allele = ea_vec,
        other_allele  = oa_vec,
        effect_weight = if (!is.null(c_weight)) suppressWarnings(as.numeric(df[[c_weight]]))
                        else                    rep(NA_real_, nrow(df)),
        chr           = if (!is.null(c_chr))    as.character(df[[c_chr]])    else "",
        pos           = if (!is.null(c_pos))    as.character(df[[c_pos]])    else "",
        matched       = TRUE,
        allele_status = "",   # filled in by .buildDosageMatrix QC step
        strand_flipped = FALSE,
        extra_cols    = extra_str,
        n_missing     = NA_integer_,
        pct_missing   = NA_real_,
        effect_af     = NA_real_,
        hwe_p         = NA_real_,
        stringsAsFactors = FALSE
      )


      # Keep ALL catalog SNPs, plus add placeholder for any selected SNPs not in catalog
      result <- catalog

      # Mark which SNPs are selected (present in snpCols)
      result$selected_flag <- result$rsid %in% snpCols

      # Add placeholder rows for selected SNPs not in catalog
      missing_from_catalog <- setdiff(snpCols, catalog$rsid)
      if (length(missing_from_catalog) > 0) {
        extra_rows <- data.frame(
          rsid = missing_from_catalog, 
          effect_allele = "", 
          other_allele = "",
          effect_weight = NA_real_, 
          chr = "", 
          pos = "",
          matched = FALSE, 
          allele_status = "\u274c not in weights file",
          strand_flipped = FALSE, 
          extra_cols = "",
          n_missing = NA_integer_, 
          pct_missing = NA_real_,
          effect_af = NA_real_,
          hwe_p = NA_real_,
          selected_flag = TRUE,  # These are selected but not in catalog
          stringsAsFactors = FALSE
        )
        result <- rbind(result, extra_rows)
      }

      # Now also mark catalog SNPs that are NOT selected (but we still show them)
      result$selected_flag[!result$rsid %in% snpCols] <- FALSE

      # Sort: selected first, then unselected
      result <- result[order(-result$selected_flag, result$rsid), ]

      # Add a visual indicator in allele_status for unselected SNPs
      result$allele_status <- ifelse(
        !result$selected_flag & result$matched,
        paste0(result$allele_status, " (in catalog but not selected)"),
        result$allele_status
      )

      # Remove the temporary flag column before returning
      result$selected_flag <- NULL
      
      attr(result, "pgs_meta") <- meta
      result
    },


    # ────────────────────────────────────────────────────────────────────────
    .parseFileMetadata = function(raw_lines) {
      meta <- list(pgs_id = "", pgs_name = "", trait_reported = "",
                   weight_type = "", genome_build = "", variants_number = "")
      for (ln in raw_lines) {
        if (!grepl("^#", ln)) break
        inner <- sub("^#+", "", ln)
        eq    <- regexpr("=", inner)
        if (eq == -1) next
        k <- trimws(substr(inner, 1, eq - 1))
        v <- trimws(substr(inner, eq + 1, nchar(inner)))
        if (tolower(k) %in% names(meta)) meta[[tolower(k)]] <- v
      }
      meta
    },

    .unitWeightTable = function(snpCols) {
      data.frame(
        rsid           = snpCols,
        effect_allele  = "",
        other_allele   = "",
        effect_weight  = 1,
        chr            = "",
        pos            = "",
        matched        = FALSE,
        allele_status  = "\u2705 no weights file (unit weights)",
        strand_flipped = FALSE,
        extra_cols     = "",
        n_missing      = NA_integer_,
        pct_missing    = NA_real_,
        effect_af      = NA_real_,
        hwe_p          = NA_real_,
        stringsAsFactors = FALSE
      )
    },

    # Builds wtable from the observed data when no catalog file is loaded.
    # Allele detection logic:
    #   - Numeric columns: treated as dosage (0/1/2); alleles left blank.
    #   - Factor/character columns: genotype strings are parsed to extract
    #     unique single-base alleles. The reference (other) allele is taken
    #     from the homozygous genotype of the first factor level (e.g. "C/C"
    #     → ref = C). The effect allele is the remaining allele.
    #   QC flags: monomorphic (only one allele), non-biallelic (>2 alleles).
    .unitWeightTableFromData = function(snpCols) {
      dat  <- self$data

      rows <- lapply(snpCols, function(snp) {
        col <- dat[[snp]]
        ea  <- ""
        oa  <- ""

        if (is.numeric(col)) {
          # Numeric dosage: allele inference not possible
          qc <- "\u26A0 numeric dosage: no allele QC possible"
        } else {
          # infer alleles (best effort, no QC decisions)
          # Apply the same cleaning as .cleanGenotypeColumn so null-allele
          # codings (0/0 etc.) don't pollute allele inference.
          cleaned_inf <- private$.cleanGenotypeColumn(col)
          col_char    <- cleaned_inf$col_clean

          bases_list <- strsplit(col_char[!is.na(col_char)],
                                "[^ACGT]|(?<=.)(?=.)", perl = TRUE)
          bases_list <- lapply(bases_list, function(b) b[b %in% c("A","C","G","T")])
          alleles <- sort(unique(unlist(bases_list)))

          if (length(alleles) >= 1) {
            oa <- alleles[1]
            ea <- if (length(alleles) >= 2) alleles[2] else alleles[1]
          }

          qc <- "\u2705 ok (unit weight)"
        }        
        data.frame(
          rsid           = snp,
          effect_allele  = ea,
          other_allele   = oa,
          effect_weight  = 1,
          chr            = "",
          pos            = "",
          matched        = TRUE,
          allele_status  = qc,
          strand_flipped = FALSE,
          extra_cols     = "",
          n_missing      = NA_integer_,
          pct_missing    = NA_real_,
          effect_af      = NA_real_,
          hwe_p          = NA_real_,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    },


    # ════════════════════════════════════════════════════════════════════════
    # .complement  —  single-base DNA complement
    # ════════════════════════════════════════════════════════════════════════
    .complement = function(alleles) {
      chartr("ACGT", "TGCA", toupper(alleles))
    },

    # ════════════════════════════════════════════════════════════════════════
    # .cleanGenotypeColumn
    #
    # Pre-processes a raw genotype column before any QC or dosage conversion.
    # Handles the following solvable encoding problems in-place:
    #
    #   1. Null-allele patterns  — "0/0", "0|0", "00", "0" — coerced to NA.
    #      These appear when genotyping pipelines encode missing genotypes as
    #      zero rather than NA or "./." and are NOT valid biallelic calls.
    #
    #   2. Canonical NA strings  — "", "NA", "N/A", ".", "./.", "?/?" — → NA.
    #
    #   3. Whitespace trimming and uppercase normalisation.
    #
    #   4. Numeric columns are left untouched (handled by the dosage branch).
    #
    # Returns a character vector (same length as input) with corrections applied.
    # A named list is returned:
    #   $col_clean  : character vector, corrected values
    #   $n_null     : number of cells fixed from null-allele coding
    #   $n_na_str   : number of cells fixed from NA-string coding
    # ════════════════════════════════════════════════════════════════════════
    .cleanGenotypeColumn = function(col_raw) {

      # Numeric columns: return as-is (dosage branch does its own cleaning)
      if (is.numeric(col_raw)) {
        return(list(col_clean = col_raw, n_null = 0L, n_na_str = 0L))
      }

      col_char <- toupper(trimws(as.character(col_raw)))

      # ── Step 1: canonical NA strings ────────────────────────────────────────
      na_strings <- c("", "NA", "N/A", ".", "./.", "?/?", "?", "-/-", "-")
      is_na_str  <- col_char %in% na_strings
      n_na_str   <- sum(is_na_str & !is.na(col_raw))  # only those not already NA
      col_char[is_na_str] <- NA_character_

      # ── Step 2: null-allele patterns (0/0, 0|0, 00, lone 0) ────────────────
      # Match strings that contain ONLY zeros and separators (/ | space),
      # with no ACGT bases whatsoever.
      is_null_allele <- !is.na(col_char) &
                        grepl("^[0/|[:space:]]+$", col_char) &
                        !grepl("[ACGT]", col_char)
      n_null <- sum(is_null_allele)
      col_char[is_null_allele] <- NA_character_

      list(col_clean = col_char, n_null = n_null, n_na_str = n_na_str)
    },

    # ════════════════════════════════════════════════════════════════════════
    # .isNumericDosage
    #
    # Decides whether a column should be treated as numeric dosage (0/1/2)
    # rather than as allele strings.
    #
    # Robustness fix over the old "> 0.5 numeric fraction" heuristic:
    #   - A column is numeric-dosage ONLY if it is already stored as numeric,
    #     OR if it is a factor/character whose non-NA values ALL convert to
    #     a valid dosage integer (0, 1, or 2) with no ACGT characters present.
    #   - Mixed columns ("A/A", "1", NA) are treated as allele strings so
    #     that the allele-parsing branch can flag them correctly rather than
    #     silently dropping the letter genotypes as NA.
    # ════════════════════════════════════════════════════════════════════════
    .isNumericDosage = function(col) {
      if (is.numeric(col)) return(TRUE)

      col_str <- toupper(trimws(as.character(col)))
      obs     <- col_str[!is.na(col_str) & col_str != "" & col_str != "NA"]
      if (length(obs) == 0) return(FALSE)

      # Reject if ANY observed value contains an ACGT base character
      if (any(grepl("[ACGT]", obs))) return(FALSE)

      # Accept only if every observed value is exactly 0, 1, or 2
      num_try <- suppressWarnings(as.integer(obs))
      all(!is.na(num_try) & num_try %in% 0:2)
    },

    # ════════════════════════════════════════════════════════════════════════
    # .buildDosageMatrix
    #
    # Full pipeline:
    #   1. Per-SNP column cleaning (.cleanGenotypeColumn)
    #   2. Input-type detection (.isNumericDosage — robust version)
    #   3. Allele QC: mismatch, multiallelic, monomorphic, no-allele-info
    #   4. Dosage conversion
    #   5. Per-SNP statistics: n_missing, pct_missing, effect_af, hwe_p
    #      (computed on the CLEANED column, before imputation)
    #   6. Configurable QC filters: missingness threshold, HWE p threshold
    #      (applied AFTER step 5, so stats are always shown)
    #   7. Missing-value imputation for keep_snps
    #   8. Row exclusion if missing_st == "exclude"
    # ════════════════════════════════════════════════════════════════════════
    .buildDosageMatrix = function(snpCols, wtable, missing_st) {

      useCols <- intersect(wtable$rsid, names(self$data))

      if (length(useCols) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>None of the selected SNP columns are present in the dataset.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return(NULL)
      }

      # ── Read user QC thresholds ────────────────────────────────────────────
      miss_thresh <- self$options$qcMaxMissingPct   # numeric, e.g. 10
      hwe_thresh  <- self$options$qcHweP            # numeric, e.g. 0.001
      apply_miss  <- isTRUE(self$options$qcFilterMissing) && !is.na(miss_thresh)
      apply_hwe   <- isTRUE(self$options$qcFilterHwe)    && !is.na(hwe_thresh)

      # ── Resolve HWE control vector ────────────────────────────────────────
      # When a binary response variable is selected, HWE is automatically
      # tested in controls only (first/lowest level), the standard approach
      # in case-control GWAS QC.  No separate column is needed.
      ctrl_vec  <- NULL
      ctrl_note <- ""
      if (apply_hwe) {
        resp_col_nm <- self$options$responseCol
        if (!is.null(resp_col_nm) && nchar(trimws(resp_col_nm)) > 0 &&
            resp_col_nm %in% names(self$data)) {
          rv      <- self$data[[resp_col_nm]]
          rv_vals <- unique(rv[!is.na(rv)])
          if (length(rv_vals) == 2) {
            ctrl_levels <- if (is.factor(rv)) levels(rv)
                           else sort(rv_vals)
            ctrl_ref  <- ctrl_levels[1]   # lower value / first factor level = controls
            ctrl_vec  <- rv
            ctrl_note <- paste0(" (controls: ", resp_col_nm, "=", ctrl_ref, ")")
          }
        }
      }

      n_rows <- nrow(self$data)

      mat <- matrix(NA_real_,
                    nrow = n_rows,
                    ncol = length(useCols),
                    dimnames = list(NULL, useCols))

      # Track exclusion reasons separately so stats are always computed first
      qc_exclude    <- character(0)   # SNPs excluded by allele/monomorphic QC
      filter_exclude <- character(0)  # SNPs excluded by missingness/HWE filters

      # Per-SNP cleaned columns (stored so HWE/AF can use cleaned data)
      cleaned_cols  <- vector("list", length(useCols))
      names(cleaned_cols) <- useCols

      for (snp in useCols) {

        col_raw <- self$data[[snp]]
        idx     <- which(wtable$rsid == snp)[1]

        ea_cat <- toupper(as.character(wtable$effect_allele[idx]))
        oa_cat <- toupper(as.character(wtable$other_allele[idx]))
        has_allele_info <- !is.na(ea_cat) && nchar(trimws(ea_cat)) > 0 &&
                           ea_cat != "0" && oa_cat != "0"

        # ── Step 1: clean the column ─────────────────────────────────────────
        cleaned     <- private$.cleanGenotypeColumn(col_raw)
        col_clean   <- cleaned$col_clean
        n_null_fix  <- cleaned$n_null
        cleaned_cols[[snp]] <- col_clean

        # Single helper: appends null-allele note or "" (used at every status site)
        null_sfx <- function(n) {
          if (n > 0) paste0("; ", n, " null-allele genotypes set to NA") else ""
        }

        # ── Step 2: all missing after cleaning ───────────────────────────────
        if (all(is.na(col_clean))) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0("\u274c all missing", null_sfx(n_null_fix))
          next
        }

        # ── Step 3: determine input type ─────────────────────────────────────
        # Use the robust detector: numeric ONLY when values are strictly 0/1/2
        # with no ACGT characters anywhere.
        if (private$.isNumericDosage(col_clean)) {

          col_num <- suppressWarnings(as.numeric(as.character(col_clean)))
          # Out-of-range values (should not exist after cleaning, but guard anyway)
          invalid <- !is.na(col_num) & (col_num < 0 | col_num > 2)
          if (any(invalid)) col_num[invalid] <- NA_real_

          mat[, snp] <- col_num

          obs_vals <- sort(unique(col_num[!is.na(col_num)]))
          obs_str  <- paste0("dosage(", paste(obs_vals, collapse = ","), ")")

          if (length(obs_vals) <= 1) {
            qc_exclude <- c(qc_exclude, snp)
            wtable$allele_status[idx] <- paste0(
              "\u274c monomorphic numeric dosage: ", obs_str, null_sfx(n_null_fix))
          } else {
            # \u26A0: allele QC cannot be performed on numeric dosage columns
            n_oor <- sum(!is.na(col_num) & (col_num < 0 | col_num > 2))
            num_actions <- c(
              "no allele QC possible (numeric dosage)",
              if (n_null_fix > 0) paste0(n_null_fix, " null-allele genotypes set to NA"),
              if (n_oor      > 0) paste0(n_oor,      " out-of-range values set to NA")
            )
            wtable$allele_status[idx] <- paste0(
              "\u26A0 ", obs_str, ": ", paste(num_actions, collapse = "; "))
          }
          next
        }

        # ── Step 4: allele-string SNPs ────────────────────────────────────────
        # col_clean is already character, NAs already set above.
        # Parse observed bases (ACGT only).
        bases_list <- strsplit(col_clean[!is.na(col_clean)],
                               "[^ACGT]|(?<=.)(?=.)", perl = TRUE)
        bases_list <- lapply(bases_list, function(b) b[b %in% c("A","C","G","T")])
        alleles    <- sort(unique(unlist(bases_list)))

        if (length(alleles) == 0) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0(
            "\u274c no valid alleles observed", null_sfx(n_null_fix))
          next
        }

        obs_str <- paste(alleles, collapse = "/")

        # ── Multiallelic ─────────────────────────────────────────────────────
        if (length(alleles) > 2) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0(
            "\u274c multiallelic: ", obs_str, null_sfx(n_null_fix))
          next
        }

        # ── No allele info in weights file ───────────────────────────────────
        if (!has_allele_info) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0(
            "\u274c no allele info in weights file: obs=", obs_str, null_sfx(n_null_fix))
          next
        }

        # ── Allele matching ──────────────────────────────────────────────────
        ea_comp <- private$.complement(ea_cat)
        oa_comp <- private$.complement(oa_cat)

        matches_direct <-
          length(alleles) == 2 && all(alleles %in% c(ea_cat, oa_cat))
        matches_complement <-
          length(alleles) == 2 && all(alleles %in% c(ea_comp, oa_comp))

        if (!matches_direct && !matches_complement) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0(
            "\u274c allele mismatch: obs=", obs_str, ", exp=", ea_cat, "/", oa_cat,
            null_sfx(n_null_fix))
          next
        }

        # ── Strand handling ──────────────────────────────────────────────────
        # Collect action flags; icon decided after all checks on this SNP.
        strand_flipped_snp <- FALSE
        is_ambiguous_snp   <- FALSE

        if (matches_direct) {
          ea_use <- ea_cat
        } else {
          ea_use <- ea_comp
          strand_flipped_snp <- TRUE
          wtable$strand_flipped[idx] <- TRUE
        }

        if (private$.isAmbiguous(ea_cat, oa_cat))
          is_ambiguous_snp <- TRUE

        # ── Dosage computation ───────────────────────────────────────────────
        # Each genotype string is split on non-ACGT separators.
        # Genotypes with != 2 alleles → NA (handles residual oddities).
        dos <- vapply(col_clean, function(g) {
          if (is.na(g)) return(NA_real_)
          b <- strsplit(g, "[^ACGT]|(?<=.)(?=.)", perl = TRUE)[[1]]
          b <- b[b %in% c("A","C","G","T")]
          if (length(b) != 2) return(NA_real_)
          sum(b == ea_use)
        }, numeric(1))

        mat[, snp] <- dos

        # ── Monomorphism check ───────────────────────────────────────────────
        if (length(unique(dos[!is.na(dos)])) <= 1) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0(
            "\u274c monomorphic: ", obs_str, null_sfx(n_null_fix))
        } else {
          # \u2705 clean pass  /  \u26A0 kept with actions
          actions_snp <- c(
            if (strand_flipped_snp) "strand flipped",
            if (is_ambiguous_snp)   "ambiguous AT/CG",
            if (n_null_fix > 0)     paste0(n_null_fix, " null-allele genotypes set to NA")
          )
          wtable$allele_status[idx] <-
            if (length(actions_snp) == 0)
              paste0("\u2705 ", obs_str)
            else
              paste0("\u26A0 ", obs_str, ": ", paste(actions_snp, collapse = "; "))
        }
      }  # end per-SNP loop

      # ── Step 5: per-SNP statistics (on cleaned data, before imputation) ────
      # Computed for ALL SNPs (including qc_exclude) so the table always shows
      # complete statistics even for excluded SNPs.
      for (snp in colnames(mat)) {
        idx  <- which(wtable$rsid == snp)[1]
        n_na <- sum(is.na(mat[, snp]))
        wtable$n_missing[idx]   <- n_na
        wtable$pct_missing[idx] <- round(n_na / n_rows * 100, 1)

        # AF and HWE from the cleaned column.
        # For HWE in controls: subset cleaned column to control individuals.
        ea <- as.character(wtable$effect_allele[idx])
        col_for_stats <- cleaned_cols[[snp]]

        # ctrl_vec and ctrl_note were resolved above from the binary response variable.
        # When a binary response is present, HWE is computed in controls only.
        if (!is.null(ctrl_vec) && length(ctrl_vec) == n_rows) {
          ctrl_levels <- if (is.factor(ctrl_vec)) levels(ctrl_vec)
                         else sort(unique(ctrl_vec[!is.na(ctrl_vec)]))
          ctrl_ref    <- ctrl_levels[1]
          ctrl_mask   <- !is.na(ctrl_vec) & as.character(ctrl_vec) == as.character(ctrl_ref)
          col_hwe     <- col_for_stats
          col_hwe[!ctrl_mask] <- NA
        } else {
          col_hwe <- col_for_stats
        }

        res <- tryCatch(
          snp_af_hwe(col_hwe, effect_allele = ea),
          error = function(e) list(effect_af = NA_real_, hwe_p = NA_real_)
        )
        wtable$effect_af[idx] <- res$effect_af
        wtable$hwe_p[idx]     <- res$hwe_p
      }

      # ── Step 6: configurable QC filters ───────────────────────────────────
      # Applied after statistics are computed, so stats remain visible.
      # Only SNPs that passed allele QC are eligible for these filters.
      candidate_snps <- setdiff(colnames(mat), qc_exclude)

      for (snp in candidate_snps) {
        idx <- which(wtable$rsid == snp)[1]

        # Missingness filter
        if (apply_miss) {
          pct_miss <- wtable$pct_missing[idx]
          if (!is.na(pct_miss) && pct_miss > miss_thresh) {
            filter_exclude <- c(filter_exclude, snp)
            # Replace leading icon with \u274c; preserve the existing detail after it
            prev <- sub("^\\\\u[0-9a-f]{4}\\s*", "", wtable$allele_status[idx])
            wtable$allele_status[idx] <- paste0(
              sprintf("\u274c excl (missing %.1f%% > %.1f%%): ", pct_miss, miss_thresh),
              prev)
            next  # no need to check HWE
          }
        }

        # HWE filter
        if (apply_hwe) {
          hwe_p_val <- wtable$hwe_p[idx]
          if (!is.na(hwe_p_val) && hwe_p_val < hwe_thresh) {
            filter_exclude <- c(filter_exclude, snp)
            prev <- sub("^\\\\u[0-9a-f]{4}\\s*", "", wtable$allele_status[idx])
            wtable$allele_status[idx] <- paste0(
              sprintf("\u274c excl (HWE p=%.2e < %.2e%s): ", hwe_p_val, hwe_thresh, ctrl_note),
              prev)
          }
        }
      }

      all_exclude <- unique(c(qc_exclude, filter_exclude))
      keep_snps   <- setdiff(colnames(mat), all_exclude)

      if (length(keep_snps) == 0) {
        n_qc  <- length(qc_exclude)
        n_flt <- length(filter_exclude)
        msg   <- paste0(
          "<p style='color:#c0392b;'>No SNPs passed QC filters.</p>",
          if (n_qc  > 0) paste0("<p>", n_qc,  " SNP(s) excluded by allele/monomorphic QC.</p>") else "",
          if (n_flt > 0) paste0("<p>", n_flt, " SNP(s) excluded by missingness/HWE thresholds.</p>") else ""
        )
        self$results$validationMsg$setContent(msg)
        self$results$validationMsg$setVisible(TRUE)
        return(NULL)
      }

      # valid_counts: logical matrix — TRUE where dosage was observed (pre-imputation)
      valid_counts <- !is.na(mat[, keep_snps, drop = FALSE])

      # ── Step 7: missing imputation for kept SNPs ───────────────────────────
      for (snp in keep_snps) {
        na_mask <- is.na(mat[, snp])
        if (any(na_mask)) {
          mat[na_mask, snp] <- switch(missing_st,
            mean       = { m <- mean(mat[, snp], na.rm = TRUE); if (is.nan(m)) 0 else m },
            zero       = 0,
            "SNP-wise" = 0,
            exclude    = NA_real_
          )
        }
      }

      # ── Step 8: row exclusion for missing_st == "exclude" ─────────────────
      if (missing_st == "exclude") {
        keep_rows    <- complete.cases(mat[, keep_snps, drop = FALSE])
        mat          <- mat[keep_rows, , drop = FALSE]
        valid_counts <- valid_counts[keep_rows, , drop = FALSE]
        private$.keepMask <- keep_rows
      } else {
        private$.keepMask <- rep(TRUE, n_rows)
      }

      mat <- mat[, keep_snps, drop = FALSE]

      list(
        mat          = mat,
        wtable       = wtable,
        valid_counts = valid_counts,
        valid_snps   = keep_snps,
        n_qc_excl    = length(qc_exclude),
        n_flt_excl   = length(filter_exclude)
      )
    },

    # ════════════════════════════════════════════════════════════════════════
    # Output helpers
    # ════════════════════════════════════════════════════════════════════════

    .fillSnpGridTable = function(wtable, valid_snps = character(0)) {
      tbl    <- self$results$snpGridTable
      tbl$deleteRows()
      inData <- names(self$data)

      if (isTRUE(self$options$filterValidSnps)) {
        wtable <- wtable[wtable$rsid %in% valid_snps, , drop = FALSE]
      }

      # Hide columns that are entirely empty (field not present in file)
      has_chr   <- any(nchar(wtable$chr)           > 0, na.rm = TRUE)
      has_pos   <- any(nchar(wtable$pos)           > 0, na.rm = TRUE)
      has_ea    <- any(nchar(wtable$effect_allele) > 0, na.rm = TRUE)
      has_oa    <- any(nchar(wtable$other_allele)  > 0, na.rm = TRUE)
      has_extra <- any(nchar(wtable$extra_cols)    > 0, na.rm = TRUE)
      has_miss  <- any(!is.na(wtable$n_missing))
      has_af    <- any(!is.na(wtable$effect_af))
      has_hwe   <- any(!is.na(wtable$hwe_p))
      # Show QC filter column only when at least one filter is active
      has_qc_flt <- isTRUE(self$options$qcFilterMissing) ||
                    isTRUE(self$options$qcFilterHwe)

      tbl$getColumn("chr")$setVisible(has_chr)
      tbl$getColumn("pos")$setVisible(has_pos)
      tbl$getColumn("effect_allele")$setVisible(has_ea)
      tbl$getColumn("other_allele")$setVisible(has_oa)
      tbl$getColumn("extra_cols")$setVisible(has_extra)
      tbl$getColumn("n_missing")$setVisible(has_miss)
      tbl$getColumn("pct_missing")$setVisible(has_miss)
      tbl$getColumn("effect_af")$setVisible(has_af)
      tbl$getColumn("hwe_p")$setVisible(has_hwe)
      tbl$getColumn("qc_excl_reason")$setVisible(has_qc_flt)

      for (i in seq_len(nrow(wtable))) {
        r      <- wtable[i, ]
        status <- as.character(r$allele_status)
        in_ds  <- r$rsid %in% inData

        # QC filter column: derive from the leading icon in allele_status.
        # \u2705 = clean pass, \u26A0 = kept with actions, \u274c = excluded.
        qc_flt_reason <- if (grepl("^\\u274c excl", status, fixed = FALSE)) {
          sub("^\\u274c excl ([^:]+):.*", "\\u274c excl \\1", status)
        } else if (r$rsid %in% valid_snps) {
          if (startsWith(status, "\\u26A0")) "\\u26A0 pass (with actions)"
          else "\u2705 pass"
        } else {
          ""   # not in dataset or excluded before filter stage
        }

        tbl$addRow(rowKey = i, values = list(
          rsid           = as.character(r$rsid),
          chr            = as.character(r$chr),
          pos            = as.character(r$pos),
          effect_allele  = as.character(r$effect_allele),
          other_allele   = as.character(r$other_allele),
          effect_weight  = if (is.na(r$effect_weight)) '' else r$effect_weight,
          matched        = if (in_ds) "\u2713" else "\u2717",
          allele_status  = status,
          extra_cols     = as.character(r$extra_cols),
          n_missing      = if (is.na(r$n_missing))   '' else as.integer(r$n_missing),
          pct_missing    = if (is.na(r$pct_missing)) '' else r$pct_missing,
          effect_af      = if (is.na(r$effect_af))   '' else r$effect_af,
          hwe_p          = if (is.na(r$hwe_p))       '' else r$hwe_p,
          qc_excl_reason = qc_flt_reason
        ))
      }
    },

    .fillCoverageTable = function(snpCols, wtable, missing_st, valid_snps,
                                   valid_counts, n_total, has_file = FALSE,
                                   n_qc_excl = 0L, n_flt_excl = 0L) {
      inData    <- names(self$data)
      matched   <- intersect(wtable$rsid[wtable$matched], inData)
      n_weights <- sum(wtable$matched)
      n_indata  <- length(intersect(snpCols, inData))
      n_matched <- length(matched)
      pct       <- if (n_weights > 0) round(n_matched / n_weights * 100, 1) else 0
      ambiguous <- sum(private$.isAmbiguous(wtable$effect_allele, wtable$other_allele))
      flipped   <- sum(wtable$strand_flipped == TRUE, na.rm = TRUE)
      mismatch  <- sum(grepl("mismatch", wtable$allele_status, ignore.case = TRUE))
      null_fix  <- sum(grepl("coded set to NA", wtable$allele_status, fixed = TRUE))
      # Complete cases: individuals with observed genotypes for ALL kept SNPs,
      # computed from valid_counts (pre-imputation observation flags) before
      # any row-exclusion.  valid_counts may already be row-filtered when
      # missing_st == "exclude"; n_total is the original full-data row count.
      complete  <- if (!is.null(valid_counts) && ncol(valid_counts) > 0)
                     sum(rowSums(!valid_counts) == 0)
                   else n_total


      meta <- attr(wtable, "pgs_meta") %||%
              list(pgs_id = "", pgs_name = "", trait_reported = "",
                   weight_type = "", genome_build = "")

      tbl <- self$results$coverageTable
      tbl$deleteRows()

      # Helper: only add row if value is non-empty
      add <- function(field, value) {
        v <- as.character(value)
        if (nchar(trimws(v)) > 0)
          tbl$addRow(rowKey = field, values = list(field = field, value = v))
      }

      # Score metadata (from file header)
      add("PGS ID",          meta$pgs_id         %||% "")
      add("Score name",      meta$pgs_name       %||% "")
      add("Trait",           meta$trait_reported %||% "")
      add("Weight type",     meta$weight_type    %||% "")
      add("Genome build",    meta$genome_build   %||% "")

      # SNP counts
      add("SNPs in weights file",       n_weights)
      add("SNPs in dataset",            n_indata)
      add("SNPs matched",               paste0(n_matched, " (", pct, "%)"))
      add("Ambiguous SNPs (AT/CG)",            ambiguous)
      add("Strand flipped (corrected)",        flipped)
      add("Allele mismatch (excluded)",        mismatch)
      if (null_fix > 0)
        add("Null-allele genotypes fixed (0/0 \u2192 NA)", null_fix)
      add("Excluded by allele/monomorphic QC", n_qc_excl)
      if (n_flt_excl > 0) {
        add("Excluded by missingness filter",  sum(grepl("\u274c excl (missing", wtable$allele_status, fixed = TRUE)))
        add("Excluded by HWE filter",          sum(grepl("\u274c excl (HWE", wtable$allele_status, fixed = TRUE)))
      }
      add("SNPs used in score",                length(valid_snps))
      add("Missing genotype strategy",         missing_st)
      run_w  <- self$options$weightingMode %in% c("weighted", "both") && has_file
      run_uw <- self$options$weightingMode %in% c("unweighted", "both") || !has_file
      if (run_w)  add("Score scale (Weighted)",   private$.scaleLabel("Weighted"))
      if (run_uw) add("Score scale (Unweighted)", private$.scaleLabel("Unweighted"))
      add("Total sample size",   n_total)
      add("Complete cases (no missing SNPs)", paste0(complete, " (", round(complete / n_total * 100, 1), "%)"))
    },

    .isAmbiguous = function(ea, oa) {
      pairs <- paste0(toupper(ea), toupper(oa))
      pairs %in% c("AT", "TA", "CG", "GC")
    },

    .fillSummaryTable = function(scores, resp = NULL, score_type = "") {
      tbl <- self$results$summaryTable

      add_row <- function(grp_label, sc) {
        sc <- sc[!is.na(sc)]
        n  <- length(sc)
        if (n == 0) return()
        mu  <- mean(sc)
        s   <- sd(sc)
        ci  <- if (n >= 2) {
          err <- qt(0.975, df = n - 1) * s / sqrt(n)
          c(mu - err, mu + err)
        } else c(NA_real_, NA_real_)
        row_key <- paste0(score_type, "_", grp_label)
        tbl$addRow(rowKey = row_key, values = list(
          score_type = score_type,
          group      = grp_label,
          n          = n,
          mean       = mu,
          sd         = s,
          ci_low     = ci[1],
          ci_high    = ci[2],
          min        = min(sc),
          max        = max(sc),
          skew       = skewness(sc)
        ))
      }

      is_binary <- !is.null(resp) && (is.factor(resp) ||
                    length(unique(resp[!is.na(resp)])) == 2)

      if (is_binary) {
        lvls <- levels(droplevels(factor(resp[!is.na(resp)])))
        for (lv in lvls)
          add_row(as.character(lv),
                  scores[!is.na(resp) & as.character(resp) == lv])
        add_row("Overall", scores)
      } else {
        add_row("Overall", scores)
      }
    },

    .saveScoresToData = function(all_scores) {
      out   <- self$results$saveScores
      n_tot <- nrow(self$data)
      keep  <- private$.keepMask   # NULL or logical[n_tot]

      keys   <- character(0)
      titles <- character(0)
      descs  <- character(0)
      types  <- character(0)

      for (mode_label in names(all_scores)) {
        key   <- paste0("PGS_", mode_label)
        title <- paste0("PGS (", mode_label, ")")
        desc  <- paste0("Polygenic score — ", tolower(mode_label), " weighting")
        keys   <- c(keys,   key)
        titles <- c(titles, title)
        descs  <- c(descs,  desc)
        types  <- c(types,  "continuous")
      }

      out$set(keys = keys, titles = titles,
              descriptions = descs, measureTypes = types)

      # Row numbers: map scores back to original data rows via keepMask
      all_row_nums <- seq_len(n_tot)
      score_rows   <- if (!is.null(keep)) all_row_nums[keep] else all_row_nums
      out$setRowNums(as.character(score_rows))

      for (mode_label in names(all_scores)) {
        key    <- paste0("PGS_", mode_label)
        scores <- all_scores[[mode_label]]
        # setValues receives the score-row values directly; setRowNums already
        # tells jamovi which original row numbers these correspond to.
        out$setValues(key = key, values = scores)
      }
    },

    .fillPercentileTable = function(all_scores, resp = NULL, respCol = NULL,
                                    covs = NULL) {

      # ── Parse percentile breaks ────────────────────────────────────────────
      breaks_str <- trimws(self$options$percentileBreaks)
      breaks_num <- suppressWarnings(as.numeric(strsplit(breaks_str, ",")[[1]]))
      breaks_num <- sort(unique(breaks_num[!is.na(breaks_num) &
                                           breaks_num > 0 & breaks_num < 100]))
      if (length(breaks_num) == 0) breaks_num <- c(20, 40, 60, 80, 90, 95)

      ref_opt <- self$options$pgsRefCategory   # "lowest" | "highest" | "middle"

      # ── Determine response type (mirrors .fillAssocTable logic) ───────────
      resp_type <- "none"
      if (!is.null(resp)) {
        n_lvls <- length(unique(resp[!is.na(resp)]))
        resp_type <- if (!is.factor(resp) && n_lvls > 5) "continuous"
                     else if (n_lvls == 2)               "binary"
                     else if (n_lvls > 2)                "polytomous"
                     else                                "continuous"
      }

      has_covs    <- !is.null(covs) && ncol(covs) > 0
      has_resp    <- resp_type %in% c("binary", "continuous")
      do_strat    <- resp_type %in% c("binary", "polytomous") && !is.null(resp)
      resp_levels <- if (do_strat) levels(droplevels(factor(resp[!is.na(resp)]))) else character(0)

      # ── Build category labels from breaks ─────────────────────────────────
      make_labels <- function(brks) {
        n <- length(brks)
        if (n == 0) return(">P0")
        lbl <- character(n + 1)
        lbl[1] <- paste0("<P", brks[1])
        for (i in seq_len(n - 1))
          lbl[i + 1] <- paste0("P", brks[i], "\u2013P", brks[i + 1])
        lbl[n + 1] <- paste0(">P", brks[n])
        lbl
      }

      cat_labels <- make_labels(breaks_num)
      n_cats     <- length(cat_labels)

      # ── Pre-compute cuts and cat_idx per mode (shared between both tables) ─
      mode_data <- lapply(all_scores, function(scores) {
        cuts    <- quantile(scores, breaks_num / 100, na.rm = TRUE)
        cat_idx <- findInterval(scores, cuts) + 1L
        cat_idx <- pmin(pmax(cat_idx, 1L), n_cats)
        list(scores = scores, cuts = cuts, cat_idx = cat_idx)
      })

      # ── Reference index ────────────────────────────────────────────────────
      ref_idx <- switch(ref_opt,
        lowest  = 1L,
        highest = n_cats,
        middle  = {
          s0  <- all_scores[[1]]
          med <- median(s0, na.rm = TRUE)
          idx <- findInterval(med, mode_data[[1]]$cuts) + 1L
          as.integer(pmin(pmax(idx, 1L), n_cats))
        }
      )

      # ══════════════════════════════════════════════════════════════════════
      # TABLE 1: Counts  (percentileThreshTable) — shown first
      # ══════════════════════════════════════════════════════════════════════
      thr_tbl <- self$results$percentileThreshTable
      thr_tbl$deleteRows()

      # Dynamic per-level columns (added once, before any rows)
      if (do_strat) {
        for (lv in resp_levels)
          thr_tbl$addColumn(
            name  = paste0("n_lv_", make.names(lv)),
            title = as.character(lv),
            type  = "text")
      }

      show_mode_col <- length(all_scores) > 1

      for (mode_label in names(mode_data)) {
        md      <- mode_data[[mode_label]]
        scores  <- md$scores
        cat_idx <- md$cat_idx
        n_total <- sum(!is.na(scores))
        resp_chr <- if (do_strat) as.character(resp) else NULL

        for (ci in seq_len(n_cats)) {
          lbl    <- cat_labels[ci]
          mask   <- cat_idx == ci & !is.na(scores)
          n_i    <- sum(mask)
          sc_i   <- scores[mask]
          rng    <- if (n_i > 0)
            sprintf("[%.3f, %.3f]", min(sc_i), max(sc_i)) else ""

          overall_str <- sprintf("%d (%.1f%%)",
                                 n_i,
                                 if (n_total > 0) 100 * n_i / n_total else 0)

          row_vals <- list(
            score_type = if (show_mode_col) mode_label else "",
            category   = lbl,
            score_range = rng,
            n_overall  = overall_str
          )

          # Per-level counts (binary / polytomous)
          if (do_strat) {
            for (lv in resp_levels) {
              n_lv_total <- sum(resp_chr == lv, na.rm = TRUE)
              n_lv_cat   <- sum(mask & resp_chr == lv, na.rm = TRUE)
              row_vals[[paste0("n_lv_", make.names(lv))]] <-
                sprintf("%d (%.1f%%)",
                        n_lv_cat,
                        if (n_lv_total > 0) 100 * n_lv_cat / n_lv_total else 0)
            }
          }

          thr_tbl$addRow(rowKey = paste0(mode_label, "_cat", ci),
                         values = row_vals)
        }
      }

      # ══════════════════════════════════════════════════════════════════════
      # TABLE 2: Regression  (percentileTable)
      # ══════════════════════════════════════════════════════════════════════
      cat_tbl <- self$results$percentileTable
      cat_tbl$deleteRows()

      # has_resp now includes polytomous — all three response types get a model
      has_resp <- resp_type %in% c("binary", "continuous", "polytomous")

      # Column title and notes
      ref_lbl <- cat_labels[ref_idx]

      if (!has_resp) {
        cat_tbl$getColumn("estimate")$setTitle("Estimate")
        cat_tbl$setNote("modelNote",
          "No response variable selected \u2014 counts and score ranges shown.")
        cat_tbl$setNote("refNote", NULL)
        cat_tbl$setNote("covNote", NULL)
      } else if (resp_type == "binary") {
        cat_tbl$getColumn("estimate")$setTitle("OR")
        cat_tbl$setNote("modelNote", "Logistic regression (OR, 95% CI)")
        cat_tbl$setNote("refNote",
          paste0("Reference category: ", ref_lbl, " (OR = 1)"))
        cat_tbl$setNote("covNote",
          if (has_covs) paste0("Adjusted for: ", paste(names(covs), collapse = ", "))
          else NULL)
      } else if (resp_type == "continuous") {
        cat_tbl$getColumn("estimate")$setTitle("\u03b2")
        cat_tbl$setNote("modelNote", "Linear regression (\u03b2, 95% CI)")
        cat_tbl$setNote("refNote",
          paste0("Reference category: ", ref_lbl, " (\u03b2 = 0)"))
        cat_tbl$setNote("covNote",
          if (has_covs) paste0("Adjusted for: ", paste(names(covs), collapse = ", "))
          else NULL)
      } else {
        # polytomous
        cat_tbl$getColumn("estimate")$setTitle("OR")
        cat_tbl$setNote("modelNote",
          paste0("Polytomous logistic regression (nnet::multinom); ",
                 "ORs relative to outcome reference level: ",
                 resp_levels[1]))
        cat_tbl$setNote("refNote",
          paste0("Reference score category: ", ref_lbl, " (OR = 1 per contrast)"))
        cat_tbl$setNote("covNote",
          if (has_covs) paste0("Adjusted for: ", paste(names(covs), collapse = ", "))
          else NULL)
      }

      for (mode_label in names(mode_data)) {
        md      <- mode_data[[mode_label]]
        scores  <- md$scores
        cat_idx <- md$cat_idx

        # ── Pre-compute per-category display values (shared across contrasts) ──
        cat_display <- lapply(seq_len(n_cats), function(ci) {
          mask <- cat_idx == ci & !is.na(scores)
          sc_i <- scores[mask]
          list(
            n   = sum(mask),
            rng = if (sum(mask) > 0)
                    sprintf("[% .3f, % .3f]", min(sc_i), max(sc_i))
                  else ""
          )
        })

        # ── Build the analysis data frame (common to all response types) ──────
        df <- if (has_covs)
          data.frame(cat = cat_idx, resp = resp, covs, check.names = FALSE)
        else
          data.frame(cat = cat_idx, resp = resp)
        df <- df[complete.cases(df), ]

        df$cat <- relevel(factor(df$cat, levels = seq_len(n_cats),
                                 labels = cat_labels),
                          ref = cat_labels[ref_idx])

        cov_terms <- if (has_covs)
          paste(paste0("`", names(covs), "`"), collapse = " + ") else ""

        # ── Fit model and emit rows ──────────────────────────────────────────

        # helper: emit one block of n_cats rows for a given contrast label,
        # coefficient matrix row name, and value extractor function.
        emit_rows <- function(contrast_lbl, get_est_ci_p) {
          for (ci in seq_len(n_cats)) {
            lbl    <- cat_labels[ci]
            is_ref <- (ci == ref_idx)
            disp   <- cat_display[[ci]]

            if (is_ref) {
              est <- if (resp_type == "continuous") 0 else 1
              ci_lo <- ''; ci_hi <- ''; p_v <- ''
            } else {
              res <- get_est_ci_p(lbl)   # returns list(est, ci_lo, ci_hi, p)
              est   <- res$est
              ci_lo <- res$ci_lo
              ci_hi <- res$ci_hi
              p_v   <- res$p
            }

            cat_tbl$addRow(
              rowKey = paste0(mode_label, "_", make.names(contrast_lbl), "_cat", ci),
              values = list(
                score_type  = mode_label,
                contrast    = contrast_lbl,
                category    = if (is_ref) paste0(lbl, " \u25c6") else lbl,
                n           = disp$n,
                score_range = disp$rng,
                estimate    = est,
                ci_low      = ci_lo,
                ci_high     = ci_hi,
                p           = p_v
              ))
          }
        }

        # ── Binary ──────────────────────────────────────────────────────────
        if (resp_type == "binary") {
          df$resp <- factor(df$resp)
          frm <- if (has_covs) as.formula(paste("resp ~ cat +", cov_terms))
                 else resp ~ cat
          fit <- tryCatch(glm(frm, data = df, family = binomial()),
                          error = function(e) NULL)

          cf  <- if (!is.null(fit)) coef(summary(fit))          else NULL
          cis <- if (!is.null(fit)) tryCatch(confint.default(fit),
                                             error = function(e) NULL) else NULL

          get_binary <- function(lbl) {
            coef_nm <- paste0("cat", lbl)
            if (is.null(cf) || !coef_nm %in% rownames(cf))
              return(list(est = '', ci_lo = '', ci_hi = '', p = ''))
            b  <- cf[coef_nm, 1]; se <- cf[coef_nm, 2]; p_v <- cf[coef_nm, 4]
            ci <- if (!is.null(cis) && coef_nm %in% rownames(cis))
                    cis[coef_nm, ] else c(b - 1.96*se, b + 1.96*se)
            list(est = .exp_or(b), ci_lo = .exp_or(ci[1]), ci_hi = .exp_or(ci[2]), p = p_v)
          }

          emit_rows("", get_binary)

        # ── Continuous ───────────────────────────────────────────────────────
        } else if (resp_type == "continuous") {
          df$resp <- as.numeric(df$resp)
          frm <- if (has_covs) as.formula(paste("resp ~ cat +", cov_terms))
                 else resp ~ cat
          fit <- tryCatch(lm(frm, data = df), error = function(e) NULL)

          cf  <- if (!is.null(fit)) coef(summary(fit))     else NULL
          cis <- if (!is.null(fit)) tryCatch(confint(fit),
                                             error = function(e) NULL) else NULL

          get_linear <- function(lbl) {
            coef_nm <- paste0("cat", lbl)
            if (is.null(cf) || !coef_nm %in% rownames(cf))
              return(list(est = '', ci_lo = '', ci_hi = '', p = ''))
            b  <- cf[coef_nm, 1]; se <- cf[coef_nm, 2]; p_v <- cf[coef_nm, 4]
            ci <- if (!is.null(cis) && coef_nm %in% rownames(cis))
                    cis[coef_nm, ] else c(b - 1.96*se, b + 1.96*se)
            list(est = b, ci_lo = ci[1], ci_hi = ci[2], p = p_v)
          }

          emit_rows("", get_linear)

        # ── Polytomous ───────────────────────────────────────────────────────
        } else if (resp_type == "polytomous") {
          df$resp <- factor(df$resp)
          frm <- if (has_covs) as.formula(paste("resp ~ cat +", cov_terms))
                 else resp ~ cat
          fit <- tryCatch(nnet::multinom(frm, data = df, trace = FALSE),
                          error = function(e) NULL)

          cf_mat <- if (!is.null(fit)) coef(fit) else NULL        # (K-1) × p
          vc_mat <- if (!is.null(fit))
                      tryCatch(vcov(fit), error = function(e) NULL)
                    else NULL

          for (lv in resp_levels[-1]) {
            local({
              lv_row       <- as.character(lv)
              contrast_lbl <- paste0(lv_row, " vs ", resp_levels[1])

              get_poly <- function(lbl) {
                coef_nm <- paste0("cat", lbl)
                if (is.null(cf_mat) || !lv_row %in% rownames(cf_mat) ||
                    !coef_nm %in% colnames(cf_mat))
                  return(list(est = '', ci_lo = '', ci_hi = '', p = ''))
                b     <- cf_mat[lv_row, coef_nm]
                vc_nm <- paste0(lv_row, ":", coef_nm)
                se    <- if (!is.null(vc_mat) && vc_nm %in% rownames(vc_mat))
                           sqrt(vc_mat[vc_nm, vc_nm]) else NA_real_
                p_v   <- if (!is.na(se) && se > 0) 2 * pnorm(-abs(b / se)) else NA_real_
                ci    <- if (!is.na(se)) c(b - 1.96*se, b + 1.96*se)
                         else            c(NA_real_, NA_real_)
                list(est = .exp_or(b), ci_lo = .exp_or(ci[1]), ci_hi = .exp_or(ci[2]), p = p_v)
              }

              emit_rows(contrast_lbl, get_poly)
            })
          }
        }

        # ── No response: counts-only rows ────────────────────────────────────
        if (!has_resp) {
          for (ci in seq_len(n_cats)) {
            lbl  <- cat_labels[ci]
            disp <- cat_display[[ci]]
            cat_tbl$addRow(
              rowKey = paste0(mode_label, "_cat", ci),
              values = list(
                score_type  = mode_label,
                contrast    = "",
                category    = lbl,
                n           = disp$n,
                score_range = disp$rng,
                estimate    = '',
                ci_low      = '',
                ci_high     = '',
                p           = ''
              ))
          }
        }
      }
    },

    .fillAssocTable = function(scores, resp, respCol, covs = NULL,
                              score_type = "") {
      tbl <- self$results$assocTable

      has_covs <- !is.null(covs) && ncol(covs) > 0
      df <- if (has_covs)
              data.frame(pgs = scores, resp = resp, covs, check.names = FALSE)
            else
              data.frame(pgs = scores, resp = resp)
      df <- df[complete.cases(df), ]
      if (nrow(df) < 3) return()

      lvls      <- levels(droplevels(factor(df$resp)))
      n_lvls    <- length(lvls)
      resp_type <- if (!is.factor(df$resp) && n_lvls > 5) "continuous"
                   else if (n_lvls == 2)                  "binary"
                   else if (n_lvls > 2)                   "polytomous"
                   else                                   "continuous"

      # ── Table notes ───────────────────────────────────────────────────────
      if (has_covs)
        self$results$assocTable$setNote(
          "covNote", paste0("Adjusted for: ", paste(names(covs), collapse = ", ")))
      else
        self$results$assocTable$setNote("covNote", NULL)

      resp_note <- switch(resp_type,
        binary     = paste0("Response: ", respCol,
                            " (", lvls[2], " vs ", lvls[1], ")"),
        polytomous = paste0("Response: ", respCol,
                            " (ref: ", lvls[1], "; ",
                            paste(lvls[-1], collapse = ", "), ")"),
        paste0("Response: ", respCol)
      )
      self$results$assocTable$setNote("respNote", resp_note)

      cov_terms <- if (has_covs)
        paste(paste0("`", names(covs), "`"), collapse = " + ") else ""

      add_row <- function(test, stat_label, estimate, se, ci_low, ci_high,
                          stat, df_val, p) {
        row_key <- paste0(score_type, "_", test)
        tbl$addRow(rowKey = row_key, values = list(
          score_type = score_type,
          test       = test,
          stat_label = stat_label,
          estimate   = estimate,
          se         = se,
          ci_low     = ci_low,
          ci_high    = ci_high,
          stat       = stat,
          df         = as.character(df_val),
          p          = p
        ))
      }

      # ── Binary response ───────────────────────────────────────────────────
      if (resp_type == "binary") {

        df$resp <- factor(df$resp)
        g1 <- df$pgs[df$resp == lvls[1]]
        g2 <- df$pgs[df$resp == lvls[2]]

        frm <- if (has_covs) as.formula(paste("resp ~ pgs +", cov_terms))
               else resp ~ pgs
        fit <- tryCatch(glm(frm, data = df, family = binomial()),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          cf <- coef(summary(fit))
          ci <- tryCatch(confint.default(fit)["pgs", ],
                         error = function(e) c('', ''))
          if ("pgs" %in% rownames(cf))
            add_row("Logistic regression", "OR",
                                        .exp_or(cf["pgs", 1]), '',
                    .exp_or(ci[1]), .exp_or(ci[2]),
                    cf["pgs", 3], "", cf["pgs", 4])
        }

        tt <- tryCatch(t.test(g2, g1), error = function(e) NULL)
        if (!is.null(tt))
          add_row("Welch t-test", "t",
                  diff(tt$estimate), '',
                  tt$conf.int[1], tt$conf.int[2],
                  tt$statistic, round(tt$parameter, 1), tt$p.value)

        mw <- tryCatch(wilcox.test(g2, g1, exact = FALSE, conf.int = TRUE),
                       error = function(e) NULL)
        if (!is.null(mw))
          add_row("Mann-Whitney U", "W",
                  mw$estimate, '',
                  mw$conf.int[1], mw$conf.int[2],
                  mw$statistic, "", mw$p.value)

      # ── Polytomous response (>2 groups) ───────────────────────────────────
      } else if (resp_type == "polytomous") {

        df$resp <- factor(df$resp)   # first level is reference by default

        # Polytomous logistic via nnet::multinom — one OR row per non-ref level
        fit <- tryCatch(
          nnet::multinom(
            if (has_covs) as.formula(paste("resp ~ pgs +", cov_terms))
            else resp ~ pgs,
            data = df, trace = FALSE),
          error = function(e) NULL)

        if (!is.null(fit)) {
          cf_mat <- coef(fit)                       # levels × predictors matrix
          # Wald SEs from the vcov matrix (rows ordered as coef)
          vc     <- tryCatch(vcov(fit), error = function(e) NULL)

          for (lv in lvls[-1]) {
            lv_row <- as.character(lv)
            if (!lv_row %in% rownames(cf_mat)) next
            b  <- cf_mat[lv_row, "pgs"]
            se_b <- if (!is.null(vc)) {
              nm <- paste0(lv_row, ":pgs")
              if (nm %in% rownames(vc)) sqrt(vc[nm, nm]) else NA_real_
            } else NA_real_
            z    <- if (!is.na(se_b) && se_b > 0) b / se_b else NA_real_
            p_z  <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA_real_
            ci_lo <- b - 1.96 * se_b
            ci_hi <- b + 1.96 * se_b
            lbl   <- paste0("Polytomous logistic (", lv_row, " vs ", lvls[1], ")")
            add_row(lbl, "OR",
                    .exp_or(b), NA_real_, .exp_or(ci_lo), .exp_or(ci_hi),
                    z, "", p_z)
          }

          # Overall likelihood-ratio test for pgs across all levels
          fit0 <- tryCatch(
            nnet::multinom(
              if (has_covs) as.formula(paste("resp ~", cov_terms))
              else resp ~ 1,
              data = df, trace = FALSE),
            error = function(e) NULL)
          if (!is.null(fit0)) {
            lr    <- 2 * (logLik(fit) - logLik(fit0))
            df_lr <- length(lvls) - 1
            p_lr  <- pchisq(as.numeric(lr), df = df_lr, lower.tail = FALSE)
            add_row("Polytomous logistic (overall)", "\u03c7\u00b2",
                    '', '', '', '',
                    as.numeric(lr), df_lr, p_lr)
          }
        }

        # One-way ANOVA (parametric)
        aov_fit <- tryCatch(aov(pgs ~ resp, data = df), error = function(e) NULL)
        if (!is.null(aov_fit)) {
          sm  <- summary(aov_fit)[[1]]
          f   <- sm["resp", "F value"]
          df1 <- sm["resp", "Df"]
          df2 <- sm["Residuals", "Df"]
          p_f <- sm["resp", "Pr(>F)"]
          add_row("ANOVA", "F",
                  '', '', '', '',
                  f, paste0(df1, ", ", df2), p_f)
        }

        # Kruskal-Wallis (non-parametric)
        kw <- tryCatch(kruskal.test(pgs ~ resp, data = df), error = function(e) NULL)
        if (!is.null(kw))
          add_row("Kruskal-Wallis", "\u03c7\u00b2",
                  '', '', '', '',
                  kw$statistic, round(kw$parameter, 0), kw$p.value)

      # ── Continuous response ───────────────────────────────────────────────
      } else {

        resp_num <- as.numeric(df$resp)
        df$resp  <- resp_num

        frm <- if (has_covs) as.formula(paste("resp ~ pgs +", cov_terms))
               else resp ~ pgs
        fit <- tryCatch(lm(frm, data = df), error = function(e) NULL)
        if (!is.null(fit)) {
          cf <- coef(summary(fit))
          ci <- tryCatch(confint(fit)["pgs", ],
                         error = function(e) c('', ''))
          if ("pgs" %in% rownames(cf))
            add_row("Linear regression", "\u03b2",
                    cf["pgs", 1], cf["pgs", 2],
                    ci[1], ci[2],
                    cf["pgs", 3], df.residual(fit), cf["pgs", 4])
        }

        pc <- tryCatch(cor.test(df$pgs, resp_num, method = "pearson"),
                       error = function(e) NULL)
        if (!is.null(pc))
          add_row("Pearson correlation", "r",
                  pc$estimate, '',
                  pc$conf.int[1], pc$conf.int[2],
                  pc$statistic, round(pc$parameter, 0), pc$p.value)

        sp <- tryCatch(cor.test(df$pgs, resp_num, method = "spearman", exact = FALSE),
                       error = function(e) NULL)
        if (!is.null(sp)) {
          r  <- as.numeric(sp$estimate)
          n  <- nrow(df)
          z  <- 0.5 * log((1 + r) / (1 - r))
          se <- 1 / sqrt(n - 3)
          ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
          add_row("Spearman correlation", "\u03c1",
                  r, '', ci[1], ci[2],
                  sp$statistic, "", sp$p.value)
        }
      }
    },

    .fillInteractionTable = function(scores, resp, respCol, covs, score_type = "") {
      tbl     <- self$results$interactionTable
      cov1_nm <- names(covs)[1]
      cov1    <- covs[[1]]

      # Build analysis data frame
      other_covs <- if (ncol(covs) > 1) covs[, -1, drop = FALSE] else NULL
      has_other  <- !is.null(other_covs) && ncol(other_covs) > 0

      df_cols <- list(pgs = scores, resp = resp, cov1 = cov1)
      other_col_names <- character(0)
      if (has_other) {
        for (j in seq_len(ncol(other_covs))) {
          cn <- paste0("ocov", j)
          df_cols[[cn]] <- other_covs[[j]]
          other_col_names <- c(other_col_names, cn)
        }
      }
      df <- as.data.frame(df_cols, check.names = FALSE)
      df <- df[complete.cases(df), ]
      if (nrow(df) < 5) return()

      lvls      <- levels(droplevels(factor(df$resp)))
      n_lvls    <- length(lvls)
      resp_type <- if (!is.factor(df$resp) && n_lvls > 5) "continuous"
                   else if (n_lvls == 2)                   "binary"
                   else                                    "polytomous"

      # Rename estimate column to OR or β depending on response type
      est_title <- switch(resp_type,
        binary     = "OR",
        continuous = "β",
        "Estimate"
      )
      self$results$interactionTable$getColumn("estimate")$setTitle(est_title)

      # Set table note describing the interaction being tested
      if (has_other)
        self$results$interactionTable$setNote("intNote",
                paste0("Adjusted for: ",
                paste(names(other_covs), collapse = ", ")))
      else
        self$results$interactionTable$setNote("intNote", NULL)

      other_terms <- if (length(other_col_names) > 0)
        paste(" +", paste(other_col_names, collapse = " + "))
      else ""

      # Formulas use safe internal names; cov1 is always the column name
      frm_int_str  <- paste0("resp ~ pgs * cov1", other_terms)
      frm_main_str <- paste0("resp ~ pgs + cov1", other_terms)

      add_row <- function(model_lbl, term_lbl, estimate, ci_low, ci_high, p) {
        row_key <- paste0(score_type, "_", model_lbl, "_", term_lbl)
        tbl$addRow(rowKey = row_key, values = list(
          score_type = score_type,
          model      = model_lbl,
          term       = term_lbl,
          estimate   = estimate,
          ci_low     = ci_low,
          ci_high    = ci_high,
          p          = p
        ))
      }

      # ── Polytomous: multinomial logistic regression ───────────────────────
      if (resp_type == "polytomous") {
        df$resp <- factor(df$resp)

        fit_int  <- tryCatch(
          nnet::multinom(as.formula(frm_int_str),  data = df, trace = FALSE),
          error = function(e) NULL)
        fit_main <- tryCatch(
          nnet::multinom(as.formula(frm_main_str), data = df, trace = FALSE),
          error = function(e) NULL)
        if (is.null(fit_int)) return()

        cf_mat <- coef(fit_int)          # (K-1) × p matrix
        vc_mat <- tryCatch(vcov(fit_int), error = function(e) NULL)

        se_from_vc <- function(lv_row, coef_nm) {
          nm <- paste0(lv_row, ":", coef_nm)
          if (!is.null(vc_mat) && nm %in% rownames(vc_mat))
            sqrt(vc_mat[nm, nm])
          else NA_real_
        }

        for (lv in lvls[-1]) {
          local({
            lv_row    <- as.character(lv)
            model_lbl <- paste0("Polytomous logistic (", lv_row,
                                " vs ", lvls[1], ")")

            if (!lv_row %in% rownames(cf_mat)) return()

            report_term <- function(coef_nm, display_nm) {
              if (!coef_nm %in% colnames(cf_mat)) return()
              b  <- cf_mat[lv_row, coef_nm]
              se <- se_from_vc(lv_row, coef_nm)
              p  <- if (!is.na(se) && se > 0) 2 * pnorm(-abs(b / se)) else NA_real_
              ci <- if (!is.na(se)) c(b - 1.96 * se, b + 1.96 * se)
                    else            c(NA_real_, NA_real_)
              add_row(model_lbl, display_nm,
                      .exp_or(b), .exp_or(ci[1]), .exp_or(ci[2]), p)
            }

            report_term("pgs", "PGS (main)")

            cov1_main_nms <- colnames(cf_mat)[startsWith(colnames(cf_mat), "cov1") &
                                              !startsWith(colnames(cf_mat), "pgs:")]
            for (cnm in cov1_main_nms) {
              lbl <- if (cnm == "cov1") paste0(cov1_nm, " (main)")
                     else paste0(cov1_nm, " (", sub("^cov1", "", cnm), ")")
              report_term(cnm, lbl)
            }

            cov1_int_nms <- colnames(cf_mat)[startsWith(colnames(cf_mat), "pgs:cov1")]
            for (cnm in cov1_int_nms) {
              lbl <- if (cnm == "pgs:cov1") paste0("PGS \u00d7 ", cov1_nm, " (int)")
                     else paste0("PGS \u00d7 ", cov1_nm,
                                 " (", sub("^pgs:cov1", "", cnm), ")")
              report_term(cnm, lbl)
            }
          })
        }

        # LRT for the interaction block — one row shared across all contrasts
        if (!is.null(fit_main)) {
          lrt <- tryCatch(anova(fit_main, fit_int, test = "Chisq"),
                          error = function(e) NULL)
          if (!is.null(lrt) && nrow(lrt) >= 2) {
            p_col <- grep("^Pr", colnames(lrt), value = TRUE)[1]
            p_l   <- if (!is.null(p_col)) as.numeric(lrt[2, p_col]) else NA_real_
            add_row("Polytomous logistic (int, LRT)", "LRT (interaction)",
                    '', '', '', p_l)
          }
        }

        return()
      }

      # ── Binary: logistic regression ───────────────────────────────────────
      if (resp_type == "binary") {
        df$resp <- factor(df$resp)

        frm_int  <- as.formula(frm_int_str)
        frm_main <- as.formula(frm_main_str)

        fit_int  <- tryCatch(glm(frm_int,  data = df, family = binomial()),
                             error = function(e) NULL)
        fit_main <- tryCatch(glm(frm_main, data = df, family = binomial()),
                             error = function(e) NULL)
        if (is.null(fit_int)) return()

        cf  <- coef(summary(fit_int))
        cis <- tryCatch(confint.default(fit_int), error = function(e) NULL)

        report_term <- function(coef_nm, display_nm) {
          if (!coef_nm %in% rownames(cf)) return()
          b  <- cf[coef_nm, 1]
          p  <- cf[coef_nm, 4]
          ci <- if (!is.null(cis) && coef_nm %in% rownames(cis))
                  cis[coef_nm, ] else c(NA_real_, NA_real_)
          add_row("Logistic (int)", display_nm,
                  .exp_or(b), .exp_or(ci[1]), .exp_or(ci[2]), p)
        }

        report_term("pgs", "PGS (main)")
        # For factor cov1, R appends the level label: "cov1Male", "cov11", etc.
        # Report one row per non-reference level.
        cov1_main_nms <- rownames(cf)[startsWith(rownames(cf), "cov1") &
                                      !startsWith(rownames(cf), "pgs:")]
        for (cnm in cov1_main_nms) {
          lbl <- if (cnm == "cov1") paste0(cov1_nm, " (main)")
                 else paste0(cov1_nm, " (", sub("^cov1", "", cnm), ")")
          report_term(cnm, lbl)
        }
        cov1_int_nms <- rownames(cf)[startsWith(rownames(cf), "pgs:cov1")]
        for (cnm in cov1_int_nms) {
          lbl <- if (cnm == "pgs:cov1") paste0("PGS × ", cov1_nm, " (int)")
                 else paste0("PGS × ", cov1_nm,
                             " (", sub("^pgs:cov1", "", cnm), ")")
          report_term(cnm, lbl)
        }

        # LRT: interaction model vs main-effects model
        if (!is.null(fit_main)) {
          lrt <- tryCatch(anova(fit_main, fit_int, test = "LRT"), error = function(e) NULL)
          if (!is.null(lrt) && nrow(lrt) >= 2) {
            chi2 <- lrt[2, "Deviance"]
            df_l <- lrt[2, "Df"]
            p_l  <- lrt[2, "Pr(>Chi)"]
            add_row("Logistic (int)", "LRT (interaction)",
                    '', '', '', p_l)
          }
        }

      # ── Continuous: linear regression ─────────────────────────────────────
      } else {
        df$resp <- as.numeric(df$resp)

        frm_int  <- as.formula(frm_int_str)
        frm_main <- as.formula(frm_main_str)

        fit_int  <- tryCatch(lm(frm_int,  data = df), error = function(e) NULL)
        fit_main <- tryCatch(lm(frm_main, data = df), error = function(e) NULL)
        if (is.null(fit_int)) return()

        cf  <- coef(summary(fit_int))
        cis <- tryCatch(confint(fit_int), error = function(e) NULL)

        report_term <- function(coef_nm, display_nm) {
          if (!coef_nm %in% rownames(cf)) return()
          b  <- cf[coef_nm, 1]
          p  <- cf[coef_nm, 4]
          ci <- if (!is.null(cis) && coef_nm %in% rownames(cis))
                  cis[coef_nm, ] else c('', '')
          add_row("Linear (int)", display_nm,
                  b, ci[1], ci[2], p)
        }

        report_term("pgs", "PGS (main)")
        cov1_main_nms <- rownames(cf)[startsWith(rownames(cf), "cov1") &
                                      !startsWith(rownames(cf), "pgs:")]
        for (cnm in cov1_main_nms) {
          lbl <- if (cnm == "cov1") paste0(cov1_nm, " (main)")
                 else paste0(cov1_nm, " (", sub("^cov1", "", cnm), ")")
          report_term(cnm, lbl)
        }
        cov1_int_nms <- rownames(cf)[startsWith(rownames(cf), "pgs:cov1")]
        for (cnm in cov1_int_nms) {
          lbl <- if (cnm == "pgs:cov1") paste0("PGS × ", cov1_nm, " (int)")
                 else paste0("PGS × ", cov1_nm,
                             " (", sub("^pgs:cov1", "", cnm), ")")
          report_term(cnm, lbl)
        }

        # F-test: interaction model vs main-effects model
        if (!is.null(fit_main)) {
          ftest <- tryCatch(anova(fit_main, fit_int), error = function(e) NULL)
          if (!is.null(ftest) && nrow(ftest) >= 2) {
            f_val <- ftest[2, "F"]
            df1   <- ftest[2, "Df"]
            df2   <- ftest[2, "Res.Df"]
            p_f   <- ftest[2, "Pr(>F)"]
            add_row("Linear (int)", "F-test (interaction)",
                    '', '', '', p_f)
          }
        }
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # .plotDist
    #
    # Distribution plot — one panel per scoring mode (side-by-side when both
    # Weighted and Unweighted are active).  Within each panel:
    #   - No response / continuous response: overall distribution
    #   - Binary response: one curve/bar-set per group, overlaid
    #
    # Plot type controlled by distPlotType option:
    #   "density"   → kernel density curve(s) with filled polygon(s)
    #   "histogram" → side-by-side or stacked bars (freq = FALSE)
    #   "both"      → histogram bars + density overlay
    # ════════════════════════════════════════════════════════════════════════
    .plotDist = function(image, ...) {

      all_scores <- private$.cache$all_scores
      resp       <- private$.cache$resp
      respCol    <- private$.cache$respCol

      if (is.null(all_scores) || length(all_scores) == 0) return(FALSE)

      # ── Options ────────────────────────────────────────────────────────────
      plot_type <- self$options$distPlotType   # "density" | "histogram" | "both"
      n_breaks  <- self$options$histBreaks     # integer, used when histogram active
      breaks    <- if (!is.null(n_breaks) && !is.na(n_breaks) && n_breaks >= 2)
                     as.integer(n_breaks) else "Sturges"

      is_binary <- !is.null(resp) && (is.factor(resp) ||
                    length(unique(resp[!is.na(resp)])) == 2)

      # Group palette (controls / cases / extra levels)
      grp_pal <- c("#2980B9", "#C0392B", "#27AE60", "#8E44AD")

      n_modes  <- length(all_scores)
      mode_nms <- names(all_scores)

      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))

      # Layout: one column per mode
      par(mfrow = c(1, n_modes),
          bg    = "white",
          oma   = c(0, 0, if (n_modes > 1) 2 else 0, 0))

      for (mi in seq_along(all_scores)) {
        mode_lbl <- mode_nms[mi]
        scores   <- all_scores[[mi]]

        # Build a joint non-missing mask over scores AND resp (when present)
        # so that group subsetting uses perfectly aligned indices.
        valid_mask <- !is.na(scores)
        if (!is.null(resp)) valid_mask <- valid_mask & !is.na(resp)
        scores_v <- scores[valid_mask]
        resp_v   <- if (!is.null(resp)) resp[valid_mask] else NULL

        if (length(scores_v) < 2) next

        # ── Build per-group score lists ────────────────────────────────────
        if (is_binary) {
          lvls    <- levels(factor(resp_v))
          n_grps  <- length(lvls)
          colours <- grp_pal[seq_len(min(n_grps, length(grp_pal)))]
          resp_ch <- as.character(resp_v)
          grp_sc  <- lapply(lvls, function(lv) scores_v[resp_ch == lv])
          names(grp_sc) <- lvls
        } else {
          lvls    <- "Overall"
          n_grps  <- 1L
          colours <- "#2980B9"
          grp_sc  <- list(Overall = scores_v)
        }

        # ── Compute y-axis limits across all groups ────────────────────────
        x_range <- range(scores_v, na.rm = TRUE)
        # Ensure a non-degenerate range (guards against all-identical scores)
        if (diff(x_range) < .Machine$double.eps)
          x_range <- x_range + c(-0.5, 0.5)
        x_pad   <- diff(x_range) * 0.04
        x_range <- c(x_range[1] - x_pad, x_range[2] + x_pad)

        y_max <- 0
        hist_list <- list()
        dens_list <- list()

        for (lv in lvls) {
          sc <- grp_sc[[lv]]
          if (length(sc) < 2) next
          if (plot_type %in% c("histogram", "both")) {
            h <- hist(sc, plot = FALSE, breaks = breaks)
            hist_list[[lv]] <- h
            y_max <- max(y_max, max(h$density, na.rm = TRUE))
          }
          if (plot_type %in% c("density", "both")) {
            d <- density(sc, na.rm = TRUE)
            dens_list[[lv]] <- d
            y_max <- max(y_max, max(d$y, na.rm = TRUE))
          }
        }
        y_max <- y_max * 1.20

        # ── Set up panel ───────────────────────────────────────────────────
        par(mar = c(4.5, 4.5, 3.5, 1.5))
        panel_title <- if (n_modes > 1) mode_lbl else "PGS Distribution"

        y_lab <- "Density"

        # Empty plot frame
        plot(NULL,
             xlim = x_range, ylim = c(0, y_max),
             xlab = "PGS Score", ylab = y_lab,
             main = panel_title, las = 1)

        # ── Draw histogram bars (side-by-side for multiple groups) ─────────
        if (plot_type %in% c("histogram", "both") && length(hist_list) > 0) {
          if (n_grps == 1) {
            # Single group: standard histogram bars
            h <- hist_list[[lvls[1]]]
            rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density,
                 col    = adjustcolor(colours[1], alpha.f = 0.55),
                 border = adjustcolor(colours[1], alpha.f = 0.80))
          } else {
            # Multiple groups: side-by-side bars within each common break set
            # Use a shared break sequence based on overall data
            all_br  <- hist(scores_v, plot = FALSE, breaks = breaks)$breaks
            bw_vec  <- diff(all_br)          # per-bin widths (may differ)
            n_bins  <- length(bw_vec)
            sub_bw  <- bw_vec / n_grps       # per-bin sub-bar width

            for (gi in seq_along(lvls)) {
              lv <- lvls[gi]
              sc <- grp_sc[[lv]]
              if (length(sc) < 1) next
              # Count into shared bins
              cnts  <- tabulate(findInterval(sc, all_br, rightmost.closed = TRUE),
                                nbins = n_bins)
              dens_vals <- cnts / (length(sc) * bw_vec)   # density per bin
              x_left  <- all_br[-length(all_br)] + (gi - 1) * sub_bw
              x_right <- x_left + sub_bw * 0.92   # small gap between sub-bars
              rect(x_left, 0, x_right, dens_vals,
                   col    = adjustcolor(colours[gi], alpha.f = 0.60),
                   border = adjustcolor(colours[gi], alpha.f = 0.85))
            }
          }
        }

        # ── Draw density curves ────────────────────────────────────────────
        if (plot_type %in% c("density", "both") && length(dens_list) > 0) {
          for (gi in seq_along(lvls)) {
            lv <- lvls[gi]
            d  <- dens_list[[lv]]
            if (is.null(d)) next
            # Filled polygon for density
            polygon(c(d$x, rev(d$x)),
                    c(d$y,  rep(0, length(d$y))),
                    col    = adjustcolor(colours[gi], alpha.f = 0.18),
                    border = NA)
            lines(d, col = colours[gi], lwd = 2.2)
          }
        }

        # ── Group mean lines ───────────────────────────────────────────────
        for (gi in seq_along(lvls)) {
          lv <- lvls[gi]
          sc <- grp_sc[[lv]]
          if (length(sc) < 1) next
          abline(v   = mean(sc, na.rm = TRUE),
                 col = colours[gi], lwd = 1.6, lty = 2)
        }

        # ── Legend ─────────────────────────────────────────────────────────
        leg_labels <- lvls
        if (!is.null(respCol) && is_binary)
          leg_labels <- paste0(respCol, "=", lvls)
        legend("topright",
               legend = leg_labels,
               col    = colours,
               lwd    = 2, lty = 1,
               bty    = "n", cex = 0.82)

        # ── p-value annotation for binary (Mann-Whitney, two groups only) ──
        if (is_binary && length(lvls) == 2) {
          g1  <- grp_sc[[lvls[1]]]
          g2  <- grp_sc[[lvls[2]]]
          mw  <- tryCatch(wilcox.test(g2, g1, exact = FALSE), error = function(e) NULL)
          if (!is.null(mw)) {
            p_fmt <- if (mw$p.value < 0.001) "p < 0.001"
                     else paste0("p = ", round(mw$p.value, 3))
            mtext(paste0("Mann-Whitney ", p_fmt), side = 3, line = 0.25,
                  cex = 0.80, col = "#555555")
          }
        }
      }

      # Overall title when both modes shown
      if (n_modes > 1)
        mtext("PGS Distribution", outer = TRUE, cex = 1.05, font = 2, line = 0.4)

      TRUE
    },

    # ════════════════════════════════════════════════════════════════════════
    # .plotStrat
    # Scatter plot (continuous response only).  Binary response is now
    # handled inside .plotDist with per-group curves/bars.
    # ════════════════════════════════════════════════════════════════════════
    .plotStrat = function(image, ...) {

      all_scores <- private$.cache$all_scores
      resp       <- private$.cache$resp
      respCol    <- private$.cache$respCol

      if (is.null(all_scores) || is.null(resp)) return(FALSE)

      # Only continuous response reaches here (binary is in .plotDist)
      is_binary <- is.factor(resp) || length(unique(resp[!is.na(resp)])) == 2
      if (is_binary) return(FALSE)

      n_modes  <- length(all_scores)
      mode_nms <- names(all_scores)

      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))

      par(mfrow = c(1, n_modes),
          bg    = "white",
          oma   = c(0, 0, if (n_modes > 1) 2 else 0, 0))

      for (mi in seq_along(all_scores)) {
        mode_lbl <- mode_nms[mi]
        scores   <- all_scores[[mi]]

        df <- data.frame(pgs = scores, resp = as.numeric(resp))
        df <- df[complete.cases(df), ]
        if (nrow(df) < 3) next

        par(mar = c(5, 4.5, 4, 1.5))
        panel_title <- if (n_modes > 1) mode_lbl else "PGS vs Response"
        y_lab       <- if (!is.null(respCol) && nchar(respCol) > 0) respCol else "Response"

        plot(df$pgs, df$resp,
             pch = 19, col = adjustcolor("#2980B9", alpha.f = 0.40), cex = 0.70,
             xlab = "PGS Score", ylab = y_lab,
             main = panel_title, las = 1)

        fit <- tryCatch(lm(resp ~ pgs, data = df), error = function(e) NULL)
        if (!is.null(fit)) {
          abline(fit, col = "#C0392B", lwd = 2)
          r2    <- summary(fit)$r.squared
          p_val <- coef(summary(fit))[2, 4]
          p_fmt <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
          mtext(paste0("R² = ", round(r2, 3), "  |  ", p_fmt),
                side = 3, line = 0.25, cex = 0.82, col = "#555555")
        }
      }

      if (n_modes > 1)
        mtext("PGS vs Response", outer = TRUE, cex = 1.05, font = 2, line = 0.4)

      TRUE
    }

  )  # end private
)

