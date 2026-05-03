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

snpPGSClass <- R6::R6Class(
  "snpPGSClass",
  inherit = snpPGSBase,

  private = list(

    .pgsScores = NULL,
    .idLabels  = NULL,
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
      normalize   <- self$options$normalize
      standardize <- self$options$standardize
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

      private$.fillCoverageTable(snpCols, wtable, missing_st, valid_snps)

      # ── Clear tables before multi-mode fill ─────────────────────────────
      self$results$summaryTable$deleteRows()
      self$results$assocTable$deleteRows()
      self$results$interactionTable$deleteRows()

      first_scores <- NULL
      all_scores   <- list()   # named by mode_label, for saveScores output

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

        # normalize_eff strictly follows the user's checkbox.
        # For SNP-wise missing strategy, zero-fill already happened in buildDosageMatrix,
        # so the raw sum is over observed SNPs only. When the user enables normalization,
        # valid_counts correctly provides per-individual observed-SNP denominators.
        normalize_eff <- normalize

        # Fix: pass a valid_counts slice restricted to the SNPs actually in dos_m.
        # This matters when weighted mode drops NA-weight SNPs (dos_m has fewer
        # columns than qc$valid_counts), ensuring the normalization denominator
        # is computed over exactly the same SNP set as the raw score.
        qc_mode <- qc
        qc_mode$valid_counts <- qc$valid_counts[,
          intersect(colnames(qc$valid_counts), colnames(dos_m)), drop = FALSE]

        scores <- private$.computeScores(dos_m, wvec, qc_mode, normalize_eff, standardize,
                                          mode_label == "Unweighted")
        all_scores[[mode_label]] <- scores

        if (is.null(first_scores)) {
          first_scores       <- scores
          private$.pgsScores <- scores
          private$.idLabels  <- as.character(seq_len(nrow(dosage)))
        }

        private$.fillSummaryTable(scores, resp, mode_label)

        if (self$options$showAssoc && !is.null(resp))
          private$.fillAssocTable(scores, resp, respCol, covs, mode_label)

        if (self$options$showInteraction && !is.null(resp) &&
            !is.null(covs) && ncol(covs) > 0)
          private$.fillInteractionTable(scores, resp, respCol, covs, mode_label)
      }

      if (is.null(first_scores)) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>No SNPs with valid weights — cannot compute PGS.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      show_strat <- self$options$showDistPlot && !is.null(resp)
      self$results$stratPlot$setVisible(show_strat)
      if (show_strat) private$.cache$resp <- resp

      if (self$options$showPercentiles)
        private$.fillPercentileTable(all_scores)

      if (self$results$saveScores$isNotFilled())
        private$.saveScoresToData(all_scores)
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
    # ════════════════════════════════════════════════════════════════════════
    # .computeScores
    #
    # Pipeline (in order):
    #   1. Raw score  = dosage %*% wvec  (weighted sum of dosage per individual)
    #   2. Normalize  = raw / max_possible  where max_possible is computed
    #                   per individual using only the SNPs they had observed
    #                   (pre-imputation valid_counts):
    #                     unweighted: 2 * n_valid_snps
    #                     weighted:   2 * sum(wvec for valid snps)
    #                   Result is in [0, 1].
    #   3. Standardize = score / SD(score across individuals)
    #                   SD is the empirical SD of the sample scores — this is
    #                   the standard approach (used by PLINK --score). A
    #                   theoretical SD from marginal Hardy-Weinberg colMeans
    #                   would ignore LD between SNPs and give population-expected
    #                   rather than observed SD, making it inappropriate here.
    #
    # is_unweighted: signals unit weights so we use SNP-count denominator.
    # ════════════════════════════════════════════════════════════════════════
    .buildWeightTable = function(snpCols) {
      path  <- self$options$weightsPath

      has_file <- !is.null(path) && nchar(trimws(path)) > 0 && file.exists(path)

      if (has_file)
        return(private$.parseCatalogFile(path, snpCols))
      
      private$.unitWeightTableFromData(snpCols)
    },



    .computeScores = function(dosage, wvec, qc, normalize, standardize,
                              is_unweighted = FALSE) {

      # ── Step 1: raw weighted sum ─────────────────────────────────────────
      scores <- as.numeric(dosage %*% wvec)

      # ── Step 2: normalize by maximum possible score per individual ────────
      if (normalize) {
        vc_cols <- intersect(names(wvec), colnames(qc$valid_counts))

        if (length(vc_cols) > 0) {
          vc <- qc$valid_counts[, vc_cols, drop = FALSE]   # logical matrix

          if (is_unweighted) {
            # Max unweighted = 2 alleles × number of observed SNPs
            max_score <- 2 * rowSums(vc)
          } else {
            # Max weighted = 2 × sum of positive weights for observed SNPs per individual.
            # Only positive-weight SNPs can increase the score (dosage in {0,1,2}), so
            # negative weights must be excluded from the maximum — otherwise the
            # denominator is under-estimated and scores can exceed 1.
            w_sub     <- wvec[vc_cols]
            w_pos     <- pmax(w_sub, 0)
            max_score <- 2 * as.numeric(vc %*% w_pos)
          }
        } else {
          # valid_counts unavailable — fall back to global maximum
          max_score <- if (is_unweighted) 2 * length(wvec)
                       else               2 * sum(pmax(wvec, 0), na.rm = TRUE)
          max_score <- rep(max_score, nrow(dosage))
        }

        max_score[max_score == 0] <- NA_real_
        scores <- scores / max_score
      }

      # ── Step 3: standardize to SD = 1 ────────────────────────────────────
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

      # Scope to selected SNPs; unmatched SNPs get placeholder rows
      # result       <- catalog[catalog$rsid %in% snpCols, , drop = FALSE]
      # missing_snps <- setdiff(snpCols, result$rsid)
      # if (length(missing_snps) > 0) {
      #   extra_rows <- data.frame(
      #     rsid = missing_snps, effect_allele = "", other_allele = "",
      #     effect_weight = NA_real_, chr = "", pos = "",
      #     matched = FALSE, allele_status = "not in weights file \u274c",
      #     strand_flipped = FALSE, extra_cols = "",
      #     n_missing = NA_integer_, pct_missing = NA_real_,
      #     stringsAsFactors = FALSE
      #   )
      #   result <- rbind(result, extra_rows)
      # }

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
          allele_status = "not in weights file \u274c",
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
        allele_status  = "no weights file",
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
        qc  <- "unweighted"

        # SNP data must be factor or dosage
        # numeric → no allele inference
        if (is.numeric(col)) {
          ea <- ""
          oa <- ""
          qc <- "numeric (dosage)"
        } else {
          # infer alleles (best effort, no QC decisions)
          col_char <- toupper(trimws(as.character(col)))
          col_char[col_char %in% c("", "NA")] <- NA_character_

          bases_list <- strsplit(col_char[!is.na(col_char)],
                                "[^ACGT]|(?<=.)(?=.)", perl = TRUE)
          bases_list <- lapply(bases_list, function(b) b[b %in% c("A","C","G","T")])
          alleles <- sort(unique(unlist(bases_list)))

          if (length(alleles) >= 1) {
            oa <- alleles[1]
            ea <- if (length(alleles) >= 2) alleles[2] else alleles[1]
          }

          qc <- "unweighted"
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
    # .buildDosageMatrix (Revised)
    # ════════════════════════════════════════════════════════════════════════

    .buildDosageMatrix = function(snpCols, wtable, missing_st) {

      useCols <- intersect(wtable$rsid, names(self$data))

      if (length(useCols) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>None of the selected SNP columns are present in the dataset.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return(NULL)
      }

      mat <- matrix(NA_real_,
                    nrow = nrow(self$data),
                    ncol = length(useCols),
                    dimnames = list(NULL, useCols))

      exclude <- character(0)
      snps_excluded <- setNames(rep(FALSE, length(useCols)), useCols)

      for (snp in useCols) {

        allele_ok <- NA

        col_raw <- self$data[[snp]]
        idx     <- which(wtable$rsid == snp)[1]

        ea_cat <- toupper(as.character(wtable$effect_allele[idx]))
        oa_cat <- toupper(as.character(wtable$other_allele[idx]))
        has_allele_info <- !is.na(ea_cat) && nchar(trimws(ea_cat)) > 0

        # ── All missing ─────────────────────────────────────────────────────
        if (all(is.na(col_raw))) {
          exclude <- c(exclude, snp)
          snps_excluded[snp] <- TRUE
          wtable$allele_status[idx] <- "all missing  \u274c"
          next
        }

        # ── Determine input type ─────────────────────────────────────────────
        is_numeric_like <- is.numeric(col_raw)
        if (!is_numeric_like && (is.factor(col_raw) || is.character(col_raw))) {
          tmp <- suppressWarnings(as.numeric(as.character(col_raw)))
          if (mean(!is.na(tmp)) > 0.5) is_numeric_like <- TRUE
        }

        # ── Numeric dosage SNPs (NO allele QC) ───────────────────────────────
        if (is_numeric_like) {

          col_num <- suppressWarnings(as.numeric(as.character(col_raw)))
          invalid <- !is.na(col_num) & (col_num < 0 | col_num > 2)

          if (any(invalid)) {
            col_num[invalid] <- NA_real_
          }

          mat[, snp] <- col_num

          obs_vals <- sort(unique(col_num[!is.na(col_num)]))
          obs_str  <- paste0("dosage(", paste(obs_vals, collapse = ","), ")")

          wtable$allele_status[idx] <- paste0("numeric dosage (no allele QC): ", obs_str)

          if (length(obs_vals) <= 1) {
            exclude <- c(exclude, snp)
            snps_excluded[snp] <- TRUE
            wtable$allele_status[idx] <-
              paste0("constant numeric dosage \u274c: ", obs_str)
          }

          next
        }

        # ── Genotype / character SNPs ────────────────────────────────────────
        col_char <- toupper(trimws(as.character(col_raw)))
        col_char[col_char %in% c("", "NA")] <- NA_character_

        bases_list <- strsplit(col_char[!is.na(col_char)],
                              "[^ACGT]|(?<=.)(?=.)", perl = TRUE)
        bases_list <- lapply(bases_list, function(b) b[b %in% c("A","C","G","T")])
        alleles <- sort(unique(unlist(bases_list)))

        if (length(alleles) == 0) {
          exclude <- c(exclude, snp)
          snps_excluded[snp] <- TRUE
          wtable$allele_status[idx] <- "no valid alleles observed \u274c "
          next
        }

        obs_str <- paste(alleles, collapse = "/")

        # ── Multiallelic ─────────────────────────────────────────────────────
        if (length(alleles) > 2) {
          exclude <- c(exclude, snp)
          snps_excluded[snp] <- TRUE
          wtable$allele_status[idx] <-
            paste0("multiallelic \u274c : ", obs_str)
          next
        }

        # ── No allele info ───────────────────────────────────────────────────
        if (!has_allele_info) {
          exclude <- c(exclude, snp)
          snps_excluded[snp] <- TRUE
          wtable$allele_status[idx] <-
            paste0("no allele info \u274c : ", obs_str)
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
          allele_ok <- FALSE
          exclude <- c(exclude, snp)
          snps_excluded[snp] <- TRUE
          wtable$allele_status[idx] <-
            paste0("allele mismatch \u274c obs: ", obs_str,
                  "; exp: ", ea_cat, "/", oa_cat)
          next
        }

        allele_ok <- TRUE

        # ── Strand handling ──────────────────────────────────────────────────
        if (matches_direct) {
          ea_use <- ea_cat
          status <- "\u2705 ok"
        } else {
          ea_use <- ea_comp
          status <- "\u2705 strand flip \u2757"
          wtable$strand_flipped[idx] <- TRUE
        }

        if (private$.isAmbiguous(ea_cat, oa_cat)) {
          status <- paste0(status, "; ambiguous AT/CG \u2757")
        }

        # ── Dosage computation ───────────────────────────────────────────────
        dos <- vapply(col_char, function(g) {
          if (is.na(g)) return(NA_real_)
          b <- strsplit(g, "[^ACGT]|(?<=.)(?=.)", perl = TRUE)[[1]]
          b <- b[b %in% c("A","C","G","T")]
          if (length(b) != 2) return(NA_real_)
          sum(b == ea_use)
        }, numeric(1))

        mat[, snp] <- dos

        base_status <- paste0(status, ": ", obs_str)

        # ── Monomorphism ONLY if allele_ok == TRUE ───────────────────────────
        if (isTRUE(allele_ok)) {
          if (length(unique(dos[!is.na(dos)])) <= 1) {
            exclude <- c(exclude, snp)
            snps_excluded[snp] <- TRUE
            wtable$allele_status[idx] <-
              paste0("monomorphic \u274c : ", obs_str)
          } else {
            wtable$allele_status[idx] <- base_status
          }
        }
      }

      # ── Post‑QC summary ────────────────────────────────────────────────────
      keep_snps <- setdiff(colnames(mat), exclude)

      if (length(keep_snps) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>No SNPs passed QC filters.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return(NULL)
      }

      valid_counts <- !is.na(mat[, keep_snps, drop = FALSE])

      for (snp in colnames(mat)) {
        idx <- which(wtable$rsid == snp)[1]
        n_na <- sum(is.na(mat[, snp]))
        wtable$n_missing[idx]   <- n_na
        wtable$pct_missing[idx] <- round(n_na / nrow(mat) * 100, 1)
      }

      # ── Allele frequency + HWE for all matched SNPs ─────────────────────────
      # Uses snp_af_hwe() from snp_helpers.R — same parse_genotype / HWE.exact
      # pipeline as snpDesc, ensuring consistency across modules.
      for (snp in intersect(wtable$rsid, names(self$data))) {
        idx <- which(wtable$rsid == snp)[1]
        ea  <- as.character(wtable$effect_allele[idx])
        res <- snp_af_hwe(self$data[[snp]], effect_allele = ea)
        wtable$effect_af[idx] <- res$effect_af
        wtable$hwe_p[idx]     <- res$hwe_p
      }

      for (snp in keep_snps) {
        na_mask <- is.na(mat[, snp])
        if (any(na_mask)) {
          mat[na_mask, snp] <- switch(missing_st,
            mean     = { m <- mean(mat[, snp], na.rm = TRUE); if (is.nan(m)) 0 else m },
            zero     = 0,
            "SNP-wise" = 0,
            exclude  = NA_real_
          )
        }
      }
      for (snp in exclude) {
        idx <- which(wtable$rsid == snp)[1]
        if (!is.na(idx)) {
          wtable$n_missing[idx] <- NA_integer_
          wtable$pct_missing[idx] <- NA_real_
        }
      }

      if (missing_st == "exclude") {
        keep_rows <- complete.cases(mat[, keep_snps, drop = FALSE])
        mat <- mat[keep_rows, , drop = FALSE]
        valid_counts <- valid_counts[keep_rows, , drop = FALSE]
        private$.keepMask <- keep_rows
      } else {
        private$.keepMask <- rep(TRUE, nrow(self$data))
      }

      mat <- mat[, keep_snps, drop = FALSE]

      list(mat = mat, wtable = wtable, valid_counts = valid_counts, valid_snps = keep_snps)
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

      tbl$getColumn("chr")$setVisible(has_chr)
      tbl$getColumn("pos")$setVisible(has_pos)
      tbl$getColumn("effect_allele")$setVisible(has_ea)
      tbl$getColumn("other_allele")$setVisible(has_oa)
      tbl$getColumn("extra_cols")$setVisible(has_extra)
      tbl$getColumn("n_missing")$setVisible(has_miss)
      tbl$getColumn("pct_missing")$setVisible(has_miss)
      tbl$getColumn("effect_af")$setVisible(has_af)
      tbl$getColumn("hwe_p")$setVisible(has_hwe)

      for (i in seq_len(nrow(wtable))) {
        r      <- wtable[i, ]
        status <- as.character(r$allele_status)
        in_ds  <- r$rsid %in% inData

        tbl$addRow(rowKey = i, values = list(
          rsid          = as.character(r$rsid),
          chr           = as.character(r$chr),
          pos           = as.character(r$pos),
          effect_allele = as.character(r$effect_allele),
          other_allele  = as.character(r$other_allele),
          effect_weight = if (is.na(r$effect_weight)) '' else r$effect_weight,
          matched       = if (in_ds) "\u2713" else "\u2717",
          allele_status = status,
          extra_cols    = as.character(r$extra_cols),
          n_missing     = if (is.na(r$n_missing))   '' else as.integer(r$n_missing),
          pct_missing   = if (is.na(r$pct_missing)) '' else r$pct_missing,
          effect_af     = if (is.na(r$effect_af))   '' else r$effect_af,
          hwe_p         = if (is.na(r$hwe_p))       '' else r$hwe_p
        ))
      }
    },

    .fillCoverageTable = function(snpCols, wtable, missing_st, valid_snps) {
      inData    <- names(self$data)
      matched   <- intersect(wtable$rsid[wtable$matched], inData)
      n_weights <- sum(wtable$matched)
      n_indata  <- length(intersect(snpCols, inData))
      n_matched <- length(matched)
      pct       <- if (n_weights > 0) round(n_matched / n_weights * 100, 1) else 0
      ambiguous <- sum(private$.isAmbiguous(wtable$effect_allele, wtable$other_allele))
      flipped   <- sum(wtable$strand_flipped == TRUE, na.rm = TRUE)
      mismatch  <- sum(grepl("mismatch", wtable$allele_status, ignore.case = TRUE))

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
      add("Ambiguous SNPs (AT/CG)",     ambiguous)
      add("Strand flipped (corrected)", flipped)
      add("Allele mismatch (excluded)", mismatch)
      add("SNPs used in score",         length(valid_snps))
      add("Missing genotype strategy",  missing_st)
    },

    .isAmbiguous = function(ea, oa) {
      pairs <- paste0(toupper(ea), toupper(oa))
      pairs %in% c("AT", "TA", "CG", "GC")
    },

    .fillSummaryTable = function(scores, resp = NULL, score_type = "") {
      tbl <- self$results$summaryTable

      skewness <- function(x) {
        x <- x[!is.na(x)]; n <- length(x)
        if (n < 3) return(NA_real_)
        x <- x - mean(x)
        y <- sqrt(n) * sum(x ^ 3) / (sum(x ^ 2) ^ (3/2))      
        y * ((1 - 1 / n)) ^ (3/2)
      }

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
        lvls <- levels(factor(resp[!is.na(resp)]))
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
        # Expand to full-data length: NAs for excluded rows
        full_vals <- rep(NA_real_, n_tot)
        full_vals[score_rows] <- scores
        out$setValues(key = key, values = full_vals[score_rows])
      }
    },

    .fillPercentileTable = function(all_scores) {
      breaks_str <- trimws(self$options$percentileBreaks)
      breaks_num <- suppressWarnings(as.numeric(strsplit(breaks_str, ",")[[1]]))
      breaks_num <- breaks_num[!is.na(breaks_num) & breaks_num >= 0 & breaks_num <= 100]
      if (length(breaks_num) == 0) breaks_num <- c(20, 40, 60, 80, 90, 95)

      tbl <- self$results$percentileTable
      tbl$deleteRows()

      # Add one score column per active mode
      show_label <- length(all_scores) > 1
      for (mode_label in names(all_scores)) {
        col_name  <- paste0("score_", mode_label)
        col_title <- if (show_label) paste0("PGS (", mode_label, ")") else "PGS Score"
        tbl$addColumn(name = col_name, title = col_title,
                      type = "number", format = "zto")
      }

      for (b in breaks_num) {
        vals <- list(threshold = paste0("P", b))
        for (mode_label in names(all_scores)) {
          col_name       <- paste0("score_", mode_label)
          vals[[col_name]] <- quantile(all_scores[[mode_label]], b / 100, na.rm = TRUE)
        }
        tbl$addRow(rowKey = b, values = vals)
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

      lvls      <- levels(factor(df$resp))
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
                    exp(cf["pgs", 1]), '',
                    exp(ci[1]), exp(ci[2]),
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
                    exp(b), NA_real_, exp(ci_lo), exp(ci_hi),
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

      lvls      <- levels(factor(df$resp))
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

      # ── Polytomous: not supported ─────────────────────────────────────────
      if (resp_type == "polytomous") {
        add_row("(not yet supported for polytomous response)", "",
                '', '', '', '')
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
                  exp(b), exp(ci[1]), exp(ci[2]), p)
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

    .plotDist = function(image, ...) {
      scores <- private$.pgsScores
      if (is.null(scores) || length(scores) < 2) return(FALSE)

      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))
      par(mar = c(4.5, 4.5, 3, 1.5), bg = "white")

      dens  <- density(scores, na.rm = TRUE)
      h     <- hist(scores, plot = FALSE, breaks = "Sturges")
      y_max <- max(max(h$density), max(dens$y)) * 1.15

      hist(scores, freq = FALSE, breaks = "Sturges",
           col = "#AED6F1", border = "#2980B9",
           main = "PGS Distribution (Overall)",
           xlab = "PGS Score", ylab = "Density",
           ylim = c(0, y_max), las = 1)
      lines(dens, col = "#C0392B", lwd = 2.5)
      abline(v = mean(scores, na.rm = TRUE), col = "#27AE60", lwd = 2, lty = 2)
      legend("topright",
             legend = c("Density", "Mean"),
             col    = c("#C0392B", "#27AE60"),
             lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
      TRUE
    },

    .plotStrat = function(image, ...) {
      scores <- private$.pgsScores
      resp   <- private$.cache$resp
      if (is.null(scores) || length(scores) < 2 || is.null(resp)) return(FALSE)

      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))

      is_binary <- is.factor(resp) || length(unique(resp[!is.na(resp)])) == 2

      df <- data.frame(pgs = scores, resp = resp)
      df <- df[complete.cases(df), ]
      if (nrow(df) < 3) return(FALSE)

      if (is_binary) {
        # ── Overlapping density curves + boxplot per group ─────────────────
        lvls   <- levels(factor(df$resp))
        grp_colours <- c("#2980B9", "#C0392B", "#27AE60", "#8E44AD")
        grp_colours <- grp_colours[seq_along(lvls)]

        # Compute densities per group
        dens_list <- lapply(lvls, function(lv)
          density(df$pgs[as.character(df$resp) == lv], na.rm = TRUE))

        y_max <- max(sapply(dens_list, function(d) max(d$y))) * 1.2

        par(mar = c(5, 4.5, 4, 1.5), bg = "white")
        plot(NULL, xlim = range(df$pgs), ylim = c(0, y_max),
             xlab = "PGS Score", ylab = "Density",
             main = "PGS Distribution by Group", las = 1)

        for (i in seq_along(lvls)) {
          polygon(c(dens_list[[i]]$x, rev(dens_list[[i]]$x)),
                  c(dens_list[[i]]$y, rep(0, length(dens_list[[i]]$y))),
                  col = adjustcolor(grp_colours[i], alpha.f = 0.25), border = NA)
          lines(dens_list[[i]], col = grp_colours[i], lwd = 2.5)
          abline(v = mean(df$pgs[as.character(df$resp) == lvls[i]], na.rm = TRUE),
                 col = grp_colours[i], lwd = 1.5, lty = 2)
        }

        legend("topright", legend = lvls,
               col = grp_colours, lwd = 2, bty = "n", cex = 0.9)

        # Add p-value annotation (Mann-Whitney)
        g1  <- df$pgs[as.character(df$resp) == lvls[1]]
        g2  <- df$pgs[as.character(df$resp) == lvls[2]]
        mw  <- tryCatch(wilcox.test(g2, g1, exact = FALSE), error = function(e) NULL)
        if (!is.null(mw)) {
          p_fmt <- if (mw$p.value < 0.001) "p < 0.001"
                   else paste0("p = ", round(mw$p.value, 3))
          mtext(paste0("Mann-Whitney ", p_fmt), side = 3, line = 0.3,
                cex = 0.85, col = "#555555")
        }

      } else {
        # ── Scatter plot with regression line ──────────────────────────────
        par(mar = c(5, 4.5, 4, 1.5), bg = "white")
        resp_num <- as.numeric(df$resp)
        plot(df$pgs, resp_num,
             pch = 19, col = adjustcolor("#2980B9", alpha.f = 0.4), cex = 0.7,
             xlab = "PGS Score", ylab = "Response",
             main = "PGS vs Response", las = 1)

        fit <- tryCatch(lm(resp_num ~ pgs, data = df), error = function(e) NULL)
        if (!is.null(fit)) {
          abline(fit, col = "#C0392B", lwd = 2)
          r2    <- summary(fit)$r.squared
          p_val <- coef(summary(fit))[2, 4]
          p_fmt <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
          mtext(paste0("R² = ", round(r2, 3), "  |  ", p_fmt),
                side = 3, line = 0.3, cex = 0.85, col = "#555555")
        }
      }
      TRUE
    }

  )  # end private
)

`%||%` <- function(a, b) if (!is.null(a)) a else b