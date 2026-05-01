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
      idCol       <- self$options$idCol
      respCol     <- self$options$responseCol
      missing_st  <- self$options$missingStrategy
      normalize   <- self$options$normalize
      standardize <- self$options$standardize
      wmode       <- self$options$weightingMode

      if (is.null(snpCols) || length(snpCols) == 0) return()

      self$results$validationMsg$setContent(
        "<span style='color:#888;font-size:0.85em;'>snpPGS v0.5 — allele QC active</span>")
      self$results$validationMsg$setVisible(TRUE)

      # ── Build weight table from file (or unit weights) ───────────────────
      wtable <- private$.buildWeightTable(snpCols)
      if (is.null(wtable)) return()

      # ── Dosage matrix + allele QC ────────────────────────────────────────
      qc     <- private$.buildDosageMatrix(snpCols, wtable, missing_st)
      if (is.null(qc)) return()
      dosage <- qc$mat
      wtable <- qc$wtable

      private$.fillSnpGridTable(wtable)

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
      show_type_col <- length(modes) > 1

      # ── Individual IDs ───────────────────────────────────────────────────
      if (!is.null(idCol) && nchar(idCol) > 0 && idCol %in% names(self$data)) {
        ids_full <- as.character(self$data[[idCol]])
        ids <- if (!is.null(private$.keepMask)) ids_full[private$.keepMask] else ids_full
      } else {
        ids <- as.character(seq_len(nrow(dosage)))
      }

      # ── Resolve response vector ──────────────────────────────────────────
      resp <- NULL
      if (!is.null(respCol) && nchar(respCol) > 0 && respCol %in% names(self$data)) {
        resp_full <- self$data[[respCol]]
        resp <- if (!is.null(private$.keepMask) &&
                    length(private$.keepMask) == length(resp_full))
                  resp_full[private$.keepMask]
                else resp_full
        if (is.numeric(resp) && length(unique(resp[!is.na(resp)])) == 2 &&
            all(resp[!is.na(resp)] %in% c(0, 1)))
          resp <- factor(resp)
      }

      private$.fillCoverageTable(snpCols, wtable, missing_st)

      # ── Clear tables before multi-mode fill ─────────────────────────────
      self$results$summaryTable$deleteRows()
      self$results$assocTable$deleteRows()
      self$results$summaryTable$getColumn("score_type")$setVisible(show_type_col)
      self$results$assocTable$getColumn("score_type")$setVisible(show_type_col)

      first_scores <- NULL

      for (mode_label in names(modes)) {
        wvec <- modes[[mode_label]]

        # Drop NA-weight SNPs (only applies to catalog weights)
        na_snps <- names(wvec)[is.na(wvec)]
        if (length(na_snps) > 0) {
          keep_cols <- names(wvec)[!is.na(wvec)]
          wvec  <- wvec[keep_cols]
          dos_m <- dosage[, keep_cols, drop = FALSE]
        } else {
          dos_m <- dosage
        }

        if (length(wvec) == 0) next

        scores <- private$.computeScores(dos_m, wvec, qc, normalize, standardize,
                                          missing_st, mode_label == "Unweighted")
        if (is.null(first_scores)) {
          first_scores <- scores
          private$.pgsScores <- scores
          private$.idLabels  <- ids
        }

        private$.fillSummaryTable(scores, resp, mode_label, show_type_col)

        if (self$options$showAssoc && !is.null(resp))
          private$.fillAssocTable(scores, resp, respCol, mode_label, show_type_col)
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

      if (self$options$showScoreTable || self$options$showPercentiles)
        private$.fillScoreTable(first_scores, ids, resp)

      if (self$options$showPercentiles)
        private$.fillPercentileTable(first_scores)
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
    .computeScores = function(dosage, wvec, qc, normalize, standardize, missing_st,
                              is_unweighted = FALSE) {
      scores_raw <- as.numeric(dosage %*% wvec)
      if (normalize || missing_st == "snp_wise") {
        if (is_unweighted) {
          # For unweighted scoring, divide by the number of SNPs in the score.
          # valid_counts is unreliable here because the unweighted dosage bypass
          # may leave mat as NA for unparseable genotypes, making rowSums = 0.
          n_valid <- rep(length(wvec), nrow(dosage))
        } else {
          vc_cols <- intersect(names(wvec), colnames(qc$valid_counts))
          if (length(vc_cols) > 0) {
            vc      <- qc$valid_counts[, vc_cols, drop = FALSE]
            n_valid <- rowSums(vc)
          } else {
            n_valid <- rep(length(wvec), nrow(dosage))
          }
        }
        n_valid[n_valid == 0] <- NA_real_
        scores <- scores_raw / n_valid
      } else {
        scores <- scores_raw
      }
      if (standardize) {
        s <- sd(scores, na.rm = TRUE)
        if (!is.na(s) && s > 0) scores <- scores / s
      }
      scores
    },

    .buildWeightTable = function(snpCols) {
      path  <- self$options$weightsPath
      wmode <- self$options$weightingMode

      has_file <- !is.null(path) && nchar(trimws(path)) > 0 && file.exists(path)

      if (has_file && wmode != "unweighted")
        return(private$.parseCatalogFile(path, snpCols))

      # No file or unweighted mode: unit weights populated from actual data
      private$.unitWeightTableFromData(snpCols)
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

      c_rsid   <- find_col("rsid", "variant_id", "snp", "snp_id", "marker_name")
      c_ea     <- find_col("effect_allele")
      c_oa     <- find_col("other_allele", "ref_allele", "non_effect_allele", "reference_allele")
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
        stringsAsFactors = FALSE
      )

      # Scope to selected SNPs; unmatched SNPs get placeholder rows
      result       <- catalog[catalog$rsid %in% snpCols, , drop = FALSE]
      missing_snps <- setdiff(snpCols, result$rsid)
      if (length(missing_snps) > 0) {
        extra_rows <- data.frame(
          rsid = missing_snps, effect_allele = "", other_allele = "",
          effect_weight = NA_real_, chr = "", pos = "",
          matched = FALSE, allele_status = "not in weights file",
          strand_flipped = FALSE, extra_cols = "",
          n_missing = NA_integer_, pct_missing = NA_real_,
          stringsAsFactors = FALSE
        )
        result <- rbind(result, extra_rows)
      }

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

        if (!is.numeric(col)) {
          # Normalise to uppercase character vector, respecting factor level order
          if (is.factor(col)) {
            lvls     <- levels(col)
            col_char <- as.character(col)
          } else {
            lvls     <- character(0)
            col_char <- as.character(col)
          }
          col_char <- toupper(trimws(col_char))
          col_char[col_char %in% c("", "NA")] <- NA_character_

          # Split every genotype into single-base alleles (handles / | concatenated)
          bases_list <- strsplit(col_char[!is.na(col_char)], "[^ACGT]|(?<=.)(?=.)",
                                 perl = TRUE)
          bases_list <- lapply(bases_list, function(b) b[b %in% c("A","C","G","T")])
          all_alleles <- sort(unique(unlist(bases_list)))

          if (length(all_alleles) == 0) {
            qc <- "no valid genotypes \u274c"
          } else if (length(all_alleles) > 2) {
            qc <- "non-biallelic \u26a0"
            oa <- all_alleles[1]; ea <- all_alleles[2]
          } else if (length(all_alleles) == 1) {
            qc <- "monomorphic \u26a0"
            oa <- all_alleles[1]; ea <- all_alleles[1]
          } else {
            # Biallelic: identify reference from first homozygous level or
            # the homozygous genotype appearing at the lowest factor level.
            ref_allele <- NA_character_

            # Try factor levels first: find first level whose genotype is homozygous
            if (length(lvls) > 0) {
              for (lv in lvls) {
                lv_up   <- toupper(trimws(lv))
                lv_base <- unique(strsplit(lv_up, "[^ACGT]|(?<=.)(?=.)",
                                           perl = TRUE)[[1]])
                lv_base <- lv_base[lv_base %in% c("A","C","G","T")]
                if (length(lv_base) == 1) {
                  ref_allele <- lv_base
                  break
                }
              }
            }

            # Fallback: homozygous genotype with highest frequency
            if (is.na(ref_allele)) {
              homozyg <- sapply(all_alleles, function(a) {
                sum(sapply(bases_list, function(b) length(b) >= 2 && all(b == a)),
                    na.rm = TRUE)
              })
              if (any(homozyg > 0))
                ref_allele <- all_alleles[which.max(homozyg)]
              else
                ref_allele <- all_alleles[1]
            }

            oa <- ref_allele
            ea <- setdiff(all_alleles, ref_allele)
            qc <- "unweighted"
          }
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
    # .buildDosageMatrix
    #
    # For each SNP column:
    #   • Numeric columns  — trusted as dosage (0/1/2); flag as "dosage (unverified)"
    #     when allele info is available (we cannot check orientation).
    #   • Character columns — observed alleles are extracted from the data,
    #     compared against the catalog, and one of these statuses is assigned:
    #       "ok"                  alleles match catalog, dosage counts EA
    #       "strand flip"         complement alleles match; dosage corrected
    #       "allele mismatch"     alleles do not match catalog; SNP excluded
    #       "ambiguous"           AT/CG SNP — no strand resolution possible;
    #                             dosage still attempted (counts EA as given)
    #       "no allele info"      catalog has no EA — dosage set to NA
    #
    # wtable is annotated in-place (allele_status, strand_flipped columns).
    # SNPs with "allele mismatch" status are dropped from the returned matrix.
    # ════════════════════════════════════════════════════════════════════════
    .buildDosageMatrix = function(snpCols, wtable, missing_st) {

      useCols <- intersect(wtable$rsid, names(self$data))

      if (length(useCols) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>None of the selected SNP columns are present in the dataset.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return(NULL)
      }

      mat      <- matrix(NA_real_, nrow = nrow(self$data), ncol = length(useCols),
                         dimnames = list(NULL, useCols))
      exclude  <- character(0)   # SNPs to drop due to allele mismatch

      for (snp in useCols) {
        col_raw <- self$data[[snp]]
        idx     <- which(wtable$rsid == snp)[1]   # row index in wtable
        ea_cat  <- wtable$effect_allele[idx]       # already uppercase from parse
        oa_cat  <- wtable$other_allele[idx]

        has_allele_info <- !is.na(ea_cat) && nchar(ea_cat) > 0
        is_ambiguous    <- private$.isAmbiguous(ea_cat, oa_cat)[1]

        # ── Numeric column (dosage already encoded) ──────────────────────
        if (is.numeric(col_raw)) {
          mat[, snp] <- as.numeric(col_raw)
          if (wtable$allele_status[idx] == "") {
            wtable$allele_status[idx] <- if (has_allele_info)
              "dosage (orientation unverified)" else "ok (numeric dosage)"
          }
          next
        }

        # Skip allele QC for unweighted SNPs that had no valid genotypes or
        # are numeric — allele_status already set by unitWeightTableFromData.
        if (wtable$allele_status[idx] %in%
            c("no valid genotypes \u274c", "monomorphic \u26a0",
              "non-biallelic \u26a0")) {
          mat[, snp] <- NA_real_
          exclude    <- c(exclude, snp)
          next
        }

        if (wtable$allele_status[idx] == "unweighted") {
          # effect_allele and other_allele are populated from the data;
          # count copies of effect_allele per individual.
          col_char <- toupper(trimws(as.character(col_raw)))
          col_char[col_char %in% c("", "NA")] <- NA_character_
          ea_use   <- wtable$effect_allele[idx]

          if (nchar(ea_use) == 0) {
            # Numeric dosage column — already handled above by is.numeric;
            # should not reach here, but guard anyway.
            mat[, snp] <- NA_real_
            exclude    <- c(exclude, snp)
          } else {
            mat[, snp] <- vapply(col_char, function(g) {
              if (is.na(g)) return(NA_real_)
              bases <- strsplit(g, "[^ACGT]|(?<=.)(?=.)", perl = TRUE)[[1]]
              bases <- bases[bases %in% c("A","C","G","T")]
              as.numeric(sum(bases == ea_use))
            }, numeric(1))
          }
          next
        }

        # ── Character/factor column — allele strings ─────────────────────
        col_char <- toupper(trimws(as.character(col_raw)))
        col_char[col_char == "" | col_char == "NA"] <- NA_character_

        # Extract unique observed alleles: split each genotype string and keep
        # only single ACGT characters — this handles separators like "/" or "|"
        # (e.g. "G/T", "G|T", "GT") as well as plain concatenated strings.
        obs_alleles <- unique(unlist(strsplit(col_char[!is.na(col_char)], "")))
        obs_alleles <- obs_alleles[obs_alleles %in% c("A","C","G","T")]

        if (!has_allele_info) {
          mat[, snp] <- NA_real_
          wtable$allele_status[idx] <- "no allele info in weights file"
          exclude <- c(exclude, snp)
          next
        }

        if (length(obs_alleles) == 0) {
          wtable$allele_status[idx] <- "no valid genotypes observed"
          exclude <- c(exclude, snp)
          next
        }

        ea_comp <- private$.complement(ea_cat)
        oa_comp <- private$.complement(oa_cat)

        matches_direct     <- all(obs_alleles %in% c(ea_cat, oa_cat))
        matches_complement <- all(obs_alleles %in% c(ea_comp, oa_comp))

        # Determine effect allele to count and whether strand flip needed.
        # Ambiguous (AT/CG) catalog pairs are checked the same way — if the
        # data alleles don't match either orientation they are a mismatch, not
        # merely ambiguous. We only apply the ambiguous warning when the data
        # alleles DO match, because in that case strand resolution is uncertain.
        ea_use <- ea_cat
        status <- ""

        if (matches_direct) {
          status <- if (is_ambiguous) "ambiguous (AT/CG) \u26a0" else "ok"
        } else if (matches_complement) {
          ea_use <- ea_comp
          wtable$strand_flipped[idx] <- TRUE
          status <- if (is_ambiguous)
            "ambiguous (AT/CG) \u26a0"          # can't trust flip for ambiguous
          else
            "strand flip \u2194 (corrected)"
        } else {
          wtable$allele_status[idx] <- "allele mismatch \u274c"
          exclude <- c(exclude, snp)
          next
        }

        mat[, snp] <- vapply(col_char, function(g) {
          if (is.na(g)) return(NA_real_)
          bases <- strsplit(g, "")[[1]]
          bases <- bases[bases %in% c("A","C","G","T")]
          sum(bases == ea_use)
        }, numeric(1))

        wtable$allele_status[idx] <- status
      }

      # Drop mismatched SNPs
      mat <- mat[, setdiff(colnames(mat), exclude), drop = FALSE]

      # Per-SNP missingness (before imputation)
      for (snp in colnames(mat)) {
        idx <- which(wtable$rsid == snp)[1]
        n_na <- sum(is.na(mat[, snp]))
        wtable$n_missing[idx]   <- n_na
        wtable$pct_missing[idx] <- if (nrow(mat) > 0) round(n_na / nrow(mat) * 100, 1) else NA_real_
      }

      # Snapshot valid cells before any imputation/exclusion
      valid_before_impute <- !is.na(mat)

      # Handle within-column missingness
      for (snp in colnames(mat)) {
        na_mask <- is.na(mat[, snp])
        if (any(na_mask)) {
          mat[, snp] <- switch(missing_st,
            mean    = { m <- mean(mat[, snp], na.rm = TRUE)
                        mat[na_mask, snp] <- if (is.nan(m)) 0 else m
                        mat[, snp] },
            zero    = { mat[na_mask, snp] <- 0; mat[, snp] },
            exclude  = mat[, snp],
            snp_wise = { mat[na_mask, snp] <- 0; mat[, snp] }
          )
        }
      }

      if (missing_st == "exclude") {
        keep <- complete.cases(mat)
        # Apply keep mask to both mat and the valid_counts snapshot so they stay aligned
        valid_before_impute <- valid_before_impute[keep, , drop = FALSE]
        mat  <- mat[keep, , drop = FALSE]
        private$.keepMask <- keep
      } else {
        private$.keepMask <- rep(TRUE, nrow(self$data))
      }

      list(mat = mat, wtable = wtable, valid_counts = valid_before_impute)
    },


    # ════════════════════════════════════════════════════════════════════════
    # Output helpers
    # ════════════════════════════════════════════════════════════════════════

    .fillSnpGridTable = function(wtable) {
      tbl    <- self$results$snpGridTable
      tbl$deleteRows()
      inData <- names(self$data)

      # Hide columns that are entirely empty (field not present in file)
      has_chr   <- any(nchar(wtable$chr)           > 0, na.rm = TRUE)
      has_pos   <- any(nchar(wtable$pos)           > 0, na.rm = TRUE)
      has_ea    <- any(nchar(wtable$effect_allele) > 0, na.rm = TRUE)
      has_oa    <- any(nchar(wtable$other_allele)  > 0, na.rm = TRUE)
      has_extra <- any(nchar(wtable$extra_cols)    > 0, na.rm = TRUE)
      has_miss  <- any(!is.na(wtable$n_missing))

      tbl$getColumn("chr")$setVisible(has_chr)
      tbl$getColumn("pos")$setVisible(has_pos)
      tbl$getColumn("effect_allele")$setVisible(has_ea)
      tbl$getColumn("other_allele")$setVisible(has_oa)
      tbl$getColumn("extra_cols")$setVisible(has_extra)
      tbl$getColumn("n_missing")$setVisible(has_miss)
      tbl$getColumn("pct_missing")$setVisible(has_miss)

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
          effect_weight = if (is.na(r$effect_weight)) NA_real_ else r$effect_weight,
          matched       = if (in_ds) "\u2713" else "\u2717",
          allele_status = status,
          extra_cols    = as.character(r$extra_cols),
          n_missing     = if (is.na(r$n_missing))   NA_integer_ else as.integer(r$n_missing),
          pct_missing   = if (is.na(r$pct_missing)) NA_real_    else r$pct_missing
        ))
      }
    },

    .fillCoverageTable = function(snpCols, wtable, missing_st) {
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
      add("Missing genotype strategy",  missing_st)
    },

    .isAmbiguous = function(ea, oa) {
      pairs <- paste0(toupper(ea), toupper(oa))
      pairs %in% c("AT", "TA", "CG", "GC")
    },

    .fillSummaryTable = function(scores, resp = NULL, score_type = "", show_type_col = FALSE) {
      tbl <- self$results$summaryTable

      skewness <- function(x) {
        x <- x[!is.na(x)]; n <- length(x)
        if (n < 3) return(NA_real_)
        mu <- mean(x); s <- sd(x)
        if (s == 0) return(NA_real_)
        tryCatch(e1071::skewness(x),
          error = function(e)
            sum(((x - mu) / s)^3) * n / ((n - 1) * (n - 2)))
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
        row_key <- if (show_type_col) paste0(score_type, "_", grp_label) else grp_label
        tbl$addRow(rowKey = row_key, values = list(
          score_type = if (show_type_col) score_type else "",
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

    .fillScoreTable = function(scores, ids, resp = NULL) {
      tbl <- self$results$scoreTable
      tbl$deleteRows()
      mu  <- mean(scores, na.rm = TRUE)
      sig <- sd(scores,   na.rm = TRUE)
      for (i in seq_along(scores)) {
        z   <- if (!is.na(sig) && sig > 0) (scores[i] - mu) / sig else NA_real_
        pct <- if (!is.na(scores[i]))
                 mean(scores <= scores[i], na.rm = TRUE) * 100
               else NA_real_
        tbl$addRow(rowKey = i, values = list(
          individual = ids[i],
          pgs        = scores[i],
          pgs_z      = z,
          percentile = pct
        ))
      }
    },

    .fillPercentileTable = function(scores) {
      breaks_str <- trimws(self$options$percentileBreaks)
      breaks_num <- suppressWarnings(as.numeric(strsplit(breaks_str, ",")[[1]]))
      breaks_num <- breaks_num[!is.na(breaks_num) & breaks_num >= 0 & breaks_num <= 100]
      if (length(breaks_num) == 0) breaks_num <- c(20, 40, 60, 80, 90, 95)

      tbl <- self$results$percentileTable
      tbl$deleteRows()
      for (b in breaks_num) {
        tbl$addRow(rowKey = b, values = list(
          threshold = paste0("P", b),
          score     = quantile(scores, b / 100, na.rm = TRUE)
        ))
      }
    },

    .fillAssocTable = function(scores, resp, respCol, score_type = "", show_type_col = FALSE) {
      tbl <- self$results$assocTable

      df <- data.frame(pgs = scores, resp = resp)
      df <- df[complete.cases(df), ]
      if (nrow(df) < 3) return()

      is_binary <- is.factor(df$resp) || length(unique(df$resp)) == 2

      add_row <- function(test, stat_label, estimate, se, ci_low, ci_high,
                          stat, df_val, p) {
        row_key <- if (show_type_col) paste0(score_type, "_", test) else test
        tbl$addRow(rowKey = row_key, values = list(
          score_type = if (show_type_col) score_type else "",
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

      if (is_binary) {
        lvls <- levels(factor(df$resp))
        g1   <- df$pgs[as.character(df$resp) == lvls[1]]
        g2   <- df$pgs[as.character(df$resp) == lvls[2]]
        lbl  <- paste0(lvls[2], " vs ", lvls[1])

        # ── Logistic regression (Wald 95% CI on log-OR scale) ────────────
        fit <- tryCatch(
          glm(resp ~ pgs, data = df, family = binomial()),
          error = function(e) NULL)
        if (!is.null(fit)) {
          cf  <- coef(summary(fit))
          ci  <- tryCatch(confint.default(fit)["pgs", ],
                          error = function(e) c(NA_real_, NA_real_))
          if (nrow(cf) >= 2)
            add_row(paste0("Logistic regression (", respCol, ")"),
                    "log-OR", cf[2,1], cf[2,2], ci[1], ci[2],
                    cf[2,3], "", cf[2,4])
        }

        # ── Welch t-test (95% CI of mean difference) ─────────────────────
        tt <- tryCatch(t.test(g2, g1), error = function(e) NULL)
        if (!is.null(tt))
          add_row(paste0("Welch t-test (", lbl, ")"),
                  "t", diff(tt$estimate), NA_real_,
                  tt$conf.int[1], tt$conf.int[2],
                  tt$statistic, round(tt$parameter, 1), tt$p.value)

        # ── Mann-Whitney U (Hodges-Lehmann 95% CI) ────────────────────────
        mw <- tryCatch(wilcox.test(g2, g1, exact = FALSE, conf.int = TRUE),
                       error = function(e) NULL)
        if (!is.null(mw))
          add_row(paste0("Mann-Whitney U (", lbl, ")"),
                  "W", mw$estimate, NA_real_,
                  mw$conf.int[1], mw$conf.int[2],
                  mw$statistic, "", mw$p.value)

      } else {
        resp_num <- as.numeric(df$resp)

        # ── Linear regression (95% CI of slope) ──────────────────────────
        fit <- tryCatch(lm(resp_num ~ pgs, data = df), error = function(e) NULL)
        if (!is.null(fit)) {
          cf  <- coef(summary(fit))
          ci  <- tryCatch(confint(fit)["pgs", ],
                          error = function(e) c(NA_real_, NA_real_))
          dff <- df.residual(fit)
          if (nrow(cf) >= 2)
            add_row(paste0("Linear regression (", respCol, ")"),
                    "\u03b2", cf[2,1], cf[2,2], ci[1], ci[2],
                    cf[2,3], dff, cf[2,4])
        }

        # ── Pearson correlation (Fisher-z 95% CI via cor.test) ────────────
        pc <- tryCatch(
          cor.test(df$pgs, resp_num, method = "pearson"),
          error = function(e) NULL)
        if (!is.null(pc))
          add_row("Pearson correlation",
                  "r", pc$estimate, NA_real_,
                  pc$conf.int[1], pc$conf.int[2],
                  pc$statistic, round(pc$parameter, 0), pc$p.value)

        # ── Spearman correlation (Fisher-z approximation for CI) ──────────
        sp <- tryCatch(
          cor.test(df$pgs, resp_num, method = "spearman", exact = FALSE),
          error = function(e) NULL)
        if (!is.null(sp)) {
          r  <- as.numeric(sp$estimate)
          n  <- nrow(df)
          # Fisher-z transform CI
          z  <- 0.5 * log((1 + r) / (1 - r))
          se <- 1 / sqrt(n - 3)
          ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
          add_row("Spearman correlation",
                  "\u03c1", r, NA_real_,
                  ci[1], ci[2],
                  sp$statistic, "", sp$p.value)
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
