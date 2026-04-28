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

    # ════════════════════════════════════════════════════════════════════════
    # .run
    # ════════════════════════════════════════════════════════════════════════
    .run = function() {

      snpCols    <- self$options$snpCols
      idCol      <- self$options$idCol
      respCol    <- self$options$responseCol
      missing_st <- self$options$missingStrategy
      normalize  <- self$options$normalize

      if (is.null(snpCols) || length(snpCols) == 0) return()

      # ── Build weight table from file (or unit weights) ───────────────────
      wtable <- private$.buildWeightTable(snpCols)
      if (is.null(wtable)) return()

      # ── Dosage matrix + allele QC (annotates wtable in-place) ───────────
      dosage <- private$.buildDosageMatrix(snpCols, wtable, missing_st)
      if (is.null(dosage)) return()

      # ── Show SNP weights table (after QC so allele_status is complete) ──
      private$.fillSnpGridTable(wtable)

      # ── Weights vector aligned to dosage columns ─────────────────────────
      matched_rows <- wtable[wtable$rsid %in% colnames(dosage), , drop = FALSE]
      weights_vec  <- setNames(matched_rows$effect_weight, matched_rows$rsid)
      weights_vec  <- weights_vec[colnames(dosage)]

      # Drop SNPs with no catalog weight (NA) when using catalog mode
      na_snps <- names(weights_vec)[is.na(weights_vec)]
      if (length(na_snps) > 0) {
        keep_cols   <- names(weights_vec)[!is.na(weights_vec)]
        dosage      <- dosage[, keep_cols, drop = FALSE]
        weights_vec <- weights_vec[keep_cols]
      }

      if (length(weights_vec) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>No SNPs with valid weights — cannot compute PGS.<br/>
           Load a weights file or use SNP columns with numeric dosage data.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      scores <- as.numeric(dosage %*% weights_vec)
      if (normalize) scores <- scores / length(weights_vec)
      private$.pgsScores <- scores

      # ── Individual IDs ───────────────────────────────────────────────────
      if (!is.null(idCol) && nchar(idCol) > 0 && idCol %in% names(self$data)) {
        ids_full <- as.character(self$data[[idCol]])
        private$.idLabels <- if (!is.null(private$.keepMask))
                               ids_full[private$.keepMask]
                             else ids_full
      } else {
        private$.idLabels <- as.character(seq_along(scores))
      }

      private$.fillCoverageTable(snpCols, wtable, missing_st)
      private$.fillSummaryTable(scores)

      if (self$options$showScoreTable || self$options$showPercentiles)
        private$.fillScoreTable(scores, private$.idLabels)

      if (self$options$showPercentiles)
        private$.fillPercentileTable(scores)

      if (self$options$showAssoc && !is.null(respCol) && nchar(respCol) > 0)
        private$.fillAssocTable(scores, respCol)
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

      path   <- self$options$weightsPath
      sepOpt <- self$options$weightsSep

      if (!is.null(path) && nchar(trimws(path)) > 0 && file.exists(path))
        return(private$.parseCatalogFile(path, sepOpt, snpCols))

      # No file — unit weights for all selected columns
      private$.unitWeightTable(snpCols)
    },


    # ════════════════════════════════════════════════════════════════════════
    # .parseCatalogFile
    # ════════════════════════════════════════════════════════════════════════
    .parseCatalogFile = function(path, sepOpt, snpCols) {

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

      sep <- switch(sepOpt,
        comma = ",",
        tab   = "\t",
        {
          first   <- dataLines[1]
          n_tabs  <- lengths(regmatches(first, gregexpr("\t", first)))
          n_commas <- lengths(regmatches(first, gregexpr(",",  first)))
          if (n_tabs >= n_commas) "\t" else ","
        }
      )

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

      # Pre-flag ambiguous SNPs from catalog alleles alone
      ambig_flag <- private$.isAmbiguous(ea_vec, oa_vec)

      catalog <- data.frame(
        rsid          = as.character(df[[c_rsid]]),
        effect_allele = ea_vec,
        other_allele  = oa_vec,
        effect_weight = if (!is.null(c_weight)) suppressWarnings(as.numeric(df[[c_weight]]))
                        else                    rep(NA_real_, nrow(df)),
        chr           = if (!is.null(c_chr))    as.character(df[[c_chr]])    else "",
        pos           = if (!is.null(c_pos))    as.character(df[[c_pos]])    else "",
        matched       = TRUE,
        allele_status = ifelse(ambig_flag, "ambiguous", ""),   # filled in by QC step
        strand_flipped = FALSE,
        extra_cols    = extra_str,
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
        stringsAsFactors = FALSE
      )
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
          if (has_allele_info && wtable$allele_status[idx] == "") {
            wtable$allele_status[idx] <<- "dosage (orientation unverified)"
          } else if (!has_allele_info && wtable$allele_status[idx] == "") {
            wtable$allele_status[idx] <<- "ok (numeric dosage)"
          }
          next
        }

        # ── Character/factor column — allele strings ─────────────────────
        col_char <- toupper(trimws(as.character(col_raw)))
        col_char[col_char == "" | col_char == "NA"] <- NA_character_

        # Extract unique observed alleles (ignoring NA)
        obs_alleles <- unique(unlist(strsplit(
          col_char[!is.na(col_char)], "")))
        obs_alleles <- obs_alleles[nchar(obs_alleles) == 1]

        if (!has_allele_info) {
          # No catalog allele info — cannot count, mark NA
          mat[, snp] <- NA_real_
          wtable$allele_status[idx] <<- "no allele info in weights file"
          exclude <- c(exclude, snp)
          next
        }

        # Determine which allele to count as effect allele dosage,
        # and whether a strand flip is needed.
        ea_use    <- ea_cat
        flipped   <- FALSE
        status    <- ""

        if (is_ambiguous) {
          # AT/CG: cannot resolve strand — count EA as given, warn
          status <- "ambiguous (AT/CG) \u26a0"
        } else if (length(obs_alleles) == 0) {
          status <- "no valid genotypes observed"
        } else {
          ea_comp <- private$.complement(ea_cat)
          oa_comp <- private$.complement(oa_cat)

          alleles_match_direct     <- all(obs_alleles %in% c(ea_cat, oa_cat))
          alleles_match_complement <- all(obs_alleles %in% c(ea_comp, oa_comp))

          if (alleles_match_direct) {
            status  <- "ok"
          } else if (alleles_match_complement) {
            ea_use  <- ea_comp   # count complement of EA
            flipped <- TRUE
            status  <- "strand flip \u2194 (corrected)"
            wtable$strand_flipped[idx] <<- TRUE
          } else {
            status  <- "allele mismatch \u274c"
            wtable$allele_status[idx] <<- status
            exclude <- c(exclude, snp)
            next
          }
        }

        # Count copies of ea_use per genotype string
        dosage_col <- vapply(col_char, function(g) {
          if (is.na(g)) return(NA_real_)
          sum(strsplit(g, "")[[1]] == ea_use)
        }, numeric(1))

        mat[, snp] <- dosage_col
        wtable$allele_status[idx] <<- status
      }

      # Drop mismatched SNPs
      keep_cols <- setdiff(colnames(mat), exclude)
      mat       <- mat[, keep_cols, drop = FALSE]

      # Handle within-column missingness
      for (snp in colnames(mat)) {
        na_mask <- is.na(mat[, snp])
        if (any(na_mask)) {
          mat[, snp] <- switch(missing_st,
            mean    = { m <- mean(mat[, snp], na.rm = TRUE)
                        mat[na_mask, snp] <- if (is.nan(m)) 0 else m
                        mat[, snp] },
            zero    = { mat[na_mask, snp] <- 0; mat[, snp] },
            exclude = mat[, snp]
          )
        }
      }

      if (missing_st == "exclude") {
        keep <- complete.cases(mat)
        mat  <- mat[keep, , drop = FALSE]
        private$.keepMask <- keep
      } else {
        private$.keepMask <- rep(TRUE, nrow(self$data))
      }

      mat
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

      tbl$getColumn("chr")$setVisible(has_chr)
      tbl$getColumn("pos")$setVisible(has_pos)
      tbl$getColumn("effect_allele")$setVisible(has_ea)
      tbl$getColumn("other_allele")$setVisible(has_oa)
      tbl$getColumn("extra_cols")$setVisible(has_extra)

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
          extra_cols    = as.character(r$extra_cols)
        ))
      }
    },

    .fillCoverageTable = function(snpCols, wtable, missing_st) {
      inData    <- names(self$data)
      matched   <- intersect(wtable$rsid[wtable$matched], inData)
      ambiguous <- sum(private$.isAmbiguous(wtable$effect_allele, wtable$other_allele))
      flipped   <- sum(wtable$strand_flipped == TRUE, na.rm = TRUE)
      mismatch  <- sum(grepl("mismatch", wtable$allele_status, ignore.case = TRUE))

      meta <- attr(wtable, "pgs_meta") %||%
              list(pgs_id = "", pgs_name = "", trait_reported = "",
                   weight_type = "", genome_build = "", variants_number = "")

      self$results$coverageTable$setRow(rowNo = 1, values = list(
        pgs_id         = meta$pgs_id         %||% "",
        pgs_name       = meta$pgs_name       %||% "",
        trait_reported = meta$trait_reported %||% "",
        weight_type    = meta$weight_type    %||% "",
        genome_build   = meta$genome_build   %||% "",
        snpsInWeights  = sum(wtable$matched),
        snpsInData     = length(intersect(snpCols, inData)),
        snpsMatched    = length(matched),
        pctMatched     = if (sum(wtable$matched) > 0)
                           length(matched) / sum(wtable$matched)
                         else 0,
        snpsAmbiguous  = ambiguous,
        snpsFlipped    = flipped,
        snpsMismatch   = mismatch,
        missingStrategy = missing_st
      ))
    },

    .isAmbiguous = function(ea, oa) {
      pairs <- paste0(toupper(ea), toupper(oa))
      pairs %in% c("AT", "TA", "CG", "GC")
    },

    .fillSummaryTable = function(scores) {
      e1071_skew <- tryCatch(
        e1071::skewness(scores, na.rm = TRUE),
        error = function(e) {
          n  <- sum(!is.na(scores))
          mu <- mean(scores, na.rm = TRUE)
          s  <- sd(scores, na.rm = TRUE)
          if (s == 0 || n < 3) return(NA_real_)
          sum(((scores[!is.na(scores)] - mu) / s)^3) * n / ((n - 1) * (n - 2))
        }
      )
      self$results$summaryTable$setRow(rowNo = 1, values = list(
        n    = sum(!is.na(scores)),
        mean = mean(scores, na.rm = TRUE),
        sd   = sd(scores,   na.rm = TRUE),
        min  = min(scores,  na.rm = TRUE),
        max  = max(scores,  na.rm = TRUE),
        skew = e1071_skew
      ))
    },

    .fillScoreTable = function(scores, ids) {
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

    .fillAssocTable = function(scores, respCol) {
      if (!(respCol %in% names(self$data))) return()

      resp <- self$data[[respCol]]
      if (!is.null(private$.keepMask) &&
          length(private$.keepMask) == length(resp))
        resp <- resp[private$.keepMask]

      df_assoc <- data.frame(pgs = scores, resp = resp)
      df_assoc <- df_assoc[complete.cases(df_assoc), ]
      if (nrow(df_assoc) < 3) return()

      is_binary <- is.factor(df_assoc$resp) ||
                   (length(unique(df_assoc$resp)) == 2 &&
                    all(df_assoc$resp %in% c(0, 1, NA)))

      fit <- tryCatch({
        if (is_binary) glm(resp ~ pgs, data = df_assoc, family = binomial())
        else            lm(resp ~ pgs, data = df_assoc)
      }, error = function(e) NULL)

      if (is.null(fit)) return()
      coefs <- coef(summary(fit))
      if (nrow(coefs) < 2) return()

      self$results$assocTable$setRow(rowNo = 1, values = list(
        responseVar = respCol,
        model       = if (is_binary) "Logistic" else "Linear",
        beta        = coefs[2, 1],
        se          = coefs[2, 2],
        stat        = coefs[2, 3],
        p           = coefs[2, 4]
      ))
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
           main = "PGS Distribution",
           xlab = "PGS Score", ylab = "Density",
           ylim = c(0, y_max), las = 1)

      lines(dens, col = "#C0392B", lwd = 2.5)
      abline(v = mean(scores, na.rm = TRUE), col = "#27AE60", lwd = 2, lty = 2)
      legend("topright",
             legend = c("Density", "Mean"),
             col    = c("#C0392B", "#27AE60"),
             lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
      TRUE
    }

  )  # end private
)

`%||%` <- function(a, b) if (!is.null(a)) a else b
