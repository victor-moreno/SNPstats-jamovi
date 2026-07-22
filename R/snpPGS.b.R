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

    .keepMask  = NULL,
    .cache     = new.env(parent = emptyenv()),
    .plotVis   = NULL,   # plot visibility set in .init (guards .run re-touch)
    .assoc_acc = NULL,   # per-run row accumulator for assocTable (across modes)
    .inter_acc = NULL,   # per-run row accumulator for interactionTable (across modes)

    # ── Refresh helpers (ported from snpStats.b.R) ─────────────────────────────
    # jamovi rebuilds the analysis object on every option click and restores a
    # table's cells (Table$fromProtoBuf) only into rows that already exist by the
    # end of .init(). Rows added in .run() via addRow() are never restored, so a
    # table filled that way comes back empty (rowCount 0) and is forced to rebuild
    # on every click — the "flash". Fix: pre-create the predicted rows in .init()
    # (positional keys "1".."n"; labels live in columns), and in .run() write by
    # position with setRow into those rows, rebuilding only when the predicted
    # count is wrong. A wrong prediction costs a one-off blank, never wrong data.

    # Full data during .init: jamovi hands .init header-only data, so read the
    # real values (one dataset read per option click) to size the tables exactly.
    .init_data = function() {
      d <- tryCatch(self$data, error = function(e) NULL)
      if (!is.null(d) && nrow(d) > 0) return(d)          # provided (R / tests)
      tryCatch(self$readDataset(FALSE), error = function(e) d)
    },

    # Pre-create n empty positional rows ("1".."n") unless the table already has
    # rows (restore may have refilled them). Row labels live in columns, so .init
    # only predicts the COUNT; .run rewrites every cell.
    .pre_rows = function(tbl, n) {
      if (tbl$rowCount == 0L && n > 0L)
        for (i in seq_len(n)) tbl$addRow(rowKey = as.character(i))
    },

    # True when the table's existing rows can be setRow'd into rather than rebuilt.
    .reuse_rows = function(tbl, n) tbl$rowCount == n && n > 0L,

    # Clear a table's run-phase note keys as init=FALSE. A note set only in .run()
    # (via .fill*) is otherwise recreated by jmvcore as an empty init=TRUE
    # placeholder on the next cycle; because the fill is gated behind a fill flag,
    # a restored run does not re-enter it, so the placeholder renders as a blank
    # footnote. Declaring the keys init=FALSE here removes the placeholder; the
    # protobuf restore then refills only the notes the previous run actually saved.
    .clearRunNotes = function(tbl, keys) {
      for (k in keys) tbl$setNote(k, NULL, init = FALSE)
    },

    # A small link to the online tutorial (docs/ is not bundled into the
    # installed module). Shown only in the "Getting started" state (no SNP
    # columns assigned yet); hidden once the guidance disappears.
    .setHelpBanner = function() {
      has_snps <- !is.null(self$options$snpCols) && length(self$options$snpCols) > 0
      self$results$helpBanner$setVisible(!has_snps)
      self$results$helpBanner$setContent(paste0(
        "<div style=\"font-size:0.85em; color:#666; padding:2px 0 6px;\">",
        "\U0001F4D6 <a href=\"https://victor-moreno.github.io/SNPstats-jamovi/TUTORIAL.html\" ",
        "target=\"_blank\" rel=\"noopener\">SNPstats tutorial &amp; help</a></div>"))
    },

    # Column keys for the grid's dynamic "extra" catalog fields. Shared by .init
    # (which pre-creates them so restore keeps them) and .fillSnpGridTable.
    .gridExtraKeys = function(extra_names) {
      if (length(extra_names) == 0) return(character(0))
      make.unique(paste0("xc_", gsub("[^A-Za-z0-9]", "_", extra_names)))
    },

    # The weight table, memoised per (snpCols, weightsPath, weightingMode) so
    # .init can predict the grid's row/column shape without re-parsing the file
    # on every cycle when the analysis object is reused. .run builds its own via
    # the two-level cache; this only serves .init's pre-creation.
    .gridWeightTable = function(snpCols) {
      key <- list(snpCols = sort(snpCols),
                  weightsPath = self$options$weightsPath %||% "",
                  wmode = self$options$weightingMode)
      if (!identical(key, private$.cache$grid_wt_key)) {
        private$.cache$grid_wt_key <- key
        private$.cache$grid_wt <-
          tryCatch(private$.buildWeightTable(snpCols), error = function(e) NULL)
      }
      private$.cache$grid_wt
    },

    # Write a list of row value-lists positionally. Reuses the rows .init created
    # (setRow) when the count matches so an unrelated click updates in place; else
    # rebuilds once (deleteRows + addRow). setRow/addRow leave unlisted columns
    # untouched, so dynamic columns created in .init keep their restored cells.
    .writeRows = function(tbl, rows) {
      reuse <- private$.reuse_rows(tbl, length(rows))
      if (!reuse) tbl$deleteRows()
      for (i in seq_along(rows))
        (if (reuse) tbl$setRow else tbl$addRow)(
          rowKey = as.character(i), values = rows[[i]])
    },

    # ── .init row-count predictors ─────────────────────────────────────────────
    # Each mirrors the corresponding .fill*() so .init pre-creates exactly the
    # rows .run will produce. Best-effort: an over/under-count only triggers a
    # one-off rebuild (.reuse_rows falls back), never wrong data.

    # Response vector resolved the way .run does (binary-numeric -> factor,
    # caseLevel relevel to reference = first level), read from .init_data();
    # NULL when unusable.
    .initResp = function() {
      respCol <- self$options$responseCol
      if (is.null(respCol) || !nzchar(trimws(respCol %||% ""))) return(NULL)
      dat <- private$.init_data()
      if (is.null(dat) || !(respCol %in% names(dat))) return(NULL)
      resp <- dat[[respCol]]
      if (is.numeric(resp) && length(unique(resp[!is.na(resp)])) == 2 &&
          all(resp[!is.na(resp)] %in% c(0, 1)))
        resp <- factor(resp)
      # Character response: use data order of appearance, not alphabetical (see .run).
      if (is.character(resp))
        resp <- factor(resp, levels = unique(resp[!is.na(resp)]))
      cl <- trimws(self$options$caseLevel %||% "")
      if (nchar(cl) > 0) {
        # Factor level order (not appearance) so other levels keep their order.
        lvls_obs <- if (is.factor(resp)) levels(droplevels(resp))
                    else unique(as.character(resp[!is.na(resp)]))
        if (cl %in% lvls_obs)
          resp <- factor(resp, levels = c(cl, setdiff(lvls_obs, cl)))
      }
      resp
    },

    # Response type as .fillAssocTable / .fillPercentileTable classify it.
    .respType = function(resp) {
      if (is.null(resp)) return("none")
      n_lvls <- length(unique(resp[!is.na(resp)]))
      if (!is.factor(resp) && n_lvls > 5) "continuous"
      else if (n_lvls == 2)               "binary"
      else if (n_lvls > 2)                "polytomous"
      else                                "continuous"
    },

    .respLevels = function(resp) {
      if (is.null(resp)) return(character(0))
      levels(droplevels(factor(resp[!is.na(resp)])))
    },

    # Number of percentile categories from the percentileBreaks option.
    .pctNCats = function() {
      b <- suppressWarnings(as.numeric(
        strsplit(trimws(self$options$percentileBreaks %||% ""), ",")[[1]]))
      b <- sort(unique(b[!is.na(b) & b > 0 & b < 100]))
      if (length(b) == 0) b <- c(20, 40, 60, 80, 90, 95)
      length(b) + 1L
    },

    # assocTable rows across all modes (see .fillAssocTable).
    .assocNRows = function(resp, n_modes) {
      per <- switch(private$.respType(resp),
        binary     = 3L,                                          # logistic, t, MW
        continuous = 3L,                                          # linear, Pearson, Spearman
        polytomous = (length(private$.respLevels(resp)) - 1L) + 3L,  # OR/lvl + LRT + ANOVA + KW
        0L)
      n_modes * per
    },

    # percentileTable (regression) rows across all modes (see .fillPercentileTable).
    .pctCatNRows = function(resp, n_modes) {
      contrasts <- if (identical(private$.respType(resp), "polytomous"))
                     max(1L, length(private$.respLevels(resp)) - 1L) else 1L
      n_modes * contrasts * private$.pctNCats()
    },

    # interactionTable rows across all modes (see .fillInteractionTable).
    .interNRows = function(resp, n_modes) {
      covCols <- self$options$covCols
      dat     <- private$.init_data()
      if (is.null(covCols) || length(covCols) == 0 || is.null(dat)) return(0L)
      cov1_nm <- covCols[[1]]
      if (!(cov1_nm %in% names(dat))) return(0L)
      cov1    <- dat[[cov1_nm]]
      is_fac  <- is.factor(cov1) || is.character(cov1)
      k        <- if (is_fac) length(unique(as.character(cov1[!is.na(cov1)]))) else 0L
      cov_rows <- if (is_fac) max(1L, k - 1L) else 1L      # cov1 main / interaction coef rows
      nl <- length(private$.respLevels(resp))
      per <- if (identical(private$.respType(resp), "polytomous"))
               max(1L, nl - 1L) * (1L + cov_rows + cov_rows) + 1L   # PGS + cov main + cov int; +LRT
             else
               1L + cov_rows + cov_rows + 1L                         # PGS + cov main + cov int + LRT/F
      n_modes * per
    },

    # ════════════════════════════════════════════════════════════════════════
    # .init  — build table skeletons from options alone (no data access)
    #
    # Jamovi calls .init() immediately on every option change so the UI can
    # show the correct table structure before the heavy computation in .run()
    # has finished.  Nothing here touches self$data.
    # ════════════════════════════════════════════════════════════════════════
    .init = function() {

      private$.setHelpBanner()

      snpCols  <- self$options$snpCols
      respCol  <- self$options$responseCol
      wmode    <- self$options$weightingMode
      has_file <- !is.null(self$options$weightsPath) &&
                  nchar(trimws(self$options$weightsPath)) > 0

      has_snps <- !is.null(snpCols) && length(snpCols) > 0
      has_resp <- !is.null(respCol) && nchar(trimws(respCol)) > 0

      # ── Determine which scoring modes will be run ──────────────────────
      # (mirrors the logic in .run so the skeleton matches the output)
      run_weighted   <- wmode %in% c("weighted", "both") && has_file
      run_unweighted <- wmode %in% c("unweighted", "both") || !has_file

      mode_labels <- character(0)
      if (run_weighted)   mode_labels <- c(mode_labels, "Weighted")
      if (run_unweighted) mode_labels <- c(mode_labels, "Unweighted")
      if (length(mode_labels) == 0) mode_labels <- "Unweighted"

      show_mode_col <- length(mode_labels) > 1

      # ── snpGridTable ───────────────────────────────────────────────────
      # Refresh-safe like the other tables: pre-create the rows (positional
      # integer keys) and the dynamic "extra" columns predicted from the weights
      # file, so protobuf restore refills the cells and an unrelated option click
      # leaves the grid untouched (.fillSnpGridTable writes by position via
      # .writeRows). The grid row set is the full weight table (catalog SNPs may
      # exceed the selected columns), so predict it from .gridWeightTable.
      snpGrid <- self$results$snpGridTable
      snpGrid$setVisible(has_snps && isTRUE(self$options$showSnpGrid))

      if (has_snps) {
        gwt         <- private$.gridWeightTable(snpCols)
        n_grid      <- if (!is.null(gwt)) nrow(gwt) else length(snpCols)
        extra_names <- if (!is.null(gwt))
                         attr(gwt, "extra_col_names") %||% character(0)
                       else character(0)
        # Dynamic extra columns must exist by end of .init (columns added only in
        # .run are never restored, like rows).
        xc_keys <- private$.gridExtraKeys(extra_names)
        for (i in seq_along(xc_keys)) {
          if (is.null(tryCatch(snpGrid$getColumn(xc_keys[i]), error = function(e) NULL)))
            snpGrid$addColumn(name = xc_keys[i], title = extra_names[i], type = "text")
        }
        private$.pre_rows(snpGrid, n_grid)

        # Column visibility is option-driven; hide columns that cannot have
        # values when no weights file is loaded.
        snpGrid$getColumn("effect_weight")$setVisible(TRUE)   # always show: 1 when unweighted
        snpGrid$getColumn("effect_allele")$setVisible(has_file)
        snpGrid$getColumn("other_allele")$setVisible(has_file)
        snpGrid$getColumn("chr")$setVisible(has_file)
        snpGrid$getColumn("pos")$setVisible(has_file)
        # QC filter columns visible only when at least one QC filter is active
        has_qc_flt <- isTRUE(self$options$qcFilterMissing) ||
                      isTRUE(self$options$qcFilterHwe)
        snpGrid$getColumn("qc_excl_reason")$setVisible(has_qc_flt)
      }

      # Data-derived predictions below need the response resolved as .run does;
      # read it once (one dataset read) and reuse for summary/assoc/percentile/
      # interaction row-count predictions.
      covCols   <- self$options$covCols
      has_covs  <- !is.null(covCols) && length(covCols) > 0
      init_resp <- if (has_resp) private$.initResp() else NULL
      n_modes   <- length(mode_labels)

      # ── coverageTable — key/value summary rows ─────────────────────────
      # The row SET is option-determined (see .coverageNRows); pre-create it so
      # restore refills the cells and an unrelated click does not rebuild.
      covTbl <- self$results$coverageTable
      covTbl$setVisible(has_snps && isTRUE(self$options$showCoverage))
      if (has_snps && isTRUE(self$options$showCoverage))
        private$.pre_rows(covTbl, private$.coverageNRows(
          has_file, private$.coverageMetaNonEmpty(self$options$weightsPath)))

      # ── summaryTable — one row per scoring mode × response group ──────
      # We can pre-build skeleton rows: one "Overall" row per mode.
      # If a binary response is selected a second stratified row per mode
      # will be added in .run(), but the skeleton still gives instant feedback.
      summaryTbl <- self$results$summaryTable
      summaryTbl$setVisible(has_snps && isTRUE(self$options$showSummary))

      if (has_snps && summaryTbl$rowCount == 0) {
        # Pre-create EVERY row .run() will produce (one per group per mode), not
        # just the Overall rows. A binary/categorical response also yields a row
        # per level; if only the Overall rows are seeded, restore refills them
        # (their keys match), the isNotFilled() gate then reads "filled" and the
        # stratified rows are never re-added — the table silently drops to Overall
        # only on any unrelated click. Seeding all rows lets restore refill them
        # and the gate correctly skips (no rebuild, no data loss).
        groups <- private$.predictSummaryGroups(init_resp)
        for (ml in mode_labels)
          for (g in groups)
            summaryTbl$addRow(
              rowKey = paste0(ml, "_", g),
              values = list(score_type = ml, group = g))
      }

      # Pre-create the predicted assoc rows (one block per mode; see .assocNRows).
      # Rows added only in .run() are never restored, so restore refills these and
      # .run writes them by position.
      assocTbl <- self$results$assocTable
      assocTbl$setVisible(has_snps && has_resp && isTRUE(self$options$showAssoc))
      if (has_snps && has_resp && isTRUE(self$options$showAssoc))
        private$.pre_rows(assocTbl, private$.assocNRows(init_resp, n_modes))
      # Re-declare the run-phase notes as init=FALSE so jmvcore does not leave an
      # empty init placeholder that a gated .run() cannot clear (protobuf restore
      # then repopulates any note the previous run actually saved). See
      # .fillAssocTable / .clearRunNotes.
      private$.clearRunNotes(assocTbl,
        c("covNote", "respNote", "fitWarning_Weighted", "fitWarning_Unweighted"))

      # ── percentileTable / percentileThreshTable ────────────────────────
      show_pct <- has_snps && isTRUE(self$options$showPercentiles)
      pctThr <- self$results$percentileThreshTable
      pctCat <- self$results$percentileTable
      pctThr$setVisible(show_pct)
      pctCat$setVisible(show_pct)
      private$.clearRunNotes(pctCat, c("modelNote", "refNote", "covNote"))
      if (show_pct) {
        # Dynamic per-level columns for the counts table must exist by end of
        # .init (columns added in .run are never restored, like rows).
        if (private$.respType(init_resp) %in% c("binary", "polytomous")) {
          for (lv in private$.respLevels(init_resp)) {
            nm <- paste0("n_lv_", make.names(lv))
            if (is.null(tryCatch(pctThr$getColumn(nm), error = function(e) NULL)))
              pctThr$addColumn(name = nm, title = as.character(lv), type = "text")
          }
        }
        private$.pre_rows(pctThr, n_modes * private$.pctNCats())
        private$.pre_rows(pctCat, private$.pctCatNRows(init_resp, n_modes))
      }

      # ── interactionTable ───────────────────────────────────────────────
      interTbl <- self$results$interactionTable
      interTbl$setVisible(
        has_snps && has_resp && has_covs && isTRUE(self$options$showInteraction))
      if (has_snps && has_resp && has_covs && isTRUE(self$options$showInteraction))
        private$.pre_rows(interTbl, private$.interNRows(init_resp, n_modes))
      private$.clearRunNotes(interTbl, "intNote")



      # ── Plot visibility ──────────────────────────────────────────────────
      # Set here (not in .run) from the predicted response type: an image touched
      # during .run() is re-rendered by the engine every run, so setting the
      # response-dependent plots' visibility in .run made them re-render on every
      # unrelated click (e.g. toggling another plot) while distPlot — set only
      # here — did not. .run keeps a guarded correction (.setPlotVis) for the
      # header-only-.init case, which only touches an image on a real change.
      rtype   <- private$.respType(init_resp)
      has_resp2 <- has_snps && !identical(rtype, "none")
      init_vis <- list(
        distPlot   = has_snps   && isTRUE(self$options$showDistPlot),
        stratPlot  = has_resp2  && isTRUE(self$options$showDistPlot) &&
                       rtype %in% c("continuous", "polytomous"),
        forestPlot = has_resp2  && isTRUE(self$options$showForestPlot),
        rocPlot    = has_resp2  && isTRUE(self$options$showRocPlot) &&
                       rtype %in% c("binary", "polytomous"),
        calibPlot  = has_resp2  && isTRUE(self$options$showCalibPlot) &&
                       rtype %in% c("binary", "polytomous"))
      for (nm in names(init_vis)) self$results$get(nm)$setVisible(init_vis[[nm]])
      private$.plotVis <- init_vis

      # ── Plot sizes: set here (before the first render) so the image is not
      # rendered at the default size and then upscaled/resized. Same formula and
      # facet layout the render functions use (.pgs_plot_size). n_modes columns
      # for the row plots; ROC/calibration add a comparison grid when polytomous.
      n_comps <- if (identical(rtype, "polytomous"))
                   max(1L, length(private$.respLevels(init_resp)) - 1L) else 1L
      set_sz <- function(img, nc, nr) {
        sz <- .pgs_plot_size(self$options, nc, nr); img$setSize(sz[["w"]], sz[["h"]])
      }
      set_sz(self$results$distPlot,   n_modes, 1L)
      set_sz(self$results$stratPlot,  n_modes, 1L)
      # forest facets by contrast (mode × comparison) when polytomous, like ROC.
      if (n_comps == 1L) {
        set_sz(self$results$forestPlot, n_modes, 1L)
        set_sz(self$results$rocPlot,    n_modes, 1L)
        set_sz(self$results$calibPlot,  n_modes, 1L)
      } else {
        set_sz(self$results$forestPlot, n_comps, n_modes)
        set_sz(self$results$rocPlot,    n_comps, n_modes)
        set_sz(self$results$calibPlot,  n_comps, n_modes)
      }
    },

    # ════════════════════════════════════════════════════════════════════════
    # .run
    # ════════════════════════════════════════════════════════════════════════
    .run = function() {

      snpCols     <- self$options$snpCols
      covCols     <- self$options$covCols
      respCol     <- self$options$responseCol
      missing_st  <- self$options$missingStrategy
      wmode       <- self$options$weightingMode

      # ── Nothing selected yet: show guidance and return ───────────────────
      # Table/plot visibility was already set to FALSE by .init().
      if (is.null(snpCols) || length(snpCols) == 0) {
        self$results$validationMsg$setContent(
          "<div style='color:#555; padding:6px 0;'>
             <b>Getting started:</b><br>
             \u2022 Drag one or more SNP columns into <i>SNP columns</i>.<br>
             \u2022 Optionally load a PGS Catalog weights file (.csv / .tsv) for weighted scoring. An example (CRCgenet-PGS.txt) is provided in the data folder of the package<br>
             \u2022 Optionally select a response variable to test association / interaction.
           </div>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      # ── Two-level cache for the SNP-QC pipeline ────────────────────────────
      # Full key: reuse the filtered result wholesale when nothing SNP/QC/missing
      # related changed (e.g. only a response/covariate/plot option moved).
      snp_qc_key <- list(
        snpCols     = sort(snpCols),
        weightsPath = self$options$weightsPath %||% "",
        wmode       = wmode,
        missing_st  = missing_st,
        qcMiss      = isTRUE(self$options$qcFilterMissing),
        qcMissPct   = self$options$qcMaxMissingPct,
        qcHwe       = isTRUE(self$options$qcFilterHwe),
        qcHweP      = self$options$qcHweP,
        hweResp     = if (isTRUE(self$options$qcFilterHwe)) respCol else NULL
      )

      if (identical(snp_qc_key, private$.cache$snp_qc_key) &&
          !is.null(private$.cache$qc)) {
        qc                <- private$.cache$qc
        private$.keepMask <- private$.cache$keepMask
      } else {
        # Core key: the heavy per-SNP compute (cleaning, dosage, AF/HWE stats) is
        # independent of the QC thresholds and the missing strategy, so a filter
        # or imputation change reuses it and only re-runs the light filter step.
        # HWE-in-controls makes it depend on the HWE flag + response, not the p.
        apply_hwe <- isTRUE(self$options$qcFilterHwe) && !is.na(self$options$qcHweP)
        core_key  <- list(
          snpCols     = sort(snpCols),
          weightsPath = self$options$weightsPath %||% "",
          wmode       = wmode,
          applyHwe    = apply_hwe,
          hweResp     = if (apply_hwe) respCol else NULL
        )
        if (identical(core_key, private$.cache$core_key) &&
            !is.null(private$.cache$core)) {
          core <- private$.cache$core
        } else {
          wtable <- private$.buildWeightTable(snpCols)
          if (is.null(wtable)) return()
          core <- private$.buildDosageCore(snpCols, wtable)
          if (is.null(core)) return()
          private$.cache$core_key <- core_key
          private$.cache$core     <- core
        }

        # Reset keepMask before the filter step sets it (exclude strategy).
        private$.keepMask <- rep(TRUE, nrow(self$data))
        qc <- private$.applyDosageFilter(core, missing_st)
        if (is.null(qc)) return()

        private$.cache$snp_qc_key <- snp_qc_key
        private$.cache$qc         <- qc
        private$.cache$keepMask   <- private$.keepMask
      }

      dosage     <- qc$mat
      wtable     <- qc$wtable
      valid_snps <- qc$valid_snps

      # Note: qc$valid_counts is already row-filtered inside buildDosageMatrix
      # for the 'exclude' strategy, so it always matches qc$mat row count on
      # return. No additional filtering is needed here.

      # Grid refill gate — the same contract as every other table now that the
      # grid is refresh-safe (positional keys + .writeRows + .init pre-creation).
      # Everything that changes the grid's content, order or visible rows — the
      # scored SNP set, the QC filters, and the sort/filter display options — is
      # in the grid's clearWith, so a change sets isNotFilled() (which survives
      # protobuf restore) and an unrelated click leaves it FALSE → the restored
      # rows are shown without a rebuild.
      if (self$results$snpGridTable$isNotFilled() ||
          self$results$snpGridTable$rowCount == 0) {
        private$.fillSnpGridTable(wtable, valid_snps)
      }

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
        # A character response carries no stored level order, so R's factor()
        # would sort it alphabetically. Use order of appearance in the data
        # instead, so the reference is the first observed level. Existing factors
        # (ordinal or user-ordered) already carry their intended order and are
        # left untouched.
        if (is.character(resp))
          resp <- factor(resp, levels = unique(resp[!is.na(resp)]))
        # Apply caseLevel: relevel so the selected level is FIRST (= reference).
        # Every model/plot treats lvls[1] as the baseline — binary glm compares
        # lvls[2] vs lvls[1], multinom uses lvls[1] as the reference category, and
        # the plots set ref_lv <- lvls[1] — so a single relevel propagates to all
        # tables and plots. Reference (not event) is used because it generalizes to
        # the polytomous case, where there is no single event level.
        cl <- trimws(self$options$caseLevel %||% "")
        if (nchar(cl) > 0) {
          # Use the factor's level order (not row order of appearance) so moving
          # the reference to the front preserves the order of the other levels.
          lvls_obs <- if (is.factor(resp)) levels(droplevels(resp))
                      else unique(as.character(resp[!is.na(resp)]))
          if (cl %in% lvls_obs) {
            other <- setdiff(lvls_obs, cl)
            resp  <- factor(resp, levels = c(cl, other))
          }
        }
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

      # isNotFilled() fires after clearWith; rowCount==0 catches first run
      if (self$results$coverageTable$isNotFilled() ||
          self$results$coverageTable$rowCount == 0)
        private$.fillCoverageTable(snpCols, wtable, missing_st, valid_snps,
                                    valid_counts  = qc$valid_counts,
                                    n_total       = nrow(self$data),
                                    has_file      = has_file,
                                    n_qc_excl     = qc$n_qc_excl %||% 0L,
                                    n_flt_excl    = qc$n_flt_excl %||% 0L)

      # ── Fill tables only when clearWith has marked them as not filled ──────
      # isNotFilled() is TRUE after clearWith fires or on first run, FALSE once a
      # table's pre-created rows (from .init) have been restored — so an unrelated
      # click skips the fill and the restored rows are shown without a rebuild.
      # rowCount == 0 is the belt-and-suspenders first-run guard.
      need_summary <- self$results$summaryTable$isNotFilled()
      need_assoc   <- self$results$assocTable$isNotFilled() ||
                      self$results$assocTable$rowCount == 0
      need_inter   <- self$results$interactionTable$isNotFilled() ||
                      self$results$interactionTable$rowCount == 0
      if (need_summary) self$results$summaryTable$deleteRows()

      do_assoc <- need_assoc && self$options$showAssoc && !is.null(resp)
      do_inter <- need_inter && self$options$showInteraction && !is.null(resp) &&
                  !is.null(covs) && ncol(covs) > 0

      # assocTable and interactionTable are one table each fed from every scoring
      # mode; accumulate their rows across the loop, then write once positionally
      # into the rows .init pre-created (see .writeRows).
      private$.assoc_acc <- list()
      private$.inter_acc <- list()

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

        if (need_summary)
          private$.fillSummaryTable(scores, resp, mode_label)

        if (do_assoc)
          private$.fillAssocTable(scores, resp, respCol, covs, mode_label)

        if (do_inter)
          private$.fillInteractionTable(scores, resp, respCol, covs, mode_label)
      }

      # Write the accumulated rows in one positional pass (reuse .init's rows).
      if (do_assoc) private$.writeRows(self$results$assocTable, private$.assoc_acc)
      if (do_inter) private$.writeRows(self$results$interactionTable, private$.inter_acc)

      if (length(all_scores) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>No SNPs with valid weights — cannot compute PGS.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      # ── Cache everything the plot render functions need ─────────────────
      # Always populated so any plot can render independently of the others.
      private$.cache$all_scores <- all_scores
      private$.cache$resp       <- resp
      private$.cache$respCol    <- respCol
      private$.cache$covs       <- covs

      # ── Reconcile plot visibility with the ACTUAL response type ───────────
      # .init already set visibility from the predicted type. Only correct here
      # when the actual differs (e.g. .init ran header-only), and touch each image
      # only on a real change — an image touched in .run is re-rendered by the
      # engine, so re-asserting the same visibility would refresh it on every
      # unrelated click. Numerical responses stay off ROC/calibration (huge grids).
      n_resp_lvls    <- if (!is.null(resp)) length(unique(resp[!is.na(resp)])) else 0L
      has_resp       <- !is.null(resp)
      is_binary_resp <- has_resp && (is.factor(resp) || n_resp_lvls == 2)
      is_cat_resp    <- has_resp && (is.factor(resp) || n_resp_lvls <= 5)
      run_vis <- list(
        stratPlot  = has_resp && self$options$showDistPlot && !is_binary_resp,
        forestPlot = has_resp && isTRUE(self$options$showForestPlot),
        rocPlot    = is_cat_resp && isTRUE(self$options$showRocPlot),
        calibPlot  = is_cat_resp && isTRUE(self$options$showCalibPlot))
      for (nm in names(run_vis)) {
        if (!identical(run_vis[[nm]], private$.plotVis[[nm]])) {
          self$results$get(nm)$setVisible(run_vis[[nm]])
          private$.plotVis[[nm]] <- run_vis[[nm]]
        }
      }

      if (self$options$showPercentiles &&
          (self$results$percentileTable$isNotFilled()      || self$results$percentileTable$rowCount == 0 ||
           self$results$percentileThreshTable$isNotFilled() || self$results$percentileThreshTable$rowCount == 0))
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
      extra_names <- setdiff(orig_names, known_cols)
      # safe R column names for extra fields (xc_ prefix)
      xc_keys <- if (length(extra_names) > 0)
                   make.unique(paste0("xc_", gsub("[^A-Za-z0-9]", "_", extra_names)))
                 else character(0)

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
        pgs_row_order  = seq_len(nrow(df)),
        n_missing     = NA_integer_,
        pct_missing   = NA_real_,
        effect_af     = NA_real_,
        hwe_p         = NA_real_,
        stringsAsFactors = FALSE
      )
      # Add individual extra columns (one per extra PGS field)
      for (i in seq_along(extra_names))
        catalog[[xc_keys[i]]] <- as.character(df[[extra_names[i]]])


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
          pgs_row_order  = NA_integer_,
          n_missing = NA_integer_,
          pct_missing = NA_real_,
          effect_af = NA_real_,
          hwe_p = NA_real_,
          selected_flag = TRUE,
          stringsAsFactors = FALSE
        )
        for (xk in xc_keys) extra_rows[[xk]] <- ""
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

      # Restore attributes dropped by rbind/subsetting
      attr(result, "pgs_meta")        <- meta
      attr(result, "extra_col_names") <- extra_names
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
        pgs_row_order  = seq_along(snpCols),
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

      rows <- lapply(seq_along(snpCols), function(idx) {
        snp <- snpCols[idx]
        col <- dat[[snp]]
        ea  <- ""
        oa  <- ""

        if (is.numeric(col)) {
          # Numeric dosage: allele inference not possible
          qc <- "\u26A0 numeric dosage"
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

          qc <- "\u2705 unit weight"
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
          pgs_row_order  = idx,
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
    # .buildDosageCore  (heavy, threshold-independent — cached across QC-filter
    # and missing-strategy changes; see .run's two-level cache)
    #
    #   1. Per-SNP column cleaning (.cleanGenotypeColumn)
    #   2. Input-type detection (.isNumericDosage — robust version)
    #   3. Allele QC: mismatch, multiallelic, monomorphic, no-allele-info
    #   4. Dosage conversion
    #   5. Per-SNP statistics: n_missing, pct_missing, effect_af, hwe_p
    #      (computed on the CLEANED column, before imputation)
    #
    # Returns the pre-imputation dosage matrix (all SNP columns), the wtable with
    # stats + base allele_status, the allele-QC exclusion set and the HWE control
    # note. The value depends only on snpCols / weightsPath / weightingMode and —
    # because HWE is tested in controls when the HWE filter is on — on the HWE
    # filter flag and response variable. It is INDEPENDENT of the QC thresholds
    # and the missing-value strategy, which .applyDosageFilter re-applies cheaply.
    # ════════════════════════════════════════════════════════════════════════
    .buildDosageCore = function(snpCols, wtable) {

      useCols <- intersect(wtable$rsid, names(self$data))

      if (length(useCols) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:#c0392b;'>None of the selected SNP columns are present in the dataset.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return(NULL)
      }

      # HWE is computed in controls only when the HWE filter is active; that makes
      # the per-SNP hwe_p depend on the flag + response, but not on the threshold.
      apply_hwe   <- isTRUE(self$options$qcFilterHwe) && !is.na(self$options$qcHweP)

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

      # Allele/monomorphic exclusions (threshold-independent); missingness/HWE
      # filter exclusions are applied later in .applyDosageFilter.
      qc_exclude    <- character(0)

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
          if (n > 0) paste0("; ", n, "×0/0→NA") else ""
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
              "\u274c mono dosage: ", obs_str, null_sfx(n_null_fix))
          } else {
            # \u26A0: allele QC cannot be performed on numeric dosage columns
            n_oor <- sum(!is.na(col_num) & (col_num < 0 | col_num > 2))
            num_actions <- c(
              "numeric dosage",
              if (n_null_fix > 0) paste0(n_null_fix, "×0/0→NA"),
              if (n_oor      > 0) paste0(n_oor,      "×OOR→NA")
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
            "\u274c no ACGT alleles", null_sfx(n_null_fix))
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
            "\u274c no allele info: ", obs_str, null_sfx(n_null_fix))
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
            "\u274c mismatch: ", obs_str, "\u2260", ea_cat, "/", oa_cat,
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

        # ── Dosage computation (vectorised) ──────────────────────────────────
        # Strip every non-ACGT character from all genotypes at once, then count
        # how many of the two remaining bases equal the effect allele. Replaces
        # a per-row strsplit (perl look-around) that was ~half the PGS runtime.
        # A cleaned genotype must reduce to exactly two ACGT bases; anything else
        # (NA, malformed, wrong length) → NA, matching the previous behaviour.
        bases    <- gsub("[^ACGT]", "", col_clean)
        two_base <- !is.na(col_clean) & nchar(bases) == 2L
        dos      <- rep(NA_real_, length(col_clean))
        dos[two_base] <- (substr(bases[two_base], 1L, 1L) == ea_use) +
                         (substr(bases[two_base], 2L, 2L) == ea_use)

        mat[, snp] <- dos

        # ── Monomorphism check ───────────────────────────────────────────────
        if (length(unique(dos[!is.na(dos)])) <= 1) {
          qc_exclude <- c(qc_exclude, snp)
          wtable$allele_status[idx] <- paste0(
            "\u274c mono: ", obs_str, null_sfx(n_null_fix))
        } else {
          # \u2705 clean pass  /  \u26A0 kept with actions
          actions_snp <- c(
            if (strand_flipped_snp) "flip",
            if (is_ambiguous_snp)   "ambig",
            if (n_null_fix > 0)     paste0(n_null_fix, "×0/0→NA")
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

      list(
        mat        = mat,          # pre-imputation, all SNP columns
        wtable     = wtable,       # stats + base allele_status
        qc_exclude = qc_exclude,   # allele/monomorphic exclusions
        ctrl_note  = ctrl_note,
        n_rows     = n_rows
      )
    },

    # ════════════════════════════════════════════════════════════════════════
    # .applyDosageFilter  (light — QC-threshold + missing-strategy dependent)
    #
    #   6. Configurable QC filters: missingness threshold, HWE p threshold
    #      (applied AFTER stats, so stats are always shown)
    #   7. Missing-value imputation for keep_snps
    #   8. Row exclusion if missing_st == "exclude"  (sets private$.keepMask)
    #
    # Operates on a cached .buildDosageCore() result; copies its mat / wtable so
    # the cache is never mutated. Returns the same list .buildDosageMatrix used
    # to (mat, wtable, valid_counts, valid_snps, n_qc_excl, n_flt_excl).
    # ════════════════════════════════════════════════════════════════════════
    .applyDosageFilter = function(core, missing_st) {

      mat        <- core$mat
      wtable     <- core$wtable
      qc_exclude <- core$qc_exclude
      ctrl_note  <- core$ctrl_note
      n_rows     <- core$n_rows

      miss_thresh <- self$options$qcMaxMissingPct
      hwe_thresh  <- self$options$qcHweP
      apply_miss  <- isTRUE(self$options$qcFilterMissing) && !is.na(miss_thresh)
      apply_hwe   <- isTRUE(self$options$qcFilterHwe)      && !is.na(hwe_thresh)

      filter_exclude <- character(0)

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
            prev <- sub("^[\u2705\u26A0\u274c]\\s*", "", wtable$allele_status[idx])
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
            prev <- sub("^[\u2705\u26A0\u274c]\\s*", "", wtable$allele_status[idx])
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
      inData <- names(self$data)

      sort_by    <- self$options$snpGridSortBy  %||% "as_selected"
      sort_field <- trimws(self$options$snpGridSortField %||% "")

      extra_names <- attr(wtable, "extra_col_names") %||% character(0)
      xc_keys     <- private$.gridExtraKeys(extra_names)

      if (isTRUE(self$options$filterValidSnps))
        wtable <- wtable[wtable$rsid %in% valid_snps, , drop = FALSE]

      # \u2500\u2500 Sort \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      if (nchar(sort_field) > 0 && length(extra_names) > 0) {
        # TextBox takes priority: case-insensitive match on original extra name
        m <- match(tolower(sort_field), tolower(extra_names))
        if (!is.na(m)) {
          vals     <- wtable[[xc_keys[m]]]
          num_vals <- suppressWarnings(as.numeric(vals))
          wtable   <- if (!all(is.na(num_vals[!is.na(vals)])))
                        wtable[order(num_vals, na.last = TRUE), , drop = FALSE]
                      else
                        wtable[order(vals,     na.last = TRUE), , drop = FALSE]
        }
      } else {
        ord <- switch(sort_by,
          pgs_order     = order(wtable$pgs_row_order, na.last = TRUE),
          effect_weight = order(wtable$effect_weight, decreasing = TRUE, na.last = TRUE),
          effect_af     = order(wtable$effect_af,     decreasing = TRUE, na.last = TRUE),
          pct_missing   = order(wtable$pct_missing,   decreasing = TRUE, na.last = TRUE),
          hwe_p         = order(wtable$hwe_p,                            na.last = TRUE),
          allele_status = order(wtable$allele_status,                    na.last = TRUE),
          {
            # as_selected: order by position in snpCols
            snp_order <- self$options$snpCols %||% character(0)
            sel_pos   <- match(wtable$rsid, snp_order)
            order(is.na(sel_pos), sel_pos, na.last = TRUE)
          }
        )
        wtable <- wtable[ord, , drop = FALSE]
      }

      # \u2500\u2500 Manage dynamic extra columns \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      # Add new xc_ columns; hide stale ones left over from a different file.
      existing <- sapply(tbl$columns, function(c) c$name)
      for (i in seq_along(extra_names)) {
        if (!(xc_keys[i] %in% existing)) {
          tbl$addColumn(name = xc_keys[i], title = extra_names[i], type = "text")
        } else {
          tbl$getColumn(xc_keys[i])$setTitle(extra_names[i])
        }
      }
      all_cols <- sapply(tbl$columns, function(c) c$name)
      for (cn in all_cols[startsWith(all_cols, "xc_")])
        tbl$getColumn(cn)$setVisible(cn %in% xc_keys)

      # \u2500\u2500 Fixed column visibility \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      tbl$getColumn("chr")$setVisible(any(nchar(wtable$chr) > 0, na.rm = TRUE))
      tbl$getColumn("pos")$setVisible(any(nchar(wtable$pos) > 0, na.rm = TRUE))
      tbl$getColumn("effect_allele")$setVisible(any(nchar(wtable$effect_allele) > 0, na.rm = TRUE))
      tbl$getColumn("other_allele")$setVisible(any(nchar(wtable$other_allele)  > 0, na.rm = TRUE))
      tbl$getColumn("n_missing")$setVisible(any(!is.na(wtable$n_missing)))
      tbl$getColumn("pct_missing")$setVisible(any(!is.na(wtable$n_missing)))
      tbl$getColumn("effect_af")$setVisible(any(!is.na(wtable$effect_af)))
      tbl$getColumn("hwe_p")$setVisible(any(!is.na(wtable$hwe_p)))
      has_qc_flt <- isTRUE(self$options$qcFilterMissing) || isTRUE(self$options$qcFilterHwe)
      tbl$getColumn("qc_excl_reason")$setVisible(has_qc_flt)

      # \u2500\u2500 Add rows \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      # Write rows positionally: .writeRows setRow's into the .init pre-created
      # rows when the count matches (refresh-safe, no rebuild) and only rebuilds
      # when it differs. Keys are "1".."n", matching .pre_rows.
      rows <- lapply(seq_len(nrow(wtable)), function(i) {
        r      <- wtable[i, ]
        status <- as.character(r$allele_status)
        in_ds  <- r$rsid %in% inData

        # QC filter column: derive from actual Unicode icons in allele_status
        qc_flt_reason <- if (grepl("^\u274c excl", status)) {
          sub("^\u274c excl ([^:]+):.*", "\u274c excl \\1", status)
        } else if (r$rsid %in% valid_snps) {
          if (startsWith(status, "\u26A0")) "\u26A0 pass (with actions)"
          else "\u2705 pass"
        } else ""

        extra_vals <- setNames(
          lapply(xc_keys, function(xk) {
            v <- r[[xk]]; if (is.null(v) || is.na(v)) "" else as.character(v)
          }),
          xc_keys
        )

        c(
          list(
            rsid          = as.character(r$rsid),
            chr           = as.character(r$chr),
            pos           = as.character(r$pos),
            effect_allele = as.character(r$effect_allele),
            other_allele  = as.character(r$other_allele),
            effect_weight = if (is.na(r$effect_weight)) '' else r$effect_weight,
            matched       = if (in_ds) "\u2713" else "\u2717",
            allele_status = status,
            n_missing     = if (is.na(r$n_missing))   '' else as.integer(r$n_missing),
            pct_missing   = if (is.na(r$pct_missing)) '' else r$pct_missing,
            effect_af     = if (is.na(r$effect_af))   '' else r$effect_af,
            hwe_p         = if (is.na(r$hwe_p))       '' else r$hwe_p,
            qc_excl_reason = qc_flt_reason
          ),
          extra_vals
        )
      })
      private$.writeRows(tbl, rows)
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
      null_fix  <- sum(grepl("0/0→NA", wtable$allele_status, fixed = TRUE))
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

      # Build an ORDERED row list, then write it positionally into the rows .init
      # pre-created (see .writeRows / .coverageNRows). The row SET is a pure
      # function of the options {has_file, weightingMode, qcFilterMissing,
      # qcFilterHwe} plus the count of non-empty header metadata fields, so .init
      # predicts it exactly and an unrelated click does not rebuild the table.
      rows <- list()
      add <- function(field, value) {
        rows[[length(rows) + 1L]] <<-
          list(field = field, value = as.character(value))
      }
      # Metadata rows are skipped when the file omits the field (same non-empty
      # test .coverageMetaNonEmpty uses in .init, so the counts agree).
      add_meta <- function(field, value) {
        v <- as.character(value)
        if (nchar(trimws(v)) > 0) add(field, v)
      }

      # Score metadata (from file header)
      add_meta("PGS ID",       meta$pgs_id         %||% "")
      add_meta("Score name",   meta$pgs_name       %||% "")
      add_meta("Trait",        meta$trait_reported %||% "")
      add_meta("Weight type",  meta$weight_type    %||% "")
      add_meta("Genome build", meta$genome_build   %||% "")

      # SNP counts \u2014 layout differs between weighted (file present) and unweighted
      if (has_file) {
        add("SNPs in weights file",       n_weights)
        add("SNPs in dataset",            n_indata)
        add("SNPs matched",               paste0(n_matched, " (", pct, "%)"))
        add("Ambiguous SNPs (AT/CG)",     ambiguous)
        add("Strand flipped (corrected)", flipped)
        add("Allele mismatch (excluded)", mismatch)
      } else {
        add("SNPs selected",              n_indata)
      }
      # QC-count rows are always emitted when their filter is active (showing 0
      # rather than appearing/disappearing) so the row set is option-determined.
      add("Null-allele genotypes fixed (0/0 \u2192 NA)", null_fix)
      add("Excluded by allele/monomorphic QC", n_qc_excl)
      if (isTRUE(self$options$qcFilterMissing))
        add("Excluded by missingness filter",
            sum(grepl("\u274c excl (missing", wtable$allele_status, fixed = TRUE)))
      if (isTRUE(self$options$qcFilterHwe))
        add("Excluded by HWE filter",
            sum(grepl("\u274c excl (HWE", wtable$allele_status, fixed = TRUE)))
      add("SNPs used in score",                length(valid_snps))
      add("Missing genotype strategy",         missing_st)
      run_w  <- self$options$weightingMode %in% c("weighted", "both") && has_file
      run_uw <- self$options$weightingMode %in% c("unweighted", "both") || !has_file
      if (run_w)  add("Score scale (Weighted)",   private$.scaleLabel("Weighted"))
      if (run_uw) add("Score scale (Unweighted)", private$.scaleLabel("Unweighted"))
      add("Total sample size",   n_total)
      add("Complete cases (no missing SNPs)",
          paste0(complete, " (", round(complete / n_total * 100, 1), "%)"))

      private$.writeRows(tbl, rows)
    },

    # Predicted coverageTable row count for .init(), matching the option-determined
    # row set of .fillCoverageTable exactly. meta_nonempty is the number of
    # non-empty header metadata fields (from .coverageMetaNonEmpty).
    .coverageNRows = function(has_file, meta_nonempty) {
      wmode  <- self$options$weightingMode
      run_w  <- wmode %in% c("weighted", "both") && has_file
      run_uw <- wmode %in% c("unweighted", "both") || !has_file
      meta_nonempty +
        (if (has_file) 6L else 1L) +
        1L +                                      # null-allele fixed
        1L +                                      # excluded by allele/mono QC
        (if (isTRUE(self$options$qcFilterMissing)) 1L else 0L) +
        (if (isTRUE(self$options$qcFilterHwe))     1L else 0L) +
        1L +                                      # SNPs used in score
        1L +                                      # missing strategy
        (if (run_w)  1L else 0L) +
        (if (run_uw) 1L else 0L) +
        1L +                                      # total sample size
        1L                                        # complete cases
    },

    # Non-empty header metadata field count without a full file parse.
    .coverageMetaNonEmpty = function(path) {
      if (is.null(path) || !nzchar(trimws(path %||% "")) || !file.exists(path))
        return(0L)
      raw <- tryCatch(readLines(path, n = 200, warn = FALSE),
                      error = function(e) NULL)
      if (is.null(raw)) return(0L)
      meta <- private$.parseFileMetadata(raw)
      sum(vapply(meta, function(v) nchar(trimws(v %||% "")) > 0, logical(1)))
    },

    .isAmbiguous = function(ea, oa) {
      pairs <- paste0(toupper(ea), toupper(oa))
      pairs %in% c("AT", "TA", "CG", "GC")
    },

    # Surface model-fit convergence/separation issues (from fit_diagnostics())
    # as a note on the association table, namespaced per scoring mode so the two
    # modes' notes don't overwrite each other. `issues` is empty for a clean fit.
    .setFitNote = function(score_type, model_lbl, issues) {
      key <- paste0("fitWarning_", score_type)
      if (length(issues) == 0) {
        self$results$assocTable$setNote(key, NULL, init = FALSE)
      } else {
        self$results$assocTable$setNote(key, paste0(
          "Model fit warning(s) — ", score_type, " ", model_lbl, ": ",
          paste(issues, collapse = "; "), ". Interpret the estimate with caution."),
          init = FALSE)
      }
    },

    # Group labels .fillSummaryTable will emit per mode, predicted from the raw
    # response at .init() time (before dosage/keepMask exist). Mirrors .run()'s
    # response resolution (binary numeric -> factor, caseLevel relevel) and
    # .fillSummaryTable's grouping (a row per level then Overall for a binary/
    # categorical response, else just Overall). Uses the unmasked column, so it
    # is a superset of the groups .run() produces — never fewer, so the seeded
    # rows are never short; a rare extra row (a level QC drops entirely) simply
    # stays empty and triggers a one-off rebuild.
    .predictSummaryGroups = function(resp) {
      # Stratification is an explicit user choice; when off, one Overall row/mode.
      if (!isTRUE(self$options$summaryStratify)) return("Overall")
      is_binary <- !is.null(resp) && (is.factor(resp) ||
                    length(unique(resp[!is.na(resp)])) == 2)
      if (is_binary) c(private$.respLevels(resp), "Overall") else "Overall"
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
      # Only stratify when the user asked for it (mirrors .predictSummaryGroups).
      stratify <- isTRUE(self$options$summaryStratify) && is_binary

      if (stratify) {
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

      # Dynamic per-level columns: response *levels* are data-derived. These are
      # pre-created in .init() so restore refills them (columns added only in
      # .run() are never restored); add here too, guarded, in case .init's
      # prediction missed a level (e.g. header-only .init data).
      if (do_strat) {
        for (lv in resp_levels) {
          nm <- paste0("n_lv_", make.names(lv))
          if (is.null(tryCatch(thr_tbl$getColumn(nm), error = function(e) NULL)))
            thr_tbl$addColumn(name = nm, title = as.character(lv), type = "text")
        }
      }

      show_mode_col <- length(all_scores) > 1

      thr_rows <- list()
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

          thr_rows[[length(thr_rows) + 1L]] <- row_vals
        }
      }
      private$.writeRows(thr_tbl, thr_rows)

      # ══════════════════════════════════════════════════════════════════════
      # TABLE 2: Regression  (percentileTable)
      # ══════════════════════════════════════════════════════════════════════
      cat_tbl <- self$results$percentileTable
      cat_rows <- list()

      # has_resp now includes polytomous — all three response types get a model
      has_resp <- resp_type %in% c("binary", "continuous", "polytomous")

      # Column title and notes
      ref_lbl <- cat_labels[ref_idx]

      # init=FALSE on every note so the text survives protobuf restore (see the
      # note in .fillAssocTable \u2014 this fill is gated, so it will not re-set them).
      if (!has_resp) {
        cat_tbl$getColumn("estimate")$setTitle("Estimate")
        cat_tbl$setNote("modelNote",
          "No response variable selected \u2014 counts and score ranges shown.",
          init = FALSE)
        cat_tbl$setNote("refNote", NULL, init = FALSE)
        cat_tbl$setNote("covNote", NULL, init = FALSE)
      } else if (resp_type == "binary") {
        cat_tbl$getColumn("estimate")$setTitle("OR")
        cat_tbl$setNote("modelNote", "Logistic regression (OR, 95% CI)", init = FALSE)
        cat_tbl$setNote("refNote",
          paste0("Reference category: ", ref_lbl, " (OR = 1)"), init = FALSE)
        cat_tbl$setNote("covNote",
          if (has_covs) paste0("Adjusted for: ", paste(names(covs), collapse = ", "))
          else NULL, init = FALSE)
      } else if (resp_type == "continuous") {
        cat_tbl$getColumn("estimate")$setTitle("\u03b2")
        cat_tbl$setNote("modelNote", "Linear regression (\u03b2, 95% CI)", init = FALSE)
        cat_tbl$setNote("refNote",
          paste0("Reference category: ", ref_lbl, " (\u03b2 = 0)"), init = FALSE)
        cat_tbl$setNote("covNote",
          if (has_covs) paste0("Adjusted for: ", paste(names(covs), collapse = ", "))
          else NULL, init = FALSE)
      } else {
        # polytomous
        cat_tbl$getColumn("estimate")$setTitle("OR")
        cat_tbl$setNote("modelNote",
          paste0("Polytomous logistic regression (nnet::multinom); ",
                 "ORs relative to outcome reference level: ",
                 resp_levels[1]), init = FALSE)
        cat_tbl$setNote("refNote",
          paste0("Reference score category: ", ref_lbl, " (OR = 1 per contrast)"),
          init = FALSE)
        cat_tbl$setNote("covNote",
          if (has_covs) paste0("Adjusted for: ", paste(names(covs), collapse = ", "))
          else NULL, init = FALSE)
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
          safe_rhs(names(covs)) else ""

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

            cat_rows[[length(cat_rows) + 1L]] <<- list(
                score_type  = mode_label,
                contrast    = contrast_lbl,
                category    = if (is_ref) paste0(lbl, " \u25c6") else lbl,
                n           = disp$n,
                score_range = disp$rng,
                estimate    = est,
                ci_low      = ci_lo,
                ci_high     = ci_hi,
                p           = p_v
              )
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
            cat_rows[[length(cat_rows) + 1L]] <- list(
                score_type  = mode_label,
                contrast    = "",
                category    = lbl,
                n           = disp$n,
                score_range = disp$rng,
                estimate    = '',
                ci_low      = '',
                ci_high     = '',
                p           = ''
              )
          }
        }
      }
      private$.writeRows(cat_tbl, cat_rows)
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
      # init=FALSE marks these as run-phase results so protobuf restore keeps
      # their text across an unrelated option click (an init=TRUE note is treated
      # as an init placeholder, regenerated empty by .init and never restored —
      # and this fill is gated behind do_assoc, so it does not re-set them).
      if (has_covs)
        self$results$assocTable$setNote(
          "covNote", paste0("Adjusted for: ", paste(names(covs), collapse = ", ")),
          init = FALSE)
      else
        self$results$assocTable$setNote("covNote", NULL, init = FALSE)

      resp_note <- switch(resp_type,
        binary     = paste0("Response: ", respCol,
                            " (", lvls[2], " vs ", lvls[1], ")"),
        polytomous = paste0("Response: ", respCol,
                            " (ref: ", lvls[1], "; ",
                            paste(lvls[-1], collapse = ", "), ")"),
        paste0("Response: ", respCol)
      )
      self$results$assocTable$setNote("respNote", resp_note, init = FALSE)

      cov_terms <- if (has_covs)
        safe_rhs(names(covs)) else ""

      # Accumulate rows across scoring modes; .run writes them positionally once
      # the loop is done (see .writeRows). Labels live in columns, so row order —
      # not a semantic key — identifies each row.
      add_row <- function(test, stat_label, estimate, se, ci_low, ci_high,
                          stat, df_val, p) {
        private$.assoc_acc[[length(private$.assoc_acc) + 1L]] <- list(
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
        )
      }

      # ── Binary response ───────────────────────────────────────────────────
      if (resp_type == "binary") {

        df$resp <- factor(df$resp)
        g1 <- df$pgs[df$resp == lvls[1]]
        g2 <- df$pgs[df$resp == lvls[2]]

        frm <- if (has_covs) as.formula(paste("resp ~ pgs +", cov_terms))
               else resp ~ pgs
        ff  <- with_warnings(glm(frm, data = df, family = binomial()))
        fit <- ff$value
        private$.setFitNote(score_type, "logistic regression",
                            fit_diagnostics(fit, ff$warnings))
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
        ff  <- with_warnings(
          nnet::multinom(
            if (has_covs) as.formula(paste("resp ~ pgs +", cov_terms))
            else resp ~ pgs,
            data = df, trace = FALSE))
        fit <- ff$value
        private$.setFitNote(score_type, "polytomous logistic",
                            fit_diagnostics(fit, ff$warnings))

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
      # init=FALSE so the text survives protobuf restore (see .fillAssocTable).
      if (has_other)
        self$results$interactionTable$setNote("intNote",
                paste0("Adjusted for: ",
                paste(names(other_covs), collapse = ", ")), init = FALSE)
      else
        self$results$interactionTable$setNote("intNote", NULL, init = FALSE)

      other_terms <- if (length(other_col_names) > 0)
        paste(" +", paste(other_col_names, collapse = " + "))
      else ""

      # Formulas use safe internal names; cov1 is always the column name
      frm_int_str  <- paste0("resp ~ pgs * cov1", other_terms)
      frm_main_str <- paste0("resp ~ pgs + cov1", other_terms)

      # Accumulate across scoring modes; .run writes positionally (see .writeRows).
      add_row <- function(model_lbl, term_lbl, estimate, ci_low, ci_high, p) {
        private$.inter_acc[[length(private$.inter_acc) + 1L]] <- list(
          score_type = score_type,
          model      = model_lbl,
          term       = term_lbl,
          estimate   = estimate,
          ci_low     = ci_low,
          ci_high    = ci_high,
          p          = p
        )
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
  .plotDist        = function(image, ...) plotDist(image,        private$.cache, self$options),
  .plotStrat       = function(image, ...) plotStrat(image,       private$.cache, self$options),
  .plotForest      = function(image, ...) plotForest(image,      private$.cache, self$options),
  .plotROC         = function(image, ...) plotROC(image,         private$.cache, self$options),
  .plotCalibration = function(image, ...) plotCalibration(image, private$.cache, self$options)
  )  # end private
)

