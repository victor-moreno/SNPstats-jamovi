#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype allele HWE.exact LD
#' @importFrom haplo.stats setupGeno haplo.em haplo.glm haplo.glm.control



# ══════════════════════════════════════════════════════════════════════════════
# snpStatsClass
# ══════════════════════════════════════════════════════════════════════════════


snpStatsClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "snpStatsClass",
  inherit = snpStatsBase,
  private = list(

    # ── Private state ────────────────────────────────────────────────────────
    .miss_cache    = NULL,   # populated during descriptive run; read by .plotMissingness
    .prep_key      = NULL,   # signature of prep inputs at last snp_prepare() call
    .prep_cache    = NULL,   # cached snp_prepare() result (see .run)
    .snpsum_rowkeys = NULL,  # stratified snpSummaryTable rowKeys from .init
    .covdesc_rowkeys = NULL, # covDescTable rowKeys pre-created by .init; .run
                             # setRows into them when its own keys match (.init,
                             # .load and .run share one object, so this carries)

    # ── Refresh guards ────────────────────────────────────────────────────────
    # These say whether a table still needs its ROWS written. They cannot say
    # whether it needs RECOMPUTING: jamovi rebuilds the analysis object on every
    # option click and Table$fromProtoBuf only restores cells into rows that
    # already exist by the end of .init(). Rows added by .run() via addRow() are
    # therefore never restored, so rowCount is always 0 here and .need_fill is
    # always TRUE. Use .cached() to skip recomputation — not these.
    .need_fill     = function(tbl) tbl$rowCount == 0 || tbl$isNotFilled(),

    # Full data during .init. jamovi hands .init header-only data (columns and
    # factor levels, no values), which cannot tell a numeric 0/1 response from a
    # continuous one, nor which variables have missing values — both decide the
    # row structure. Reading the values here buys an exact prediction, at the
    # cost of one dataset read per option click.
    .init_data = function() {
      d <- tryCatch(self$data, error = function(e) NULL)
      if (!is.null(d) && nrow(d) > 0) return(d)          # provided (R / tests)
      tryCatch(self$readDataset(FALSE), error = function(e) d)
    },

    # Predict covDescTable's rowKeys. Mirrors the row order of compute_cov_desc()
    # exactly, so .run can setRow into the rows .init created rather than
    # rebuilding them. Returns NULL when the structure cannot be determined.
    .covdesc_keys = function(opts, data) {
      if (is.null(data)) return(list())
      rows <- list()
      add  <- function(key, variable, level)
        rows[[length(rows) + 1L]] <<- list(key = key, variable = variable, level = level)

      resp <- opts$response
      if (!is.null(resp) && length(resp) == 1L && resp %in% names(data)) {
        col   <- data[[resp]]
        rtype <- detect_response_type(col, opts$responseType)
        if (identical(rtype, "quantitative"))
          add(paste0(resp, "|mean"), resp, "Mean ± SD")
        add(paste0(resp, "|valid"), resp, "Valid")
        if (anyNA(col)) add(paste0(resp, "|miss"), resp, "Missing")
      }

      for (v in opts$covariates) {
        if (!v %in% names(data)) next
        col <- data[[v]]
        if (is.factor(col) || is.character(col)) {
          for (lvl in levels(as.factor(col))) add(paste0(v, "|lvl|", lvl), v, lvl)
        } else {
          add(paste0(v, "|mean"), v, "Mean ± SD")
        }
        if (anyNA(col)) add(paste0(v, "|miss"), v, "Missing")
      }
      rows
    },

    # Response levels used as the stratified group columns, or NULL when not
    # stratified. compute_cov_desc derives these from the response after any
    # mask; a mask that removed an entire response level would differ, which is
    # why .write_cov_desc only adds columns that are not already present.
    .strat_levels = function(opts, data) {
      if (!isTRUE(opts$subpop) || is.null(data)) return(NULL)
      resp <- opts$response
      if (is.null(resp) || length(resp) != 1L || !resp %in% names(data)) return(NULL)
      col   <- data[[resp]]
      rtype <- detect_response_type(col, opts$responseType)
      if (!rtype %in% c("binary", "categorical")) return(NULL)
      levels(as.factor(col))
    },

    # Add an Array item unless it is already there, and return it.
    # Array$addItem() appends unconditionally — it is NOT idempotent. Calling it
    # in both .init and .run therefore produced two items per key; restore
    # indexes the saved items by name, so the LAST duplicate won — the empty one
    # .run appended — and every per-SNP table came back blank.
    .ensure_item = function(arr, key) {
      if (!key %in% unlist(arr$itemKeys)) arr$addItem(key = key)
      tryCatch(arr$get(key = key), error = function(e) NULL)
    },

    # Pre-create n empty rows with positional keys. The per-SNP tables carry
    # their labels in a column (allele / genotype / group), so position is the
    # only thing a row identifies — .run rewrites every cell anyway. That means
    # .init only has to predict the row COUNT, not reproduce the ref-allele
    # ordering. Writers reuse these rows when the count matches (see .reuse_rows)
    # and rebuild otherwise, so a wrong prediction costs a blank, never a wrong
    # number.
    .pre_rows = function(tbl, n) {
      if (tbl$rowCount == 0L && n > 0L)
        for (i in seq_len(n)) tbl$addRow(rowKey = as.character(i))
    },

    # True when the table's existing rows can be setRow'd into rather than
    # rebuilt. Rows added in .run are never restored, so rebuilding blanks.
    .reuse_rows = function(tbl, n) tbl$rowCount == n && n > 0L,

    # Stratified group columns must exist by the end of .init or their cells are
    # never restored — addColumn has the same defect as addRow.
    .pre_strat_cols = function(tbl, grp_levels, pval = FALSE) {
      if (is.null(grp_levels) || length(grp_levels) == 0L) return(invisible())
      for (j in seq_along(grp_levels)) {
        nm <- paste0("stat_g", j - 1L)
        if (is.null(tryCatch(tbl$getColumn(nm), error = function(e) NULL)))
          tbl$addColumn(name = nm, title = grp_levels[j], type = "text")
      }
      if (pval && is.null(tryCatch(tbl$getColumn("pval"), error = function(e) NULL)))
        tbl$addColumn(name = "pval", title = "P-value", type = "text")
      invisible()
    },

    # Pre-create one LD group's rows/columns from the SNP names, mirroring the
    # keys .run_ld uses (pair keys "a___b" for the pairwise table, "row_i" rows
    # and one column per SNP for the matrix). Restore then refills them, so each
    # table's per-table .need_fill gate skips recompute on an unrelated click.
    # Predicting from the selected SNPs can only over-create (validation drops
    # SNPs), which just falls back to a rebuild — never wrong data.
    .pre_ld_item = function(item, nms, opts) {
      n <- length(nms)
      if (n < 2L) return(invisible())
      if (isTRUE(opts$ldAnalysis) && item$ldTable$rowCount == 0L)
        for (pair in combn(nms, 2, simplify = FALSE))
          item$ldTable$addRow(rowKey = paste(pair, collapse = "___"))
      if (isTRUE(opts$ldMatrix)) {
        mtbl <- item$ldMatrixTable
        for (snp in nms) {
          safe_nm <- gsub("[^A-Za-z0-9_]", "_", snp)
          if (is.null(tryCatch(mtbl$getColumn(safe_nm), error = function(e) NULL)))
            mtbl$addColumn(name = safe_nm, title = snp, type = "text")
        }
        if (mtbl$rowCount == 0L)
          for (i in seq_len(n)) mtbl$addRow(rowKey = paste0("row_", i))
      }
      invisible()
    },

    # Whether a SNP column has genotypes that snp_prepare will treat as missing
    # (NA or a null-allele code). Decides the presence of the Missing row.
    .snp_any_missing = function(col) {
      if (is.null(col)) return(FALSE)
      ch <- as.character(col)
      anyNA(ch) || any(grepl(.NULL_ALLELE_PAT, ch, ignore.case = TRUE), na.rm = TRUE)
    },

    # Number of assocTable rows a SNP will produce: one block per selected
    # model. Verified against compute_assoc — the count depends only on the
    # models and the genotype count, not on covariates or response type.
    .assoc_nrows = function(opts, col) {
      n_genos <- length(private$.snp_geno_levels(col))
      if (n_genos == 0L) return(0L)
      sum(vapply(private$.get_interaction_models(opts), function(m) switch(m,
        codominant   = n_genos,   # one row per genotype
        logadditive  = 1L,        # single per-allele row
        2L),                      # dominant / recessive / overdominant
        integer(1)))
    },

    # Observed genotypes of a SNP, excluding missing and null-allele codes.
    # NOT get_snp_level_order(): that returns NULL for alphabetical levels
    # (meaning "the user set no custom order"), which says nothing about which
    # genotypes exist. Only the count is used, to size the tables.
    .snp_geno_levels = function(col) {
      if (is.null(col)) return(character(0))
      ch <- as.character(col)
      ch <- ch[!is.na(ch) & nzchar(ch) & !grepl(.NULL_ALLELE_PAT, ch, ignore.case = TRUE)]
      unique(ch)
    },

    # Predict the row count of each interaction table from options + data, so
    # .init can pre-create exactly those rows and restore refills them (no
    # rebuild flash). Exact for the common case (multiplicative interaction,
    # categorical covariate); other cases just cost a one-off redraw since the
    # writers reconcile any mismatch with row reuse. Returns a 4-element list.
    .int_nrows = function(opts, snp_col, cov_col, response_type) {
      z <- list(main = 0L, cov = 0L, geno = 0L, cross = 0L)
      if (identical(response_type, "categorical")) return(z)
      models  <- private$.get_interaction_models(opts)
      n_genos <- length(private$.snp_geno_levels(snp_col))
      if (length(models) == 0L || n_genos == 0L || is.null(cov_col)) return(z)
      is_fac <- is.factor(cov_col) || is.character(cov_col)
      n_cat  <- if (is.factor(cov_col)) nlevels(cov_col)
                else length(unique(cov_col[!is.na(cov_col)]))
      iv_d   <- if (is_fac) max(1L, n_cat - 1L) else 1L
      geno_rows <- function(m) switch(m, codominant = n_genos, logadditive = 1L, 2L)
      snp_d     <- function(m) switch(m, codominant = max(1L, n_genos - 1L), 1L)
      seps      <- function(k) max(0L, k - 1L)
      mg    <- setdiff(models, "logadditive")            # geno / cross drop log-additive
      cov_g <- if (!is_fac && n_cat > 6) 1L else n_cat   # numeric covariate collapses to "Overall"
      list(
        main  = sum(vapply(models, function(m) snp_d(m) + iv_d + snp_d(m) * iv_d, 0L)) +
                seps(length(models)),
        cov   = if (n_cat > 6) 0L else
                sum(vapply(models, function(m) n_cat * geno_rows(m), 0L)) + seps(length(models)),
        geno  = if (length(mg) == 0L) 0L else
                sum(vapply(mg, function(m) geno_rows(m) * cov_g, 0L)) + seps(length(mg)),
        cross = if (n_cat > 6 || length(mg) == 0L) 0L else
                sum(vapply(mg, function(m) n_cat * geno_rows(m), 0L)) + seps(length(mg)))
    },

    # Cache a computed result across option clicks.
    # Unlike rowCount, a table's state IS gated by its clearWith: jamovi restores
    # it only when no option in that clearWith changed. So a surviving state
    # truthfully means "nothing this table depends on changed" — return it and
    # skip the model fitting. Wrapped in list(v=) so a NULL result still caches.
    .cached = function(tbl, compute) {
      st <- tbl$state
      if (!is.null(st)) return(st$v)
      v <- compute()
      tbl$setState(list(v = v))
      v
    },

    # ══════════════════════════════════════════════════════════════════════════
    # .init  — build table skeletons from options alone (no data access)
    #
    # Jamovi calls this on every option change. Each section guards itself so
    # it only runs when the options it depends on actually changed, leaving
    # unrelated tables untouched and preserving any data already computed.
    # Nothing here touches self$data.
    # ══════════════════════════════════════════════════════════════════════════
    .init = function() {

      opts     <- self$options
      snps     <- opts$snps
      response <- opts$response
      covs     <- opts$covariates

      has_snps <- !is.null(snps)     && length(snps)     > 0

      # Link to the online tutorial (docs/ is not bundled into the installed
      # module). Shown only in the "Getting started" state (no SNPs assigned
      # yet); hidden once the guidance disappears.
      self$results$helpBanner$setVisible(!has_snps)
      self$results$helpBanner$setContent(paste0(
        "<div style=\"font-size:0.85em; color:#666; padding:2px 0 6px;\">",
        "\U0001F4D6 <a href=\"https://victor-moreno.github.io/SNPstats-jamovi/TUTORIAL.html\" ",
        "target=\"_blank\" rel=\"noopener\">SNPstats tutorial &amp; help</a></div>"))
      has_resp <- !is.null(response) && nchar(trimws(response)) > 0
      has_covs <- !is.null(covs)     && length(covs) > 0

      res <- self$results

      # ── Descriptive: covDescTable — pre-create rows so restore can refill ──
      # Rows must exist by the end of .init or jamovi never restores their cells
      # (Table$fromProtoBuf only copies into rows that already exist), which is
      # why this table used to blank on every option click. Creating them here —
      # empty, so a first run still fills them — lets the previous results
      # reappear before .run starts. Keys are remembered for .write_cov_desc.
      if (isTRUE(opts$covDesc)) {
        tbl  <- res$descGroup$covDescGroup$covDescTable
        dat  <- private$.init_data()
        keys <- private$.covdesc_keys(opts, dat)
        private$.covdesc_rowkeys <- vapply(keys, `[[`, "", "key")
        if (tbl$rowCount == 0)
          for (r in keys) tbl$addRow(rowKey = r$key)

        # Group/pval columns must exist here too: addColumn in .run has the same
        # defect as addRow — a column absent at restore never gets its cells
        # back, so the stratified columns blanked on every click.
        grp <- private$.strat_levels(opts, dat)
        if (!is.null(grp)) {
          for (j in seq_along(grp)) {
            nm <- paste0("stat_g", j - 1L)
            if (is.null(tryCatch(tbl$getColumn(nm), error = function(e) NULL)))
              tbl$addColumn(name = nm, title = grp[j], type = "text")
          }
          if (is.null(tryCatch(tbl$getColumn("pval"), error = function(e) NULL)))
            tbl$addColumn(name = "pval", title = "P-value", type = "text")
        }
      }

      # ── Descriptive: snpSummaryTable ──────────────────────────────────────
      # rows:(snps) binding in snpStats.r.yaml auto-creates one row per SNP
      # before .init is called — no addRow needed here when unstratified.
      # Hide the Missing column when completeCases is on (no missing by definition).
      if (has_snps) {
        stbl <- res$descGroup$snpSummaryTablesGroup$snpSummaryTable
        miss_col <- tryCatch(stbl$getColumn("missing"), error = function(e) NULL)
        if (!is.null(miss_col))
          miss_col$setVisible(!isTRUE(opts$completeCases))

        # Group visibility follows the subpop option alone, so set it here: a
        # column that only becomes visible in .run makes the table change shape
        # between restore and run.
        grp_col <- tryCatch(stbl$getColumn("group"), error = function(e) NULL)
        if (!is.null(grp_col)) grp_col$setVisible(isTRUE(opts$subpop))

        # Stratified needs snp x (All + response levels) rows, which the
        # rows:(snps) binding cannot express. Build them here rather than in
        # .run: rows created in .run are never restored, so the table blanked on
        # every click whenever subpop was on.
        sgrp <- private$.strat_levels(opts, private$.init_data())
        private$.snpsum_rowkeys <- NULL
        if (!is.null(sgrp)) {
          keys <- as.vector(t(outer(opts$snps, c("All", sgrp),
                                    function(s, g) paste0(s, "|", g))))
          private$.snpsum_rowkeys <- keys
          stbl$deleteRows()                       # drop the rows:(snps) rows
          for (k in keys) stbl$addRow(rowKey = k)
        }
      }

      # ── Descriptive: per-SNP tables — pre-create items, rows and columns ──
      # .run used to addItem/addRow/addColumn here, none of which survive the
      # restore, so every per-SNP table blanked on each click. Everything these
      # need is derivable now that .init reads the data: the genotype levels give
      # the row counts, the response levels the group columns.
      if (has_snps && (isTRUE(opts$allFreq) || isTRUE(opts$genoFreq) || isTRUE(opts$hweTest))) {
        dat       <- private$.init_data()
        arr       <- res$descGroup$descSnpResults
        sgrp      <- private$.strat_levels(opts, dat)
        # showMissing is ignored under completeCases (nothing is missing then),
        # so the Missing row follows the same combination .run uses.
        show_miss <- isTRUE(opts$showMissing) && !isTRUE(opts$completeCases)

        for (nm in snps) {
          item <- private$.ensure_item(arr, nm)
          if (is.null(item)) next

          col <- if (!is.null(dat) && nm %in% names(dat)) dat[[nm]] else NULL
          lv  <- private$.snp_geno_levels(col)
          n_miss_row <- if (show_miss && private$.snp_any_missing(col)) 1L else 0L

          if (isTRUE(opts$allFreq)) {
            n_alleles <- length(unique(unlist(strsplit(lv, "/", fixed = TRUE))))
            private$.pre_rows(item$allFreqTable, n_alleles + n_miss_row)
            private$.pre_strat_cols(item$allFreqTable, sgrp)
          }
          if (isTRUE(opts$genoFreq)) {
            private$.pre_rows(item$genoFreqTable, length(lv) + n_miss_row)
            private$.pre_strat_cols(item$genoFreqTable, sgrp)
          }
          if (isTRUE(opts$hweTest)) {
            # One row for all subjects, plus one per response group.
            private$.pre_rows(item$hweTable, 1L + length(sgrp))
          }
        }
      }

      # ── Descriptive: missingnessPlot visibility ────────────────────────────
      res$descGroup$missingnessPlot$setVisible(
        has_snps && isTRUE(opts$showMissingnessPlot) && !isTRUE(opts$completeCases))

      ld_ok        <- has_snps && length(snps) >= 2
      haplo_ok     <- ld_ok   && has_resp
      haplo_int_ok <- haplo_ok && has_covs

      # ── Association: pre-seed one Array item per SNP so tables appear ──────
      # immediately (empty) while .run() computes. Items must exist by the end
      # of .init or the restore cannot refill them. Always go through
      # .ensure_item: addItem() appends unconditionally, so adding the same key
      # here and again in .run silently duplicates the item.
      if (isTRUE(opts$snpAssoc) && has_snps && has_resp) {
        show_aic  <- isTRUE(opts$showAIC)
        assoc_arr <- res$assocGroup$assocSnpResults
        adat      <- private$.init_data()
        for (snp in snps) {
          private$.ensure_item(assoc_arr, snp)
          tbl <- tryCatch(assoc_arr$get(key = snp)$assocTable, error = function(e) NULL)
          if (!is.null(tbl)) {
            tbl$getColumn("AIC")$setVisible(show_aic)
            tbl$getColumn("BIC")$setVisible(show_aic)
            # Rows here too, or .run's addRow leaves nothing for restore to
            # refill and the table blanks on every click.
            acol <- if (!is.null(adat) && snp %in% names(adat)) adat[[snp]] else NULL
            private$.pre_rows(tbl, private$.assoc_nrows(opts, acol))
          }
        }
      }

      # ── Interaction: seed one SNP item per SNP and pre-create the predicted
      # rows of its stacked tables, so restore refills them before .run and the
      # tables do not rebuild on an unrelated click (see .int_nrows).
      if (isTRUE(opts$snpInteraction) && has_snps && has_resp && has_covs) {
        show_aic  <- isTRUE(opts$showAIC)
        assoc_arr <- res$assocGroup$assocSnpResults
        adat      <- private$.init_data()
        rcol      <- if (!is.null(adat) && length(opts$response) == 1L &&
                         opts$response %in% names(adat)) adat[[opts$response]] else NULL
        rtype     <- if (!is.null(rcol)) detect_response_type(rcol, opts$responseType) else NULL
        cov1      <- opts$covariates[[1]]
        cov_col   <- if (!is.null(adat) && !is.null(cov1) && cov1 %in% names(adat))
                       adat[[cov1]] else NULL
        for (snp in snps) {
          snpItem <- private$.ensure_item(assoc_arr, snp)
          if (is.null(snpItem)) next
          snpItem$interactionTable$getColumn("AIC")$setVisible(show_aic)
          snpItem$interactionTable$getColumn("BIC")$setVisible(show_aic)
          scol <- if (!is.null(adat) && snp %in% names(adat)) adat[[snp]] else NULL
          n    <- private$.int_nrows(opts, scol, cov_col, rtype)
          if (isTRUE(opts$showInteractionTable)) private$.pre_rows(snpItem$interactionTable, n$main)
          if (isTRUE(opts$showStratByCovariate)) private$.pre_rows(snpItem$stratByCovariate, n$cov)
          if (isTRUE(opts$showStratByGenotype))  private$.pre_rows(snpItem$stratByGenotype,  n$geno)
          if (isTRUE(opts$showCrossClassTable))  private$.pre_rows(snpItem$crossClassTable,  n$cross)
        }
      }

      # ── LD: hide immediately when no sub-option is active ─────────────────
      if (!isTRUE(opts$ldAnalysis) && !isTRUE(opts$ldMatrix) && !isTRUE(opts$ldPlot))
        res$ldHaploGroup$ldGroup$setVisible(FALSE)

      # Pre-seed the LD group items + their rows/columns so restore refills them
      # and the ldResults array is not recomputed/rebuilt on an unrelated click.
      if (ld_ok && (isTRUE(opts$ldAnalysis) || isTRUE(opts$ldMatrix))) {
        ldg_arr <- res$ldHaploGroup$ldGroup$ldResults
        grps    <- private$.strat_levels(opts, private$.init_data())  # NULL unless subpop
        for (g in c("Overall", if (!is.null(grps)) grps)) {
          private$.ensure_item(ldg_arr, g)
          it <- tryCatch(ldg_arr$get(key = g), error = function(e) NULL)
          if (!is.null(it)) private$.pre_ld_item(it, snps, opts)
        }
      }

      # ── Haplotype: explicit show/hide so tables appear empty immediately ───
      # (only hiding in .init is unreliable — show explicitly too)
      hg <- res$ldHaploGroup$haploGroup
      hg$haploFreqTable$setVisible(       isTRUE(opts$haploFreq)        && ld_ok)
      haplo_int_vis <- isTRUE(opts$haploInteraction) && haplo_int_ok
      # Assoc and interaction not implemented for categorical — show message instead
      is_cat_haplo_blocked <- identical(opts$responseType, "categorical")
      if (is_cat_haplo_blocked) {
        hg$haploAssocTable$setVisible(FALSE)
        hg$haploInteractionTable$setVisible(FALSE)
        hg$haploCondCovarTable$setVisible(FALSE)
        hg$haploCondHaploTable$setVisible(FALSE)
        if (isTRUE(opts$haploAssoc) || isTRUE(opts$haploInteraction)) {
          hg$haploNotImplMsg$setContent(
            "<p style='color:orange;'>Haplotype association and interaction analyses are only implemented for binary and quantitative responses.</p>")
          hg$haploNotImplMsg$setVisible(TRUE)
        }
      } else {
        hg$haploAssocTable$setVisible(      isTRUE(opts$haploAssoc)       && haplo_ok)
        hg$haploInteractionTable$setVisible(haplo_int_vis)
        hg$haploCondCovarTable$setVisible(  haplo_int_vis)
        hg$haploCondHaploTable$setVisible(  haplo_int_vis)
        hg$haploNotImplMsg$setVisible(FALSE)
      }

      # ── Hide sections when required variables are absent ───────────────────
      if (!has_snps) {
        res$descGroup$snpSummaryTablesGroup$setVisible(FALSE)
        res$descGroup$descSnpResults$setVisible(FALSE)
      }
      if (!has_resp && !has_covs)
        res$descGroup$covDescGroup$covDescTable$setVisible(FALSE)
      if (!has_snps || !has_resp)
        res$assocGroup$assocSnpResults$setVisible(FALSE)
      if (!ld_ok)
        res$ldHaploGroup$ldGroup$setVisible(FALSE)
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
      #
      # snp_prepare() parses every SNP genotype string and is the most expensive
      # step; cache it keyed on the inputs it depends on so a pure display-toggle
      # (which fills no table) does not re-parse. nrow guards against row filters.
      prep_key <- list(
        snps          = opts$snps,
        response      = opts$response,
        covariates    = opts$covariates,
        responseType  = opts$responseType %||% "auto",
        completeCases = isTRUE(opts$completeCases),
        n             = nrow(self$data)
      )
      if (identical(prep_key, private$.prep_key) && !is.null(private$.prep_cache)) {
        prep <- private$.prep_cache
      } else {
        prep <- snp_prepare(
          data           = self$data,
          snps           = opts$snps,
          response       = opts$response,
          covariates     = opts$covariates,
          response_type  = opts$responseType %||% "auto",
          # completeCases forces joint SNP-complete-case mask globally;
          # rmSnpMissing only affects covDesc (handled in .run_descriptive)
          rm_snp_missing = isTRUE(opts$completeCases)
        )
        private$.prep_key   <- prep_key
        private$.prep_cache <- prep
      }

      # Validation messages
      if (nchar(prep$warnings) > 0) {
        self$results$validationMsg$setContent(prep$warnings)
        self$results$validationMsg$setVisible(TRUE)
      } else {
        self$results$validationMsg$setVisible(FALSE)
      }

      # ── Descriptive ───────────────────────────────────────────────────────
      any_desc <- isTRUE(opts$covDesc) || isTRUE(opts$snpSummary) || isTRUE(opts$allFreq) || isTRUE(opts$genoFreq) || isTRUE(opts$hweTest) ||
          isTRUE(opts$subpop) || isTRUE(opts$showMissingnessPlot)
      if (any_desc) 
        private$.run_descriptive(prep, opts)
      
     
      # ── Association ───────────────────────────────────────────────────────
      any_assoc <- isTRUE(opts$snpAssoc) || isTRUE(opts$snpInteraction)
      if (any_assoc) 
        private$.run_association(prep, opts)
      

      # ── LD / Haplotype ────────────────────────────────────────────────────
      any_ld <- isTRUE(opts$ldAnalysis) || isTRUE(opts$ldMatrix) || isTRUE(opts$ldPlot) ||
                isTRUE(opts$haploFreq)  || isTRUE(opts$haploAssoc) || isTRUE(opts$haploInteraction)
      if (any_ld) {
        if (length(prep$snp_vars) < 2) {
          self$results$validationMsg$setContent(
            "<p style='color:red;'>LD and haplotype analyses require at least 2 SNPs.</p>")
          self$results$validationMsg$setVisible(TRUE)
          # Options may be on but there are not enough validated SNPs — hide the
          # tables so empty frames don't show. .init already hides when snps is
          # absent from options; this catches the post-validation case. The LD
          # tables live inside the ldResults array (not direct children of
          # ldGroup), so hide the whole ldGroup rather than reaching for
          # non-existent ldGroup$ldTable etc.
          self$results$ldHaploGroup$ldGroup$setVisible(FALSE)
          hg  <- self$results$ldHaploGroup$haploGroup
          if (isTRUE(opts$haploFreq))        hg$haploFreqTable$setVisible(FALSE)
          if (isTRUE(opts$haploAssoc))       hg$haploAssocTable$setVisible(FALSE)
          if (isTRUE(opts$haploInteraction)) {
            hg$haploInteractionTable$setVisible(FALSE)
            hg$haploCondCovarTable$setVisible(FALSE)
            hg$haploCondHaploTable$setVisible(FALSE)
          }
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
      # completeCases suppresses all missing-value display rows
      run_showMissing <- isTRUE(opts$showMissing) && !isTRUE(opts$completeCases)
      run_showMissingnessPlot <- isTRUE(opts$showMissingnessPlot)

      private$.miss_cache <- list()

      res <- self$results$descGroup

      # ── Covariate descriptives ─────────────────────────────────────────────
      if (run_covDesc) {

        # cc_mask for covDesc:
        # - completeCases: use global complete_mask (already includes joint SNP filter)
        # - rmSnpMissing: compute joint-SNP mask locally, only for covDesc; does not
        #   change the global complete_mask used by SNP frequency/association analyses
        cc_mask <- if (isTRUE(opts$completeCases)) {
          prep$complete_mask
        } else if (isTRUE(opts$rmSnpMissing) && length(prep$snp_vars) > 0) {
          all_snps_ok <- Reduce(`&`, lapply(prep$snp_vars,
                                function(nm) !is.na(prep$snp_data[[nm]]$clean)))
          prep$complete_mask & all_snps_ok
        } else {
          NULL
        }

        tbl <- res$covDescGroup$covDescTable
        if (!is.null(prep$cov_df) || !is.null(prep$response_raw)) {
          # Recompute only when covDescTable's clearWith fired; otherwise reuse
          # the cached descriptives. The rows themselves are data-derived, so
          # they must be rewritten either way — but from cached values.
          result <- private$.cached(tbl, function()
            compute_cov_desc(prep, subpop = run_subpop, mask = cc_mask))
          # Always write: .init's rows start empty, so on a first run there is
          # nothing to restore. Writing cached values into the existing rows is
          # cheap and produces exactly what restore already put there, so it
          # costs nothing visually — it is the recompute above that is skipped.
          if (!is.null(result))
            private$.write_cov_desc(tbl, result, prep)
        }
        tbl$setVisible(tbl$rowCount > 0)
      }

      # ── SNP summary table ──────────────────────────────────────────────────
      if (! length(prep$snp_vars) > 0) return() # no SNPs, nothing to do here

      if (run_snpSummary) {
        tbl <- res$snpSummaryTablesGroup$snpSummaryTable
        tbl$getColumn("group")$setVisible(run_subpop)

        if (!run_subpop) {
          # rows:(snps) binding pre-creates one row per SNP.
          # Only compute rows that are not yet filled — new SNPs get filled
          # without disturbing existing results.
          for (snp in prep$snp_vars) {
            if (tryCatch(!tbl$isFilled(rowKey = snp), error = function(e) TRUE)) {
              sp     <- prep; sp$snp_vars <- snp
              result <- compute_snp_summary(sp, subpop = FALSE)
              if (!is.null(result) && nrow(result) > 0)
                private$.fill_snp_summary_row(tbl, result[1L, ])
            }
          }
        } else {
          # Stratified needs snp x (All + groups) rows, which .init pre-creates
          # (the rows:(snps) binding cannot express them). Recompute only when
          # the table's clearWith fired; writing cached values into the existing
          # rows is cheap and matches what restore already put there.
          result <- private$.cached(tbl, function() compute_snp_summary(prep, subpop = TRUE))
          if (!is.null(result) && nrow(result) > 0)
            private$.write_snp_summary_strat(tbl, result)
        }

        # Null-allele note
        total_null_across_snps <- sum(sapply(prep$snp_vars, function(nm) {
          raw <- as.character(prep$data[[nm]])
          sum(!is.na(raw) & grepl(.NULL_ALLELE_PAT, raw, ignore.case = TRUE))
        }))
        tbl$setNote(key = "null_allele",
          note = if (total_null_across_snps > 0)
            paste0(total_null_across_snps,
                   " genotype(s) coded as 0/0 were treated as missing (NA).")
          else NULL)

        if (isTRUE(opts$completeCases)) {
          n_analyzed <- sum(prep$complete_mask)
          n_removed  <- prep$n_rows - n_analyzed
          tbl$setNote(key = "missing_resp_cov",
            note = sprintf(
              "Complete-case analysis: %d of %d observations used (%d removed due to missing SNP, response, or covariate).",
              n_analyzed, prep$n_rows, n_removed))
        } else {
          parts <- c(if (!is.null(prep$cov_df) && ncol(prep$cov_df) > 0) "covariates",
                     if (!is.null(prep$response_raw)) "response")
          if (length(parts) > 0)
            tbl$setNote(key = "missing_resp_cov",
              note = paste0("Rows missing any ",
                            paste(parts, collapse = " or "), " or SNP value are excluded."))
        }
      }

      # ── Per-SNP descriptives ───────────────────────────────────────────────
      if (run_allFreq || run_genoFreq || run_hweTest || run_showMissingnessPlot) {
        arr <- res$descSnpResults

        # .init already seeded one item per requested SNP; add only those it
        # could not know about (a SNP that passes validation but was absent from
        # the data at .init). .ensure_item, never addItem: a duplicate key makes
        # restore pick the empty item and the whole SNP comes back blank.
        for (nm in prep$snp_vars) {
          sd   <- prep$snp_data[[nm]]
          item <- private$.ensure_item(arr, nm)
          if (is.null(item)) next

          n_typed          <- sd$n_typed
          n_total_eligible <- sum(prep$complete_mask)
          total_missing    <- n_total_eligible - n_typed

          # Typing-rate HTML
          typing_html <- sprintf("<b>Typed samples:</b> %d / %d (%.1f%%)",
            n_typed, n_total_eligible,
            if (n_total_eligible > 0) n_typed / n_total_eligible * 100 else 0)
          if (total_missing > 0 && !isTRUE(opts$completeCases))
            typing_html <- paste0(typing_html, sprintf(
              " ── <b>Missing SNP:</b> %d (%.1f%%)",
              total_missing,
              if (n_total_eligible > 0) total_missing / n_total_eligible * 100 else 0))
          if (isTRUE(opts$completeCases)) {
            n_removed <- prep$n_rows - n_total_eligible
            if (n_removed > 0)
              typing_html <- paste0(typing_html, sprintf(
                " ── <b>Complete-case analysis:</b> %d removed from dataset of %d.",
                n_removed, prep$n_rows))
          }
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

          # Allele frequency — recompute only when the table's clearWith fired;
          # .init pre-created the rows, so writing cached values into them
          # reproduces what restore already showed.
          if (run_allFreq) {
            result <- private$.cached(item$allFreqTable, function()
              compute_allele_freq(nm, prep, subpop = run_subpop,
                                  show_missing = run_showMissing))
            if (!is.null(result))
              private$.write_allele_freq(item$allFreqTable, result, run_subpop, prep)
          }

          # Genotype frequency
          if (run_genoFreq) {
            result <- private$.cached(item$genoFreqTable, function()
              compute_geno_freq(nm, prep, subpop = run_subpop,
                                show_missing = run_showMissing))
            if (!is.null(result))
              private$.write_geno_freq(item$genoFreqTable, result, run_subpop,
                                       prep, prep$response_type)
          }

          # HWE
          if (run_hweTest) {
            result <- private$.cached(item$hweTable, function()
              compute_hwe(nm, prep, subpop = run_subpop,
                          show_missing = run_showMissing))
            if (!is.null(result))
              private$.write_hwe(item$hweTable, result)
          }
        }
      }

      # Missingness plot: hide when completeCases is on (no missingness by definition)
      if (run_showMissingnessPlot) {
        res$missingnessPlot$setVisible(
          !isTRUE(opts$completeCases) && length(private$.miss_cache) > 0)
      }

    },

    # ══════════════════════════════════════════════════════════════════════════
    # Table writers — thin jamovi-specific adapters over compute_* results
    # ══════════════════════════════════════════════════════════════════════════

    .write_cov_desc = function(tbl, result, prep) {
      tbl_df     <- result$table
      grp_levels <- result$grp_levels   # NULL when not stratified
      do_strat   <- !is.null(grp_levels) && length(grp_levels) > 0
      keys       <- as.character(tbl_df[["key"]])

      # Reuse the rows .init created when they match, so the table is updated in
      # place and never structurally churns. They can only mismatch when .init's
      # header-only guess of the response type was wrong; rebuild from scratch
      # then — correct, at the cost of the blank this change exists to avoid.
      reuse <- identical(keys, private$.covdesc_rowkeys)
      if (!reuse) {
        tbl$deleteRows()
        private$.covdesc_rowkeys <- keys
      }

      # Group columns and pval are data-derived — all added dynamically so
      # that pval stays last. clearWith resets columns between runs.
      if (do_strat) {
        for (j in seq_along(grp_levels)) {
          nm <- paste0("stat_g", j - 1L)
          if (is.null(tryCatch(tbl$getColumn(nm), error = function(e) NULL)))
            tbl$addColumn(name = nm, title = grp_levels[j], type = "text")
        }
        if (is.null(tryCatch(tbl$getColumn("pval"), error = function(e) NULL)))
          tbl$addColumn(name = "pval", title = "P-value", type = "text")
      }

      for (i in seq_len(nrow(tbl_df))) {
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
        if (reuse) tbl$setRow(rowKey = keys[i], values = row)
        else       tbl$addRow(rowKey = keys[i], values = row)
      }

      if (length(result$notes) > 0)
        tbl$setNote(key = "cov_desc_note",
                    note = paste(result$notes, collapse = " "))
    },

    # Non-stratified: table has rows:(snps) pre-binding; fill the single row
    # keyed by SNP name. content:($key) populates the snp column automatically.
    .fill_snp_summary_row = function(tbl, row) {
      tbl$setRow(rowKey = as.character(row[["snp"]]), values = list(
        alleles    = as.character(row[["alleles"]]),
        n          = as.integer(row[["n"]]),
        missing    = if (is.na(row[["missing"]])) NA_integer_
                     else as.integer(row[["missing"]]),
        maf        = as.numeric(row[["maf"]]),
        genoCounts = as.character(row[["geno_counts"]]),
        hwePval    = as.numeric(row[["hwe_pval"]])))
    },

    # Stratified: multiple rows per SNP — deleteRows() then addRow() loop.
    .write_snp_summary_strat = function(tbl, result) {
      # Rows are pre-created in .init (snp|group keys) so restore can refill
      # them; setRow into those when they match, else rebuild.
      keys  <- paste0(as.character(result[["snp"]]), "|", as.character(result[["group"]]))
      reuse <- identical(keys, private$.snpsum_rowkeys)
      if (!reuse) {
        tbl$deleteRows()
        private$.snpsum_rowkeys <- keys
      }
      for (i in seq_len(nrow(result)))
        (if (reuse) tbl$setRow else tbl$addRow)(rowKey = keys[i], values = list(
          snp        = as.character(result[["snp"]][i]),
          alleles    = as.character(result[["alleles"]][i]),
          group      = as.character(result[["group"]][i]),
          n          = as.integer(result[["n"]][i]),
          missing    = if (is.na(result[["missing"]][i])) NA_integer_
                       else as.integer(result[["missing"]][i]),
          maf        = as.numeric(result[["maf"]][i]),
          genoCounts = as.character(result[["geno_counts"]][i]),
          hwePval    = as.numeric(result[["hwe_pval"]][i])))
    },

    .write_allele_freq = function(tbl, result, do_strat, prep) {
      reuse <- private$.reuse_rows(tbl, nrow(result))
      if (!reuse) tbl$deleteRows()
      grp_levels <- if (do_strat && !is.null(prep$response_raw))
                      levels(as.factor(prep$response_raw)) else NULL
      # .init adds these too; only add what is missing, so a restored column
      # keeps its cells.
      private$.pre_strat_cols(tbl, grp_levels)
      for (i in seq_len(nrow(result))) {
        row <- list(
          allele = as.character(result[["allele"]][i]),
          stat   = as.character(result[["overall"]][i]))
        if (do_strat && !is.null(grp_levels))
          for (j in seq_along(grp_levels))
            row[[paste0("stat_g", j-1L)]] <-
              as.character(result[[paste0("stat_g", j-1L)]][i] %||% "")
        (if (reuse) tbl$setRow else tbl$addRow)(rowKey = as.character(i), values = row)
      }
    },

    .write_geno_freq = function(tbl, result, do_strat, prep, response_type) {
      reuse <- private$.reuse_rows(tbl, nrow(result))
      if (!reuse) tbl$deleteRows()
      grp_levels <- if (do_strat && !is.null(prep$response_raw))
                      levels(as.factor(prep$response_raw)) else NULL
      if (response_type == "quantitative")
        tbl$getColumn("responseStat")$setVisible(TRUE)
      private$.pre_strat_cols(tbl, grp_levels)
      for (i in seq_len(nrow(result))) {
        row <- list(
          genotype     = as.character(result[["genotype"]][i]),
          stat         = as.character(result[["overall"]][i]),
          responseStat = as.character(result[["response_stat"]][i] %||% ""))
        if (do_strat && !is.null(grp_levels))
          for (j in seq_along(grp_levels))
            row[[paste0("stat_g", j-1L)]] <-
              as.character(result[[paste0("stat_g", j-1L)]][i] %||% "")
        (if (reuse) tbl$setRow else tbl$addRow)(rowKey = as.character(i), values = row)
      }
    },

    .write_hwe = function(tbl, result) {
      reuse <- private$.reuse_rows(tbl, nrow(result$rows))
      if (!reuse) tbl$deleteRows()
      labels <- result$col_labels
      if (length(labels) == 3) {
        tbl$getColumn("n11")$setTitle(labels[1])
        tbl$getColumn("n12")$setTitle(labels[2])
        tbl$getColumn("n22")$setTitle(labels[3])
      }
      rows <- result$rows
      for (i in seq_len(nrow(rows))) {
        (if (reuse) tbl$setRow else tbl$addRow)(rowKey = as.character(i), values = list(
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
      # Interaction is only implemented for binary responses
      is_cat_blocked <- run_snpInteraction && prep$response_type == "categorical"
      if (is_cat_blocked) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>Interaction analyses are only implemented for binary response.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_snpInteraction <- FALSE
      }

      int_models   <- private$.get_interaction_models(opts)
      assoc_arr    <- self$results$assocGroup$assocSnpResults

      # Association array items: count depends on validated SNPs (data-derived),
      # so addItem() must stay in .run() (Level 3). Association group heading
      # only appears once content is actually being written.
      for (nm in prep$snp_vars) {
        private$.ensure_item(assoc_arr, nm)
        # Categorical response: interaction is not supported. The global
        # validationMsg already explains this; add a note by the empty table.
        if (is_cat_blocked && length(int_models) > 0) {
          assoc_arr$get(key = nm)$interactionTable$setNote(key = "catBlocked",
            note = "Interaction analyses are only implemented for binary response.")
        }
      }

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

        # Main association table — recompute only when its clearWith fired.
        # .init pre-created the rows, so writing cached values into them
        # reproduces exactly what restore already showed.
        if (run_snpAssoc) {
          result <- private$.cached(item$assocTable, function()
            compute_assoc(nm, prep,
                          models   = int_models,
                          ci_width = opts$ciWidth %||% 95))
          if (!is.null(result$rows) && nrow(result$rows) > 0)
            private$.write_assoc(item$assocTable, result, opts, nm)
        }

        # Interaction tables
        if (run_snpInteraction && !is.null(prep$cov_df) && ncol(prep$cov_df) > 0) {
          interaction_var <- names(prep$cov_df)[1]  # first covariate as interaction var
          response_cc <- prep$response_enc[sd$snp_mask]
          cov_df_cc   <- prep$cov_df[sd$snp_mask, , drop = FALSE]
          ref         <- sd$ref
          user_levels <- sd$user_levels
          # Ref-allele-first orientation (as compute_assoc uses), not raw
          # clean_cc: genotype counts in .compute_stats match display labels
          # only when het alleles are ordered like user_levels. Passing clean_cc
          # made the heterozygous row show 0 in stratified / cross-class tables.
          snp_geno_cc <- as.character(sd$geno_cc)

          # One table per SNP, all models stacked (blank separator between them).
          # The writers loop int_models and stack with a model column. They are
          # called unconditionally (not gated by .need_fill): each caches its fit
          # in the table state and renders with row reuse, so an unrelated click
          # neither refits nor rebuilds — .init pre-created the rows for restore.
          if (isTRUE(opts$showInteractionTable))
            private$.fill_interaction(
              item$interactionTable, snp_geno_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              prep$response_type, opts, int_models, user_levels,
              prep$response_raw[sd$snp_mask], nm, separators = TRUE)
          if (isTRUE(opts$showStratByCovariate))
            private$.fill_strat_by_covariate(
              item$stratByCovariate, snp_geno_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              prep$response_type, opts, int_models, user_levels,
              prep$response_raw[sd$snp_mask], nm)
          if (isTRUE(opts$showStratByGenotype))
            private$.fill_strat_by_genotype(
              item$stratByGenotype, snp_geno_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              prep$response_type, opts, int_models, user_levels,
              prep$response_raw[sd$snp_mask], nm)
          if (isTRUE(opts$showCrossClassTable))
            private$.fill_cross_class(
              item$crossClassTable, snp_geno_cc, ref,
              response_cc, cov_df_cc, interaction_var,
              prep$response_type, opts, int_models, user_levels,
              prep$response_raw[sd$snp_mask], nm)
        }
      }
    },

    # Write compute_assoc result to a jamovi assocTable
    .write_assoc = function(tbl, result, opts, snp_nm) {
      reuse <- private$.reuse_rows(tbl, nrow(result$rows))
      if (!reuse) tbl$deleteRows()
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
        (if (reuse) tbl$setRow else tbl$addRow)(rowKey = as.character(i), values = list(
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
    # .run_ldhaplo  — runs LD (pairwise, direct) and haplotype analyses
    # ══════════════════════════════════════════════════════════════════════════
    .run_ldhaplo = function(prep, opts) {
      run_ldAnalysis       <- isTRUE(opts$ldAnalysis)
      run_ldMatrix         <- isTRUE(opts$ldMatrix)
      run_ldPlot           <- isTRUE(opts$ldPlot)
      run_haploFreq        <- isTRUE(opts$haploFreq)
      run_haploAssoc       <- isTRUE(opts$haploAssoc)
      run_haploInteraction <- isTRUE(opts$haploInteraction)
      run_subpop           <- isTRUE(opts$subpop) &&
                              !is.null(prep$response_raw) &&
                              prep$response_type %in% c("binary", "categorical")

      # Assoc/interaction not implemented for categorical (auto-detected case)
      hg <- self$results$ldHaploGroup$haploGroup
      if (prep$response_type == "categorical" && (run_haploAssoc || run_haploInteraction)) {
        hg$haploAssocTable$setVisible(FALSE)
        hg$haploInteractionTable$setVisible(FALSE)
        hg$haploCondCovarTable$setVisible(FALSE)
        hg$haploCondHaploTable$setVisible(FALSE)
        hg$haploNotImplMsg$setContent(
          "<p style='color:orange;'>Haplotype association and interaction analyses are only implemented for binary and quantitative responses.</p>")
        hg$haploNotImplMsg$setVisible(TRUE)
        run_haploAssoc       <- FALSE
        run_haploInteraction <- FALSE
      }

      # ── LD ────────────────────────────────────────────────────────────────
      # Each pair uses only individuals where BOTH of that pair's SNPs are
      # non-missing (pairwise complete cases), so a missing third SNP never
      # shrinks the sample for an unrelated pair.
      # completeCases additionally requires response/covariate non-missingness.
      ldg_arr <- self$results$ldHaploGroup$ldGroup$ldResults
      if (run_ldAnalysis || run_ldMatrix || run_ldPlot) {
        complete_cases <- isTRUE(opts$completeCases)
        masks <- list(Overall = NULL)
        if (run_subpop) {
          grp_levels <- if (is.factor(prep$response_raw)) levels(prep$response_raw)
                        else sort(unique(na.omit(as.character(prep$response_raw))))
          for (grp in grp_levels)
            masks[[grp]] <- !is.na(prep$response_raw) & as.character(prep$response_raw) == grp
        }
        for (g in names(masks)) {
          private$.ensure_item(ldg_arr, g)
          item <- ldg_arr$get(key = g)
          # Fill each output only when it is shown AND still needs filling, so
          # toggling one (e.g. the matrix) never recomputes/rebuilds a sibling
          # that is already computed. .run_ld computes the shared LD store once
          # if any of the three needs it.
          need_table  <- run_ldAnalysis && private$.need_fill(item$ldTable)
          need_matrix <- run_ldMatrix   && private$.need_fill(item$ldMatrixTable)
          need_plot   <- run_ldPlot     && is.null(item$ldPlotImage$state)
          private$.run_ld(item, prep, opts, complete_cases,
                          need_table, need_matrix, need_plot, group_mask = masks[[g]])
        }
      }

      # ── Haplotype ─────────────────────────────────────────────────────────
      # haplo.em / haplo.glm estimate missing genotypes via EM, so by default we
      # keep every subject typed at >= 1 SNP and let the EM handle the rest;
      # for association/interaction, na.geno.keep then drops rows missing the
      # response or covariates. completeCases overrides this with a strict
      # all-SNP complete-case mask.
      if (run_haploFreq || run_haploAssoc || run_haploInteraction) {
        jg         <- .make_haplo_geno_list(prep, complete_case = isTRUE(opts$completeCases))
        geno_list  <- jg$geno_list
        hap_mask   <- jg$mask
        response_jm     <- if (!is.null(prep$response_enc)) prep$response_enc[hap_mask] else NULL
        response_raw_jm <- if (!is.null(prep$response_raw)) prep$response_raw[hap_mask] else NULL
        cov_df_jm       <- if (!is.null(prep$cov_df) && nrow(prep$cov_df) > 0)
                             prep$cov_df[hap_mask, , drop = FALSE] else prep$cov_df
        private$.run_haplo(
          geno_list     = geno_list,
          response      = response_jm,
          response_raw  = response_raw_jm,
          response_type = prep$response_type,
          cov_df        = cov_df_jm,
          opts          = opts,
          run_haploFreq        = run_haploFreq,
          run_haploAssoc       = run_haploAssoc,
          run_haploInteraction = run_haploInteraction,
          run_subpop           = run_subpop
        )
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

    .fill_interaction = function(tbl, snp_raw, ref, response, cov_df,
                                  interaction_var, response_type, opts,
                                  int_models, user_levels = NULL, response_raw, snp_lbl,
                                  separators = FALSE) {
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "β")
      adj_vars <- setdiff(names(cov_df), interaction_var)
      int_lbl  <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      int_type <- if (is.null(opts$interactionType)) "multiplicative" else opts$interactionType
      formula_token <- switch(int_type,
        multiplicative       = paste0(snp_lbl, " × ", int_lbl),
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
      model_labels <- .MODEL_LABELS
      .label_term <- function(term, mdl, geno_labels) {
        lbl <- gsub("^snp([^:]+)", paste0(snp_lbl, "(\\1)"), term)
        lbl <- gsub(":snp([^:]+)", paste0(":", snp_lbl, "(\\1)"), lbl)
        lbl
      }
      # Fitting the per-model interaction models is the expensive step; cache the
      # rows it produces so an unrelated click reuses them (see .fill_strat_*).
      built <- private$.cached(tbl, function() {
        out <- list()
        any_clamped <- FALSE
        wrote_any   <- FALSE
        inter_pvals <- character(0)
        fit_diags   <- character(0)
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
          if (!is.null(p_inter) && !is.na(p_inter))
            inter_pvals[[model_labels[[mdl]]]] <- fmt_pval(p_inter)
          dg <- attr(res_list, "diagnostics")
          if (length(dg) > 0)
            fit_diags <- c(fit_diags, paste0(model_labels[[mdl]], ": ", paste(dg, collapse = "; ")))
          snp_char_l    <- as.character(snp_raw)
          all_genos_l   <- .all_genos_for_snp(user_levels, snp_char_l, ref)
          geno_labels_l <- private$.geno_labels_for_model(mdl, all_genos_l, ref)
          first_row   <- TRUE
          for (res in res_list) {
            rtype <- if (is.null(res$row_type)) "snp" else res$row_type
            if (rtype == "adjustment" && !show_adj) next
            if (separators && first_row && wrote_any)
              out[[length(out) + 1L]] <- list(
                model = "", term = "", effect = "", ciLow = "", ciHigh = "",
                pval = "", AIC = "", BIC = "")
            if (is.na(res$effect)) any_clamped <- TRUE
            term_label <- .label_term(res$term, mdl, geno_labels_l)
            aic_val <- if (first_row && !is.nan(res$aic)) fmt3(res$aic) else ""
            bic_val <- if (first_row && !is.null(res$bic) && !is.nan(res$bic))
              fmt3(res$bic) else ""
            out[[length(out) + 1L]] <- list(
              model = if (first_row) model_labels[mdl] else "",
              term  = term_label,
              effect = fmt3(res$effect), ciLow = fmt3(res$ci_low), ciHigh = fmt3(res$ci_high),
              pval   = fmt_pval(res$pval),
              AIC = aic_val, BIC = bic_val)
            first_row <- FALSE
            wrote_any <- TRUE
          }
        }
        list(rows = out, pvals = inter_pvals, any_clamped = any_clamped,
             fit_diags = fit_diags)
      })
      if (length(built$pvals) == 1L)
        tbl$setNote(key = "interactionPval",
          note = paste0("Interaction p-value (LRT): ", built$pvals[[1L]]))
      else if (length(built$pvals) > 1L)
        tbl$setNote(key = "interactionPval",
          note = paste0("Interaction p-value (LRT) — ",
                        paste(names(built$pvals), built$pvals, sep = ": ", collapse = "; ")))
      if (isTRUE(built$any_clamped))
        tbl$setNote(key = "separation",
                    note = "One or more OR/CI suppressed (shown as blank) due to complete or quasi-complete separation.")
      if (length(built$fit_diags) > 0)
        tbl$setNote(key = "fitWarning",
                    note = paste0("Model fit warning(s) — ",
                                  paste(built$fit_diags, collapse = " | "),
                                  ". Interpret the affected estimates with caution."))
      private$.render_flat(tbl, built$rows)
    },
    # Flat-row BUILDER for the stratified tables: collects rows into a list
    # instead of writing to the table, so the (expensive) fit that produces them
    # can be cached and the list rendered with row reuse. model+grp1 are shown
    # once per block (combineBelow blanks repeats) with a blank separator row
    # between models. Returns list(emit = fn(model,grp1,grp2,vals), rows = fn()).
    .mk_flat_rows = function() {
      out <- list(); cur_model <- NULL; cur_grp1 <- NULL
      blank <- list(stat0 = "", stat1 = "", effect = "", ciLow = "", ciHigh = "", pval = "")
      emit <- function(model_lbl, grp1_lbl, grp2, vals) {
        if (!is.null(cur_model) && !identical(model_lbl, cur_model))
          out[[length(out) + 1L]] <<- c(list(model = "", grp1 = "", grp2 = ""), blank)
        show_m <- is.null(cur_model) || !identical(model_lbl, cur_model)
        show_g <- show_m || is.null(cur_grp1) || !identical(grp1_lbl, cur_grp1)
        out[[length(out) + 1L]] <<- c(list(
          model = if (show_m) model_lbl else "",
          grp1  = if (show_g) grp1_lbl else "",
          grp2  = grp2), vals)
        cur_model <<- model_lbl; cur_grp1 <<- grp1_lbl
      }
      list(emit = emit, rows = function() out)
    },

    # Render a pre-built list of row value-lists into a table, reusing existing
    # rows (setRow) when the count matches so an unrelated option click updates
    # cells in place instead of a delete/rebuild that flashes. Row keys are the
    # positional "1".."N" that .init pre-creates, so restore refills them. Works
    # for any column set (missing columns are left untouched by setRow/addRow).
    .render_flat = function(tbl, rows) {
      reuse <- private$.reuse_rows(tbl, length(rows))
      if (!reuse) tbl$deleteRows()
      for (i in seq_along(rows)) {
        if (reuse) tbl$setRow(rowNo = i, values = rows[[i]])
        else       tbl$addRow(rowKey = as.character(i), values = rows[[i]])
      }
      if (length(rows) == 0L && tbl$rowCount > 0L) tbl$deleteRows()
      invisible()
    },

    # One table per view, all models stacked: flat output (Model + Stratum
    # columns), blank separator row between models. The (expensive) per-model
    # fit that produces the rows is wrapped in .cached, so an unrelated option
    # click reuses the stored rows instead of refitting; .render_flat then
    # updates cells in place when the row count is unchanged (no rebuild flash).
    .fill_strat_by_covariate = function(tbl, snp_raw, ref, response, cov_df,
                                                 interaction_var, response_type, opts,
                                                 int_models, user_levels = NULL, response_raw, snp_lbl) {
      int_var_data <- cov_df[[interaction_var]]
      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      snp_char     <- as.character(snp_raw)
      all_genos    <- .all_genos_for_snp(user_levels, snp_char, ref)
      model_labels <- .MODEL_LABELS
      adj_vars     <- setdiff(names(cov_df), interaction_var)
      adj_cov_df   <- if (length(adj_vars) > 0) cov_df[, adj_vars, drop = FALSE] else NULL
      cov_levels   <- if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(int_var_data[!is.na(int_var_data)]))
      tbl$getColumn("grp1")$setTitle(int_lbl)
      tbl$getColumn("grp2")$setTitle(snp_lbl)
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "β")
      if (response_type == "binary") {
        resp_lv <- levels(as.factor(response_raw))
        tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
        tbl$getColumn("stat0")$setVisible(TRUE);     tbl$getColumn("stat1")$setVisible(TRUE)
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }
      built <- private$.cached(tbl, function() {
        if (length(table(int_var_data)) > 6) return(list(rows = list(), pvals = character(0)))
        b <- private$.mk_flat_rows(); emit <- b$emit
        inter_pvals <- character(0)
        for (mdl in int_models) {
          snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
          res_list <- fit_interaction_model(snp_enc, response, cov_df, interaction_var, mdl,
                                            response_type, opts$ciWidth, conditional = TRUE)
          if (is.null(res_list)) next
          p_i <- attr(res_list, "pval_interaction")
          if (!is.null(p_i) && !is.na(p_i)) inter_pvals[[model_labels[[mdl]]]] <- fmt_pval(p_i)
          inter_only  <- res_list[sapply(res_list, function(r)
            is.null(r$row_type) || r$row_type == "interaction")]
          geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
          for (cl in cov_levels) {
            cl_label  <- as.character(cl)
            level_res <- inter_only[grepl(cl_label, sapply(inter_only, `[[`, "term"), fixed = TRUE)]
            if (length(level_res) == 0) next
            mask_k <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
            if (!is.null(adj_cov_df) && ncol(adj_cov_df) > 0) mask_k <- mask_k & complete.cases(adj_cov_df)
            st       <- private$.compute_stats(geno_labels, snp_char[mask_k], response[mask_k], response_type)
            grp1_lbl <- paste0(int_lbl, ": ", cl_label)
            if (mdl == "logadditive") {
              emit(model_labels[[mdl]], grp1_lbl, "per allele", list(
                stat0 = "", stat1 = "",
                effect = fmt3(level_res[[1]]$effect), ciLow = fmt3(level_res[[1]]$ci_low),
                ciHigh = fmt3(level_res[[1]]$ci_high), pval = fmt_pval(level_res[[1]]$pval)))
              next
            }
            emit(model_labels[[mdl]], grp1_lbl, geno_labels[1], list(
              stat0 = st$s0[1], stat1 = st$s1[1],
              effect = if (response_type == "binary") fmt3(1.0) else fmt3(0.0),
              ciLow = "", ciHigh = "", pval = ""))
            for (i in seq_along(level_res)) {
              res <- level_res[[i]]
              gl  <- if ((i + 1) <= length(geno_labels)) geno_labels[i + 1] else sub("snp", "", res$term)
              emit(model_labels[[mdl]], grp1_lbl, gl, list(
                stat0 = if ((i+1) <= length(st$s0)) st$s0[i+1] else "-",
                stat1 = if ((i+1) <= length(st$s1)) st$s1[i+1] else " ",
                effect = fmt3(res$effect), ciLow = fmt3(res$ci_low),
                ciHigh = fmt3(res$ci_high), pval = fmt_pval(res$pval)))
            }
          }
        }
        list(rows = b$rows(), pvals = inter_pvals)
      })
      tbl$setNote(key = "interStratCov",
        note = paste0("Reference genotype is the first row of each stratum. ",
                      "Stratified by ", int_lbl, "."))
      if (length(built$pvals) > 0)
        tbl$setNote(key = "interStratCovPval",
          note = paste0("Interaction p-value — ",
                        paste(names(built$pvals), built$pvals, sep = ": ", collapse = "; ")))
      private$.render_flat(tbl, built$rows)
    },

    # One table per view, all models stacked (see .fill_strat_by_covariate).
    .fill_strat_by_genotype = function(tbl, snp_raw, ref, response, cov_df,
                                                interaction_var, response_type, opts,
                                                int_models, user_levels = NULL, response_raw, snp_lbl) {
      snp_char     <- as.character(snp_raw)
      all_genos    <- .all_genos_for_snp(user_levels, snp_char, ref)
      int_var_data <- cov_df[[interaction_var]]
      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      model_labels <- .MODEL_LABELS
      resp_lv      <- levels(as.factor(response_raw))
      is_numerical <- length(unique(int_var_data)) > 6 && sum(is.na(as.numeric(int_var_data))) == 0
      cov_levels   <- if (!is_numerical) {
        if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(as.character(int_var_data[!is.na(int_var_data)])))
      } else interaction_var
      tbl$getColumn("grp1")$setTitle(snp_lbl)
      tbl$getColumn("grp2")$setTitle(int_lbl)
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "β")
      if (response_type == "binary") {
        tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
        tbl$getColumn("stat0")$setVisible(TRUE);     tbl$getColumn("stat1")$setVisible(TRUE)
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }
      built <- private$.cached(tbl, function() {
        b <- private$.mk_flat_rows(); emit <- b$emit
        inter_pvals <- character(0)
        for (mdl in int_models) {
          if (mdl == "logadditive") next
          geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
          snp_enc_m   <- encode_model(snp_char, ref, mdl, user_levels)
          res_list    <- fit_interaction_model(snp_enc_m, response, cov_df, interaction_var, mdl,
                                               response_type, opts$ciWidth, conditional = TRUE, cond_var = "snp")
          if (is.null(res_list)) next
          p_i <- attr(res_list, "pval_interaction")
          if (!is.null(p_i) && !is.na(p_i)) inter_pvals[[model_labels[[mdl]]]] <- fmt_pval(p_i)
          n_cov_contrasts <- if (is_numerical) 1L else max(1L, length(cov_levels) - 1L)
          for (gl in geno_labels) {
            gl_idx   <- match(gl, geno_labels)
            grp1_lbl <- paste0(snp_lbl, ": ", gl)
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
            if (!is_numerical) {
              cl_ref <- cov_levels[1]
              stat0  <- if (response_type == "binary") fmt_cat(counts[cl_ref, 1], totals[1]) else { vals <- resp_g[as.character(int_g) == cl_ref]; fmt_cont(vals) }
              stat1  <- if (response_type == "binary") fmt_cat(counts[cl_ref, 2], totals[2]) else ""
              emit(model_labels[[mdl]], grp1_lbl, cl_ref, list(
                stat0 = stat0, stat1 = stat1,
                effect = if (response_type == "binary") fmt3(1.0) else fmt3(0.0), ciLow = "", ciHigh = "", pval = ""))
              for (i in seq_along(cov_levels[-1])) {
                cl    <- cov_levels[-1][i]
                res   <- if (i <= length(gl_res)) gl_res[[i]] else NULL
                stat0 <- if (response_type == "binary") fmt_cat(counts[cl, 1], totals[1]) else { vals <- resp_g[as.character(int_g) == cl]; fmt_cont(vals) }
                stat1 <- if (response_type == "binary") fmt_cat(counts[cl, 2], totals[2]) else ""
                emit(model_labels[[mdl]], grp1_lbl, cl, list(
                  stat0 = stat0, stat1 = stat1,
                  effect = if (!is.null(res)) fmt3(res$effect) else if (response_type == "binary") fmt3(1.0) else fmt3(0.0),
                  ciLow  = if (!is.null(res)) fmt3(res$ci_low)  else "",
                  ciHigh = if (!is.null(res)) fmt3(res$ci_high) else "",
                  pval   = if (!is.null(res)) fmt_pval(res$pval) else ""))
              }
            } else {
              stat0 <- if (response_type == "binary") fmt_cont(int_g[resp_raw_g == resp_lv[1]]) else fmt_cont(resp_g)
              stat1 <- if (response_type == "binary") fmt_cont(int_g[resp_raw_g == resp_lv[2]]) else ""
              res   <- if (length(gl_res) > 0) gl_res[[1]] else NULL
              emit(model_labels[[mdl]], grp1_lbl, "Overall", list(
                stat0 = stat0, stat1 = stat1,
                effect = if (!is.null(res)) fmt3(res$effect) else if (response_type == "binary") fmt3(1.0) else fmt3(0.0),
                ciLow  = if (!is.null(res)) fmt3(res$ci_low)  else "",
                ciHigh = if (!is.null(res)) fmt3(res$ci_high) else "",
                pval   = if (!is.null(res)) fmt_pval(res$pval) else ""))
            }
          }
        }
        list(rows = b$rows(), pvals = inter_pvals)
      })
      tbl$setNote(key = "interStratGeno",
        note = paste0("Reference is the first ", int_lbl, " level in each genotype stratum. ",
                      "Log-additive is omitted (no discrete genotype strata)."))
      if (length(built$pvals) > 0)
        tbl$setNote(key = "interStratGenoPval",
          note = paste0("Interaction p-value — ",
                        paste(names(built$pvals), built$pvals, sep = ": ", collapse = "; ")))
      private$.render_flat(tbl, built$rows)
    },

    # One table per view, all models stacked (see .fill_strat_by_covariate).
    .fill_cross_class = function(tbl, snp_raw, ref, response, cov_df,
                                          interaction_var, response_type, opts,
                                          int_models, user_levels = NULL, response_raw, snp_lbl) {
      int_var_data <- cov_df[[interaction_var]]
      int_lbl      <- attr(self$data[[interaction_var]], "label") %||% interaction_var
      snp_char     <- as.character(snp_raw)
      all_genos    <- .all_genos_for_snp(user_levels, snp_char, ref)
      model_labels <- .MODEL_LABELS
      cov_levels   <- if (is.factor(int_var_data)) levels(int_var_data) else sort(unique(int_var_data[!is.na(int_var_data)]))
      adj_vars     <- setdiff(names(cov_df), interaction_var)
      resp_lv      <- levels(as.factor(response_raw))
      tbl$getColumn("grp1")$setTitle(int_lbl)
      tbl$getColumn("grp2")$setTitle(snp_lbl)
      tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "β")
      if (response_type == "binary") {
        tbl$getColumn("stat0")$setTitle(resp_lv[1]); tbl$getColumn("stat1")$setTitle(resp_lv[2])
        tbl$getColumn("stat0")$setVisible(TRUE);     tbl$getColumn("stat1")$setVisible(response_type == "binary")
      } else {
        tbl$getColumn("stat0")$setTitle("Mean (SD)"); tbl$getColumn("stat0")$setVisible(TRUE)
        tbl$getColumn("stat1")$setVisible(FALSE)
      }
      built <- private$.cached(tbl, function() {
        if (length(table(int_var_data)) > 6) return(list(rows = list(), pvals = character(0)))
        b <- private$.mk_flat_rows(); emit <- b$emit
        inter_pvals <- character(0)
        for (mdl in int_models) {
          if (mdl == "logadditive") next
          snp_enc          <- encode_model(snp_char, ref, mdl, user_levels)
          df_fit           <- data.frame(resp = response, snp = snp_enc, interaction_var = int_var_data)
          if (length(adj_vars) > 0) df_fit <- cbind(df_fit, cov_df[, adj_vars, drop = FALSE])
          adj_part         <- if (length(adj_vars) > 0) paste("+", safe_rhs(adj_vars)) else ""
          formula_str      <- paste("resp ~ snp * interaction_var", adj_part)
          formula_main_str <- paste("resp ~ snp + interaction_var", adj_part)
          fit <- if (response_type == "binary")
            glm(as.formula(formula_str), data = df_fit, family = binomial())
          else lm(as.formula(formula_str), data = df_fit)
          fit_main_cc <- if (response_type == "binary")
            glm(as.formula(formula_main_str), data = df_fit, family = binomial())
          else lm(as.formula(formula_main_str), data = df_fit)
          lrtest       <- if (response_type == "binary") "Chisq" else "F"
          lrtest_label <- if (response_type == "binary") "Pr(>Chi)" else "Pr(>F)"
          lrt_cc       <- tryCatch(anova(fit_main_cc, fit, test = lrtest), error = function(e) NULL)
          p_inter_cc   <- if (!is.null(lrt_cc)) lrt_cc[2, lrtest_label] else NA_real_
          if (!is.na(p_inter_cc)) inter_pvals[[model_labels[[mdl]]]] <- fmt_pval(p_inter_cc)
          betas       <- coef(fit)
          v_cov       <- vcov(fit)
          ci_z        <- qnorm(1 - (1 - opts$ciWidth/100)/2)
          geno_labels <- private$.geno_labels_for_model(mdl, all_genos, ref)
          for (j in seq_along(cov_levels)) {
            cl       <- cov_levels[j]
            grp1_lbl <- paste0(int_lbl, ": ", as.character(cl))
            mask_k <- !is.na(int_var_data) & int_var_data == cl & !is.na(snp_raw)
            st <- private$.compute_stats(geno_labels, snp_char[mask_k], response[mask_k], response_type)
            for (i in seq_along(geno_labels)) {
              gl <- geno_labels[i]
              if (i == 1 && j == 1) {
                emit(model_labels[[mdl]], grp1_lbl, gl, list(
                  stat0 = st$s0[i], stat1 = st$s1[i],
                  effect = if (response_type == "binary") fmt3(1.0) else fmt3(0.0), ciLow = "", ciHigh = "", pval = ""))
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
              emit(model_labels[[mdl]], grp1_lbl, gl, list(
                stat0 = st$s0[i], stat1 = st$s1[i],
                effect = fmt3(if (response_type == "binary") .exp_or(combined_beta) else combined_beta),
                ciLow  = fmt3(if (response_type == "binary") .exp_or(lo_beta) else lo_beta),
                ciHigh = fmt3(if (response_type == "binary") .exp_or(hi_beta) else hi_beta),
                pval   = fmt_pval(p_val)))
            }
          }
        }
        list(rows = b$rows(), pvals = inter_pvals)
      })
      tbl$setNote(key = "interCrossClass",
        note = paste0("Reference cell is ", int_lbl, ": ", cov_levels[1],
                      " with the first genotype. Log-additive is omitted."))
      if (length(built$pvals) > 0)
        tbl$setNote(key = "interCrossClassPval",
          note = paste0("Interaction p-value — ",
                        paste(names(built$pvals), built$pvals, sep = ": ", collapse = "; ")))
      private$.render_flat(tbl, built$rows)
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
    # item          : one Group from ldResults Array (ldTable, ldMatrixTable, ldPlotImage)
    # complete_cases: when TRUE also exclude rows missing response/covariates
    # group_mask    : optional logical (n_rows) to restrict to one response group
    #
    # Each pair gets its OWN mask (both SNPs non-missing), so individuals missing
    # an unrelated third SNP are not excluded from the pair's computation.
    .run_ld = function(item, prep, opts, complete_cases,
                       need_table, need_matrix, need_plot, group_mask = NULL) {
      if (!(need_table || need_matrix || need_plot)) return()   # nothing to (re)fill
      nms   <- names(prep$snp_data)
      n     <- length(nms)
      if (n < 2) return()
      pairs <- combn(nms, 2, simplify = FALSE)
      ld_store <- list()
      for (pair in pairs) {
        sd1 <- prep$snp_data[[pair[1]]]; sd2 <- prep$snp_data[[pair[2]]]
        # Pairwise mask: only this pair needs to be non-missing
        pair_mask <- !is.na(sd1$clean) & !is.na(sd2$clean)
        if (complete_cases)         pair_mask <- pair_mask & prep$complete_mask
        if (!is.null(group_mask))   pair_mask <- pair_mask & group_mask
        g1 <- parse_genotype(sd1$clean[pair_mask], sd1$user_levels)
        g2 <- parse_genotype(sd2$clean[pair_mask], sd2$user_levels)
        if (is.null(g1) || is.null(g2)) next
        key    <- paste(pair, collapse = "___")
        ld_res <- tryCatch(genetics::LD(g1, g2), error = function(e) NULL)
        if (!is.null(ld_res)) ld_store[[key]] <- ld_res
      }
      if (need_table) {
        tbl <- item$ldTable
        tbl$deleteRows()                    # drop any .init-seeded rows before rebuild
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          tbl$addRow(rowKey = key, values = list(
            snp1   = pair[1], snp2 = pair[2],
            r2     = fmt3(ld_res$`r`^2),
            Dprime = fmt3(ld_res$`D'`),
            D      = fmt3(ld_res$`D`),
            pval   = fmt_pval(ld_res$`P-value`)))
        }
      }
      if (need_matrix) {
        mtbl   <- item$ldMatrixTable
        metric <- opts$ldMetric %||% "r2"
        # Add SNP columns dynamically, but only if .init did not already seed
        # them (re-adding an existing column would duplicate it).
        for (snp in nms) {
          safe_nm <- gsub("[^A-Za-z0-9_]", "_", snp)
          if (is.null(tryCatch(mtbl$getColumn(safe_nm), error = function(e) NULL)))
            mtbl$addColumn(name = safe_nm, title = snp, type = "text")
        }
        mtbl$deleteRows()                   # drop any .init-seeded rows before rebuild
        upper_mat <- matrix("", n, n)
        lower_mat <- matrix("", n, n)
        diag(upper_mat) <- nms
        diag(lower_mat) <- nms
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          i <- match(pair[1], nms)
          j <- match(pair[2], nms)
          p_str  <- fmt_pval(ld_res$`P-value`)
          up_val <- switch(metric,
            Dprime = fmt3(ld_res$`D'`),
            r2     = fmt3(ld_res$`r`^2),
            D      = fmt3(ld_res$`D`),
            fmt3(ld_res$`r`^2))
          upper_mat[i, j] <- up_val
          lower_mat[j, i] <- p_str
        }
        for (i in seq_len(n)) {
          row_vals <- list(snp = nms[i])
          for (j in seq_len(n)) {
            safe_nm <- gsub("[^A-Za-z0-9_]", "_", nms[j])
            row_vals[[safe_nm]] <- if (i==j) nms[i] else if (j>i) upper_mat[i,j] else lower_mat[i,j]
          }
          mtbl$addRow(rowKey = paste0("row_", i), values = row_vals)
        }
        metric_label <- switch(metric, Dprime="D'", r2="r²", D="D")
        mtbl$setNote(key = "layout",
                     note = paste0("Upper triangle: ", metric_label,
                                   ". Lower triangle: P-value. Diagonal: SNP name."))
      }
      if (need_plot) {
        item$ldPlotImage$setState(list(ld_store = ld_store, nms = nms, metric = opts$ldMetric))
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

    .run_haplo = function(geno_list, response, response_raw, response_type,
                           cov_df, opts, run_haploFreq, run_haploAssoc,
                           run_haploInteraction, run_subpop) {
      hg <- self$results$ldHaploGroup$haploGroup
      # Each haplotype table's clearWith blanks it when its content options change;
      # only recompute (haplo.em / haplo.glm are expensive) the ones that need it.
      need_freq  <- run_haploFreq && private$.need_fill(hg$haploFreqTable)
      need_assoc <- run_haploAssoc && !is.null(response) &&
                    private$.need_fill(hg$haploAssocTable)
      need_inter <- run_haploInteraction && !is.null(cov_df) && ncol(cov_df) > 0 &&
                    !is.null(response) && private$.need_fill(hg$haploInteractionTable)
      if (!need_freq && !need_assoc && !need_inter) return()

      snp_names   <- names(geno_list)
      allele_mat  <- do.call(cbind, lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]])))
      geno_setup  <- tryCatch(haplo.stats::setupGeno(allele_mat, locus.label = snp_names),
                              error = function(e) NULL)
      if (is.null(geno_setup)) return()
      u_alleles     <- attr(geno_setup, "unique.alleles")
      snp_miss_mask <- apply(is.na(allele_mat), 1, all)
      # geno_list is already on the joint mask, so complete_mask is all-TRUE here.
      keep          <- !snp_miss_mask
      n_miss        <- sum(snp_miss_mask)
      if (need_freq)
        private$.compute_haplo_freqs(geno_setup, response_raw, response_type, keep,
                                     n_miss, opts, run_subpop, snp_names, u_alleles)
      # Note the haplotype genetic model used for the GLM fits. The recessive
      # effect only makes sense when the region has many common haplotypes (few
      # homozygotes otherwise), so append that caveat when it is selected.
      model_note <- paste0("Haplotype genetic model: ", opts$haploEffect)
      if (identical(opts$haploEffect, "recessive"))
        model_note <- paste(model_note,
                            "— only meaningful when the region has many common haplotypes.")
      if (need_assoc) {
        hg$haploAssocTable$setNote("hmodel", model_note)
        private$.compute_haplo_assoc(geno_setup, response, response_type, cov_df, keep,
                                     n_miss, opts, snp_names, u_alleles)
      }
      if (need_inter) {
        hg$haploInteractionTable$setNote("hmodel", model_note)
        private$.compute_haplo_interaction(geno_setup, response, response_type, cov_df, keep,
                                           n_miss, opts, snp_names, u_alleles)
      }
    },

    .compute_haplo_freqs = function(geno_setup, response_raw, response_type, keep,
                                    n_miss, opts, run_subpop, snp_names, u_alleles) {
      tbl <- self$results$ldHaploGroup$haploGroup$haploFreqTable
      tbl$deleteRows()
      tbl$setTitle("<b>Haplotype Frequencies</b>")
      # categorical included: k-group stratification, columns added dynamically
      do_strat_haplo   <- isTRUE(run_subpop) && !is.null(response_raw) &&
                          response_type %in% c("binary", "categorical")
      grp_levels_haplo <- if (do_strat_haplo) levels(as.factor(response_raw[keep])) else character(0)
      if (do_strat_haplo && length(grp_levels_haplo) >= 2) {
        for (j in seq_along(grp_levels_haplo)) {
          tbl$addColumn(name = paste0("freq_g", j - 1L),
                        title = as.character(grp_levels_haplo[j]),
                        type = "number", format = "zto")
        }
      }
      # Cache the (expensive) EM fits in the table state so an unrelated click
      # reuses them instead of re-running haplo.em; the frequency rows are then
      # rebuilt cheaply. State survives restore gated by the table's clearWith.
      ems <- private$.cached(tbl, function() {
        # keep only hap.prob + haplotype: the full haplo.em object is huge
        # (per-subject posteriors) and must not be serialised into the state.
        slim <- function(e) if (is.null(e)) NULL else list(hap.prob = e$hap.prob, haplotype = e$haplotype)
        ea <- tryCatch(haplo.stats::haplo.em(subset_geno(geno_setup, keep), locus.label=snp_names),
                       error=function(e) NULL)
        eg <- list()
        if (do_strat_haplo && length(grp_levels_haplo) >= 2 && !is.null(ea)) {
          for (lvl in grp_levels_haplo) {
            keep_lvl <- keep & !is.na(response_raw) & as.character(response_raw)==lvl
            if (sum(keep_lvl) < 5) next
            eg[[lvl]] <- tryCatch(
              haplo.stats::haplo.em(subset_geno(geno_setup, keep_lvl), locus.label=snp_names),
              error=function(e) NULL)
          }
        }
        list(em_all = slim(ea), em_grp = lapply(eg, slim))
      })
      em_all <- ems$em_all; em_grp <- ems$em_grp
      if (!is.null(em_all)) {
        freqs    <- em_all$hap.prob
        rare_sum <- 0
        grp_freq <- list()
        if (do_strat_haplo && length(grp_levels_haplo) >= 2) {
          grp_freq <- lapply(em_grp, function(em_g) {
            if (is.null(em_g)) return(list())
            setNames(as.list(lapply(em_g$hap.prob, fmt3)),
                     sapply(seq_len(nrow(em_g$haplotype)), function(j)
                       decode_haplo_row(as.numeric(em_g$haplotype[j,]), u_alleles)))
          })
          for (j in seq_along(grp_levels_haplo)) {
            col_nm <- paste0("freq_g", j - 1L)
            tbl$getColumn(col_nm)$setVisible(TRUE)
            tbl$getColumn(col_nm)$setTitle(as.character(grp_levels_haplo[j]))
          }
        } else {
          tbl$getColumn('freq_g0')$setVisible(FALSE)
          tbl$getColumn('freq_g1')$setVisible(FALSE)
        }
        # Rare-haplotype pooling follows the selected criterion (haploRareCriterion):
        # a frequency cutoff, or a count converted to an expected frequency on the
        # EM fit's subject count so the display matches the haplo.glm control.
        freq_cut  <- haplo_rare_freq_cut(opts, sum(keep))
        rare_lbl  <- haplo_rare_label(opts)
        sorted_idx <- order(freqs, decreasing = TRUE)
        for (i in sorted_idx) {
          if (freqs[i] < freq_cut) { rare_sum <- rare_sum + freqs[i]; next }
          label    <- decode_haplo_row(as.numeric(em_all$haplotype[i,]), u_alleles)
          row_vals <- list(haplotype=label, freq=fmt3(freqs[i]))
          if (do_strat_haplo && length(grp_levels_haplo) >= 2) {
            for (j in seq_along(grp_levels_haplo)) {
              row_vals[[paste0("freq_g", j - 1L)]] <- grp_freq[[grp_levels_haplo[j]]][[label]] %||% NA_real_
            }
          }
          tbl$addRow(rowKey=paste0("f",i), values=row_vals)
        }
        if (rare_sum > 0) {
          row_vals <- list(haplotype=rare_lbl, freq=fmt3(rare_sum))
          if (do_strat_haplo && length(grp_levels_haplo) >= 2) {
            for (j in seq_along(grp_levels_haplo)) {
              em_g <- em_grp[[grp_levels_haplo[j]]]
              rare_val <- if (!is.null(em_g)) fmt3(sum(em_g$hap.prob[em_g$hap.prob < freq_cut])) else NA_real_
              row_vals[[paste0("freq_g", j - 1L)]] <- if (!is.na(rare_val) && rare_val > 0) rare_val else NA_real_
            }
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
      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        m_model    <- cbind(m_model, cov_df[keep, , drop = FALSE])
        formula_str <- paste("y ~ geno +", safe_rhs(names(cov_df)))
      } else {
        formula_str <- "y ~ geno"
      }
      null_formula_str <- if (!is.null(cov_df) && ncol(cov_df) > 0)
        paste("y ~", safe_rhs(names(cov_df))) else "y ~ 1"
      tbl <- self$results$ldHaploGroup$haploGroup$haploAssocTable
      # Cache the (small) products of the expensive haplo.glm \u2014 coefficient
      # table, CI matrix, the LRT p-value and the haplotype metadata \u2014 NOT the
      # fit object itself (it holds per-subject arrays and would bloat the saved
      # state). An unrelated click then reuses these instead of refitting.
      fitres <- private$.cached(tbl, function() {
        hf <- tryCatch(
          with_fixed_seed(haplo.stats::haplo.glm(as.formula(formula_str), family = family, data = m_model,
                                 na.action = na.geno.keep,
                                 control = haplo_glm_control(opts))),
          error = function(e) {
            self$results$validationMsg$setContent(paste0("<b>Haplotype GLM error:</b> ", e$message)); NULL
          })
        if (is.null(hf)) return(list(ok = FALSE))
        nf <- tryCatch(
          if (family == "binomial") glm(as.formula(null_formula_str), family = binomial(), data = m_model)
          else lm(as.formula(null_formula_str), data = m_model),
          error = function(e) NULL)
        p_lrt <- NA_real_
        if (!is.null(nf)) {
          # df is the number of haplotype (geno) terms in the fit. haplo.glm's EM
          # row-expansion makes hf$df.residual unusable and lm null models lack
          # $df.null (the old subtraction then yielded numeric(0) -> empty P for
          # quantitative responses), so count the coefficients directly, as the
          # interaction LRT does. For gaussian the deviance is on the variance
          # scale, so divide by the model dispersion to get a chi-square statistic
          # (binomial dispersion = 1); this matches haplo.glm's own $lrt.
          dev_diff <- tryCatch(deviance(nf) - hf$deviance, error = function(e) NA_real_)
          df_geno  <- sum(grepl("^geno", names(stats::coef(hf))))
          disp     <- if (family == "binomial") 1
                      else tryCatch(summary(hf)$dispersion, error = function(e) NA_real_)
          if (length(dev_diff) == 1L && !is.na(dev_diff) && df_geno > 0 &&
              length(disp) == 1L && !is.na(disp) && disp > 0)
            p_lrt <- pchisq(dev_diff / disp, df = df_geno, lower.tail = FALSE)
        }
        list(ok = TRUE, p_lrt = p_lrt,
             coef_sum = tryCatch(summary(hf)$coefficients, error = function(e) NULL),
             ci_mat   = tryCatch(confint(hf, level = opts$ciWidth / 100), error = function(e) NULL),
             haplo.base = hf$haplo.base, haplo.unique = hf$haplo.unique, haplo.freq = hf$haplo.freq,
             haplo.common = hf$haplo.common, haplo.rare = hf$haplo.rare,
             haplo.rare.term = hf$haplo.rare.term)
      })
      if (isTRUE(fitres$ok)) {
        tbl$deleteRows()
        tbl$setTitle("<b>Haplotype Association</b>")
        tbl$getColumn("effect")$setTitle(if (response_type == "binary") "OR" else "\u03B2")
        tbl$setNote(note = paste0("Likelihood ratio test for overall haplotype association: P = ",
                                  fmt_pval(fitres$p_lrt)), key = "lrt_assoc")
        label_from_unique_row <- function(row_vec) paste(as.character(row_vec), collapse = "-")
        coef_sum  <- fitres$coef_sum
        ci_mat    <- fitres$ci_mat
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
        base_idx   <- fitres$haplo.base
        base_label <- label_from_unique_row(fitres$haplo.unique[base_idx, ])
        base_freq  <- fitres$haplo.freq[base_idx]
        tbl$addRow(rowKey = "base", values = list(
          haplotype = paste0(base_label, " (Ref)"), freq = fmt3(base_freq),
          effect = if (response_type == "binary") 1.0 else 0.0, ciLow = '', ciHigh = '', pval = ''))
        common_idx   <- fitres$haplo.common
        common_freqs <- fitres$haplo.freq[common_idx]
        sorted_j     <- order(common_freqs, decreasing = TRUE)
        for (j in sorted_j) {
          h_idx   <- common_idx[j]
          h_label <- label_from_unique_row(fitres$haplo.unique[h_idx, ])
          h_freq  <- fitres$haplo.freq[h_idx]
          tbl$addRow(rowKey = paste0("h", j), values = make_row(h_label, h_freq, get_stats(j)))
        }
        has_rare <- isTRUE(fitres$haplo.rare.term) || (length(fitres$haplo.rare) > 0)
        if (has_rare) {
          rare_freq <- sum(fitres$haplo.freq[fitres$haplo.rare])
          tbl$addRow(rowKey = "rare",
            values = make_row(haplo_rare_label(opts), rare_freq,
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
      # Cache the (small) products of the two expensive haplo.glm fits — the
      # coefficient table, vcov, haplotype metadata and the LRT deviances — NOT
      # the fit objects (they hold per-subject arrays and would bloat the saved
      # state). An unrelated click then reuses these instead of refitting.
      hfits <- private$.cached(self$results$ldHaploGroup$haploGroup$haploInteractionTable, function() {
        fm <- tryCatch(
          with_fixed_seed(haplo.stats::haplo.glm(as.formula(formula_mult_str), family = family_int, data = m_int,
                                 na.action = na.geno.keep,
                                 control = haplo_glm_control(opts))),
          error = function(e) {
            self$results$validationMsg$setContent(paste0("<b>Haplotype interaction GLM error:</b> ", e$message)); NULL
          })
        if (is.null(fm)) return(list(ok = FALSE))
        fa <- tryCatch(
          with_fixed_seed(haplo.stats::haplo.glm(as.formula(formula_add_str), family = family_int, data = m_int,
                                 na.action = na.geno.keep,
                                 control = haplo_glm_control(opts))),
          error = function(e) NULL)
        list(ok = TRUE,
             coef_sum = tryCatch(summary(fm)$coefficients, error = function(e) NULL),
             vcov_mat = tryCatch(vcov(fm), error = function(e) NULL),
             haplo.base = fm$haplo.base, haplo.unique = fm$haplo.unique, haplo.freq = fm$haplo.freq,
             haplo.common = fm$haplo.common, haplo.rare = fm$haplo.rare, haplo.rare.term = fm$haplo.rare.term,
             mult_deviance = fm$deviance, mult_df.residual = fm$df.residual,
             add_deviance = if (!is.null(fa)) fa$deviance else NULL,
             add_df.residual = if (!is.null(fa)) fa$df.residual else NULL,
             has_add = !is.null(fa))
      })
      if (!isTRUE(hfits$ok)) return()
      coef_sum <- hfits$coef_sum
      if (is.null(coef_sum)) return()
      se_col   <- if ("SE" %in% colnames(coef_sum)) "SE" else "se"
      vcov_mat <- hfits$vcov_mat
      decode_haplo_label <- function(row_vec) paste(as.character(row_vec), collapse = "-")
      rare_label  <- haplo_rare_label(opts)
      z_crit      <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
      base_idx    <- hfits$haplo.base
      base_label  <- decode_haplo_label(hfits$haplo.unique[base_idx, ])
      int_var_factor <- as.factor(m_int[[int_var]])
      covar_levels   <- levels(int_var_factor)
      ref_covar_lvl  <- covar_levels[1]
      common_idx   <- hfits$haplo.common
      common_freqs <- hfits$haplo.freq[common_idx]
      has_rare     <- isTRUE(hfits$haplo.rare.term) || (length(hfits$haplo.rare) > 0)
      rare_freq    <- if (has_rare) sum(hfits$haplo.freq[hfits$haplo.rare]) else 0
      all_haplo_entries <- list()
      all_haplo_entries[["base"]] <- list(label = base_label, freq = hfits$haplo.freq[base_idx], coef_pos = NA_integer_)
      sorted_common_j <- order(common_freqs, decreasing = TRUE)
      for (j in sorted_common_j) {
        h_idx   <- common_idx[j]
        h_label <- decode_haplo_label(hfits$haplo.unique[h_idx, ])
        all_haplo_entries[[h_label]] <- list(label = h_label, freq = hfits$haplo.freq[h_idx], coef_pos = j)
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
        # .exp_or clamps separation-driven extremes to NA; suppress the cell
        # (shown blank) rather than printing "NA (NA-NA)".
        or <- .exp_or(beta)
        if (is.na(or)) return(NA_character_)
        sprintf("%.2f (%.2f\u2013%.2f)", or, .exp_or(beta - z_crit*se), .exp_or(beta + z_crit*se))
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
      if (isTRUE(hfits$has_add)) {
        dev_diff <- hfits$add_deviance - hfits$mult_deviance
        # df = number of haplotype x covariate interaction coefficients. haplo.glm
        # expands each subject into EM-weighted haplotype-pair rows and that
        # expansion differs between the additive and multiplicative fits, so their
        # df.residual difference is NOT the interaction parameter count; count the
        # geno:covariate terms directly instead.
        inter_terms <- rownames(coef_sum)[grepl(":", rownames(coef_sum)) &
                                          grepl("geno", rownames(coef_sum))]
        df_diff  <- length(inter_terms)
        p_inter  <- if (!is.na(dev_diff) && df_diff > 0)
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
      tbl_cross$deleteRows()
      tbl_cross$setTitle(paste0("<b>Haplotype \u00D7 ", int_var, " (cross-classification)</b>"))
      build_notes(tbl_cross)
      if (!is.na(p_inter)) tbl_cross$setNote(note = paste0("Interaction p-value (LRT): ", fmt_pval(p_inter)), key = "lrt_inter")
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
      tbl_cond_covar$deleteRows()
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
      if (!is.na(p_inter)) tbl_cond_covar$setNote(note = paste0("Interaction p-value (LRT): ", fmt_pval(p_inter)), key = "lrt_inter2")
      # Table 3: Covariate effect conditional on haplotype
      tbl_cond_haplo <- self$results$ldHaploGroup$haploGroup$haploCondHaploTable
      tbl_cond_haplo$deleteRows()
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
      if (!is.na(p_inter)) tbl_cond_haplo$setNote(note = paste0("Interaction p-value (LRT): ", fmt_pval(p_inter)), key = "lrt_inter2")
    }

  )  # end private
)
