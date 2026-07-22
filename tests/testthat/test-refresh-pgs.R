# Refresh / restore regression tests for snpPGS.
#
# Replays jamovi's option-click cycle (fresh object -> .init -> .load restore ->
# .run) the way the engine drives it: with the FULL option set. A sparse option
# set makes jmvcore's compProtoBuf report every unsent option as changed, which
# clears every table on restore and hides the real behaviour (see the note in
# test-refresh.R). Requires RProtoBuf.

skip_if_not_installed("RProtoBuf")

jmvcore:::initProtoBuf()

# ── Inverse of jmvcore:::parseOptionPB — build an option protobuf ────────────
.pgs_optPB <- function(v) {
  pb <- RProtoBuf::new(jamovi.coms.AnalysisOption)
  if (is.logical(v) && length(v) == 1L)        pb$o <- if (isTRUE(v)) 1L else 0L
  else if (is.character(v) && length(v) == 1L) pb$s <- v
  else if (is.numeric(v) && length(v) == 1L)   pb$d <- as.numeric(v)
  else if (is.null(v))                         pb$o <- 2L
  else {
    inner <- RProtoBuf::new(jamovi.coms.AnalysisOptions)
    inner$options  <- lapply(v, .pgs_optPB)
    inner$hasNames <- FALSE
    pb$c <- inner
  }
  pb
}

.pgs_optionsPB <- function(vals) {
  pb <- RProtoBuf::new(jamovi.coms.AnalysisOptions)
  pb$names    <- names(vals)
  pb$hasNames <- TRUE
  pb$options  <- lapply(vals, .pgs_optPB)
  pb
}

# Full defaults, minus saveScores (an Output-type option that cannot be set as a
# plain value through fromProtoBuf); the engine handles it separately.
.pgs_defaults <- local({
  o <- snpPGSOptions$new()
  d <- as.list(o$values())[o$names]
  d[setdiff(names(d), "saveScores")]
})

.pgs_mk <- function(over = list()) {
  vals <- .pgs_defaults
  vals[names(over)] <- over
  o <- snpPGSOptions$new()
  o$fromProtoBuf(.pgs_optionsPB(vals))
  o
}

# ── Weights fixture (four dataset SNPs) ──────────────────────────────────────
.pgs_snps4    <- c("rs12080929", "rs10911251", "rs10936599", "rs6691170")
.pgs_wfile <- local({
  f <- tempfile(fileext = ".tsv")
  writeLines("# test PGS weights file", f)
  df <- data.frame(rsID = .pgs_snps4,
                   effect_allele = c("C", "C", "T", "T"),
                   other_allele  = c("T", "A", "C", "G"),
                   effect_weight = c(0.5, -0.3, 0.8, 0.2),
                   chr_name = seq_along(.pgs_snps4),
                   chr_position = seq_along(.pgs_snps4) * 100L)
  suppressWarnings(write.table(df, f, sep = "\t", row.names = FALSE,
                               quote = FALSE, append = TRUE))
  f
})

.pgs_base <- list(snpCols = as.list(.pgs_snps4), weightsPath = .pgs_wfile,
                  weightingMode = "both", responseCol = "phenotype",
                  covCols = list("sex"), showCoverage = TRUE, showSnpGrid = TRUE,
                  showAssoc = TRUE, showInteraction = TRUE, showPercentiles = TRUE)

.summary_of <- function(a) a$results$summaryTable

test_that("summaryTable keeps its stratified rows after an unrelated click", {
  # Regression: .init used to seed only the {mode}_Overall rows. On a binary
  # response .run adds {mode}_Case / {mode}_Control too, but after restore only
  # the Overall rows came back filled, the isNotFilled() gate read "done", and
  # the stratified rows were dropped — the table silently fell from 6 rows to 2
  # on any unrelated click. .init now seeds every group row.
  statefile <- tempfile(fileext = ".pb")

  a0 <- snpPGSClass$new(options = .pgs_mk(.pgs_base), data = .test_data,
                        analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base <- .summary_of(a0)$asDF
  expect_equal(nrow(base), 6L)                       # 2 modes x (Case, Control, Overall)

  # plotWidth is not in summaryTable's clearWith — an unrelated click.
  a1 <- snpPGSClass$new(options = .pgs_mk(c(.pgs_base, list(plotWidth = 700))),
                        data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  expect_equal(.summary_of(a1)$rowCount, 6L)         # all rows restored before .run
  a1$run()
  expect_equal(.summary_of(a1)$asDF, base)           # full stratified content, unchanged
})

test_that("snpGridTable survives an unrelated click without a rebuild", {
  statefile <- tempfile(fileext = ".pb")
  a0 <- snpPGSClass$new(options = .pgs_mk(.pgs_base), data = .test_data,
                        analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base <- a0$results$snpGridTable$asDF
  expect_equal(nrow(base), length(.pgs_snps4))

  a1 <- snpPGSClass$new(options = .pgs_mk(c(.pgs_base, list(plotWidth = 700))),
                        data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  expect_equal(a1$results$snpGridTable$rowCount, nrow(base))
  a1$run()
  expect_equal(a1$results$snpGridTable$asDF, base)
})

test_that("snpGridTable is NOT refilled on an unrelated click (incl. active sort)", {
  # The grid used to rebuild on every option click because its refill gate keyed
  # on instance state (lost on restore). It is now refresh-safe (positional keys
  # + setRow), so a plot toggle must leave the restored rows untouched — proven
  # by a sentinel poked into a restored cell that a refill would overwrite.
  check_no_refill <- function(extra0) {
    statefile <- tempfile(fileext = ".pb")
    a0 <- snpPGSClass$new(options = .pgs_mk(c(.pgs_base, extra0)),
                          data = .test_data, analysisId = 1, revision = 1)
    a0$.setStatePathSource(function() statefile)
    a0$init(); a0$run(); a0$.save()

    a1 <- snpPGSClass$new(options = .pgs_mk(c(.pgs_base, extra0,
                            list(showRocPlot = TRUE))),   # unrelated toggle
                          data = .test_data, analysisId = 1, revision = 2)
    a1$.setStatePathSource(function() statefile)
    a1$init(); a1$.load()
    a1$results$snpGridTable$setRow(rowKey = "1", values = list(pct_missing = -999))
    a1$run()
    expect_equal(a1$results$snpGridTable$asDF$pct_missing[1], -999)  # sentinel survived
  }
  check_no_refill(list())                                # default sort
  check_no_refill(list(snpGridSortBy = "effect_weight")) # active sort
})

# ── The tables that used to blank/rebuild on every option click ──────────────
# Each was fixed by pre-creating its predicted rows in .init() (positional keys)
# and writing them by position in .run() (see .writeRows). plotWidth is in none
# of their clearWith lists, so after an unrelated click restore must refill the
# rows and .run must leave them untouched — same rowCount, identical asDF.
.pgs_survives_click <- function(get_tbl, extra = list()) {
  statefile <- tempfile(fileext = ".pb")
  opts0 <- c(.pgs_base, extra)
  a0 <- snpPGSClass$new(options = .pgs_mk(opts0), data = .test_data,
                        analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base <- get_tbl(a0)$asDF
  expect_gt(nrow(base), 0L)

  a1 <- snpPGSClass$new(options = .pgs_mk(c(opts0, list(plotWidth = 701))),
                        data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  # rows restored before .run (no rebuild needed)
  expect_equal(get_tbl(a1)$rowCount, nrow(base))
  a1$run()
  expect_equal(get_tbl(a1)$asDF, base)
}

test_that("coverageTable survives an unrelated click", {
  .pgs_survives_click(function(a) a$results$coverageTable)
})

test_that("assocTable survives an unrelated click", {
  .pgs_survives_click(function(a) a$results$assocTable)
})

test_that("interactionTable survives an unrelated click", {
  .pgs_survives_click(function(a) a$results$interactionTable)
})

test_that("percentileTable and percentileThreshTable survive an unrelated click", {
  .pgs_survives_click(function(a) a$results$percentileTable)
  .pgs_survives_click(function(a) a$results$percentileThreshTable)
})

test_that("continuous-response assoc/percentile survive an unrelated click", {
  # exercises the continuous branch (different predicted row counts)
  ext <- list(responseCol = "age")   # numeric response in the test data
  .pgs_survives_click(function(a) a$results$assocTable, ext)
  .pgs_survives_click(function(a) a$results$percentileTable, ext)
})

# ── Table notes survive restore (no empty footnotes) ─────────────────────────
# Regression: run-phase notes were set init=TRUE, so protobuf restore blanked
# their text and the gated .run() fill never re-set them — the assoc table
# showed two empty footnotes. Notes are now init=FALSE and their keys are
# cleared in .init so no empty placeholder survives.
.note_text <- function(tbl, key) {
  n <- tbl$notes[[key]]
  if (is.null(n)) NA_character_ else as.character(n$note)
}

test_that("assocTable covNote/respNote keep their text across an unrelated click", {
  statefile <- tempfile(fileext = ".pb")
  a0 <- snpPGSClass$new(options = .pgs_mk(.pgs_base), data = .test_data,
                        analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  cov0 <- .note_text(a0$results$assocTable, "covNote")
  resp0 <- .note_text(a0$results$assocTable, "respNote")
  expect_true(nzchar(cov0) && nzchar(resp0))

  a1 <- snpPGSClass$new(options = .pgs_mk(c(.pgs_base, list(plotWidth = 702))),
                        data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load(); a1$run()
  expect_equal(.note_text(a1$results$assocTable, "covNote"), cov0)
  expect_equal(.note_text(a1$results$assocTable, "respNote"), resp0)
})

# ── QC-filter changes must refresh the tables (not show stale content) ───────
# Regression: the QC-filter options were in no table's clearWith / refill key, so
# changing a threshold left the grid, coverage and score tables showing the old
# exclusions. They are now in clearWith (coverage/summary/...) and the grid's
# refill key, so a QC change updates them to match a fresh recompute.
test_that("changing a QC threshold refreshes grid/coverage/summary", {
  # Two of the four fixture SNPs are ~0.39% missing, the others < 0.1%, so
  # pct=0.1 excludes two SNPs and pct=0.5 excludes none — a discriminating pair.
  qcbase <- c(.pgs_base[setdiff(names(.pgs_base), "covCols")],
              list(qcFilterMissing = TRUE, qcMaxMissingPct = 0.1))
  statefile <- tempfile(fileext = ".pb")
  a0 <- snpPGSClass$new(options = .pgs_mk(qcbase), data = .test_data,
                        analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  strict <- list(grid = a0$results$snpGridTable$asDF,
                 cov  = a0$results$coverageTable$asDF,
                 sm   = a0$results$summaryTable$asDF)

  # relax the missingness threshold (was not in any clearWith / refill key before)
  a1 <- snpPGSClass$new(options = .pgs_mk(c(qcbase, list(qcMaxMissingPct = 0.5))),
                        data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load(); a1$run()
  relaxed <- list(grid = a1$results$snpGridTable$asDF,
                  cov  = a1$results$coverageTable$asDF,
                  sm   = a1$results$summaryTable$asDF)

  # fresh full recompute at the relaxed threshold
  b <- snpPGSClass$new(options = .pgs_mk(c(qcbase, list(qcMaxMissingPct = 0.5))),
                       data = .test_data, analysisId = 1, revision = 1)
  b$init(); b$run()
  fresh <- list(grid = b$results$snpGridTable$asDF,
                cov  = b$results$coverageTable$asDF,
                sm   = b$results$summaryTable$asDF)

  # updated (not stale) and identical to the fresh recompute
  expect_false(isTRUE(all.equal(relaxed$grid, strict$grid)))
  expect_equal(relaxed$grid, fresh$grid)
  expect_equal(relaxed$cov,  fresh$cov)
  expect_equal(relaxed$sm,   fresh$sm)
})

test_that("no empty note placeholder survives when the note does not apply", {
  # No covariates -> assoc covNote must be absent (not an empty footnote); the
  # single-covariate interaction likewise has no intNote.
  base_nocov <- .pgs_base[setdiff(names(.pgs_base), "covCols")]
  base_nocov$showInteraction <- FALSE
  statefile <- tempfile(fileext = ".pb")
  a0 <- snpPGSClass$new(options = .pgs_mk(base_nocov), data = .test_data,
                        analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()

  a1 <- snpPGSClass$new(options = .pgs_mk(c(base_nocov, list(plotWidth = 703))),
                        data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load(); a1$run()
  expect_true(is.na(.note_text(a1$results$assocTable, "covNote")))
  expect_true(nzchar(.note_text(a1$results$assocTable, "respNote")))
})
