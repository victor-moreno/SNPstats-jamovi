
# openNew is the only way this analysis produces data: option$perform()
# writes .buildColumns()'s columns to a temp .omv on disk and hands it to
# jamovi's file-open pipeline. See the comments on openNewState in
# snpimport.r.yaml / .performOpen in snpimport.b.R for the re-fire guard this
# exercises. What .buildColumns() actually builds -- correct column names and
# values for every format, filter and covariate combination -- is covered by
# the rest of the suite through helper-fixtures.R's out_values(), which reads
# it directly (see the comment there for why that does not need perform()).
#
# jmvReadWrite -- what jmvcore's Action$perform() calls to write the temp
# .omv -- ships with jamovi's own engine image, not with this project's local
# R library (CLAUDE.md), so it is not installed here. Tests that need it to
# actually run skip when it is absent; what can be checked without it (the
# columns built, the re-fire guard) is checked unconditionally.

test_that("openNew is a plain constructor argument, not an Output needing private state", {
  # Unlike the Output option this replaced (not a parameter of the generated
  # function -- see helper-fixtures.R's run_import), openNew is an ordinary
  # Action option and takes its value straight from the constructor.
  ns   <- asNamespace("SNPstats")
  opts <- ns$snpImportOptions$new(openNew = TRUE)
  expect_true(opts$openNew)

  opts2 <- ns$snpImportOptions$new()
  expect_false(opts2$openNew)
})

test_that("openNewState's row exists from .init, empty until openNew fires", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  p   <- make_payloads("small", ids)

  ns   <- asNamespace("SNPstats")
  opts <- ns$snpImportOptions$new(
    genoContent = p$bed, variantContent = p$bim, sampleContent = p$fam,
    snpListText = paste(ids, collapse = "\n"), showSummary = FALSE)
  a <- ns$snpImportClass$new(options = opts, data = data.frame())
  a$run()

  expect_equal(a$results$openNewState$rowCount, 1)
  expect_identical(a$results$openNewState$getCell(col = "fired", rowNo = 1)$value, "")
})

# perform()'s write_omv call needs analysis$.getSessionTemp() to resolve to a
# real, writable directory. jmvcore's .getSessionTemp() (analysis.R) reads it
# from private$.resourcesPathSource, which the engine wires in via
# .setResourcesPathSource() and a headless analysis object never receives --
# left unset it is NA, and .getSessionTemp() fails with "attempt to apply
# non-function". No other test in this suite needs this stub, because none of
# them call an Action option's perform().
with_session_temp <- function(a) {
  tmp <- tempfile("snpimport-sessiontemp-")
  dir.create(file.path(tmp, "temp"), recursive = TRUE)
  a$.setResourcesPathSource(function(name, ext)
    list(rootPath = file.path(tmp, "x")))
  invisible(tmp)
}

test_that("openNew performs at most once per analysis instance", {
  skip_without_fixture("small")
  skip_if_not_installed("jmvReadWrite")

  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  p   <- make_payloads("small", ids)

  ns   <- asNamespace("SNPstats")
  opts <- ns$snpImportOptions$new(
    genoContent = p$bed, variantContent = p$bim, sampleContent = p$fam,
    snpListText = paste(ids, collapse = "\n"), showSummary = FALSE,
    openNew = TRUE)
  a <- ns$snpImportClass$new(options = opts, data = data.frame())
  with_session_temp(a)
  a$run()

  expect_identical(a$results$openNewState$getCell(col = "fired", rowNo = 1)$value, "yes")
  # 'openNew' the results item exists only as whatever perform() dynamically
  # added to results$.items -- it is not a declared item (that crashes, see
  # the r.yaml comment on openNewState), so it is reachable through the
  # generic Group$get()/itemNames methods but not as a$results$openNew, which
  # errors "does not exist in this results element".
  expect_true("openNew" %in% a$results$itemNames)
  openNewResults <- a$results$get("openNew")
  expect_length(openNewResults$itemNames, 1)
  # perform() (jmvcore's, options.R) only sets status on the error path; a
  # successful write_omv leaves it unset and reports the temp .omv's path
  # instead, with the data.frame itself stripped back out of the result.
  result <- openNewResults$get(1)$result
  expect_null(result$status)
  expect_true(nzchar(result$path))
  expect_null(result$data)

  # What a rebuilt instance sees on the next option click: openNew is still
  # TRUE (jamovi never resets it) and openNewState's row is restored with
  # 'fired' already set -- .r.yaml's justification for using that table at
  # all. .performOpen must not call perform() again: a second call is
  # observable because perform() requires self$value (it does, openNew is
  # still TRUE) but would add a *second* item to results$openNew if it ran.
  opts2 <- ns$snpImportOptions$new(
    genoContent = p$bed, variantContent = p$bim, sampleContent = p$fam,
    snpListText = paste(ids, collapse = "\n"), showSummary = FALSE,
    openNew = TRUE)
  a2 <- ns$snpImportClass$new(options = opts2, data = data.frame())
  with_session_temp(a2)
  # .init() has to run first to create openNewState's row (as it does inside
  # $run(), before .run()) -- poking setRow any earlier is poking a table with
  # no rows yet, which is not what a restore onto an already-initialised
  # results tree looks like.
  a2$init()
  a2$results$openNewState$setRow(rowNo = 1, values = list(fired = "yes"))
  a2$run()

  expect_false("openNew" %in% a2$results$itemNames)
})

test_that("openNew with nothing to write reports an error result, not a silent no-op", {
  skip_if_not_installed("jmvReadWrite")

  ns <- asNamespace("SNPstats")
  # loaded data, but filtered down to zero SNPs and sample columns off: a
  # different route to an empty .buildColumns() than "nothing loaded", which
  # .haveData() catches earlier and .performOpen never even reaches.
  tfam <- c("f1 s1 0 0 1 1", "f1 s2 0 0 2 1")
  geno <- as_payload("1 snpA 0 100 A A A G")
  opts <- ns$snpImportOptions$new(
    genoContent = geno, sampleContent = as_payload(tfam),
    sourceFormat = "tped", snpListText = "snpA", showSummary = FALSE,
    emitSamples = FALSE, applyFilters = TRUE, minMaf = 0.5, openNew = TRUE)
  a <- ns$snpImportClass$new(options = opts, data = data.frame())
  with_session_temp(a)
  a$run()

  item <- a$results$get("openNew")$get(1)
  expect_identical(item$result$status, "error")
})
