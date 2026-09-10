# What jamovi does between clicks, which run_import() does not: a fresh analysis
# object per request, .init() on empty tables, then the previous results
# restored into it, then .run().
#
# The 2026-09 jamovi review read .init()'s unconditional addRows as a
# duplicate-row bug on any option change outside a table's clearWith. It is not
# one -- the engine builds the object fresh (enginer.cpp: create(), then init(),
# then .load(), then run()) and Table$fromProtoBuf only copies cells into rows
# that already exist -- but nothing in the suite said so, because every other
# test constructs one analysis and runs it once. This is that test.
#
# asProtoBuf is used directly rather than through .save()/.load(): the file
# round trip needs an analysisid and revision the engine assigns, and the bytes
# in between are not what is under test.

ns <- asNamespace("SNPstats")

# One request. `restore` is the results element of the previous cycle, or NULL
# for the first, and `changed` the option names jamovi says changed.
cycle <- function(restore = NULL, changed = character(), ...) {
  opts <- do.call(ns$snpImportOptions$new, list(...))
  a <- ns$snpImportClass$new(options = opts, data = data.frame())
  a$init(noThrow = TRUE)
  if (!is.null(restore))
    a$results$fromProtoBuf(restore, changed, character())
  a$run(noThrow = TRUE)
  a
}

test_that("rows do not accumulate across option changes outside clearWith", {
  skip_without_fixture("small")
  jmvcore:::initProtoBuf()

  ids <- read_bim_file("small")$id[c(1, 5, 10, 20, 40)]
  p   <- make_payloads("small", ids)
  base <- list(genoContent = p$geno %||% p$bed, variantContent = p$bim,
               sampleContent = p$fam, snpListText = paste(ids, collapse = "\n"),
               showSummary = TRUE, showProvenance = TRUE)

  a1 <- do.call(cycle, base)
  n_prov <- a1$results$provenance$rowCount
  expect_equal(a1$results$summary$rowCount, length(ids))
  expect_equal(a1$results$openNewState$rowCount, 1L)

  # dosage is in neither table's clearWith, so jamovi hands the old cells back.
  a2 <- do.call(cycle, c(list(restore = a1$results$asProtoBuf(),
                              changed = "dosage", dosage = TRUE), base))
  expect_equal(a2$results$provenance$rowCount, n_prov)
  expect_equal(a2$results$summary$rowCount, length(ids))
  expect_equal(a2$results$openNewState$rowCount, 1L)

  # And again, so a second restore on top of a restored state is covered too.
  a3 <- do.call(cycle, c(list(restore = a2$results$asProtoBuf(),
                              changed = "dosage", dosage = FALSE), base))
  expect_equal(a3$results$provenance$rowCount, n_prov)
  expect_equal(a3$results$summary$rowCount, length(ids))

  # The keys, not just the count: a duplicated key would still be a duplicate
  # if something else had gone missing.
  expect_equal(unlist(a3$results$summary$rowKeys), which(rep(TRUE, length(ids))))
  expect_false(any(duplicated(unlist(a3$results$provenance$rowKeys))))
})

test_that("a table whose row set changes is rebuilt, not appended to", {
  skip_without_fixture("small")
  jmvcore:::initProtoBuf()

  ids <- read_bim_file("small")$id[c(1, 5, 10, 20, 40)]
  p   <- make_payloads("small", ids)
  base <- list(genoContent = p$geno %||% p$bed, variantContent = p$bim,
               sampleContent = p$fam, snpListText = paste(ids, collapse = "\n"),
               showSummary = TRUE, showProvenance = TRUE)

  a1 <- do.call(cycle, base)
  n1 <- a1$results$provenance$rowCount

  # applyFilters adds "SNPs kept after filters" to .provKeys() and is in
  # neither clearWith, so this is the case where the restored state and the
  # wanted row set genuinely disagree.
  a2 <- do.call(cycle, c(list(restore = a1$results$asProtoBuf(),
                              changed = "applyFilters", applyFilters = TRUE),
                         base))
  expect_equal(a2$results$provenance$rowCount, n1 + 1L)
  expect_true("SNPs kept after filters" %in%
                unlist(a2$results$provenance$rowKeys))
  expect_false(any(duplicated(unlist(a2$results$provenance$rowKeys))))

  # Every cell is filled, i.e. .fillProvenance's key lookup still matched.
  vals <- as.data.frame(a2$results$provenance)$value
  expect_true(all(nzchar(vals)))
})

test_that(".syncRows rebuilds a table that already holds the wrong rows", {
  # The guard itself, exercised directly: the engine never hands .init() a
  # populated table, so nothing above can reach the deleteRows branch.
  tbl <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "probe")
  tbl$addColumn(name = "item", type = "text")
  tbl$addRow(rowKey = "a", values = list(item = "a"))
  tbl$addRow(rowKey = "b", values = list(item = "b"))

  a <- ns$snpImportClass$new(options = ns$snpImportOptions$new(),
                             data = data.frame())
  sync <- a$.__enclos_env__$private$.syncRows

  expect_true(sync(tbl, list("a", "c"), list(item = c("a", "c"))))
  expect_equal(unlist(tbl$rowKeys), c("a", "c"))
  expect_equal(as.data.frame(tbl)$item, c("a", "c"))

  # Same keys: left alone, so restored cells are not thrown away.
  expect_false(sync(tbl, list("a", "c"), list(item = c("x", "y"))))
  expect_equal(as.data.frame(tbl)$item, c("a", "c"))
})
