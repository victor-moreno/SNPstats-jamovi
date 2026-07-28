# fast_rows() swaps fixed copies of jmvcore's Table$addRow / $deleteRows onto
# the module's own tables (see R/jmvcore_fastrows.R). The copies must be
# behaviourally indistinguishable from the stock methods — including the
# .rowNames vector, which is what protobuf restore matches saved cells on. Every
# test here compares a patched table against a stock one built the same way.

fast_rows    <- getFromNamespace("fast_rows", "SNPstats")
.rownames_of <- function(tbl) tbl$.__enclos_env__$private$.rowNames

.new_table <- function(patch) {
  tbl <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "t")
  tbl$addColumn(name = "a", type = "text")
  tbl$addColumn(name = "b", type = "number")
  if (patch) expect_true(fast_rows(tbl))
  tbl
}

.fill <- function(tbl, keys) {
  for (k in keys) tbl$addRow(rowKey = k, values = list(a = paste0("v", k), b = nchar(k)))
  tbl
}

test_that("the patched addRow is indistinguishable from the stock one", {
  keys <- c("k1", "rs1___rs2", "a b", "quote\"inside", "10")
  stock   <- .fill(.new_table(FALSE), keys)
  patched <- .fill(.new_table(TRUE),  keys)

  expect_equal(patched$rowCount, stock$rowCount)
  expect_equal(as_df(patched), as_df(stock))
  expect_equal(patched$rowKeys, stock$rowKeys)
  expect_equal(.rownames_of(patched), .rownames_of(stock))
})

test_that("the patched deleteRows leaves a table that refills identically", {
  # The stock deleteRows leaves .rowNames behind; harmless only because the
  # stock addRow overwrites the whole vector. An appending addRow would index
  # into the stale one, so the fixed deleteRows must clear it.
  stock   <- .fill(.new_table(FALSE), c("x", "y", "z"))
  patched <- .fill(.new_table(TRUE),  c("x", "y", "z"))
  stock$deleteRows();   .fill(stock,   c("p", "q"))
  patched$deleteRows(); .fill(patched, c("p", "q"))

  expect_equal(patched$rowCount, 2L)
  expect_equal(as_df(patched), as_df(stock))
  expect_equal(.rownames_of(patched), .rownames_of(stock))
})

test_that("fast_rows is idempotent and only touches the table it is given", {
  a <- .new_table(TRUE)
  expect_true(fast_rows(a))                   # second call is a no-op, still TRUE
  b <- .new_table(FALSE)
  expect_true(grepl("sapply", paste(deparse(body(b$addRow)), collapse = " ")))
  expect_false(fast_rows(NULL))
})

test_that("row-name maintenance survives the real restore path", {
  # The point of matching .rowNames exactly: fromProtoBuf looks saved cells up by
  # row name. Save a patched table and restore into a fresh one; the cells must
  # come back.
  skip_if_not_installed("RProtoBuf")
  jmvcore:::initProtoBuf()

  src <- .fill(.new_table(TRUE), c("k1", "k2", "k3"))
  pb  <- src$asProtoBuf()
  dst <- .fill(.new_table(TRUE), c("k1", "k2", "k3"))
  for (i in seq_len(dst$rowCount)) dst$setRow(rowNo = i, values = list(a = "", b = 0))
  dst$fromProtoBuf(pb, character(), character())

  expect_equal(as_df(dst), as_df(src))
})

test_that("the LD table really is patched in a live analysis", {
  opts <- snpStatsOptions$new(snps = as.list(.snps4), covDesc = FALSE,
                              snpSummary = FALSE, ldAnalysis = TRUE)
  a <- snpStatsClass$new(options = opts, data = .test_data,
                         analysisId = 1, revision = 1)
  a$init()

  tbl <- a$results$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldTable
  expect_equal(tbl$rowCount, choose(length(.snps4), 2))
  expect_false(grepl("sapply", paste(deparse(body(tbl$addRow)), collapse = " ")))
  # .init's pre-created keys are what .run_ld setRows into
  expect_equal(.rownames_of(tbl),
               vapply(combn(.snps4, 2, simplify = FALSE),
                      function(p) paste0("\"", paste(p, collapse = "___"), "\""), character(1)))
})
