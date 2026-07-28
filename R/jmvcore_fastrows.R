# ══════════════════════════════════════════════════════════════════════════════
# Fixed copies of two jmvcore Table methods
#
# jmvcore's Table$addRow recomputes EVERY row name on EVERY call:
#
#     private$.rowNames <- sapply(private$.rowKeys, toJSON, USE.NAMES = FALSE)
#
# so filling a table is quadratic in its row count — 2016 rows (the pairwise LD
# table for 64 SNPs) costs ~2 million toJSON calls, 41 s, against 1.5 s for all
# 2016 LD computations. Appending the one new name instead is 0.14 s.
#
# Table$deleteRows is the other half: it clears .rowKeys but leaves .rowNames
# behind. That is harmless only because the stock addRow overwrites the whole
# vector next time; an appending addRow would index into the stale one.
#
# Both fixes are one line each and have been reported upstream (see
# docs/jmvcore-addrow-pr/). Until they land, fast_rows() swaps these copies onto
# the individual Table objects this module owns — nothing global, no other
# module's tables, and no change to jmvcore's namespace.
#
# The copies are otherwise character-for-character the stock methods, so their
# observable behaviour — .rowNames content included, which protobuf restore
# matches on — is identical. That is asserted in tests/testthat/test-fastrows.R
# against the real jmvcore.
# ══════════════════════════════════════════════════════════════════════════════

# jmvcore's addRow with the row-name recompute replaced by an append.
# Unqualified names (isValue, reject, toJSON, private) resolve because
# fast_rows() sets this function's environment to the table's own R6 enclosure,
# whose parent is the jmvcore namespace — exactly where the stock method runs.
# They are fetched through get() so R CMD check does not see undefined globals.
.fixed_addRow <- function(rowKey, values = list()) {
  .jmv     <- parent.env(environment())
  isValue  <- get("isValue", envir = .jmv)
  toJSON   <- get("toJSON",  envir = .jmv)
  for (value in values)
    if (!isValue(value))
      get("reject", envir = .jmv)("Table$addRow(): value is not atomic", code = "error")
  private$.rowKeys[length(private$.rowKeys) + 1] <- list(rowKey)
  private$.rowCount <- private$.rowCount + 1
  private$.rowNames[[private$.rowCount]] <- toJSON(rowKey)   # was: recompute all
  valueNames <- names(values)
  for (column in private$.columns) {
    if (column$name %in% valueNames)
      column$addCell(values[[column$name]], .key = rowKey, .index = private$.rowCount)
    else
      column$addCell(.key = rowKey, .index = private$.rowCount)
  }
}

# jmvcore's deleteRows, additionally dropping the row names it leaves stale.
.fixed_deleteRows <- function() {
  private$.rowKeys  <- list()
  private$.rowNames <- character()                           # was: left behind
  for (column in private$.columns) column$clear()
  private$.rowCount <- 0
}

# The exact method bodies these copies were forked from. The guard below is an
# equality test against them, not a search for the bad line: if jmvcore changes
# either method AT ALL — the fix landing, or any unrelated edit — the match
# fails, the table keeps the stock methods, and the module is slow instead of
# subtly wrong. test-fastrows.R asserts the patch DOES apply, so a jmvcore
# upgrade that ends the match fails the suite loudly rather than silently
# reverting the speed-up.
.JMV_STOCK_ADDROW <- paste(
  '{ for (value in values) { if (!isValue(value)) reject("Table$addRow(): value',
  'is not atomic", code = "error") } private$.rowKeys[length(private$.rowKeys) +',
  '1] <- list(rowKey) private$.rowCount <- private$.rowCount + 1',
  'private$.rowNames <- sapply(private$.rowKeys, toJSON, USE.NAMES = FALSE)',
  'valueNames <- names(values) for (column in private$.columns) { if (column$name',
  '%in% valueNames) column$addCell(values[[column$name]], .key = rowKey, .index =',
  'private$.rowCount) else column$addCell(.key = rowKey, .index =',
  'private$.rowCount) } }')

.JMV_STOCK_DELETEROWS <- paste(
  '{ private$.rowKeys <- list() for (column in private$.columns) column$clear()',
  'private$.rowCount <- 0 }')

.squish_body <- function(f) gsub("[[:space:]]+", " ", paste(deparse(body(f)), collapse = " "))

# Swap the fixed methods onto one Table, and report whether it took.
fast_rows <- function(tbl) {
  if (is.null(tbl) || !inherits(tbl, "Table")) return(invisible(FALSE))
  if (identical(body(tbl$addRow), body(.fixed_addRow))) return(invisible(TRUE))  # already done

  ee   <- tryCatch(tbl$.__enclos_env__, error = function(e) NULL)
  priv <- if (is.null(ee)) NULL else ee$private
  if (is.null(priv)) return(invisible(FALSE))
  if (!all(c(".rowKeys", ".rowCount", ".rowNames", ".columns") %in% ls(priv, all.names = TRUE)))
    return(invisible(FALSE))
  if (!identical(.squish_body(tbl$addRow),     .JMV_STOCK_ADDROW) ||
      !identical(.squish_body(tbl$deleteRows), .JMV_STOCK_DELETEROWS))
    return(invisible(FALSE))

  ok <- tryCatch({
    for (nm in c("addRow", "deleteRows")) {
      fixed <- if (nm == "addRow") .fixed_addRow else .fixed_deleteRows
      environment(fixed) <- ee
      unlockBinding(nm, tbl)
      assign(nm, fixed, envir = tbl)
      lockBinding(nm, tbl)
    }
    TRUE
  }, error = function(e) FALSE)
  invisible(ok)
}
