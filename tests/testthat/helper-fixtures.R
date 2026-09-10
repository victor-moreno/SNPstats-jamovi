
# Test helpers.
#
# The package itself never opens a file — `grep -rn "file(" R/`
# finds nothing, and that is a property worth keeping. So the seeking that the
# browser does in production is done here instead, reusing the package's own
# decoder for everything below the seek.

fixture_dir <- function(tier = "small") {
  for (p in c(file.path("..", "..", "data-raw", "fixtures", tier),
              file.path("data-raw", "fixtures", tier),
              file.path("..", "..", "..", "data-raw", "fixtures", tier)))
    if (dir.exists(p)) return(normalizePath(p))
  NULL
}

have_fixture <- function(tier = "small") {
  d <- fixture_dir(tier)
  !is.null(d) && file.exists(file.path(d, paste0(tier, ".bed")))
}

skip_without_fixture <- function(tier = "small") {
  if (!have_fixture(tier))
    testthat::skip(paste0("fixture '", tier,
                          "' not generated — run data-raw/make_fixtures.sh"))
}

fx <- function(tier, ext) file.path(fixture_dir(tier), paste0(tier, ext))

read_bim_file <- function(tier) read_bim(readLines(fx(tier, ".bim")))
read_fam_file <- function(tier) read_fam(readLines(fx(tier, ".fam")))

#' Do what the browser does: seek to each variant's offset and slice.
#'
#' Returns the concatenated raw bytes, i.e. exactly what lands in genoContent
#' before base64, so tests exercise the same path production does.
slice_bed <- function(tier, indices, n_samples) {
  bpv <- bed_bytes_per_variant(n_samples)
  con <- file(fx(tier, ".bed"), "rb")
  on.exit(close(con))
  bed_check_magic(readBin(con, "raw", 3L))
  out <- vector("list", length(indices))
  for (j in seq_along(indices)) {
    seek(con, where = 3 + (as.double(indices[j]) - 1) * bpv, origin = "start")
    out[[j]] <- readBin(con, "raw", bpv)
  }
  unlist(out)
}

#' plink --recode A output for a tier, as a matrix plus what it counted.
read_raw_oracle <- function(tier, file = paste0(tier, ".raw")) {
  raw <- utils::read.table(file.path(fixture_dir(tier), file), header = TRUE,
                           check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(raw[, -(1:6), drop = FALSE])
  storage.mode(m) <- "integer"
  list(values  = m,
       iid     = as.character(raw$IID),
       ids     = sub("_[^_]*$", "", colnames(m)),
       counted = sub("^.*_", "", colnames(m)))
}

#' base64 a raw vector or a character vector, as the browser would.
as_payload <- function(x) {
  if (is.character(x)) x <- charToRaw(paste(x, collapse = "\n"))
  base64enc::base64encode(x)
}


# ── driving the analysis ─────────────────────────────────────────────────────
# Shared by every test file: testthat sources helper-*.R globally but does not
# share definitions between test files.

# Build the three payloads the browser would send for a selection of SNP IDs.
make_payloads <- function(tier, ids) {
  bim <- read_bim_file(tier)
  fam <- read_fam_file(tier)
  idx <- match(ids, bim$id)
  idx <- idx[!is.na(idx)]

  bim_lines <- readLines(fx(tier, ".bim"))[idx]
  fam_lines <- readLines(fx(tier, ".fam"))

  list(bed = as_payload(slice_bed(tier, idx, length(fam$iid))),
       bim = as_payload(bim_lines),
       fam = as_payload(fam_lines),
       n   = length(fam$iid),
       idx = idx)
}

# Run snpImport and hand back the results. openNew needs no special-casing
# the way saveCols (OptionOutput) used to -- it is an ordinary constructor
# argument -- but this analysis never writes into the current sheet at all
# (see snpimport.b.R), so run_import does not set it: out_values below reads
# what would have been opened straight from .buildColumns(), the function
# .performOpen() itself calls, without going through option$perform() (which
# needs jmvReadWrite -- jamovi's engine only, not this project's local R
# library -- and a real session-temp directory; test-option-action.R covers
# that path where jmvReadWrite is available).
.LAST_ANALYSIS <- new.env(parent = emptyenv())

run_import <- function(p, ids, ...) {
  ns   <- asNamespace("SNPstats")
  args <- utils::modifyList(
    list(genoContent = p$geno %||% p$bed,
         variantContent = p$bim, sampleContent = p$fam,
         snpListText = paste(ids, collapse = "\n")),
    list(...))

  if (is.null(args$showSummary)) args$showSummary <- TRUE
  opts <- do.call(ns$snpImportOptions$new, args)

  a <- ns$snpImportClass$new(options = opts, data = data.frame())
  a$run()
  .LAST_ANALYSIS$value <- a
  a$results
}

# Ignores its argument's identity and reads .LAST_ANALYSIS instead, because
# the columns .performOpen would build are private state on the *analysis*,
# not on anything reachable from `results` (there is no more Output object to
# read off). Safe because every test calls this immediately after the
# matching run_import(), same as the rest of this file assumes.
#
# force(results) is load-bearing, not decoration: every caller in this suite
# nests the two, out_values(run_import(...)), and R's lazy evaluation means an
# argument nothing in the body reads is an argument never evaluated. Without
# it run_import()'s call — and the .LAST_ANALYSIS$value <- a inside it — never
# happens, and this silently reads whatever the *previous* test left there.
# Measured: every out_values(run_import(...)) call after the first in a file
# returned the previous test's columns until this was added.
out_values <- function(results) {
  force(results)
  a  <- .LAST_ANALYSIS$value
  pr <- a$.__enclos_env__$private
  d  <- pr$.cache
  if (is.null(d)) return(list())
  pr$.buildColumns(d)$values
}



`%||%` <- function(a, b) if (is.null(a)) b else a

#' Table notes as plain strings. `$notes` holds R6 Note objects, so unlist()
#' gives their methods rather than their text.
note_text <- function(tbl) vapply(tbl$notes, function(n) n$note, "")
