
# ══════════════════════════════════════════════════════════════════════════════
# pgs_weights() — read a PGS-Catalog weights file for snpPGS()
#
# weightsFile is a native jamovi File option: jamovi's own FileSelector picks
# it and carries it as a resource inside the saved .omv, so a reopened
# analysis keeps working without the original file. This helper is the
# scripting equivalent of picking the file in the UI.
#
# Usage:
#   snpPGS(data = mydata, snpCols = c("rs1", "rs2"),
#          weightsFile = pgs_weights("pgs.csv")$weightsFile)
# or, more usually:
#   w <- pgs_weights("pgs.csv")
#   do.call(snpPGS, c(list(data = mydata, snpCols = c("rs1", "rs2")), w))
# ══════════════════════════════════════════════════════════════════════════════

# Upper bound on a weights file, applied both here and after gunzipping
# (.weightsRawLines). Without it a small crafted .gz expands without bound in
# the engine process.
PGS_MAX_WEIGHTS_BYTES <- 64 * 1024^2        # 64 MB

#' Read a PGS-Catalog weights file for use with snpPGS
#'
#' Validates \code{path} and returns the \code{snpPGS()} argument that carries
#' a weights file. jmvcore's File option accepts a plain path string from R
#' directly; this helper exists mainly to give scripted callers the same
#' size/existence checks the UI's file picker gets for free.
#'
#' A \code{.gz} file is passed through still compressed; \code{snpPGS()}
#' decompresses it based on the \code{.gz} extension in the file name.
#'
#' @param path Path to a PGS-Catalog format file (\code{.csv}, \code{.tsv},
#'   \code{.txt} or their \code{.gz} forms).
#' @return A named list with \code{weightsFile} (the path), suitable for
#'   splicing into a \code{snpPGS()} call.
#' @export
pgs_weights <- function(path) {

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("pgs_weights(): 'path' must be a single non-empty file path.")
  if (!file.exists(path))
    stop("pgs_weights(): file not found: ", path)

  size <- file.info(path)$size
  # Same ceiling the backend applies to decompressed content; a PGS catalog of
  # this size is already far larger than any real scoring file.
  if (!is.na(size) && size > PGS_MAX_WEIGHTS_BYTES)
    stop("pgs_weights(): file is larger than ",
         round(PGS_MAX_WEIGHTS_BYTES / 1024^2), " MB: ", path)

  list(weightsFile = path)
}
