
# ══════════════════════════════════════════════════════════════════════════════
# pgs_weights() — read a PGS-Catalog weights file for snpPGS()
#
# snpPGS has no file-path option on purpose. Analysis options are serialised
# into the .omv file and re-run on whoever opens it, so a path option would let
# a crafted document read an arbitrary file from the opener's machine. The
# weights therefore travel as *content*: the file-browse button embeds the bytes
# into weightsContent, and this helper does the same thing for R scripting.
#
# Usage:
#   snpPGS(data = mydata, snpCols = c("rs1", "rs2"),
#          weightsContent  = pgs_weights("pgs.csv")$weightsContent,
#          weightsFilename = pgs_weights("pgs.csv")$weightsFilename)
# or, more usually:
#   w <- pgs_weights("pgs.csv")
#   do.call(snpPGS, c(list(data = mydata, snpCols = c("rs1", "rs2")), w))
# ══════════════════════════════════════════════════════════════════════════════

# Upper bound on a weights file, applied both when embedding (pgs_weights) and
# after gunzipping embedded content (.weightsRawLines). Without it a small
# crafted .gz inside a saved .omv expands without bound in the engine process.
PGS_MAX_WEIGHTS_BYTES <- 64 * 1024^2        # 64 MB

#' Read a PGS-Catalog weights file for use with snpPGS
#'
#' Reads \code{path} and returns the pair of \code{snpPGS()} arguments that
#' carry a weights file: the base64-encoded file contents and its name. The
#' file is read once, here, in the caller's own session -- nothing about the
#' path is stored in the analysis, so a saved \code{.omv} carries the weights
#' themselves rather than a path that would be re-resolved on another machine.
#'
#' A \code{.gz} file is passed through still compressed; \code{snpPGS()}
#' decompresses it based on the \code{.gz} extension in the file name.
#'
#' @param path Path to a PGS-Catalog format file (\code{.csv}, \code{.tsv},
#'   \code{.txt} or their \code{.gz} forms).
#' @return A named list with \code{weightsContent} (base64 string) and
#'   \code{weightsFilename}, suitable for splicing into a \code{snpPGS()} call.
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

  list(weightsContent  = base64enc::base64encode(path),
       weightsFilename = basename(path))
}
