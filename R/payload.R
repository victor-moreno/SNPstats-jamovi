
# The transport layer: turning what the browser sent into usable data.
#
# Files reach the analysis only as base64 in hidden String options -- the only
# channel available, and the widest one jamovi has. This file decodes that
# channel and nothing else; it is separate from plink_read.R so the parsers stay
# free of any notion of how the bytes arrived.
#
# Everything here is defensive. A payload is attacker-controlled input in a
# process that also holds other users' sessions in a cloud deployment.

# Refuse to decode more than this from one option, whatever it claims to be.
# The transport cannot deliver more than 4 MiB anyway, so anything larger means
# something is wrong rather than large.
MAX_PAYLOAD_BYTES <- 8 * 1024^2

#' Decode one base64 option into raw bytes, transparently gunzipping.
#'
#' The browser may gzip before base64 (CompressionStream), which is the only
#' lever on the transport ceiling. Detection is by magic bytes rather than by
#' trusting a flag, so a mismatch between what the UI thinks it sent and what
#' arrived surfaces here instead of as corrupt genotypes later.
payload_bytes <- function(b64, what = "payload") {
  if (!is.character(b64) || length(b64) != 1L || !nzchar(b64))
    return(raw(0))

  n <- nchar(b64)
  if (n %% 4L != 0L)
    stop(what, " is truncated (base64 length ", n, " is not a multiple of 4)")
  if (n / 4 * 3 > MAX_PAYLOAD_BYTES)
    stop(what, " is ", round(n / 4 * 3 / 1024^2, 1),
         " MB \u2014 refusing to decode more than ", MAX_PAYLOAD_BYTES / 1024^2, " MB")

  # base64enc::base64decode does NOT reject characters outside the alphabet --
  # it silently returns garbage -- so the alphabet is checked here. Without this
  # a corrupted option becomes plausible-looking wrong genotypes.
  if (!grepl("^[A-Za-z0-9+/\n\r]*={0,2}$", b64))
    stop(what, " is not valid base64")

  raw_bytes <- tryCatch(base64enc::base64decode(b64),
                        error = function(e)
                          stop(what, " is not valid base64"))

  # The browser's CompressionStream('gzip') emits a real gzip stream (1f 8b),
  # and that is the only compressed form this transport carries.
  is_gzip <- length(raw_bytes) >= 2 &&
             raw_bytes[1] == as.raw(0x1f) && raw_bytes[2] == as.raw(0x8b)

  # A zlib stream (78 xx) is what R's own memCompress(type = "gzip") writes, and
  # it used to be accepted for symmetry. It is refused now: it cannot be inflated
  # under a bound (gzcon() does not recognise the header and passes the bytes
  # through untouched, silently), so accepting it means memDecompress() and an
  # allocation sized by the attacker. Nothing on the wire is ever zlib, so this
  # costs no caller anything.
  if (length(raw_bytes) >= 2 && raw_bytes[1] == as.raw(0x78) &&
      as.integer(raw_bytes[2]) %in% c(0x01, 0x5e, 0x9c, 0xda))
    stop(what, " is zlib-compressed; only gzip is accepted")

  if (is_gzip)
    raw_bytes <- .gunzip_bounded(raw_bytes, MAX_PAYLOAD_BYTES, what)

  raw_bytes
}

# Inflate a gzip member, never allocating more than `cap` bytes of output.
#
# memDecompress() materialises the whole output in one allocation, so testing
# its length afterwards tests memory that has already been taken: 424 KB of
# base64 -- comfortably inside every transport ceiling -- expanded to 400 MB
# before the length could be objected to, and five payloads arrive per request.
# In a shared engine that is a one-message denial of service.
#
# gzcon() over a rawConnection() is base R's streaming inflater, so readBin()
# with n = cap + 1 stops there: the peak allocation is the cap, whatever the
# stream claims or contains. Reading one byte past the cap is what distinguishes
# "exactly at the limit" from "over it".
.gunzip_bounded <- function(b, cap, what) {
  rc <- rawConnection(b, "rb")
  on.exit(try(suppressWarnings(close(rc)), silent = TRUE), add = TRUE)
  con <- gzcon(rc)
  on.exit(try(suppressWarnings(close(con)), silent = TRUE), add = TRUE)

  out <- tryCatch(readBin(con, "raw", n = cap + 1),
                  error = function(e)
                    stop(what, " looks compressed but will not decompress"))

  if (length(out) > cap)
    stop(what, " expands to more than ", cap / 1024^2,
         " MB \u2014 refusing (possible decompression bomb)")

  # A corrupt stream is not an error to gzcon(): it stops where it stops and the
  # partial contents would parse as a file simply missing its last variants. The
  # checksum mismatch it notices goes to stderr through REprintf, so no handler
  # can see it and the read returns as if nothing happened.
  #
  # ISIZE is the evidence that is left. It is a number in the payload, so it is
  # no use as a size guard -- the cap above is that, and is what an ISIZE that
  # lies upward runs into -- but comparing it against what actually came out
  # catches the short read.
  claimed <- .gzip_isize(b)
  if (!is.na(claimed) && claimed != length(out))
    stop(what, " is truncated or corrupt (declares ", claimed,
         " bytes of contents, ", length(out), " decompressed)")

  out
}

# How large a gzip member claims its contents are: ISIZE, the last four bytes,
# little-endian and taken mod 2^32. NA when the bytes are not a gzip member or
# are too short to carry one.
#
# Untrusted by construction -- it is a number in the payload -- so it never
# sizes an allocation; it is compared against what was actually inflated.
.gzip_isize <- function(b) {
  n <- length(b)
  if (n < 18 || b[1] != as.raw(0x1f) || b[2] != as.raw(0x8b)) return(NA_real_)
  sum(as.double(as.integer(b[(n - 3):n])) * c(1, 256, 65536, 16777216))
}

#' Decode a base64 option into text lines.
#'
#' Handles the three line endings and drops a UTF-8 BOM, both of which appear
#' in real files often enough that failing on them would be a support burden.
payload_lines <- function(b64, what = "file") {
  b <- payload_bytes(b64, what)
  if (length(b) == 0) return(character(0))

  if (length(b) >= 3 && all(b[1:3] == as.raw(c(0xef, 0xbb, 0xbf))))
    b <- b[-(1:3)]

  txt <- rawToChar(b)
  Encoding(txt) <- "UTF-8"
  txt <- gsub("\r\n", "\n", txt, fixed = TRUE)
  txt <- gsub("\r",   "\n", txt, fixed = TRUE)
  strsplit(txt, "\n", fixed = TRUE)[[1]]
}


# -- the SNP selection list ---------------------------------------------------

#' Header names that can hold the variant ID, best first.
#'
#' The order is the contract, not just a list: snpimport.js reads the same names
#' out of the same file to decide which variants to slice, and it walks them in
#' this order. Matching by position instead -- which this did -- made the two
#' disagree on a header like `id, rsid, effect_allele`: the browser sliced the
#' rsIDs while R reported the contents of `id` as requested and every one of
#' them as not found. The import worked and the report was nonsense.
SNP_ID_NAMES <- c("rsid", "variant_id", "snpid", "snp", "id", "name")

#' First of `names` that appears in `hdr`, by preference rather than by column
#' position. Returns integer(0) when none does.
.pick_col <- function(hdr, names) {
  for (nm in names) {
    at <- match(nm, hdr)
    if (!is.na(at)) return(at)
  }
  integer(0)
}

#' Extract variant IDs from a pasted list, a plain ID file, or a weights file.
#'
#' A PGS-Catalog weights file is a first-class way in, because scoring a
#' published PGS is a main use case, and it carries the effect allele, which
#' enables the orientation check that plain ID lists cannot support. The three
#' shapes are told apart by looking at the content rather than asking the user,
#' since getting it wrong is silent.
#'
#' Returns a list(ids, effect_allele, weight) -- the last two NULL unless a
#' weights file was recognised.
parse_snp_selection <- function(lines) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0)
    return(list(ids = character(0), effect_allele = NULL, weight = NULL))

  # PGS Catalog files carry '#' comment headers; drop them before looking for
  # the column header.
  body <- lines[!startsWith(lines, "#")]
  if (length(body) == 0)
    return(list(ids = character(0), effect_allele = NULL, weight = NULL))

  hdr    <- tolower(strsplit(body[1], "[\t,;]")[[1]])
  hdr    <- trimws(hdr)
  id_col <- .pick_col(hdr, SNP_ID_NAMES)
  ea_col <- .pick_col(hdr, c("effect_allele", "ea", "a1", "risk_allele"))

  if (length(id_col) == 0) {
    # No recognisable header: treat every whitespace-separated token as an ID.
    ids <- unlist(strsplit(body, "[ \t,;]+"))
    return(list(ids = ids[nzchar(ids)], effect_allele = NULL, weight = NULL))
  }

  w_col <- .pick_col(hdr, c("effect_weight", "weight", "beta"))
  rows  <- strsplit(body[-1], "[\t,;]")
  get   <- function(k) if (length(k)) trimws(vapply(rows, function(r)
                          if (length(r) >= k[1]) r[k[1]] else NA_character_, "")) else NULL

  ids <- get(id_col)
  keep <- !is.na(ids) & nzchar(ids)
  wt <- get(w_col)

  list(ids           = ids[keep],
       effect_allele = if (!is.null(get(ea_col))) toupper(get(ea_col))[keep] else NULL,
       weight        = if (!is.null(wt)) suppressWarnings(as.numeric(wt))[keep] else NULL)
}


# -- allele orientation -------------------------------------------------------

#' Compare a weights file's effect allele against the .bim alleles.
#'
#' Silent allele flips are the commonest source of wrong PGS results, and the
#' check is nearly free once the effect allele is in hand, so it belongs in the
#' import summary rather than in the user's head. Never auto-flips: a wrong
#' flip is unrecoverable downstream, so this reports and lets the user decide.
#'
#' Returns one of: "match" (effect allele is a2, the dosage we count),
#' "swapped" (it is a1 -- dosage counts the other allele), "flipped" (matches
#' only after complementing), "ambiguous" (A/T or C/G, where strand cannot be
#' resolved from alleles alone), "absent".
#' Frequency of the effect allele, or NA when it is not one of the two observed.
#'
#' A PGS weight is signed for a named allele, so "MAF = 0.12" is the wrong
#' number to check a score against: what matters is how common *that* allele is,
#' which is 0.88 as often as it is 0.12. Complements count -- a weights file on
#' the opposite strand still names the same allele -- and the caller reports the
#' strand question separately through allele_match(), so this only has to answer
#' the frequency.
effect_freq <- function(effect, a1, a2, freq_a2) {
  if (is.na(effect) || is.na(freq_a2)) return(NA_real_)
  e <- toupper(trimws(effect)); x1 <- toupper(a1); x2 <- toupper(a2)
  if (!nzchar(e)) return(NA_real_)

  if (e == x2) return(freq_a2)
  if (e == x1) return(1 - freq_a2)

  ce <- unname(c(A = "T", T = "A", C = "G", G = "C")[e])
  if (!is.na(ce)) {
    if (ce == x2) return(freq_a2)
    if (ce == x1) return(1 - freq_a2)
  }
  NA_real_
}

allele_match <- function(effect, a1, a2) {
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  e <- toupper(effect); x1 <- toupper(a1); x2 <- toupper(a2)

  if (is.na(e) || !nzchar(e)) return("absent")

  palindromic <- (x1 == "A" && x2 == "T") || (x1 == "T" && x2 == "A") ||
                 (x1 == "C" && x2 == "G") || (x1 == "G" && x2 == "C")

  if (e == x2) return(if (palindromic) "ambiguous" else "match")
  if (e == x1) return(if (palindromic) "ambiguous" else "swapped")

  ce <- comp[e]
  if (!is.na(ce)) {
    if (ce == x2) return("flipped")
    if (ce == x1) return("flipped")
  }
  "absent"
}
