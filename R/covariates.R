
# Covariate files and the merge onto genotype samples.
#
# Pure functions: text in, plain lists out. No jmvcore, no file paths.
#
# The genotype sample order is authoritative. A covariate file is matched onto
# it, never the other way round -- the .fam defines which rows exist, and once
# SNP columns have been written those rows cannot be reordered or extended
# without misaligning genotypes against phenotypes.

# Names a covariate file might use for the sample identifier. IID is the PLINK
# convention (.fam is FID/IID, and --covar expects FID IID then the covariates;
# plink2 accepts a '#IID' header alone), so it wins over the generic ones.
COV_ID_NAMES  <- c("iid", "id", "sample", "sample_id", "sampleid", "subject")
COV_FID_NAMES <- c("fid", "family", "family_id", "familyid")

#' Guess the delimiter of one line. Tab, comma and semicolon all occur in files
#' people actually have; whichever yields the most fields wins, ties to tab.
.guess_delim <- function(line) {
  best <- "\t"; best_n <- length(strsplit(line, "\t", fixed = TRUE)[[1]])
  for (d in c(",", ";")) {
    n <- length(strsplit(line, d, fixed = TRUE)[[1]])
    if (n > best_n) { best <- d; best_n <- n }
  }
  list(delim = best, n = best_n)
}

#' Drop comment lines and work out which line is the header.
#'
#' Two conventions both appear, and they disagree about what a leading '#'
#' means:
#'
#'   PGS Catalog and friends   '# free text' comments, then a plain header
#'   plink2                    the header *is* the '#' line, e.g. '#IID<TAB>age'
#'
#' They are told apart by field count rather than by guessing: a header has the
#' same number of fields as the data, a comment almost never does. The
#' delimiter is guessed from the first data line, which is the one line certain
#' to carry every field -- guessing from a free-text comment is what made a
#' perfectly good file report "only one column".
.cov_prepare <- function(lines) {
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("the covariate file is empty")

  lines <- lines[!startsWith(trimws(lines), "##")]      # VCF-style meta lines
  hashed <- startsWith(trimws(lines), "#")
  data_at <- which(!hashed)
  if (length(data_at) == 0) stop("the covariate file has no data rows")
  first <- data_at[1]

  g <- .guess_delim(lines[first])
  if (g$n < 2)
    stop("the covariate file has only one column \u2014 is the delimiter unusual? ",
         "tab, comma and semicolon are recognised")
  nf <- function(x) length(strsplit(x, g$delim, fixed = TRUE)[[1]])

  if (first > 1L) {
    cand <- sub("^[[:space:]]*#", "", lines[first - 1L])
    if (nf(cand) == g$n)
      return(list(delim = g$delim, lines = c(cand, lines[first:length(lines)])))
  }
  list(delim = g$delim, lines = lines[first:length(lines)])
}

#' A PGS weights file is per-SNP, not per-sample, so it can never be a
#' covariate table. It is an easy slot to mistake, and falling back to "column 1
#' is the ID" would match rsIDs against sample IIDs and quietly produce a table
#' of NAs, so it is named and refused instead.
.looks_like_weights <- function(header) {
  h <- tolower(trimws(header))
  any(h %in% c("effect_weight", "effect_allele", "beta", "or", "hm_effect_allele"))
}

#' Read a covariate table.
#'
#' `id_col` names the identifier column explicitly; empty means detect it. A
#' file with no recognisable header falls back to the first column, which is
#' what a headerless PLINK-style covariate file has.
#'
#' Returns list(id, fid, vars, id_name, fid_name) where `vars` is a named list
#' of raw character vectors -- typing happens in cov_typed(), so the caller can
#' report what it did.
read_covariate_table <- function(lines, id_col = "") {

  too_long <- which(nchar(lines) > 1048576L)
  if (length(too_long))
    stop("covariate file line ", too_long[1], " is too long to parse")

  prep  <- .cov_prepare(lines)
  if (length(prep$lines) < 2)
    stop("the covariate file has no data rows")

  parts <- lapply(strsplit(prep$lines, prep$delim, fixed = TRUE), trimws)
  hdr   <- parts[[1]]
  low   <- tolower(hdr)

  if (.looks_like_weights(hdr))
    stop("this looks like a PGS weights file, which lists SNPs rather than ",
         "samples \u2014 load it under 'SNPs to import' instead of as covariates")

  if (nzchar(id_col)) {
    i <- match(tolower(trimws(id_col)), low)
    if (is.na(i))
      stop("no column named '", id_col, "' in the covariate file (found: ",
           paste(hdr, collapse = ", "), ")")
  } else {
    i <- NA_integer_
    for (nm in COV_ID_NAMES) { j <- match(nm, low); if (!is.na(j)) { i <- j; break } }
    if (is.na(i)) i <- 1L                     # no header we recognise
  }

  fi <- NA_integer_
  for (nm in COV_FID_NAMES) { j <- match(nm, low); if (!is.na(j)) { fi <- j; break } }

  body <- parts[-1]
  n_f  <- length(hdr)
  col  <- function(k) vapply(body, function(r)
    if (length(r) >= k) r[k] else NA_character_, "")

  keep <- setdiff(seq_len(n_f), c(i, fi))
  if (length(keep) == 0)
    stop("the covariate file has an identifier column but no covariates")

  vars <- stats::setNames(lapply(keep, col), hdr[keep])

  list(id       = col(i),
       fid      = if (is.na(fi)) NULL else col(fi),
       vars     = vars,
       id_name  = hdr[i],
       fid_name = if (is.na(fi)) NULL else hdr[fi])
}

#' Type a raw covariate column: numeric when every non-blank value parses.
#'
#' Blanks and the usual missing markers become NA either way. A column that is
#' *mostly* numeric but has a stray label stays categorical rather than
#' silently dropping the odd value.
cov_typed <- function(x) {
  x[x %in% c("", "NA", "na", "N/A", "-9", ".")] <- NA_character_
  num <- suppressWarnings(as.numeric(x))
  if (all(is.na(num) == is.na(x))) num else factor(x)
}

#' Align a covariate table onto the genotype samples.
#'
#' A left join with the genotype order authoritative: every genotype sample
#' keeps its row, covariate rows that match nothing are dropped. Both
#' directions are reported, because an ID convention mismatch (leading zeros,
#' FID_IID versus IID) shows up as "nothing matched" and needs to be
#' diagnosable from the report alone.
#'
#' Matching uses FID+IID when the .fam's IIDs are not unique and the covariate
#' file supplies a family column; IID alone otherwise.
merge_covariates <- function(cov, fam) {

  use_pair <- !isTRUE(fam$iid_unique) && !is.null(cov$fid)

  geno_key <- if (use_pair) paste(fam$fid, fam$iid, sep = "\r") else fam$iid
  cov_key  <- if (use_pair) paste(cov$fid, cov$id, sep = "\r") else cov$id

  dup_cov <- anyDuplicated(cov_key) > 0
  idx <- match(geno_key, cov_key)

  values <- lapply(cov$vars, function(v) cov_typed(v[idx]))

  list(values      = values,
       matched     = sum(!is.na(idx)),
       n_geno      = length(geno_key),
       n_cov       = length(cov_key),
       unmatched_geno = fam$iid[is.na(idx)],
       unmatched_cov  = cov$id[!(seq_along(cov_key) %in% idx)],
       matched_on  = if (use_pair) "FID + IID" else cov$id_name,
       dup_cov     = dup_cov)
}

#' Make covariate names safe and unique against names already being emitted.
#'
#' A clash would otherwise become "age (2)" in the spreadsheet, which is silent
#' and easy to miss. Renaming here means the import report can say what it did.
#' The first of a duplicate pair keeps its name and the second is renamed, which
#' is the opposite of what comparing each name against the whole vector does.
#' The straightforward way to get that is a loop with a growing `used` vector,
#' and it was quadratic: `%in%` rebuilds its hash table on every iteration, so
#' 30 000 variants -- which the cell budget allows at 100 samples -- took 2.8 s,
#' on every option click, because this runs inside .prepared().
#'
#' Instead the whole thing is settled in one pass. `taken` names are seeded as
#' already-seen so a variant can never take one, then each name's occurrence
#' number within its own group gives the suffix. The only part that still needs
#' iterating is the rare case where the suffixed name is itself taken (a file
#' that really does contain `age` and `age_2`), and that loop runs over
#' collisions rather than over names.
cov_unique_names <- function(names, taken) {
  out <- sanitise_id(names)
  if (length(out) == 0) return(out)

  # Occurrence number of each name among the ones settled before it, counting
  # the reserved names as having been settled first.
  held <- if (length(taken)) as.character(sanitise_id(taken)) else character(0)
  all_names <- c(held, out)
  seen <- stats::ave(seq_along(all_names), all_names, FUN = seq_along)
  seen <- seen[seq_along(out) + length(held)]

  dup <- seen > 1L
  if (any(dup)) out[dup] <- paste0(out[dup], "_", seen[dup])

  # A suffix can collide with a name that was already there. Rare, so resolved
  # one at a time -- but only for the names that actually collide.
  pool <- new.env(hash = TRUE, parent = emptyenv())
  for (nm in c(held, out[!dup])) assign(nm, TRUE, envir = pool)
  for (i in which(dup)) {
    k <- seen[i]
    base <- sanitise_id(names[i])
    while (exists(out[i], envir = pool, inherits = FALSE)) {
      k <- k + 1L
      out[i] <- paste0(base, "_", k)
    }
    assign(out[i], TRUE, envir = pool)
  }
  out
}
