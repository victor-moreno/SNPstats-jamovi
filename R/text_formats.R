
# Variant-major text genotype formats: .tped/.tfam and VCF.
#
# Both are one line per variant, which is what makes them workable: the browser
# streams the file, keeps only the lines whose ID matches, and sends those. R
# never sees the variants that were not asked for.
#
# Each reader returns the same shape as the .bed path -- fam, a .bim-like
# variant table, and an n_samples x n_variants matrix of allele-2 dosages -- so
# everything downstream is format-blind.
#
# Pure functions: text in, plain lists out. No jmvcore, no file paths.

MAX_TEXT_FIELDS <- 4e6      # a single line's field count; a guard, not a policy

.split_line <- function(line) {
  f <- strsplit(trimws(line), "[ \t]+")[[1]]
  if (length(f) > MAX_TEXT_FIELDS)
    stop("a line has ", length(f), " fields \u2014 refusing to parse")
  f
}

#' Assign A1/A2 to a set of observed allele calls.
#'
#' A .tped carries no .bim, so the roles have to be derived. PLINK's convention
#' when it builds a .bim from text is A1 = the minor allele, so that is used
#' here; ties break on sort order for determinism. Missing codes never name an
#' allele (see MISSING_ALLELE).
#'
#' Where a site carries more than two allele codes this keeps the two commonest
#' and leaves the rest to be read as no-calls by the caller. It used to stop,
#' which took the whole import down over one stray letter in one sample -- and
#' the third code is usually exactly that, since a genuinely triallelic SNP is
#' rare and a typo, a lower-case call or a converter's own missing marker is
#' not.
tped_alleles <- function(calls) {
  obs <- calls[!is_missing_allele(calls)]
  if (length(obs) == 0) return(c("0", "0"))
  tab <- sort(table(obs), decreasing = TRUE)
  al  <- names(tab)
  if (length(al) == 1) return(c(al, al))       # monomorphic
  c(al[2], al[1])                              # A1 = minor, A2 = major
}

#' Two haplotype columns to an allele-2 count, anything unexpected missing.
#'
#' A call counts only when both its alleles are ones this variant actually has.
#' Everything else -- a missing code, a half call, a third allele, a stray
#' letter -- is a no-call. The alternative, comparing against A2 alone, reads an
#' unrecognised code as "not A2" and so reports it as a homozygous A1: a
#' fabricated genotype, and the one failure mode here that no downstream check
#' would ever catch.
.call_dose <- function(h1, h2, a1, a2) {
  ok <- (h1 == a1 | h1 == a2) & (h2 == a1 | h2 == a2)
  ok[is.na(ok)] <- FALSE
  d <- as.integer(h1 == a2) + as.integer(h2 == a2)
  d[!ok] <- NA_integer_
  list(dose = d,
       # Missing codes are the file saying "not typed" and are not worth
       # reporting; anything else is the file saying something we could not
       # read, which is.
       unreadable = sum(!ok & !(is_missing_allele(h1) | is_missing_allele(h2))))
}

#' Parse .tped lines with their .tfam.
#'
#' .tped: 2N+4 fields -- chr, id, cM, bp, then two allele calls per sample.
#' .tfam is byte-identical to a .fam, so read_fam handles it.
read_tped <- function(tped_lines, tfam_lines) {

  fam <- read_fam(tfam_lines)
  n   <- length(fam$iid)

  lines <- tped_lines[nzchar(trimws(tped_lines))]
  if (length(lines) == 0) stop("no .tped data lines")

  k    <- length(lines)
  dose <- matrix(NA_integer_, nrow = n, ncol = k)
  chr <- id <- a1 <- a2 <- character(k)
  bp  <- numeric(k)
  unreadable <- 0L

  for (j in seq_len(k)) {
    f <- .split_line(lines[j])
    if (length(f) != 2 * n + 4)
      stop(sprintf(".tped line %d has %d fields, expected %d (4 + 2 x %d samples)",
                   j, length(f), 2 * n + 4, n))

    chr[j] <- f[1]; id[j] <- f[2]
    p <- suppressWarnings(as.numeric(f[4]))
    if (is.na(p) || p < 0) stop(".tped line ", j, " has a bad base-pair position")
    bp[j] <- p

    calls <- f[-(1:4)]
    al <- tped_alleles(calls)
    a1[j] <- al[1]; a2[j] <- al[2]

    cd <- .call_dose(calls[c(TRUE, FALSE)], calls[c(FALSE, TRUE)], al[1], al[2])
    dose[, j] <- cd$dose
    unreadable <- unreadable + cd$unreadable
  }

  list(fam = fam,
       bim = list(chr = chr, id = id, bp = bp, a1 = a1, a2 = a2),
       dose = dose, unreadable = unreadable)
}


# -- VCF ----------------------------------------------------------------------

#' Parse VCF lines: the #CHROM header plus the variant lines that were selected.
#'
#' A VCF carries its own sample list, so no companion file is needed. Sex and
#' phenotype are not part of the format, so they come back as NA -- the .fam of a
#' PLINK trio has them, a VCF does not.
#'
#' REF becomes allele 1 and ALT allele 2, so the dosage this returns is the ALT
#' count, which is what a GT of 0/1 or 1/1 says directly.
read_vcf <- function(lines) {

  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!startsWith(lines, "##")]            # meta lines carry no data
  hdr_i <- which(startsWith(lines, "#CHROM"))
  if (length(hdr_i) == 0)
    stop("no #CHROM header line \u2014 the VCF header must be included")

  hdr <- .split_line(sub("^#", "", lines[hdr_i[1]]))
  if (length(hdr) < 10)
    stop("the VCF header lists no samples")
  iid <- hdr[10:length(hdr)]
  n   <- length(iid)

  body <- lines[-seq_len(hdr_i[1])]
  body <- body[!startsWith(body, "#")]
  if (length(body) == 0) stop("no VCF variant lines")

  keep <- rep(TRUE, length(body))
  unreadable <- 0L
  chr <- id <- a1 <- a2 <- character(length(body))
  bp  <- numeric(length(body))
  dose <- matrix(NA_integer_, nrow = n, ncol = length(body))

  for (j in seq_along(body)) {
    f <- .split_line(body[j])
    if (length(f) != n + 9)
      stop(sprintf("VCF line %d has %d fields, expected %d (9 + %d samples)",
                   j, length(f), n + 9, n))

    # Multi-allelic sites and indels do not fit a biallelic genotype; they are
    # skipped and counted rather than silently mis-decoded.
    if (grepl(",", f[5], fixed = TRUE) ||
        nchar(f[4]) != 1 || nchar(f[5]) != 1) {
      keep[j] <- FALSE
      next
    }

    chr[j] <- f[1]; id[j] <- f[3]; a1[j] <- f[4]; a2[j] <- f[5]
    p <- suppressWarnings(as.numeric(f[2]))
    if (is.na(p) || p < 0) stop("VCF line ", j, " has a bad POS")
    bp[j] <- p

    # GT is not reliably the first sub-field: its position comes from FORMAT,
    # per line.
    fmt <- strsplit(f[9], ":", fixed = TRUE)[[1]]
    gt_at <- match("GT", fmt)
    if (is.na(gt_at))
      stop("VCF line ", j, " has no GT in its FORMAT field")

    smp <- f[10:(n + 9)]
    gt <- if (length(fmt) == 1) smp else
      vapply(strsplit(smp, ":", fixed = TRUE),
             function(x) if (length(x) >= gt_at) x[gt_at] else NA_character_, "")

    dose[, j] <- .vcf_gt_dose(gt)
    # A GT of './.' is the file saying "not called"; anything else that came
    # back NA is the file saying something this reader could not read, and the
    # import report says how much of that there was.
    unreadable <- unreadable +
      sum(is.na(dose[, j]) & !is.na(gt) & !grepl(".", gt, fixed = TRUE))
  }

  if (!any(keep)) stop("no biallelic SNVs among the selected VCF lines")

  list(fam = list(fid = iid, iid = iid,
                  sex   = factor(rep(NA_character_, n), levels = c("male", "female")),
                  pheno = factor(rep(NA_character_, n), levels = c("control", "case")),
                  fid_informative = FALSE, iid_unique = !anyDuplicated(iid)),
       bim = list(chr = chr[keep], id = id[keep], bp = bp[keep],
                  a1 = a1[keep], a2 = a2[keep]),
       dose = dose[, keep, drop = FALSE],
       skipped = sum(!keep), unreadable = unreadable)
}

#' GT strings to an ALT-allele count. Handles phased and unphased separators,
#' and haploid calls, which appear on the sex chromosomes.
#'
#' Only 0 and 1 are alleles here: the line has already been checked to be
#' biallelic, so an index of 2 or more does not belong to it. Comparing against
#' "1" alone counted such an index as a REF copy, which turned a call this
#' reader cannot represent into a confident homozygous reference.
.vcf_gt_dose <- function(gt) {
  a   <- substr(gt, 1L, 1L)
  sep <- substr(gt, 2L, 2L)
  b   <- substr(gt, 3L, 3L)                  # '' when the call is haploid
  hap <- !is.na(b) & !nzchar(b)

  ok <- !is.na(a) & (a == "0" | a == "1") &
        (hap | (!is.na(b) & (b == "0" | b == "1") &
                (sep == "/" | sep == "|")))
  ok[is.na(ok)] <- FALSE

  d <- as.integer(a == "1") + ifelse(hap, 0L, as.integer(b == "1"))
  d[!ok] <- NA_integer_
  d
}


# -- .ped / .map --------------------------------------------------------------

#' Parse .map lines: chr, id, cM (optional), bp.
#'
#' Three or four columns -- PLINK writes the cM column but it is optional, and
#' a three-column .map puts the base-pair position where cM would be.
read_map <- function(lines) {
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("no .map lines")
  p <- lapply(lines, .split_line)
  nf <- lengths(p)
  if (any(nf < 3)) stop(".map line ", which(nf < 3)[1], " has fewer than 3 fields")

  bp_at <- if (all(nf >= 4)) 4L else 3L
  bp <- suppressWarnings(as.numeric(vapply(p, `[`, "", bp_at)))
  if (any(is.na(bp) | bp < 0))
    stop(".map contains a non-numeric or negative base-pair position")

  list(chr = vapply(p, `[`, "", 1L), id = vapply(p, `[`, "", 2L), bp = bp)
}

#' Parse .ped lines against their .map.
#'
#' A .ped is sample-major: one line per sample, 2V+6 fields, the first six of
#' which are exactly a .fam row. That is why no companion sample file is
#' needed -- unlike .tped, the .ped carries its own IDs, sex and phenotype.
#'
#' The lines are expected to be *already reduced* to the selected variants, in
#' the order of `map_lines`, because selecting from a sample-major format means
#' dropping columns from every line and that is done before the data is sent.
read_ped <- function(ped_lines, map_lines) {

  map <- read_map(map_lines)
  k   <- length(map$id)

  lines <- ped_lines[nzchar(trimws(ped_lines))]
  if (length(lines) == 0) stop("no .ped data lines")
  n <- length(lines)

  parts <- lapply(lines, .split_line)
  nf <- lengths(parts)
  bad <- which(nf != 2 * k + 6)
  if (length(bad))
    stop(sprintf(".ped line %d has %d fields, expected %d (6 + 2 x %d variants)",
                 bad[1], nf[bad[1]], 2 * k + 6, k))

  # One matrix rather than per-variant list indexing: a variant is then two
  # whole columns, which keeps this linear in cells instead of quadratic.
  m <- matrix(unlist(parts), nrow = n, byrow = TRUE)

  fam <- read_fam(vapply(seq_len(n), function(i)
    paste(m[i, 1:6], collapse = " "), ""))

  dose <- matrix(NA_integer_, nrow = n, ncol = k)
  a1 <- a2 <- character(k)
  unreadable <- 0L

  for (j in seq_len(k)) {
    h1 <- m[, 6 + 2 * j - 1]
    h2 <- m[, 6 + 2 * j]
    al <- tped_alleles(c(h1, h2))
    a1[j] <- al[1]; a2[j] <- al[2]

    cd <- .call_dose(h1, h2, al[1], al[2])
    dose[, j] <- cd$dose
    unreadable <- unreadable + cd$unreadable
  }

  list(fam = fam,
       bim = list(chr = map$chr, id = map$id, bp = map$bp, a1 = a1, a2 = a2),
       dose = dose, unreadable = unreadable)
}
