
# PLINK format readers -- pure functions, no jmvcore, no R6, no file paths.
#
# Everything here takes bytes or character vectors and returns plain lists, so
# it is testable without jamovi and auditable on its own. Two deliberate
# absences:
#
#   * no `file()`, `readLines(<path>)` or any other filesystem access. The
#     module never opens a path, by design; bytes arrive from the browser.
#     Tests read fixtures through their own helper, which calls the decoders
#     here but does its own seeking.
#   * no `eval`, `parse`, `source` or `system`.
#
# Every function treats its input as hostile: sizes are checked against each
# other before anything is allocated, and nothing is sized from a number that
# came out of the data without being reconciled with the actual byte count.


# -- limits -------------------------------------------------------------------
# Structural sanity bounds. These are not the transport budget (that is
# enforced in the browser, see snpImport.js) -- they exist so a corrupt or
# hostile file cannot make us allocate absurdly before we notice.

MAX_SAMPLES  <- 5e6
MAX_VARIANTS <- 1e6
MAX_CELLS    <- 5e7      # far above the ~3M transport cap; a backstop only


# -- .bed decoding ------------------------------------------------------------

# .bed packs four samples per byte, two bits each, sample 1 in the low-order
# bits. Codes (PLINK 1.9 spec): 00 = hom allele 1, 01 = MISSING, 10 = het,
# 11 = hom allele 2. Note 01 = missing, which is not the intuitive ordering and
# is the classic off-by-one in home-grown readers.
#
# lut[b + 1, j] is the dosage of allele 2 for the j-th sample in a byte of
# value b, NA for missing. Built once; decoding is then a single vectorised
# index over the raw bytes.

.plink_lut <- function() {
  b <- 0:255
  code <- function(j) bitwAnd(bitwShiftR(b, 2L * j), 3L)
  dose <- function(cd) {
    out <- integer(length(cd))
    out[cd == 0L] <- 0L    # hom allele 1
    out[cd == 1L] <- NA    # missing
    out[cd == 2L] <- 1L    # het
    out[cd == 3L] <- 2L    # hom allele 2
    out
  }
  matrix(c(dose(code(0L)), dose(code(1L)), dose(code(2L)), dose(code(3L))),
         nrow = 256L, ncol = 4L)
}

PLINK_LUT <- .plink_lut()

#' Bytes per variant block for n samples.
bed_bytes_per_variant <- function(n_samples) ceiling(n_samples / 4)

#' Decode one variant block into n_samples allele-2 dosages (0/1/2, NA).
bed_decode_block <- function(block, n_samples) {
  bpv <- bed_bytes_per_variant(n_samples)
  if (length(block) != bpv)
    stop(sprintf("genotype block is %d bytes, expected %d for %d samples",
                 length(block), bpv, n_samples))
  g <- PLINK_LUT[as.integer(block) + 1L, , drop = FALSE]
  as.vector(t(g))[seq_len(n_samples)]     # drop the padding in the last byte
}

#' Decode concatenated variant blocks into an n_samples x n_variants matrix.
#'
#' `raw_bytes` is the browser's concatenation of the selected slices, in the
#' same order as the .bim lines it sent. The length check is the structural
#' guarantee that the two agree -- a mismatch means the trio does not belong
#' together, or the selection and the slices got out of step.
bed_decode_slices <- function(raw_bytes, n_samples, n_variants) {
  stopifnot(is.raw(raw_bytes))
  if (n_samples < 1 || n_samples > MAX_SAMPLES)
    stop("implausible sample count: ", n_samples)
  if (n_variants < 1 || n_variants > MAX_VARIANTS)
    stop("implausible variant count: ", n_variants)
  if (as.double(n_samples) * n_variants > MAX_CELLS)
    stop("refusing to decode ", format(as.double(n_samples) * n_variants,
                                       big.mark = ","), " genotypes")

  bpv      <- bed_bytes_per_variant(n_samples)
  expected <- as.double(n_variants) * bpv     # double: overflows int on real data
  if (length(raw_bytes) != expected)
    stop(sprintf(paste("genotype data is %s bytes, expected %s",
                       "(%d variants x %d bytes for %d samples)"),
                 format(length(raw_bytes), big.mark = ","),
                 format(expected, big.mark = ","), n_variants, bpv, n_samples))

  out <- matrix(NA_integer_, nrow = n_samples, ncol = n_variants)
  for (j in seq_len(n_variants)) {
    from <- (j - 1) * bpv + 1
    out[, j] <- bed_decode_block(raw_bytes[from:(from + bpv - 1)], n_samples)
  }
  out
}

#' Re-check that a .bed, its .bim and its .fam describe each other.
#'
#' `dims` is "bedBytes,bimLines,famLines" as the browser measured the three
#' files, and this is the only place their mutual consistency can be tested. The
#' payload cannot show it: the browser computes variant i's offset as
#' 3 + i * ceil(n/4) from the .fam's sample count, and bed_decode_slices()
#' re-derives its expected length from that same count, so a .fam from a
#' different build of the dataset agrees with itself and simply decodes the
#' wrong bytes. A whole .bed states its own size, and nothing in the other two
#' files can fake it.
#'
#' Empty `dims` means the payload predates this check (an .omv saved by an
#' earlier version). That is reported by the caller, not refused: the file is not
#' necessarily wrong, it is only unverified.
bed_check_dims <- function(dims, n_fam) {
  if (!is.character(dims) || length(dims) != 1L || !nzchar(dims))
    return(invisible(FALSE))

  p <- suppressWarnings(as.numeric(strsplit(dims, ",", fixed = TRUE)[[1]]))
  if (length(p) != 3L || any(is.na(p)) || any(p < 0))
    stop("the genotype dimensions are unreadable")

  bed_bytes <- p[1]; bim_lines <- p[2]; fam_lines <- p[3]

  # The .fam that was measured has to be the .fam that arrived, or the rest of
  # the arithmetic is about a different file.
  if (fam_lines != n_fam)
    stop(sprintf(paste("the .fam carries %d samples but the browser measured %d",
                       "\u2014 the sample file changed between reading and sending"),
                 n_fam, fam_lines))

  expected <- 3 + bim_lines * bed_bytes_per_variant(n_fam)
  if (bed_bytes != expected)
    stop(sprintf(paste("the .bed is %s bytes, but its .bim (%s variants) and",
                       ".fam (%s samples) describe %s \u2014 these files are not",
                       "from the same dataset"),
                 format(bed_bytes, big.mark = ",", scientific = FALSE),
                 format(bim_lines, big.mark = ",", scientific = FALSE),
                 format(n_fam, big.mark = ",", scientific = FALSE),
                 format(expected, big.mark = ",", scientific = FALSE)))
  invisible(TRUE)
}

#' Validate a .bed header and report whether the rest is a whole file.
#'
#' Not reached in production: the browser sends variant blocks, never the header,
#' so the magic bytes are checked there (_magic in snpimport.js) and R never sees
#' them. It is kept because the test helper reads whole fixture files through it,
#' and because it is the one part of "R re-validates everything it receives" that
#' R structurally cannot do for the .bed path -- which is why bed_check_dims()
#' exists.
bed_check_magic <- function(first3) {
  if (length(first3) < 3)
    stop("not a .bed file: fewer than 3 bytes")
  m <- as.integer(first3[1:3])
  if (!identical(m[1:2], c(0x6cL, 0x1bL)))
    stop("not a .bed file: wrong magic bytes")
  if (m[3] == 0x00L)
    stop(paste("this .bed is in the legacy sample-major layout.",
               "Convert it first:  plink --bfile <stem> --make-bed --out <stem>"))
  if (m[3] != 0x01L)
    stop("unrecognised .bed layout byte: ", m[3])
  invisible(TRUE)
}


# -- text sidecar files -------------------------------------------------------

# Split on any run of whitespace, which covers both the space-delimited .fam
# and the tab-delimited .bim without having to know which we were given.
.split_fields <- function(lines, expect, what) {
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) stop("empty ", what)

  # Bound the line length before splitting: one pathological line is a denial
  # of service against the engine, and no legitimate .bim or .fam line is
  # anywhere near this.
  too_long <- which(nchar(lines) > 65536L)
  if (length(too_long))
    stop(what, " line ", too_long[1], " is ", nchar(lines[too_long[1]]),
         " characters \u2014 refusing to parse")

  parts <- strsplit(trimws(lines), "[ \t]+")
  nf    <- lengths(parts)
  bad   <- which(nf < expect)
  if (length(bad))
    stop(what, " line ", bad[1], " has ", nf[bad[1]],
         " fields, expected at least ", expect)
  parts
}

#' Parse .bim lines: chr, variant ID, cM, bp, allele 1, allele 2.
read_bim <- function(lines) {
  p <- .split_fields(lines, 6L, ".bim")
  bp <- suppressWarnings(as.numeric(vapply(p, `[`, "", 4L)))
  if (any(is.na(bp) | bp < 0))
    stop(".bim contains a non-numeric or negative base-pair position")

  list(
    chr = vapply(p, `[`, "", 1L),
    id  = vapply(p, `[`, "", 2L),
    bp  = bp,
    a1  = vapply(p, `[`, "", 5L),
    a2  = vapply(p, `[`, "", 6L)
  )
}

#' Parse .fam lines: FID, IID, father, mother, sex, phenotype.
#'
#' Sex and phenotype are mapped to factors here rather than left numeric,
#' because that is what makes them usable in snpStats without recoding
#' without recoding. Sex: 1 = male, 2 = female, 0/other = NA. Phenotype:
#' 1 = control, 2 = case, -9/0/non-numeric = NA.
read_fam <- function(lines) {
  p <- .split_fields(lines, 6L, ".fam")

  sex_raw <- suppressWarnings(as.integer(vapply(p, `[`, "", 5L)))
  sex <- factor(ifelse(sex_raw == 1L, "male",
                ifelse(sex_raw == 2L, "female", NA_character_)),
                levels = c("male", "female"))

  ph_raw <- suppressWarnings(as.numeric(vapply(p, `[`, "", 6L)))
  pheno <- factor(ifelse(!is.na(ph_raw) & ph_raw == 1, "control",
                  ifelse(!is.na(ph_raw) & ph_raw == 2, "case", NA_character_)),
                  levels = c("control", "case"))

  fid <- vapply(p, `[`, "", 1L)
  iid <- vapply(p, `[`, "", 2L)

  # A .fam identifies a sample by the *pair* (FID, IID) -- the spec only
  # requires IID to be unique within a family. So the pair is what identifies a
  # row, and whether FID carries information is worth knowing: when every FID
  # is the same it is noise, when IIDs repeat it is essential.
  list(
    fid       = fid,
    iid       = iid,
    sex       = sex,
    pheno     = pheno,
    fid_informative = length(unique(fid)) > 1L,
    iid_unique      = !anyDuplicated(iid)
  )
}


# -- genotypes and per-SNP statistics -----------------------------------------

#' Allele codes that name no allele.
#'
#' PLINK's own missing code is `0`; VCF uses `.`; converters and hand-edited
#' files use `-`, `N` (IUPAC "unknown base") and the odd literal `NA`. None of
#' them is a nucleotide, so treating them as one would invent a genotype that
#' the file does not claim. Case is ignored: `n` occurs in files written by
#' tools that lower-case their sequence.
MISSING_ALLELE <- c("0", ".", "-", "n", "na", "")

is_missing_allele <- function(x) is.na(x) | tolower(trimws(x)) %in% MISSING_ALLELE

#' Dosages to a genotype factor.
#'
#' dose counts allele 2, so 0 = a1/a1, 1 = a1/a2, 2 = a2/a2. Levels are always
#' all three, present or not, so a table of counts never silently loses a
#' genotype class.
#'
#' The two alleles are labelled in sorted order rather than in the order the
#' source file happens to name them. Which allele a format calls "first" is
#' arbitrary and differs between them -- a VCF names REF first, a .bim picks by
#' frequency -- so without this the *same* heterozygote reads as "A/G" from one
#' file and "G/A" from another, and the two imports would not compare.
geno_factor <- function(dose, a1, a2) {
  # A genotype that needs an allele the variant file does not name cannot be
  # written down. It happens in a .bim, whose A1 is `0` when nothing carries the
  # minor allele: a well-formed one then has no sample with a copy of it, but a
  # malformed one does, and "0/0" or "0/A" would read as a real genotype rather
  # than as the absence of one. The dosage is left alone -- the number of copies
  # of allele 2 is a fact whatever the other allele turns out to be called.
  # !is.na first: a logical index carrying NA is an error in an assignment, and
  # a missing call is exactly what dose holds NA for.
  if (is_missing_allele(a1)) dose[!is.na(dose) & dose < 2L] <- NA_integer_
  if (is_missing_allele(a2)) dose[!is.na(dose) & dose > 0L] <- NA_integer_

  if (a1 <= a2) {
    lo <- a1; hi <- a2; d <- dose
  } else {
    lo <- a2; hi <- a1; d <- 2L - dose      # relabelled, so the count flips
  }
  lv <- c(paste0(lo, "/", lo), paste0(lo, "/", hi), paste0(hi, "/", hi))
  # unique(): a monomorphic variant read from text has one observed allele, so
  # tped_alleles() names it as both, and all three labels collapse to the same
  # string. Duplicated levels are an error in factor(), which took the whole
  # import down over a SNP that simply did not vary. A .bim writes '0' for the
  # absent allele and so never hits it.
  factor(lv[d + 1L], levels = unique(lv))
}

#' Per-SNP statistics for every column at once.
#'
#' The same quantities as snp_stats(), computed with four column-wise passes
#' over the matrix instead of one call per SNP. Measured at 250 SNPs x 10000
#' samples: 0.076 s down to 0.014 s. Parallelism was considered and rejected --
#' the whole computation is under 0.1 s, less than the startup cost of a PSOCK
#' cluster, and mclapply does not work on Windows. Vectorising is faster, has no
#' dependency and behaves identically on every platform.
#'
#' test-plink-read.R checks this agrees with snp_stats() column by column, so
#' the two cannot drift apart.
#' `hwe_rows` restricts the Hardy-Weinberg test to a subset of samples while
#' every other statistic still describes all of them. The usual subset is the
#' controls: departure from equilibrium in cases can be a sign of association
#' rather than of a genotyping problem, which is why PLINK's --hwe also looks
#' at controls only by default for case/control phenotypes. NULL uses everyone.
snp_stats_all <- function(dose, hwe_rows = NULL) {
  n_miss <- colSums(is.na(dose))
  n11    <- colSums(dose == 0L, na.rm = TRUE)
  n12    <- colSums(dose == 1L, na.rm = TRUE)
  n22    <- colSums(dose == 2L, na.rm = TRUE)
  n      <- n11 + n12 + n22

  freq2 <- ifelse(n > 0, (n12 + 2 * n22) / (2 * n), NA_real_)
  maf   <- pmin(freq2, 1 - freq2)

  if (is.null(hwe_rows)) {
    h11 <- n11; h12 <- n12; h22 <- n22
    hwe_n <- nrow(dose)
  } else {
    sub <- dose[hwe_rows, , drop = FALSE]
    h11 <- colSums(sub == 0L, na.rm = TRUE)
    h12 <- colSums(sub == 1L, na.rm = TRUE)
    h22 <- colSums(sub == 2L, na.rm = TRUE)
    hwe_n <- nrow(sub)
  }

  hwe <- vapply(seq_along(h11), function(j)
    if (h11[j] + h12[j] + h22[j] > 0)
      tryCatch(hwe_exact_p(h11[j], h12[j], h22[j]),
               error = function(e) NA_real_) else NA_real_,
    numeric(1))

  list(n = n, missing = n_miss, n11 = n11, n12 = n12, n22 = n22,
       freq_a2 = freq2, maf = maf, hwe_p = hwe, hwe_n = hwe_n)
}

#' Per-individual missingness and heterozygosity, for sample-level QC.
#'
#' Both are computed over the *selected* SNPs only, which is what we have.
#' That makes them differ from a PLINK run over a full panel, and heterozygosity
#' in particular is noisy on a small selection -- the caller warns about it.
sample_stats <- function(dose) {
  n_snp  <- ncol(dose)
  n_miss <- rowSums(is.na(dose))
  n_obs  <- n_snp - n_miss
  n_het  <- rowSums(dose == 1L, na.rm = TRUE)

  list(missing_rate = if (n_snp > 0) n_miss / n_snp else rep(0, nrow(dose)),
       het_rate     = ifelse(n_obs > 0, n_het / n_obs, NA_real_),
       n_obs        = n_obs)
}

#' Per-SNP descriptive statistics from a dosage vector.
#'
#' MAF is the minor-allele frequency, so it is reported for whichever allele is
#' rarer, and `minor` names it. HWE comes from the vendored exact test.
snp_stats <- function(dose) {
  ok <- !is.na(dose)
  n  <- sum(ok)
  d  <- dose[ok]

  n11 <- sum(d == 0L)    # hom allele 1
  n12 <- sum(d == 1L)
  n22 <- sum(d == 2L)    # hom allele 2

  freq2 <- if (n > 0) (n12 + 2 * n22) / (2 * n) else NA_real_
  maf   <- if (is.na(freq2)) NA_real_ else min(freq2, 1 - freq2)

  hwe <- if (n > 0) tryCatch(hwe_exact_p(n11, n12, n22),
                             error = function(e) NA_real_) else NA_real_

  list(n = n, missing = sum(!ok), n11 = n11, n12 = n12, n22 = n22,
       freq_a2 = freq2, maf = maf, minor_is_a2 = !is.na(freq2) && freq2 <= 0.5,
       hwe_p = hwe)
}


# -- identifiers --------------------------------------------------------------

#' Make a variant ID safe to use as a jamovi column name.
#'
#' jamovi has its own `fix_column_names`, but doing it ourselves means the name
#' in the summary table and the name in the spreadsheet are the same string,
#' and that control characters never reach either.
sanitise_id <- function(x) {
  # ifelse() on an empty vector returns logical(0), not character(0), and the
  # callers go on to c() the result with names -- so an empty set of covariates
  # or variants would change the type of the whole vector. Type-stable in, type-
  # stable out.
  if (length(x) == 0) return(character(0))
  x <- gsub("[[:cntrl:]]", "", x)
  x <- gsub("`", "_", x)
  x <- sub("^\\.", "_", x)
  x <- trimws(x)
  ifelse(nzchar(x), x, "unnamed")
}

#' Column names for a set of variants: one name each, and each one meaningful.
#'
#' Converting a VCF to PLINK without --set-missing-var-ids leaves every variant
#' named '.', so the IDs are neither unique nor informative. Those become
#' chr:bp, which is what plink2's own @:# template produces and is unique per
#' variant. Anything still duplicated after that gets a numbered suffix:
#' emitting one name twice would have jamovi silently merge the two into
#' 'rs1' and 'rs1 (2)', with no way to tell which variant is which.
variant_ids <- function(bim) {
  id <- trimws(bim$id)
  unnamed <- !nzchar(id) | id == "."
  id[unnamed] <- paste0(bim$chr[unnamed], ":",
                        format(bim$bp[unnamed], scientific = FALSE, trim = TRUE))
  cov_unique_names(id, character(0))
}
