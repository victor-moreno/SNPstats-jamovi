
# .tped and VCF, and the cross-format equivalence that is their best test.
#
# The fixture holds the same data in every format, so "the same SNPs read from
# .bed, .tped and .vcf must give the same genotypes" needs no oracle at all —
# and it is stronger than checking each against `--recode A` separately,
# because a format-specific misread would have to be made identically by two
# independent parsers to slip through.
#
# What must NOT be compared is the dosage. Which allele is called 1 and which 2
# differs legitimately by format: a VCF names REF first, a .bim picks by
# frequency. The genotypes are the same either way, so the comparison is on the
# emitted genotype strings.

# The genotypes as the module actually emits them. Written through
# geno_factor() rather than reconstructed by hand: a hand-rolled version has to
# re-derive the label/dose relationship, and getting that subtly wrong makes a
# correct parser look broken — which it did once while writing these.
geno_strings <- function(dose, bim) {
  vapply(seq_len(ncol(dose)), function(j)
    paste(as.character(geno_factor(dose[, j], bim$a1[j], bim$a2[j])),
          collapse = "|"), "")
}

bed_all <- function(tier) {
  bim <- read_bim_file(tier); fam <- read_fam_file(tier)
  list(bim = bim, fam = fam,
       dose = bed_decode_slices(slice_bed(tier, seq_along(bim$id),
                                          length(fam$iid)),
                                length(fam$iid), length(bim$id)))
}


test_that("allele roles are derived the way plink assigns them", {
  # no .bim in a .tped, so A1 = the minor allele, A2 = the major one
  expect_identical(tped_alleles(c("A", "A", "A", "G")), c("G", "A"))
  expect_identical(tped_alleles(c("C", "C")), c("C", "C"))   # monomorphic
  expect_identical(tped_alleles(c("0", "0")), c("0", "0"))   # nothing called
  expect_identical(tped_alleles(c("A", "0", "G", "G")), c("A", "G"))

  # Every missing marker in circulation, not just plink's own 0
  expect_identical(tped_alleles(c("A", ".", "-", "N", "G", "G")), c("A", "G"))

  # A third code no longer takes the import down with it: the two commonest
  # win and the rest are read as no-calls by .call_dose
  expect_identical(tped_alleles(c("A", "A", "G", "G", "C")), c("G", "A"))
})

test_that("a call in a code the variant does not have is missing, not a genotype", {
  # The dangerous shape: comparing against A2 alone counts an unrecognised
  # code as "not A2" and reports it as homozygous A1 — an invented genotype.
  cd <- SNPstats:::.call_dose(c("A", "A", "G", "C", "0", "A"),
                   c("A", "G", "G", "C", "0", "0"), "A", "G")
  expect_identical(cd$dose, c(0L, 1L, 2L, NA, NA, NA))
  expect_identical(cd$unreadable, 1L)      # the C/C, not the 0s

  # a half call is no call
  expect_identical(SNPstats:::.call_dose("A", "0", "A", "G")$dose, NA_integer_)
})

test_that("a stray allele code is imported as missing, not refused", {
  tfam <- c("f1 s1 0 0 1 1", "f1 s2 0 0 2 1", "f1 s3 0 0 1 2")
  tped <- c("1 snpA 0 100 A A A G X Y",
            "1 snpB 0 200 G G A G A A")
  t <- read_tped(tped, tfam)
  # A is the major allele of snpA, so it is A2 and the dosage counts it
  expect_identical(t$dose[, 1], c(2L, 1L, NA))
  expect_identical(t$unreadable, 1L)
  expect_identical(unname(t$bim$a1), c("G", "G"))
})

test_that("GT strings become ALT counts, phased or not", {
  expect_identical(SNPstats:::.vcf_gt_dose(
    c("0/0", "0/1", "1/1", "0|1", "1|0", "./.", ".")),
    c(0L, 1L, 2L, 1L, 1L, NA, NA))
  # haploid calls appear on the sex chromosomes
  expect_identical(SNPstats:::.vcf_gt_dose(c("0", "1")), c(0L, 1L))

  # The line is biallelic by the time this runs, so an allele index of 2 is
  # not one of its alleles. It used to count as a REF copy.
  expect_identical(SNPstats:::.vcf_gt_dose(c("0/2", "2/2", "2")),
                   c(NA_integer_, NA_integer_, NA_integer_))
  # and anything that is not a GT at all
  expect_identical(SNPstats:::.vcf_gt_dose(c("", "A/A", "0x1", NA)),
                   rep(NA_integer_, 4))
})

test_that(".tped decodes to the same genotypes as the .bed", {
  skip_without_fixture("small")

  b <- bed_all("small")
  t <- read_tped(readLines(fx("small", ".tped")), readLines(fx("small", ".tfam")))

  expect_identical(t$fam$iid, b$fam$iid)          # sample order preserved
  expect_identical(t$bim$id, b$bim$id)            # variant order preserved
  expect_identical(dim(t$dose), dim(b$dose))
  expect_identical(is.na(t$dose), is.na(b$dose))  # missingness identical

  expect_identical(geno_strings(t$dose, t$bim), geno_strings(b$dose, b$bim))
})

test_that("VCF decodes to the same genotypes as the .bed", {
  skip_without_fixture("small")

  b <- bed_all("small")
  v <- read_vcf(readLines(fx("small", ".vcf")))

  expect_identical(v$fam$iid, b$fam$iid)
  expect_equal(v$skipped, 0)                      # the fixture is all SNVs
  expect_identical(is.na(v$dose), is.na(b$dose))

  i <- match(b$bim$id, v$bim$id)
  expect_false(anyNA(i))
  vb <- list(a1 = v$bim$a1[i], a2 = v$bim$a2[i])
  expect_identical(geno_strings(v$dose[, i, drop = FALSE], vb),
                   geno_strings(b$dose, b$bim))
})

test_that("the allele labels legitimately differ even though genotypes agree", {
  skip_without_fixture("small")

  # Worth pinning: a future change that made the labels match would probably
  # have flipped a dosage somewhere. VCF names REF first, so its A1/A2 are the
  # other way round from a .bim that picks by frequency.
  b <- bed_all("small")
  v <- read_vcf(readLines(fx("small", ".vcf")))
  i <- match(b$bim$id, v$bim$id)
  expect_true(sum(v$bim$a1[i] == b$bim$a1) < length(i))
})

test_that("a VCF supplies sample IDs but no sex or phenotype", {
  skip_without_fixture("small")

  v <- read_vcf(readLines(fx("small", ".vcf")))
  expect_true(all(is.na(v$fam$sex)))
  expect_true(all(is.na(v$fam$pheno)))
  expect_true(v$fam$iid_unique)
})

test_that("multi-allelic sites and indels are skipped and counted", {
  hdr <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"
  v <- read_vcf(c(
    "##fileformat=VCFv4.2", hdr,
    "1\t100\tsnv\tA\tG\t.\t.\t.\tGT\t0/0\t1/1",
    "1\t200\tmulti\tA\tG,T\t.\t.\t.\tGT\t0/0\t0/1",     # multi-allelic
    "1\t300\tindel\tAT\tA\t.\t.\t.\tGT\t0/0\t0/1"))     # indel
  expect_identical(v$bim$id, "snv")
  expect_equal(v$skipped, 2)
  expect_identical(as.integer(v$dose[, 1]), c(0L, 2L))
})

test_that("GT is located from FORMAT, not assumed to be first", {
  hdr <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"
  v <- read_vcf(c(hdr,
    "1\t100\trs1\tA\tG\t.\t.\t.\tDP:GT:GQ\t30:0/1:99\t25:1/1:80"))
  expect_identical(as.integer(v$dose[, 1]), c(1L, 2L))

  expect_error(read_vcf(c(hdr,
    "1\t100\trs1\tA\tG\t.\t.\t.\tDP:GQ\t30:99\t25:80")), "no GT")
})

test_that("malformed text is refused rather than half-read", {
  expect_error(read_tped("1 rs1 0 100 A A", c("F1 S1 0 0 1 1", "F2 S2 0 0 1 1")),
               "expected 8")
  expect_error(read_tped(character(0), "F1 S1 0 0 1 1"), "no .tped data")

  expect_error(read_vcf("1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/0"), "no #CHROM")
  expect_error(read_vcf(c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")),
               "no samples")
})

test_that("a variant that does not vary imports instead of aborting the run", {
  # A .tped or .ped names its alleles from what it observes, so a monomorphic
  # SNP names the same allele twice and all three genotype labels collapse to
  # one string. factor() rejects duplicated levels, which took the whole import
  # down over a SNP that was merely uninformative — and a small sample makes
  # monomorphic SNPs common.
  expect_identical(levels(geno_factor(c(0L, 0L), "C", "C")), "C/C")
  expect_identical(as.character(geno_factor(c(2L, 2L, NA), "C", "C")),
                   c("C/C", "C/C", NA))

  tfam <- c("F1 S1 0 0 1 1", "F2 S2 0 0 2 2", "F3 S3 0 0 1 1")
  t <- read_tped(c("1 rsMono 0 100 C C C C C C",
                   "1 rsPoly 0 200 A A A G G G"), tfam)
  expect_identical(t$bim$a1, c("C", "G"))
  expect_identical(t$bim$a2, c("C", "A"))
  expect_identical(as.character(geno_factor(t$dose[, 1], t$bim$a1[1], t$bim$a2[1])),
                   rep("C/C", 3))

  # and end to end, where the failure actually showed
  v <- out_values(run_import(
    list(geno = as_payload(c("1 rsMono 0 100 C C C C C C",
                             "1 rsPoly 0 200 A A A G G G")),
         bim = as_payload(character(0)), fam = as_payload(tfam)),
    c("rsMono", "rsPoly"), sourceFormat = "tped", showSummary = FALSE))
  expect_identical(as.character(v$rsMono), rep("C/C", 3))
  expect_identical(as.character(v$rsPoly), c("A/A", "A/G", "G/G"))
})

test_that("a .tped import produces the same columns as the .bed import", {
  skip_without_fixture("small")

  # end to end through the analysis, not just the parser
  bim <- read_bim_file("small")
  ids <- bim$id[1:8]
  tped <- readLines(fx("small", ".tped"))
  keep <- tped[match(ids, vapply(strsplit(trimws(tped), "[ \t]+"),
                                 `[`, "", 2L))]

  v_bed <- out_values(run_import(make_payloads("small", ids), ids))
  v_tpd <- out_values(run_import(
    list(geno = as_payload(keep),
         bim  = as_payload(character(0)),
         fam  = as_payload(readLines(fx("small", ".tfam"))),
         n    = length(read_fam_file("small")$iid)),
    ids, sourceFormat = "tped"))

  for (id in ids)
    expect_identical(as.character(v_tpd[[id]]), as.character(v_bed[[id]]),
                     info = id)
  expect_identical(as.character(v_tpd$IID), as.character(v_bed$IID))
})

test_that("a VCF import reports its format and needs no companion files", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:5]
  vcf <- readLines(fx("small", ".vcf"))
  hdr <- vcf[startsWith(vcf, "#CHROM")]
  body <- vcf[!startsWith(vcf, "#")]
  keep <- body[match(ids, vapply(strsplit(body, "\t"), `[`, "", 3L))]

  r <- run_import(list(geno = as_payload(c(hdr, keep)),
                       bim = as_payload(character(0)),
                       fam = as_payload(character(0)), n = 200),
                  ids, sourceFormat = "vcf")

  expect_equal(nrow(r$summary$asDF), 5)
  pv <- r$provenance$asDF
  expect_match(pv$value[pv$item == "Format"], "VCF")
})


# ── .ped / .map ──────────────────────────────────────────────────────────────

test_that(".map is read with or without the centimorgan column", {
  four  <- read_map(c("1\trs1\t0\t100", "1\trs2\t0\t200"))
  three <- read_map(c("1\trs1\t100", "1\trs2\t200"))
  expect_equal(four$bp, c(100, 200))
  expect_equal(three$bp, c(100, 200))          # bp moves up when cM is absent
  expect_identical(four$id, three$id)

  expect_error(read_map("1\trs1"), "fewer than 3")
  expect_error(read_map(c("1\trs1\t0\tx")), "non-numeric")
})

test_that(".ped decodes to the same genotypes as the .bed", {
  skip_without_fixture("small")

  b <- bed_all("small")
  p <- read_ped(readLines(fx("small", ".ped")), readLines(fx("small", ".map")))

  expect_identical(p$bim$id, b$bim$id)
  expect_identical(is.na(p$dose), is.na(b$dose))
  expect_identical(geno_strings(p$dose, p$bim), geno_strings(b$dose, b$bim))
})

test_that("a .ped supplies its own samples, sex and phenotype", {
  skip_without_fixture("small")

  # the first six fields of a .ped line are exactly a .fam row, which is why
  # no companion sample file is needed
  f <- read_fam_file("small")
  p <- read_ped(readLines(fx("small", ".ped")), readLines(fx("small", ".map")))

  expect_identical(p$fam$iid, f$iid)
  expect_identical(as.character(p$fam$sex), as.character(f$sex))
  expect_identical(as.character(p$fam$pheno), as.character(f$pheno))
})

test_that("a .ped is checked against its .map rather than half-read", {
  map <- c("1\trs1\t0\t100", "1\trs2\t0\t200")
  expect_silent(read_ped("F1 S1 0 0 1 1 A A G G", map))   # 6 + 2*2 fields
  expect_error(read_ped("F1 S1 0 0 1 1 A A G", map), "expected 10")
  expect_error(read_ped(character(0), map), "no .ped data")
})

test_that("a reduced .ped imports end to end and matches the .bed import", {
  skip_without_fixture("small")

  # This is what the browser sends: the six sample fields plus only the
  # selected variants' allele pairs, in .map order.
  bim <- read_bim_file("small")
  ids <- bim$id[c(2, 7, 11)]
  mapl <- readLines(fx("small", ".map"))
  mid  <- vapply(strsplit(trimws(mapl), "[ \t]+"), `[`, "", 2L)
  j    <- match(ids, mid)

  cols <- as.vector(rbind(6 + 2 * j - 1, 6 + 2 * j))
  ped  <- vapply(strsplit(trimws(readLines(fx("small", ".ped"))), "[ \t]+"),
                 function(f) paste(c(f[1:6], f[cols]), collapse = " "), "")

  v_ped <- out_values(run_import(
    list(geno = as_payload(ped), bim = as_payload(mapl[j]),
         fam = as_payload(character(0)), n = 200),
    ids, sourceFormat = "ped"))
  v_bed <- out_values(run_import(make_payloads("small", ids), ids))

  for (id in ids)
    expect_identical(as.character(v_ped[[id]]), as.character(v_bed[[id]]), info = id)
  expect_identical(as.character(v_ped$IID), as.character(v_bed$IID))
})
