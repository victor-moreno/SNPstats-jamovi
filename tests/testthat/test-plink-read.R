
test_that("the lookup table encodes the PLINK 1.9 genotype codes", {
  # 00 = hom allele 1, 01 = MISSING, 10 = het, 11 = hom allele 2, low bits
  # first. 01 = missing is not the intuitive ordering and is the classic bug.
  expect_identical(PLINK_LUT[0x00 + 1L, ], c(0L, 0L, 0L, 0L))
  expect_identical(PLINK_LUT[0xff + 1L, ], c(2L, 2L, 2L, 2L))

  # one byte holding, in sample order: hom1, missing, het, hom2
  b <- 0L + 0L * 1L + 1L * 4L + 2L * 16L + 3L * 64L
  expect_identical(PLINK_LUT[b + 1L, ], c(0L, NA_integer_, 1L, 2L))
})

test_that("block decoding drops the padding in the last byte", {
  # 6 samples -> 2 bytes, the last holding 2 real genotypes and 2 pad slots
  block <- as.raw(c(0x00, 0x0f))
  g <- bed_decode_block(block, 6L)
  expect_length(g, 6L)
  expect_identical(g[1:4], c(0L, 0L, 0L, 0L))
})

test_that("a block of the wrong size is rejected", {
  expect_error(bed_decode_block(as.raw(c(0, 0)), 100L), "expected")
})

test_that("magic bytes are checked, and sample-major is refused by name", {
  expect_true(bed_check_magic(as.raw(c(0x6c, 0x1b, 0x01))))
  expect_error(bed_check_magic(as.raw(c(0x00, 0x00, 0x01))), "wrong magic")
  expect_error(bed_check_magic(as.raw(c(0x6c, 0x1b, 0x00))), "sample-major")
  expect_error(bed_check_magic(as.raw(c(0x6c, 0x1b, 0x09))), "layout byte")
  expect_error(bed_check_magic(as.raw(0x6c)), "fewer than 3")
})

test_that("slice decoding insists the byte count matches the variant count", {
  n <- 10L; bpv <- bed_bytes_per_variant(n)
  ok <- as.raw(rep(0, 3 * bpv))
  expect_equal(dim(bed_decode_slices(ok, n, 3L)), c(10L, 3L))

  # the structural check that the .bed slices and the .bim lines agree
  expect_error(bed_decode_slices(ok, n, 4L), "expected")
  expect_error(bed_decode_slices(ok[-1], n, 3L), "expected")
})

test_that("implausible dimensions are refused before allocating", {
  expect_error(bed_decode_slices(raw(0), 0L, 1L), "implausible sample")
  expect_error(bed_decode_slices(raw(0), 10L, 0L), "implausible variant")
  expect_error(bed_decode_slices(raw(0), 4e6, 1e5), "refusing to decode")
})

test_that(".bim and .fam parse, and .fam codes become usable factors", {
  bim <- read_bim(c("1\trs1\t0\t1000\tA\tG", "2\trs2\t0\t2000\tC\tT"))
  expect_identical(bim$id, c("rs1", "rs2"))
  expect_identical(bim$a1, c("A", "C"))
  expect_identical(bim$bp, c(1000, 2000))

  fam <- read_fam(c("F1 S1 0 0 1 2", "F2 S2 0 0 2 1", "F3 S3 0 0 0 -9"))
  expect_identical(fam$iid, c("S1", "S2", "S3"))
  expect_identical(as.character(fam$sex),   c("male", "female", NA))
  expect_identical(as.character(fam$pheno), c("case", "control", NA))
})

test_that("sex and phenotype codes outside the spec are missing, not guessed", {
  # A .fam and a .tfam are the same file, and a .ped's first six fields are a
  # .fam row, so this rule covers every format that carries either field.
  fam <- read_fam(c("F1 S1 0 0 M 2.5",     # letters, quantitative phenotype
                    "F2 S2 0 0 3 0",       # out of range, 0 = missing
                    "F3 S3 0 0 NA NA",
                    "F4 S4 0 0 -9 -9",
                    "F5 S5 0 0 1 1"))
  expect_identical(as.character(fam$sex),
                   c(NA, NA, NA, NA, "male"))
  expect_identical(as.character(fam$pheno),
                   c(NA, NA, NA, NA, "control"))
})

test_that("a genotype needing an allele the .bim does not name is missing", {
  # PLINK writes A1 = 0 where nothing carries the minor allele. A well-formed
  # file then has no sample with a copy of it; a malformed one does, and
  # "0/0" would read as a genotype rather than as the absence of one.
  g <- geno_factor(c(0L, 1L, 2L, NA), "0", "A")
  expect_identical(as.character(g), c(NA, NA, "A/A", NA))

  # the dosage itself is untouched — the number of copies of allele 2 is a
  # fact whatever the other allele turns out to be called
  expect_identical(snp_stats(c(0L, 1L, 2L, NA))$missing, 1L)
})

test_that("malformed text is refused rather than half-parsed", {
  expect_error(read_bim("1\trs1\t0"), "expected at least")
  expect_error(read_bim("1\trs1\t0\tnotanumber\tA\tG"), "non-numeric")
  expect_error(read_bim(character(0)), "empty")

  # a single enormous line is a denial of service against the engine
  expect_error(read_bim(strrep("x", 70000)), "refusing to parse")
})

test_that("genotype factors follow the .bim allele order and keep all levels", {
  g <- geno_factor(c(0L, 1L, 2L, NA), "A", "G")
  expect_identical(levels(g), c("A/A", "A/G", "G/G"))
  expect_identical(as.character(g), c("A/A", "A/G", "G/G", NA))

  # a monomorphic SNP still reports three levels, so counts never lose a class
  expect_identical(levels(geno_factor(c(0L, 0L), "C", "T")),
                   c("C/C", "C/T", "T/T"))

  # Labels are canonical, not source order: a1/a2 reversed must give the same
  # string for the same genotype, or two formats of one dataset would not
  # compare. dose counts allele 2, so reversing the roles flips the count.
  fwd <- geno_factor(c(0L, 1L, 2L), "A", "G")
  rev <- geno_factor(c(2L, 1L, 0L), "G", "A")
  expect_identical(as.character(fwd), as.character(rev))
  expect_identical(levels(fwd), levels(rev))
})

test_that("per-SNP statistics are counted the way the summary table reports", {
  s <- snp_stats(c(0L, 0L, 1L, 2L, NA))
  expect_equal(s$n, 4); expect_equal(s$missing, 1)
  expect_equal(c(s$n11, s$n12, s$n22), c(2, 1, 1))
  expect_equal(s$freq_a2, 3 / 8)
  expect_equal(s$maf, 3 / 8)          # a2 is the minor allele here

  # MAF is the *minor* allele frequency whichever allele that is
  expect_equal(snp_stats(c(2L, 2L, 2L, 1L))$maf, 1 / 8)
  expect_true(is.na(snp_stats(rep(NA_integer_, 4))$maf))
})

test_that("column names are made safe without becoming ambiguous", {
  expect_identical(sanitise_id(c("rs1", "a`b", ".hidden", " pad ")),
                   c("rs1", "a_b", "_hidden", "pad"))
  expect_identical(sanitise_id(""), "unnamed")
})

test_that("vectorised statistics agree with the per-SNP version", {
  # snp_stats_all() replaced a per-column loop for speed. The two must not
  # drift apart, so every quantity is compared column by column.
  set.seed(42)
  dose <- matrix(sample(c(0L, 1L, 2L, NA), 400 * 25, TRUE, c(.4, .35, .2, .05)),
                 nrow = 400, ncol = 25)

  all_st <- snp_stats_all(dose)
  for (j in seq_len(ncol(dose))) {
    one <- snp_stats(dose[, j])
    expect_equal(all_st$n[j],       one$n)
    expect_equal(all_st$missing[j], one$missing)
    expect_equal(all_st$n11[j],     one$n11)
    expect_equal(all_st$n12[j],     one$n12)
    expect_equal(all_st$n22[j],     one$n22)
    expect_equal(all_st$maf[j],     one$maf)
    expect_equal(all_st$hwe_p[j],   one$hwe_p)
  }
})

test_that("a monomorphic or all-missing column does not break the vectorised path", {
  dose <- cbind(rep(0L, 50), rep(NA_integer_, 50), rep(2L, 50))
  st <- snp_stats_all(dose)
  expect_equal(st$n, c(50, 0, 50))
  expect_equal(st$maf, c(0, NA, 0))
  expect_true(all(is.na(st$hwe_p)))       # undefined without both alleles
})

test_that("the trio is checked against the .bed's own size, not against itself", {
  # 200 samples is 50 bytes a variant, so a 40-variant .bed is 3 + 40 * 50.
  expect_true(bed_check_dims("2003,40,200", 200))

  # The failure the slice length check cannot see: a .fam from a different build
  # of the dataset. Both sides then compute ceil(196/4) = 49 and agree with each
  # other about a payload that decodes the wrong variants.
  expect_error(bed_check_dims("2003,40,196", 196), "not from the same dataset")

  # ...and the .fam that was measured has to be the one that arrived
  expect_error(bed_check_dims("2003,40,200", 199), "sample file changed")

  # a stale .bim: fewer lines than the .bed has blocks
  expect_error(bed_check_dims("2003,30,200", 200), "not from the same dataset")

  # absent means an import written before the check existed, which the caller
  # reports rather than refuses
  expect_false(bed_check_dims("", 200))
  expect_error(bed_check_dims("2003,40", 200), "unreadable")
  expect_error(bed_check_dims("a,b,c", 200), "unreadable")
})

test_that("per-sample statistics count missingness and heterozygosity", {
  #            snp1 snp2 snp3 snp4
  dose <- rbind(c(0L,  1L,  1L,  2L),    # no missing, 2 of 4 het
                c(NA,  NA,  0L,  0L),    # half missing, no het
                c(1L,  1L,  1L,  1L))    # all het
  ss <- sample_stats(dose)
  expect_equal(ss$missing_rate, c(0, 0.5, 0))
  expect_equal(ss$het_rate, c(0.5, 0, 1))
  expect_equal(ss$n_obs, c(4, 2, 4))

  # an all-missing sample has no defined heterozygosity rather than zero
  expect_true(is.na(sample_stats(matrix(NA_integer_, 1, 4))$het_rate))
})
