
# The correctness layer: every genotype we decode is compared to what plink
# wrote for the same data. This is what catches a wrong index, a wrong byte
# offset, or a wrong allele orientation — each of which otherwise produces
# perfectly plausible genotypes.

test_that("decoded genotypes match plink --recode A on every cell", {
  skip_without_fixture("small")

  tier <- "small"
  bim <- read_bim_file(tier); fam <- read_fam_file(tier)
  n <- length(fam$iid); v <- length(bim$id)

  raw  <- slice_bed(tier, seq_len(v), n)
  dose <- bed_decode_slices(raw, n, v)

  o <- read_raw_oracle(tier)
  expect_identical(o$iid, fam$iid)      # .raw rows follow the .fam
  expect_identical(o$ids, bim$id)       # and columns follow the .bim

  # plink --recode A counts the *minor* allele, which is allele 1 for only
  # about half the variants in real data, so the comparison is per variant.
  is_a1 <- o$counted == bim$a1
  expect_true(all(is_a1 | o$counted == bim$a2))

  oriented <- dose
  oriented[, is_a1] <- 2L - dose[, is_a1]

  expect_identical(unname(is.na(oriented)), unname(is.na(o$values)))
  both <- !is.na(oriented) & !is.na(o$values)
  expect_equal(sum(oriented[both] != o$values[both]), 0)
})

test_that("a selected subset decodes identically to the same variants in full", {
  skip_without_fixture("small")

  tier <- "small"
  bim <- read_bim_file(tier); fam <- read_fam_file(tier)
  n <- length(fam$iid)

  set.seed(1)
  idx <- sort(sample.int(length(bim$id), 25))

  full <- bed_decode_slices(slice_bed(tier, seq_along(bim$id), n),
                            n, length(bim$id))
  sub  <- bed_decode_slices(slice_bed(tier, idx, n), n, length(idx))

  # Selecting must not change what a variant decodes to — this is the property
  # that makes seeking to offsets equivalent to reading the whole file.
  expect_identical(sub, full[, idx, drop = FALSE])
})

test_that("MAF and genotype counts agree with the oracle", {
  skip_without_fixture("small")

  tier <- "small"
  bim <- read_bim_file(tier); fam <- read_fam_file(tier)
  n <- length(fam$iid); v <- length(bim$id)
  dose <- bed_decode_slices(slice_bed(tier, seq_len(v), n), n, v)

  o <- read_raw_oracle(tier)
  is_a1 <- o$counted == bim$a1

  for (j in seq_len(min(v, 50))) {
    s <- snp_stats(dose[, j])
    col <- o$values[, j]

    expect_equal(s$n, sum(!is.na(col)))
    expect_equal(s$missing, sum(is.na(col)))

    # the oracle counts whichever allele plink named; ours counts allele 2
    cnt <- if (is_a1[j]) 2 - dose[, j] else dose[, j]
    expect_equal(sum(cnt == 0, na.rm = TRUE) + sum(cnt == 1, na.rm = TRUE) +
                 sum(cnt == 2, na.rm = TRUE), s$n)

    freq <- sum(col, na.rm = TRUE) / (2 * sum(!is.na(col)))
    expect_equal(s$maf, min(freq, 1 - freq), tolerance = 1e-9)
  }
})

test_that("HWE p-values are probabilities and behave at the extremes", {
  # a locus in perfect equilibrium is not rejected
  expect_gt(hwe_exact_p(25, 50, 25), 0.9)
  # complete absence of heterozygotes at equal frequencies is
  expect_lt(hwe_exact_p(50, 0, 50), 1e-10)
  # monomorphic is undefined, not zero
  expect_true(is.na(hwe_exact_p(100, 0, 0)))

  for (h in c(0, 10, 40, 50, 60))
    expect_true(hwe_exact_p(30, h, 30) >= 0 && hwe_exact_p(30, h, 30) <= 1)
})
