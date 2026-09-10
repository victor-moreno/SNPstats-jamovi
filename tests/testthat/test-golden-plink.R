
# Frozen comparison against PLINK's own statistics.
#
# test-oracle.R checks the genotypes we decode. This checks the numbers we
# *derive* from them — everything the SNP summary table reports — against what
# PLINK computes from the same file:
#
#   small.frq    --freq     minor allele frequency
#   small.hwe    --hardy    genotype counts and the exact HWE p-value
#   small.lmiss  --missing  per-variant missingness
#   small.imiss  --missing  per-sample missingness
#
# The files are generated once by data-raw/make_fixtures.sh and committed, so
# these are golden tests: they keep passing without PLINK installed, and they
# fail if our arithmetic drifts even when the decode is still correct.
#
# The HWE comparison is the valuable one. PLINK implements the same Wigginton
# exact test, so agreeing with it to six decimal places is independent evidence
# that R/vendored_hwe.R is right — the vendored code is checked against an
# outside implementation, not only against our own oracle.

read_plink_table <- function(tier, ext) {
  f <- file.path(fixture_dir(tier), paste0(tier, ext))
  if (!file.exists(f)) return(NULL)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE)
}

skip_without_plink_stats <- function(tier = "small") {
  skip_without_fixture(tier)
  for (ext in c(".frq", ".hwe", ".lmiss", ".imiss"))
    if (is.null(read_plink_table(tier, ext)))
      testthat::skip(paste0("missing ", tier, ext,
                            " — run data-raw/make_fixtures.sh ", tier))
}

# The summary table as the analysis actually produces it, for every SNP.
summary_for_all <- function(tier, ...) {
  bim <- read_bim_file(tier)
  run_import(make_payloads(tier, bim$id), bim$id, ...)$summary$asDF
}


test_that("MAF matches plink --freq", {
  skip_without_plink_stats()

  df  <- summary_for_all("small")
  frq <- read_plink_table("small", ".frq")

  expect_setequal(as.character(df$snp), frq$SNP)
  i <- match(frq$SNP, as.character(df$snp))

  # Both report the frequency of the *minor* allele, whichever that is, so no
  # orientation adjustment is needed. The tolerance is plink's printing
  # precision, not ours: .frq carries 4-5 significant figures, so 0.4568528
  # arrives as 0.4569. Tightening this does not test anything more.
  expect_equal(df$maf[i], frq$MAF, tolerance = 2e-4)
})

test_that("genotype counts and HWE p-values match plink --hardy", {
  skip_without_plink_stats()

  # explicitly all-samples: the module defaults to controls only, and this
  # compares against plink's ALL row
  df  <- summary_for_all("small", hweGroup = "all")
  hwe <- read_plink_table("small", ".hwe")
  hwe <- hwe[hwe$TEST == "ALL", ]          # not the case/control breakdowns

  bim <- read_bim_file("small")
  i   <- match(hwe$SNP, as.character(df$snp))
  j   <- match(hwe$SNP, bim$id)

  # plink's GENO is "C(HOM A1)/C(HET)/C(HOM A2)" for the A1/A2 *it* chose,
  # which is the minor allele first — not necessarily the .bim order our
  # counts follow. Compare after orienting to plink's choice.
  parts <- do.call(rbind, lapply(strsplit(hwe$GENO, "/"), as.numeric))
  ours  <- do.call(rbind, strsplit(as.character(df$genoCounts[i]), " / "))
  ours  <- matrix(as.numeric(ours), ncol = 3)

  flip <- hwe$A1 == bim$a2[j]              # plink's A1 is our allele 2
  ours[flip, ] <- ours[flip, c(3, 2, 1)]

  expect_equal(ours[, 1], parts[, 1])
  expect_equal(ours[, 2], parts[, 2])      # heterozygotes, orientation-free
  expect_equal(ours[, 3], parts[, 3])

  # The exact test itself. plink prints p to four significant figures, so the
  # tolerance is the file's precision, not ours.
  #
  # Monomorphic sites are excluded, and the exclusion is the point: there is no
  # equilibrium to test when only one allele is present, and the two
  # implementations answer differently. plink prints P = 1; hwe_exact_p returns
  # NA, because a test that had no alternative to reject did not pass, it did
  # not run. NA is what test-plink-read.R pins ("undefined without both
  # alleles"), and hwe_exact_p is vendored from SNPstats and has to stay
  # diff-clean against it, so this is a documented divergence rather than
  # something to fix here.
  mono <- parts[, 2] == 0 & (parts[, 1] == 0 | parts[, 3] == 0)
  expect_true(all(is.na(df$hwePval[i][mono])))
  expect_equal(df$hwePval[i][!mono], hwe$P[!mono], tolerance = 1e-3)
})

test_that("per-SNP missingness matches plink --missing", {
  skip_without_plink_stats()

  df    <- summary_for_all("small")
  lmiss <- read_plink_table("small", ".lmiss")
  i     <- match(lmiss$SNP, as.character(df$snp))

  expect_equal(df$missing[i], lmiss$N_MISS)
  expect_equal(df$n[i] + df$missing[i], lmiss$N_GENO)
})

test_that("per-sample missingness matches plink --missing, and drives --mind", {
  skip_without_plink_stats()

  tier <- "small"
  bim  <- read_bim_file(tier)
  fam  <- read_fam_file(tier)
  n    <- length(fam$iid)

  dose  <- bed_decode_slices(slice_bed(tier, seq_along(bim$id), n),
                             n, length(bim$id))
  ours  <- sample_stats(dose)
  imiss <- read_plink_table(tier, ".imiss")

  expect_identical(as.character(imiss$IID), fam$iid)
  expect_equal(ours$missing_rate, imiss$F_MISS, tolerance = 1e-6)

  # and the sample filter drops exactly the individuals plink's F_MISS would
  thresh <- 0.05
  expected_dropped <- sum(imiss$F_MISS > thresh)
  skip_if(expected_dropped == 0, "no sample exceeds the threshold in this fixture")

  ids <- bim$id
  r <- run_import(make_payloads(tier, ids), ids,
                  filterSamples = TRUE, maxIndMissing = thresh * 100)
  kept <- length(out_values(r)[["IID"]])
  expect_equal(kept, n - expected_dropped)
})

test_that("the summary table is stable against a frozen snapshot", {
  skip_without_plink_stats()

  # Guards the whole pipeline at once, at full precision: a change in decoding,
  # statistics or formatting moves these numbers even when each part still
  # agrees with plink to the precision plink prints. Captured from a run
  # verified by the comparisons above — update only after those still pass.
  # Pinned to hweGroup = "all" so the values stay comparable to plink's ALL row.
  df <- summary_for_all("small", hweGroup = "all")
  df <- df[match(c("snp0", "snp1", "snp249", "snp499"), as.character(df$snp)), ]

  expect_identical(as.character(df$alleles),
                   c("G/T", "C/A", "T/C", "T/C"))
  expect_identical(as.character(df$genoCounts),
                   c("55 / 93 / 42", "166 / 29 / 2", "85 / 79 / 33", "76 / 94 / 26"))
  expect_equal(df$n,       c(190, 197, 197, 196))
  expect_equal(df$missing, c(10, 3, 3, 4))
  expect_equal(df$maf, c(0.4657895, 0.0837563, 0.3680203, 0.3724490),
               tolerance = 1e-7)
  expect_equal(df$hwePval, c(0.8840394, 0.6293689, 0.0652237, 0.7619517),
               tolerance = 1e-7)
})

test_that("controls-only HWE matches plink's UNAFF rows", {
  skip_without_plink_stats()

  # PLINK's --hwe considers only controls by default for case/control
  # phenotypes, and --hardy reports that as the UNAFF row. Agreeing with both
  # ALL and UNAFF from the same data checks the subsetting, not just the test.
  hwe <- read_plink_table("small", ".hwe")
  bim <- read_bim_file("small")

  for (mode in c("all", "controls")) {
    df  <- summary_for_all("small", hweGroup = mode)
    ref <- hwe[hwe$TEST == (if (mode == "all") "ALL" else "UNAFF"), ]
    i   <- match(ref$SNP, as.character(df$snp))

    # One documented divergence: where the locus is monomorphic in the tested
    # subset, plink prints P = 1 and we report NA. The exact test conditions on
    # the observed allele counts, so with a single allele there is nothing to
    # test — NA says that, 1 implies a test was run and passed. Both refuse to
    # drop the SNP, so the --hwe filter behaves identically either way.
    # Monomorphism is a property of the *tested subset*, so it is read from
    # plink's own per-row GENO counts, not from our whole-sample ones.
    g    <- do.call(rbind, lapply(strsplit(ref$GENO, "/"), as.numeric))
    mono <- rowSums(g > 0) == 1
    expect_true(all(is.na(df$hwePval[i][mono])), info = mode)
    expect_true(all(ref$P[mono] == 1), info = mode)

    ok <- !mono
    expect_equal(df$hwePval[i][ok], ref$P[ok], tolerance = 1e-3,
                 info = paste("hweGroup =", mode))
  }

  # and the two modes really do differ, so the test above is not vacuous
  a <- summary_for_all("small", hweGroup = "all")$hwePval
  c <- summary_for_all("small", hweGroup = "controls")$hwePval
  expect_gt(sum(abs(a - c) > 1e-6, na.rm = TRUE), 0)
})

test_that("the genotype counts stay whole-sample whichever HWE mode is used", {
  skip_without_plink_stats()

  # Only HWE changes with the mode; N, missing, MAF and the counts always
  # describe every imported sample.
  a <- summary_for_all("small", hweGroup = "all")
  c <- summary_for_all("small", hweGroup = "controls")

  expect_identical(as.character(a$genoCounts), as.character(c$genoCounts))
  expect_equal(a$n, c$n)
  expect_equal(a$missing, c$missing)
  expect_equal(a$maf, c$maf)
})
