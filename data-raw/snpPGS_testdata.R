# snpPGS_testdata.R
# ─────────────────────────────────────────────────────────────────────────────
# Generates a synthetic test dataset and accompanying weights files that
# exercise every distinct code path in the SNPstats PGS module:
#
#   SNP FORMAT VARIANTS
#     slash_clean      A/G, C/T  — clean slash-separated, no issues
#     slash_nullallele A/G       — contains 0/0 null-allele codings
#     pipe_sep         A|G       — pipe separator (normalised to slash)
#     concat_sep       AG        — concatenated two-char (e.g. "AG", "GG")
#     strand_flip      C/A       — weights file lists G/T; complement matches
#     ambiguous_at     A/T       — ambiguous A/T SNP (kept with warning)
#     ambiguous_cg     C/G       — ambiguous C/G SNP (kept with warning)
#     numeric_dos      0/1/2     — numeric dosage column (no allele QC)
#     multiallelic     A/G/C     — three alleles; should be excluded
#     all_missing      NA        — all values NA; should be excluded
#     monomorphic      A/A only  — single genotype; should be excluded
#     mismatch         A/G       — weights file lists C/T; mismatch, excluded
#     no_allele_info            — weight file row has blank effect/other allele
#
#   RESPONSE VARIANTS (one column each)
#     resp_binary      factor 0/1                — binary, two levels
#     resp_polytomous  factor 1/2/3              — polytomous, three levels
#     resp_continuous  numeric                   — continuous linear
#     resp_smalln      binary, last 190 rows NA  — triggers near-separation
#
#   COVARIATES
#     cov_factor       factor (Male/Female)      — two-level factor
#     cov_numeric      numeric                   — continuous covariate
#     cov_multicategory factor (A/B/C)           — three-level factor
#
#   WEIGHTS FILES (written to disk as CSV)
#     weights_standard.csv  — normal PGS Catalog format with all SNPs
#     weights_partial.csv   — only a subset of SNPs (tests coverage < 100%)
#     weights_negweights.csv— includes negative weights (tests max_possible)
#     weights_noalleles.csv — effect_allele / other_allele columns blank
#     weights_header.csv    — PGS Catalog-style # header metadata rows
#
# Usage:
#   source("snpPGS_testdata.R")
#   # Datasets are returned by make_pgs_testdata(); weights files are written
#   # to the directory given by `weights_dir`.
#
# ─────────────────────────────────────────────────────────────────────────────

make_pgs_testdata <- function(
    n            = 300,        # number of individuals
    seed         = 42,
    weights_dir  = tempdir()   # where to write weights CSV files
) {

  set.seed(seed)

  # ── helpers ────────────────────────────────────────────────────────────────

  # Random biallelic genotypes as "A/B" strings, with specified frequencies
  rsnp <- function(n, a1, a2, p1 = 0.5, na_frac = 0) {
    g <- sample(
      c(paste0(a1,"/",a1), paste0(a1,"/",a2), paste0(a2,"/",a2)),
      n, replace = TRUE,
      prob = c(p1^2, 2*p1*(1-p1), (1-p1)^2)
    )
    if (na_frac > 0)
      g[sample.int(n, size = round(n * na_frac))] <- NA_character_
    g
  }

  # ── SNP columns ────────────────────────────────────────────────────────────

  # 1. Clean slash-separated A/G — perfect match to weights file
  slash_clean      <- rsnp(n, "A", "G", p1 = 0.4)

  # 2. Slash-separated with ~8% null-allele codings (0/0 → NA after cleaning)
  slash_nullallele <- rsnp(n, "C", "T", p1 = 0.35, na_frac = 0)
  null_idx <- sample.int(n, size = round(n * 0.08))
  slash_nullallele[null_idx] <- "0/0"

  # 3. Pipe separator (A|G) — normalised by the cleaning step
  pipe_sep <- rsnp(n, "A", "C", p1 = 0.45)
  pipe_sep[!is.na(pipe_sep)] <- gsub("/", "|", pipe_sep[!is.na(pipe_sep)], fixed = TRUE)

  # 4. Concatenated two-char ("AG", "GG") — no separator
  concat_snp <- rsnp(n, "A", "G", p1 = 0.6)
  concat_sep <- gsub("/", "", concat_snp, fixed = TRUE)

  # 5. Strand-flip: data has C/A, weights file will say G/T  (complement match)
  strand_flip <- rsnp(n, "C", "A", p1 = 0.5)

  # 6. Ambiguous A/T SNP — kept with warning flag
  ambiguous_at <- rsnp(n, "A", "T", p1 = 0.45)

  # 7. Ambiguous C/G SNP — kept with warning flag
  ambiguous_cg <- rsnp(n, "C", "G", p1 = 0.55)

  # 8. Numeric dosage 0/1/2 — no allele QC possible
  numeric_dos <- sample(c(0L, 1L, 2L), n, replace = TRUE, prob = c(0.25, 0.50, 0.25))

  # 9. Multiallelic — should be excluded (three distinct alleles)
  multiallelic <- sample(c("A/G", "A/C", "G/C", NA), n, replace = TRUE,
                         prob = c(0.3, 0.3, 0.3, 0.1))

  # 10. All-missing — should be excluded
  all_missing <- rep(NA_character_, n)

  # 11. Monomorphic — all homozygous reference, should be excluded
  monomorphic <- rep("A/A", n)
  monomorphic[sample.int(n, 5)] <- NA_character_   # a few NAs, still monomorphic

  # 12. Allele mismatch — data has A/G but weights file says C/T; excluded
  mismatch <- rsnp(n, "A", "G", p1 = 0.5)

  # 13. Missing > threshold (15% NA) — survives allele QC but caught by miss filter
  high_missing <- rsnp(n, "A", "G", p1 = 0.4, na_frac = 0.18)

  # 14. SNP with ~30% missing — tests mean imputation / SNP-wise strategies
  partial_missing <- rsnp(n, "C", "T", p1 = 0.5, na_frac = 0.30)

  # 15. SNP used only in partial weights file (tests coverage < 100%)
  partial_only    <- rsnp(n, "A", "G", p1 = 0.38)

  # ── Response columns ───────────────────────────────────────────────────────

  # Binary (0/1); correlated weakly with slash_clean dosage
  dos_clean <- (slash_clean == "A/A") * 2 +
               (slash_clean %in% c("A/G","G/A")) * 1
  lp_binary <- -0.5 + 0.3 * ifelse(is.na(dos_clean), 0, dos_clean)
  resp_binary <- factor(rbinom(n, 1, plogis(lp_binary)), levels = c(0, 1))

  # Polytomous (three ordered levels)
  lp2 <- lp_binary + rnorm(n, 0, 0.5)
  resp_polytomous <- factor(
    ifelse(lp2 < -0.3, 1L, ifelse(lp2 < 0.3, 2L, 3L)),
    levels = c(1, 2, 3)
  )

  # Continuous
  resp_continuous <- 10 + 2 * ifelse(is.na(dos_clean), 0, dos_clean) + rnorm(n, 0, 3)

  # Small-n binary: only 10 cases — likely to trigger near-separation / wild ORs
  resp_smalln <- resp_binary
  resp_smalln[sample(which(resp_binary == 1), size = sum(resp_binary == 1) - 10)] <- NA

  # ── Covariate columns ──────────────────────────────────────────────────────

  cov_factor        <- factor(sample(c("Male", "Female"), n, replace = TRUE))
  cov_numeric       <- rnorm(n, mean = 50, sd = 10)
  cov_multicategory <- factor(sample(c("A", "B", "C"), n, replace = TRUE,
                                     prob = c(0.4, 0.35, 0.25)))

  # ── Assemble data frame ────────────────────────────────────────────────────

  dat <- data.frame(
    # SNP columns
    snp_slash_clean      = slash_clean,
    snp_slash_null       = slash_nullallele,
    snp_pipe             = pipe_sep,
    snp_concat           = concat_sep,
    snp_strand_flip      = strand_flip,
    snp_ambiguous_at     = ambiguous_at,
    snp_ambiguous_cg     = ambiguous_cg,
    snp_numeric_dos      = numeric_dos,
    snp_multiallelic     = multiallelic,
    snp_all_missing      = all_missing,
    snp_monomorphic      = monomorphic,
    snp_mismatch         = mismatch,
    snp_high_missing     = high_missing,
    snp_partial_missing  = partial_missing,
    snp_partial_only     = partial_only,

    # Response columns
    resp_binary          = resp_binary,
    resp_polytomous      = resp_polytomous,
    resp_continuous      = resp_continuous,
    resp_smalln          = resp_smalln,

    # Covariates
    cov_factor           = cov_factor,
    cov_numeric          = cov_numeric,
    cov_multicategory    = cov_multicategory,

    stringsAsFactors     = FALSE
  )

  # ── Weights files ──────────────────────────────────────────────────────────

  # Helper: write a data.frame as CSV with optional # header lines
  write_weights <- function(df, filename, header_lines = NULL) {
    path <- file.path(weights_dir, filename)
    if (!is.null(header_lines)) {
      writeLines(header_lines, path)
      write.table(df, path, sep = ",", row.names = FALSE,
                  col.names = TRUE, quote = FALSE, append = TRUE)
    } else {
      write.csv(df, path, row.names = FALSE, quote = FALSE)
    }
    invisible(path)
  }

  # Standard weights — all SNP columns present, full allele info
  # Note: strand_flip entry uses G/T so complement matching is required
  # Note: mismatch entry uses C/T so allele mismatch triggers exclusion
  # Note: no_allele_info entry has blank alleles → excluded (no allele info)
  weights_standard <- data.frame(
    rsid            = c("snp_slash_clean", "snp_slash_null",
                        "snp_pipe",        "snp_concat",
                        "snp_strand_flip", "snp_ambiguous_at",
                        "snp_ambiguous_cg","snp_numeric_dos",
                        "snp_multiallelic","snp_all_missing",
                        "snp_monomorphic", "snp_mismatch",
                        "snp_high_missing","snp_partial_missing"),
    effect_allele   = c("A",  "C",  "A",  "A",  "G",  "A",
                        "C",  "",   "A",  "A",  "A",  "C",
                        "A",  "C"),
    other_allele    = c("G",  "T",  "C",  "G",  "T",  "T",
                        "G",  "",   "G",  "G",  "A",  "T",
                        "G",  "T"),
    effect_weight   = c( 0.25, 0.18, 0.30, 0.22, 0.15, 0.12,
                         0.08, 0.20, 0.10, 0.05, 0.11, 0.19,
                         0.14, 0.09),
    chr_name        = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
    chr_position    = c(100000, 200000, 300000, 400000, 500000,
                        600000, 700000, 800000, 900000, 1000000,
                        1100000, 1200000, 1300000, 1400000),
    stringsAsFactors = FALSE
  )
  write_weights(weights_standard, "weights_standard.csv")

  # Partial weights — only the clean SNPs; snp_partial_only is in the file
  # but NOT in the standard set, to test the "in catalog but not selected" path.
  # snp_slash_clean and snp_slash_null are omitted to test coverage < 100%.
  weights_partial <- data.frame(
    rsid          = c("snp_pipe", "snp_concat", "snp_strand_flip",
                      "snp_partial_only", "snp_not_in_data"),
    effect_allele = c("A", "A", "G", "A", "C"),
    other_allele  = c("C", "G", "T", "G", "T"),
    effect_weight = c(0.30, 0.22, 0.15, 0.18, 0.25),
    stringsAsFactors = FALSE
  )
  write_weights(weights_partial, "weights_partial.csv")

  # Negative weights — tests max_possible denominator with mixed-sign weights
  weights_negweights <- data.frame(
    rsid          = c("snp_slash_clean", "snp_slash_null",
                      "snp_pipe",        "snp_concat",
                      "snp_strand_flip", "snp_ambiguous_at"),
    effect_allele = c("A",  "C",  "A",  "A",  "G",  "A"),
    other_allele  = c("G",  "T",  "C",  "G",  "T",  "T"),
    effect_weight = c( 0.30, -0.15,  0.22, -0.10,  0.18,  0.05),
    stringsAsFactors = FALSE
  )
  write_weights(weights_negweights, "weights_negweights.csv")

  # No-alleles weights — effect_allele and other_allele are blank
  # All SNPs should be excluded with "no allele info" status
  weights_noalleles <- data.frame(
    rsid          = c("snp_slash_clean", "snp_slash_null", "snp_pipe"),
    effect_allele = c("", "", ""),
    other_allele  = c("", "", ""),
    effect_weight = c(0.25, 0.18, 0.30),
    stringsAsFactors = FALSE
  )
  write_weights(weights_noalleles, "weights_noalleles.csv")

  # PGS Catalog format with # header metadata lines
  header_lines <- c(
    "#pgs_id=PGS000001",
    "#pgs_name=TestScore",
    "#trait_reported=Test trait",
    "#weight_type=beta",
    "#genome_build=GRCh38",
    "#variants_number=3"
  )
  weights_header <- data.frame(
    rsid          = c("snp_slash_clean", "snp_pipe", "snp_concat"),
    effect_allele = c("A", "A", "A"),
    other_allele  = c("G", "C", "G"),
    effect_weight = c(0.25, 0.30, 0.22),
    chr_name      = c(1, 3, 4),
    chr_position  = c(100000, 300000, 400000),
    stringsAsFactors = FALSE
  )
  write_weights(weights_header, "weights_header.csv",
                header_lines = header_lines)

  # ── Print summary ──────────────────────────────────────────────────────────

  message("─── SNPstats PGS test dataset created (n = ", n, ") ───")
  message("")
  message("SNP columns and expected QC outcomes:")
  snp_outcomes <- data.frame(
    column = c("snp_slash_clean", "snp_slash_null", "snp_pipe",
               "snp_concat", "snp_strand_flip", "snp_ambiguous_at",
               "snp_ambiguous_cg", "snp_numeric_dos", "snp_multiallelic",
               "snp_all_missing", "snp_monomorphic", "snp_mismatch",
               "snp_high_missing", "snp_partial_missing", "snp_partial_only"),
    format   = c("A/G slash", "A/G + 0/0 nulls", "A|C pipe",
                 "AG concat", "C/A (flip of G/T)", "A/T ambiguous",
                 "C/G ambiguous", "0/1/2 numeric", "A/G/C multiallelic",
                 "all NA", "A/A monomorphic", "A/G (file: C/T)",
                 "A/G 18% NA", "C/T 30% NA", "A/G"),
    expected = c("\u2705 clean pass", "\u26A0 null-allele fix", "\u26A0 strand ok",
                 "\u2705 concat ok", "\u26A0 strand flip", "\u26A0 ambiguous AT",
                 "\u26A0 ambiguous CG", "\u26A0 no allele QC", "\u274c multiallelic",
                 "\u274c all missing", "\u274c monomorphic", "\u274c allele mismatch",
                 "\u274c miss filter (>10%)", "\u26A0 high missingness", "\u274c not in weights"),
    stringsAsFactors = FALSE
  )
  print(snp_outcomes, row.names = FALSE)

  message("")
  message("Response columns:  resp_binary, resp_polytomous, resp_continuous, resp_smalln")
  message("Covariate columns: cov_factor (M/F), cov_numeric, cov_multicategory (A/B/C)")
  message("")
  message("Weights files written to: ", weights_dir)
  message("  weights_standard.csv   — full, all alleles")
  message("  weights_partial.csv    — subset + one phantom SNP not in data")
  message("  weights_negweights.csv — includes negative effect weights")
  message("  weights_noalleles.csv  — blank effect/other allele columns")
  message("  weights_header.csv     — PGS Catalog # header metadata")

  invisible(dat)
}


# ─────────────────────────────────────────────────────────────────────────────
# Suggested test scenarios
# (run manually after sourcing this file)
# ─────────────────────────────────────────────────────────────────────────────

pgs_test_scenarios <- list(

  # 1. Basic weighted scoring with all SNP formats — verifies format detection,
  #    allele QC routing, and coverage table completeness.
  basic_weighted = list(
    desc       = "All SNP formats, standard weights, binary response",
    snpCols    = c("snp_slash_clean","snp_slash_null","snp_pipe",
                   "snp_concat","snp_strand_flip","snp_ambiguous_at",
                   "snp_ambiguous_cg","snp_numeric_dos","snp_multiallelic",
                   "snp_all_missing","snp_monomorphic","snp_mismatch",
                   "snp_high_missing","snp_partial_missing"),
    weightsFile = "weights_standard.csv",
    mode        = "both",
    response    = "resp_binary",
    covariates  = "cov_factor",
    expected_qc = list(pass = 7, warn = 5, fail = 7)
  ),

  # 2. Negative weights — verifies that max_possible denominator excludes
  #    negatives when computing per-individual correction.
  negative_weights = list(
    desc        = "Negative effect weights — correction denominator",
    snpCols     = c("snp_slash_clean","snp_slash_null","snp_pipe",
                    "snp_concat","snp_strand_flip","snp_ambiguous_at"),
    weightsFile = "weights_negweights.csv",
    mode        = "weighted",
    response    = "resp_continuous",
    covariates  = NULL
  ),

  # 3. Partial weights — some SNPs in data have no weights file entry,
  #    one entry in the file has no match in data.
  partial_coverage = list(
    desc        = "Partial weights file — coverage < 100%, phantom SNP",
    snpCols     = c("snp_pipe","snp_concat","snp_strand_flip",
                    "snp_partial_only"),
    weightsFile = "weights_partial.csv",
    mode        = "both",
    response    = "resp_binary",
    covariates  = NULL
  ),

  # 4. No allele info in weights — all SNPs should be excluded
  no_allele_info = list(
    desc        = "Blank allele columns — all SNPs excluded with status message",
    snpCols     = c("snp_slash_clean","snp_slash_null","snp_pipe"),
    weightsFile = "weights_noalleles.csv",
    mode        = "weighted",
    response    = NULL,
    covariates  = NULL
  ),

  # 5. PGS Catalog header — verifies metadata parsing from # comment lines
  catalog_header = list(
    desc        = "PGS Catalog # header metadata parsing",
    snpCols     = c("snp_slash_clean","snp_pipe","snp_concat"),
    weightsFile = "weights_header.csv",
    mode        = "weighted",
    response    = "resp_binary",
    covariates  = NULL
  ),

  # 6. Polytomous response — exercises the multinomial branches in assocTable,
  #    percentileTable, and interactionTable.
  polytomous = list(
    desc        = "Polytomous response (3 levels) — multinom in all analysis tables",
    snpCols     = c("snp_slash_clean","snp_pipe","snp_concat",
                    "snp_strand_flip","snp_partial_missing"),
    weightsFile = "weights_standard.csv",
    mode        = "both",
    response    = "resp_polytomous",
    covariates  = c("cov_factor", "cov_numeric")
  ),

  # 7. Small n binary — triggers near-separation, exercises .exp_or clamping
  small_n_separation = list(
    desc        = "Near-separation binary (10 cases) — OR clamping to Inf",
    snpCols     = c("snp_slash_clean","snp_pipe","snp_concat"),
    weightsFile = "weights_standard.csv",
    mode        = "weighted",
    response    = "resp_smalln",
    covariates  = c("cov_factor", "cov_numeric", "cov_multicategory")
  ),

  # 8. Missing genotype strategies — all four strategies with same data
  missing_strategies = list(
    desc        = "All four missing-genotype strategies with 30%-missing SNP",
    snpCols     = c("snp_slash_clean","snp_partial_missing",
                    "snp_slash_null"),
    weightsFile = "weights_standard.csv",
    mode        = "both",
    response    = "resp_continuous",
    covariates  = "cov_numeric",
    strategies  = c("SNP-wise", "zero", "mean", "exclude")
  ),

  # 9. Continuous response with all scale/rescale options
  scale_methods = list(
    desc        = "All rescaling methods and standardise flag",
    snpCols     = c("snp_slash_clean","snp_pipe","snp_concat",
                    "snp_strand_flip"),
    weightsFile = "weights_standard.csv",
    mode        = "both",
    response    = "resp_continuous",
    covariates  = "cov_numeric",
    scale_methods = c("none","proportion","percent","multiply","perNAlleles"),
    standardize   = c(TRUE, FALSE)
  ),

  # 10. Interaction with multicategory covariate — first covariate has 3 levels
  multicategory_interaction = list(
    desc        = "Interaction with 3-level factor covariate (cov_multicategory)",
    snpCols     = c("snp_slash_clean","snp_pipe","snp_concat"),
    weightsFile = "weights_standard.csv",
    mode        = "both",
    response    = "resp_binary",
    covariates  = c("cov_multicategory", "cov_numeric")
  ),

  # 11. HWE filter in controls (binary response) vs overall (no response)
  hwe_filter = list(
    desc        = "HWE filter — tested in controls only when binary response given",
    snpCols     = c("snp_slash_clean","snp_slash_null","snp_pipe",
                    "snp_concat","snp_strand_flip"),
    weightsFile = "weights_standard.csv",
    mode        = "weighted",
    response    = "resp_binary",
    covariates  = NULL,
    qc_hwe      = TRUE,
    hwe_thresh  = 0.10   # intentionally lenient to not filter too aggressively
  ),

  # 12. Scale weights to unit mean — verifies weighted ≈ unweighted numerically
  scale_weights = list(
    desc        = "Scale weights to unit mean — weighted score comparable to unweighted",
    snpCols     = c("snp_slash_clean","snp_pipe","snp_concat",
                    "snp_strand_flip","snp_ambiguous_at"),
    weightsFile = "weights_standard.csv",
    mode        = "both",
    response    = "resp_continuous",
    covariates  = NULL,
    scaleWeights = TRUE
  )
)


# ─────────────────────────────────────────────────────────────────────────────
# Quick smoke-test runner (no jamovi; tests helpers directly)
# Call:  smoke_test_helpers(dat)
# ─────────────────────────────────────────────────────────────────────────────

smoke_test_helpers <- function(dat) {

  source("R/snp_helpers.R")
  pass <- 0L; fail <- 0L

  chk <- function(desc, expr) {
    ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
    if (ok) {
      pass <<- pass + 1L
      message("  \u2705  ", desc)
    } else {
      fail <<- fail + 1L
      message("  \u274c  FAIL: ", desc)
    }
  }

  message("\n── detect_snp_sep ────────────────────────────────────────")
  chk("slash_clean detected",     !is.null(detect_snp_sep(dat$snp_slash_clean)))
  chk("pipe_sep detected",        !is.null(detect_snp_sep(dat$snp_pipe)))
  chk("concat detected",          !is.null(detect_snp_sep(dat$snp_concat)))
  chk("numeric_dos NOT a SNP",    is.null(detect_snp_sep(dat$snp_numeric_dos)))
  chk("all_missing returns NULL", is.null(detect_snp_sep(dat$snp_all_missing)))

  message("\n── is_snp_column ─────────────────────────────────────────")
  chk("slash_clean is SNP",     is_snp_column(dat$snp_slash_clean))
  chk("ambiguous_at is SNP",    is_snp_column(dat$snp_ambiguous_at))
  chk("cov_factor is not SNP",  !is_snp_column(dat$cov_factor))
  chk("cov_numeric is not SNP", !is_snp_column(dat$cov_numeric))

  message("\n── snp_af_hwe ────────────────────────────────────────────")
  r1 <- snp_af_hwe(dat$snp_slash_clean, effect_allele = "A")
  chk("AF in [0,1] for slash_clean",    r1$effect_af >= 0 && r1$effect_af <= 1)
  chk("HWE p in [0,1] for slash_clean", r1$hwe_p    >= 0 && r1$hwe_p    <= 1)

  r2 <- snp_af_hwe(dat$snp_numeric_dos, effect_allele = NULL)
  chk("AF computable from dosage",      r2$effect_af >= 0 && r2$effect_af <= 1)

  r3 <- snp_af_hwe(dat$snp_all_missing)
  chk("all_missing returns NA AF",  is.na(r3$effect_af))
  chk("all_missing returns NA HWE", is.na(r3$hwe_p))

  message("\n── .exp_or ───────────────────────────────────────────────")
  chk("normal OR preserved",    abs(.exp_or(log(2)) - 2) < 0.001)
  chk("large OR clamped to Inf", .exp_or(20) == Inf)
  chk("tiny OR clamped to 0",    .exp_or(-20) == 0)
  chk("NA preserved",            is.na(.exp_or(NA_real_)))

  message("\n── skewness ──────────────────────────────────────────────")
  chk("right-skewed > 0",  skewness(c(1,1,1,1,2,2,3,10)) > 0)
  chk("symmetric ≈ 0",     abs(skewness(rnorm(1000))) < 0.3)
  chk("n < 3 returns NA",  is.na(skewness(c(1, 2))))

  message("\n── Summary ───────────────────────────────────────────────")
  message("  Passed: ", pass, "  Failed: ", fail)
  invisible(list(pass = pass, fail = fail))
}


smoke_test_helpers()
