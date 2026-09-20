
# Covariate files and the merge onto genotype samples.
#
# The merge is the part with teeth: a covariate silently attached to the wrong
# sample is worse than one that fails to attach at all, so most of these check
# alignment rather than parsing.

test_that("the ID column is detected, with IID winning over generic names", {
  # PLINK's convention, and plink2 writes it with a leading '#'
  t1 <- read_covariate_table(c("FID\tIID\tage", "F1\tS1\t50", "F2\tS2\t60"))
  expect_identical(t1$id_name, "IID")
  expect_identical(t1$fid_name, "FID")
  expect_identical(t1$id, c("S1", "S2"))
  expect_identical(names(t1$vars), "age")

  expect_identical(read_covariate_table(c("#IID\tage", "S1\t50"))$id_name, "IID")

  # IID wins even when a generic name comes first
  t2 <- read_covariate_table(c("id\tIID\tage", "x\tS1\t50"))
  expect_identical(t2$id_name, "IID")

  # and the generic ones are accepted when IID is absent
  for (nm in c("id", "sample", "subject"))
    expect_identical(
      read_covariate_table(c(paste0(nm, "\tage"), "S1\t50"))$id_name, nm)
})

test_that("a file with no recognisable header falls back to the first column", {
  t <- read_covariate_table(c("S0\t1\t2", "S1\t3\t4"))
  expect_identical(t$id, c("S1"))        # row 1 is treated as the header
  expect_length(t$vars, 2)
})

test_that("comment lines are skipped before the header is found", {
  # The delimiter used to be guessed from line 1. A '# free text' comment has
  # no tabs, so a perfectly good file reported "only one column" — which is
  # exactly what a PGS Catalog-style file looks like.
  t <- read_covariate_table(c("# exported by our pipeline",
                              "# second comment line",
                              "IID\tage\tbmi",
                              "S1\t50\t24.0",
                              "S2\t60\t27.5"))
  expect_identical(t$id_name, "IID")
  expect_identical(names(t$vars), c("age", "bmi"))
  expect_identical(t$id, c("S1", "S2"))
})

test_that("a '#' header is kept while '#' comments are dropped", {
  # plink2 writes the header itself with a leading '#', so the two conventions
  # disagree about what '#' means. They are told apart by field count.
  p2 <- read_covariate_table(c("#IID\tage", "S1\t50", "S2\t60"))
  expect_identical(p2$id_name, "IID")
  expect_identical(p2$id, c("S1", "S2"))

  # a '#' line that is *not* header-shaped is a comment, so the next line is
  both <- read_covariate_table(c("## meta", "# note", "#FID\tIID\tage",
                                 "F1\tS1\t50"))
  expect_identical(both$id_name, "IID")
  expect_identical(both$fid_name, "FID")
  expect_identical(both$id, "S1")
})

test_that("a PGS weights file in the covariate slot is named, not mis-parsed", {
  # It lists SNPs, not samples. Falling back to "column 1 is the ID" would
  # match rsIDs against IIDs and quietly produce a table of NAs.
  w <- c("# PGS Catalog scoring file",
         "rsID\teffect_allele\tother_allele\teffect_weight",
         "rs1\tA\tG\t0.1")
  expect_error(read_covariate_table(w), "PGS weights file")
  expect_error(read_covariate_table(w), "SNPs to import")

  expect_error(read_covariate_table(c("variant_id,beta", "rs1,0.2")), "weights")
})

test_that("the real weights fixture is refused as a covariate file", {
  skip_without_fixture("small")
  expect_error(
    read_covariate_table(readLines(file.path(fixture_dir("small"), "weights.tsv"))),
    "PGS weights file")
})

test_that("the delimiter is guessed, and a one-column file is refused", {
  for (d in c("\t", ",", ";")) {
    t <- read_covariate_table(c(paste("IID", "age", sep = d),
                                paste("S1", "50", sep = d)))
    expect_identical(t$id, "S1")
    expect_identical(names(t$vars), "age")
  }
  expect_error(read_covariate_table(c("IID", "S1")), "one column")
  expect_error(read_covariate_table("IID\tage"), "no data rows")
})

test_that("an explicit ID column overrides detection, and a wrong one errors", {
  lines <- c("IID\talt_id\tage", "S1\tX1\t50", "S2\tX2\t60")
  expect_identical(read_covariate_table(lines, "alt_id")$id, c("X1", "X2"))
  expect_error(read_covariate_table(lines, "nope"), "no column named")
})

test_that("columns are typed numeric only when every value parses", {
  expect_type(cov_typed(c("1", "2.5", "")), "double")
  expect_true(is.na(cov_typed(c("1", "NA", "-9"))[2]))

  # one stray label keeps the whole column categorical rather than quietly
  # discarding the value
  expect_s3_class(cov_typed(c("1", "2", "high")), "factor")
  expect_s3_class(cov_typed(c("yes", "no")), "factor")
})

test_that("the merge keeps the genotype order and reports both directions", {
  fam <- read_fam(c("F1 S1 0 0 1 1", "F2 S2 0 0 2 2", "F3 S3 0 0 1 1"))
  # S2 is absent, and there is a covariate row for a sample that is not genotyped
  tab <- read_covariate_table(c("IID\tage", "S3\t30", "S1\t10", "GHOST\t99"))
  m   <- merge_covariates(tab, fam)

  # values follow the .fam, not the covariate file
  expect_equal(as.numeric(m$values$age), c(10, NA, 30))
  expect_equal(m$matched, 2)
  expect_equal(m$n_geno, 3)
  expect_identical(m$unmatched_geno, "S2")
  expect_identical(m$unmatched_cov, "GHOST")
  expect_identical(m$matched_on, "IID")
})

test_that("FID and IID together identify a sample when IIDs repeat", {
  # the .fam spec only requires IID to be unique within a family
  fam <- read_fam(c("F1 S1 0 0 1 1", "F2 S1 0 0 2 2"))
  expect_false(fam$iid_unique)

  tab <- read_covariate_table(c("FID\tIID\tage", "F2\tS1\t20", "F1\tS1\t10"))
  m   <- merge_covariates(tab, fam)

  expect_identical(m$matched_on, "FID + IID")
  expect_equal(as.numeric(m$values$age), c(10, 20))   # not c(20, 10)
})

test_that("nothing matching is reported rather than silently producing NAs", {
  fam <- read_fam(c("F1 S1 0 0 1 1", "F2 S2 0 0 2 2"))
  m <- merge_covariates(read_covariate_table(c("IID\tage", "0001\t10")), fam)
  expect_equal(m$matched, 0)
  expect_length(m$unmatched_geno, 2)
  expect_true(all(is.na(m$values$age)))
})

test_that("duplicate covariate IDs are flagged", {
  fam <- read_fam(c("F1 S1 0 0 1 1"))
  m <- merge_covariates(
    read_covariate_table(c("IID\tage", "S1\t10", "S1\t20")), fam)
  expect_true(m$dup_cov)
  expect_equal(as.numeric(m$values$age), 10)     # first match wins
})

test_that("covariate names are made unique against what is already emitted", {
  expect_identical(cov_unique_names(c("age", "sex"), c("IID", "sex")),
                   c("age", "sex_2"))
  expect_identical(cov_unique_names(c("a", "a"), character(0)), c("a", "a_2"))
  expect_identical(cov_unique_names("a`b", character(0)), "a_b")

  # the first of a pair keeps its name; only the later ones move
  expect_identical(cov_unique_names(c("a", "a", "a"), character(0)),
                   c("a", "a_2", "a_3"))
  # a suffix that is itself taken keeps counting
  expect_identical(cov_unique_names("a", c("a", "a_2")), "a_3")
  # type-stable when there is nothing to name
  expect_identical(cov_unique_names(character(0), "a"), character(0))
  expect_identical(sanitise_id(character(0)), character(0))
})

test_that("naming a large selection is linear, not quadratic", {
  # This runs inside .prepared(), so it is paid on every option click, and the
  # cell budget allows ~30 000 variants at 100 samples. The loop it replaced
  # rebuilt a hash table per name and took 2.8 s at that size.
  nm <- paste0("rs", seq_len(30000))
  t <- system.time(out <- cov_unique_names(nm, c("FID", "IID")))[["elapsed"]]
  expect_identical(out, nm)
  expect_lt(t, 1)
})


# ── through the analysis ─────────────────────────────────────────────────────

# covFile is a native File option now, so the "payload" is a real path.
cov_payload <- function(lines) {
  f <- tempfile(fileext = ".tsv")
  writeLines(lines, f)
  f
}

test_that("covariates are emitted between the phenotype and the SNPs", {
  skip_without_fixture("small")

  # Position matters: jamovi fixes a column's place when it is first created,
  # so the emit order is what puts covariates in the middle on a first load.
  bim <- read_bim_file("small")
  fam <- read_fam_file("small")
  ids <- bim$id[1:4]
  p   <- make_payloads("small", ids)

  cov <- c("IID\tage\tbmi",
           paste(fam$iid, seq_along(fam$iid), 20 + seq_along(fam$iid),
                 sep = "\t"))

  v <- out_values(run_import(p, ids, covFile = cov_payload(cov)))
  nms <- names(v)

  expect_true(all(c("age", "bmi") %in% nms))
  expect_lt(max(match(c("age", "bmi"), nms)), min(match(ids, nms)))
  expect_gt(min(match(c("age", "bmi"), nms)), match("phenotype", nms))

  # and they are aligned to the genotype samples
  expect_length(v$age, length(fam$iid))
  expect_equal(as.numeric(v$age), seq_along(fam$iid))
})

test_that("covariate types survive into the emitted columns", {
  skip_without_fixture("small")

  bim <- read_bim_file("small"); fam <- read_fam_file("small")
  ids <- bim$id[1:2]
  cov <- c("IID\tage\tsmoker",
           paste(fam$iid, seq_along(fam$iid),
                 ifelse(seq_along(fam$iid) %% 2 == 0, "yes", "no"), sep = "\t"))

  v <- out_values(run_import(make_payloads("small", ids), ids,
                             covFile = cov_payload(cov)))
  expect_type(v$age, "double")
  expect_s3_class(v$smoker, "factor")
  expect_setequal(levels(v$smoker), c("yes", "no"))
})

test_that("a partly matching covariate file reports both directions", {
  skip_without_fixture("small")

  bim <- read_bim_file("small"); fam <- read_fam_file("small")
  ids <- bim$id[1:2]
  keep <- fam$iid[1:50]
  cov <- c("IID\tage", paste(c(keep, "GHOST1", "GHOST2"),
                             c(seq_along(keep), 1, 2), sep = "\t"))

  r <- run_import(make_payloads("small", ids), ids,
                  covFile = cov_payload(cov))
  pv <- r$provenance$asDF
  val <- pv$value[pv$item == "Covariates matched"]

  expect_match(val, "50 of 200")
  expect_match(val, "no covariates for 150")
  expect_match(val, "2 covariate rows unused")

  v <- out_values(r)
  expect_equal(sum(!is.na(v$age)), 50)
})

test_that("a name clash with an emitted column is renamed, not silently doubled", {
  skip_without_fixture("small")

  bim <- read_bim_file("small"); fam <- read_fam_file("small")
  ids <- bim$id[1:2]
  # 'sex' and 'phenotype' already come from the .fam
  cov <- c("IID\tsex\tage", paste(fam$iid, "M", 1, sep = "\t"))

  v <- out_values(run_import(make_payloads("small", ids), ids,
                             covFile = cov_payload(cov)))
  expect_true("sex" %in% names(v))       # the .fam's
  expect_true("sex_2" %in% names(v))     # the covariate file's
  expect_equal(sum(names(v) == "sex"), 1)
})

test_that("a covariate named after an imported SNP does not orphan a column", {
  skip_without_fixture("small")

  # Two columns named snp1 cannot both survive into the opened dataset -- the
  # second would silently overwrite or shadow the first depending on how it
  # is built. The only cure is not to emit the name twice.
  bim <- read_bim_file("small"); fam <- read_fam_file("small")
  ids <- bim$id[1:3]
  cov <- c(paste("IID", ids[2], "age", sep = "\t"),
           paste(fam$iid, seq_along(fam$iid), 40 + seq_along(fam$iid), sep = "\t"))

  r <- run_import(make_payloads("small", ids), ids, covFile = cov_payload(cov))
  v <- out_values(r)
  expect_equal(anyDuplicated(names(v)), 0)
  expect_true(all(vapply(v, length, 0L) > 0))     # none left unwritten

  expect_true(ids[2] %in% names(v))               # the covariate keeps the name
  expect_true(paste0(ids[2], "_2") %in% names(v)) # the SNP moves aside

  # and the summary table names the column the user will actually find
  expect_true(paste0(ids[2], "_2") %in% as.character(r$summary$asDF$snp))
})

test_that("covariates attach to the samples that survived filtering", {
  skip_without_fixture("small")

  # A dropped sample must not come back through the covariate file, and the
  # remaining covariate values must still line up with the right rows.
  bim <- read_bim_file("small"); fam <- read_fam_file("small")
  ids <- bim$id[1:20]
  cov <- c("IID\tage", paste(fam$iid, seq_along(fam$iid), sep = "\t"))

  r <- run_import(make_payloads("small", ids), ids,
                  covFile = cov_payload(cov),
                  filterSamples = TRUE, maxIndMissing = 0)
  v <- out_values(r)

  expect_lt(length(v$age), length(fam$iid))
  expect_equal(length(unique(vapply(v, length, 0L))), 1L)

  # every emitted age must be the one belonging to its IID
  expect_equal(as.numeric(v$age), match(as.character(v$IID), fam$iid))
})

test_that("an unreadable covariate file fails the import loudly", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:2]
  r <- run_import(make_payloads("small", ids), ids,
                  covFile = cov_payload("IID"))        # header only
  expect_true(r$notice$visible)
  expect_match(r$notice$content, "Import failed")
})
