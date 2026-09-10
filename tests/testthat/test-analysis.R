
# The analysis end to end, driven exactly as jamovi drives it: options carrying
# base64 payloads, output columns enabled, results read back.
#
# Two jamovi-specific workarounds, both found by running against a real jamovi,
# without which this cannot run at all:
#
#   * data = data.frame(). An analysis with a Data option and no variable
#     options cannot be constructed on a non-empty frame: jmvcore's select()
#     builds a 0-column frame from an empty varsRequired and then assigns the
#     original row names to it.
#   * OptionOutput is not a parameter of the generated function, so a scripted
#     call always leaves the output disabled and emits nothing. Enabling it
#     means reaching into the private option object.

test_that("an import produces genotype columns matching plink", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[c(1, 5, 10, 20, 40)]
  p   <- make_payloads("small", ids)

  r <- run_import(p, ids)
  v <- out_values(r)

  # sample columns first, then one per SNP. No FID here: plink2 --dummy gives
  # every sample FID "0", so the column would carry no information and is
  # suppressed — see the FID rule test below.
  expect_true(all(c("IID", "sex", "phenotype") %in% names(v)))
  expect_false("FID" %in% names(v))
  expect_true(all(ids %in% names(v)))
  expect_length(v[[ids[1]]], p$n)

  o <- read_raw_oracle("small")
  for (id in ids) {
    j     <- match(id, bim$id)
    oj    <- match(id, o$ids)
    ours  <- v[[id]]
    expect_s3_class(ours, "factor")

    # Labels are canonical (alleles sorted), not in .bim order, so the level
    # set is checked rather than the level sequence.
    al <- sort(c(bim$a1[j], bim$a2[j]))
    expect_setequal(levels(ours),
                    c(paste0(al[1], "/", al[1]),
                      paste0(al[1], "/", al[2]),
                      paste0(al[2], "/", al[2])))

    # Count the allele plink counted, straight out of the genotype string.
    # Reading it off the level *position* would assume a label order that is
    # deliberately not guaranteed.
    lab  <- as.character(ours)
    dose <- vapply(strsplit(lab, "/", fixed = TRUE),
                   function(p) sum(p == o$counted[oj]), 0L)
    dose[is.na(lab)] <- NA_integer_

    theirs <- o$values[, oj]
    expect_identical(is.na(dose), unname(is.na(theirs)))
    both <- !is.na(dose) & !is.na(theirs)
    expect_equal(sum(dose[both] != theirs[both]), 0)
  }
})

test_that("the summary table reports one row per requested SNP", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:12]
  r <- run_import(make_payloads("small", ids), ids)

  df <- r$summary$asDF
  expect_equal(nrow(df), 12)
  expect_identical(as.character(df$snp), ids)
  expect_true(all(df$n + df$missing == length(read_fam_file("small")$iid)))
  expect_true(all(df$maf >= 0 & df$maf <= 0.5, na.rm = TRUE))
  expect_true(all(df$status == "ok"))
})

test_that("sample columns can be suppressed", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  v <- out_values(run_import(make_payloads("small", ids), ids,
                             emitSamples = FALSE))

  expect_false(any(c("FID", "IID", "sex", "phenotype") %in% names(v)))
  expect_identical(sort(names(v)), sort(ids))
})

test_that("dosage mode emits numeric columns", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  v <- out_values(run_import(make_payloads("small", ids), ids, dosage = TRUE))

  expect_type(v[[ids[1]]], "double")
  expect_true(all(v[[ids[1]]] %in% c(0, 1, 2, NA)))
})

test_that("the allele-label note is not shown for dosage columns", {
  skip_without_fixture("small")

  # The note explains that the spreadsheet writes A/A, A/G, G/G where this
  # table counts G/G, G/A, A/A. In dosage mode the spreadsheet holds numbers
  # and no labels at all, so the note describes something that is not there.
  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  p   <- make_payloads("small", ids)

  geno <- note_text(run_import(p, ids, dosage = FALSE)$summary)
  expect_true(any(grepl("alphabetical order", geno)))

  dose <- note_text(run_import(p, ids, dosage = TRUE)$summary)
  expect_false(any(grepl("alphabetical order", dose)))
  # the notes that still apply are untouched
  expect_true(any(grepl("HWE is the exact test", dose)))
})

test_that("QC filters drop SNPs and say why", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:30]
  p   <- make_payloads("small", ids)

  # a MAF floor of 0.5 can keep almost nothing
  r <- run_import(p, ids, applyFilters = TRUE, minMaf = 0.5,
                  maxMissing = 100, hweP = 0)
  df <- r$summary$asDF

  # The summary lists what was imported, not what was asked for: a table
  # showing SNPs that are not in the spreadsheet reads as a bug. The count that
  # was dropped goes in the import report instead.
  v <- out_values(r)
  snp_cols <- setdiff(names(v), c("FID", "IID", "sex", "phenotype"))
  expect_equal(nrow(df), length(snp_cols))
  expect_lt(nrow(df), 30)
  expect_setequal(as.character(df$snp), snp_cols)

  pv <- r$provenance$asDF
  expect_true(any(grepl("SNPs kept after filters", pv$item)))
  expect_match(pv$value[pv$item == "SNPs kept after filters"],
               sprintf("^%d of 30$", length(snp_cols)))
})

test_that("unmatched IDs are reported rather than silently ignored", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- c(bim$id[1:3], "rs_does_not_exist_1", "rs_does_not_exist_2")
  p   <- make_payloads("small", ids)          # only the 3 real ones get sliced

  r  <- run_import(p, ids)
  pv <- r$provenance$asDF
  expect_true(any(grepl("not found", pv$item)))
  expect_true(any(grepl("rs_does_not_exist", pv$value)))
})

test_that("a weights file supplies effect alleles and flags orientation", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:6]
  # effect allele = allele 1 for every SNP, i.e. the opposite of what our
  # dosage counts, so all six must come back flagged as swapped
  w <- c("rsID\teffect_allele\tother_allele\teffect_weight",
         paste(ids, bim$a1[1:6], bim$a2[1:6], "0.1", sep = "\t"))

  r <- run_import(make_payloads("small", ids), ids,
                  snpListContent = as_payload(w), snpListText = "")
  df <- r$summary$asDF
  expect_true(all(nzchar(as.character(df$allele))))
  expect_true(all(grepl("swapped|ambiguous", df$status)))
})

test_that("with a weights file the frequency column reports the effect allele", {
  skip_without_fixture("small")

  # A PGS weight belongs to a named allele, so the number worth checking is how
  # common that allele is — which is 1 - MAF exactly when the effect allele is
  # the major one, and reading 0.12 for an allele carried by 88% of samples is
  # how a score gets checked against the wrong figure.
  bim <- read_bim_file("small")
  ids <- bim$id[1:6]
  plain <- run_import(make_payloads("small", ids), ids)
  maf   <- plain$summary$asDF$maf

  # effect allele = allele 1, the one the dosage does not count
  w <- c("rsID\teffect_allele\tother_allele\teffect_weight",
         paste(ids, bim$a1[1:6], bim$a2[1:6], "0.1", sep = "\t"))
  r  <- run_import(make_payloads("small", ids), ids,
                   snpListContent = as_payload(w), snpListText = "")
  df <- r$summary$asDF

  expect_identical(r$summary$getColumn("maf")$title, "EAF")
  expect_identical(plain$summary$getColumn("maf")$title, "MAF")
  expect_true(any(grepl("effect allele named in the weights file",
                        note_text(r$summary), fixed = TRUE)))
  expect_true(any(grepl("whichever allele is rarer",
                        note_text(plain$summary), fixed = TRUE)))

  # and the same list with the effect allele the other way round
  w2 <- c("rsID\teffect_allele\tother_allele\teffect_weight",
          paste(ids, bim$a2[1:6], bim$a1[1:6], "0.1", sep = "\t"))
  r2 <- run_import(make_payloads("small", ids), ids,
                   snpListContent = as_payload(w2), snpListText = "")
  d2 <- r2$summary$asDF

  # the two orientations are the two allele frequencies of the same site, and
  # the rarer of them is the MAF the plain run reported
  expect_equal(df$maf + d2$maf, rep(1, length(ids)), tolerance = 1e-12)
  expect_equal(pmin(df$maf, d2$maf), maf, tolerance = 1e-12)
  # ...and it is not simply the same number relabelled
  expect_true(any(abs(df$maf - maf) > 1e-8))
})

test_that("an effect allele on the other strand is still the same allele", {
  expect_equal(effect_freq("A", "A", "G", 0.3), 0.7)   # effect = allele 1
  expect_equal(effect_freq("G", "A", "G", 0.3), 0.3)   # effect = allele 2
  expect_equal(effect_freq("C", "A", "G", 0.3), 0.3)   # complement of G
  expect_equal(effect_freq("T", "A", "G", 0.3), 0.7)   # complement of A
  expect_true(is.na(effect_freq("A", "C", "G", 0.3)))  # neither, either strand
  expect_true(is.na(effect_freq(NA, "A", "G", 0.3)))
  expect_true(is.na(effect_freq("A", "A", "G", NA)))
})

test_that("a mismatched trio is refused, not half-read", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:5]
  p   <- make_payloads("small", ids)

  # .bim says five variants, the .bed carries four blocks: the structural check
  # that stops a mismatched trio producing plausible wrong genotypes
  short <- base64enc::base64decode(p$bed)
  bpv   <- bed_bytes_per_variant(p$n)
  p$bed <- as_payload(short[seq_len(4 * bpv)])

  r <- run_import(p, ids)
  expect_true(r$notice$visible)
  expect_match(r$notice$content, "expected", ignore.case = TRUE)
  # .init creates the SNP rows from the .bim before .run reads the .bed, so the
  # table is hidden rather than left showing SNPs with no statistics
  expect_false(r$summary$visible)
  expect_false(r$provenance$visible)
})

test_that("a .fam from another build of the dataset is refused, not decoded", {
  skip_without_fixture("small")

  # The failure the byte-count check is blind to. Nothing about the payload is
  # wrong: the browser computed ceil(196/4) = 49 bytes a variant from a stale
  # .fam, sliced 4 whole blocks at that stride, and R re-derives the same
  # expectation. The bytes are simply the wrong bytes — measured at ~30% of
  # calls disagreeing with plink, with no error anywhere before this check.
  bim <- read_bim_file("small")
  ids <- bim$id[1:4]
  fam <- readLines(fx("small", ".fam"))[1:196]

  p <- list(bed = as_payload(slice_bed("small", match(ids, bim$id), 196)),
            bim = as_payload(readLines(fx("small", ".bim"))[match(ids, bim$id)]),
            fam = as_payload(fam), n = 196)

  # without the attestation it decodes, which is why the browser sends one
  r <- run_import(p, ids)
  expect_true(length(out_values(r)) > 0)

  # with it, the .bed's own size settles the argument
  bed_bytes <- file.size(fx("small", ".bed"))
  r <- run_import(p, ids,
                  sourceDims = paste(bed_bytes, length(bim$id), 196, sep = ","))
  expect_true(r$notice$visible)
  expect_match(r$notice$content, "not from the same dataset")
  expect_length(out_values(r), 0)
})

test_that("a .bed loaded before the trio check says it is unverified", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:4]
  r <- run_import(make_payloads("small", ids), ids, showProvenance = TRUE)
  prov <- r$provenance$asDF
  expect_true("File consistency" %in% prov$item)
  expect_match(prov$value[prov$item == "File consistency"], "not checked")

  # and stops saying so once the browser has attested to it
  bed_bytes <- file.size(fx("small", ".bed"))
  r <- run_import(make_payloads("small", ids), ids, showProvenance = TRUE,
                  sourceDims = paste(bed_bytes, length(bim$id),
                                     length(read_fam_file("small")$iid), sep = ","))
  expect_false("File consistency" %in% r$provenance$asDF$item)
})

test_that("the fast row builder is proved against jmvcore's own, not assumed", {
  ns <- asNamespace("SNPstats")

  # .addRows writes into a Table's private state. Testing that the fields exist
  # proves nothing about what they mean now, so the fast path is compared
  # against addRow() once per session and demoted if they disagree. If this ever
  # returns FALSE against a supported jmvcore, the tables are still correct --
  # they are just built the slow way -- which is the point of having it.
  expect_true(ns$.addRowsUsable())

  # and the two really do agree on a table with mixed column types
  mk <- function() {
    t <- jmvcore::Table$new(options = jmvcore::Options$new(), name = "t")
    t$addColumn(name = "a", type = "text")
    t$addColumn(name = "b", type = "number")
    t
  }
  vals <- list(a = c("x", "y", "z"), b = c(1.5, 2.5, 3.5))
  slow <- mk()
  for (i in 1:3) slow$addRow(rowKey = i, values = lapply(vals, function(v) v[i]))
  fast <- mk()
  ns$.add_rows_fast(fast, as.list(1:3), vals)
  expect_identical(as.data.frame(slow), as.data.frame(fast))

  # a column with no values supplied still gets its cells
  bare <- mk()
  ns$.add_rows_fast(bare, as.list(1:3), list(a = c("x", "y", "z")))
  expect_equal(bare$rowCount, 3)
})

test_that("a sample filter that keeps nobody is refused, not written empty", {
  skip_without_fixture("small")

  # Every call missing, so --mind at 0% removes every sample. It used to emit
  # nine columns of length zero with no notice at all, and a summary table
  # listing every SNP at N = 0 -- an import that looks like it half-worked.
  bim <- read_bim_file("small")
  ids <- bim$id[1:4]
  p   <- make_payloads("small", ids)
  p$bed <- as_payload(rep(as.raw(0x55), bed_bytes_per_variant(p$n) * 4))

  r <- run_import(p, ids, filterSamples = TRUE, maxIndMissing = 0)
  expect_true(r$notice$visible)
  expect_match(r$notice$content, "Every sample was dropped")
  expect_false(r$summary$visible)
  expect_length(out_values(r), 0)

  # the report still explains itself
  pv <- r$provenance$asDF
  expect_match(pv$value[pv$item == "Samples"], "0 kept of 200")
})

test_that("calls in an unexpected code are counted in the import report", {
  # Silence would leave them looking like ordinary no-calls, and the difference
  # decides whether the file is worth re-exporting.
  tfam <- c("f1 s1 0 0 1 1", "f1 s2 0 0 2 1", "f1 s3 0 0 1 2")
  r <- run_import(
    list(geno = as_payload(c("1 snpA 0 100 A A A G X Y",
                             "1 snpB 0 200 G G A G A A")),
         bim = as_payload(character(0)), fam = as_payload(tfam)),
    c("snpA", "snpB"), sourceFormat = "tped", showSummary = FALSE)

  df <- r$provenance$asDF
  row <- df$value[df$item == "Unreadable calls"]
  expect_length(row, 1)
  expect_match(row, "^1 set to missing")

  # and it says nothing at all when every call is readable
  clean <- run_import(
    list(geno = as_payload(c("1 snpA 0 100 A A A G G G")),
         bim = as_payload(character(0)), fam = as_payload(tfam)),
    "snpA", sourceFormat = "tped", showSummary = FALSE)
  expect_false("Unreadable calls" %in% clean$provenance$asDF$item)
})

test_that("a browser refusal is restated in the results, not just the status box", {
  # The panel refuses an oversized selection and leaves the carriers empty,
  # which from R is indistinguishable from an import nobody has started. Only
  # loadProblem tells the two apart, and without it the reason lived solely in
  # a one-line text box that shows about four words of it.
  msg <- paste("50000 SNPs x 2000 samples is too much to transfer -- jamovi",
               "would drop the columns without an error. Use at most 1500",
               "SNPs at this sample size.")
  r <- run_import(list(geno = "", bim = "", fam = ""), character(0),
                  loadStatus = msg, loadProblem = "error")
  expect_true(r$notice$visible)
  expect_match(r$notice$content, "too much to transfer", fixed = TRUE)
  expect_match(r$notice$content, "Use at most 1500", fixed = TRUE)
  expect_match(r$notice$content, "#d93025", fixed = TRUE)   # the red border
  expect_true(r$instructions$visible)

  # ...and a load nobody has attempted yet stays quiet.
  q <- run_import(list(geno = "", bim = "", fam = ""), character(0),
                  loadStatus = "ready - press Load genotypes", loadProblem = "")
  expect_false(isTRUE(q$notice$visible))
})

test_that("FID is emitted when the .fam actually distinguishes families", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  p   <- make_payloads("small", ids)

  # rewrite the .fam with real family IDs, as a pedigree file would have
  fam_lines <- readLines(fx("small", ".fam"))
  parts <- strsplit(trimws(fam_lines), "[ \t]+")
  fam_lines <- vapply(seq_along(parts), function(i) {
    f <- parts[[i]]; f[1] <- paste0("FAM", i %% 7)
    paste(f, collapse = " ")
  }, "")
  p$fam <- as_payload(fam_lines)

  v <- out_values(run_import(p, ids))
  expect_true("FID" %in% names(v))
  expect_true(nlevels(v[["FID"]]) > 1)
})

test_that("FID is emitted only when it carries information", {
  # the .fam identifies a sample by the (FID, IID) pair, but an all-identical
  # FID column is noise, so it is emitted only when it distinguishes samples
  informative <- read_fam(c("F1 S1 0 0 1 1", "F2 S2 0 0 2 2"))
  expect_true(informative$fid_informative)
  expect_true(informative$iid_unique)

  flat <- read_fam(c("FAM S1 0 0 1 1", "FAM S2 0 0 2 2"))
  expect_false(flat$fid_informative)

  # IIDs need only be unique *within* a family, so repeats are legal and mean
  # the pair is what identifies a row
  dup <- read_fam(c("F1 S1 0 0 1 1", "F2 S1 0 0 2 2"))
  expect_false(dup$iid_unique)
  expect_true(dup$fid_informative)
})

test_that("variants with no usable ID get one column each, named by position", {
  # A .bim built from a VCF without --set-missing-var-ids names every variant
  # '.'. Sanitising alone made them all '_', so every SNP wanted the same
  # column; jamovi would have merged them into '_' and '_ (2)'.
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:4]
  p   <- make_payloads("small", ids)

  lines <- readLines(fx("small", ".bim"))[p$idx]
  parts <- strsplit(trimws(lines), "[ \t]+")
  p$bim <- as_payload(vapply(parts, function(f) {
    f[2] <- "."; paste(f, collapse = "\t") }, ""))

  r <- run_import(p, ids)
  v <- out_values(r)
  want <- paste0(vapply(parts, `[`, "", 1L), ":",
                 vapply(parts, `[`, "", 4L))
  expect_true(all(want %in% names(v)))
  expect_equal(length(unique(want)), 4)

  pv <- r$provenance$asDF
  expect_match(as.character(pv$value[pv$item == "SNPs received"]), "chr:bp")

  # a genuine repeat of a real ID is numbered rather than merged
  expect_identical(variant_ids(list(chr = c("1", "1"), bp = c(10, 20),
                                    id = c("rs1", "rs1"))),
                   c("rs1", "rs1_2"))
})

test_that("repeated IIDs are called out in the report", {
  # The row was being written under a key .provKeys() never created, so
  # setRow skipped it and the warning silently never appeared.
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  p   <- make_payloads("small", ids)

  parts <- strsplit(trimws(readLines(fx("small", ".fam"))), "[ \t]+")
  p$fam <- as_payload(vapply(seq_along(parts), function(i) {
    f <- parts[[i]]
    f[1] <- paste0("FAM", i %% 2)         # two families
    f[2] <- paste0("S", ceiling(i / 2))   # IIDs repeat across them
    paste(f, collapse = " ")
  }, ""))

  pv <- run_import(p, ids)$provenance$asDF
  expect_true("Sample IDs" %in% as.character(pv$item))
  expect_match(as.character(pv$value[pv$item == "Sample IDs"]), "not unique")

  # and the unique case must not grow a row that says nothing
  expect_false("Sample IDs" %in%
               as.character(run_import(make_payloads("small", ids),
                                       ids)$provenance$asDF$item))
})

test_that("an empty SNP list imports every variant supplied", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:8]
  p   <- make_payloads("small", ids)

  # the browser sends whatever it sliced; with no list, that is everything
  r <- run_import(p, ids, snpListText = "")
  expect_equal(nrow(r$summary$asDF), 8)

  pv <- r$provenance$asDF
  expect_true(any(grepl("no list given", pv$value)))
})

test_that("changing a QC threshold does not need the file again", {
  skip_without_fixture("small")

  # The point of applying filters in R: the payload is already in the options,
  # so a threshold change is pure recomputation. Two runs differing only in the
  # threshold must give different results from identical inputs.
  bim <- read_bim_file("small")
  ids <- bim$id[1:30]
  p   <- make_payloads("small", ids)

  loose <- run_import(p, ids, applyFilters = TRUE, minMaf = 0,
                      maxMissing = 100, hweP = 0)
  tight <- run_import(p, ids, applyFilters = TRUE, minMaf = 0.45,
                      maxMissing = 100, hweP = 0)

  n_loose <- sum(loose$summary$asDF$status == "ok")
  n_tight <- sum(tight$summary$asDF$status == "ok")
  expect_gt(n_loose, n_tight)
})

test_that("sample filters drop individuals and shorten every column together", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:20]
  p   <- make_payloads("small", ids)

  full <- out_values(run_import(p, ids))
  n_all <- length(full[["IID"]])

  # a 0% missingness ceiling drops every sample with any missing genotype
  r <- run_import(p, ids, filterSamples = TRUE, maxIndMissing = 0)
  v <- out_values(r)

  expect_lt(length(v[["IID"]]), n_all)
  # every emitted column must be the same length, or genotypes and phenotypes
  # would be misaligned against each other
  expect_equal(length(unique(vapply(v, length, 0L))), 1L)

  pv <- r$provenance$asDF
  expect_true(any(grepl("Samples dropped", pv$item)))
  expect_true(any(grepl("missing", pv$value[pv$item == "Samples dropped"])))
  expect_true(any(grepl("kept of", pv$value)))
})

test_that("the report keeps its shape so option clicks do not blank it", {
  skip_without_fixture("small")

  # The rows are created in .init, which is what lets jamovi restore a populated
  # table on an option click instead of rendering an empty one and refilling.
  bim <- read_bim_file("small")
  ids <- bim$id[1:6]
  p   <- make_payloads("small", ids)

  base <- run_import(p, ids)
  expect_equal(nrow(base$summary$asDF), 6)
  expect_true(all(nzchar(as.character(base$provenance$asDF$value))))

  # the same options must always give the same row set, whatever the data says
  no_missing <- run_import(p, c(ids, "rs_not_here"))
  expect_identical(as.character(base$provenance$asDF$item),
                   as.character(no_missing$provenance$asDF$item))

  # switching a filter on adds its row, and never leaves a value blank
  filt <- run_import(p, ids, applyFilters = TRUE, filterSamples = TRUE)
  expect_true("SNPs kept after filters" %in% as.character(filt$provenance$asDF$item))
  expect_true("Samples dropped" %in% as.character(filt$provenance$asDF$item))
  expect_true(all(nzchar(as.character(filt$provenance$asDF$value))))
})

test_that("a loaded import not yet opened says to press Open as new dataset", {
  skip_without_fixture("small")

  # Nothing is produced until the button is pressed (openNew defaults FALSE,
  # like every option), so this is the one way to stop a user loading data and
  # seeing nothing happen without explanation.
  bim <- read_bim_file("small")
  ids <- bim$id[1:3]
  p   <- make_payloads("small", ids)

  ns   <- asNamespace("SNPstats")
  opts <- ns$snpImportOptions$new(
    genoContent = p$bed, variantContent = p$bim, sampleContent = p$fam,
    snpListText = paste(ids, collapse = "\n"))
  a <- ns$snpImportClass$new(options = opts, data = data.frame())
  a$run()

  expect_true(a$results$notice$visible)
  expect_match(a$results$notice$content, "Open as new dataset")

  # and once it has fired (openNewState's cell, not the option -- see the
  # r.yaml comment on why), the reminder stops
  a$results$openNewState$setRow(rowNo = 1, values = list(fired = "yes"))
  a$run()
  expect_false(isTRUE(a$results$notice$visible)
               && grepl("Open as new dataset", a$results$notice$content))
})

test_that("heterozygosity outlier filtering runs and is flagged as noisy", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:50]
  p   <- make_payloads("small", ids)

  r <- run_import(p, ids, filterSamples = TRUE, maxIndMissing = 100, hetSd = 1)
  v <- out_values(r)
  expect_equal(length(unique(vapply(v, length, 0L))), 1L)

  pv <- r$provenance$asDF
  if (any(grepl("heterozygosity", pv$item)))
    expect_true(any(grepl("unreliable", pv$value)))
})

test_that("SNP statistics describe the samples that survived filtering", {
  skip_without_fixture("small")

  # PLINK's order, and the reason it matters: a MAF computed over samples that
  # were then dropped would describe a set that no longer exists.
  bim <- read_bim_file("small")
  ids <- bim$id[1:15]
  p   <- make_payloads("small", ids)

  r <- run_import(p, ids, filterSamples = TRUE, maxIndMissing = 0)
  df <- r$summary$asDF
  v  <- out_values(r)

  expect_true(all(df$n + df$missing == length(v[["IID"]])))
  expect_true(all(df$missing == 0))     # every remaining sample is complete
})


test_that("relaxing a SNP filter brings the column back", {
  skip_without_fixture("small")

  # Each run builds its columns from scratch (.buildColumns, off the current
  # d$keep), so there is no stale state from a stricter previous run to worry
  # about -- unlike the growing-Output design this replaced, where a SNP a
  # filter had removed had to be deliberately re-emitted to come back.
  bim <- read_bim_file("small")
  ids <- bim$id[1:20]
  p   <- make_payloads("small", ids)

  strict <- run_import(p, ids, applyFilters = TRUE, minMaf = 0.45,
                       maxMissing = 100, hweP = 0)
  n_strict <- length(setdiff(names(out_values(strict)),
                             c("FID", "IID", "sex", "phenotype")))

  relaxed <- run_import(p, ids, applyFilters = TRUE, minMaf = 0,
                        maxMissing = 100, hweP = 0)
  kept <- setdiff(names(out_values(relaxed)), c("FID", "IID", "sex", "phenotype"))

  expect_gt(length(kept), n_strict)
  expect_setequal(kept, ids)          # everything comes back
})


test_that("HWE falls back to all samples when the .fam has no phenotype", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:5]
  p   <- make_payloads("small", ids)

  # blank out every phenotype: asking for controls only would otherwise test
  # nobody, so it has to fall back and say so
  fam_lines <- readLines(fx("small", ".fam"))
  parts <- strsplit(trimws(fam_lines), "[ \t]+")
  p$fam <- as_payload(vapply(parts, function(f) {
    f[6] <- "-9"; paste(f, collapse = " ")
  }, ""))

  r <- run_import(p, ids, hweGroup = "controls")
  note <- note_text(r$summary)
  expect_true(any(grepl("all samples", note)))
  expect_true(any(grepl("no phenotype", note)))
  expect_true(all(!is.na(r$summary$asDF$hwePval)))
})

test_that("the summary note says which samples HWE used", {
  skip_without_fixture("small")

  bim <- read_bim_file("small")
  ids <- bim$id[1:5]
  p   <- make_payloads("small", ids)

  note <- note_text(run_import(p, ids, hweGroup = "controls")$summary)
  expect_true(any(grepl("computed in controls", note)))
  expect_true(any(grepl("describe all 200 imported samples", note)))

  note <- note_text(run_import(p, ids, hweGroup = "all")$summary)
  expect_true(any(grepl("computed in all samples", note)))
})
