
# The two halves, joined.
#
# Every other test in this suite builds its own payloads, so the browser could
# slice the wrong bytes — or send them under the wrong sourceFormat — and the R
# suite would still be green. These options are the ones snpImport.js actually
# produced, dumped by tests/js/dump-payloads.js, and they are run through the
# analysis unchanged.
#
# The claim being checked is the one that matters for the formats: whatever the
# user picked, the spreadsheet ends up with the same genotypes.

payload_file <- function(tier) {
  for (p in c(file.path("..", "..", ".tmp", paste0("browser-payloads-", tier, ".tsv")),
              file.path(".tmp", paste0("browser-payloads-", tier, ".tsv")),
              file.path("..", "..", "..", ".tmp",
                        paste0("browser-payloads-", tier, ".tsv"))))
    if (file.exists(p)) return(p)
  NULL
}

read_payloads <- function(tier) {
  f <- payload_file(tier)
  if (is.null(f))
    testthat::skip(paste0("no browser payloads for '", tier,
                          "' — run node tests/js/dump-payloads.js"))
  utils::read.delim(f, colClasses = "character", quote = "")
}

# The genotypes as strings, in sample order, for every SNP column the import
# emitted. Sample and covariate columns are excluded by name: they are the same
# in every format only when the format carries them, which .tped does not.
snp_columns <- function(results, ids) {
  v <- out_values(results)
  v[names(v) %in% ids]
}

run_payload <- function(row, ids, ...) {
  run_import(list(geno = row$genoContent, bim = row$variantContent,
                  fam = row$sampleContent),
             ids, sourceFormat = row$sourceFormat, ...)
}

IDS <- c("snp2", "snp5", "snp9", "snp100")


test_that("every format the browser can send is loaded", {
  skip_without_fixture("small")
  p <- read_payloads("small")

  expect_setequal(p$set, c("bed", "tped", "ped", "vcf", "vcfgz"))
  for (i in seq_len(nrow(p)))
    expect_true(nzchar(p$genoContent[i]),
                info = paste(p$set[i], "sent nothing:", p$loadStatus[i]))
})


test_that("the browser's own payloads decode to the same genotypes in every format", {
  skip_without_fixture("small")
  p <- read_payloads("small")

  ref <- NULL
  for (i in seq_len(nrow(p))) {
    r <- run_payload(as.list(p[i, ]), IDS)
    g <- snp_columns(r, IDS)

    expect_setequal(names(g), IDS)
    got <- lapply(g[IDS], as.character)
    if (is.null(ref)) ref <- got
    else expect_identical(got, ref, info = paste("format", p$set[i]))
  }
})


test_that("a bgzipped VCF gives exactly what the plain one gives", {
  # bgzip writes one gzip member per block. Reading only the first is a partial
  # load that still reports success, so this is checked on the decoded values
  # rather than on the status line.
  skip_without_fixture("small")
  p <- read_payloads("small")
  a <- p[p$set == "vcf", ]
  b <- p[p$set == "vcfgz", ]
  skip_if(nrow(a) == 0 || nrow(b) == 0, "no VCF payloads")

  expect_identical(b$genoContent, a$genoContent)
  expect_identical(
    lapply(snp_columns(run_payload(as.list(b), IDS), IDS), as.character),
    lapply(snp_columns(run_payload(as.list(a), IDS), IDS), as.character))
})


test_that("the samples come back in the same order whatever the format carried", {
  # .tped and .bed name their samples in a .tfam/.fam, a .ped carries them
  # inline, and a VCF names them in its header. All three must land on the same
  # rows, or the genotypes would be right and attached to the wrong people.
  skip_without_fixture("small")
  p <- read_payloads("small")

  fam <- read_fam_file("small")
  for (i in seq_len(nrow(p))) {
    r <- run_payload(as.list(p[i, ]), IDS, emitSamples = TRUE)
    v <- out_values(r)
    expect_identical(as.character(v$IID), fam$iid, info = p$set[i])
  }
})
