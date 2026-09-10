
# The transport layer, including the hostile-input cases. A payload is
# attacker-controlled input in a process that also holds other users' sessions,
# so these are regression tests for the hardening, not edge-case trivia.

test_that("base64 round-trips raw bytes and text", {
  b <- as.raw(c(0x6c, 0x1b, 0x01, 0xff, 0x00))
  expect_identical(payload_bytes(as_payload(b)), b)
  expect_identical(payload_lines(as_payload(c("a", "b", "c"))), c("a", "b", "c"))
  expect_identical(payload_bytes(""), raw(0))
  expect_identical(payload_lines(""), character(0))
})

test_that("compression is detected by magic bytes, not by being told", {
  txt <- paste(rep("rs12345\tA\tG", 200), collapse = "\n")

  # The browser's CompressionStream('gzip') emits a real gzip stream (1f 8b),
  # and that is the only compressed form the transport carries — a payload that
  # is not recognised passes through undecompressed and parses as binary noise.
  gzf <- tempfile(); con <- gzfile(gzf, "wb")
  writeBin(charToRaw(txt), con); close(con)
  real_gzip <- readBin(gzf, "raw", file.info(gzf)$size)

  expect_identical(as.integer(real_gzip[1:2]), c(0x1fL, 0x8bL))
  expect_lt(length(real_gzip), nchar(txt))
  expect_identical(rawToChar(payload_bytes(as_payload(real_gzip))), txt)
  expect_length(payload_lines(as_payload(real_gzip)), 200)

  # zlib (78 xx) — what R's own memCompress(type = "gzip") writes — is refused
  # rather than inflated, because it cannot be inflated under a bound.
  zlib <- memCompress(charToRaw(txt), type = "gzip")
  expect_identical(as.integer(zlib[1]), 0x78L)
  expect_error(payload_bytes(as_payload(zlib)), "only gzip is accepted")
})

test_that("truncated, oversized and invalid payloads are refused", {
  expect_error(payload_bytes("abcde"), "truncated")
  # over MAX_PAYLOAD_BYTES (8 MB) once decoded: 8 MB * 4/3 of base64
  expect_error(payload_bytes(strrep("A", 12 * 1024 * 1024)), "refusing to decode")

  # base64enc::base64decode does not validate its input — it silently returns
  # garbage — so this checks our own alphabet guard, not the decoder's.
  expect_error(payload_bytes("!!!!"), "not valid base64")
  expect_error(payload_bytes("ab*d"), "not valid base64")
})

test_that("a decompression bomb is refused without being expanded", {
  # 512 MB of zeros compresses to a few hundred KB; expanding it would take the
  # engine down and with it every session sharing the process.
  bf <- tempfile(); con <- gzfile(bf, "wb")
  writeBin(raw(512 * 1024^2), con); close(con)
  bomb <- readBin(bf, "raw", file.info(bf)$size)
  expect_lt(length(bomb), 1024^2)

  # ISIZE makes a bomb cheap to spot, but it is a number in the payload and can
  # say anything, so the refusal must not rest on it: this copy declares 1 KB of
  # contents and must be refused on the same terms as the honest one.
  liar <- bomb
  liar[(length(liar) - 3):length(liar)] <- as.raw(c(0x00, 0x04, 0x00, 0x00))
  expect_equal(.gzip_isize(liar), 1024)

  for (v in list(bomb, liar)) {
    b64 <- as_payload(v)
    # "without being expanded" is the assertion, so it is measured rather than
    # implied: the whole refusal must cost less memory than the bomb inflates
    # to. Raw vectors live in Vcells, 8 bytes each, and gc() reports their peak.
    gc(reset = TRUE)
    expect_error(payload_bytes(b64), "expands to more than")
    expect_lt(gc()[2L, "max used"] * 8 / 1024^2, 100)
  }

  # a payload that is merely large, not a bomb, still decodes
  okf <- tempfile(); con <- gzfile(okf, "wb")
  writeBin(charToRaw(strrep("rs1\trs2\n", 1000)), con); close(con)
  ok <- readBin(okf, "raw", file.info(okf)$size)
  expect_gt(length(payload_bytes(as_payload(ok))), 1000)
})

test_that("a truncated gzip payload is refused, not silently shortened", {
  # gzcon() treats the end of a corrupt stream as the end of the data, so a
  # payload cut short inflates to a prefix and would read as a file simply
  # missing its last variants. ISIZE is the only end-of-stream evidence a gzip
  # member carries, so a short read has to be caught against it.
  txt <- paste(rep("rs12345\tA\tG", 2000), collapse = "\n")
  gzf <- tempfile(); con <- gzfile(gzf, "wb")
  writeBin(charToRaw(txt), con); close(con)
  gz <- readBin(gzf, "raw", file.info(gzf)$size)

  cut <- c(gz[1:(length(gz) - 40L)], gz[(length(gz) - 7L):length(gz)])
  expect_error(payload_bytes(as_payload(cut)), "truncated or corrupt")
})

test_that("line endings and a BOM are handled", {
  crlf <- charToRaw("a\r\nb\r\nc")
  expect_identical(payload_lines(as_payload(crlf)), c("a", "b", "c"))

  bom <- c(as.raw(c(0xef, 0xbb, 0xbf)), charToRaw("rs1\nrs2"))
  expect_identical(payload_lines(as_payload(bom)), c("rs1", "rs2"))
})

test_that("a plain ID list is read whatever shape it arrives in", {
  expect_identical(parse_snp_selection(c("rs1", "rs2", "rs3"))$ids,
                   c("rs1", "rs2", "rs3"))
  expect_identical(parse_snp_selection("rs1 rs2  rs3")$ids,
                   c("rs1", "rs2", "rs3"))
  expect_identical(parse_snp_selection(c("rs1", "", "  ", "rs2"))$ids,
                   c("rs1", "rs2"))
  expect_length(parse_snp_selection(character(0))$ids, 0)
})

test_that("the ID column is chosen by preference, as the browser chooses it", {
  # The browser walks SNP_ID_NAMES in order and slices the variants it finds;
  # R used to take the leftmost column matching any of them. On this header the
  # two disagreed: the payload held rs1/rs2 while the report said x1/x2 had been
  # requested and neither was found. Same rule, same order, both halves.
  s <- parse_snp_selection(c("id\trsid\teffect_allele",
                            "x1\trs1\tA", "x2\trs2\tG"))
  expect_identical(s$ids, c("rs1", "rs2"))
  expect_identical(s$effect_allele, c("A", "G"))

  # preference, not position, in the other direction too
  expect_identical(parse_snp_selection(c("rsid\tid", "rs1\tx1"))$ids, "rs1")
  # and the later names still work when the preferred ones are absent
  expect_identical(parse_snp_selection(c("name\tchr", "rs9\t1"))$ids, "rs9")
})

test_that("a PGS-Catalog weights file yields IDs and effect alleles", {
  w <- c("# PGS Catalog scoring file",
         "# format_version=2.0",
         "rsID\teffect_allele\tother_allele\teffect_weight",
         "rs100\tA\tG\t0.15",
         "rs200\tT\tC\t-0.30")
  s <- parse_snp_selection(w)
  expect_identical(s$ids, c("rs100", "rs200"))
  expect_identical(s$effect_allele, c("A", "T"))
  expect_equal(s$weight, c(0.15, -0.30))
})

test_that("weights files are recognised under their alternative column names", {
  s <- parse_snp_selection(c("variant_id,effect_allele,beta",
                             "rs1,C,0.2", "rs2,G,0.4"))
  expect_identical(s$ids, c("rs1", "rs2"))
  expect_identical(s$effect_allele, c("C", "G"))
  expect_equal(s$weight, c(0.2, 0.4))
})

test_that("allele orientation is classified, and palindromes are not guessed", {
  # our dosage counts allele 2, so an effect allele equal to a2 needs no change
  expect_identical(allele_match("G", "A", "G"), "match")
  expect_identical(allele_match("A", "A", "G"), "swapped")
  expect_identical(allele_match("C", "A", "G"), "flipped")   # complement of G
  expect_identical(allele_match("T", "A", "G"), "flipped")   # complement of A
  expect_identical(allele_match("X", "A", "G"), "absent")
  expect_identical(allele_match(NA,  "A", "G"), "absent")

  # A/T and C/G cannot be resolved from alleles alone, either way round
  expect_identical(allele_match("A", "A", "T"), "ambiguous")
  expect_identical(allele_match("T", "A", "T"), "ambiguous")
  expect_identical(allele_match("C", "C", "G"), "ambiguous")
  expect_identical(allele_match("G", "G", "C"), "ambiguous")
})
