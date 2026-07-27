# Verify the proposed jmvcore::State fix end to end:
#   1. unpatched  -> analysis fails to construct
#   2. patched    -> constructs, and state survives the save/restore cycle
suppressMessages({ library(SNPstats); library(RProtoBuf) })
jmvcore:::initProtoBuf()

d <- read.delim("data/CRCgenet-SNPs.tsv", stringsAsFactors = TRUE)
snps <- c("rs12080929", "rs10911251", "rs10936599", "rs6691170")

.optPB <- function(v) {
  pb <- RProtoBuf::new(jamovi.coms.AnalysisOption)
  if (is.logical(v) && length(v) == 1L)        pb$o <- if (isTRUE(v)) 1L else 0L
  else if (is.character(v) && length(v) == 1L) pb$s <- v
  else if (is.numeric(v) && length(v) == 1L)   pb$d <- as.numeric(v)
  else if (is.null(v))                         pb$o <- 2L
  else { inner <- RProtoBuf::new(jamovi.coms.AnalysisOptions)
         inner$options <- lapply(v, .optPB); inner$hasNames <- FALSE; pb$c <- inner }
  pb
}
.optionsPB <- function(vals) {
  pb <- RProtoBuf::new(jamovi.coms.AnalysisOptions)
  pb$names <- names(vals); pb$hasNames <- TRUE
  pb$options <- lapply(vals, .optPB); pb
}
.defaults <- local({ o <- snpStatsOptions$new(); as.list(o$values())[o$names] })
mk <- function(over = list()) {
  vals <- .defaults; vals[names(over)] <- over
  o <- snpStatsOptions$new(); o$fromProtoBuf(.optionsPB(vals)); o
}
base <- list(snps = as.list(snps), response = "phenotype", covariates = list("sex"),
             haploFreq = TRUE, haploAssoc = TRUE, haploInteraction = TRUE)

cat("=== 1. UNPATCHED jmvcore", as.character(packageVersion("jmvcore")), "===\n")
r <- tryCatch({ snpStatsClass$new(options = mk(base), data = d,
                                  analysisId = 1, revision = 1); "constructed" },
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  analysis with a `type: State` element ->", r, "\n\n")

cat("=== 2. APPLYING THE PROPOSED FIX ===\n")
ns <- asNamespace("jmvcore")
St <- get("State", envir = ns)
St$set("public", "initialize",
       function(options, name = "", title = "", visible = FALSE,
                clearWith = "*", refs = character()) {
         super$initialize(options = options, name = name, title = "",
                          visible = FALSE, clearWith = clearWith, refs = refs)
       },
       overwrite = TRUE)
cat("  State$initialize replaced\n\n")

cat("=== 3. PATCHED ===\n")
statefile <- tempfile(fileext = ".pb")
a0 <- tryCatch(snpStatsClass$new(options = mk(base), data = d,
                                 analysisId = 1, revision = 1),
               error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(a0)) {
  cat("  analysis constructed OK\n")
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run()
  sx <- a0$results$ldHaploGroup$haploGroup$haploInterState
  payload <- list(rows = a0$results$ldHaploGroup$haploGroup$haploInteractionTable$asDF)
  sx$setState(payload)
  a0$.save()
  cat("  state set (", nrow(payload$rows), "rows ) and saved\n")

  # unrelated click: ciWidth is not in haploInterState's clearWith
  a1 <- snpStatsClass$new(options = mk(c(base, list(ciWidth = 90))),
                          data = d, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  s1 <- a1$results$ldHaploGroup$haploGroup$haploInterState$state
  cat("  after unrelated click + restore: state survived =", !is.null(s1),
      " rows =", if (is.null(s1)) NA else nrow(s1$rows), "\n")

  # a change that IS in clearWith must invalidate it
  a2 <- snpStatsClass$new(options = mk(c(base, list(haploEffect = "dominant"))),
                          data = d, analysisId = 1, revision = 3)
  a2$.setStatePathSource(function() statefile)
  a2$init(); a2$.load()
  s2 <- a2$results$ldHaploGroup$haploGroup$haploInterState$state
  cat("  after a clearWith-triggering change: state cleared =", is.null(s2), "\n")
}
