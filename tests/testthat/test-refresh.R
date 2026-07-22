# Refresh / recomputation regression tests.
#
# These replay jamovi's real option-click cycle, which the rest of the suite
# never exercises: jamovi rebuilds the analysis object on every click, restores
# the previous results from a protobuf state file, then runs. Only that cycle
# reveals whether an option click recomputes a table it cannot affect.
#
# Requires RProtoBuf (jmvcore's save/load is protobuf-based) — skipped without it.

skip_if_not_installed("RProtoBuf")

jmvcore:::initProtoBuf()

# ── Build an AnalysisOptions protobuf ────────────────────────────────────────
# snpStatsOptions$new() leaves the options protobuf empty, so jmvcore computes
# no option changes and clearWith can never fire. Options must arrive the way
# the engine sends them — as a protobuf — for this cycle to be realistic.
# Inverse of jmvcore:::parseOptionPB.
.optPB <- function(v) {
  pb <- RProtoBuf::new(jamovi.coms.AnalysisOption)
  if (is.logical(v) && length(v) == 1L)        pb$o <- if (isTRUE(v)) 1L else 0L
  else if (is.character(v) && length(v) == 1L) pb$s <- v
  else if (is.numeric(v) && length(v) == 1L)   pb$d <- as.numeric(v)
  else if (is.null(v))                         pb$o <- 2L
  else {
    inner <- RProtoBuf::new(jamovi.coms.AnalysisOptions)
    inner$options  <- lapply(v, .optPB)
    inner$hasNames <- FALSE
    pb$c <- inner
  }
  pb
}

.optionsPB <- function(vals) {
  pb <- RProtoBuf::new(jamovi.coms.AnalysisOptions)
  pb$names    <- names(vals)
  pb$hasNames <- TRUE
  pb$options  <- lapply(vals, .optPB)
  pb
}

.base_opts <- list(response = .resp, snps = as.list(.snps2), covariates = list("sex"),
                   covDesc = TRUE, snpSummary = TRUE, showAIC = FALSE, hweTest = FALSE)

# Every option, at its declared default. The jamovi engine always sends the FULL
# option set on every run, so the saved state records all of them and compProtoBuf
# reports only genuinely-changed options. A sparse option set (only a handful set)
# makes compProtoBuf treat every unsent option as "added/changed" — clearing every
# table on restore — which is not how jamovi behaves and masks the real refresh
# logic. Merge overrides onto the full defaults to reproduce the engine faithfully.
.all_defaults <- local({ o <- snpStatsOptions$new(); as.list(o$values())[o$names] })

.mk_opts <- function(over = list()) {
  vals <- modifyList(.all_defaults, .base_opts)
  vals[names(over)] <- over
  o <- snpStatsOptions$new()
  o$fromProtoBuf(.optionsPB(vals))
  o
}

.covdesc_of <- function(a) a$results$descGroup$covDescGroup$covDescTable

# Run one click: fresh object + restore + run, counting real computations.
# Returns the table's contents and how many times compute_cov_desc ran.
.click <- function(statefile, over) {
  n <- 0L
  bump <- function() n <<- n + 1L
  # The tracer runs in the traced function's frame, so a local name would not
  # resolve; bquote inlines the counting function itself.
  suppressMessages(trace("compute_cov_desc", tracer = bquote(.(bump)()),
                         print = FALSE, where = asNamespace("SNPstats")))
  on.exit(suppressMessages(untrace("compute_cov_desc", where = asNamespace("SNPstats"))))

  a <- snpStatsClass$new(options = .mk_opts(over), data = .test_data,
                         analysisId = 1, revision = 2)
  a$.setStatePathSource(function() statefile)
  a$init(); a$.load(); a$run()
  list(computes = n, df = .covdesc_of(a)$asDF, rows = .covdesc_of(a)$rowCount)
}

.first_run <- function(statefile) {
  a <- snpStatsClass$new(options = .mk_opts(), data = .test_data,
                         analysisId = 1, revision = 1)
  a$.setStatePathSource(function() statefile)
  a$init(); a$run(); a$.save()
  .covdesc_of(a)$asDF
}

test_that("covDescTable is not recomputed when an unrelated option is clicked", {
  statefile <- tempfile(fileext = ".pb")
  base <- .first_run(statefile)

  # showAIC and hweTest are absent from covDescTable's clearWith, so its state
  # survives the restore and the cached descriptives are reused.
  for (opt in list(list(showAIC = TRUE), list(hweTest = TRUE))) {
    got <- .click(statefile, opt)
    expect_equal(got$computes, 0L, label = paste(names(opt), "computes"))
    expect_equal(got$df, base, label = paste(names(opt), "contents"))
  }
})

test_that("covDescTable IS recomputed when an option it depends on changes", {
  statefile <- tempfile(fileext = ".pb")
  .first_run(statefile)

  # Both are in covDescTable's clearWith -> state is cleared -> must recompute.
  got <- .click(statefile, list(covariates = list("sex", "age")))
  expect_equal(got$computes, 1L)
  expect_gt(got$rows, 3L)

  got <- .click(statefile, list(response = "bmi"))
  expect_equal(got$computes, 1L)
})

test_that("addRow-built rows never survive jamovi's restore (why rowCount cannot gate)", {
  # Pins the constraint the whole fix rests on: Table$fromProtoBuf only restores
  # cells into rows that exist by the end of .init(). covDescTable's rows are
  # pre-created there, so its contents come back BEFORE .run() — which is what
  # stops the table blanking on every click. Were the rows added in .run()
  # instead, nothing would be restored and rowCount here would be 0.
  statefile <- tempfile(fileext = ".pb")
  base <- .first_run(statefile)

  a <- snpStatsClass$new(options = .mk_opts(), data = .test_data,
                         analysisId = 1, revision = 2)
  a$.setStatePathSource(function() statefile)
  a$init()
  expect_equal(.covdesc_of(a)$rowCount, nrow(base))   # rows exist, still empty
  expect_true(all(is.na(.covdesc_of(a)$asDF$stat_overall)))

  a$.load()
  expect_equal(.covdesc_of(a)$asDF, base)             # restored, before run()

  # snpSummaryTable gets the same property from its rows:(snps) binding.
  expect_equal(a$results$descGroup$snpSummaryTablesGroup$snpSummaryTable$rowCount,
               length(.snps2))
})

# No column of CRCgenet-SNPs.tsv contains NA, so missingness is injected here.
.data_with_na <- local({
  d <- .test_data
  d$age[1:20] <- NA
  d
})

test_that("covDescTable emits a Missing row only for variables that have missing", {
  # sex has no missing values, age (injected) does. A Missing row must appear
  # for the latter only — .init reads the data, so it predicts exactly this.
  expect_false(anyNA(.data_with_na$sex))
  expect_true(anyNA(.data_with_na$age))

  res <- run_snp(data = .data_with_na, snps = .snps2, response = .resp,
                 covariates = c("sex", "age"), covDesc = TRUE)
  df  <- as_df(res$descGroup$covDescGroup$covDescTable)
  expect_equal(nrow(df[df$variable == "sex" & df$level == "Missing", ]), 0L)
  expect_equal(nrow(df[df$variable == "age" & df$level == "Missing", ]), 1L)
})

test_that("a Missing row survives completeCases, which masks its values away", {
  # The row's presence is decided from the raw column, not the masked one, so
  # the structure stays predictable at .init. Masking empties the count rather
  # than removing the row, which .init could not have foreseen.
  res <- run_snp(data = .data_with_na, snps = .snps2, response = .resp,
                 covariates = "age", covDesc = TRUE, completeCases = TRUE)
  df  <- as_df(res$descGroup$covDescGroup$covDescTable)
  miss <- df[df$variable == "age" & df$level == "Missing", ]
  expect_equal(nrow(miss), 1L)
  expect_match(miss$stat_overall, "^0 ")
})

test_that("per-SNP descriptive tables survive an unrelated click", {
  statefile <- tempfile(fileext = ".pb")
  desc <- list(allFreq = TRUE, genoFreq = TRUE, hweTest = TRUE)
  tbls <- c("allFreqTable", "genoFreqTable", "hweTable")
  item_of <- function(a) a$results$descGroup$descSnpResults$get(key = .snps2[[1]])

  a1 <- snpStatsClass$new(options = .mk_opts(desc), data = .test_data,
                          analysisId = 1, revision = 1)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$run(); a1$.save()
  base <- lapply(tbls, function(t) item_of(a1)[[t]]$asDF)
  expect_true(all(vapply(base, nrow, 0L) > 0L))

  a2 <- snpStatsClass$new(options = .mk_opts(c(desc, list(showAIC = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a2$.setStatePathSource(function() statefile)
  a2$init(); a2$.load()
  for (i in seq_along(tbls))
    expect_equal(item_of(a2)[[tbls[i]]]$asDF, base[[i]], label = tbls[i])
})

test_that("Array items are never duplicated by .init and .run both adding them", {
  # Array$addItem() appends unconditionally. A duplicate key makes restore
  # index the LAST item of that name — the empty one — so the SNP comes back
  # blank. Pins that each key appears exactly once after a full run.
  res <- run_snp(data = .test_data, snps = .snps2, response = .resp,
                 covariates = "sex", allFreq = TRUE, genoFreq = TRUE,
                 hweTest = TRUE, snpAssoc = TRUE, snpInteraction = TRUE,
                 modelCodominant = TRUE)
  for (arr in list(res$descGroup$descSnpResults, res$assocGroup$assocSnpResults)) {
    keys <- unlist(arr$itemKeys)
    expect_equal(anyDuplicated(keys), 0L)
    expect_setequal(keys, .snps2)
  }
})

test_that("assocTable rows survive an unrelated click", {
  statefile <- tempfile(fileext = ".pb")
  assoc <- list(snpAssoc = TRUE, modelCodominant = TRUE, modelLogAdditive = TRUE)
  at <- function(a) a$results$assocGroup$assocSnpResults$get(key = .snps2[[1]])$assocTable

  a1 <- snpStatsClass$new(options = .mk_opts(assoc), data = .test_data,
                          analysisId = 1, revision = 1)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$run(); a1$.save()
  base <- at(a1)$asDF
  expect_gt(nrow(base), 0L)

  # showAIC is display-only: absent from assocTable's clearWith, it just
  # reveals the AIC/BIC columns.
  a2 <- snpStatsClass$new(options = .mk_opts(c(assoc, list(showAIC = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a2$.setStatePathSource(function() statefile)
  a2$init(); a2$.load()
  got <- at(a2)$asDF
  expect_equal(nrow(got), nrow(base))
  # showAIC reveals the AIC/BIC columns, so compare the ones common to both.
  common <- intersect(names(base), names(got))
  expect_equal(got[common], base[common])
})

test_that("stratified snpSummaryTable rows and covDesc group columns survive a click", {
  statefile <- tempfile(fileext = ".pb")
  strat <- list(subpop = TRUE)

  a1 <- snpStatsClass$new(options = .mk_opts(strat), data = .test_data,
                          analysisId = 1, revision = 1)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$run(); a1$.save()
  base_sum <- a1$results$descGroup$snpSummaryTablesGroup$snpSummaryTable$asDF
  base_cov <- .covdesc_of(a1)$asDF
  expect_true("stat_g0" %in% names(base_cov))   # group columns exist

  over <- c(strat, list(showAIC = TRUE))
  a2 <- snpStatsClass$new(options = .mk_opts(over), data = .test_data,
                          analysisId = 1, revision = 2)
  a2$.setStatePathSource(function() statefile)
  a2$init(); a2$.load()
  # Restored before run(): stratified rows AND the dynamically added columns.
  expect_equal(a2$results$descGroup$snpSummaryTablesGroup$snpSummaryTable$asDF, base_sum)
  expect_equal(.covdesc_of(a2)$asDF, base_cov)
})

test_that("codominant heterozygous rows show non-zero counts in stratified tables", {
  res  <- run_snp(
    data = .test_data, snps = .snps2[[1]], response = .resp, covariates = "sex",
    snpInteraction = TRUE, modelCodominant = TRUE,
    showStratByCovariate = TRUE, showCrossClassTable = TRUE)
  item <- res$assocGroup$assocSnpResults$get(key = .snps2[[1]])
  het_n <- function(tbl) {
    df  <- as_df(tbl)
    het <- grepl("^[A-Z]/[A-Z]$", df$grp2) &
           vapply(strsplit(df$grp2, "/"), function(p) p[[1]] != p[[2]], logical(1))
    as.integer(sub(" .*", "", df$stat0[het]))
  }
  expect_true(all(het_n(item$stratByCovariate) > 0))
  expect_true(all(het_n(item$crossClassTable) > 0))
})

test_that("interaction tables are neither refit nor rebuilt on an unrelated click", {
  statefile <- tempfile(fileext = ".pb")
  snp1     <- .snps2[[1]]
  int_over <- list(covDesc = FALSE, snpSummary = FALSE, snpInteraction = TRUE,
                   modelCodominant = TRUE, modelDominant = TRUE,
                   showInteractionTable = TRUE, showStratByCovariate = TRUE,
                   showStratByGenotype = TRUE, showCrossClassTable = TRUE)
  four <- function(a) {
    it <- a$results$assocGroup$assocSnpResults$get(key = snp1)
    list(it$interactionTable, it$stratByCovariate, it$stratByGenotype, it$crossClassTable)
  }
  dfs <- function(a) lapply(four(a), as_df)

  a0 <- snpStatsClass$new(options = .mk_opts(int_over), data = .test_data,
                          analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base   <- dfs(a0)
  counts <- vapply(four(a0), function(t) t$rowCount, numeric(1))

  n <- 0L; bump <- function() n <<- n + 1L
  suppressMessages(trace("fit_interaction_model", tracer = bquote(.(bump)()),
                         print = FALSE, where = asNamespace("SNPstats")))
  on.exit(suppressMessages(untrace("fit_interaction_model", where = asNamespace("SNPstats"))))

  a1 <- snpStatsClass$new(options = .mk_opts(c(int_over, list(hweTest = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  restored <- vapply(four(a1), function(t) t$rowCount, numeric(1))
  a1$run()

  expect_equal(n, 0L)                    # hweTest is not in any interaction clearWith → no refit
  expect_equal(restored, counts)         # .init pre-created rows → restore refilled them (no rebuild)
  expect_equal(dfs(a1), base)            # content identical
})

test_that(".int_nrows predicts interaction row counts exactly across common configs", {
  # Pins that .init pre-creates the same number of rows .run produces, so the
  # reuse path (no rebuild) is actually taken. If a predictor drifts from a
  # writer, pre != post here (harmless in production — just a redraw — but a bug
  # to fix). Exercised without restore: init count vs run count on one object.
  snp1  <- .snps2[[1]]
  tabs  <- c("interactionTable", "stratByCovariate", "stratByGenotype", "crossClassTable")
  base_int <- list(covDesc = FALSE, snpSummary = FALSE, snpInteraction = TRUE,
                   showInteractionTable = TRUE, showStratByCovariate = TRUE,
                   showStratByGenotype = TRUE, showCrossClassTable = TRUE,
                   modelCodominant = FALSE)
  cfgs <- list(
    codominant   = list(modelCodominant = TRUE),
    all_models   = list(modelCodominant = TRUE, modelDominant = TRUE, modelRecessive = TRUE,
                        modelOverdominant = TRUE, modelLogAdditive = TRUE),
    logadd_only  = list(modelLogAdditive = TRUE),
    quantitative = list(modelCodominant = TRUE, modelDominant = TRUE,
                        response = "age", responseType = "quantitative"))
  for (nm in names(cfgs)) {
    a <- snpStatsClass$new(options = .mk_opts(modifyList(base_int, cfgs[[nm]])),
                           data = .test_data, analysisId = 1, revision = 1)
    a$init()
    it   <- a$results$assocGroup$assocSnpResults$get(key = snp1)
    pre  <- vapply(tabs, function(t) it[[t]]$rowCount, numeric(1))
    a$run()
    post <- vapply(tabs, function(t) it[[t]]$rowCount, numeric(1))
    expect_equal(pre, post, label = nm)
  }
})

test_that("LD tables are neither recomputed nor rebuilt on an unrelated click", {
  statefile <- tempfile(fileext = ".pb")
  ld_over <- list(snps = as.list(.snps4), covDesc = FALSE, snpSummary = FALSE,
                  ldAnalysis = TRUE, ldMatrix = TRUE)
  ld_of <- function(a) {
    it <- a$results$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")
    list(it$ldTable, it$ldMatrixTable)
  }

  a0 <- snpStatsClass$new(options = .mk_opts(ld_over), data = .test_data,
                          analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base   <- lapply(ld_of(a0), function(t) t$asDF)
  counts <- vapply(ld_of(a0), function(t) t$rowCount, numeric(1))
  expect_true(all(counts > 0))

  n <- 0L; bump <- function() n <<- n + 1L
  suppressMessages(trace("LD", tracer = bquote(.(bump)()), print = FALSE,
                         where = asNamespace("genetics")))
  on.exit(suppressMessages(untrace("LD", where = asNamespace("genetics"))))

  a1 <- snpStatsClass$new(options = .mk_opts(c(ld_over, list(hweTest = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  restored <- vapply(ld_of(a1), function(t) t$rowCount, numeric(1))
  a1$run()

  expect_equal(n, 0L)                    # ldResults clearWith excludes hweTest → no genetics::LD
  expect_equal(restored, counts)         # .init pre-created rows → restore refilled them
  expect_equal(lapply(ld_of(a1), function(t) t$asDF), base)
})

test_that("adding the LD matrix does not rebuild the already-computed pairwise table", {
  statefile <- tempfile(fileext = ".pb")
  a_over <- list(snps = as.list(.snps4), covDesc = FALSE, snpSummary = FALSE,
                 ldAnalysis = TRUE)
  tbl_of <- function(a) a$results$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldTable

  a0 <- snpStatsClass$new(options = .mk_opts(a_over), data = .test_data,
                          analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base_df <- tbl_of(a0)$asDF
  base_n  <- tbl_of(a0)$rowCount
  expect_gt(base_n, 0)

  a1 <- snpStatsClass$new(options = .mk_opts(c(a_over, list(ldMatrix = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  # Toggling ldMatrix must NOT clear the pairwise table: it is restored, filled,
  # before .run — so its per-table .need_fill gate skips and it is not rebuilt.
  expect_false(tbl_of(a1)$isNotFilled())
  expect_equal(tbl_of(a1)$rowCount, base_n)
  a1$run()
  expect_equal(tbl_of(a1)$asDF, base_df)
  mtx <- a1$results$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")$ldMatrixTable
  expect_gt(mtx$rowCount, 0)            # the matrix itself did get built
})

test_that("toggling a haplotype option does not recompute or rebuild the LD tables", {
  # The LD and haplotype analyses share the ldHaploGroup parent but no clearWith
  # option, so enabling a haplotype table must leave the already-computed LD
  # tables untouched. genetics::LD (the LD cost) must not run again.
  statefile <- tempfile(fileext = ".pb")
  ld_over <- list(snps = as.list(.snps4), covDesc = FALSE, snpSummary = FALSE,
                  ldAnalysis = TRUE, ldMatrix = TRUE)
  ld_of <- function(a) {
    it <- a$results$ldHaploGroup$ldGroup$ldResults$get(key = "Overall")
    list(it$ldTable, it$ldMatrixTable)
  }

  a0 <- snpStatsClass$new(options = .mk_opts(ld_over), data = .test_data,
                          analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base   <- lapply(ld_of(a0), function(t) t$asDF)
  counts <- vapply(ld_of(a0), function(t) t$rowCount, numeric(1))
  expect_true(all(counts > 0))

  n <- 0L; bump <- function() n <<- n + 1L
  suppressMessages(trace("LD", tracer = bquote(.(bump)()), print = FALSE,
                         where = asNamespace("genetics")))
  on.exit(suppressMessages(untrace("LD", where = asNamespace("genetics"))))

  a1 <- snpStatsClass$new(options = .mk_opts(c(ld_over, list(haploFreq = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load()
  restored <- vapply(ld_of(a1), function(t) t$rowCount, numeric(1))
  a1$run()

  expect_equal(n, 0L)                    # haploFreq is not in any LD clearWith → no genetics::LD
  expect_equal(restored, counts)         # LD rows restored before .run → not rebuilt
  expect_equal(lapply(ld_of(a1), function(t) t$asDF), base)
})

test_that("haplotype EM/GLM are not refit on an unrelated click", {
  statefile <- tempfile(fileext = ".pb")
  h_over <- list(snps = as.list(.snps4), covDesc = FALSE, snpSummary = FALSE,
                 haploFreq = TRUE, haploAssoc = TRUE, haploInteraction = TRUE)
  hg_of <- function(a) a$results$ldHaploGroup$haploGroup
  dfs <- function(a) {
    hg <- hg_of(a)
    lapply(list(hg$haploFreqTable, hg$haploAssocTable, hg$haploInteractionTable), as_df)
  }

  a0 <- snpStatsClass$new(options = .mk_opts(h_over), data = .test_data,
                          analysisId = 1, revision = 1)
  a0$.setStatePathSource(function() statefile)
  a0$init(); a0$run(); a0$.save()
  base <- dfs(a0)
  expect_true(all(vapply(base, nrow, 0L) > 0L))

  nem <- 0L; nglm <- 0L
  suppressMessages(trace("haplo.em",  tracer = bquote(.(function() nem  <<- nem  + 1L)()),
                         print = FALSE, where = asNamespace("haplo.stats")))
  suppressMessages(trace("haplo.glm", tracer = bquote(.(function() nglm <<- nglm + 1L)()),
                         print = FALSE, where = asNamespace("haplo.stats")))
  on.exit({
    suppressMessages(untrace("haplo.em",  where = asNamespace("haplo.stats")))
    suppressMessages(untrace("haplo.glm", where = asNamespace("haplo.stats")))
  })

  a1 <- snpStatsClass$new(options = .mk_opts(c(h_over, list(hweTest = TRUE))),
                          data = .test_data, analysisId = 1, revision = 2)
  a1$.setStatePathSource(function() statefile)
  a1$init(); a1$.load(); a1$run()

  expect_equal(nem, 0L)                 # clicking hweTest must not re-run haplo.em
  expect_equal(nglm, 0L)                # nor haplo.glm
  expect_equal(dfs(a1), base)           # tables rebuilt from cache, content identical
})
