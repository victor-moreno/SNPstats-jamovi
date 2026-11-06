# Check one SNPstats .jmo against a jamovi install, loading it the way jamovi does
# (the module's own R library, then jamovi's base library), then run the test suite.
#
#   Rscript --vanilla check_jmo.R <file.jmo> <env.conf, or "none" for plain R> <suite dir>
#
# Runs in the R that ships inside jamovi, so it can only use what jamovi provides.
args <- commandArgs(TRUE)
jmo <- normalizePath(args[1], winslash = "/")
conf <- args[2]
suite <- normalizePath(args[3], winslash = "/")

report <- function(label, value) {
    line <- sprintf("- **%s**: %s\n", label, value)
    cat(line)
    summary_file <- Sys.getenv("GITHUB_STEP_SUMMARY")
    if (nzchar(summary_file))
        cat(line, file = summary_file, append = TRUE)
}
fail <- function(...) {
    report("FAILED", paste0(...))
    quit(save = "no", status = 1)
}
field <- function(lines, key, sep) {
    hit <- grep(paste0("^", key, sep), lines, value = TRUE)
    if (length(hit) == 0) return(NA_character_)
    trimws(sub(paste0("^", key, sep), "", hit[1]))
}
app_path <- function(value) {
    if (grepl("^(/|[A-Za-z]:)", value)) return(value)
    file.path(dirname(conf), value)
}

report("File", basename(jmo))
extracted <- tempfile("jmo")
dir.create(extracted)
utils::unzip(jmo, exdir = extracted)
module_yaml <- readLines(file.path(extracted, "SNPstats", "jamovi.yaml"), warn = FALSE)
module_stamp <- field(module_yaml, "rVersion", ": *")
report("File stamp", module_stamp)

libs <- file.path(extracted, "SNPstats", "R")
if (conf == "none") {
    report("jamovi", "none available for this platform, so plain R and CRAN packages")
    libs <- c(libs, .libPaths())
} else {
    env <- readLines(conf, warn = FALSE)
    app_stamp <- field(env, "JAMOVI_R_VERSION", "=")
    version_path <- field(env, "JAMOVI_VERSION_PATH", "=")
    report("jamovi", readLines(app_path(version_path), n = 1, warn = FALSE))
    report("jamovi stamp", app_stamp)
    # jamovi's only runtime compatibility check (modules.py, _parse_module_defn)
    if (!identical(module_stamp, app_stamp))
        fail("jamovi would mark this file incompatible: ", module_stamp, " vs ", app_stamp)
    libs <- c(libs, file.path(app_path(field(env, "JAMOVI_MODULES_PATH", "=")), "base", "R"))

    # jamovi 28.4 renaming a panel control is what broke 1.2.0, and no R test can see it
    client <- list.files(app_path(field(env, "JAMOVI_CLIENT_PATH", "=")),
                         pattern = "^analysisui.*[.]js$", recursive = TRUE, full.names = TRUE)
    panel <- readLines(file.path(extracted, "SNPstats", "ui", "snpPGS.js"), warn = FALSE)
    control <- regmatches(panel, regexpr("DefaultControls[.]File[A-Za-z]*", panel))
    control <- sub("DefaultControls[.]", "", unique(control))
    if (length(client) == 0) {
        report("File control", "client bundle not found, not checked")
    } else if (!any(grepl(paste0(control, ":class"), readLines(client[1], warn = FALSE), fixed = TRUE))) {
        fail("jamovi's client has no ", control, " control, so the snpPGS panel cannot be built")
    } else {
        report("File control", paste(control, "present in the client"))
    }
}
.libPaths(libs, include.site = FALSE)

report("R", paste(getRversion(), R.version$platform))
loaded <- tryCatch({ loadNamespace("SNPstats"); TRUE }, error = function(e) conditionMessage(e))
if (!isTRUE(loaded)) fail("module does not load: ", loaded)
report("SNPstats", as.character(utils::packageVersion("SNPstats")))
report("jmvcore", as.character(utils::packageVersion("jmvcore")))
if (!is.loaded("haplo_em_pin", PACKAGE = "haplo.stats"))
    fail("the compiled haplo.stats library did not load")
report("haplo.stats library", dirname(getLoadedDLLs()[["haplo.stats"]][["path"]]))

ns <- asNamespace("SNPstats")
if (utils::packageVersion("SNPstats") > "1.2.1") {
    if (!is.null(formals(get("snpPGS", envir = ns))$caseLevel))
        fail("snpPGS() has no default for caseLevel")
} else {
    # 1.2.1 shipped without that default. It only affects calling snpPGS() from R and
    # would hide everything else here, so correct it in memory for this old release.
    report("caseLevel", "known defect in this release, corrected in memory for the suite")
    corrected <- get("snpPGS", envir = ns)
    formals(corrected)["caseLevel"] <- list(NULL)
    suppressMessages(library(SNPstats))
    for (where in list(ns, as.environment("package:SNPstats"))) {
        unlockBinding("snpPGS", where)
        assign("snpPGS", corrected, envir = where)
        lockBinding("snpPGS", where)
    }
}

suppressMessages(library(SNPstats))
suppressMessages(library(testthat))
setwd(suite)
results <- test_dir("tests/testthat", package = "SNPstats", reporter = "silent",
                    stop_on_failure = FALSE)
for (test in results) {
    for (result in test$results) {
        if (inherits(result, c("expectation_failure", "expectation_error")))
            cat("\n[", test$file, "] ", test$test, "\n  ",
                substr(conditionMessage(result), 1, 500), "\n", sep = "")
    }
}
counts <- as.data.frame(results)
report("Test suite", sprintf("%d pass, %d fail, %d error, %d skipped",
                             sum(counts$passed), sum(counts$failed),
                             sum(counts$error), sum(counts$skipped)))
if (sum(counts$failed) > 0 || any(counts$error))
    quit(save = "no", status = 1)
