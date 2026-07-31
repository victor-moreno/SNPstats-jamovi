# Local test runner.
#
# devtools is not available in this sandbox and the user package library is not
# readable, so dependencies live in a project-local library (.Rlib-arm /
# .Rlib-x64). Run with:
#
#   PD="$(pwd)/.Rlib-arm"   # or .Rlib-x64 on the x86_64 Mac
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER="$PD" \
#     Rscript --vanilla tests/run_tests.R
#
# --vanilla + R_LIBS_USER avoid the (TCC-blocked) ~/Library R user library.

suppressMessages(library(testthat))
suppressMessages(library(SNPstats))

# The summary reporter prints per-file dots and "DONE" but never a totals line —
# only the progress reporter does that, and its per-file format is far noisier.
# So keep this reporter and tally the results object by hand.
#
# stop_on_failure is off for the same reason: it aborts inside test_dir(), which
# would suppress the counts on exactly the runs where they matter most. The
# stop() below restores the non-zero exit that CI and run_tests.sh depend on.
res <- testthat::test_dir("tests/testthat", reporter = "summary",
                          stop_on_failure = FALSE)

df <- as.data.frame(res)
cat(sprintf("\n[ PASS %d | FAIL %d | WARN %d | SKIP %d | ERROR %d ]\n",
            sum(df$passed), sum(df$failed), sum(df$warning),
            sum(df$skipped), sum(df$error)))

if (sum(df$failed) > 0L || any(df$error))
    stop("Test failures.", call. = FALSE)
