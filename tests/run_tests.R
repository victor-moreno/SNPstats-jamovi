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

testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
