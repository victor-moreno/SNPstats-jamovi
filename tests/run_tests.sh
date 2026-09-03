#!/usr/bin/env bash
# Run the SNPstats test suite against the active R's default library (managed
# by rig). Requires: bash tests/setup_test_env.sh  (run once, or after adding
# a new dependency).
set -euo pipefail

# A bare jmvtools::prepare()/install() regenerates R/snpPGS.h.R and drops the
# caseLevel default, which breaks every snpPGS() call in the suite with an
# unrelated-looking error. Re-apply before installing.
bash tools/patch_h.sh R/snpPGS.h.R || [ $? -eq 10 ]

# Reinstall the package so source changes are picked up, then run testthat.
R CMD INSTALL --no-byte-compile . >/dev/null

Rscript tests/run_tests.R
