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

# The snpImport half of the module lives as much in the browser as in R: format
# dispatch, the file dialog and the Load button are all in jamovi/js, and no R
# test can reach any of them. Run them first, then dump the payloads the panel
# actually produces so test-browser-payloads.R checks the two halves against
# each other instead of the R suite checking its own idea of what JS sends.
if command -v node >/dev/null 2>&1; then
  echo ">> js tests"
  node tests/js/test-groupfiles.js || exit 1
  node tests/js/test-panel.js || exit 1

  mkdir -p .tmp
  node tests/js/dump-payloads.js small .tmp/browser-payloads-small.tsv \
       snp2 snp5 snp9 snp100
else
  echo ">> js tests skipped (no node)"
fi

Rscript tests/run_tests.R
