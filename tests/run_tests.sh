#!/usr/bin/env bash
# Run the SNPstats test suite against the project-local library
# (.Rlib-arm / .Rlib-x64, picked by `uname -m` like install_jamovi.sh).
# Requires: bash tests/setup_test_env.sh  (run once, or after changing R code).
#
# --vanilla + R_LIBS_USER keep R off the (possibly unreadable) user library;
# R_ENVIRON_USER/R_PROFILE_USER=/dev/null stop startup files re-adding it.
set -euo pipefail

ARCH="$(uname -m)"
case "$ARCH" in
  arm64)   PD="$(pwd)/.Rlib-arm" ;;
  x86_64)  PD="$(pwd)/.Rlib-x64" ;;
  *)       echo "Unsupported architecture: $ARCH" >&2; exit 1 ;;
esac
if [ ! -d "$PD/SNPstats" ]; then
  echo "SNPstats not installed in $PD — run: bash tests/setup_test_env.sh" >&2
  exit 1
fi

# Reinstall the package so source changes are picked up, then run testthat.
R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER="$PD" \
  R CMD INSTALL --no-byte-compile --library="$PD" . >/dev/null

R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER="$PD" \
  Rscript --vanilla tests/run_tests.R
