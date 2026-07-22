#!/usr/bin/env bash
# Create a self-contained, project-local R library (.Rlib-arm / .Rlib-x64) to
# build and test SNPstats without touching the user/site R libraries.
#
# Why: in restricted/CI sandboxes the personal R library (~/Library/R/...) may be
# unreadable. Installing every dependency into the project-local library and
# pointing R at it with --vanilla + R_LIBS_USER avoids that path entirely. The
# library is architecture-specific and git-ignored (see docs/ENVIRONMENT.md);
# this picks .Rlib-arm or .Rlib-x64 from `uname -m` like install_jamovi.sh does.
#
# Usage:   bash tests/setup_test_env.sh        # from the package root
# Then:    bash tests/run_tests.sh             # run the suite
set -euo pipefail

CRAN="https://cloud.r-project.org"
ARCH="$(uname -m)"
case "$ARCH" in
  arm64)   PD="$(pwd)/.Rlib-arm" ;;
  x86_64)  PD="$(pwd)/.Rlib-x64" ;;
  *)       echo "Unsupported architecture: $ARCH" >&2; exit 1 ;;
esac
mkdir -p "$PD"

# Direct + recursive dependencies (genetics/haplo.stats pull combinat, gdata,
# gtools, mvtnorm, MASS, arsenal; jmvcore pulls rlang, jsonlite, base64enc).
PKGS='c("R6","jmvcore","nnet","genetics","haplo.stats","testthat")'

R_LIBS_USER="$PD" Rscript --vanilla -e "
  install.packages($PKGS, lib='$PD', repos='$CRAN', dependencies=c('Depends','Imports','LinkingTo'))
  miss <- Filter(function(p) !requireNamespace(p, quietly=TRUE),
                 c('R6','jmvcore','nnet','genetics','haplo.stats','testthat'))
  if (length(miss)) stop('missing after install: ', paste(miss, collapse=', '))
  cat('all dependencies available in', '$PD', '\n')
"

R_LIBS_USER="$PD" Rscript --vanilla -e "
  install.packages('jmvtools', repos='https://repo.jamovi.org')
"

# Install SNPstats itself into the local library. R_ENVIRON_USER / R_PROFILE_USER
# are nulled so the install subprocess does not re-add the blocked user library.
R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER="$PD" \
  R CMD INSTALL --no-byte-compile --library="$PD" .


echo "Environment ready. Run: bash tests/run_tests.sh"
