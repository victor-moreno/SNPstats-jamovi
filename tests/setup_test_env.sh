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

# Direct + recursive dependencies (haplo.stats pulls arsenal, MASS; jmvcore
# pulls rlang, jsonlite). ggplot2 and base64enc are declared Imports (plot
# rendering / embedded weights decoding), so they must be present for
# R CMD INSTALL to succeed. RProtoBuf is test-only: jmvcore's save/load is
# protobuf-based and the refresh suites replay jamovi's option-click cycle
# through it. Without RProtoBuf both refresh files skip entirely — the restore
# path then goes completely untested, which is how it stayed untested for a
# while. The `genetics` package is deliberately absent: it is neither a
# dependency nor an oracle any more (see helper-data.R).

# RProtoBuf goes in FIRST, on its own. jmvcore/R/protobuf.R resolves its three
# RProtoBuf wrappers at jmvcore's own install time:
#
#   RProtoBuf_serialize <- if (requireNamespace('RProtoBuf', quietly=TRUE)) RProtoBuf::serialize
#
# RProtoBuf is only a Suggests of jmvcore, so listing them in one vector installs
# jmvcore first and binds all three to NULL. Analysis$.save()/.load() then fail
# with "could not find function RProtoBuf_serialize" inside their own try(),
# which swallows it — the state file is never written and both refresh suites
# fail against tables that are entirely NA. This is the same defect the CI
# workflow works around; here the ordering is enough.
R_LIBS_USER="$PD" Rscript --vanilla -e "
  install.packages('RProtoBuf', lib='$PD', repos='$CRAN', dependencies=c('Depends','Imports','LinkingTo'))
"

PKGS='c("R6","jmvcore","nnet","haplo.stats","ggplot2","base64enc","testthat")'

R_LIBS_USER="$PD" Rscript --vanilla -e "
  install.packages($PKGS, lib='$PD', repos='$CRAN', dependencies=c('Depends','Imports','LinkingTo'))
  miss <- Filter(function(p) !requireNamespace(p, quietly=TRUE),
                 c('R6','jmvcore','nnet','haplo.stats','ggplot2','base64enc','testthat','RProtoBuf'))
  if (length(miss)) stop('missing after install: ', paste(miss, collapse=', '))
  # the wrappers above are the whole point of the ordering — prove they resolved
  stopifnot(is.function(jmvcore:::RProtoBuf_serialize),
            is.function(jmvcore:::RProtoBuf_read),
            is.function(jmvcore:::RProtoBuf_new))
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
