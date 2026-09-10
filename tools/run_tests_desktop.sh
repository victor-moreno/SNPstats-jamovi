#!/usr/bin/env bash
# Run the SNPstats suite against the copy of the module that jamovi.app loads,
# using jamovi's own bundled R.
#
# Distinct from tests/run_tests.sh, which uses CRAN R and the project-local
# library: jamovi ships its own R and its own dependency set, and those are what
# a user actually runs. A stale module dir here (old haplo.stats, packages left
# behind by a previous release) produces wrong numbers with the suite still
# green on the host, so test the installed copy directly.
#
#   Usage: tools/run_tests_desktop.sh
#          SKIP_INSTALL=1 tools/run_tests_desktop.sh     test what is installed
set -euo pipefail

SRC="$(cd "$(dirname "$0")/.." && pwd)"
APP=/Applications/jamovi.app
APP_R="$APP/Contents/Frameworks/R.framework/Versions/Current/Resources/bin/R"
APP_BASE="$APP/Contents/Resources/modules/base/R"
MODDIR="$HOME/Library/Application Support/jamovi/modules/SNPstats"

[ -x "$APP_R" ] || { echo "error: no R inside $APP" >&2; exit 1; }

cd "$SRC"
if [ "${SKIP_INSTALL:-0}" = 1 ]; then
  echo ">> reusing the module already installed in jamovi.app"
else
  bash tools/install.sh desktop
fi

[ -d "$MODDIR/R/SNPstats" ] \
  || { echo "error: SNPstats not installed at $MODDIR" >&2; exit 1; }

echo ">> running suite (jamovi.app R)"
R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=/dev/null \
R_LIBS="$MODDIR/R:$APP_BASE" \
  "$APP_R" --vanilla -q --no-echo \
    -e 'cat(R.version.string, "| jmvcore", as.character(packageVersion("jmvcore")),
            "| haplo.stats", as.character(packageVersion("haplo.stats")), "\n")' \
    -e 'source("tests/run_tests.R")'
