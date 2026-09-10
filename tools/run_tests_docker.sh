#!/usr/bin/env bash
# Run the SNPstats suite inside a running jamovi container, against the module
# exactly as jamovi loads it (container R, container jmvcore, container deps).
#
# This is the check that catches a jamovi/R upgrade shifting results: the host
# suite (tests/run_tests.sh) runs against CRAN R and the project-local library,
# which is a different stack.
#
#   Usage: tools/run_tests_docker.sh [container_name]     (default: jamovi)
#          SKIP_INSTALL=1 tools/run_tests_docker.sh       test what is installed
set -euo pipefail

CONTAINER="${1:-jamovi}"
SRC="$(cd "$(dirname "$0")/.." && pwd)"

docker exec "$CONTAINER" true 2>/dev/null \
  || { echo "error: container '$CONTAINER' is not running" >&2; exit 1; }

if [ "${SKIP_INSTALL:-0}" = 1 ]; then
  echo ">> reusing the module already installed in $CONTAINER"
else
  bash "$SRC/tools/install.sh" docker "$CONTAINER"
  # install.sh restarts the container; wait for it back
  for _ in $(seq 1 30); do
    docker exec "$CONTAINER" true 2>/dev/null && break
    sleep 1
  done
fi

# helper-data.R resolves the fixtures relative to the package root, so tests/
# and data/ have to sit together in one directory.
echo ">> copying suite -> $CONTAINER:/tmp/snpstats-test"
tar --no-mac-metadata --no-xattrs -C "$SRC" -cf - tests data \
  | docker exec -i "$CONTAINER" sh -c \
      'rm -rf /tmp/snpstats-test && mkdir -p /tmp/snpstats-test && tar -C /tmp/snpstats-test -xf -'

echo ">> running suite (container R)"
docker exec "$CONTAINER" bash -c '
set -euo pipefail
cd /tmp/snpstats-test
# R appends its own library to .libPaths(), so only the module dirs are needed
export R_LIBS="/usr/lib/jamovi/modules/SNPstats/R:/usr/lib/jamovi/modules/base/R"
# tests/run_tests.R came across in the tarball; use it rather than a second copy
# of the invocation, so the pass/fail tally is defined in exactly one place
Rscript --vanilla \
  -e "cat(R.version.string, \"| jmvcore\", as.character(packageVersion(\"jmvcore\")),
          \"| haplo.stats\", as.character(packageVersion(\"haplo.stats\")), \"\n\")" \
  -e "source(\"tests/run_tests.R\")"'
