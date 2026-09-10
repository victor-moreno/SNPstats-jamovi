#!/usr/bin/env bash
# Re-apply the hand-patches that the jamovi UI compiler wipes out.
#
# jmvtools::install() / jmvtools::prepare() regenerate R/snpPGS.h.R from
# snpPGS.a.yaml. A `type: Level` option cannot carry a yaml `default:` (the
# compiler rejects it), so the compiler emits a bare `caseLevel,` formal on the
# public snpPGS() function — and every call that omits caseLevel then dies with
# "argument caseLevel is missing" (tests, R scripts, the documented usage).
#
# jmc --install does the same thing in the container, so this applies there too.
#
# Scripted rather than hand-edited so it survives every rebuild. Idempotent.
# Called by install_jamovi.sh and install_snpstats_docker.sh (after the
# compiler) and by run_tests.sh (which installs whatever is on disk, including a
# tree left behind by a bare jmvtools::prepare()).
set -euo pipefail

HFILE="${1:-R/snpPGS.h.R}"

# `sed -i` takes a mandatory suffix on BSD and rejects one on GNU, and this now
# also runs inside the linux container (install_snpstats_docker.sh) — so write
# through a temp file, which behaves the same everywhere.
if grep -q '^    caseLevel,$' "$HFILE"; then
  sed 's/^    caseLevel,$/    caseLevel = NULL,/' "$HFILE" > "$HFILE.patched"
  mv "$HFILE.patched" "$HFILE"
  echo "patched caseLevel default in $HFILE"
  exit 10          # signals "changed" to callers that need to reinstall
fi
exit 0
