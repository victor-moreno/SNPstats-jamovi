#!/usr/bin/env bash
# Build and install SNPstats into jamovi desktop and/or a running jamovi
# Docker container. Same shape as the other modules' tools/install.sh.
#
#   bash tools/install.sh              both targets, whichever are available
#   bash tools/install.sh desktop
#   bash tools/install.sh docker [container]     (default container: jamovi)
#
# The desktop target uses whichever R `Rscript` resolves to (respecting
# ~/.Rprofile, which appends jamovi.app's bundled module library) -- R here is
# managed by rig, which already isolates by version and architecture. Pick the
# active version with `rig default <version>` first if needed; it must match the
# R version jamovi.app bundles, or jmvcore segfaults on load.
#
# Replaces the old install_jamovi.sh (desktop) and install_snpstats_docker.sh
# (docker); both are folded in here unchanged in substance.
set -euo pipefail

TARGET="${1:-both}"
CONTAINER="${2:-jamovi}"

# Unlike the other modules, the R package sits at the repo root, not in a
# subdirectory named after itself.
HERE="$(cd "$(dirname "$0")/.." && pwd)"
MODULE=SNPstats
VERSION="$(awk -F': *' '$1 == "Version" { print $2; exit }' "$HERE/DESCRIPTION")"
ARTIFACT="$HERE/${MODULE}_${VERSION}.jmo"

# ── desktop ──────────────────────────────────────────────────────────────────
install_desktop() {
  local APP APP_R APP_BASE MODDIR LOG APP_MAJOR JMVT_MAJOR
  APP=/Applications/jamovi.app
  APP_R="$APP/Contents/Frameworks/R.framework/Versions/Current/Resources/bin/R"
  APP_BASE="$APP/Contents/Resources/modules/base/R"
  MODDIR="$HOME/Library/Application Support/jamovi/modules/$MODULE"
  [ -x "$APP_R" ] || { echo "error: no R inside $APP" >&2; return 1; }

  echo ">> desktop: building $MODULE $VERSION with $(Rscript -e 'cat(R.version.string)')"
  cd "$HERE"

  # jmvtools bundles the jamovi-compiler, which hard-refuses a jamovi whose
  # major version it predates ("a newer version of the jamovi-compiler (or
  # jmvtools) is required", installer.js). jamovi 2.7 -> 28 tripped exactly
  # that. Majors track each other, so mirror jamovi's when they diverge.
  APP_MAJOR="$(cut -d. -f1 "$APP/Contents/Resources/jamovi/version")"
  JMVT_MAJOR="$(Rscript -e 'cat(strsplit(as.character(packageVersion("jmvtools")), "[.]")[[1]][1])' 2>/dev/null || echo 0)"
  if [ "$JMVT_MAJOR" != "$APP_MAJOR" ]; then
    echo "   jmvtools major $JMVT_MAJOR != jamovi major $APP_MAJOR — updating jmvtools"
    Rscript -e "install.packages('jmvtools', repos='https://repo.jamovi.org')"
  fi

  LOG="$(mktemp)"
  Rscript -e 'jmvtools::install()' 2>&1 | tee "$LOG" | grep -vE '^\s*$' || true

  # jmvtools::install() can report errors on stdout while exiting successfully.
  # It can also claim installation succeeded after a SingletonLock failure.
  [ -f "$ARTIFACT" ] || {
    echo "error: jmvtools did not produce $ARTIFACT" >&2
    rm -f "$LOG"; return 1
  }
  # jmvtools::install() has already regenerated R/snpPGS.h.R from snpPGS.a.yaml
  # and wiped the caseLevel default. Put it back FIRST: the compiler ran no
  # matter what happened to jamovi.app afterwards, so every exit path below --
  # including the SingletonLock one -- would otherwise leave the working tree
  # holding a header that breaks every scripted snpPGS() call. See NEWS.md.
  bash "$HERE/tools/patch_h.sh" R/snpPGS.h.R || [ $? -eq 10 ]
  R CMD INSTALL --no-byte-compile . >/dev/null
  echo ">> desktop: reinstalled locally (patched)"

  if grep -q 'SingletonLock' "$LOG"; then
    echo
    echo "!! jamovi.app could not be driven (SingletonLock denied) -- quit jamovi"
    echo "!! and re-run, or install the .jmo that was still built by hand:"
    echo "!!   jamovi -> Modules -> Install from file -> $ARTIFACT"
    rm -f "$LOG"
    return 0
  fi
  if ! grep -q 'Module installed successfully' "$LOG"; then
    echo "error: jmvtools::install() did not install the module (see above)" >&2
    rm -f "$LOG"; return 1
  fi
  rm -f "$LOG"

  # jamovi.app writes the module directory asynchronously, so it can still be
  # missing for a second or two after jmvtools has reported success -- don't
  # mistake that for a failed install.
  local i
  for i in $(seq 1 20); do
    [ -d "$MODDIR" ] && break
    sleep 1
  done
  if [ ! -d "$MODDIR" ]; then
    echo "!! desktop: install reported success but $MODDIR never appeared."
    echo "!! Install $ARTIFACT by hand (Modules -> Install from file)."
    return 1
  fi
  echo ">> desktop: installed at $MODDIR"

  # The copy jamovi loads was built from the unpatched header, so scripted
  # snpPGS() calls (tests, Rj) still die on the missing caseLevel default.
  # Rebuild it from the patched source with jamovi's own R.
  R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=/dev/null \
    R_LIBS="$MODDIR/R:$APP_BASE" \
    "$APP_R" CMD INSTALL --no-byte-compile --library="$MODDIR/R" . >/dev/null
  echo ">> desktop: reinstalled into $MODDIR/R (patched)"
}

# ── docker ───────────────────────────────────────────────────────────────────
# jamovi's own route (jamovi >= 28.4): the jamovi compiler on this machine
# builds the module in a throwaway container off the running container's image
# (haplo.stats included, compiled against the image's R) and hands the .jmo to
# the server over its stdin, which installs it into $HOME/.jamovi/modules in the
# container. Nothing is added to the image. Needs here: docker, node, and a
# jamovi-src checkout at the image's version (tools/jamovi-src, a sibling
# ../jamovi-src, or JMC_COMPILER=/path/to/jamovi-compiler). The container must
# run upstream's docker-compose.yaml (stdin_open and --stdin-slave), and the
# Docker VM must see this repo read-write. To keep the module across
# 'down'/'up', mount a volume at /root/.jamovi (jamovi-skill's
# templates/docker/docker-compose.override.yaml).
find_compiler() {
  local c
  for c in "${JMC_COMPILER:-}" "$HERE/tools/jamovi-src/jamovi-compiler" "$HERE/../jamovi-src/jamovi-compiler"; do
    [ -n "$c" ] && [ -f "$c/index.js" ] && [ -f "$c/docker.js" ] && { (cd "$c" && pwd); return 0; }
  done
  return 1
}

install_docker() {
  local JMC CHOME MODR OLD NEW i
  if ! docker ps --format '{{.Names}}' 2>/dev/null | grep -qx "$CONTAINER"; then
    echo "!! docker: container '$CONTAINER' is not running — skipping"
    return 0
  fi
  JMC="$(find_compiler)" || {
    echo "!! docker: no jamovi-compiler with docker support (jamovi >= 28.4) found." >&2
    echo "!! Link a jamovi-src checkout at tools/jamovi-src, or set JMC_COMPILER." >&2
    return 1; }
  command -v node >/dev/null || { echo "!! docker: node is needed on this machine" >&2; return 1; }
  if [ ! -d "$JMC/node_modules" ]; then
    echo ">> docker: installing the compiler's npm dependencies (once)"
    ( cd "$JMC" && npm install --no-audit --no-fund >/dev/null ) || return 1
  fi
  if [ "$(docker inspect -f '{{.Config.OpenStdin}}' "$CONTAINER")" != true ]; then
    echo "!! docker: '$CONTAINER' was not started with stdin open (--stdin-slave);" >&2
    echo "!! start it with jamovi-src's docker-compose.yaml" >&2
    return 1
  fi

  CHOME="$(docker exec "$CONTAINER" sh -c 'echo $HOME')"
  MODR="$CHOME/.jamovi/modules/$MODULE/R"
  OLD="$(docker exec "$CONTAINER" sh -c "grep '^build-time' '$CHOME/.jamovi/modules/$MODULE/jamovi.yaml' 2>/dev/null" || true)"

  echo ">> docker: building $MODULE $VERSION in $(docker inspect -f '{{.Config.Image}}' "$CONTAINER")"
  node "$JMC/index.js" --install "$HERE" --home "docker:$CONTAINER" || return 1
  rm -f "$HERE/.jmc-docker.jmo"   # the build container's artifact, already handed over

  # the server installs it asynchronously after reading its stdin
  NEW=""
  for i in $(seq 1 60); do
    NEW="$(docker exec "$CONTAINER" sh -c "grep '^build-time' '$CHOME/.jamovi/modules/$MODULE/jamovi.yaml' 2>/dev/null" || true)"
    [ -n "$NEW" ] && [ "$NEW" != "$OLD" ] && break
    sleep 1
  done
  if [ -z "$NEW" ] || [ "$NEW" = "$OLD" ]; then
    echo "!! docker: jamovi has not installed the new $MODULE after 60 s (docker logs $CONTAINER)" >&2
    return 1
  fi

  # jmc regenerated R/snpPGS.h.R in this tree and dropped the caseLevel default,
  # so the installed package has a snpPGS() no scripted call can use (the jamovi
  # UI is unaffected). Re-apply the patch and rebuild the R package over the
  # installed copy with the container's R; the yaml/ui jamovi installed stay.
  bash "$HERE/tools/patch_h.sh" "$HERE/R/snpPGS.h.R" || [ $? -eq 10 ]
  tar --no-mac-metadata --no-xattrs -C "$HERE" -cf - DESCRIPTION NAMESPACE NEWS.md R data \
    | docker exec -i "$CONTAINER" sh -c \
        "rm -rf /tmp/${MODULE}-src && mkdir -p /tmp/${MODULE}-src && tar -C /tmp/${MODULE}-src -xf -"
  docker exec "$CONTAINER" bash -c "source /usr/lib/jamovi/bin/env.conf 2>/dev/null || true; \
    R_LIBS='$MODR:/usr/lib/jamovi/modules/base/R' R CMD INSTALL --no-byte-compile \
      --library='$MODR' /tmp/${MODULE}-src >/dev/null" || return 1
  echo ">> docker: installed at $CHOME/.jamovi/modules/$MODULE (patched)"

  if ! docker inspect -f '{{range .Mounts}}{{println .Destination}}{{end}}' "$CONTAINER" | grep -qx "$CHOME/.jamovi"; then
    echo "!! docker: $CHOME/.jamovi is not a volume: the module is lost when the container is"
    echo "!! recreated (down/up). See jamovi-skill's templates/docker/docker-compose.override.yaml"
  fi
  echo ">> docker: installed $MODULE. Open jamovi in the browser (Analyses menu)."
}

case "$TARGET" in
  desktop) install_desktop ;;
  docker)  install_docker ;;
  both)    install_desktop || true; echo; install_docker || true ;;
  *)       echo "usage: install.sh [desktop|docker|both] [container]" >&2; exit 1 ;;
esac
