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
# Runtime (ephemeral) install into a running container. Idempotent: Node, jmc
# and the R deps are bootstrapped only once (they survive `docker restart` but
# not `compose down/up`); every run recompiles the module and reloads jamovi.
install_docker() {
  local JSRC NODE_VER
  if ! docker ps --format '{{.Names}}' 2>/dev/null | grep -qx "$CONTAINER"; then
    echo "!! docker: container '$CONTAINER' is not running — skipping"
    return 0
  fi

  # jamovi-src is symlinked into tools/; older checkouts had it at the root.
  # Only needed to bootstrap jmc and to read the pinned Node version.
  if   [ -d "$HERE/tools/jamovi-src" ]; then JSRC="$HERE/tools/jamovi-src"
  elif [ -d "$HERE/jamovi-src" ];       then JSRC="$HERE/jamovi-src"
  else JSRC=""
  fi

  NODE_VER=""
  if [ -n "$JSRC" ]; then
    # Track whatever the Dockerfile pins rather than restating it here — they
    # drifted apart once already (v22 hardcoded while the image had moved).
    NODE_VER="$(grep -oE 'nodejs\.org/dist/v[0-9]+\.[0-9]+\.[0-9]+' "$JSRC/docker/jamovi-Dockerfile" \
                | head -1 | grep -oE 'v[0-9.]+' || true)"
  fi
  if ! docker exec "$CONTAINER" sh -c 'command -v jmc >/dev/null 2>&1'; then
    [ -n "$JSRC" ] && [ -n "$NODE_VER" ] || {
      echo "!! docker: jmc is not in the container and there is no jamovi-src" >&2
      echo "!! checkout to bootstrap it from. Install the compiler in the image." >&2
      return 1
    }
  fi

  echo ">> docker: copying source into $CONTAINER"
  # --no-mac-metadata/--no-xattrs: AppleDouble ._ files otherwise land in the
  # container and jmc tries to compile them.
  tar --no-mac-metadata --no-xattrs -C "$HERE" -cf - DESCRIPTION NAMESPACE NEWS.md R jamovi data \
    | docker exec -i "$CONTAINER" sh -c \
        "rm -rf /tmp/${MODULE}-src && mkdir -p /tmp/${MODULE}-src && tar -C /tmp/${MODULE}-src -xf -"
  tar --no-mac-metadata --no-xattrs -C "$HERE/tools" -cf - patch_h.sh \
    | docker exec -i "$CONTAINER" tar -C /tmp -xf -

  if docker exec "$CONTAINER" sh -c 'command -v jmc >/dev/null 2>&1'; then
    echo ">> docker: jmc already present (baked image)"
  else
    echo ">> docker: copying compiler source for the jmc bootstrap"
    tar --no-mac-metadata --no-xattrs -C "$JSRC" -cf - jamovi-compiler \
      | docker exec -i "$CONTAINER" sh -c \
          'rm -rf /tmp/jamovi-compiler && tar -C /tmp -xf -'
  fi

  echo ">> docker: bootstrap toolchain (once) + jmc --install"
  docker exec -i "$CONTAINER" bash -s "$MODULE" "$NODE_VER" <<'INCONTAINER'
set -euo pipefail
MODULE="$1"
NODE_VER="${2:-}"
source /usr/lib/jamovi/bin/env.conf 2>/dev/null || true
# ask R where it lives rather than pinning a version that moves with the image
RHOME="${R_HOME:-$(R RHOME 2>/dev/null || true)}"
[ -n "$RHOME" ] || { echo "   error: no R in the container" >&2; exit 1; }
RLIBS=/usr/lib/jamovi/modules/base/R

if ! command -v node >/dev/null 2>&1; then
  [ -n "$NODE_VER" ] || { echo "   error: node missing and no version to fetch" >&2; exit 1; }
  echo "   installing Node ${NODE_VER}"
  case "$(uname -m)" in
    aarch64|arm64) NODEARCH=arm64 ;;
    *)             NODEARCH=x64   ;;
  esac
  curl -L -f -o /tmp/node.tar.gz \
    "https://nodejs.org/dist/${NODE_VER}/node-${NODE_VER}-linux-${NODEARCH}.tar.gz"
  mkdir -p /opt/node && tar -xzf /tmp/node.tar.gz -C /opt/node --strip-components=1
  ln -sf /opt/node/bin/node /usr/local/bin/node
  ln -sf /opt/node/bin/npm  /usr/local/bin/npm
fi
export PATH=/opt/node/bin:$PATH

if ! command -v jmc >/dev/null 2>&1; then
  echo "   installing jamovi-compiler (jmc)"
  ( cd /tmp/jamovi-compiler && npm install --no-audit --no-fund && npm install -g )
  # only symlink our own /opt/node jmc; never clobber a jmc baked at /usr/local
  [ -e /opt/node/bin/jmc ] && ln -sf /opt/node/bin/jmc /usr/local/bin/jmc
fi

# haplo.stats is compiled and is not in the base image; the rest of SNPstats'
# dependencies already resolve from RLIBS.
if [ ! -d "$RLIBS/haplo.stats" ]; then
  echo "   installing R dep haplo.stats (this compiles, ~minutes)"
  "$RHOME/bin/R" --vanilla -q -e \
    "install.packages('haplo.stats', lib='$RLIBS', repos='https://cloud.r-project.org')"
fi

echo "   jmc --install"
jmc --install "/tmp/${MODULE}-src" \
    --to /usr/lib/jamovi/modules \
    --rhome "$RHOME" \
    --rlibs "$RLIBS" \
    --patch-version --skip-deps

[ -f "/usr/lib/jamovi/modules/${MODULE}/jamovi.yaml" ] || {
  echo "   error: jmc did not install ${MODULE}" >&2; exit 1; }

# jmc regenerates R/snpPGS.h.R in place and drops the caseLevel default, so the
# module it just installed has a snpPGS() no scripted call can use (the jamovi
# UI is unaffected — it goes through the options class). Re-apply the patch and
# rebuild the R package over jmc's copy; the yaml/ui it wrote stay.
bash /tmp/patch_h.sh "/tmp/${MODULE}-src/R/snpPGS.h.R" || [ $? -eq 10 ]
echo "   re-installing patched R package"
R_LIBS="$RLIBS" "$RHOME/bin/R" CMD INSTALL --no-byte-compile \
    --library="/usr/lib/jamovi/modules/${MODULE}/R" "/tmp/${MODULE}-src" >/dev/null
INCONTAINER

  echo ">> docker: restarting $CONTAINER to load the module"
  docker restart "$CONTAINER" >/dev/null
  echo ">> docker: installed $MODULE. Open http://127.0.0.1:41337 (Analyses menu)."
}

case "$TARGET" in
  desktop) install_desktop ;;
  docker)  install_docker ;;
  both)    install_desktop || true; echo; install_docker || true ;;
  *)       echo "usage: install.sh [desktop|docker|both] [container]" >&2; exit 1 ;;
esac
