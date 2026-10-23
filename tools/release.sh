#!/usr/bin/env bash
# Build one SNPstats .jmo per OS x CPU, and optionally publish them as a
# GitHub release.
#
#   bash tools/release.sh              # write dist/SNPstats_<version>_<os>-<cpu>.jmo
#   bash tools/release.sh --publish    # ... and create (or replace) release v<version>
#
# Why not the template's two-artifact release (tools/ template, prepare-jmo.sh):
# that repack is only valid for modules without compiled code. SNPstats bundles
# haplo.stats, which has a native library, and jamovi's install gate compares
# only rVersion ("4.6.0-x64") -- no OS -- so a macOS .jmo would install on
# Windows and then fail to load. One artifact per OS x CPU, named so users can
# tell them apart.
#
# How: build once with jmc against the local jamovi.app (the module itself is
# pure R, so its installed copy is the same on every OS), then swap the bundled
# dependencies for the binaries of each target, taken from the same snapshot
# jamovi's own compiler uses for that R (jamovi-compiler/snapshots.js). The
# Posit macOS binaries are byte-identical to what jmc bundles natively.
# There is no Windows-on-ARM jamovi; those users run the x64 build.
set -euo pipefail

PUBLISH=no
case "${1:-}" in
  "") ;;
  --publish) PUBLISH=yes ;;
  *) echo "usage: release.sh [--publish]" >&2; exit 1 ;;
esac

HERE="$(cd "$(dirname "$0")/.." && pwd)"
MODULE=SNPstats
VERSION="$(awk -F': *' '$1 == "Version" { print $2; exit }' "$HERE/DESCRIPTION")"
APP="${JAMOVI_APP:-/Applications/jamovi.app}"

STAMP="$(sed -n 's/^JAMOVI_R_VERSION=//p' "$APP/Contents/Resources/env.conf" | tr -d '"[:space:]')"
R_VERSION="${STAMP%-*}"
R_MINOR="${R_VERSION%.*}"
# the snapshot date jmc pins for this R (jamovi-compiler/snapshots.js)
case "$R_VERSION" in
  4.6.0) SNAPSHOT=2026-05-11 ;;
  *) echo "error: no snapshot known for R $R_VERSION -- add it from snapshots.js" >&2; exit 1 ;;
esac
POSIT="https://packagemanager.posit.co/cran/$SNAPSHOT"

WORK="$HERE/.tmp/release"
DIST="$HERE/dist"
rm -rf "$WORK"; mkdir -p "$WORK" "$DIST"
rm -f "$DIST/${MODULE}_${VERSION}_"*.jmo

# ── 1. one base build ────────────────────────────────────────────────────────
NODE="$(Rscript --vanilla -e 'cat(system.file("node-darwin", "bin", "node", package = "node"))')"
JMC="$(Rscript --vanilla -e 'cat(system.file("node_modules", "jamovi-compiler", "index.js", package = "jmvtools"))')"
[ -x "$NODE" ] && [ -f "$JMC" ] || { echo "error: jmvtools/node not available to Rscript" >&2; exit 1; }

BASE="$WORK/base.jmo"
echo ">> building $MODULE $VERSION against $(basename "$APP") (R $R_VERSION)"
# jmc --build writes the .jmo relative to the cwd; give it an absolute path
"$NODE" "$JMC" --build "$HERE" --home "$APP" --jmo "$BASE"
[ -f "$BASE" ] || { echo "error: jmc did not produce $BASE" >&2; exit 1; }
# jmc regenerated R/snpPGS.h.R; restore the caseLevel default (CLAUDE.md)
bash "$HERE/tools/patch_h.sh" "$HERE/R/snpPGS.h.R" || [ $? -eq 10 ]

# the dependencies jmc decided to bundle = everything under R/ but the module
DEPS="$(unzip -Z1 "$BASE" | sed -n "s#^$MODULE/R/\([^/]*\)/DESCRIPTION\$#\1#p" | grep -vx "$MODULE" | sort -u | tr '\n' ' ')"
echo ">> bundled dependencies: $DEPS"

# ── 2. per-target dependency binaries ────────────────────────────────────────
# fetch PKG into DEST from the repo for one target; check it is a binary built
# for that platform (Posit silently serves source when it has no binary)
fetch() {
  local target="$1" pkg="$2" dest="$3" repo ua="" ext ver url file built
  case "$target" in
    windows-x64) repo="$POSIT/bin/windows/contrib/$R_MINOR"; ext=zip ;;
    macos-x64)   repo="$POSIT/bin/macosx/big-sur-x86_64/contrib/$R_MINOR"; ext=tgz ;;
    macos-arm64) repo="$POSIT/bin/macosx/sonoma-arm64/contrib/$R_MINOR"; ext=tgz ;;
    linux-x64)   repo="https://packagemanager.posit.co/cran/__linux__/noble/$SNAPSHOT/src/contrib"; ext=tar.gz
                 ua="R ($R_VERSION x86_64-pc-linux-gnu x86_64 linux-gnu)" ;;
    linux-arm64) repo="https://packagemanager.posit.co/cran/__linux__/noble/$SNAPSHOT/src/contrib"; ext=tar.gz
                 ua="R ($R_VERSION aarch64-unknown-linux-gnu aarch64 linux-gnu)" ;;
  esac
  ver="$(curl -fsL ${ua:+-A "$ua"} "$repo/PACKAGES" | awk -v RS= -v p="$pkg" '$0 ~ "(^|\n)Package: " p "\n"' \
         | sed -n 's/^Version: //p')"
  [ -n "$ver" ] || { echo "error: $pkg not in $repo" >&2; return 1; }
  file="$WORK/dl/$target/${pkg}_$ver.$ext"
  mkdir -p "$(dirname "$file")"
  curl -fsL ${ua:+-A "$ua"} -o "$file" "$repo/${pkg}_$ver.$ext"
  if [ "$ext" = zip ]; then unzip -q "$file" -d "$dest"; else tar xzf "$file" -C "$dest"; fi
  built="$(sed -n 's/^Built: //p' "$dest/$pkg/DESCRIPTION")"
  [ -n "$built" ] || { echo "error: $target $pkg $ver came as source, not a binary" >&2; return 1; }
  echo "   $pkg $ver  ($built)"
}

# expected platform triple of any native library, for the check below
platform_of() {
  case "$1" in
    windows-x64) echo x86_64-w64-mingw32 ;;
    macos-x64)   echo x86_64-apple-darwin ;;
    macos-arm64) echo aarch64-apple-darwin ;;
    linux-x64)   echo x86_64-pc-linux-gnu ;;
    linux-arm64) echo aarch64-unknown-linux-gnu ;;
  esac
}

for TARGET in windows-x64 macos-arm64 macos-x64 linux-x64 linux-arm64; do
  echo ">> $TARGET"
  STAGE="$WORK/$TARGET"
  mkdir -p "$STAGE"
  unzip -q "$BASE" -d "$STAGE"
  for PKG in $DEPS; do
    rm -rf "$STAGE/$MODULE/R/$PKG"
    fetch "$TARGET" "$PKG" "$STAGE/$MODULE/R"
    # a compiled package must carry this target's platform in Built:
    if [ -d "$STAGE/$MODULE/R/$PKG/libs" ]; then
      grep -q "^Built: R $R_VERSION; $(platform_of "$TARGET")" "$STAGE/$MODULE/R/$PKG/DESCRIPTION" || {
        echo "error: $PKG is not an R $R_VERSION $(platform_of "$TARGET") build" >&2; exit 1; }
    fi
  done
  # rVersion is jamovi's only compatibility gate: R version + CPU
  CPU="${TARGET#*-}"
  for META in "$STAGE/$MODULE/jamovi.yaml" "$STAGE/$MODULE/jamovi-full.yaml"; do
    perl -0pi -e "s/^rVersion: .*\$/rVersion: $R_VERSION-$CPU/m" "$META"
    grep -qx "rVersion: $R_VERSION-$CPU" "$META" || { echo "error: rVersion not set in $META" >&2; exit 1; }
  done
  OUT="$DIST/${MODULE}_${VERSION}_${TARGET}.jmo"
  ( cd "$STAGE" && zip -q -r -X "$OUT" "$MODULE" )
  echo "   wrote dist/$(basename "$OUT")"
done

[ "$PUBLISH" = yes ] || { echo ">> built; re-run with --publish to release"; exit 0; }

# ── 3. publish ───────────────────────────────────────────────────────────────
TAG="v$VERSION"
REPO="$(cd "$HERE" && gh repo view --json nameWithOwner -q .nameWithOwner)"
f() { echo "\`${MODULE}_${VERSION}_$1.jmo\`"; }
NOTES="$(cat <<NOTESEOF
$MODULE $VERSION for **jamovi 28.4 or newer** (R $R_VERSION). Pick the file for your computer:

| your computer | file |
| --- | --- |
| Windows (Intel/AMD, and Windows on ARM) | $(f windows-x64) |
| macOS, Apple silicon (M1 and later) | $(f macos-arm64) |
| macOS, Intel | $(f macos-x64) |
| Linux (flatpak), x86_64 | $(f linux-x64) |
| Linux (flatpak), ARM64 | $(f linux-arm64) |

The files differ only in the compiled haplo.stats library they bundle. jamovi
checks the R version and the CPU but not the operating system, so it would
accept a file for the wrong OS and the analyses would then fail to load: use
the file named for your system.

Install: **Modules (+) -> Sideload**, then select the \`.jmo\`.

Bundled: $DEPS(binaries from the Posit CRAN snapshot $SNAPSHOT that jamovi uses for R $R_VERSION).
NOTESEOF
)"

echo ">> releasing $TAG on $REPO"
gh release delete "$TAG" --repo "$REPO" --cleanup-tag --yes >/dev/null 2>&1 || true
gh release create "$TAG" --repo "$REPO" --title "$MODULE $VERSION" --notes "$NOTES" \
  "$DIST/${MODULE}_${VERSION}_"*.jmo
echo ">> done: https://github.com/$REPO/releases/tag/$TAG"
