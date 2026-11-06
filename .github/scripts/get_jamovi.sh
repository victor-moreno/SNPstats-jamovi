#!/usr/bin/env bash
# Install the real jamovi for this runner and record how to run its R.
#
#   get_jamovi.sh <windows|macos-arm64|macos-x64|flatpak> <jamovi version>
#
# Appends R_CMD and ENV_CONF to $GITHUB_ENV. On Linux the version argument is
# ignored and the jamovi tested is whatever Flathub has (check_jmo.R reports it).
set -euo pipefail

APP="$1"
VERSION="$2"
WORK="${RUNNER_TEMP:-.tmp}/jamovi"
mkdir -p "$WORK"

fetch() {
  curl -fsSL --retry 3 --retry-all-errors -C - -o "$2" "https://archives.jamovi.org/$1"
}

case "$APP" in
  windows)
    fetch "jamovi-$VERSION-win-x64.exe" "$WORK/jamovi.exe"
    # PowerShell, because Git Bash would rewrite the installer's /S and /D= as paths
    powershell -NoProfile -Command "Start-Process -Wait -FilePath '$(cygpath -w "$WORK/jamovi.exe")' -ArgumentList '/S','/D=C:\jamovi'"
    CONF="$(find /c/jamovi -name env.conf | head -1)"
    RS="$(find /c/jamovi -name Rscript.exe -path '*x64*' | head -1)"
    CONF="$(cygpath -m "$CONF")"
    R_CMD="$(cygpath -m "$RS")"
    ;;
  macos-*)
    fetch "jamovi-$VERSION-$APP.dmg" "$WORK/jamovi.dmg"
    hdiutil attach -quiet -nobrowse -mountpoint "$WORK/mnt" "$WORK/jamovi.dmg"
    # jamovi's R expects to live where a normal drag-install puts it
    cp -R "$WORK/mnt/jamovi.app" /Applications/
    hdiutil detach -quiet "$WORK/mnt"
    CONF=/Applications/jamovi.app/Contents/Resources/env.conf
    RH=/Applications/jamovi.app/Contents/Frameworks/R.framework/Versions/Current/Resources
    R_CMD="$RH/bin/Rscript"
    # jamovi's Rscript has /Library/Frameworks/R.framework built in: without RHOME it
    # fails, or on a machine with its own R, silently runs that R instead
    echo "RHOME=$RH" >> "$GITHUB_ENV"
    ;;
  flatpak)
    sudo apt-get update -qq
    sudo apt-get install -y -qq flatpak
    flatpak remote-add --user --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
    flatpak install --user -y --noninteractive flathub org.jamovi.jamovi
    # the check runs inside the flatpak sandbox, so these are sandbox paths
    flatpak run --command=sh org.jamovi.jamovi -c 'find /app -name env.conf -o -name Rscript | head; cat "$(find /app -name env.conf | head -1)"'
    CONF="$(flatpak run --command=sh org.jamovi.jamovi -c 'find /app -name env.conf | head -1')"
    RS="$(flatpak run --command=sh org.jamovi.jamovi -c 'find /app -name Rscript -type f | head -1')"
    R_CMD="flatpak run --filesystem=$GITHUB_WORKSPACE --filesystem=$RUNNER_TEMP --command=$RS org.jamovi.jamovi"
    ;;
  *)
    echo "unknown app: $APP" >&2
    exit 1
    ;;
esac

[ -n "$CONF" ] && [ -n "$R_CMD" ] || { echo "error: jamovi layout not recognised" >&2; exit 1; }
{
  echo "ENV_CONF=$CONF"
  echo "R_CMD=$R_CMD"
} >> "$GITHUB_ENV"
