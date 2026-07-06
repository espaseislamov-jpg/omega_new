#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-extracted_desktop}"
ARCHIVE="${2:-Desktop.part1.rar}"

if ! command -v unar >/dev/null 2>&1; then
  echo "The 'unar' extractor is required for these RAR5 split archives." >&2
  echo "Install it with: apt-get update && apt-get install -y unar" >&2
  exit 1
fi

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR.tmp"
unar -o "$OUT_DIR.tmp" "$ARCHIVE"

# unar creates a folder named after the archived Desktop directory.
if [ -d "$OUT_DIR.tmp/Desktop" ]; then
  mv "$OUT_DIR.tmp/Desktop" "$OUT_DIR"
else
  mkdir -p "$OUT_DIR"
  find "$OUT_DIR.tmp" -mindepth 1 -maxdepth 1 -exec mv {} "$OUT_DIR" \;
fi
rm -rf "$OUT_DIR.tmp"
find "$OUT_DIR" -maxdepth 1 -type f -print | sort
