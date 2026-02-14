#!/bin/sh
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
SRC_MD="$ROOT_DIR/examples/visual_walkthroughs.md"
OUT_DIR="${1:-$ROOT_DIR/examples/diagrams}"
TMP_DIR="$OUT_DIR/.tmp_mermaid"

mkdir -p "$OUT_DIR" "$TMP_DIR"
rm -f "$TMP_DIR"/*.mmd "$OUT_DIR"/*.svg

awk -v outdir="$TMP_DIR" '
  BEGIN {
    in_mermaid = 0;
    idx = 0;
  }
  /^```mermaid[[:space:]]*$/ {
    in_mermaid = 1;
    idx++;
    file = sprintf("%s/diagram_%02d.mmd", outdir, idx);
    next;
  }
  /^```[[:space:]]*$/ {
    if (in_mermaid == 1) {
      in_mermaid = 0;
      close(file);
    }
    next;
  }
  {
    if (in_mermaid == 1) {
      print $0 >> file;
    }
  }
' "$SRC_MD"

count="$(find "$TMP_DIR" -name 'diagram_*.mmd' | wc -l | tr -d ' ')"
if [ "$count" = "0" ]; then
  echo "No mermaid blocks found in $SRC_MD"
  exit 1
fi

if ! command -v mmdc >/dev/null 2>&1; then
  echo "Mermaid CLI (mmdc) not found."
  echo "Install with: npm i -g @mermaid-js/mermaid-cli"
  echo "Extracted .mmd files are in: $TMP_DIR"
  exit 2
fi

for f in "$TMP_DIR"/diagram_*.mmd; do
  base="$(basename "$f" .mmd)"
  out="$OUT_DIR/$base.svg"
  mmdc -i "$f" -o "$out" >/dev/null
  echo "Rendered: $out"
done

rm -f "$TMP_DIR"/*.mmd
rmdir "$TMP_DIR" 2>/dev/null || true
echo "Rendered $count diagram(s) into $OUT_DIR"
