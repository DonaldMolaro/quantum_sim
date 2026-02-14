#!/bin/sh
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
BUILD_DIR="${1:-$ROOT_DIR/build}"
OUT_CSV="${2:-$ROOT_DIR/benchmarks/grover_bench.csv}"

mkdir -p "$(dirname "$OUT_CSV")"

if [ ! -x "$BUILD_DIR/all_tests" ]; then
  echo "Building all_tests in $BUILD_DIR..."
  cmake -S "$ROOT_DIR" -B "$BUILD_DIR" >/dev/null
  cmake --build "$BUILD_DIR" -j >/dev/null
fi

TMP_OUT="$(mktemp)"
trap 'rm -f "$TMP_OUT"' EXIT

(cd "$BUILD_DIR" && QSIM_GROVER_BENCH=1 ./all_tests) >"$TMP_OUT" 2>&1

{
  echo "n_qubits,marked_count,iterations,expected_success,empirical_success"
  awk '
    /^Grover bench n=/ {
      n=""; m=""; r=""; e=""; p="";
      for (i=1; i<=NF; ++i) {
        if ($i ~ /^n=/) { n=substr($i,3); }
        else if ($i ~ /^M=/) { m=substr($i,3); }
        else if ($i ~ /^R=/) { r=substr($i,3); }
        else if ($i ~ /^expected=/) { e=substr($i,10); }
        else if ($i ~ /^empirical=/) { p=substr($i,11); }
      }
      if (n != "" && m != "" && r != "" && e != "" && p != "") {
        print n "," m "," r "," e "," p;
      }
    }
  ' "$TMP_OUT"
} > "$OUT_CSV"

echo "Wrote benchmark CSV: $OUT_CSV"
