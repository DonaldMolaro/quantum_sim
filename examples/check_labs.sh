#!/bin/sh
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$ROOT_DIR"

make clean >/dev/null
make quantum_sim >/dev/null

run_case() {
  name="$1"
  expected="$2"
  cmds="$3"
  printf "[LAB] %s\n" "$name"
  out="$(printf "%s\nQUIT\n" "$cmds" | ./quantum_sim 2>&1 || true)"
  printf "%s" "$out" | grep -q "$expected" || {
    echo "[FAIL] $name"
    echo "Expected to find: $expected"
    echo "$out"
    exit 1
  }
  echo "[OK] $name"
}

run_case "Bell prep" "CX(0, 1) applied." \
"SEED 7
TUTOR ON
INIT 2
H 0
CX 0 1
DISPLAY"

run_case "Grover + Counting" "Grover AUTO iterations used:" \
"SEED 11
TUTOR ON
QCOUNT RUN 3 6 3 5
GROVER AUTO 3 6 3 5"

run_case "QAOA lab" "QAOA result:" \
"SEED 99
TUTOR ON
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
QAOA QUBO 3 1 0 20 0.25 -2 0 2 0 1 0 2 0 -3"

run_case "VQE lab" "VQE result:" \
"SEED 42
TUTOR ON
VQE RUN 1 1 30 0.3 0 1 1.0 1 Z 0"

printf "[LAB] Shor lab\n"
shor_out="$(printf "SEED 1234\nTUTOR ON\nSHOR 15\nQUIT\n" | \
  QSIM_SHOR_MAX_ATTEMPTS=1 QSIM_SHOR_FORCE_A=7 QSIM_SHOR_FORCE_X=1 QSIM_SHOR_FORCE_NC=4 QSIM_SHOR_FORCE_R=4 ./quantum_sim 2>&1 || true)"
printf "%s" "$shor_out" | grep -q "Non-trivial factors found:" || {
  echo "[FAIL] Shor lab"
  echo "$shor_out"
  exit 1
}
echo "[OK] Shor lab"

echo "All guided lab checks passed."
