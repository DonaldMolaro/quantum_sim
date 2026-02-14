#!/bin/sh
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$ROOT_DIR"

if [ ! -x ./quantum_sim ]; then
  echo "Building quantum_sim..."
  make clean >/dev/null
  make quantum_sim >/dev/null
fi

run_cli() {
  name="$1"
  cmds="$2"
  printf "[RUN] %s\n" "$name"
  printf "%s\nQUIT\n" "$cmds" | ./quantum_sim >/tmp/qsim_example.out 2>/tmp/qsim_example.err || {
    echo "[FAIL] $name"
    cat /tmp/qsim_example.err
    exit 1
  }
  echo "[OK] $name"
}

run_cli "Grover basic" "GROVER 3 5"
run_cli "Grover auto-tuned" "GROVER AUTO 3 6 3 5"
run_cli "Quantum counting" "QCOUNT RUN 3 6 1 6"
run_cli "Simon" "SIMON 4 10 8"
run_cli "VQA QAOA" "VQA QAOA 3 1 0 10 0.25 -2 0 2 0 1 0 2 0 -3"
run_cli "VQE" "VQE RUN 1 1 10 0.3 0 1 1.0 1 Z 0"

echo "All example smoke runs completed."
