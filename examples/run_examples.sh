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
run_cli "QPE demo" "QPE DEMO"
run_cli "QPE T gate (pi/4, m=4)" "QPE 4 PI/4"
run_cli "QPE S gate (pi/2, m=4)" "QPE 4 PI/2"
run_cli "Sdg/Tdg gates" "INIT 1 0
S 0
SDG 0
X 0
T 0
TDG 0"
run_cli "P gate" "INIT 1 0
X 0
P 0 PI"
run_cli "CP gate" "INIT 2 0
X 0
X 1
CP 0 1 PI/2"
run_cli "MCX 3-control gate" "INIT 4 0
X 0
X 1
X 2
MCX 0 1 2 3"
run_cli "CSWAP gate" "INIT 3 0
X 0
X 2
CSWAP 2 0 1"
run_cli "NOISE + Bell state" "INIT 2 0
NOISE 0
H 0
CX 0 1
CHECK BELL PHI+"
run_cli "LOAD bell_state.qsim" "LOAD examples/bell_state.qsim"
run_cli "LOAD qpe example" "LOAD examples/qpe_t_gate.qsim"
run_cli "LOAD MCX oracle" "LOAD examples/mcx_oracle.qsim"

echo "All example smoke runs completed."
