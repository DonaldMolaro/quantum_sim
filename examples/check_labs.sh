#!/bin/sh
set -eu

ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
cd "$ROOT_DIR"

for arg in "$@"; do
  case "$arg" in
    all|bell|grover|qaoa|vqe|shor)
      ;;
    *)
      echo "Unknown lab selector: $arg" >&2
      echo "Expected one of: all bell grover qaoa vqe shor" >&2
      exit 1
      ;;
  esac
done

if [ ! -x ./quantum_sim ]; then
  make clean >/dev/null
  make quantum_sim >/dev/null
fi

want_lab() {
  if [ "$#" -eq 0 ]; then
    return 0
  fi

  lab="$1"
  shift
  for arg in "$@"; do
    if [ "$arg" = "$lab" ] || [ "$arg" = "all" ]; then
      return 0
    fi
  done
  return 1
}

require_output() {
  name="$1"
  out="$2"
  needle="$3"
  printf "%s" "$out" | grep -F -q "$needle" || {
    echo "[FAIL] $name"
    echo "Expected to find: $needle"
    echo "$out"
    exit 1
  }
}

run_case() {
  name="$1"
  cmds="$2"
  shift 2
  printf "[LAB] %s\n" "$name"
  out="$(printf "%s\nQUIT\n" "$cmds" | ./quantum_sim 2>&1 || true)"
  for needle in "$@"; do
    require_output "$name" "$out" "$needle"
  done
  echo "[OK] $name"
}

run_case_env() {
  name="$1"
  env_prefix="$2"
  cmds="$3"
  shift 3
  printf "[LAB] %s\n" "$name"
  out="$(printf "%s\nQUIT\n" "$cmds" | env $env_prefix ./quantum_sim 2>&1 || true)"
  for needle in "$@"; do
    require_output "$name" "$out" "$needle"
  done
  echo "[OK] $name"
}

run_bell_lab() {
  run_case "Bell states lab" \
"TUTOR ON
INIT 2
H 0
CX 0 1
CHECK BELL PHI+
INIT 2
H 0
CX 0 1
Z 0
CHECK BELL PHI-
INIT 2
H 0
CX 0 1
X 1
CHECK BELL PSI+
INIT 2
H 0
CX 0 1
X 1
Z 0
CHECK BELL PSI-" \
    "CHECK BELL PHI+: PASS" \
    "CHECK BELL PHI-: PASS" \
    "CHECK BELL PSI+: PASS" \
    "CHECK BELL PSI-: PASS"
}

run_grover_lab() {
  run_case "Grover + Counting lab" \
"SEED 11
TUTOR ON
QCOUNT RUN 3 6 3 5
GROVER 3 5
GROVER AUTO 3 6 3 5" \
    "estimated_targets=2" \
    "Grover iterations used: 1" \
    "Grover AUTO iterations used: 1"
}

run_qaoa_lab() {
  run_case "QAOA lab" \
"SEED 99
TUTOR ON
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
QAOA QUBO 3 1 0 20 0.25 -2 0 2 0 1 0 2 0 -3
QAOA QUBO 3 2 0 60 0.20 -2 0 2 0 1 0 2 0 -3" \
    "QUBO exact result: min=-3 at 0b100 (4)" \
    "QAOA result:" \
    "Reference exact minimum:" \
    "min_state=0b100 (4)"
}

run_vqe_lab() {
  run_case "VQE lab" \
"SEED 42
TUTOR ON
VQE RUN 1 1 30 0.3 0 1 1.0 1 Z 0
VQE RUN 2 1 30 0.3 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
VQE RUN 2 2 80 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1" \
    "best_energy=-0.999" \
    "best_energy=-2.06048" \
    "best_energy=-2.06155"
}

run_shor_lab() {
  run_case_env "Shor lab" \
    "QSIM_SHOR_MAX_ATTEMPTS=1 QSIM_SHOR_FORCE_A=7 QSIM_SHOR_FORCE_X=1 QSIM_SHOR_FORCE_NC=4 QSIM_SHOR_FORCE_R=4" \
"SEED 1234
TUTOR ON
SHOR 15" \
    "Non-trivial factors found:"
}

if want_lab bell "$@"; then
  run_bell_lab
fi
if want_lab grover "$@"; then
  run_grover_lab
fi
if want_lab qaoa "$@"; then
  run_qaoa_lab
fi
if want_lab vqe "$@"; then
  run_vqe_lab
fi
if want_lab shor "$@"; then
  run_shor_lab
fi

echo "All requested guided lab checks passed."
