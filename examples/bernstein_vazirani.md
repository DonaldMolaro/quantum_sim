# Bernstein-Vazirani Usage (CLI Demo)

This example shows how to recover a hidden bitstring `s` from a linear oracle
`f(x)=s.x xor b` with one query.

Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

## 1) Start the CLI

```
./quantum_sim
```

## 2) Run BV with a hidden string

```
BV 5 22
```

`22` is binary `10110`, so the measured secret should be `22`.

## 3) Run BV with a bias bit

```
BV 5 22 1
```

Expected: measured secret is still `22` (global bias does not change recovery).

## 4) By hand (gate-level), n=4, secret s=11 (binary 1011)

This reproduces BV manually with input qubits `q0..q3` and ancilla `q4`.

`s=1011` means CNOT controls on `q0`, `q1`, and `q3`.

```
INIT 5 5
X 4
H 0
H 1
H 2
H 3
H 4
CX 0 4
CX 1 4
CX 3 4
H 0
H 1
H 2
H 3
MEASURE 0 0
MEASURE 1 1
MEASURE 2 2
MEASURE 3 3
DISPLAY
```

Expected input-register bits: `1011` (LSB-first qubit indexing in the simulator).
