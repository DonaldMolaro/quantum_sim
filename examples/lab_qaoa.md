# Guided Lab: QAOA on QUBO

This lab compares exact QUBO optimization with variational QAOA approximations.

## Goal

Understand approximation quality versus depth and optimization budget.

## Setup

```bash
./quantum_sim
```

Recommended:

```
SEED 99
TUTOR ON
```

If you want to validate the expected checkpoints after working through the lab:

```bash
./examples/check_labs.sh qaoa
```

Use this 3-variable QUBO matrix in all steps:

`-2 0 2 0 1 0 2 0 -3`

## Part A: Exact reference

Predict first:
- Which result from this section should every later variational run be compared
  against?

```
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
```

Checkpoint:
- Record the exact minimum energy and assignment.
- In the current seeded lab path, the exact minimum is `-3` at `0b100`.

## Part B: Shallow QAOA

Predict first:
- Should a shallow variational run exactly match the brute-force optimum every
  time?

```
QAOA QUBO 3 1 0 20 0.25 -2 0 2 0 1 0 2 0 -3
```

Checkpoint:
- Record `best_energy`.
- Compare it with the exact minimum from Part A.

## Part C: Deeper/longer QAOA

Predict first:
- If we increase both depth and optimization budget, what should usually happen
  to `best_energy`?

```
QAOA QUBO 3 2 0 60 0.20 -2 0 2 0 1 0 2 0 -3
```

Checkpoint:
- Compare with Part B and exact reference.
- You should still see the exact reference block so the approximation gap is easy to inspect.

## Extension

Run with finite shots (`shots=128`) and compare noise effects against `shots=0`.
