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

Use this 3-variable QUBO matrix in all steps:

`-2 0 2 0 1 0 2 0 -3`

## Part A: Exact reference

```
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
```

Checkpoint:
- Record the exact minimum energy and assignment.

## Part B: Shallow QAOA

```
QAOA QUBO 3 1 0 20 0.25 -2 0 2 0 1 0 2 0 -3
```

Checkpoint:
- Record `best_energy`.

## Part C: Deeper/longer QAOA

```
QAOA QUBO 3 2 0 60 0.20 -2 0 2 0 1 0 2 0 -3
```

Checkpoint:
- Compare with Part B and exact reference.

## Extension

Run with finite shots (`shots=128`) and compare noise effects against `shots=0`.
