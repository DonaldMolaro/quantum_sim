# Guided Lab: VQE

This lab walks through manual VQE tuning from simple to richer ansatz.

## Goal

Observe how optimization budget and ansatz depth affect `best_energy`.

## Setup

```bash
./quantum_sim
```

Recommended:

```
SEED 42
TUTOR ON
```

## Part A: Single-qubit baseline

```
VQE RUN 1 1 30 0.3 0 1 1.0 1 Z 0
```

Checkpoint:
- `best_energy` should be close to `-1`.

## Part B: Two-qubit toy Hamiltonian

```
VQE RUN 2 1 30 0.3 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

Checkpoint:
- Record `best_energy`.

## Part C: Increase expressivity and budget

```
VQE RUN 2 2 80 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

Checkpoint:
- Compare with Part B; deeper ansatz often improves minimum energy found.

## Extension

Sweep `layers` in `{1,2,3}` and plot `best_energy` versus runtime/evaluations.
