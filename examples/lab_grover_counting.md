# Guided Lab: Grover + Quantum Counting

This lab shows how counting informs Grover iteration choice.

## Goal

For a 3-qubit search space with targets `{3, 5}`:
- Estimate marked-count with `QCOUNT`.
- Compare manual Grover run vs `GROVER AUTO`.

## Setup

```bash
./quantum_sim
```

Enable guided narration:

```
TUTOR ON
```

## Part A: Estimate marked-count

Run fixed-window counting:

```
QCOUNT RUN 3 6 3 5
```

Checkpoint:
- `estimated_targets` should be near `2`.

## Part B: Standard Grover run

```
GROVER 3 5
```

Checkpoint:
- Output reports iterations used and target amplification behavior.

## Part C: Auto-tuned Grover run

```
GROVER AUTO 3 6 3 5
```

Checkpoint:
- Auto mode should report counting estimate and selected iteration count.

## Part D: Compare behavior

Compare:
- iteration count from `GROVER`
- iteration count from `GROVER AUTO`
- whether both produce strong amplification on targets

## Extension

Repeat with one target and then three targets to see how iteration choice changes
with marked-count.
