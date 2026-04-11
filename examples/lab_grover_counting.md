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

If you want to validate the expected checkpoints after working through the lab:

```bash
./examples/check_labs.sh grover
```

## Part A: Estimate marked-count

Predict first:
- There are two targets in an eight-state space. What should `estimated_targets`
  be close to?

Run fixed-window counting:

```
QCOUNT RUN 3 6 3 5
```

Checkpoint:
- `estimated_targets` should be near `2`.

## Part B: Standard Grover run

Predict first:
- With two marked states out of eight, should Grover need more than one
  iteration here?

```
GROVER 3 5
```

Checkpoint:
- Output reports iterations used and target amplification behavior.
- For this setup, `Grover iterations used:` should be `1`.

## Part C: Auto-tuned Grover run

Predict first:
- If counting reports about two targets, should auto mode choose a different
  iteration count from Part B?

```
GROVER AUTO 3 6 3 5
```

Checkpoint:
- Auto mode should report counting estimate and selected iteration count.
- For this setup, `Grover AUTO iterations used:` should also be `1`.

## Part D: Compare behavior

Compare:
- iteration count from `GROVER`
- iteration count from `GROVER AUTO`
- whether both produce strong amplification on targets

## Extension

Repeat with one target and then three targets to see how iteration choice changes
with marked-count.
