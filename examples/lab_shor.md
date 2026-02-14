# Guided Lab: Shor (Educational Scale)

This lab demonstrates period-finding flow and factor extraction for a small
composite.

## Goal

Factor `N=15` and connect output to the period-finding concept.

## Setup

```bash
./quantum_sim
```

Enable deterministic behavior and guided narration:

```
SEED 1234
TUTOR ON
```

## Part A: Run Shor demo

```
SHOR 15
```

Checkpoint:
- Output should report non-trivial factors for `15` (e.g., `3 x 5`).

## Part B: Try another small semiprime

```
SHOR 21
```

Checkpoint:
- If factors are found, verify they multiply back to `21`.
- If a run fails to find factors, re-run once with same seed for reproducibility checks,
  then change seed and compare behavior.

## Discussion prompts

- Why does order-finding drive factorization?
- Why are repeated attempts needed in practical/educational runs?
