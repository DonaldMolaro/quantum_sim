# Quantum Counting (CLI Demo)

Quantum Counting estimates how many marked states exist in a search space,
without enumerating every state.

In this simulator, `QCOUNT` estimates the number of marked targets using
Grover-iteration success probabilities and a fitted angle `theta`.

## Why it is useful

- Grover Search answers: "find one marked state."
- Quantum Counting answers: "how many marked states are there?"

This is useful when you do not know target count in advance, or you want to
choose better Grover iteration counts.

## Command syntax

```
QCOUNT <n_qubits> <target1> [target2 ...]
QCOUNT RUN <n_qubits> <fit_iterations> <target1> [target2 ...]
```

Where:
- `n_qubits` defines search space size `N = 2^n`
- each target is a marked basis-state index in `[0, N-1]`

## 1) Built-in demo

```
QCOUNT DEMO
```

This runs a fixed example with `n=3`, `targets={1,6}`.

## 2) Manual examples

Single marked state in 3 qubits (`N=8`):

```
QCOUNT 3 5
```

Two marked states in 3 qubits:

```
QCOUNT 3 1 6
```

Expected estimate near `2`.

Explicit fitting window (advanced mode):

```
QCOUNT RUN 3 6 1 6
```

This forces the estimator to use exactly 6 Grover-iteration samples (`k=0..5`).

## 3) "By hand" counting intuition (2 qubits)

For `n=2`, let the only marked state be `|11>`. We can run one Grover iteration
manually and inspect the amplitude concentration.

```
INIT 2
H 0
H 1
CZ 0 1
H 0
H 1
X 0
X 1
CZ 0 1
X 0
X 1
H 0
H 1
DISPLAY
```

If probability mass concentrates heavily onto one basis state after one iteration,
that behavior is consistent with `M=1` marked state out of `N=4`.

## Output fields

- `estimated_targets`: rounded integer estimate
- `real`: non-rounded estimate
- `estimated_theta`: fitted Grover angle
- `iterations_used`: number of Grover iterations sampled by the estimator

## Notes / limits

- Targets must be in range for `n_qubits`.
- This is an educational estimator in simulator style, not a fault-tolerant
  phase-estimation implementation.
