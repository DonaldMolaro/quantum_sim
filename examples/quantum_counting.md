# Quantum Counting (CLI Demo)

Quantum Counting estimates how many marked states exist in a search space,
without enumerating every state.

In this simulator, `QCOUNT` estimates the number of marked targets using
Grover-iteration success probabilities and a fitted angle `theta`.

## Classical Computing View

This is cardinality estimation for a predicate.

- Problem: estimate `M = |{x : is_target(x)=true}|` over `N=2^n`.
- Classical exact method: evaluate all `N` states.
- Classical approximate method: random sampling, then extrapolate.
- Quantum counting idea: infer `M` from Grover dynamics (the rotation angle
  depends on `M/N`) instead of direct enumeration.

In this simulator, counting is implemented with an educational fit over
Grover-success traces rather than full phase-estimation hardware logic.

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

## 4) Hand-tooled workflow for `QCOUNT RUN`

Use this sequence when you want to compare manual intuition with the built-in
fitter over a fixed number of samples.

1. Build a small marked-state experiment by hand (as above) and inspect with `DISPLAY`.
2. Run the estimator with explicit fitting samples:

```
QCOUNT RUN 2 4 3
```

3. Compare:
- Manual observation: one marked state appears strongly after one Grover step.
- Estimator output: `estimated_targets` should be close to `1`.

## Output fields

- `estimated_targets`: rounded integer estimate
- `real`: non-rounded estimate
- `estimated_theta`: fitted Grover angle
- `iterations_used`: number of Grover iterations sampled by the estimator

## Notes / limits

- Targets must be in range for `n_qubits`.
- This is an educational estimator in simulator style, not a fault-tolerant
  phase-estimation implementation.
