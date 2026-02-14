# VQA (QAOA) Usage for QUBO

This example shows how to run the new VQA/QAOA workflow against a QUBO matrix.

Optional: enable algorithm logs:

```
VERBOSE VERBOSE
```

## 1) Built-in VQA demo

```
VQA DEMO
```

This runs QAOA on a fixed 3-variable QUBO and prints:
- VQA best energy/state
- exact minimum (reference from brute force)

## 2) Run QAOA directly from matrix entries

Syntax:

```
VQA QAOA <n> <p> <shots> <iters> <step> <n*n matrix entries>
```

- `n`: number of binary variables
- `p`: QAOA layers
- `shots`: measurement shots (`0` = exact expectation)
- `iters`: optimizer iterations
- `step`: initial step size

Example (`n=3`, same matrix as QUBO example):

```
VQA QAOA 3 1 0 40 0.25 -2 0 2 0 1 0 2 0 -3
```

## 3) Shot-based (sampling) mode

Use `shots > 0` to simulate finite-shot expectation estimates:

```
VQA QAOA 3 1 128 40 0.25 -2 0 2 0 1 0 2 0 -3
```

Use larger shot counts for less noisy estimates.

## 4) Hand-tooled QAOA workflow

This mirrors what the built-in optimizer does, but in a controlled sequence:

1. Solve exact baseline first:

```
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
```

2. Run shallow QAOA:

```
VQA QAOA 3 1 0 20 0.25 -2 0 2 0 1 0 2 0 -3
```

3. Increase depth and iterations:

```
VQA QAOA 3 2 0 60 0.20 -2 0 2 0 1 0 2 0 -3
```

4. Compare best energy from each run against the exact baseline.

This gives a manual, reproducible tuning loop for `p`, step size, and optimizer
budget.
