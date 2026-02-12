# QUBO Usage (CLI Demo)

This example shows how to evaluate and solve a small QUBO in the simulator.

Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

## 1) Built-in demo

```
QUBO DEMO
```

This runs a fixed 3-variable QUBO, prints exact minimum, then runs Grover-threshold search.

## 2) Exact solve from matrix entries

Matrix is row-major (`n*n` entries):

```
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
```

Expected minimum at bitstring `0b100` (decimal `4`).

## 3) Grover threshold solve from matrix entries

Syntax:

```
QUBO GROVER <n> <threshold> <iterations> <n*n matrix entries>
```

Example:

```
QUBO GROVER 3 -3 1 -2 0 2 0 1 0 2 0 -3
```

This marks assignments with objective value `<= -3` and performs one Grover iteration.
