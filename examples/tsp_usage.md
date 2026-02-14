# TSP Usage (CLI Demo)

This example shows how the simulator solves a small Traveling Salesman Problem
by encoding it as QUBO (fixed start city `0`) and running exact search.

Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

## Classical Computing View

TSP is a constrained combinatorial search problem.

- Objective: minimize route cost visiting each city exactly once.
- Classical exact methods scale poorly (factorial growth or exponential DP).
- QUBO reformulation converts constraints + objective into one binary energy.

In this simulator, the value is educational: see how a familiar NP-hard problem
becomes a binary optimization instance that can be solved by multiple backends.

## 1) Built-in TSP demo

```
TSP DEMO
```

This runs a built-in 4-city distance matrix and prints:
- decoded route (starting/ending at city `0`)
- route cost
- QUBO energy

## 2) Exact solve from a distance matrix

Syntax:

```
TSP EXACT <n_cities> <penalty> <n*n distance entries row-major>
```

Notes:
- use diagonal zeros (`d[i][i] = 0`)
- set `penalty` to `-1` to auto-select a penalty

Example 4-city matrix:

```
TSP EXACT 4 -1 0 1 2 1 1 0 1 2 2 1 0 1 1 2 1 0
```

One optimal tour for this instance is cost `4`.

## 3) Practical size limits

The current fixed-start QUBO encoding uses `(n-1)^2` variables.
With the simulator bitstring cap, that means:
- hard cap around `n <= 8`
- practical demo range is smaller (typically `4-5` cities)
