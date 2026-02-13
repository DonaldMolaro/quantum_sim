# Max-Cut Example (via QUBO)

This example shows how to solve a small Max-Cut instance using the simulator's QUBO tools.

## Problem

Use a triangle graph with 3 vertices and unit edge weights:
- edges: `(0,1)`, `(1,2)`, `(0,2)`

Let each bit `x_i` indicate which side of the cut vertex `i` belongs to (`0` or `1`).

For this graph, maximizing cut size is equivalent to minimizing this QUBO:

```
Q = [ -2   1   1
       1  -2   1
       1   1  -2 ]
```

Row-major entries for CLI:

```
-2 1 1 1 -2 1 1 1 -2
```

## 1) Solve exactly

```
QUBO EXACT 3 -2 1 1 1 -2 1 1 1 -2
```

Expected:
- minimum objective value is `-2`
- multiple optimal bitstrings (all non-trivial partitions)

## 2) Solve with VQA (QAOA)

```
VQA QAOA 3 1 0 40 0.25 -2 1 1 1 -2 1 1 1 -2
```

## 3) Solve with annealing

```
ANNEAL QUBO SQA 3 80 20 0.1 6.0 8 -2 1 1 1 -2 1 1 1 -2
```

## 4) Interpret a result

If the solver returns bitstring `101`:
- side A: vertices with bit `1` -> `{0,2}`
- side B: vertices with bit `0` -> `{1}`

Edges crossing the partition define the cut.
