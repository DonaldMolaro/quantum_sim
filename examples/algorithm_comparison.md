# Algorithm Comparison Guide

Use this quick map to decide which algorithm to run first in the simulator.

## Decision table

- `DEUTSCH_JOZSA`: decide if oracle is constant or balanced in one query model.
- `BV`: recover hidden bit string from linear oracle.
- `SIMON`: recover xor-mask secret from 2-to-1 periodic oracle.
- `GROVER`: find one marked item in an unstructured search space.
- `QCOUNT`: estimate how many marked items exist.
- `SHOR`: factor composites (small educational instances here).
- `QAOA` / `VQE`: variational optimization workflows for objective/Hamiltonian problems.

## Recommended workflow for search problems

1. Start with `QCOUNT` when marked-count is unknown.
2. Use `GROVER AUTO` to choose iteration count from the estimate.
3. Use `GROVER` with explicit targets when marked-count is known.

Example:

```
QCOUNT RUN 3 6 3 5
GROVER AUTO 3 6 3 5
GROVER 3 5
```

## Recommended workflow for optimization problems

1. Run exact baseline (small instances):

```
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
```

2. Run variational approximation:

```
VQA QAOA 3 1 0 40 0.25 -2 0 2 0 1 0 2 0 -3
VQE RUN 2 2 50 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

3. Compare to the baseline and tune depth/iterations.
