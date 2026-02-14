# Performance and Scaling Guide

This page is for engineers deciding when to use each algorithm in this simulator.

## Practical framing

- State simulation cost is dominated by `2^n` state size.
- Most operations are linear scans over non-zero amplitudes.
- Educational demos are reliable for small/medium `n`; large `n` quickly becomes
  memory and time bound.

## Complexity Cheat Sheet

`N = 2^n` basis states, `M` marked states.

- Core gates (`X`, `H`, rotations, controlled gates):
  - Time: typically `O(N)` per gate
  - Memory: `O(N)` in dense behavior
- `GROVER`:
  - Iterations: `R ~= floor((pi/4)*sqrt(N/M))`
  - Total: `O(R * N)` simulator work
- `QCOUNT`:
  - Runs multiple Grover depths and fits theta
  - Total: roughly `O(k * N)` where `k = fit_iterations`
- `SHOR` (educational):
  - Dominated by quantum-register simulation and QFT/mod-exp steps
  - Practical only for small composites in this implementation
- `QUBO EXACT`:
  - `O(2^n)` objective evaluations
- `QAOA`/`VQE`:
  - `O(iterations * layers * N)`-style simulator scaling (plus objective terms)

## Memory rule of thumb

Rough lower bound for dense complex-state storage:

- `N` amplitudes × 16 bytes (complex<double>) ≈ `16 * 2^n` bytes
- Real simulator overhead is higher due to containers/metadata.

Approximate floor:

- `n=20`: ~16 MB raw amplitudes
- `n=24`: ~256 MB raw amplitudes
- `n=28`: ~4 GB raw amplitudes

## Recommended operating ranges (education/prototyping)

- Gate demos / Bell / DJ / BV / Simon: `n <= 12` very comfortable
- Grover / counting demos: typically `n <= 16` (depending on machine)
- Exact QUBO/TSP encodings: keep variable count low (`n <= ~20` bits for exact)
- Shor demo: small composite numbers only

## Decision checklist

1. Need guaranteed optimum and small search space? Use `QUBO EXACT`.
2. Need scalable approximation on larger binary objective? Use `QAOA`/`ANNEAL`.
3. Need unstructured search demo? Use `GROVER`, optionally `GROVER AUTO`.
4. Need target cardinality estimate? Use `QCOUNT` first.
