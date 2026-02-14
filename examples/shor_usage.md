# Shor's Algorithm Usage (CLI Demo)

This example shows how to use the built-in Shor demo from the CLI. The simulator
implements a small-scale, educational version of Shor's algorithm and is intended
for small composite numbers.
Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

## Classical Computing View

Factoring is easy to verify but hard to solve quickly in general.

- Input: composite `N`.
- Goal: find non-trivial factors `p, q` such that `N = p*q`.
- Classical baseline in this project size range: trial division / heuristics.
- Shor structure: reduce factoring to order-finding:
  - choose `a` with `gcd(a,N)=1`
  - find period `r` where `a^r mod N = 1`
  - use `gcd(a^(r/2) ± 1, N)` to extract factors

The quantum part accelerates the period-finding subproblem; the final factor
extraction is classical number theory.

## 1) Start the CLI

```
./quantum_sim
```

## 2) Run the Shor demo

```
SHOR 15
```

Expected: the demo prints the quantum order-finding process and suggests factors.

## 3) Try other small composites

```
SHOR 21
SHOR 33
SHOR 35
```

## Notes

- The demo is best suited for **small** semiprimes (product of two primes).
- Results may vary depending on random choices made during the run.
- For repeatability, you can set an environment seed when running tests; the
  interactive CLI does not force a seed.

## Suggested Workflow

If you want to inspect more detail:

1. Run `SHOR N` for a small composite.
2. Read the order-finding output and the measured phase.
3. Verify the factors manually (e.g., 15 = 3 * 5).
