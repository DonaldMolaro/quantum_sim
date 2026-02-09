# Shor's Algorithm Usage (CLI Demo)

This example shows how to use the built-in Shor demo from the CLI. The simulator
implements a small-scale, educational version of Shor's algorithm and is intended
for small composite numbers.

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
