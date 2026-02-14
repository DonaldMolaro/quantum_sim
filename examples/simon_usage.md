# Simon's Algorithm Usage (CLI Demo)

Simon solves a hidden-xor-mask problem: find `s` such that
`f(x) = f(x xor s)` for all `x`.

## Commands

Run built-in demo:

```
SIMON DEMO
```

Run with explicit parameters:

```
SIMON <n_inputs> <secret> [shots]
```

Example:

```
SIMON 4 10 8
```

This asks the simulator to recover secret `10` (`1010b`) using 8 sampled
equations.

## Notes

- `n_inputs` must be in `[1,16]`.
- `secret` must be in `[0, 2^n - 1]`.
- `shots` defaults to an internal value when omitted.

## Hand-tooled Simon (equation method)

The simulator command automates equation sampling and solving. A hand-tooled
equivalent is to collect linear equations `y·s = 0 (mod 2)` and solve for `s`.

Example for `n=4`, secret `s=1010`:

- `y1=0101` gives `0* s3 + 1* s2 + 0* s1 + 1* s0 = 0`
- `y2=1010` gives `1* s3 + 0* s2 + 1* s1 + 0* s0 = 0`
- `y3=1111` gives `s3 + s2 + s1 + s0 = 0`

One non-zero solution satisfying all equations is `s=1010`.

Use CLI output for comparison:

```
SIMON 4 10 8
```

This demonstrates the same recovery path as the manual GF(2) solve.
