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
