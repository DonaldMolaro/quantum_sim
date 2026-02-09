# Quantum Random Number Generator (QRNG) Demo

This example uses the `QRNG` command in the CLI to generate random bits using
Hadamard-superposed qubits and measurement.

## 1) Generate a single random bit

```
QRNG 1
```

Expected: prints a single random bit.

## 2) Generate multiple random bits (count)

```
QRNG 1 16
```

Expected: prints 16 random bits.

## 3) Generate small random integers

With `n` qubits, `QRNG n` samples a number in [0, 2^n - 1].

```
QRNG 4 8
```

Expected: prints 8 random integers between 0 and 15.

## 4) Compare against deterministic seed (tests)

The test suite uses a fixed seed for reproducible output. In the interactive CLI,
no seed is forced, so results are non-deterministic.

## Notes

- `QRNG n count` uses `n` qubits and repeats `count` times.
- You can use `DISPLAY` after manual preparation to inspect amplitude states,
  but `QRNG` handles setup internally.
