# Deutsch–Jozsa Algorithm (CLI Demo)

This example uses the built-in Deutsch–Jozsa demo to distinguish between
constant and balanced oracles with a single query.

## 1) Start the CLI

```
./quantum_sim
```

## 2) Constant Oracles

```
DEUTSCH_JOZSA 3 CONST0
DEUTSCH_JOZSA 3 CONST1
```

Expected: the measured input register is all zeros, and the result reports CONSTANT.

## 3) Balanced Oracles

```
DEUTSCH_JOZSA 3 BALANCED_XOR0
DEUTSCH_JOZSA 3 BALANCED_PARITY
```

Expected: the measured input register contains at least one 1, and the result
reports BALANCED.

## Notes

- `n` is the number of input qubits; the demo allocates one extra ancilla.
- `BALANCED_XOR0` uses f(x) = x0 (first input bit).
- `BALANCED_PARITY` uses f(x) = parity of all input bits.
