# Deutsch–Jozsa Algorithm (CLI Demo)

This example uses the built-in Deutsch–Jozsa demo to distinguish between
constant and balanced oracles with a single query.
Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

## Classical Computing View

Think of this as a black-box function classification problem.

- Input: a function `f(x)` promised to be either:
  - `CONSTANT`: same output for all inputs
  - `BALANCED`: output `1` for exactly half the inputs
- Classical deterministic strategy: in the worst case, you need `2^(n-1)+1`
  queries to be certain.
- Deutsch-Jozsa strategy: uses interference so one oracle query is enough to
  distinguish the two promised classes.

If you are coming from classical testing/property checking, this is like a
promise-problem accelerator: we are not learning the full truth table; we only
decide which class the function belongs to.

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

---

## By Hand (n = 2)

Below is a manual circuit for n=2 (inputs q0, q1; ancilla q2). We show one
constant and one balanced oracle so you can see the difference without
calling the built-in demo.

### A) Constant oracle f(x) = 0 (no-op)

```
INIT 3
X 2
H 0
H 1
H 2
# oracle: f(x)=0, do nothing
H 0
H 1
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

Expected: measured input register = `00` (constant).

### B) Balanced oracle f(x) = x0 (flip ancilla when q0=1)

```
INIT 3
X 2
H 0
H 1
H 2
# oracle: f(x)=x0
CX 0 2
H 0
H 1
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

Expected: measured input register has at least one 1 (balanced).
