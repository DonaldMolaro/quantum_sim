# Bell States: Create and Verify

These examples use the `quantum_sim` CLI to prepare and verify all four Bell states.
Each state is shown with a short preparation circuit and a simple measurement-based
verification strategy. Run the commands in order inside the simulator shell.

## Common Setup

All examples start from a clean 2-qubit state:

```
INIT 2
```

`DISPLAY` shows the current amplitudes if you want to inspect the state vector.

---

## |Φ+> = (|00> + |11>) / sqrt(2)

**Create**
```
INIT 2
H 0
CNOT 0 1
DISPLAY
```

**Verify (correlated measurements)**
```
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: measurements are always the same (00 or 11).

---

## |Φ-> = (|00> - |11>) / sqrt(2)

**Create**
```
INIT 2
H 0
CNOT 0 1
Z 0
DISPLAY
```

**Verify (phase-sensitive check)**
Apply `H` to both qubits and measure. The relative phase flips which outcome is certain.

```
H 0
H 1
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: outcomes are always 01 or 10 (anti-correlated).

---

## |Ψ+> = (|01> + |10>) / sqrt(2)

**Create**
```
INIT 2
H 0
CNOT 0 1
X 1
DISPLAY
```

**Verify (anti-correlated measurements)**
```
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: measurements are always different (01 or 10).

---

## |Ψ-> = (|01> - |10>) / sqrt(2)

**Create**
```
INIT 2
H 0
CNOT 0 1
X 1
Z 0
DISPLAY
```

**Verify (phase-sensitive check)**
Apply `H` to both qubits and measure.

```
H 0
H 1
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: outcomes are always 00 or 11 (correlated).

---

## Notes

- The `H`/`CNOT` core prepares |Φ+>. Local `X` and `Z` gates transform to the other Bell states.
- Phase-sensitive verification uses a basis change (`H` on both qubits) to convert relative
  phase into measurable correlation changes.
