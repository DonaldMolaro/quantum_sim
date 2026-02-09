# Quantum Teleportation (State Transfer)

This example shows how to teleport a single-qubit state from Alice to Bob using
an entangled Bell pair and two classical bits. We use 3 qubits:
Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

- `q0`: Alice's *message* qubit (the state to teleport)
- `q1`: Alice's entangled qubit
- `q2`: Bob's entangled qubit (destination)

We will teleport a simple test state: `|ψ> = H|0> = (|0> + |1>)/sqrt(2)`.

## 1) Initialize and prepare |ψ> on q0

```
INIT 3
H 0
```

## 2) Create Bell pair between q1 and q2

```
H 1
CNOT 1 2
```

## 3) Bell-measure q0 and q1 (Alice)

```
CNOT 0 1
H 0
MEASURE 0 0
MEASURE 1 1
```

The two classical bits (c0, c1) are Alice's measurement results.

## 4) Bob applies corrections to q2

- If c1 == 1, apply X to q2
- If c0 == 1, apply Z to q2

Example correction sequence (apply both, then the invalid ones are no-ops if bit = 0):

```
X 2
Z 2
```

## 5) Verify (optional)

To verify, apply `H` to q2 and measure — it should match the original `|ψ>`
statistics (50/50 for 0 and 1).

```
H 2
MEASURE 2 2
DISPLAY
```

## Notes

- Teleportation transfers the **state**, not the physical qubit.
- The original state on q0 is destroyed by measurement.
- You can replace the initial `H 0` with any single-qubit preparation (e.g., `RY`, `RX`).
