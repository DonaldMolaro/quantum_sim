# Superdense Coding (2 classical bits via 1 qubit)

This example uses the `quantum_sim` CLI to demonstrate superdense coding. Alice and Bob
share an entangled Bell pair. Alice encodes two classical bits by applying a local gate
to her qubit, then sends that qubit to Bob. Bob decodes with `CNOT` and `H` and measures.
Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

Qubits:
- `q0` = Alice
- `q1` = Bob

## Shared Bell Pair Setup (|Φ+>)

```
INIT 2
H 0
CNOT 0 1
```

At this point, Alice and Bob share |Φ+> = (|00> + |11>)/sqrt(2).

---

## Encoding Map (Alice on q0)

Alice encodes two classical bits by applying a gate to `q0`:

- `00` : I (do nothing)
- `01` : X
- `10` : Z
- `11` : X then Z (or Z then X, up to a global phase)

---

## Decoding (Bob)

After Alice applies her gate, Bob performs:

```
CNOT 0 1
H 0
```

Then Bob measures both qubits to recover the 2-bit message.

---

## Full CLI Transcripts (All Four Messages)

### Message = 00

```
INIT 2
H 0
CNOT 0 1
# Alice encodes 00 (I) -> no gate
CNOT 0 1
H 0
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: classical bits = 00.

### Message = 01

```
INIT 2
H 0
CNOT 0 1
X 0
CNOT 0 1
H 0
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: classical bits = 01.

### Message = 10

```
INIT 2
H 0
CNOT 0 1
Z 0
CNOT 0 1
H 0
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: classical bits = 10.

### Message = 11

```
INIT 2
H 0
CNOT 0 1
X 0
Z 0
CNOT 0 1
H 0
MEASURE 0 0
MEASURE 1 1
DISPLAY
```
Expected: classical bits = 11.

---

## Notes

- The protocol transmits **2 classical bits** using **1 qubit**, plus shared entanglement.
- Any global phase does not affect measurement outcomes.
- You can add `DISPLAY` after the encoding step to inspect the Bell basis states.
