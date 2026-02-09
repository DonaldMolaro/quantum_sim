# BB84 Quantum Key Distribution (Sketch with CLI)

This is a **conceptual walkthrough** using the `quantum_sim` CLI to demonstrate
BB84. The simulator does not model a noisy channel or explicit eavesdropping,
but we can show the prepare/measure steps and how sifting works.
Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

We use one qubit at a time. Alice chooses a random bit and random basis, prepares
that qubit, sends it to Bob, and Bob measures in a random basis.

## Basis and Encoding

- Z-basis (computational):
  - bit 0 => |0>
  - bit 1 => |1>
- X-basis (Hadamard):
  - bit 0 => |+> = H|0>
  - bit 1 => |-> = H|1> = H X |0>

## Alice’s Preparation Commands

Assume fresh qubit `q0` for each round:

- Basis Z, bit 0:
  ```
  INIT 1
  # |0>
  ```

- Basis Z, bit 1:
  ```
  INIT 1
  X 0
  ```

- Basis X, bit 0:
  ```
  INIT 1
  H 0
  ```

- Basis X, bit 1:
  ```
  INIT 1
  X 0
  H 0
  ```

## Bob’s Measurement Commands

Bob chooses a measurement basis:

- Z-basis measurement:
  ```
  MEASURE 0 0
  ```

- X-basis measurement (apply H first, then measure Z):
  ```
  H 0
  MEASURE 0 0
  ```

## One Example Round (Same Basis)

Alice chooses basis X, bit 1; Bob chooses basis X:

```
INIT 1
X 0
H 0
H 0
MEASURE 0 0
DISPLAY
```
Expected: Bob recovers bit 1 with certainty.

## One Example Round (Different Bases)

Alice chooses basis X, bit 1; Bob chooses basis Z:

```
INIT 1
X 0
H 0
MEASURE 0 0
DISPLAY
```
Expected: Bob’s result is random (0 or 1 with 50/50 probability).

## Sifting (Classical Step)

After many rounds:
1. Alice and Bob publicly compare bases (not bits).
2. Keep only rounds where bases match.
3. The remaining bits form the raw key.

## Notes

- This example is single-qubit-at-a-time for clarity.
- To simulate an eavesdropper (Eve), insert a measurement in a random basis
  before Bob’s measurement; mismatched bases introduce errors.
