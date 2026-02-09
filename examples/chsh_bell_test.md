# CHSH / Bell Inequality Test (Sketch with CLI)

This example demonstrates the CHSH test using the simulator. We prepare an entangled
Bell pair and measure correlations along different bases. With ideal measurements,
quantum correlations violate the classical CHSH bound of 2.

We use two qubits:
- `q0`: Alice
- `q1`: Bob

## 1) Prepare |Φ+> = (|00> + |11>)/sqrt(2)

```
INIT 2
H 0
CNOT 0 1
```

## 2) Measurement Settings

CHSH uses two settings per side:
- Alice: A0, A1
- Bob:   B0, B1

For a simple idealized CHSH test, we can approximate rotations to set measurement
bases. A standard choice uses:
- Alice: A0 = Z, A1 = X
- Bob:   B0 = (Z + X)/sqrt(2), B1 = (Z - X)/sqrt(2)

We implement rotated measurement bases by applying a rotation before Z measurement.
A convenient choice:
- Measuring in basis (Z+X)/sqrt(2): apply `RY -45DEG` then measure Z
- Measuring in basis (Z-X)/sqrt(2): apply `RY 45DEG` then measure Z

## 3) Example Correlation Runs

For each setting pair (A?, B?), do the following:

### A0, B0
Alice measures Z, Bob measures (Z+X)/sqrt(2):
```
INIT 2
H 0
CNOT 0 1
# Alice A0: Z basis
# Bob B0: rotate then Z measure
RY 1 -45DEG
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

### A0, B1
```
INIT 2
H 0
CNOT 0 1
RY 1 45DEG
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

### A1, B0
Alice measures X, Bob measures (Z+X)/sqrt(2):
```
INIT 2
H 0
CNOT 0 1
H 0
RY 1 -45DEG
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

### A1, B1
```
INIT 2
H 0
CNOT 0 1
H 0
RY 1 45DEG
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

## 4) Interpreting Results

Repeat each setting many times, collect correlation values E(Ai, Bj), and compute:

S = E(A0,B0) + E(A0,B1) + E(A1,B0) - E(A1,B1)

Classical bound: |S| <= 2
Quantum ideal:   |S| can reach 2*sqrt(2) ≈ 2.828

## Notes

- This is a *conceptual* CHSH sketch. To get stable statistics, repeat each
  setting many times and compute averages externally.
- You can add a small script to automate the repeated trials if desired.
