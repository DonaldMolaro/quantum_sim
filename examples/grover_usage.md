# Grover's Algorithm Usage (CLI Demo)

This example shows how to use the built-in Grover demo and how to build a small
Grover oracle by hand in the CLI.
Optional: enable step-by-step logs while you run commands:

```
VERBOSE VERBOSE
```

## 1) Run the Grover demo

Start the CLI:

```
./quantum_sim
```

Then run a demo search. Example: search for targets 3 and 5 with 3 qubits:

```
GROVER 3 5
```

The demo prints the marked items and the amplified probabilities.

If target count is unknown, auto-tune iterations with counting first:

```
GROVER AUTO 3 6 3 5
```

Format: `GROVER AUTO <n_qubits> <count_iterations> <targets...>`.

---

## 2) Build a 2-qubit Grover oracle by hand (target = |11>)

We will amplify the state |11> out of 4 possibilities.

### Step A: Initialize and create uniform superposition

```
INIT 2
H 0
H 1
DISPLAY
```

### Step B: Oracle for |11>

To phase-flip |11>, we can use a controlled-Z between the qubits:

```
CZ 0 1
```

### Step C: Diffusion operator (2 qubits)

The diffusion operator is:

1. Apply H to all qubits
2. Apply X to all qubits
3. Apply CZ (phase flip |11>)
4. Apply X to all qubits
5. Apply H to all qubits

Commands:

```
H 0
H 1
X 0
X 1
CZ 0 1
X 0
X 1
H 0
H 1
```

### Step D: Inspect result

```
DISPLAY
```

Expected: amplitude for |11> is significantly amplified.

---

## 3) Repeat for a different target (example: |10>)

To mark |10>, we can flip qubit 0 before and after the CZ so that |10> maps to |11>:

```
INIT 2
H 0
H 1
X 0
CZ 0 1
X 0
# diffusion
H 0
H 1
X 0
X 1
CZ 0 1
X 0
X 1
H 0
H 1
DISPLAY
```

---

## Notes

- For 2 qubits, one Grover iteration is enough to maximize success.
- The built-in `GROVER` command handles all of this automatically for larger N.
- This manual approach is useful for learning or verifying oracle behavior.
