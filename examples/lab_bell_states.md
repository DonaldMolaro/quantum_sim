# Guided Lab: Bell States

This lab walks through preparing and verifying all four Bell states.

## Goal

Prepare `|Phi+>`, `|Phi->`, `|Psi+>`, `|Psi->` and verify the expected basis
support/probabilities.

## Setup

```bash
./quantum_sim
```

Enable guided narration:

```
TUTOR ON
```

If you want to validate your results after working through the lab:

```bash
./examples/check_labs.sh bell
```

## Part A: Prepare `|Phi+> = (|00> + |11>)/sqrt(2)`

Predict first:
- Which two basis states should survive?
- Should the amplitudes have the same sign or opposite sign?

```
INIT 2
H 0
CX 0 1
DISPLAY
CHECK BELL PHI+
```

Checkpoint:
- Non-zero amplitudes should be at `|00>` and `|11>`.
- `CHECK BELL PHI+` should report `PASS`.

## Part B: Prepare `|Phi-> = (|00> - |11>)/sqrt(2)`

Predict first:
- Which basis states stay non-zero?
- What changed relative to `|Phi+>`: support or phase?

```
INIT 2
H 0
CX 0 1
Z 0
DISPLAY
CHECK BELL PHI-
```

Checkpoint:
- Same basis support as `|Phi+>`, but opposite phase on one component.
- `CHECK BELL PHI-` should report `PASS`.

## Part C: Prepare `|Psi+> = (|01> + |10>)/sqrt(2)`

Predict first:
- Which computational basis states should now carry the amplitude?

```
INIT 2
H 0
CX 0 1
X 1
DISPLAY
CHECK BELL PSI+
```

Checkpoint:
- Non-zero amplitudes should be at `|01>` and `|10>`.
- `CHECK BELL PSI+` should report `PASS`.

## Part D: Prepare `|Psi-> = (|01> - |10>)/sqrt(2)`

Predict first:
- Does `Z 0` change basis support here, or only the relative sign?

```
INIT 2
H 0
CX 0 1
X 1
Z 0
DISPLAY
CHECK BELL PSI-
```

Checkpoint:
- Same support as `|Psi+>`, with relative sign flip.
- `CHECK BELL PSI-` should report `PASS`.

## Extension

Measure both qubits repeatedly and compare observed correlations for each Bell state.
