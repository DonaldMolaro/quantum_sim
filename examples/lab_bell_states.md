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

## Part A: Prepare `|Phi+> = (|00> + |11>)/sqrt(2)`

```
INIT 2
H 0
CX 0 1
DISPLAY
```

Checkpoint:
- Non-zero amplitudes should be at `|00>` and `|11>`.

## Part B: Prepare `|Phi-> = (|00> - |11>)/sqrt(2)`

```
INIT 2
H 0
CX 0 1
Z 0
DISPLAY
```

Checkpoint:
- Same basis support as `|Phi+>`, but opposite phase on one component.

## Part C: Prepare `|Psi+> = (|01> + |10>)/sqrt(2)`

```
INIT 2
H 0
CX 0 1
X 1
DISPLAY
```

Checkpoint:
- Non-zero amplitudes should be at `|01>` and `|10>`.

## Part D: Prepare `|Psi-> = (|01> - |10>)/sqrt(2)`

```
INIT 2
H 0
CX 0 1
X 1
Z 0
DISPLAY
```

Checkpoint:
- Same support as `|Psi+>`, with relative sign flip.

## Extension

Measure both qubits repeatedly and compare observed correlations for each Bell state.
