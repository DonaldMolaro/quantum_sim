# Lesson Plan Index

This index maps teaching topics to simulator commands, expected outcomes, and
recommended labs.

## Module 1: Core State Manipulation

- Topic: superposition and phase
- Commands: `INIT`, `H`, `X`, `Z`, `DISPLAY`, `MEASURE`
- Expected output: amplitude redistribution and phase-only changes
- Lab: `lab_bell_states.md`

## Module 2: Entanglement

- Topic: Bell-state preparation and verification
- Commands: `H`, `CX`, `Z`, `DISPLAY`
- Expected output: two-basis-state support with correlated measurements
- Lab: `lab_bell_states.md`

## Module 3: Oracle Algorithms

- Topic: constant vs balanced / hidden string / hidden xor mask
- Commands: `DEUTSCH_JOZSA`, `BV`, `SIMON`
- Expected output: algorithm-specific classification/recovery output
- Lab pairing: `deutsch_jozsa.md`, `bernstein_vazirani.md`, `simon_usage.md`

## Module 4: Search and Counting

- Topic: amplitude amplification and marked-count estimation
- Commands: `GROVER`, `QCOUNT`, `GROVER AUTO`
- Expected output: iteration counts and amplified target probability
- Lab: `lab_grover_counting.md`

## Module 5: Factoring (Educational Scale)

- Topic: period-finding workflow and classical post-processing
- Commands: `SHOR`
- Expected output: non-trivial factors on small composites
- Lab: `lab_shor.md`

## Module 6: Variational Methods

- Topic: objective minimization with parameterized circuits
- Commands: `QAOA QUBO`, `VQE RUN`, `QUBO EXACT`
- Expected output: approximate minimum energy plus exact reference comparison
- Labs: `lab_qaoa.md`, `lab_vqe.md`

## Teaching notes

- Use `SEED <n>` at the beginning of each session for reproducibility.
- Use `TUTOR ON` for step-wise narration during live demos.
- Use `make test-examples` and `./examples/check_labs.sh` to verify scripted
  command paths before class.
