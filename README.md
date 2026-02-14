# quantum_sim
Quantum Computing Simulator in C++

This is a quantum computing simulator based on the paper "Quantum Computing without the Linear Algebra" by Aws Albarghouthi.

Aws has an implementation in Python, but this one is in C++.

Much of the code was created by NotebookLM, with subsequent additions and refactors by Codex (OpenAI).

## Build (CMake)

```bash
cmake -S . -B build
cmake --build build -j
```

## Tests (CMake)

```bash
ctest --test-dir build --output-on-failure
ctest --test-dir build -R test-slow --output-on-failure
ctest --test-dir build -R test-demo --output-on-failure

# Example smoke/lab checks
cmake --build build --target test-examples
cmake --build build --target test-labs
```

## Coverage Gate (CMake)

```bash
cmake -S . -B build-cov -DQSIM_ENABLE_COVERAGE=ON
cmake --build build-cov -j
cmake --build build-cov --target coverage-check
```

## Project Layout

- `state.hh`, `state.cc`: core quantum state and gates
- `display.cc`: state visualization
- `math/`: bit operations and modular arithmetic helpers
- `algorithms/`: core Grover/Shor routines (no demo output)
- `algorithms/shor_classical.*`: continued fractions + order extraction helpers
- `algorithms/shor_quantum.*`: quantum order-finding routines
- `demos/`: CLI/demo wrappers for Shor and Latin square output
- `tests/`: unit tests and Grover tests
- `cli/`: interactive CLI driver program
- `examples/`: worked CLI transcripts (see `examples/README.md` for the index)

## Style Notes

- Include order: project headers first, then standard library headers.

## Conventions

- Bit ordering uses little-endian indexing: qubit 0 is the least-significant bit.
- Registers are laid out as contiguous ranges; for Shor, target comes first, then control.
- Register ranges are represented by `QubitRange` (`start`, `end`, `size()`), used in `RegisterLayout`.
- QFT has two equivalent simulator paths selected with `QFTMODE`:
  - `DIRECT` (default): direct matrix-style transform, usually faster here.
  - `GATE`: explicit gate decomposition with controlled phase rotations.

## Gate Index (CLI)

Single-qubit:
`H`, `X`, `Y`, `Z`, `S`, `T`, `RX`, `RY`, `RZ`, `RU`

Two-qubit (controlled):
`CX` (`CNOT`), `CZ`, `CY`, `CH`, `CRX`, `CRY`, `CRZ`, `CU`

Three-qubit:
`CCX` (`TOFFOLI`)

Other:
`SWAP`
