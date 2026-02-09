# quantum_sim
Quantum Computing Simulator in C++

This is a quantum computing simulator based on the paper "Quantum Computing without the Linear Algebra" by Aws Albarghouthi.

Aws has an implementation in Python, but this one is in C++.

Much of the code was created by NotebookLM, with subsequent additions and refactors by Codex (OpenAI).

## Build

```bash
make
```

## Tests

```bash
make test          # fast suite + 100% coverage (library)
make test-slow     # includes slow Shor coverage paths
./all_tests        # direct runner (no coverage gate)

# Optional Grover benchmark sweep
QSIM_GROVER_BENCH=1 ./all_tests
```

Notes:
- `make test` and `make test-slow` enforce 100% coverage using `gcov`.
- `make clean` removes build artifacts and `gcov` files.

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

## Style Notes

- Include order: project headers first, then standard library headers.

## Conventions

- Bit ordering uses little-endian indexing: qubit 0 is the least-significant bit.
- Registers are laid out as contiguous ranges; for Shor, target comes first, then control.
- Register ranges are represented by `QubitRange` (`start`, `end`, `size()`), used in `RegisterLayout`.

## Gate Index (CLI)

Single-qubit:
`H`, `X`, `Y`, `Z`, `S`, `T`, `RX`, `RY`, `RZ`, `RU`

Two-qubit (controlled):
`CX` (`CNOT`), `CZ`, `CY`, `CH`, `CRX`, `CRY`, `CRZ`, `CU`

Three-qubit:
`CCX` (`TOFFOLI`)

Other:
`SWAP`
