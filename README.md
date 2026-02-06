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
make unit_tests
./unit_tests

make test   # runs grover_test
```

## Project Layout

- `state.hh`, `state.cc`: core quantum state and gates
- `display.cc`: state visualization
- `math/`: bit operations and modular arithmetic helpers
- `algorithms/`: Grover and Shor implementations
- `tests/`: unit tests and Grover tests
- `cli/`: interactive CLI driver program

## Style Notes

- Include order: project headers first, then standard library headers.

## Conventions

- Bit ordering uses little-endian indexing: qubit 0 is the least-significant bit.
- Registers are laid out as contiguous ranges; for Shor, target comes first, then control.
