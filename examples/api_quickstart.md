# C++ API Quickstart

This quickstart is aimed at engineers using the simulator as a C++ library.

## 1) Build library artifacts

```bash
cmake -S . -B build
cmake --build build -j
```

Library outputs are in `build/` (for example `libquantum_sim.a`).

## 2) Minimal gate usage

```cpp
#include "include/quantum_sim.hh"
#include <iostream>

int main() {
  State s(2, 0);
  s.h(0).cx(0, 1); // Bell pair
  s.display();
  return 0;
}
```

Compile/link (example):

```bash
c++ -std=c++11 -I. my_program.cc build/libquantum_sim.a -o my_program
```

## 3) Grover API (library path)

```cpp
#include "algorithms/api/grover_api.hh"
#include <vector>

int main() {
  std::vector<Bitstring> targets{3, 5};
  GroverResult r = run_grover(3, targets, -1);
  if (!r.ok) return 1;
  return 0;
}
```

## 4) Quantum counting + auto-tuned Grover

```cpp
#include "algorithms/quantum_counting.hh"
#include "algorithms/api/grover_api.hh"
#include <vector>

int main() {
  std::vector<Bitstring> targets{3, 5};
  QuantumCountingResult c = run_quantum_counting(3, targets, 6);
  if (!c.ok) return 1;

  GroverResult g = run_grover_auto_tuned(3, targets, 6);
  if (!g.ok) return 1;
  return 0;
}
```

## 5) Variational API sketch

`VQA/QAOA` and `VQE` are available through algorithm headers in `algorithms/`.
Use the CLI examples as parameter templates, then map those fields into
`VqaQaoaOptions`, `VqeHamiltonian`, and related types.

## 6) Reproducibility

- For deterministic behavior in demos/CLI: set `SEED <n>`.
- For test-like runs: set deterministic RNG inputs where APIs allow it.
