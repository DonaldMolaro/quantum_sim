# Glossary (Classical-First)

This glossary maps common quantum terms to familiar classical concepts.

- `Basis state`:
  A concrete bit pattern (like a single classical register value), e.g. `|101>`.
- `Amplitude`:
  Complex weight attached to a basis state. Probability is `|amplitude|^2`.
- `Superposition`:
  One state object representing weighted many-state possibilities at once.
- `Phase`:
  Sign/angle component of amplitude; does not directly change probability, but
  changes interference outcomes.
- `Interference`:
  Amplitudes add/cancel, analogous to constructive/destructive signal addition.
- `Oracle`:
  A black-box predicate/function used by an algorithm (e.g., Grover marking).
- `Entanglement`:
  Joint state cannot be factored into independent per-qubit states.
- `Measurement`:
  Sampling step that returns classical bits and collapses state.
- `QFT`:
  Fourier transform over amplitudes; similar role to DFT/FFT in signal analysis.
- `Hamiltonian`:
  Objective/operator defining energy in VQE.
- `Ansatz`:
  Parameterized circuit family used as search space in variational algorithms.
- `QUBO`:
  Binary quadratic objective `x^T Q x`, widely used for combinatorial optimization.
