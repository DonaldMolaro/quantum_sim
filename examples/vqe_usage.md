# VQE Usage Examples

This file shows two illustrative Variational Quantum Eigensolver (VQE) examples.

## Classical Computing View

VQE is an optimization loop around a costly objective evaluation.

- Objective: minimize expected energy of a Hamiltonian.
- Classical loop:
  1. propose parameters
  2. evaluate objective (here via quantum-state simulation)
  3. update parameters
- This is comparable to derivative-free numerical optimization where objective
  calls are expensive and noisy.

For teaching/prototyping, VQE is best viewed as a hybrid optimizer with
physics-structured objective functions.

## Command format

```
VQE RUN <n> <layers> <iters> <step> <shots> <terms> <coeff op_count [op qubit]...>
```

Notes:
- `shots` is currently exact-only in this simulator, so use `0`.
- Each term is encoded as:
  - `<coeff> <op_count> [op qubit]...`
  - Example for `0.5 * X0 X1`: `0.5 2 X 0 X 1`

## 1) Single-qubit ground state (H = Z0)

Hamiltonian:
- `H = Z0`
- Ground state: `|1>`
- Ground energy: `-1`

Run:

```
VQE RUN 1 1 30 0.3 0 1 1.0 1 Z 0
```

You should see `best_energy` close to `-1`.

You can also run the built-in demo:

```
VQE DEMO
```

## 2) Two-qubit toy Hamiltonian (H = Z0 + Z1 + 0.5 X0 X1)

Hamiltonian:
- `H = Z0 + Z1 + 0.5 X0 X1`
- This is a compact example where entangling structure matters.

Run:

```
VQE RUN 2 2 50 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

Try adjusting:
- `layers` (e.g. `1`, `2`, `3`)
- `iters`
- `step`

to see how optimization quality changes.

## 3) Hand-tooled VQE workflow

A practical manual loop:

1. Start with a very small ansatz:

```
VQE RUN 2 1 20 0.3 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

2. Increase optimization budget:

```
VQE RUN 2 1 80 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

3. Increase ansatz expressivity:

```
VQE RUN 2 2 80 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

4. Compare the final `best_energy` values. Lower is better.

This is the recommended hand-tuned path before automating sweeps.
