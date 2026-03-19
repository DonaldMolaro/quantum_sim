# Quantum Bomb Detection (Hand-Tooled)

This is a hand-tooled Elitzur–Vaidman style walkthrough using the simulator as
an educational model.

## Classical Computing View

Classically, checking whether a bomb sensor is live usually means interacting
with it and risking detonation.

Quantum interaction-free measurement changes that tradeoff:
- Some runs still explode.
- Some runs prove a live bomb without explosion (`dark-port` click).

## Toy model used here

- One qubit `q0` models path:
  - `|0>`: upper path
  - `|1>`: lower path (bomb path)
- `H` acts like a beamsplitter.
- A live bomb on `|1>` can absorb the photon (explosion branch).
- A `dark-port` click is modeled as measuring `1` at output.

## Case A: No bomb (reference interference)

```text
INIT 1 1
H 0
H 0
MEASURE 0 0
DISPLAY
```

Expected:
- Always measure `0` (bright port only).
- Intuition: second beamsplitter recombines paths coherently.

## Case B: Live bomb, no-explosion branch (post-selected)

After first beamsplitter:
- state is `( |0> + |1> ) / sqrt(2)`
- if path is `|1>`, bomb explodes (discard run)
- conditioning on "no explosion" collapses to `|0>`

Hand-tool the surviving branch:

```text
INIT 1 1
H 0
# live bomb on |1| path
# condition on no explosion => surviving state behaves as |0>
INIT 1 1
H 0
MEASURE 0 0
DISPLAY
```

Expected in no-explosion branch:
- `P(0) ≈ 0.5`, `P(1) ≈ 0.5`
- Any observed `1` (dark-port click) certifies a live bomb without detonation.

## Case C: Dud bomb (no absorption)

```text
INIT 1 1
H 0
# dud bomb => no interaction
H 0
MEASURE 0 0
DISPLAY
```

Expected:
- Always `0`, same as no-bomb case.

## Why this demonstrates interaction-free detection

- No-bomb and dud-bomb preserve interference: dark port never clicks.
- Live-bomb possibility breaks interference via which-path disturbance.
- Therefore dark-port clicks become possible and indicate a live bomb, even
  though that specific run did not explode.
