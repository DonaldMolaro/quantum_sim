# quantum_sim White Paper
## System Overview, Usage, and Self-Training Path

## 1. Executive Summary

`quantum_sim` is a quantum computing simulator written in C++ with two goals:

- provide an understandable implementation of core quantum ideas and algorithms
- provide a practical training environment for engineers who already think in classical systems terms

The project is based on the paper *Quantum Computing without the Linear Algebra* by Aws Albarghouthi, then extended into a larger educational and experimental system. It supports direct state inspection, an interactive CLI, multiple algorithm demonstrations, guided examples, and a substantial self-study course aimed at classical engineers.

This white paper explains:

- what the simulator is and how it is structured
- how to build and use it
- how to train yourself with the included course and examples
- what the system is good at educationally
- where its practical limits are

## 2. Problem Statement

Many introductions to quantum computing fail engineers in one of two ways:

- they stay too abstract, leaving learners with vocabulary but little operational intuition
- they hide too much behind large frameworks, so the learner never sees what the system is doing

`quantum_sim` takes the opposite approach. It makes the internal state inspectable, keeps the command surface direct, and lets the learner move from small gates to full algorithms while staying close to the underlying mechanics.

The intended outcome is not only "I ran Grover" or "I saw Shor factor 15." The intended outcome is: "I understand what changed in the state, why it changed, and how the command sequence maps to the algorithm."

## 3. System Goals

### 3.1 Educational goals

- make superposition, interference, entanglement, and measurement visible
- let users verify results with simple built-in checks
- support progressive learning from one qubit to full algorithms
- give classical engineers a path from intuition to quantum reasoning

### 3.2 Engineering goals

- keep the codebase practical to build and run locally
- provide both Make and CMake build paths
- maintain test coverage and example/lab validation
- keep CLI behavior consistent with course and examples

## 4. What the System Includes

The simulator combines several layers:

- a core quantum state engine
- a CLI for interactive exploration
- algorithm implementations and demos
- examples, labs, and a full self-study course
- automated tests, example smoke checks, and lab validators

At a high level, this is not just a library and not just a tutorial set. It is a small educational platform built around a quantum state simulator.

## 5. Architecture Overview

### 5.1 Core state engine

The state engine is centered on `State`, implemented in [state.hh](state.hh) and [state.cc](state.cc).

It is responsible for:

- sparse quantum state representation
- gate application
- measurement and classical-register updates
- analysis helpers like Bloch vectors, expectation values, and entropy
- state-display support

Related visualization logic lives in [display.cc](display.cc).

### 5.2 CLI layer

The CLI lives in [cli/](cli/) and provides:

- direct gate entry
- algorithm commands
- layered help
- tutor mode
- file-based workflows with `LOAD` and `SAVE`
- multi-shot measurement workflows
- terminal-native visual tools such as `DISPLAY TOP <k>`, `DISPLAY PROBS`, `BLOCH`, and `CIRCUIT`

This layer is intentionally central to the learning experience. The project is designed to be explored interactively.

### 5.3 Algorithms and demos

Core algorithms live in [algorithms/](algorithms/), with demo wrappers in [demos/](demos/).

The supported material includes:

- Deutsch-Jozsa
- Bernstein-Vazirani
- Grover
- Quantum counting
- Simon
- QPE
- Shor
- QAOA / VQA
- VQE
- annealing-style QUBO search
- simple QEC

### 5.4 Training material

The training material lives in [examples/](examples/).

The key entry points are:

- [examples/course.md](examples/course.md): the full self-study course
- [examples/README.md](examples/README.md): examples index
- guided labs such as Bell, Grover, Shor, QAOA, and VQE
- helper scripts such as `run_examples.sh` and `check_labs.sh`

## 6. Build and Run

### 6.1 Requirements

At minimum, you need:

- a C++ compiler with C++17 support
- `make` for the Make-based build path
- optionally `cmake` for the CMake path

### 6.2 Build with Make

```bash
make
./quantum_sim
```

### 6.3 Build with CMake

```bash
cmake -S . -B build
cmake --build build
./build/quantum_sim
```

## 7. How to Use the Simulator

The normal learning loop is:

1. initialize a state
2. apply gates or run an algorithm
3. inspect the result
4. validate your interpretation

Example:

```text
INIT 2 2
H 0
CX 0 1
DISPLAY
CHECK BELL PHI+
```

Useful inspection and learning commands include:

- `DISPLAY`
- `DISPLAY TOP <k>`
- `DISPLAY PROBS`
- `BLOCH <q>`
- `EXPECT ...`
- `ENTROPY ...`
- `CHECK ...`
- `SHOTS <n> <file>`
- `CIRCUIT`
- `TUTOR ON`
- `VERBOSE VERBOSE`

### 7.1 Terminal-native visualization

The current visualization approach is intentionally terminal-first. Rather than introducing a separate GUI, the simulator exposes information through compact textual views:

- amplitude tables
- probability-only views
- histogram bars
- tutor-mode state deltas
- Bloch summaries
- ASCII circuit views from command history

This keeps the tool portable and lowers the cost of experimenting.

## 8. How to Train Yourself with the Course

The recommended self-training path is:

### Stage 1: First contact

Start with [examples/course.md](examples/course.md), especially the "First 60 Minutes" path near the top. This gives you:

- one-qubit superposition
- phase/interference
- entanglement
- one complete algorithm workflow

This is the right starting point if you want quick momentum before committing to the full course.

### Stage 2: Full course

Then work through the full module sequence in [examples/course.md](examples/course.md).

The course is built for classical engineers. It starts with what a qubit is, then builds through:

- phase and interference
- Bell states and entanglement
- oracle-based algorithms
- period finding and factoring
- optimization and variational methods
- error correction

### Stage 3: Guided labs

Once you are comfortable with the basic flow, move into the lab material:

- [examples/lab_bell_states.md](examples/lab_bell_states.md)
- [examples/lab_grover_counting.md](examples/lab_grover_counting.md)
- [examples/lab_shor.md](examples/lab_shor.md)
- [examples/lab_qaoa.md](examples/lab_qaoa.md)
- [examples/lab_vqe.md](examples/lab_vqe.md)

These labs are structured around prediction, execution, and validation.

### Stage 4: Self-check and automation

Use the built-in scripts to close the loop:

```bash
./examples/run_examples.sh
./examples/check_labs.sh
./examples/check_labs.sh bell
./examples/check_labs.sh grover
```

This gives the learner the kind of feedback loop engineers usually expect from tests or exercises.

## 9. Recommended Learning Workflow

For each topic:

1. read the short section in the course
2. run the command block manually
3. predict the outcome before inspecting it
4. inspect with `DISPLAY`, `DISPLAY PROBS`, `BLOCH`, or `CIRCUIT`
5. validate with `CHECK`, `SHOTS`, or the lab script

That pattern matters. The simulator is strongest when used as an interactive reasoning tool, not just as a demo launcher.

## 10. Why This Works Well for Classical Engineers

The system is especially well-suited to engineers who already think in terms of:

- state transitions
- representations
- probability distributions
- debugging and observability
- experiments with expected outcomes
- verification loops

The simulator maps quantum ideas onto those familiar habits:

- `DISPLAY` exposes state
- `CHECK` acts like assertions
- `SHOTS` exposes repeated-run behavior
- `CIRCUIT` connects commands to structure
- `TUTOR ON` explains intent and effect

That combination reduces the usual abstraction shock.

## 11. Verification and Quality Controls

The repo includes several quality mechanisms:

- unit tests
- full coverage gate
- example smoke tests
- guided lab validators
- both Make and CMake build paths

Typical verification commands:

```bash
make test
./examples/run_examples.sh
./examples/check_labs.sh
```

And with CMake:

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## 12. Limits and Scope

This is an educational simulator, not a production quantum SDK and not a high-scale distributed simulator.

Its strengths are:

- clarity
- inspectability
- local experimentation
- algorithm intuition

Its limits are the usual simulator limits:

- state growth is exponential in qubit count
- some workflows are illustrative rather than hardware-realistic
- the tool is optimized for understanding, not for scaling to very large systems

That tradeoff is deliberate.

## 13. Best Practices for Users

- Start with the course, not with Shor.
- Turn on `TUTOR ON` early.
- Use `DISPLAY PROBS` when you care about outcomes more than amplitudes.
- Use `DISPLAY TOP <k>` when the state is getting busy.
- Use `CIRCUIT` after a command sequence to connect behavior to circuit structure.
- Use `SHOTS` when you want the probability distribution, not just one measurement.
- Use the lab checks as progress gates.

## 14. Recommended Next Improvements

The current visualization system is already useful, but the next high-value steps are clear:

- richer tutor-mode deltas for algorithms, not just gates
- tighter circuit rendering for more command types
- step-indexed history playback
- optional export of circuit/state snapshots for static diagrams or HTML views

Those would extend the same terminal-native teaching model without changing the system’s educational character.

## 15. Conclusion

`quantum_sim` is most effective when understood as a training system wrapped around a simulator.

It gives engineers:

- a buildable local environment
- a visible quantum state model
- practical CLI exploration tools
- a structured self-study course
- examples and labs with validation

If your goal is to understand what quantum programs are doing, not just invoke them, this is the right way to use the repository:

start in the course, inspect everything, validate often, and use the simulator as an instrument rather than a black box.
