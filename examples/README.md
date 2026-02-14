# Examples Index

This folder contains runnable CLI transcripts and walkthroughs for common quantum
protocols and algorithms supported by the simulator.

## Contents

- `bell_states.md`: Prepare and verify all four Bell states
- `superdense_coding.md`: Superdense coding (2 classical bits via 1 qubit)
- `quantum_teleportation.md`: Teleport a single-qubit state
- `bb84_qkd.md`: BB84 quantum key distribution (conceptual walkthrough)
- `chsh_bell_test.md`: CHSH/Bell inequality test sketch
- `qrng_demo.md`: Quantum random number generator demo
- `shor_usage.md`: Shor demo usage (small composites)
- `grover_usage.md`: Grover demo + manual oracle by hand
- `deutsch_jozsa.md`: Deutsch–Jozsa demo (constant vs balanced)
- `bernstein_vazirani.md`: Bernstein-Vazirani hidden-string recovery
- `qubo_usage.md`: QUBO exact and Grover-threshold solving
- `tsp_usage.md`: Traveling Salesman as QUBO (fixed-start encoding)
- `quantum_counting.md`: Estimate number of marked states (`QCOUNT`)
- `simon_usage.md`: Simon's hidden-mask recovery (`SIMON`)
- `algorithm_comparison.md`: quick algorithm-selection and workflow guide
- `vqa_qaoa.md`: VQA/QAOA optimization of QUBO instances
- `anneal_usage.md`: Simulated annealing and SQA-style QUBO optimization
- `max_cut.md`: Max-Cut encoded and solved as QUBO
- `vqe_usage.md`: VQE examples (`H=Z0` and a two-qubit toy Hamiltonian)

## Running the CLI

From the project root:

```
./quantum_sim
```

Then copy/paste the command blocks from any example.

Tip: set verbosity to see algorithm steps as they run:

```
VERBOSE VERBOSE
```

## Smoke-run examples

You can run a lightweight examples smoke test:

```bash
./examples/run_examples.sh
```
