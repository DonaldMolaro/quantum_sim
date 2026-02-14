# Annealing Usage (QUBO)

This example shows the pragmatic annealing workflow for QUBO in the simulator.

Optional verbose logs:

```
VERBOSE VERBOSE
```

## Classical Computing View

Annealing here is a stochastic local-search family for binary optimization.

- SA: one replica, temperature schedule controls hill-climbing vs exploration.
- SQA-style: multiple coupled replicas to improve basin traversal.
- Baseline comparison: exact brute force for small `n`.

Use this when you want scalable approximate solutions and a familiar
meta-heuristic workflow.

## 1) Built-in anneal demo

```
ANNEAL DEMO
```

Runs SQA on the built-in 3-variable QUBO and compares to the exact minimum.

## 2) Run simulated annealing (SA)

Syntax:

```
ANNEAL QUBO SA <n> <steps> <sweeps> <beta_start> <beta_end> <replicas> <n*n matrix entries>
```

Example:

```
ANNEAL QUBO SA 3 80 20 0.1 6.0 1 -2 0 2 0 1 0 2 0 -3
```

## 3) Run simulated quantum annealing (SQA-style replicas)

Syntax:

```
ANNEAL QUBO SQA <n> <steps> <sweeps> <beta_start> <beta_end> <replicas> <n*n matrix entries>
```

Example:

```
ANNEAL QUBO SQA 3 80 20 0.1 6.0 8 -2 0 2 0 1 0 2 0 -3
```

Note: `replicas` is meaningful for `SQA` and should be at least `2`.
