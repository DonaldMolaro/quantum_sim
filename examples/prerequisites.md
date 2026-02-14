# Prerequisites Map

Use this as a learning dependency graph before running each algorithm.

## Foundation

Required first:
- `INIT`, `DISPLAY`, `MEASURE`
- Single-qubit gates: `H`, `X`, `Z`, `RZ`
- Controlled gate intuition: `CX`, `CZ`
- Bit ordering rule in this simulator: `q0` is least-significant bit

## Oracle Algorithms

Before `DEUTSCH_JOZSA`, `BV`, `SIMON`:
- Phase kickback idea (`X` ancilla + `H`)
- Boolean function/oracle modeling
- For Simon: linear algebra over GF(2)

## Search and Counting

Before `GROVER`, `QCOUNT`:
- Oracle marking as phase flip
- Diffusion operator intuition (reflection about average)
- Relationship between target count `M` and iterations

## Factoring

Before `SHOR`:
- Modular arithmetic basics (`a^r mod N`)
- GCD and coprime checks
- Continued-fraction/order interpretation

## Variational Optimization

Before `QAOA`, `VQE`:
- Objective minimization and local search behavior
- Hyperparameter tuning (`layers`, `iters`, `step`, `shots`)
- Compare against exact baseline when possible (`QUBO EXACT`)
