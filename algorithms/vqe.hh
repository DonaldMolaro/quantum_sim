#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct VqePauliOp {
  char op = 'I'; // Supported: X, Y, Z
  int qubit = 0;
  VqePauliOp() {}
  VqePauliOp(char op_in, int qubit_in) : op(op_in), qubit(qubit_in) {}
};

struct VqePauliTerm {
  double coeff = 0.0;
  std::vector<VqePauliOp> ops;
};

struct VqeHamiltonian {
  int n_qubits = 0;
  std::vector<VqePauliTerm> terms;
};

struct VqeOptions {
  int layers = 1;
  int max_iters = 40;
  double step_size = 0.25;
  int shots = 0; // 0 = exact expectation (deterministic)
  unsigned seed = 7u;
};

struct VqeResult {
  bool ok = false;
  std::string error;
  double best_energy = 0.0;
  std::vector<double> best_params;
  std::vector<double> energy_history;
  int evaluations = 0;
};

bool vqe_hamiltonian_valid(const VqeHamiltonian& h, std::string& error);
double vqe_expectation_exact(const State& state, const VqeHamiltonian& h);
VqeResult run_vqe(const VqeHamiltonian& h, const VqeOptions& options);
