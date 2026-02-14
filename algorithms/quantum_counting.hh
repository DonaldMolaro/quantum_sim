#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct QuantumCountingResult {
  bool ok = false;
  std::string error;
  int n_qubits = 0;
  Bitstring state_size = 0ULL;
  Bitstring estimated_targets = 0ULL;
  double estimated_targets_real = 0.0;
  double estimated_theta = 0.0;
  int iterations_used = 0;
  std::vector<double> success_probabilities;
};

QuantumCountingResult run_quantum_counting(int n_qubits,
                                           const std::vector<Bitstring>& targets,
                                           int max_iterations = -1);
