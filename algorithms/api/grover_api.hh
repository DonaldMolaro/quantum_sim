#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct GroverResult {
  bool ok = false;
  std::string error;
  int iterations = 0;
  double expected_success = 0.0;
  bool auto_tuned = false;
  Bitstring estimated_targets = 0ULL;
  double estimated_targets_real = 0.0;
  int counting_iterations = 0;
};

GroverResult run_grover(State& s, const std::vector<Bitstring>& targets, int iterations = -1);
GroverResult run_grover(int n_qubits, const std::vector<Bitstring>& targets, int iterations = -1);
GroverResult run_grover_auto_tuned(int n_qubits,
                                   const std::vector<Bitstring>& targets,
                                   int counting_iterations = -1);
