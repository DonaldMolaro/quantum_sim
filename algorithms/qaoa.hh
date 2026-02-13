#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct QaoaOptions {
  int p_layers = 1;
  int max_iters = 40;
  int shots = 0; // 0 = exact expectation (deterministic)
  double step_size = 0.25;
  unsigned seed = 7u;
};

struct QaoaResult {
  bool ok = false;
  std::string error;
  std::vector<double> best_gamma;
  std::vector<double> best_beta;
  Bitstring best_bitstring = 0ULL;
  double best_energy = 0.0;
  std::vector<double> energy_history;
  int evaluations = 0;
};

QaoaResult run_qaoa_qubo(int n,
                         const std::vector<double>& q,
                         const QaoaOptions& options);
