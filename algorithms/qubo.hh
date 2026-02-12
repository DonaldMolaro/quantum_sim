#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct QuboExactResult {
  bool ok = false;
  std::string error;
  Bitstring argmin = 0ULL;
  double min_value = 0.0;
  std::vector<Bitstring> minimizers;
};

struct QuboGroverResult {
  bool ok = false;
  std::string error;
  bool found = false;
  Bitstring candidate = 0ULL;
  double candidate_value = 0.0;
  int iterations = 0;
  double expected_success = 0.0;
};

bool qubo_matrix_valid(int n, const std::vector<double>& q, std::string& error);
double qubo_evaluate(Bitstring assignment, int n, const std::vector<double>& q);
QuboExactResult qubo_solve_exact(int n, const std::vector<double>& q);
QuboGroverResult qubo_solve_grover_threshold(int n,
                                             const std::vector<double>& q,
                                             double threshold,
                                             int iterations = -1);
