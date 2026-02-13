#pragma once

#include "state.hh"
#include <string>
#include <vector>

enum class AnnealMethod {
  SA,
  SQA
};

struct AnnealOptions {
  AnnealMethod method = AnnealMethod::SQA;
  int steps = 80;
  int sweeps_per_step = 20;
  double beta_start = 0.10;
  double beta_end = 6.00;
  int replicas = 8; // SQA only
  unsigned seed = 7u;
};

struct AnnealResult {
  bool ok = false;
  std::string error;
  Bitstring best_assignment = 0ULL;
  double best_value = 0.0;
  std::vector<double> best_history;
  int accepted_moves = 0;
  int attempted_moves = 0;
};

bool anneal_options_valid(int n,
                          const std::vector<double>& q,
                          const AnnealOptions& options,
                          std::string& error);

AnnealResult anneal_qubo(int n,
                         const std::vector<double>& q,
                         const AnnealOptions& options);

const char* anneal_method_name(AnnealMethod method);
bool parse_anneal_method(const std::string& token, AnnealMethod& out);
