#pragma once
#include "state.hh"
#include <string>

struct ShorQuantumResult {
  bool ok = false;
  Bitstring measured_x = 0ULL;
  int n_c = 0;
  int n_t = 0;
  int total_qubits = 0;
  Bitstring max_n = 0ULL;
  std::string error;
};

ShorQuantumResult run_shor_algorithm_quantum_part(Bitstring N, Bitstring a);
