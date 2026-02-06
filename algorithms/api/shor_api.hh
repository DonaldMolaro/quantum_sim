#pragma once

#include "state.hh"
#include <string>

struct ShorResult {
  bool ok = false;
  std::string error;
  Bitstring measured_x = 0;
  int n_c = 0;
};

ShorResult run_shor_quantum_part(Bitstring N, Bitstring a);
