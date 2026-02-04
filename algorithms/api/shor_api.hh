#pragma once

#include "state.hh"

struct ShorResult {
  Bitstring measured_x;
  int n_c;
};

ShorResult run_shor_quantum_part(Bitstring N, Bitstring a);
