#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct BernsteinVaziraniResult {
  bool ok = false;
  std::string error;
  Bitstring measured_secret = 0ULL;
  std::vector<int> measured_bits;
};

BernsteinVaziraniResult run_bernstein_vazirani(int n_inputs, Bitstring secret, int bias = 0);
