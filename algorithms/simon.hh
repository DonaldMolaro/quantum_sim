#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct SimonResult {
  bool ok = false;
  std::string error;
  int n_inputs = 0;
  Bitstring secret = 0ULL;
  Bitstring recovered_secret = 0ULL;
  bool verified = false;
  std::vector<Bitstring> equations;
};

SimonResult run_simon(int n_inputs, Bitstring secret, int shots = -1);
