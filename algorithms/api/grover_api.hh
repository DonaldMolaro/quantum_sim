#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct GroverResult {
  bool ok = false;
  std::string error;
  int iterations = 0;
};

GroverResult run_grover(State& s, const std::vector<Bitstring>& targets, int iterations = -1);
