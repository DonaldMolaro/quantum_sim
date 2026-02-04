#pragma once

#include "state.hh"
#include <vector>

struct GroverResult {
  int iterations;
};

GroverResult run_grover(State& s, const std::vector<Bitstring>& targets, int iterations = -1);
