#pragma once

#include "state.hh"
#include <vector>

void run_grover_search(State *s, Bitstring target_w);
void run_grover_search_multi(State *s, const std::vector<Bitstring>& targets);
