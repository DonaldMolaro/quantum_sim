#pragma once

#include "state.hh"
#include <vector>

void run_quantum_counting_cli(int n_qubits, const std::vector<Bitstring>& targets, int max_iterations = -1);
void run_quantum_counting_demo();
