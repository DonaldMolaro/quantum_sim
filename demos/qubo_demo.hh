#pragma once

#include <vector>

void run_qubo_demo();
void run_qubo_exact_cli(int n, const std::vector<double>& matrix);
void run_qubo_grover_cli(int n, double threshold, int iterations, const std::vector<double>& matrix);
