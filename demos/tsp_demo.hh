#pragma once

#include <vector>

void run_tsp_demo();
void run_tsp_exact_cli(int n_cities,
                       double penalty,
                       const std::vector<double>& distance);
