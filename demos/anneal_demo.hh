#pragma once

#include <string>
#include <vector>

void run_anneal_demo();
void run_anneal_qubo_cli(const std::string& method_token,
                         int n,
                         int steps,
                         int sweeps,
                         double beta_start,
                         double beta_end,
                         int replicas,
                         const std::vector<double>& matrix);
