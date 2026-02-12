#pragma once

#include <vector>

void run_vqa_demo();
void run_vqa_qaoa_cli(int n,
                      int p_layers,
                      int shots,
                      int max_iters,
                      double step_size,
                      const std::vector<double>& matrix);
