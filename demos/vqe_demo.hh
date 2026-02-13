#pragma once

#include "algorithms/vqe.hh"

void run_vqe_demo();
void run_vqe_cli(const VqeHamiltonian& h,
                 int layers,
                 int iters,
                 double step,
                 int shots);
