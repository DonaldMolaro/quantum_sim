#include "demos/vqe_demo.hh"

#include <iostream>

void run_vqe_cli(const VqeHamiltonian& h,
                 int layers,
                 int iters,
                 double step,
                 int shots)
{
  VqeOptions options;
  options.layers = layers;
  options.max_iters = iters;
  options.step_size = step;
  options.shots = shots;

  VqeResult result = run_vqe(h, options);
  if (!result.ok) {
    std::cerr << "VQE error: " << result.error << "\n";
    return;
  }

  std::cout << "VQE result:\n";
  std::cout << "  best_energy=" << result.best_energy << "\n";
  std::cout << "  evaluations=" << result.evaluations << "\n";
  std::cout << "  params=[";
  for (size_t i = 0; i < result.best_params.size(); ++i) {
    if (i != 0) std::cout << ", ";
    std::cout << result.best_params[i];
  }
  std::cout << "]\n";
}

void run_vqe_demo()
{
  // Simple one-qubit ground-state problem: minimize H = -Z.
  VqeHamiltonian h;
  h.n_qubits = 1;
  VqePauliTerm term;
  term.coeff = -1.0;
  term.ops.push_back(VqePauliOp{'Z', 0});
  h.terms.push_back(term);

  std::cout << "VQE demo: minimize H = -Z0\n";
  run_vqe_cli(h, 1, 20, 0.25, 0);
}
