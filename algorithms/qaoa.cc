#include "algorithms/qaoa.hh"

#include "algorithms/vqa_qaoa.hh"

QaoaResult run_qaoa_qubo(int n,
                         const std::vector<double>& q,
                         const QaoaOptions& options)
{
  VqaQaoaOptions legacy_options;
  legacy_options.p_layers = options.p_layers;
  legacy_options.max_iters = options.max_iters;
  legacy_options.shots = options.shots;
  legacy_options.step_size = options.step_size;
  legacy_options.seed = options.seed;

  VqaQaoaResult legacy = run_vqa_qaoa_qubo(n, q, legacy_options);

  QaoaResult out;
  out.ok = legacy.ok;
  out.error = legacy.error;
  out.best_gamma = legacy.best_gamma;
  out.best_beta = legacy.best_beta;
  out.best_bitstring = legacy.best_bitstring;
  out.best_energy = legacy.best_energy;
  out.energy_history = legacy.energy_history;
  out.evaluations = legacy.evaluations;
  return out;
}
