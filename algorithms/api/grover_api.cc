#include "algorithms/api/grover_api.hh"
#include <cmath>
#include <unordered_set>

GroverResult run_grover(State& s, const std::vector<Bitstring>& targets, int iterations)
{
  GroverResult result{0};
  if (targets.empty()) {
    return result;
  }

  int n_qubits = s.get_num_qubits();
  if (n_qubits >= 63) {
    return result;
  }

  const double N = std::ldexp(1.0, n_qubits);
  const double PI = std::acos(-1.0);
  const double M = static_cast<double>(targets.size());
  const int R_default = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N / M)));
  const int R = (iterations >= 0) ? iterations : R_default;

  std::unordered_set<Bitstring> target_set;
  for (Bitstring t : targets) {
    target_set.insert(t);
  }

  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }

  for (int k = 0; k < R; ++k) {
    s.grover_oracle_Uf_multi(target_set);
    s.grover_diffusion_Us();
  }

  result.iterations = R;
  return result;
}
