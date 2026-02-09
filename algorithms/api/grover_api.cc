#include "algorithms/api/grover_api.hh"
#include <cmath>
#include <unordered_set>

GroverResult run_grover(State& s, const std::vector<Bitstring>& targets, int iterations)
{
  GroverResult result;
  if (targets.empty()) {
    result.error = "No targets provided.";
    return result;
  }

  int n_qubits = s.get_num_qubits();
  if (n_qubits >= 63) {
    result.error = "Grover supports up to 62 qubits.";
    return result;
  }

  const double N = std::ldexp(1.0, n_qubits);
  const double PI = std::acos(-1.0);

  std::unordered_set<Bitstring> target_set;
  for (Bitstring t : targets) {
    if (t >= static_cast<Bitstring>(N)) {
      result.error = "Target out of range for current number of qubits.";
      return result;
    }
    target_set.insert(t);
  }

  const double M = static_cast<double>(target_set.size());
  if (M == 0.0) {
    result.error = "No valid targets provided.";
    return result;
  }

  const int R_default = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N / M)));
  const int R = (iterations >= 0) ? iterations : R_default;

  const double theta = std::asin(std::sqrt(M / N));
  result.expected_success = std::sin((2.0 * R + 1.0) * theta);
  result.expected_success *= result.expected_success;

  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }

  if (R > 0) {
    if (N <= 1048576.0) {
      std::vector<uint8_t> mask(static_cast<size_t>(N), 0);
      for (Bitstring t : target_set) {
        mask[static_cast<size_t>(t)] = 1;
      }
      for (int k = 0; k < R; ++k) {
        s.grover_oracle_Uf_mask(mask);
        s.grover_diffusion_Us();
      }
    } else {
      for (int k = 0; k < R; ++k) {
        s.grover_oracle_Uf_multi(target_set);
        s.grover_diffusion_Us();
      }
    }
  }

  result.iterations = R;
  result.ok = true;
  return result;
}
