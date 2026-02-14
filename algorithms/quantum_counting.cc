#include "algorithms/quantum_counting.hh"

#include "internal/limits.hh"
#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace {

double success_probability_for_iterations(int n_qubits,
                                          const std::unordered_set<Bitstring>& targets,
                                          int iterations)
{
  State s(n_qubits, 0);
  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }
  for (int k = 0; k < iterations; ++k) {
    s.grover_oracle_Uf_multi(targets);
    s.grover_diffusion_Us();
  }

  double p = 0.0;
  const QuantumState& qs = s.get_state();
  for (size_t i = 0; i < qs.size(); ++i) {
    if (targets.find(qs[i].first) != targets.end()) {
      p += std::norm(qs[i].second);
    }
  }
  return p;
}

double fit_theta_from_probabilities(const std::vector<double>& probs)
{
  const long double PI = std::acos(-1.0L);
  const int grid = 20000;
  double best_theta = 0.0;
  long double best_loss = 1e300L;

  for (int i = 0; i <= grid; ++i) {
    const long double theta = (PI / 2.0L) * static_cast<long double>(i) / static_cast<long double>(grid);
    long double loss = 0.0L;
    for (size_t k = 0; k < probs.size(); ++k) {
      const long double pred = std::sin((2.0L * static_cast<long double>(k) + 1.0L) * theta);
      const long double pk = pred * pred;
      const long double diff = pk - static_cast<long double>(probs[k]);
      loss += diff * diff;
    }
    if (loss < best_loss) {
      best_loss = loss;
      best_theta = static_cast<double>(theta);
    }
  }
  return best_theta;
}

} // namespace

QuantumCountingResult run_quantum_counting(int n_qubits,
                                           const std::vector<Bitstring>& targets,
                                           int max_iterations)
{
  QuantumCountingResult result;
  result.n_qubits = n_qubits;
  if (!qsim::limits::valid_bitstring_qubit_count(n_qubits)) {
    result.error = "n_qubits must be in [1, 62]";
    return result;
  }

  const Bitstring N = 1ULL << n_qubits;
  result.state_size = N;
  std::unordered_set<Bitstring> target_set;
  for (size_t i = 0; i < targets.size(); ++i) {
    if (targets[i] >= N) {
      result.error = "target out of range for given n_qubits";
      return result;
    }
    target_set.insert(targets[i]);
  }
  if (target_set.empty()) {
    result.error = "at least one target is required";
    return result;
  }

  int iters = max_iterations;
  if (iters < 0) {
    iters = std::min(16, std::max(2, n_qubits + 1));
  }
  if (iters < 1) {
    result.error = "max_iterations must be >= 1";
    return result;
  }

  result.success_probabilities.reserve(static_cast<size_t>(iters));
  for (int k = 0; k < iters; ++k) {
    result.success_probabilities.push_back(success_probability_for_iterations(n_qubits, target_set, k));
  }

  const double theta = fit_theta_from_probabilities(result.success_probabilities);
  const double m_real = static_cast<double>(N) * std::pow(std::sin(theta), 2.0);
  Bitstring m_round = static_cast<Bitstring>(std::llround(m_real));
  if (m_round > N) m_round = N;

  result.estimated_theta = theta;
  result.estimated_targets_real = m_real;
  result.estimated_targets = m_round;
  result.iterations_used = iters;
  result.ok = true;
  return result;
}
