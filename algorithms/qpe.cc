#include "algorithms/qpe.hh"
#include "state.hh"
#include "internal/limits.hh"
#include <cmath>
#include <stdexcept>

QpeResult run_qpe(int m, double phase_radians)
{
  if (m < 1 || m > 20) {
    throw std::invalid_argument("QPE: m must be in [1, 20]");
  }

  // Total qubits: 1 target (qubit 0) + m precision (qubits 1..m).
  const int n_qubits = m + 1;
  const int target = 0;

  State s(n_qubits);

  // 1. Prepare eigenstate |1> on target qubit.
  s.x(target);

  // 2. Superpose precision register.
  for (int k = 1; k <= m; k++) {
    s.h(k);
  }

  // 3. Controlled-U^{2^{k-1}}: CP(k, target, phase * 2^{k-1}).
  for (int k = 1; k <= m; k++) {
    double angle = phase_radians * static_cast<double>(1LL << (k - 1));
    s.cp(k, target, angle);
  }

  // 4. Inverse QFT on precision register (qubits 1..m).
  s.iqft(1, m);

  // 5. Find most-probable bitstring in precision register.
  const QuantumState& qs = s.get_state();
  Bitstring best_bs = 0;
  double best_prob = -1.0;
  for (const auto& pair : qs) {
    double prob = std::norm(pair.second);
    if (prob > best_prob) {
      best_prob = prob;
      best_bs = pair.first;
    }
  }

  // Extract precision bits (little-endian): bits 1..m of the bitstring.
  int measured_int = 0;
  for (int k = 0; k < m; k++) {
    measured_int |= (static_cast<int>((best_bs >> (k + 1)) & 1) << k);
  }

  const double two_pow_m = static_cast<double>(1LL << m);
  QpeResult result;
  result.m = m;
  result.true_phase_frac = phase_radians / (2.0 * qsim::limits::PI);
  result.measured_int = measured_int;
  result.estimated_frac = measured_int / two_pow_m;
  result.estimated_radians = result.estimated_frac * 2.0 * qsim::limits::PI;
  result.best_prob = best_prob;
  return result;
}
