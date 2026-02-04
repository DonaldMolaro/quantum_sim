#include "state.hh"
#include <cmath>
#include <complex>
#include <iostream>

static double probability_of_state(const State& s, Bitstring target)
{
  double prob = 0.0;
  const QuantumState& qs = s.get_state();
  for (const auto& pair : qs) {
    if (pair.first == target) {
      prob += std::norm(pair.second);
    }
  }
  return prob;
}

static bool run_grover_case(int n_qubits, Bitstring target, double min_prob)
{
  State s(n_qubits, 0);
  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }

  const double N = std::ldexp(1.0, n_qubits);
  const double PI = std::acos(-1.0);
  const int R = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N)));

  for (int i = 0; i < R; ++i) {
    s.grover_oracle_Uf(target);
    s.grover_diffusion_Us();
  }

  double prob = probability_of_state(s, target);
  std::cout << "Grover test n=" << n_qubits
            << " target=" << target
            << " R=" << R
            << " prob=" << prob << "\n";

  return prob >= min_prob;
}

int main()
{
  bool ok = true;

  // For N=4 (n=2), one Grover iteration should concentrate near 1.0.
  ok = ok && run_grover_case(2, 3, 0.9);

  // For N=8 (n=3), two iterations should noticeably amplify the target.
  ok = ok && run_grover_case(3, 5, 0.5);

  if (!ok) {
    std::cerr << "Grover tests failed.\n";
    return 1;
  }

  std::cout << "Grover tests passed.\n";
  return 0;
}
