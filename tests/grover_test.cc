#include "state.hh"
#include <cmath>
#include <complex>
#include <iostream>
#include <unordered_set>
#include "algorithms/latin_square.hh"

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

static bool run_grover_multi_case(int n_qubits,
                                  const std::vector<Bitstring>& targets,
                                  double min_prob)
{
  State s(n_qubits, 0);
  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }

  const double N = std::ldexp(1.0, n_qubits);
  const double PI = std::acos(-1.0);
  const double M = static_cast<double>(targets.size());
  const int R = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N / M)));

  std::unordered_set<Bitstring> target_set;
  for (Bitstring t : targets) {
    target_set.insert(t);
  }

  for (int i = 0; i < R; ++i) {
    s.grover_oracle_Uf_multi(target_set);
    s.grover_diffusion_Us();
  }

  double min_found = 1.0;
  for (Bitstring t : targets) {
    double prob = probability_of_state(s, t);
    if (prob < min_found) {
      min_found = prob;
    }
  }

  std::cout << "Grover multi-test n=" << n_qubits
            << " targets=" << targets.size()
            << " R=" << R
            << " min_prob=" << min_found << "\n";

  return min_found >= min_prob;
}

int main()
{
  bool ok = true;

  // For N=4 (n=2), one Grover iteration should concentrate near 1.0.
  ok = ok && run_grover_case(2, 3, 0.9);

  // For N=8 (n=3), two iterations should noticeably amplify the target.
  ok = ok && run_grover_case(3, 5, 0.5);

  // Multi-solution case: N=8, targets {1,6} should both be amplified.
  ok = ok && run_grover_multi_case(3, {1, 6}, 0.2);

  // Latin-3 fixed-row has exactly 2 solutions.
  int latin_count = 0;
  for (Bitstring b = 0; b < (1ULL << 12); ++b) {
    if (is_valid_latin3_fixedrow(b)) {
      ++latin_count;
    }
  }
  std::cout << "Latin-3 fixed-row solution count: " << latin_count << "\n";
  ok = ok && (latin_count == 2);

  if (!ok) {
    std::cerr << "Grover tests failed.\n";
    return 1;
  }

  std::cout << "Grover tests passed.\n";
  return 0;
}
