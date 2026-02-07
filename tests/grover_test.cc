#include "state.hh"
#include "algorithms/api/grover_api.hh"
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

static bool run_grover_api_case(int n_qubits,
                                const std::vector<Bitstring>& targets,
                                int iterations,
                                double min_prob)
{
  State s(n_qubits, 0);
  GroverResult result = run_grover(s, targets, iterations);
  if (!result.ok) {
    std::cerr << "Grover API error: " << result.error << "\n";
    return false;
  }

  double min_found = 1.0;
  for (Bitstring t : targets) {
    double prob = probability_of_state(s, t);
    if (prob < min_found) {
      min_found = prob;
    }
  }

  std::cout << "Grover API test n=" << n_qubits
            << " targets=" << targets.size()
            << " R=" << result.iterations
            << " min_prob=" << min_found
            << " expected_success=" << result.expected_success << "\n";

  return min_found >= min_prob;
}

static bool run_grover_api_duplicates_case()
{
  State s(3, 0);
  std::vector<Bitstring> targets = {1, 1, 6};
  GroverResult result = run_grover(s, targets, -1);
  if (!result.ok) {
    std::cerr << "Grover API error: " << result.error << "\n";
    return false;
  }

  double p1 = probability_of_state(s, 1);
  double p6 = probability_of_state(s, 6);
  double min_found = std::min(p1, p6);

  std::cout << "Grover API duplicates test min_prob=" << min_found
            << " expected_success=" << result.expected_success << "\n";

  return min_found >= 0.2;
}

static bool run_grover_api_zero_iterations_case()
{
  const int n_qubits = 2;
  const double N = std::ldexp(1.0, n_qubits);
  State s(n_qubits, 0);
  std::vector<Bitstring> targets = {0, 1, 2}; // M=3 => R_default=0
  GroverResult result = run_grover(s, targets, -1);
  if (!result.ok) {
    std::cerr << "Grover API error: " << result.error << "\n";
    return false;
  }
  if (result.iterations != 0) {
    std::cerr << "Expected R=0, got " << result.iterations << "\n";
    return false;
  }

  double expected = 1.0 / N;
  for (Bitstring t : targets) {
    double prob = probability_of_state(s, t);
    if (std::abs(prob - expected) > 1e-6) {
      std::cerr << "Expected uniform probability for R=0.\n";
      return false;
    }
  }

  return true;
}

static bool run_grover_api_out_of_range_case()
{
  State s(3, 0);
  std::vector<Bitstring> targets = {8}; // out of range for n=3 (max 7)
  GroverResult result = run_grover(s, targets, -1);
  if (result.ok) {
    std::cerr << "Expected out-of-range error.\n";
    return false;
  }
  return true;
}

static bool run_grover_api_stress_case()
{
  const int n_qubits = 5; // N=32
  std::vector<Bitstring> targets = {3, 7, 12, 19}; // M=4
  return run_grover_api_case(n_qubits, targets, -1, 0.2);
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

  // Grover API single-target (explicit iterations).
  ok = ok && run_grover_api_case(3, {5}, 2, 0.5);

  // Grover API duplicate targets should de-duplicate without hurting amplification.
  ok = ok && run_grover_api_duplicates_case();

  // Grover API when R=0 should leave uniform distribution.
  ok = ok && run_grover_api_zero_iterations_case();

  // Grover API out-of-range target should fail.
  ok = ok && run_grover_api_out_of_range_case();

  // Grover API stress case (n=5, M=4).
  ok = ok && run_grover_api_stress_case();

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
