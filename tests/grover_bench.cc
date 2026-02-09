#include "algorithms/api/grover_api.hh"
#include "state.hh"
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

static bool sample_hit(int n_qubits,
                       const std::vector<Bitstring>& targets,
                       int iterations,
                       std::mt19937& rng)
{
  State s(n_qubits, n_qubits);
  GroverResult result = run_grover(s, targets, iterations);
  if (!result.ok) {
    return false;
  }

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<double> rands(n_qubits);
  for (int i = 0; i < n_qubits; ++i) {
    rands[i] = dist(rng);
  }

  std::vector<int> out;
  s.measure_all_with_rng(rands, out);
  Bitstring measured = 0;
  for (int i = 0; i < n_qubits; ++i) {
    if (out[i]) {
      measured |= (1ULL << i);
    }
  }

  for (Bitstring t : targets) {
    if (t == measured) {
      return true;
    }
  }
  return false;
}

static void run_case(int n_qubits, const std::vector<Bitstring>& targets, int samples)
{
  State s(n_qubits, 0);
  GroverResult result = run_grover(s, targets, -1);
  if (!result.ok) {
    std::cerr << "Grover error: " << result.error << "\n";
    return;
  }

  std::mt19937 rng(12345);
  int hits = 0;
  for (int i = 0; i < samples; ++i) {
    if (sample_hit(n_qubits, targets, result.iterations, rng)) {
      ++hits;
    }
  }

  double empirical = static_cast<double>(hits) / static_cast<double>(samples);
  std::cout << "Grover bench n=" << n_qubits
            << " M=" << targets.size()
            << " R=" << result.iterations
            << " expected=" << result.expected_success
            << " empirical=" << empirical
            << "\n";
}

int run_grover_bench()
{
  auto start = std::chrono::steady_clock::now();

  run_case(5, {3}, 200);
  run_case(6, {5, 27}, 200);
  run_case(7, {7, 19, 88, 101}, 200);

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Grover bench elapsed: " << elapsed.count() << "s\n";
  return 0;
}

#ifndef ALL_TESTS
int main()
{
  return run_grover_bench();
}
#endif
