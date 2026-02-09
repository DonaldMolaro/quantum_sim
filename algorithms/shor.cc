/*
 * Shor's algorithm implementation (order-finding + classical post-processing).
 * This is a pedagogical simulator implementation that directly applies
 * unitary permutations for modular exponentiation (rather than decomposing
 * into elementary reversible gates).
 */
#include "state.hh"
#include "algorithms/shor_classical.hh"
#include "algorithms/shor_quantum.hh"
#include "math/mod_arith.hh"
#include <cstdlib>
#include <ctime>
#include <iostream>

static bool run_shor_algorithm(Bitstring N)
{
  static bool seeded = false;
  if (!seeded) {
    seeded = true;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
  }
  if (N < 2) {
    std::cout << "N must be >= 2.\n";
    return false;
  }
  if (N % 2 == 0) {
    std::cout << "Trivial factor found: 2 x " << (N / 2) << "\n";
    return true;
  }

  // Pick random a in [2, N-2], unless overridden for tests.
  Bitstring a = 2 + (std::rand() % (N - 3));
  if (const char* env_a = std::getenv("QSIM_SHOR_FORCE_A")) {
    a = static_cast<Bitstring>(std::strtoull(env_a, nullptr, 10));
  }
  Bitstring g = gcd_bitstring(a, N);
  if (g > 1 && g < N) {
    std::cout << "Lucky gcd found: " << g << " x " << (N / g) << "\n";
    return true;
  }

  int n_c = 0;
  Bitstring measured_x = 0ULL;
  if (const char* env_x = std::getenv("QSIM_SHOR_FORCE_X")) {
    if (const char* env_nc = std::getenv("QSIM_SHOR_FORCE_NC")) {
      measured_x = static_cast<Bitstring>(std::strtoull(env_x, nullptr, 10));
      n_c = static_cast<int>(std::strtol(env_nc, nullptr, 10));
    } else {
      measured_x = run_shor_algorithm_quantum_part(N, a, n_c);
    }
  } else {
    measured_x = run_shor_algorithm_quantum_part(N, a, n_c);
  }
  if (measured_x == 0 && n_c == 0) {
    return false;
  }

  std::cout << "--- Classical Post-Processing ---\n";
  std::cout << "Result x/2^n = " << measured_x << " / " << (1ULL << n_c) << "\n";

  Bitstring r = 0;
  if (const char* env_r = std::getenv("QSIM_SHOR_FORCE_R")) {
    r = static_cast<Bitstring>(std::strtoull(env_r, nullptr, 10));
  } else {
    r = estimate_order(measured_x, n_c, a, N);
  }
  if (r == 0) {
    std::cout << "Failed to estimate order r from continued fractions.\n";
    return false;
  }

  std::cout << "Estimated order r = " << r << "\n";
  if (r % 2 != 0) {
    std::cout << "Order r is odd; try again.\n";
    return false;
  }

  Bitstring ar2 = mod_pow(a, r / 2, N);
  if (ar2 == N - 1) {
    std::cout << "a^(r/2) ≡ -1 (mod N); try again.\n";
    return false;
  }

  Bitstring p = gcd_bitstring(ar2 + 1, N);
  Bitstring q = gcd_bitstring((ar2 == 0 ? N : ar2 - 1), N);
  if (p == 1 || q == 1 || p == N || q == N) {
    std::cout << "Failed to extract non-trivial factors; try again.\n";
    return false;
  }

  std::cout << "Non-trivial factors found: " << p << " x " << q << "\n";
  return true;
}

void run_shor_demo(Bitstring N)
{
  int max_attempts = 5;
  if (const char* env = std::getenv("QSIM_SHOR_MAX_ATTEMPTS")) {
    int v = std::atoi(env);
    if (v >= 0) {
      max_attempts = v;
    }
  }
  // Try a few attempts to account for probabilistic outcomes.
  for (int attempt = 0; attempt < max_attempts; ++attempt) {
    std::cout << "\n=== Shor Attempt " << (attempt + 1) << " ===\n";
    if (run_shor_algorithm(N)) {
      return;
    }
  }
  std::cout << "Shor demo did not find factors after multiple attempts.\n";
}
