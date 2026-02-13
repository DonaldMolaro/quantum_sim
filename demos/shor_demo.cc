#include "demos/shor_demo.hh"
#include "algorithms/shor_classical.hh"
#include "algorithms/shor_quantum.hh"
#include "math/mod_arith.hh"
#include <cstdlib>
#include <ctime>
#include <iostream>

static bool run_shor_algorithm(Bitstring N, std::ostream& out)
{
  static bool seeded = false;
  if (const char* env = std::getenv("QSIM_RNG_SEED")) {
    unsigned int seed = static_cast<unsigned int>(std::strtoul(env, nullptr, 10));
    std::srand(seed);
  } else if (!seeded) {
    seeded = true;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
  }
  if (N < 2) {
    out << "N must be >= 2.\n";
    return false;
  }
  if (N % 2 == 0) {
    out << "Trivial factor found: 2 x " << (N / 2) << "\n";
    return true;
  }

  // Pick random a in [2, N-2], unless overridden for tests.
  Bitstring a = 2 + (std::rand() % (N - 3));
  if (const char* env_a = std::getenv("QSIM_SHOR_FORCE_A")) {
    a = static_cast<Bitstring>(std::strtoull(env_a, nullptr, 10));
  }
  Bitstring g = gcd_bitstring(a, N);
  if (g > 1 && g < N) {
    out << "Lucky gcd found: " << g << " x " << (N / g) << "\n";
    return true;
  }

  int n_c = 0;
  Bitstring measured_x = 0ULL;
  if (const char* env_x = std::getenv("QSIM_SHOR_FORCE_X")) {
    if (const char* env_nc = std::getenv("QSIM_SHOR_FORCE_NC")) {
      measured_x = static_cast<Bitstring>(std::strtoull(env_x, nullptr, 10));
      n_c = static_cast<int>(std::strtol(env_nc, nullptr, 10));
    } else {
      ShorQuantumResult q = run_shor_algorithm_quantum_part(N, a);
      if (!q.ok) {
        out << q.error << "\n";
        if (q.max_n > 0) {
          out << "Try N <= " << q.max_n << " for this simulator.\n";
        }
        return false;
      }
      measured_x = q.measured_x;
      n_c = q.n_c;
    }
  } else {
    ShorQuantumResult q = run_shor_algorithm_quantum_part(N, a);
    if (!q.ok) {
      out << q.error << "\n";
      if (q.max_n > 0) {
        out << "Try N <= " << q.max_n << " for this simulator.\n";
      }
      return false;
    }
    measured_x = q.measured_x;
    n_c = q.n_c;
  }
  if (measured_x == 0 && n_c == 0) {
    return false;
  }

  out << "--- Classical Post-Processing ---\n";
  out << "Result x/2^n = " << measured_x << " / " << (1ULL << n_c) << "\n";

  Bitstring r = 0;
  if (const char* env_r = std::getenv("QSIM_SHOR_FORCE_R")) {
    r = static_cast<Bitstring>(std::strtoull(env_r, nullptr, 10));
  } else {
    r = estimate_order(measured_x, n_c, a, N);
  }
  if (r == 0) {
    out << "Failed to estimate order r from continued fractions.\n";
    return false;
  }

  out << "Estimated order r = " << r << "\n";
  if (r % 2 != 0) {
    out << "Order r is odd; try again.\n";
    return false;
  }

  Bitstring ar2 = mod_pow(a, r / 2, N);
  if (ar2 == N - 1) {
    out << "a^(r/2) ≡ -1 (mod N); try again.\n";
    return false;
  }

  Bitstring p = gcd_bitstring(ar2 + 1, N);
  Bitstring q = gcd_bitstring((ar2 == 0 ? N : ar2 - 1), N);
  if (p == 1 || q == 1 || p == N || q == N) {
    out << "Failed to extract non-trivial factors; try again.\n";
    return false;
  }

  out << "Non-trivial factors found: " << p << " x " << q << "\n";
  return true;
}

void run_shor_demo(Bitstring N)
{
  std::ostream& out = std::cout;
  int max_attempts = 5;
  if (const char* env = std::getenv("QSIM_SHOR_MAX_ATTEMPTS")) {
    int v = std::atoi(env);
    if (v >= 0) {
      max_attempts = v;
    }
  }
  // Try a few attempts to account for probabilistic outcomes.
  for (int attempt = 0; attempt < max_attempts; ++attempt) {
    out << "\n=== Shor Attempt " << (attempt + 1) << " ===\n";
    if (run_shor_algorithm(N, out)) {
      return;
    }
  }
  out << "Shor demo did not find factors after multiple attempts.\n";
}
