#include "algorithms/qrng.hh"
#include "state.hh"
#include <cstdlib>

std::vector<int> qrng_bits(int n, const std::vector<double>* random_vals)
{
  if (n <= 0) {
    return {};
  }

  State s(n, n);
  s.set_basis_state(0, 1.0);

  for (int j = 0; j < n; ++j) {
    s.h(j);
  }

  std::vector<int> out;
  if (random_vals) {
    s.measure_all_with_rng(*random_vals, out);
  } else {
    s.measure_all(out);
  }

  return out;
}

uint64_t qrng_u64(int n, const std::vector<double>* random_vals)
{
  if (n <= 0) {
    return 0;
  }
  if (n > 63) {
    n = 63;
  }
  if (const char* env = std::getenv("QSIM_QRNG_MAX_QUBITS")) {
    int cap = std::atoi(env);
    if (cap > 0 && n > cap) {
      n = cap;
    }
  }

  std::vector<int> bits = qrng_bits(n, random_vals);
  uint64_t value = 0;
  for (size_t i = 0; i < bits.size(); ++i) {
    if (bits[i]) {
      value |= (1ULL << i);
    }
  }
  return value;
}
