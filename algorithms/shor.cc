/*
 * Shor's algorithm implementation (order-finding + classical post-processing).
 * This is a pedagogical simulator implementation that directly applies
 * unitary permutations for modular exponentiation (rather than decomposing
 * into elementary reversible gates).
 */
#include "state.hh"
#include "math/mod_arith.hh"
#include "math/register_layout.hh"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

static std::vector<Bitstring> continued_fraction(Bitstring numerator, Bitstring denominator)
{
  std::vector<Bitstring> terms;
  Bitstring n = numerator;
  Bitstring d = denominator;
  while (d != 0) {
    Bitstring a = n / d;
    terms.push_back(a);
    Bitstring r = n % d;
    n = d;
    d = r;
  }
  return terms;
}

static void convergents(const std::vector<Bitstring>& terms,
                        std::vector<Bitstring>& out_num,
                        std::vector<Bitstring>& out_den)
{
  out_num.clear();
  out_den.clear();
  Bitstring h1 = 1, h2 = 0;
  Bitstring k1 = 0, k2 = 1;
  for (size_t i = 0; i < terms.size(); ++i) {
    Bitstring a = terms[i];
    Bitstring h = a * h1 + h2;
    Bitstring k = a * k1 + k2;
    out_num.push_back(h);
    out_den.push_back(k);
    h2 = h1;
    h1 = h;
    k2 = k1;
    k1 = k;
  }
}

static Bitstring estimate_order(Bitstring measured_x, int n_c, Bitstring a, Bitstring N)
{
  Bitstring q = 1ULL << n_c;
  std::vector<Bitstring> terms = continued_fraction(measured_x, q);

  std::vector<Bitstring> nums;
  std::vector<Bitstring> dens;
  convergents(terms, nums, dens);

  for (size_t i = 0; i < dens.size(); ++i) {
    Bitstring r = dens[i];
    if (r == 0) continue;
    if (r >= N) continue;
    if (mod_pow(a, r, N) == 1ULL) {
      return r;
    }
  }
  return 0ULL;
}

static Bitstring run_shor_algorithm_quantum_part(Bitstring N, Bitstring a, int& out_n_c)
{
  int n_t = static_cast<int>(std::ceil(std::log2(static_cast<double>(N))));
  int n_c = 2 * n_t;
  int total_qubits = n_c + n_t;
  int num_cbits = n_c;

  out_n_c = n_c;
  if (n_c >= 63) {
    std::cout << "Control register too large for 64-bit demo (n_c=" << n_c << ").\n";
    return 0ULL;
  }

  State s(total_qubits, num_cbits);

  std::cout << "Running Shor's Algorithm for N=" << N << ", a=" << a << "\n";
  std::cout << "Total Qubits: " << total_qubits << " (Control=" << n_c << ", Target=" << n_t << ")\n";

  // Initialize target register to |1> (rightmost target qubit).
  s.x(total_qubits - 1);

  // QFT on control register (creates uniform superposition).
  s.qft(0, n_c - 1);

  // Controlled modular exponentiation for each control qubit.
  for (int j = 0; j < n_c; ++j) {
    int control_q = j;
    Bitstring power_of_a = 1ULL << j;
    s.controlled_modular_exponentiation(
        control_q,
        n_c, total_qubits - 1,
        a, N,
        power_of_a);
  }

  // Inverse QFT on control register.
  s.iqft(0, n_c - 1);

  // Measure control register into classical bits.
  std::cout << "\nMeasuring Control Register...\n";
  for (int j = 0; j < n_c; ++j) {
    s.measure(j, j);
  }

  Bitstring measured_x = 0;
  std::cout << "Measured result (bitstring x): ";
  for (int j = 0; j < n_c; ++j) {
    int bit = s.get_cbit(j);
    std::cout << bit;
    measured_x |= (Bitstring)bit << j;
  }
  std::cout << " (Decimal: " << measured_x << ")\n";

  return measured_x;
}

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

  // Pick random a in [2, N-2].
  Bitstring a = 2 + (std::rand() % (N - 3));
  Bitstring g = gcd_bitstring(a, N);
  if (g > 1 && g < N) {
    std::cout << "Lucky gcd found: " << g << " x " << (N / g) << "\n";
    return true;
  }

  int n_c = 0;
  Bitstring measured_x = run_shor_algorithm_quantum_part(N, a, n_c);
  if (measured_x == 0 && n_c == 0) {
    return false;
  }

  std::cout << "--- Classical Post-Processing ---\n";
  std::cout << "Result x/2^n = " << measured_x << " / " << (1ULL << n_c) << "\n";

  Bitstring r = estimate_order(measured_x, n_c, a, N);
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
  // Try a few attempts to account for probabilistic outcomes.
  for (int attempt = 0; attempt < 5; ++attempt) {
    std::cout << "\n=== Shor Attempt " << (attempt + 1) << " ===\n";
    if (run_shor_algorithm(N)) {
      return;
    }
  }
  std::cout << "Shor demo did not find factors after multiple attempts.\n";
}

State& State::run_shor_algorithm_quantum_part(Bitstring N, Bitstring a)
{
  int total_qubits = num_qubits_;
  int n_t = total_qubits / 2;
  int n_c = total_qubits - n_t;

  RegisterLayout layout = make_shor_layout(n_t, n_c);

  Bitstring modulus = N;

  // Reset to |0...0> and set target to |1>.
  set_basis_state(0ULL, ONE_COMPLEX);
  x(layout.target_start);

  // QFT on control register.
  qft(layout.control_start, layout.control_end);

  // Controlled modular exponentiation.
  for (int j = 0; j < n_c; ++j) {
    int control_q = layout.control_start + j;
    Bitstring power_of_a = 1ULL << j;
    controlled_modular_exponentiation(
        control_q,
        layout.target_start, layout.target_end,
        a, modulus,
        power_of_a);
  }

  // Inverse QFT on control register.
  iqft(layout.control_start, layout.control_end);

  return *this;
}
