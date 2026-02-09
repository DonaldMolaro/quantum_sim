#include "algorithms/shor_quantum.hh"
#include "math/register_layout.hh"
#include <cmath>
#include <iostream>

Bitstring run_shor_algorithm_quantum_part(Bitstring N, Bitstring a, int& out_n_c)
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
  const int kMaxShorQubits = 24;
  if (total_qubits > kMaxShorQubits) {
    int max_n_t = kMaxShorQubits / 3;
    Bitstring max_n = 0;
    if (max_n_t > 0 && max_n_t < 63) {
      max_n = (1ULL << max_n_t) - 1;
    }
    std::cout << "Shor demo limited to " << kMaxShorQubits
              << " qubits; N=" << N
              << " needs " << total_qubits
              << " qubits (n_t=" << n_t << ", n_c=" << n_c << ").\n";
    if (max_n > 0) {
      std::cout << "Try N <= " << max_n << " for this simulator.\n";
    }
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

State& State::run_shor_algorithm_quantum_part(Bitstring N, Bitstring a)
{
  int total_qubits = num_qubits_;
  int n_t = total_qubits / 2;
  int n_c = total_qubits - n_t;

  RegisterLayout layout = make_shor_layout(n_t, n_c);

  Bitstring modulus = N;

  // Reset to |0...0> and set target to |1>.
  set_basis_state(0ULL, ONE_COMPLEX);
  x(layout.target.start);

  // QFT on control register.
  qft(layout.control.start, layout.control.end);

  // Controlled modular exponentiation.
  for (int j = 0; j < n_c; ++j) {
    int control_q = layout.control.start + j;
    Bitstring power_of_a = 1ULL << j;
    controlled_modular_exponentiation(
        control_q,
        layout.target.start, layout.target.end,
        a, modulus,
        power_of_a);
  }

  // Inverse QFT on control register.
  iqft(layout.control.start, layout.control.end);

  return *this;
}
