#include "algorithms/shor_quantum.hh"
#include "internal/limits.hh"
#include "internal/register_layout.hh"
#include "internal/step_observer.hh"
#include <cmath>
#include <string>

ShorQuantumResult run_shor_algorithm_quantum_part(Bitstring N, Bitstring a)
{
  ShorQuantumResult result;
  result.n_t = static_cast<int>(std::ceil(std::log2(static_cast<double>(N))));
  result.n_c = 2 * result.n_t;
  result.total_qubits = result.n_c + result.n_t;
  int num_cbits = result.n_c;

  if (result.n_c >= 63) {
    result.error = "Control register too large for 64-bit demo (n_c=" + std::to_string(result.n_c) + ").";
    return result;
  }
  if (result.total_qubits > qsim::limits::kMaxShorQubits) {
    int max_n_t = qsim::limits::kMaxShorQubits / 3;
    if (max_n_t > 0 && max_n_t < 63) {
      result.max_n = (1ULL << max_n_t) - 1;
    }
    result.error = "Shor demo limited to " + std::to_string(qsim::limits::kMaxShorQubits) +
                   " qubits; N=" + std::to_string(N) +
                   " needs " + std::to_string(result.total_qubits) +
                   " qubits (n_t=" + std::to_string(result.n_t) +
                   ", n_c=" + std::to_string(result.n_c) + ").";
    return result;
  }

  State s(result.total_qubits, num_cbits);
  qsim::notify_init(s);

  // Initialize target register to |1> (rightmost target qubit).
  s.x(result.total_qubits - 1);
  qsim::notify_gate("X " + std::to_string(result.total_qubits - 1), s);

  // QFT on control register (creates uniform superposition).
  s.qft(0, result.n_c - 1);
  qsim::notify_block("Shor QFT", "QFT", "control register: uniform superposition", s, 0, result.n_c - 1);

  // Controlled modular exponentiation for each control qubit.
  for (int j = 0; j < result.n_c; ++j) {
    int control_q = j;
    Bitstring power_of_a = 1ULL << j;
    s.controlled_modular_exponentiation(
        control_q,
        result.n_c, result.total_qubits - 1,
        a, N,
        power_of_a);
    qsim::notify_block("Shor controlled-U^(2^" + std::to_string(j) + ")", "U" + std::to_string(j),
                       "control q" + std::to_string(j) + ": target <- " + std::to_string(a) + "^(2^" +
                           std::to_string(j) + ") * target mod " + std::to_string(N),
                       s, j, result.total_qubits - 1);
  }

  // Inverse QFT on control register.
  s.iqft(0, result.n_c - 1);
  qsim::notify_block("Shor inverse QFT", "iQFT", "control register: period appears as peaks", s, 0, result.n_c - 1);

  // Measure control register into classical bits.
  for (int j = 0; j < result.n_c; ++j) {
    s.measure(j, j);
    qsim::notify_gate("MEASURE " + std::to_string(j) + " " + std::to_string(j), s);
  }

  Bitstring measured_x = 0;
  for (int j = 0; j < result.n_c; ++j) {
    int bit = s.get_cbit(j);
    measured_x |= static_cast<Bitstring>(bit) << j;
  }

  result.measured_x = measured_x;
  result.ok = true;
  return result;
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
