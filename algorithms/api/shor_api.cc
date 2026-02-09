#include "algorithms/api/shor_api.hh"
#include "logging.hh"
#include <cmath>
#include <sstream>

ShorResult run_shor_quantum_part(Bitstring N, Bitstring a)
{
  ShorResult result;

  if (N < 2) {
    result.error = "N must be >= 2.";
    return result;
  }

  int n_t = static_cast<int>(std::ceil(std::log2(static_cast<double>(N))));
  int n_c = 2 * n_t;
  int total_qubits = n_c + n_t;
  int num_cbits = n_c;

  if (n_c >= 63) {
    result.error = "Control register too large for 64-bit demo.";
    return result;
  }

  result.n_c = n_c;

  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "Shor (quantum part): N=" << N << " a=" << a
        << " n_t=" << n_t << " n_c=" << n_c << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  State s(total_qubits, num_cbits);
  s.x(total_qubits - 1);
  s.qft(0, n_c - 1);

  for (int j = 0; j < n_c; ++j) {
    Bitstring power_of_a = 1ULL << j;
    s.controlled_modular_exponentiation(j, n_c, total_qubits - 1, a, N, power_of_a);
  }

  s.iqft(0, n_c - 1);
  for (int j = 0; j < n_c; ++j) {
    s.measure(j, j);
  }

  Bitstring measured_x = 0;
  for (int j = 0; j < n_c; ++j) {
    measured_x |= (static_cast<Bitstring>(s.get_cbit(j)) << j);
  }

  result.measured_x = measured_x;
  result.ok = true;
  return result;
}
