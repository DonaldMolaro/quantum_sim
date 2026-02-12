#include "algorithms/bernstein_vazirani.hh"
#include "logging.hh"

#include <sstream>

BernsteinVaziraniResult run_bernstein_vazirani(int n_inputs, Bitstring secret, int bias)
{
  BernsteinVaziraniResult result;

  if (n_inputs <= 0) {
    result.error = "n_inputs must be > 0";
    return result;
  }
  if (n_inputs + 1 >= 63) {
    result.error = "n_inputs too large for 64-bit demo";
    return result;
  }
  if (bias != 0 && bias != 1) {
    result.error = "bias must be 0 or 1";
    return result;
  }

  const Bitstring max_secret = (1ULL << n_inputs) - 1ULL;
  if (secret > max_secret) {
    result.error = "secret out of range for n_inputs";
    return result;
  }

  State state(n_inputs + 1, n_inputs + 1);
  const int ancilla = n_inputs;

  state.x(ancilla);
  for (int j = 0; j <= n_inputs; ++j) {
    state.h(j);
  }

  // Oracle for f(x)=s.x xor bias
  if (bias == 1) {
    state.x(ancilla);
  }
  for (int j = 0; j < n_inputs; ++j) {
    if ((secret >> j) & 1ULL) {
      state.cx(j, ancilla);
    }
  }

  for (int j = 0; j < n_inputs; ++j) {
    state.h(j);
  }

  std::vector<int> bits;
  state.measure_all(bits);
  result.measured_bits.assign(bits.begin(), bits.begin() + n_inputs);

  Bitstring measured = 0ULL;
  for (int j = 0; j < n_inputs; ++j) {
    if (result.measured_bits[j]) {
      measured |= (1ULL << j);
    }
  }

  result.measured_secret = measured;
  result.ok = true;

  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "Bernstein-Vazirani: n=" << n_inputs
        << " secret=" << secret
        << " measured=" << measured << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  return result;
}
