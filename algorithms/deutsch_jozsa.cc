#include "algorithms/deutsch_jozsa.hh"

#include <algorithm>
#include <cctype>
#include <string>

static void apply_oracle(State& state, int n_inputs, DeutschJozsaOracle oracle)
{
  const int ancilla = n_inputs;
  switch (oracle) {
    case DeutschJozsaOracle::ConstantZero:
      return;
    case DeutschJozsaOracle::ConstantOne:
      state.x(ancilla);
      return;
    case DeutschJozsaOracle::BalancedXor0:
      state.cx(0, ancilla);
      return;
    case DeutschJozsaOracle::BalancedParity:
      for (int j = 0; j < n_inputs; ++j) {
        state.cx(j, ancilla);
      }
      return;
  }
}

DeutschJozsaResult run_deutsch_jozsa(int n_inputs, DeutschJozsaOracle oracle)
{
  DeutschJozsaResult result;
  if (n_inputs <= 0) {
    result.error = "n_inputs must be > 0";
    return result;
  }
  if (n_inputs + 1 > 63) {
    result.error = "n_inputs too large for demo";
    return result;
  }

  State state(n_inputs + 1, n_inputs + 1);

  const int ancilla = n_inputs;
  state.x(ancilla);
  for (int j = 0; j <= n_inputs; ++j) {
    state.h(j);
  }

  apply_oracle(state, n_inputs, oracle);

  for (int j = 0; j < n_inputs; ++j) {
    state.h(j);
  }

  std::vector<int> bits;
  state.measure_all(bits);
  result.input_measurement.assign(bits.begin(), bits.begin() + n_inputs);

  result.is_constant = std::all_of(result.input_measurement.begin(),
                                   result.input_measurement.end(),
                                   [](int v) { return v == 0; });
  result.ok = true;
  return result;
}

const char* deutsch_jozsa_oracle_name(DeutschJozsaOracle oracle)
{
  switch (oracle) {
    case DeutschJozsaOracle::ConstantZero:
      return "CONST0";
    case DeutschJozsaOracle::ConstantOne:
      return "CONST1";
    case DeutschJozsaOracle::BalancedXor0:
      return "BALANCED_XOR0";
    case DeutschJozsaOracle::BalancedParity:
      return "BALANCED_PARITY";
  }
  return "UNKNOWN";
}

bool parse_deutsch_jozsa_oracle(const std::string& token, DeutschJozsaOracle& out)
{
  std::string upper = token;
  std::transform(upper.begin(), upper.end(), upper.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

  if (upper == "CONST0" || upper == "CONSTANT0" || upper == "CONSTANT_ZERO") {
    out = DeutschJozsaOracle::ConstantZero;
    return true;
  }
  if (upper == "CONST1" || upper == "CONSTANT1" || upper == "CONSTANT_ONE") {
    out = DeutschJozsaOracle::ConstantOne;
    return true;
  }
  if (upper == "BALANCED_XOR0" || upper == "BALANCED" || upper == "XOR0") {
    out = DeutschJozsaOracle::BalancedXor0;
    return true;
  }
  if (upper == "BALANCED_PARITY" || upper == "PARITY") {
    out = DeutschJozsaOracle::BalancedParity;
    return true;
  }
  return false;
}
