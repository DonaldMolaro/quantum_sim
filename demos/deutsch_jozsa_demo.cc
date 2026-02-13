#include "demos/deutsch_jozsa_demo.hh"

#include "algorithms/deutsch_jozsa.hh"
#include "internal/format_utils.hh"
#include <iostream>
#include <string>
#include <vector>

void run_deutsch_jozsa_demo(int n_inputs, const std::string& oracle_token)
{
  DeutschJozsaOracle oracle;
  if (!parse_deutsch_jozsa_oracle(oracle_token, oracle)) {
    std::cerr << "Deutsch-Jozsa: unknown oracle '" << oracle_token << "'.\n";
    std::cerr << "  Use CONST0, CONST1, BALANCED_XOR0, or BALANCED_PARITY.\n";
    return;
  }

  DeutschJozsaResult result = run_deutsch_jozsa(n_inputs, oracle);
  if (!result.ok) {
    std::cerr << "Deutsch-Jozsa error: " << result.error << "\n";
    return;
  }

  std::cout << "Deutsch-Jozsa with n=" << n_inputs
            << " oracle=" << deutsch_jozsa_oracle_name(oracle) << "\n";
  std::cout << "Measured input register: " << qsim::format::bits_to_string_msb(result.input_measurement) << "\n";
  if (result.is_constant) {
    std::cout << "Result: CONSTANT function\n";
  } else {
    std::cout << "Result: BALANCED function\n";
  }
}
