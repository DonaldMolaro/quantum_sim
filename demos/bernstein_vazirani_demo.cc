#include "demos/bernstein_vazirani_demo.hh"

#include "algorithms/bernstein_vazirani.hh"
#include "internal/format_utils.hh"
#include <iostream>
#include <string>
#include <vector>

void run_bernstein_vazirani_demo(int n_inputs, Bitstring secret, int bias)
{
  BernsteinVaziraniResult result = run_bernstein_vazirani(n_inputs, secret, bias);
  if (!result.ok) {
    std::cerr << "Bernstein-Vazirani error: " << result.error << "\n";
    return;
  }

  std::cout << "Bernstein-Vazirani with n=" << n_inputs
            << " secret=" << secret
            << " bias=" << bias << "\n";
  std::cout << "Measured bits: " << qsim::format::bits_to_string_msb(result.measured_bits) << "\n";
  std::cout << "Measured secret: " << result.measured_secret << "\n";
}
