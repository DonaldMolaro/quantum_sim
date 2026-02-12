#include "demos/bernstein_vazirani_demo.hh"

#include "algorithms/bernstein_vazirani.hh"
#include <iostream>
#include <string>
#include <vector>

static std::string bits_to_string(const std::vector<int>& bits)
{
  if (bits.empty()) return "0";
  std::string out;
  out.reserve(bits.size());
  for (int i = static_cast<int>(bits.size()) - 1; i >= 0; --i) {
    out.push_back(bits[i] ? '1' : '0');
  }
  return out;
}

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
  std::cout << "Measured bits: " << bits_to_string(result.measured_bits) << "\n";
  std::cout << "Measured secret: " << result.measured_secret << "\n";
}
