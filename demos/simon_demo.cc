#include "demos/simon_demo.hh"

#include "algorithms/simon.hh"
#include <iostream>

static void print_bits(Bitstring value, int n)
{
  for (int i = n - 1; i >= 0; --i) {
    std::cout << (((value >> i) & 1ULL) ? '1' : '0');
  }
}

void run_simon_cli(int n_inputs, Bitstring secret, int shots)
{
  SimonResult result = run_simon(n_inputs, secret, shots);
  if (!result.ok) {
    std::cerr << "SIMON error: " << result.error << "\n";
    return;
  }

  std::cout << "Simon result:\n";
  std::cout << "  n_inputs=" << result.n_inputs << "\n";
  std::cout << "  secret=";
  print_bits(result.secret, result.n_inputs);
  std::cout << "\n";
  std::cout << "  recovered=";
  print_bits(result.recovered_secret, result.n_inputs);
  std::cout << "\n";
  std::cout << "  equations=" << result.equations.size() << "\n";
  std::cout << "  verified=" << (result.verified ? "true" : "false") << "\n";
}

void run_simon_demo()
{
  std::cout << "Simon demo (n=4, secret=1010)\n";
  run_simon_cli(4, 0b1010, -1);
}
