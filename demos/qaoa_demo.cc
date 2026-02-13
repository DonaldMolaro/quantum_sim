#include "demos/qaoa_demo.hh"

#include "algorithms/qaoa.hh"
#include "algorithms/qubo.hh"
#include <iostream>
#include <vector>

static void print_bits(Bitstring x, int n)
{
  std::cout << "0b";
  for (int i = n - 1; i >= 0; --i) {
    std::cout << (((x >> i) & 1ULL) ? '1' : '0');
  }
}

void run_qaoa_qubo_cli(int n,
                       int p_layers,
                       int shots,
                       int max_iters,
                       double step_size,
                       const std::vector<double>& matrix)
{
  QaoaOptions options;
  options.p_layers = p_layers;
  options.shots = shots;
  options.max_iters = max_iters;
  options.step_size = step_size;

  QaoaResult result = run_qaoa_qubo(n, matrix, options);
  if (!result.ok) {
    std::cerr << "QAOA error: " << result.error << "\n";
    return;
  }

  QuboExactResult exact = qubo_solve_exact(n, matrix);

  std::cout << "QAOA result:\n";
  std::cout << "  best_energy=" << result.best_energy << "\n";
  std::cout << "  best_state=";
  print_bits(result.best_bitstring, n);
  std::cout << " (" << result.best_bitstring << ")\n";
  std::cout << "  evaluations=" << result.evaluations << "\n";
  std::cout << "Reference exact minimum:\n";
  std::cout << "  min_energy=" << exact.min_value << "\n";
  std::cout << "  min_state=";
  print_bits(exact.argmin, n);
  std::cout << " (" << exact.argmin << ")\n";
}

void run_qaoa_demo()
{
  // Small 3-variable QUBO with a unique minimum at x=100 (4).
  const int n = 3;
  const std::vector<double> q = {
    -2.0,  0.0,  2.0,
     0.0,  1.0,  0.0,
     2.0,  0.0, -3.0
  };

  std::cout << "QAOA demo (n=3)\n";
  run_qaoa_qubo_cli(n, 1, 0, 40, 0.25, q);
}
