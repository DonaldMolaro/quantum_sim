#include "demos/qubo_demo.hh"

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

void run_qubo_exact_cli(int n, const std::vector<double>& matrix)
{
  QuboExactResult exact = qubo_solve_exact(n, matrix);
  if (!exact.ok) {
    std::cerr << "QUBO exact error: " << exact.error << "\n";
    return;
  }

  std::cout << "QUBO exact result: min=" << exact.min_value << " at ";
  print_bits(exact.argmin, n);
  std::cout << " (" << exact.argmin << ")\n";
  std::cout << "Tie count: " << exact.minimizers.size() << "\n";
}

void run_qubo_grover_cli(int n, double threshold, int iterations, const std::vector<double>& matrix)
{
  QuboGroverResult grover = qubo_solve_grover_threshold(n, matrix, threshold, iterations);
  if (!grover.ok) {
    std::cerr << "QUBO Grover error: " << grover.error << "\n";
    return;
  }
  if (!grover.found) {
    std::cout << "QUBO Grover: no assignment <= threshold\n";
    return;
  }

  std::cout << "QUBO Grover candidate: ";
  print_bits(grover.candidate, n);
  std::cout << " (" << grover.candidate << ")"
            << " value=" << grover.candidate_value
            << " iterations=" << grover.iterations
            << " expected_success=" << grover.expected_success << "\n";
}

void run_qubo_demo()
{
  // Small 3-variable QUBO with a unique minimum at x=100 (4).
  const int n = 3;
  const std::vector<double> q = {
    -2.0,  0.0,  2.0,
     0.0,  1.0,  0.0,
     2.0,  0.0, -3.0
  };

  std::cout << "QUBO demo (n=3)\n";
  run_qubo_exact_cli(n, q);

  QuboExactResult exact = qubo_solve_exact(n, q);
  if (exact.ok) {
    run_qubo_grover_cli(n, exact.min_value, -1, q);
  }
}
