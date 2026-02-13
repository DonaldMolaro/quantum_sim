#include "demos/anneal_demo.hh"

#include "algorithms/anneal.hh"
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

void run_anneal_qubo_cli(const std::string& method_token,
                         int n,
                         int steps,
                         int sweeps,
                         double beta_start,
                         double beta_end,
                         int replicas,
                         const std::vector<double>& matrix)
{
  AnnealMethod method = AnnealMethod::SQA;
  if (!parse_anneal_method(method_token, method)) {
    std::cerr << "ANNEAL error: method must be SA or SQA.\n";
    return;
  }

  AnnealOptions options;
  options.method = method;
  options.steps = steps;
  options.sweeps_per_step = sweeps;
  options.beta_start = beta_start;
  options.beta_end = beta_end;
  options.replicas = replicas;

  AnnealResult result = anneal_qubo(n, matrix, options);
  if (!result.ok) {
    std::cerr << "ANNEAL error: " << result.error << "\n";
    return;
  }

  QuboExactResult exact = qubo_solve_exact(n, matrix);

  std::cout << "Anneal result (" << anneal_method_name(method) << "):\n";
  std::cout << "  best_energy=" << result.best_value << "\n";
  std::cout << "  best_state=";
  print_bits(result.best_assignment, n);
  std::cout << " (" << result.best_assignment << ")\n";
  std::cout << "  accepted_moves=" << result.accepted_moves
            << "/" << result.attempted_moves << "\n";
  std::cout << "Reference exact minimum:\n";
  std::cout << "  min_energy=" << exact.min_value << "\n";
  std::cout << "  min_state=";
  print_bits(exact.argmin, n);
  std::cout << " (" << exact.argmin << ")\n";
}

void run_anneal_demo()
{
  // Small 3-variable QUBO with a unique minimum at x=100 (4).
  const int n = 3;
  const std::vector<double> q = {
    -2.0,  0.0,  2.0,
     0.0,  1.0,  0.0,
     2.0,  0.0, -3.0
  };

  std::cout << "ANNEAL demo (n=3)\n";
  run_anneal_qubo_cli("SQA", n, 80, 20, 0.10, 6.00, 8, q);
}
