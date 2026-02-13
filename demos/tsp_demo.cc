#include "demos/tsp_demo.hh"

#include "algorithms/tsp.hh"
#include <iostream>

static void print_route(const std::vector<int>& route)
{
  for (size_t i = 0; i < route.size(); ++i) {
    if (i != 0) std::cout << " -> ";
    std::cout << route[i];
  }
}

void run_tsp_exact_cli(int n_cities,
                       double penalty,
                       const std::vector<double>& distance)
{
  TspSolveResult result = tsp_solve_exact(n_cities, distance, penalty);
  if (!result.ok) {
    std::cerr << "TSP error: " << result.error << "\n";
    return;
  }

  std::cout << "TSP exact result:\n";
  std::cout << "  route: ";
  print_route(result.route);
  std::cout << "\n";
  std::cout << "  route_cost=" << result.route_cost << "\n";
  std::cout << "  qubo_energy=" << result.qubo_energy << "\n";
}

void run_tsp_demo()
{
  // 4-city square-like toy instance with optimal tour cost 4.
  const int n_cities = 4;
  const std::vector<double> d = {
      0.0, 1.0, 2.0, 1.0,
      1.0, 0.0, 1.0, 2.0,
      2.0, 1.0, 0.0, 1.0,
      1.0, 2.0, 1.0, 0.0};

  std::cout << "TSP demo (n_cities=4)\n";
  run_tsp_exact_cli(n_cities, -1.0, d);
}
