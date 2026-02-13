#pragma once

#include "state.hh"
#include <string>
#include <vector>

struct TspQuboBuildResult {
  bool ok = false;
  std::string error;
  int n_cities = 0;
  int n_vars = 0;
  double penalty = 0.0;
  std::vector<double> qubo_matrix;
};

struct TspSolveResult {
  bool ok = false;
  std::string error;
  std::vector<int> route;
  double route_cost = 0.0;
  double qubo_energy = 0.0;
  Bitstring assignment = 0ULL;
};

bool tsp_distance_matrix_valid(int n_cities, const std::vector<double>& distance, std::string& error);
double tsp_route_cost(const std::vector<int>& route, const std::vector<double>& distance);
bool tsp_decode_fixed_start_assignment(int n_cities, Bitstring assignment, std::vector<int>& route);
TspQuboBuildResult tsp_build_qubo(int n_cities,
                                  const std::vector<double>& distance,
                                  double penalty = -1.0);
TspSolveResult tsp_solve_exact(int n_cities,
                               const std::vector<double>& distance,
                               double penalty = -1.0);
