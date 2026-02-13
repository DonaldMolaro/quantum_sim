#include "algorithms/tsp.hh"

#include "algorithms/qubo.hh"
#include "internal/limits.hh"
#include <algorithm>
#include <cmath>

namespace {

inline size_t matrix_index(int n, int row, int col)
{
  return static_cast<size_t>(row) * static_cast<size_t>(n) + static_cast<size_t>(col);
}

inline int var_index(int n_cities, int city, int pos)
{
  const int m = n_cities - 1;
  return (city - 1) * m + (pos - 1);
}

double auto_penalty(int n_cities, const std::vector<double>& distance)
{
  double max_abs = 0.0;
  for (size_t i = 0; i < distance.size(); ++i) {
    max_abs = std::max(max_abs, std::abs(distance[i]));
  }
  return max_abs * static_cast<double>(n_cities) + 1.0;
}

} // namespace

bool tsp_decode_fixed_start_assignment(int n_cities, Bitstring assignment, std::vector<int>& route)
{
  const int m = n_cities - 1;
  route.assign(static_cast<size_t>(n_cities + 1), 0);
  route[0] = 0;
  route[n_cities] = 0;

  for (int pos = 1; pos <= m; ++pos) {
    int chosen_city = -1;
    for (int city = 1; city <= m; ++city) {
      const int idx = var_index(n_cities, city, pos);
      if (((assignment >> idx) & 1ULL) == 0ULL) continue;
      if (chosen_city != -1) {
        return false;
      }
      chosen_city = city;
    }
    if (chosen_city == -1) {
      return false;
    }
    route[pos] = chosen_city;
  }

  std::vector<int> seen(static_cast<size_t>(n_cities), 0);
  for (int pos = 1; pos <= m; ++pos) {
    const int city = route[pos];
    if (city < 1 || city >= n_cities) return false;
    seen[city] += 1;
    if (seen[city] != 1) return false;
  }
  return true;
}

bool tsp_distance_matrix_valid(int n_cities, const std::vector<double>& distance, std::string& error)
{
  if (n_cities < 3) {
    error = "n_cities must be >= 3";
    return false;
  }
  const int n_vars = (n_cities - 1) * (n_cities - 1);
  if (!qsim::limits::valid_bitstring_qubit_count(n_vars)) {
    error = "n_cities too large for this simulator (fixed-start encoding needs <= 62 variables)";
    return false;
  }
  const size_t expected = static_cast<size_t>(n_cities) * static_cast<size_t>(n_cities);
  if (distance.size() != expected) {
    error = "distance matrix must contain n_cities*n_cities entries";
    return false;
  }
  for (int i = 0; i < n_cities; ++i) {
    if (std::abs(distance[matrix_index(n_cities, i, i)]) > 1e-12) {
      error = "distance matrix diagonal must be zero";
      return false;
    }
  }
  error.clear();
  return true;
}

double tsp_route_cost(const std::vector<int>& route, const std::vector<double>& distance)
{
  if (route.size() < 2) return 0.0;
  const int n_cities = static_cast<int>(std::sqrt(static_cast<double>(distance.size())));
  double total = 0.0;
  for (size_t i = 0; i + 1 < route.size(); ++i) {
    const int a = route[i];
    const int b = route[i + 1];
    total += distance[matrix_index(n_cities, a, b)];
  }
  return total;
}

TspQuboBuildResult tsp_build_qubo(int n_cities,
                                  const std::vector<double>& distance,
                                  double penalty)
{
  TspQuboBuildResult out;
  out.n_cities = n_cities;

  std::string error;
  if (!tsp_distance_matrix_valid(n_cities, distance, error)) {
    out.error = error;
    return out;
  }

  const int m = n_cities - 1;
  const int n_vars = m * m;
  out.n_vars = n_vars;
  out.penalty = (penalty > 0.0) ? penalty : auto_penalty(n_cities, distance);
  out.qubo_matrix.assign(static_cast<size_t>(n_vars) * static_cast<size_t>(n_vars), 0.0);

  // One-city-per-position constraints.
  for (int pos = 1; pos <= m; ++pos) {
    for (int city = 1; city <= m; ++city) {
      const int idx = var_index(n_cities, city, pos);
      out.qubo_matrix[matrix_index(n_vars, idx, idx)] += -out.penalty;
      for (int city2 = city + 1; city2 <= m; ++city2) {
        const int idx2 = var_index(n_cities, city2, pos);
        out.qubo_matrix[matrix_index(n_vars, idx, idx2)] += 2.0 * out.penalty;
      }
    }
  }

  // One-position-per-city constraints.
  for (int city = 1; city <= m; ++city) {
    for (int pos = 1; pos <= m; ++pos) {
      const int idx = var_index(n_cities, city, pos);
      out.qubo_matrix[matrix_index(n_vars, idx, idx)] += -out.penalty;
      for (int pos2 = pos + 1; pos2 <= m; ++pos2) {
        const int idx2 = var_index(n_cities, city, pos2);
        out.qubo_matrix[matrix_index(n_vars, idx, idx2)] += 2.0 * out.penalty;
      }
    }
  }

  // Travel objective:
  // 0 -> first city.
  for (int city = 1; city <= m; ++city) {
    const int idx = var_index(n_cities, city, 1);
    out.qubo_matrix[matrix_index(n_vars, idx, idx)] += distance[matrix_index(n_cities, 0, city)];
  }
  // between internal positions.
  for (int pos = 1; pos < m; ++pos) {
    for (int city_a = 1; city_a <= m; ++city_a) {
      const int idx_a = var_index(n_cities, city_a, pos);
      for (int city_b = 1; city_b <= m; ++city_b) {
        const int idx_b = var_index(n_cities, city_b, pos + 1);
        out.qubo_matrix[matrix_index(n_vars, idx_a, idx_b)] +=
            distance[matrix_index(n_cities, city_a, city_b)];
      }
    }
  }
  // last city -> 0.
  for (int city = 1; city <= m; ++city) {
    const int idx = var_index(n_cities, city, m);
    out.qubo_matrix[matrix_index(n_vars, idx, idx)] += distance[matrix_index(n_cities, city, 0)];
  }

  out.ok = true;
  return out;
}

TspSolveResult tsp_solve_exact(int n_cities,
                               const std::vector<double>& distance,
                               double penalty)
{
  TspSolveResult out;
  TspQuboBuildResult built = tsp_build_qubo(n_cities, distance, penalty);
  if (!built.ok) {
    out.error = built.error;
    return out;
  }

  QuboExactResult solved = qubo_solve_exact(built.n_vars, built.qubo_matrix);

  std::vector<int> route;
  (void)tsp_decode_fixed_start_assignment(n_cities, solved.argmin, route);

  out.ok = true;
  out.assignment = solved.argmin;
  out.qubo_energy = solved.min_value;
  out.route = route;
  out.route_cost = tsp_route_cost(route, distance);
  return out;
}
