#include "algorithms/qubo.hh"
#include "algorithms/api/grover_api.hh"
#include "logging.hh"

#include <cmath>
#include <sstream>
#include <unordered_set>

bool qubo_matrix_valid(int n, const std::vector<double>& q, std::string& error)
{
  if (n <= 0) {
    error = "n must be > 0";
    return false;
  }
  if (n >= 63) {
    error = "n must be <= 62";
    return false;
  }
  if (q.size() != static_cast<size_t>(n) * static_cast<size_t>(n)) {
    error = "matrix must contain n*n entries";
    return false;
  }
  error.clear();
  return true;
}

double qubo_evaluate(Bitstring assignment, int n, const std::vector<double>& q)
{
  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    const double xi = ((assignment >> i) & 1ULL) ? 1.0 : 0.0;
    if (xi == 0.0) continue;
    for (int j = 0; j < n; ++j) {
      const double xj = ((assignment >> j) & 1ULL) ? 1.0 : 0.0;
      if (xj == 0.0) continue;
      total += q[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)] * xi * xj;
    }
  }
  return total;
}

QuboExactResult qubo_solve_exact(int n, const std::vector<double>& q)
{
  QuboExactResult result;
  if (!qubo_matrix_valid(n, q, result.error)) {
    return result;
  }

  const Bitstring N = 1ULL << n;
  result.min_value = qubo_evaluate(0ULL, n, q);
  result.argmin = 0ULL;
  result.minimizers.clear();
  result.minimizers.push_back(0ULL);

  for (Bitstring x = 1ULL; x < N; ++x) {
    const double value = qubo_evaluate(x, n, q);
    if (value < result.min_value) {
      result.min_value = value;
      result.argmin = x;
      result.minimizers.clear();
      result.minimizers.push_back(x);
    } else if (std::abs(value - result.min_value) < 1e-12) {
      result.minimizers.push_back(x);
    }
  }

  result.ok = true;
  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "QUBO exact: n=" << n
        << " min=" << result.min_value
        << " argmin=" << result.argmin
        << " ties=" << result.minimizers.size() << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }
  return result;
}

QuboGroverResult qubo_solve_grover_threshold(int n,
                                             const std::vector<double>& q,
                                             double threshold,
                                             int iterations)
{
  QuboGroverResult result;
  std::string matrix_error;
  if (!qubo_matrix_valid(n, q, matrix_error)) {
    result.error = matrix_error;
    return result;
  }

  const Bitstring N = 1ULL << n;
  std::vector<Bitstring> targets;
  targets.reserve(static_cast<size_t>(N));
  for (Bitstring x = 0ULL; x < N; ++x) {
    if (qubo_evaluate(x, n, q) <= threshold) {
      targets.push_back(x);
    }
  }

  if (targets.empty()) {
    result.ok = true;
    result.found = false;
    result.error = "no assignments at or below threshold";
    return result;
  }

  State state(n, 0);
  GroverResult grover = run_grover(state, targets, iterations);

  Bitstring best = targets[0];
  double best_prob = -1.0;
  for (const auto& pair : state.get_state()) {
    for (size_t i = 0; i < targets.size(); ++i) {
      if (pair.first == targets[i]) {
        const double p = std::norm(pair.second);
        if (p > best_prob) {
          best_prob = p;
          best = pair.first;
        }
      }
    }
  }

  result.ok = true;
  result.found = true;
  result.candidate = best;
  result.candidate_value = qubo_evaluate(best, n, q);
  result.iterations = grover.iterations;
  result.expected_success = grover.expected_success;

  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "QUBO Grover threshold: n=" << n
        << " threshold=" << threshold
        << " candidate=" << best
        << " value=" << result.candidate_value
        << " iterations=" << result.iterations << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  return result;
}
