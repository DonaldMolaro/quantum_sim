#include "algorithms/anneal.hh"

#include "algorithms/qubo.hh"
#include "logging.hh"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

Bitstring random_assignment(int n, std::mt19937& rng)
{
  std::bernoulli_distribution bit(0.5);
  Bitstring out = 0ULL;
  for (int i = 0; i < n; ++i) {
    if (bit(rng)) {
      out |= (1ULL << i);
    }
  }
  return out;
}

double lerp(double a, double b, double t)
{
  return a + (b - a) * t;
}

AnnealResult run_sa(int n, const std::vector<double>& q, const AnnealOptions& options)
{
  AnnealResult out;
  std::mt19937 rng(options.seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_int_distribution<int> pick_bit(0, n - 1);

  Bitstring current = random_assignment(n, rng);
  double current_e = qubo_evaluate(current, n, q);

  out.best_assignment = current;
  out.best_value = current_e;
  out.best_history.reserve(static_cast<size_t>(options.steps) + 1);
  out.best_history.push_back(out.best_value);

  for (int step = 0; step < options.steps; ++step) {
    const double t = (options.steps > 1)
                       ? static_cast<double>(step) / static_cast<double>(options.steps - 1)
                       : 1.0;
    const double beta = lerp(options.beta_start, options.beta_end, t);

    for (int sweep = 0; sweep < options.sweeps_per_step; ++sweep) {
      for (int move = 0; move < n; ++move) {
        const int bit = pick_bit(rng);
        const Bitstring candidate = current ^ (1ULL << bit);
        const double cand_e = qubo_evaluate(candidate, n, q);
        const double delta = cand_e - current_e;

        ++out.attempted_moves;
        bool accept = false;
        if (delta <= 0.0) {
          accept = true;
        } else {
          const double p = std::exp(-beta * delta);
          accept = (unit(rng) < p);
        }

        if (accept) {
          ++out.accepted_moves;
          current = candidate;
          current_e = cand_e;
          if (current_e < out.best_value) {
            out.best_value = current_e;
            out.best_assignment = current;
          }
        }
      }
    }
    out.best_history.push_back(out.best_value);
  }

  out.ok = true;
  return out;
}

int bit_of(Bitstring x, int i)
{
  return ((x >> i) & 1ULL) ? 1 : 0;
}

AnnealResult run_sqa(int n, const std::vector<double>& q, const AnnealOptions& options)
{
  AnnealResult out;
  std::mt19937 rng(options.seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_int_distribution<int> pick_bit(0, n - 1);
  std::uniform_int_distribution<int> pick_replica(0, options.replicas - 1);

  std::vector<Bitstring> world(static_cast<size_t>(options.replicas), 0ULL);
  std::vector<double> world_e(static_cast<size_t>(options.replicas), 0.0);
  for (int r = 0; r < options.replicas; ++r) {
    world[static_cast<size_t>(r)] = random_assignment(n, rng);
    world_e[static_cast<size_t>(r)] = qubo_evaluate(world[static_cast<size_t>(r)], n, q);
  }

  out.best_assignment = world[0];
  out.best_value = world_e[0];
  for (int r = 1; r < options.replicas; ++r) {
    if (world_e[static_cast<size_t>(r)] < out.best_value) {
      out.best_value = world_e[static_cast<size_t>(r)];
      out.best_assignment = world[static_cast<size_t>(r)];
    }
  }

  out.best_history.reserve(static_cast<size_t>(options.steps) + 1);
  out.best_history.push_back(out.best_value);

  for (int step = 0; step < options.steps; ++step) {
    const double t = (options.steps > 1)
                       ? static_cast<double>(step) / static_cast<double>(options.steps - 1)
                       : 1.0;
    const double beta = lerp(options.beta_start, options.beta_end, t);
    const double problem_scale = t;
    const double coupling = (1.0 - t) * 1.5; // stronger transverse influence early

    for (int sweep = 0; sweep < options.sweeps_per_step; ++sweep) {
      for (int move = 0; move < n * options.replicas; ++move) {
        const int r = pick_replica(rng);
        const int bit = pick_bit(rng);
        const int r_prev = (r == 0) ? (options.replicas - 1) : (r - 1);
        const int r_next = (r + 1) % options.replicas;

        const Bitstring cur = world[static_cast<size_t>(r)];
        const Bitstring candidate = cur ^ (1ULL << bit);

        const double cur_problem = world_e[static_cast<size_t>(r)];
        const double cand_problem = qubo_evaluate(candidate, n, q);
        double delta = problem_scale * (cand_problem - cur_problem);

        const int cur_bit = bit_of(cur, bit);
        const int flip_bit = 1 - cur_bit;
        const int prev_bit = bit_of(world[static_cast<size_t>(r_prev)], bit);
        const int next_bit = bit_of(world[static_cast<size_t>(r_next)], bit);

        const int before = (cur_bit != prev_bit ? 1 : 0) + (cur_bit != next_bit ? 1 : 0);
        const int after = (flip_bit != prev_bit ? 1 : 0) + (flip_bit != next_bit ? 1 : 0);
        delta += coupling * static_cast<double>(after - before);

        ++out.attempted_moves;
        bool accept = false;
        if (delta <= 0.0) {
          accept = true;
        } else {
          const double p = std::exp(-beta * delta);
          accept = (unit(rng) < p);
        }

        if (accept) {
          ++out.accepted_moves;
          world[static_cast<size_t>(r)] = candidate;
          world_e[static_cast<size_t>(r)] = cand_problem;
          if (cand_problem < out.best_value) {
            out.best_value = cand_problem;
            out.best_assignment = candidate;
          }
        }
      }
    }

    out.best_history.push_back(out.best_value);
  }

  out.ok = true;
  return out;
}

} // namespace

bool anneal_options_valid(int n,
                          const std::vector<double>& q,
                          const AnnealOptions& options,
                          std::string& error)
{
  if (!qubo_matrix_valid(n, q, error)) {
    return false;
  }
  if (options.steps <= 0) {
    error = "steps must be > 0";
    return false;
  }
  if (options.sweeps_per_step <= 0) {
    error = "sweeps_per_step must be > 0";
    return false;
  }
  if (!(options.beta_start > 0.0) || !(options.beta_end > 0.0)) {
    error = "beta_start and beta_end must be > 0";
    return false;
  }
  if (options.beta_end < options.beta_start) {
    error = "beta_end must be >= beta_start";
    return false;
  }
  if (options.method == AnnealMethod::SQA && options.replicas < 2) {
    error = "SQA requires replicas >= 2";
    return false;
  }
  if (options.method == AnnealMethod::SA && options.replicas <= 0) {
    error = "replicas must be positive";
    return false;
  }
  error.clear();
  return true;
}

AnnealResult anneal_qubo(int n,
                         const std::vector<double>& q,
                         const AnnealOptions& options)
{
  AnnealResult out;
  if (!anneal_options_valid(n, q, options, out.error)) {
    return out;
  }

  out = (options.method == AnnealMethod::SQA)
          ? run_sqa(n, q, options)
          : run_sa(n, q, options);

  if (out.ok && qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "Anneal " << anneal_method_name(options.method)
        << ": n=" << n
        << " best=" << out.best_value
        << " assignment=" << out.best_assignment
        << " accepted=" << out.accepted_moves
        << "/" << out.attempted_moves << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  return out;
}

const char* anneal_method_name(AnnealMethod method)
{
  switch (method) {
    case AnnealMethod::SA:
      return "SA";
    case AnnealMethod::SQA:
      return "SQA";
  }
  return "UNKNOWN";
}

bool parse_anneal_method(const std::string& token, AnnealMethod& out)
{
  std::string up = token;
  std::transform(up.begin(), up.end(), up.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  if (up == "SA") {
    out = AnnealMethod::SA;
    return true;
  }
  if (up == "SQA" || up == "QUANTUM") {
    out = AnnealMethod::SQA;
    return true;
  }
  return false;
}
