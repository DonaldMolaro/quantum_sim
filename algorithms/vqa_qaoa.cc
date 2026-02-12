#include "algorithms/vqa_qaoa.hh"

#include "algorithms/qubo.hh"
#include "logging.hh"

#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

namespace {

struct PreparedQubo {
  std::vector<double> energies;
  Bitstring dim = 0ULL;
};

bool validate_options(int n,
                      const std::vector<double>& q,
                      const VqaQaoaOptions& options,
                      std::string& error)
{
  if (!qubo_matrix_valid(n, q, error)) {
    return false;
  }
  if (options.p_layers < 0) {
    error = "p_layers must be >= 0";
    return false;
  }
  if (options.max_iters < 0) {
    error = "max_iters must be >= 0";
    return false;
  }
  if (options.shots < 0) {
    error = "shots must be >= 0";
    return false;
  }
  if (!(options.step_size > 0.0)) {
    error = "step_size must be > 0";
    return false;
  }
  error.clear();
  return true;
}

PreparedQubo prepare_qubo(int n, const std::vector<double>& q)
{
  PreparedQubo prepared;
  prepared.dim = (1ULL << n);
  prepared.energies.resize(static_cast<size_t>(prepared.dim), 0.0);
  for (Bitstring x = 0ULL; x < prepared.dim; ++x) {
    prepared.energies[static_cast<size_t>(x)] = qubo_evaluate(x, n, q);
  }
  return prepared;
}

void apply_cost_layer(State& state, int n, double gamma, const PreparedQubo& prepared)
{
  (void)n;
  QuantumState next = state.get_state();
  for (size_t i = 0; i < next.size(); ++i) {
    const Bitstring basis = next[i].first;
    const double e = prepared.energies[static_cast<size_t>(basis)];
    const double phase = -gamma * e;
    const ComplexNumber rot(std::cos(phase), std::sin(phase));
    next[i].second *= rot;
  }
  state.set_superposition(next);
}

State build_qaoa_state(int n,
                       const std::vector<double>& gamma,
                       const std::vector<double>& beta,
                       const PreparedQubo& prepared)
{
  State state(n, 0);
  for (int j = 0; j < n; ++j) {
    state.h(j);
  }
  for (size_t layer = 0; layer < gamma.size(); ++layer) {
    apply_cost_layer(state, n, gamma[layer], prepared);
    for (int j = 0; j < n; ++j) {
      state.rx(j, 2.0 * beta[layer]);
    }
  }
  return state;
}

double exact_energy_expectation(const State& state, const PreparedQubo& prepared)
{
  double e = 0.0;
  const QuantumState& amps = state.get_state();
  for (size_t i = 0; i < amps.size(); ++i) {
    const Bitstring b = amps[i].first;
    const double p = std::norm(amps[i].second);
    e += p * prepared.energies[static_cast<size_t>(b)];
  }
  return e;
}

Bitstring most_likely_basis(const State& state)
{
  Bitstring best = 0ULL;
  double best_prob = -1.0;
  const QuantumState& amps = state.get_state();
  for (size_t i = 0; i < amps.size(); ++i) {
    const double p = std::norm(amps[i].second);
    if (p > best_prob) {
      best_prob = p;
      best = amps[i].first;
    }
  }
  return best;
}

double sample_energy_expectation(const State& state,
                                 const PreparedQubo& prepared,
                                 int shots,
                                 std::mt19937& rng)
{
  struct Entry {
    Bitstring basis;
    double cumulative;
  };

  std::vector<Entry> table;
  table.reserve(state.get_state().size());
  double cdf = 0.0;
  const QuantumState& amps = state.get_state();
  for (size_t i = 0; i < amps.size(); ++i) {
    const double p = std::norm(amps[i].second);
    if (p <= 0.0) continue;
    cdf += p;
    table.push_back(Entry{amps[i].first, cdf});
  }
  table.back().cumulative = 1.0;

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  double total = 0.0;
  for (int s = 0; s < shots; ++s) {
    const double r = dist(rng);
    size_t idx = 0;
    while (idx + 1 < table.size() && table[idx].cumulative < r) {
      ++idx;
    }
    total += prepared.energies[static_cast<size_t>(table[idx].basis)];
  }
  return total / static_cast<double>(shots);
}

double evaluate_energy(int n,
                       const std::vector<double>& gamma,
                       const std::vector<double>& beta,
                       const PreparedQubo& prepared,
                       int shots,
                       std::mt19937& rng,
                       Bitstring* likely_state)
{
  State state = build_qaoa_state(n, gamma, beta, prepared);
  if (likely_state != nullptr) {
    *likely_state = most_likely_basis(state);
  }
  if (shots == 0) {
    return exact_energy_expectation(state, prepared);
  }
  return sample_energy_expectation(state, prepared, shots, rng);
}

} // namespace

VqaQaoaResult run_vqa_qaoa_qubo(int n,
                                const std::vector<double>& q,
                                const VqaQaoaOptions& options)
{
  VqaQaoaResult out;
  if (!validate_options(n, q, options, out.error)) {
    return out;
  }

  const PreparedQubo prepared = prepare_qubo(n, q);

  std::vector<double> gamma(static_cast<size_t>(options.p_layers), 0.0);
  std::vector<double> beta(static_cast<size_t>(options.p_layers), 0.0);
  std::mt19937 rng(options.seed);

  Bitstring init_likely = 0ULL;
  double best = evaluate_energy(n, gamma, beta, prepared, options.shots, rng, &init_likely);
  out.evaluations = 1;
  out.energy_history.push_back(best);
  Bitstring best_likely = init_likely;

  for (int iter = 0; iter < options.max_iters; ++iter) {
    const double step = options.step_size * std::pow(0.95, static_cast<double>(iter));

    for (int idx = 0; idx < options.p_layers; ++idx) {
      std::vector<double> trial_g = gamma;
      trial_g[static_cast<size_t>(idx)] += step;
      Bitstring plus_likely = 0ULL;
      const double plus = evaluate_energy(n, trial_g, beta, prepared, options.shots, rng, &plus_likely);
      ++out.evaluations;

      trial_g[static_cast<size_t>(idx)] = gamma[static_cast<size_t>(idx)] - step;
      Bitstring minus_likely = 0ULL;
      const double minus = evaluate_energy(n, trial_g, beta, prepared, options.shots, rng, &minus_likely);
      ++out.evaluations;

      if (plus < best && plus <= minus) {
        gamma[static_cast<size_t>(idx)] += step;
        best = plus;
        best_likely = plus_likely;
      } else if (minus < best) {
        gamma[static_cast<size_t>(idx)] -= step;
        best = minus;
        best_likely = minus_likely;
      }
    }

    for (int idx = 0; idx < options.p_layers; ++idx) {
      std::vector<double> trial_b = beta;
      trial_b[static_cast<size_t>(idx)] += step;
      Bitstring plus_likely = 0ULL;
      const double plus = evaluate_energy(n, gamma, trial_b, prepared, options.shots, rng, &plus_likely);
      ++out.evaluations;

      trial_b[static_cast<size_t>(idx)] = beta[static_cast<size_t>(idx)] - step;
      Bitstring minus_likely = 0ULL;
      const double minus = evaluate_energy(n, gamma, trial_b, prepared, options.shots, rng, &minus_likely);
      ++out.evaluations;

      if (plus < best && plus <= minus) {
        beta[static_cast<size_t>(idx)] += step;
        best = plus;
        best_likely = plus_likely;
      } else if (minus < best) {
        beta[static_cast<size_t>(idx)] -= step;
        best = minus;
        best_likely = minus_likely;
      }
    }

    out.energy_history.push_back(best);
  }

  out.ok = true;
  out.best_energy = best;
  out.best_gamma = gamma;
  out.best_beta = beta;
  out.best_bitstring = best_likely;

  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "VQA QAOA: n=" << n
        << " p=" << options.p_layers
        << " best_energy=" << out.best_energy
        << " best_state=" << out.best_bitstring
        << " evals=" << out.evaluations << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  return out;
}
