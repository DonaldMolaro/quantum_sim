#include "algorithms/api/grover_api.hh"
#include "algorithms/quantum_counting.hh"
#include "internal/limits.hh"
#include "logging.hh"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <unordered_set>

GroverResult run_grover(State& s, const std::vector<Bitstring>& targets, int iterations)
{
  GroverResult result;
  if (targets.empty()) {
    result.error = "No targets provided.";
    return result;
  }

  int n_qubits = s.get_num_qubits();
  if (!qsim::limits::valid_bitstring_qubit_count(n_qubits)) {
    result.error = "Grover supports up to 62 qubits.";
    return result;
  }

  const double N = std::ldexp(1.0, n_qubits);
  const double PI = std::acos(-1.0);

  std::unordered_set<Bitstring> target_set;
  bool had_invalid = false;
  for (Bitstring t : targets) {
    if (t >= static_cast<Bitstring>(N)) {
      had_invalid = true;
      continue;
    }
    target_set.insert(t);
  }

  const double M = static_cast<double>(target_set.size());
  if (M == 0.0) {
    result.error = had_invalid
        ? "Target out of range for current number of qubits."
        : "No valid targets provided.";
    return result;
  }

  const int R_default = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N / M)));
  const int R = (iterations >= 0) ? iterations : R_default;

  const double theta = std::asin(std::sqrt(M / N));
  result.expected_success = std::sin((2.0 * R + 1.0) * theta);
  result.expected_success *= result.expected_success;

  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "Grover: N=" << static_cast<unsigned long long>(N)
        << " targets=" << target_set.size()
        << " iterations=" << R
        << " expected_success=" << result.expected_success << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }

  if (R > 0) {
    bool force_set = false;
    if (const char* env = std::getenv("QSIM_GROVER_FORCE_SET")) {
      force_set = (env[0] != '\0' && env[0] != '0');
    }
    if (!force_set && N <= 1048576.0) {
      std::vector<uint8_t> mask(static_cast<size_t>(N), 0);
      for (Bitstring t : target_set) {
        mask[static_cast<size_t>(t)] = 1;
      }
      for (int k = 0; k < R; ++k) {
        if (qsim_log::enabled(qsim_log::Level::Verbose)) {
          std::ostringstream step;
          step << "Grover: iteration " << (k + 1) << "/" << R << " (mask oracle)\n";
          qsim_log::log(qsim_log::Level::Verbose, step.str());
        }
        s.grover_oracle_Uf_mask(mask);
        s.grover_diffusion_Us();
      }
    } else {
      for (int k = 0; k < R; ++k) {
        if (qsim_log::enabled(qsim_log::Level::Verbose)) {
          std::ostringstream step;
          step << "Grover: iteration " << (k + 1) << "/" << R << " (set oracle)\n";
          qsim_log::log(qsim_log::Level::Verbose, step.str());
        }
        s.grover_oracle_Uf_multi(target_set);
        s.grover_diffusion_Us();
      }
    }
  }

  result.iterations = R;
  result.ok = true;
  return result;
}

GroverResult run_grover(int n_qubits, const std::vector<Bitstring>& targets, int iterations)
{
  GroverResult result;
  if (!qsim::limits::valid_bitstring_qubit_count(n_qubits)) {
    result.error = "Grover supports 1..62 qubits.";
    return result;
  }
  State s(n_qubits, 0);
  return run_grover(s, targets, iterations);
}

GroverResult run_grover_auto_tuned(int n_qubits,
                                   const std::vector<Bitstring>& targets,
                                   int counting_iterations)
{
  GroverResult result;
  if (!qsim::limits::valid_bitstring_qubit_count(n_qubits)) {
    result.error = "Grover supports 1..62 qubits.";
    return result;
  }

  QuantumCountingResult count = run_quantum_counting(n_qubits, targets, counting_iterations);
  if (!count.ok) {
    result.error = "Quantum counting failed: " + count.error;
    return result;
  }

  const Bitstring estimated_m = (count.estimated_targets == 0ULL) ? 1ULL : count.estimated_targets;

  const double N = std::ldexp(1.0, n_qubits);
  const int tuned_iterations =
      static_cast<int>(std::floor((std::acos(-1.0) / 4.0) * std::sqrt(N / static_cast<double>(estimated_m))));

  GroverResult grover = run_grover(n_qubits, targets, tuned_iterations);
  grover.auto_tuned = true;
  grover.estimated_targets = count.estimated_targets;
  grover.estimated_targets_real = count.estimated_targets_real;
  grover.counting_iterations = count.iterations_used;
  return grover;
}
