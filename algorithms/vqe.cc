#include "algorithms/vqe.hh"

#include "logging.hh"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cctype>
#include <sstream>
#include <unordered_map>

namespace {

bool is_supported_pauli(char op)
{
  const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(op)));
  return u == 'X' || u == 'Y' || u == 'Z';
}

State build_ansatz_state(int n_qubits, int layers, const std::vector<double>& params)
{
  State s(n_qubits, 0);
  size_t idx = 0;
  for (int l = 0; l < layers; ++l) {
    for (int q = 0; q < n_qubits; ++q) {
      s.ry(q, params[idx++]);
    }
    for (int q = 0; q + 1 < n_qubits; ++q) {
      s.cx(q, q + 1);
    }
  }
  return s;
}

void apply_pauli_to_basis(Bitstring in_basis,
                          const VqePauliTerm& term,
                          Bitstring& out_basis,
                          ComplexNumber& out_phase)
{
  out_basis = in_basis;
  out_phase = ONE_COMPLEX;

  for (size_t i = 0; i < term.ops.size(); ++i) {
    const char op = static_cast<char>(std::toupper(static_cast<unsigned char>(term.ops[i].op)));
    const int q = term.ops[i].qubit;
    const bool bit_is_one = ((out_basis >> q) & 1ULL) != 0ULL;

    if (op == 'X') {
      out_basis ^= (1ULL << q);
    } else if (op == 'Y') {
      out_phase *= bit_is_one ? ComplexNumber(0.0, -1.0) : ComplexNumber(0.0, 1.0);
      out_basis ^= (1ULL << q);
    } else if (op == 'Z') {
      if (bit_is_one) {
        out_phase *= -1.0;
      }
    }
  }
}

double evaluate_energy_exact(const VqeHamiltonian& h, const std::vector<double>& params, int layers)
{
  State s = build_ansatz_state(h.n_qubits, layers, params);
  return vqe_expectation_exact(s, h);
}

} // namespace

bool vqe_hamiltonian_valid(const VqeHamiltonian& h, std::string& error)
{
  if (h.n_qubits <= 0 || h.n_qubits >= 63) {
    error = "n_qubits must be in [1, 62]";
    return false;
  }
  if (h.terms.empty()) {
    error = "Hamiltonian must contain at least one term";
    return false;
  }
  for (size_t t = 0; t < h.terms.size(); ++t) {
    for (size_t i = 0; i < h.terms[t].ops.size(); ++i) {
      const VqePauliOp& op = h.terms[t].ops[i];
      if (!is_supported_pauli(op.op)) {
        error = "Unsupported Pauli op (use X, Y, or Z)";
        return false;
      }
      if (op.qubit < 0 || op.qubit >= h.n_qubits) {
        error = "Pauli op qubit index out of range";
        return false;
      }
    }
  }
  error.clear();
  return true;
}

double vqe_expectation_exact(const State& state, const VqeHamiltonian& h)
{
  std::unordered_map<Bitstring, ComplexNumber> amps;
  amps.reserve(state.get_state().size() * 2 + 1);
  const QuantumState& psi = state.get_state();
  for (size_t i = 0; i < psi.size(); ++i) {
    amps[psi[i].first] = psi[i].second;
  }

  double energy = 0.0;
  for (size_t t = 0; t < h.terms.size(); ++t) {
    const VqePauliTerm& term = h.terms[t];
    ComplexNumber accum = 0.0;

    for (std::unordered_map<Bitstring, ComplexNumber>::const_iterator it = amps.begin();
         it != amps.end(); ++it) {
      Bitstring transformed_basis = 0ULL;
      ComplexNumber phase = ONE_COMPLEX;
      apply_pauli_to_basis(it->first, term, transformed_basis, phase);
      std::unordered_map<Bitstring, ComplexNumber>::const_iterator jt = amps.find(transformed_basis);
      if (jt == amps.end()) continue;
      accum += std::conj(it->second) * phase * jt->second;
    }
    energy += term.coeff * accum.real();
  }
  return energy;
}

VqeResult run_vqe(const VqeHamiltonian& h, const VqeOptions& options)
{
  VqeResult out;
  if (!vqe_hamiltonian_valid(h, out.error)) {
    return out;
  }
  if (options.layers < 0) {
    out.error = "layers must be >= 0";
    return out;
  }
  if (options.max_iters < 0) {
    out.error = "max_iters must be >= 0";
    return out;
  }
  if (!(options.step_size > 0.0)) {
    out.error = "step_size must be > 0";
    return out;
  }
  if (options.shots != 0) {
    out.error = "shot-based VQE is not implemented yet; use shots=0";
    return out;
  }

  const size_t param_count = static_cast<size_t>(h.n_qubits) * static_cast<size_t>(options.layers);
  std::vector<double> params(param_count, 0.0);

  double best = evaluate_energy_exact(h, params, options.layers);
  out.evaluations = 1;
  out.energy_history.push_back(best);

  for (int iter = 0; iter < options.max_iters; ++iter) {
    const double step = options.step_size * std::pow(0.95, static_cast<double>(iter));
    for (size_t p = 0; p < param_count; ++p) {
      std::vector<double> plus_params = params;
      plus_params[p] += step;
      const double plus = evaluate_energy_exact(h, plus_params, options.layers);
      ++out.evaluations;

      std::vector<double> minus_params = params;
      minus_params[p] -= step;
      const double minus = evaluate_energy_exact(h, minus_params, options.layers);
      ++out.evaluations;

      if (plus < best && plus <= minus) {
        params[p] += step;
        best = plus;
      } else if (minus < best) {
        params[p] -= step;
        best = minus;
      }
    }
    out.energy_history.push_back(best);
  }

  out.ok = true;
  out.best_energy = best;
  out.best_params = params;

  if (qsim_log::enabled(qsim_log::Level::Normal)) {
    std::ostringstream msg;
    msg << "VQE: n=" << h.n_qubits
        << " layers=" << options.layers
        << " best_energy=" << out.best_energy
        << " evals=" << out.evaluations << "\n";
    qsim_log::log(qsim_log::Level::Normal, msg.str());
  }

  return out;
}
