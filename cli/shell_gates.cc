#include "cli/shell.hh"
#include "cli/commands.hh"
#include <iomanip>
#include <cmath>
#include <functional>
#include <iostream>
#include <unordered_map>

namespace {

using SingleQubitGate = State& (State::*)(int);
using TwoQubitGate = State& (State::*)(int, int);

} // namespace

bool QuantumShell::handle_single_qubit_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  static const std::unordered_map<std::string, SingleQubitGate> kSingleQubitGates = {
      {"H", &State::h}, {"X", &State::x}, {"Y", &State::y},
      {"Z", &State::z}, {"S", &State::s}, {"T", &State::t}};
  std::unordered_map<std::string, SingleQubitGate>::const_iterator it = kSingleQubitGates.find(cmd);
  if (it == kSingleQubitGates.end()) return false;

  int j = cli::get_arg(tokens, 1, cmd);
  if (j == -1) return true;
  try {
    (state.get()->*(it->second))(j);
    std::cout << cmd << "(" << j << ") applied.\n";
    if (cmd == "H") {
      tutor_note("H creates equal superposition from a basis state (up to phase).");
    } else if (cmd == "X") {
      tutor_note("X flips computational basis states |0> <-> |1>.");
    } else if (cmd == "Z") {
      tutor_note("Z keeps probabilities but flips phase on |1> components.");
    }
    state->display();
  } catch (const std::exception& e) {
    std::cerr << "Operation failed: " << e.what() << "\n";
  }
  return true;
}

bool QuantumShell::handle_parametric_gate_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  if (cmd == "RX" || cmd == "RY" || cmd == "RZ") {
    int j = cli::get_arg(tokens, 1, cmd);
    double theta = cli::get_angle_arg_required(tokens, 2, cmd);
    if (j == -1 || std::isnan(theta)) return true;
    try {
      if (cmd == "RX") state->rx(j, theta);
      else if (cmd == "RY") state->ry(j, theta);
      else state->rz(j, theta);
      std::cout << cmd << "(" << j << ", " << theta << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "RU") {
    int j = cli::get_arg(tokens, 1, cmd);
    double theta = cli::get_angle_arg_required(tokens, 2, cmd);
    double phi = cli::get_angle_arg_required(tokens, 3, cmd);
    double lambda = cli::get_angle_arg_required(tokens, 4, cmd);
    if (j == -1 || std::isnan(theta) || std::isnan(phi) || std::isnan(lambda)) return true;
    try {
      state->ru(j, theta, phi, lambda);
      std::cout << "RU(" << j << ", " << theta << ", " << phi << ", " << lambda << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "CRZ" || cmd == "CRX" || cmd == "CRY") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    double theta = cli::get_angle_arg_required(tokens, 3, cmd);
    if (j == -1 || k == -1 || std::isnan(theta)) return true;
    try {
      if (cmd == "CRZ") {
        state->crz(j, k, theta);
        std::cout << "CRZ(" << j << ", " << k << ", " << theta << ") applied.\n";
      } else if (cmd == "CRX") {
        state->crx(j, k, theta);
        std::cout << "CRX(" << j << ", " << k << ", " << theta << ") applied.\n";
      } else {
        state->cry(j, k, theta);
        std::cout << "CRY(" << j << ", " << k << ", " << theta << ") applied.\n";
      }
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "CU") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    double theta = cli::get_angle_arg_required(tokens, 3, cmd);
    double phi = cli::get_angle_arg_required(tokens, 4, cmd);
    double lambda = cli::get_angle_arg_required(tokens, 5, cmd);
    if (j == -1 || k == -1 || std::isnan(theta) || std::isnan(phi) || std::isnan(lambda)) return true;
    try {
      state->cu(j, k, theta, phi, lambda);
      std::cout << "CU(" << j << ", " << k << ", " << theta << ", " << phi << ", " << lambda << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }
  return false;
}

bool QuantumShell::handle_multi_qubit_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  static const std::unordered_map<std::string, TwoQubitGate> kTwoQubitGates = {
      {"CX", &State::cx}, {"CZ", &State::cz}, {"CY", &State::cy},
      {"CH", &State::ch}, {"CNOT", &State::cnot}};
  std::unordered_map<std::string, TwoQubitGate>::const_iterator it = kTwoQubitGates.find(cmd);
  if (it != kTwoQubitGates.end()) {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    if (j == -1 || k == -1) return true;
    (state.get()->*(it->second))(j, k);
    std::cout << cmd << "(" << j << ", " << k << ") applied.\n";
    state->display();
    return true;
  }

  if (cmd == "CCX" || cmd == "TOFFOLI") {
    int c1 = cli::get_arg(tokens, 1, cmd);
    int c2 = cli::get_arg(tokens, 2, cmd);
    int t = cli::get_arg(tokens, 3, cmd);
    if (c1 == -1 || c2 == -1 || t == -1) return true;
    try {
      state->ccx(c1, c2, t);
      std::cout << cmd << "(" << c1 << ", " << c2 << ", " << t << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "SWAP") {
    int j = cli::get_arg(tokens, 1, "SWAP");
    int k = cli::get_arg(tokens, 2, "SWAP");
    if (j == -1 || k == -1) return true;
    state->swap(j, k);
    std::cout << "SWAP(" << j << ", " << k << ") applied.\n";
    state->display();
    return true;
  }
  return false;
}

bool QuantumShell::handle_measurement_and_display_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  if (cmd == "MEASURE") {
    int j = cli::get_arg(tokens, 1, "MEASURE");
    int c = cli::get_arg(tokens, 2, "MEASURE");
    if (j == -1 || c == -1) return true;
    state->measure(j, c);
    std::cout << "Qubit " << j << " measured.";
    if (!state->get_cbits().empty()) {
      std::cout << " Result stored in c[" << c << "]: " << state->get_cbits()[c] << "\n";
    } else {
      std::cout << " No Register to store in \n";
    }
    tutor_note("Measurement collapses superposition on qubit j and records a classical bit.");
    state->display();
    return true;
  }
  if (cmd == "DISPLAY") {
    state->display();
    return true;
  }
  if (cmd == "CHECK") {
    if (tokens.size() < 2) {
      std::cout << "CHECK modes: NORMALIZED | BASIS <idx> | TARGETS <t1> [t2 ...] | BELL <PHI+|PHI-|PSI+|PSI->\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "NORMALIZED") {
      double total = 0.0;
      const QuantumState& qs = state->get_state();
      for (size_t i = 0; i < qs.size(); ++i) {
        total += std::norm(qs[i].second);
      }
      const double diff = std::abs(total - 1.0);
      std::cout << std::setprecision(12)
                << "CHECK NORMALIZED: total_probability=" << total
                << " diff=" << diff
                << (diff < 1e-9 ? " [PASS]\n" : " [WARN]\n");
      return true;
    }
    if (mode == "BASIS") {
      int idx = cli::get_arg(tokens, 2, "CHECK BASIS");
      if (idx == -1) return true;
      if (idx < 0) {
        std::cerr << "Error: CHECK BASIS requires idx >= 0.\n";
        return true;
      }
      const Bitstring N = (state->get_num_qubits() >= 63) ? 0ULL : (1ULL << state->get_num_qubits());
      if (N != 0ULL && static_cast<Bitstring>(idx) >= N) {
        std::cerr << "Error: basis index out of range for current qubit count.\n";
        return true;
      }
      const ComplexNumber a = state->get_amplitude(static_cast<Bitstring>(idx));
      const double p = std::norm(a);
      std::cout << "CHECK BASIS: idx=" << idx << " probability=" << p
                << (p > 0.999 ? " [PASS]\n" : " [INFO]\n");
      return true;
    }
    if (mode == "TARGETS") {
      if (tokens.size() < 3) {
        std::cerr << "Error: CHECK TARGETS requires at least one target index.\n";
        return true;
      }
      double p = 0.0;
      for (size_t i = 2; i < tokens.size(); ++i) {
        int t = cli::get_arg(tokens, i, "CHECK TARGETS");
        if (t == -1) return true;
        if (t < 0) {
          std::cerr << "Error: CHECK TARGETS entries must be >= 0.\n";
          return true;
        }
        p += std::norm(state->get_amplitude(static_cast<Bitstring>(t)));
      }
      std::cout << "CHECK TARGETS: total_probability=" << p << "\n";
      return true;
    }
    if (mode == "BELL") {
      if (state->get_num_qubits() != 2) {
        std::cerr << "Error: CHECK BELL requires exactly 2 qubits.\n";
        return true;
      }
      if (tokens.size() < 3) {
        std::cerr << "Error: CHECK BELL requires PHI+, PHI-, PSI+, or PSI-.\n";
        return true;
      }
      const std::string bell = tokens[2];
      const ComplexNumber a00 = state->get_amplitude(0);
      const ComplexNumber a01 = state->get_amplitude(1);
      const ComplexNumber a10 = state->get_amplitude(2);
      const ComplexNumber a11 = state->get_amplitude(3);
      bool pass = false;
      if (bell == "PHI+" || bell == "PHI-") {
        const double p00 = std::norm(a00);
        const double p11 = std::norm(a11);
        const double support = std::norm(a01) + std::norm(a10);
        const double rel = std::real(a00 * std::conj(a11));
        pass = (p00 > 0.49 && p11 > 0.49 && support < 1e-6);
        if (bell == "PHI+") pass = pass && (rel > 0.0);
        else pass = pass && (rel < 0.0);
      } else if (bell == "PSI+" || bell == "PSI-") {
        const double p01 = std::norm(a01);
        const double p10 = std::norm(a10);
        const double support = std::norm(a00) + std::norm(a11);
        const double rel = std::real(a01 * std::conj(a10));
        pass = (p01 > 0.49 && p10 > 0.49 && support < 1e-6);
        if (bell == "PSI+") pass = pass && (rel > 0.0);
        else pass = pass && (rel < 0.0);
      } else {
        std::cerr << "Error: CHECK BELL requires PHI+, PHI-, PSI+, or PSI-.\n";
        return true;
      }
      std::cout << "CHECK BELL " << bell << ": " << (pass ? "PASS" : "FAIL")
                << (pass ? "" : " (hint: verify qubit order and phase gate placement)") << "\n";
      return true;
    }
    std::cerr << "Error: unknown CHECK mode. Use NORMALIZED, BASIS, TARGETS, or BELL.\n";
    return true;
  }
  return false;
}
