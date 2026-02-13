#include "cli/shell.hh"
#include "cli/commands.hh"
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
    state->display();
    return true;
  }
  if (cmd == "DISPLAY") {
    state->display();
    return true;
  }
  return false;
}
