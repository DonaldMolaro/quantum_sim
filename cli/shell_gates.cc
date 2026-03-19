#include "cli/shell.hh"
#include "cli/commands.hh"
#include "internal/limits.hh"
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
      {"Z", &State::z}, {"S", &State::s}, {"T", &State::t},
      {"SDG", &State::sdg}, {"TDG", &State::tdg}};
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
  if (cmd == "P") {
    int j = cli::get_arg(tokens, 1, cmd);
    double phi = cli::get_angle_arg_required(tokens, 2, cmd);
    if (j == -1 || std::isnan(phi)) return true;
    try {
      state->p(j, phi);
      std::cout << "P(" << j << ", " << phi << ") applied.\n";
      tutor_note("P(phi) applies phase e^{i*phi} to |1> and identity to |0>.");
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "CP") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    double phi = cli::get_angle_arg_required(tokens, 3, cmd);
    if (j == -1 || k == -1 || std::isnan(phi)) return true;
    try {
      state->cp(j, k, phi);
      std::cout << "CP(" << j << ", " << k << ", " << phi << ") applied.\n";
      tutor_note("CP applies e^{i*phi} only when both control and target are |1>.");
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "ISWAP") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    if (j == -1 || k == -1) return true;
    try {
      state->iswap(j, k);
      std::cout << "iSWAP(" << j << ", " << k << ") applied.\n";
      tutor_note("iSWAP swaps j and k and multiplies the swapped amplitude by i.");
      state->display();
    } catch (const std::exception& e) { std::cerr << "Operation failed: " << e.what() << "\n"; }
    return true;
  }

  if (cmd == "XX" || cmd == "YY" || cmd == "ZZ") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    double theta = cli::get_angle_arg_required(tokens, 3, cmd);
    if (j == -1 || k == -1 || std::isnan(theta)) return true;
    try {
      if      (cmd == "XX") state->xx(j, k, theta);
      else if (cmd == "YY") state->yy(j, k, theta);
      else                  state->zz(j, k, theta);
      std::cout << cmd << "(" << j << ", " << k << ", " << theta << ") applied.\n";
      tutor_note(cmd + "(theta) is exp(-i theta/2 " + cmd + ") — native trapped-ion / superconducting gate.");
      state->display();
    } catch (const std::exception& e) { std::cerr << "Operation failed: " << e.what() << "\n"; }
    return true;
  }

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

  if (cmd == "CSWAP" || cmd == "FREDKIN") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    int l = cli::get_arg(tokens, 3, cmd);
    if (j == -1 || k == -1 || l == -1) return true;
    try {
      state->cswap(j, k, l);
      std::cout << cmd << "(" << j << ", " << k << ", " << l << ") applied.\n";
      tutor_note("CSWAP (Fredkin) swaps qubits k and l when control j is |1>.");
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return true;
  }

  if (cmd == "MCX") {
    // MCX c1 c2 ... cn t — all but last are controls, last is target
    if (tokens.size() < 3) {
      std::cerr << "Error: MCX requires at least one control and one target qubit.\n";
      return true;
    }
    std::vector<int> controls;
    for (size_t i = 1; i + 1 < tokens.size(); ++i) {
      int q = cli::get_arg(tokens, i, "MCX");
      if (q == -1) return true;
      controls.push_back(q);
    }
    int target = cli::get_arg(tokens, tokens.size() - 1, "MCX");
    if (target == -1) return true;
    try {
      state->mcx(controls, target);
      std::cout << "MCX applied (" << controls.size() << " controls -> target " << target << ").\n";
      tutor_note("MCX flips the target qubit only when all control qubits are |1>.");
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
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
                << (diff < qsim::limits::AMPLITUDE_EPSILON ? " [PASS]\n" : " [WARN]\n");
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

  // --- EXPECT <P1> <q1> [P2 q2 ...] ---
  if (cmd == "EXPECT") {
    if (tokens.size() < 3 || (tokens.size() - 1) % 2 != 0) {
      std::cerr << "Error: EXPECT requires pairs of <Pauli> <qubit> (e.g. EXPECT Z 0 Z 1).\n";
      return true;
    }
    std::vector<std::pair<char,int>> ops;
    for (size_t i = 1; i + 1 < tokens.size(); i += 2) {
      const std::string& pstr = tokens[i];
      if (pstr.size() != 1 || std::string("IXYZ").find(pstr[0]) == std::string::npos) {
        std::cerr << "Error: EXPECT Pauli must be one of I X Y Z.\n";
        return true;
      }
      int q = cli::get_arg(tokens, i + 1, "EXPECT");
      if (q == -1) return true;
      ops.push_back({pstr[0], q});
    }
    std::string obs;
    for (const auto& op : ops) obs += op.first + std::to_string(op.second) + " ";
    double ev = state->expect_pauli(ops);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "EXPECT " << obs << "= " << ev << "\n";
    tutor_note("Expectation values are real for Hermitian (Pauli) observables. Range is [-1, 1].");
    return true;
  }

  // --- FIDELITY <idx> ---
  if (cmd == "FIDELITY") {
    if (tokens.size() < 2) {
      std::cerr << "Error: FIDELITY requires a basis state index.\n";
      return true;
    }
    int idx = cli::get_arg(tokens, 1, "FIDELITY");
    if (idx == -1) return true;
    double fid = std::norm(state->get_amplitude(static_cast<Bitstring>(idx)));
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "FIDELITY with |" << idx << "> = " << fid
              << (fid > 0.999 ? "  [~1.000, nearly pure target]\n" : "\n");
    tutor_note("Fidelity |<idx|psi>|^2 = 1 means the state is exactly |idx>; 0 means orthogonal.");
    return true;
  }

  // --- BLOCH <j> ---
  if (cmd == "BLOCH") {
    int j = cli::get_arg(tokens, 1, "BLOCH");
    if (j == -1) return true;
    const State::BlochVector bv = state->bloch(j);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Bloch vector qubit " << j << ": x=" << bv.x
              << "  y=" << bv.y << "  z=" << bv.z << "\n";
    // Polar angles
    const double r = std::sqrt(bv.x*bv.x + bv.y*bv.y + bv.z*bv.z);
    if (r > 1e-9) {
      const double theta = std::acos(std::max(-1.0, std::min(1.0, bv.z / r))) * 180.0 / M_PI;
      const double phi   = std::atan2(bv.y, bv.x) * 180.0 / M_PI;
      std::cout << "  polar: r=" << r << "  theta=" << theta << "deg  phi=" << phi << "deg\n";
      std::cout << "  (r=1 => pure; r<1 => mixed due to entanglement)\n";
    } else {
      std::cout << "  (maximally mixed single-qubit state)\n";
    }
    tutor_note("Bloch z=+1 is |0>, z=-1 is |1>; equator is equal superposition; r<1 means entangled.");
    return true;
  }

  // --- ENTROPY <start_q> [end_q] ---
  if (cmd == "ENTROPY") {
    int start_q = cli::get_arg(tokens, 1, "ENTROPY");
    if (start_q == -1) return true;
    int end_q = (tokens.size() > 2) ? cli::get_arg(tokens, 2, "ENTROPY") : start_q;
    if (end_q == -1) return true;
    if (end_q - start_q >= 10) {
      std::cerr << "Error: ENTROPY subsystem can span at most 10 qubits.\n";
      return true;
    }
    double S = state->entropy(start_q, end_q);
    int n_A = std::abs(end_q - start_q) + 1;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "ENTROPY qubit" << (n_A > 1 ? "s " : " ") << start_q;
    if (n_A > 1) std::cout << ".." << end_q;
    std::cout << " = " << S << " bits  (max=" << n_A << ")\n";
    if (S < 1e-6)      std::cout << "  (product state — no entanglement)\n";
    else if (S > n_A - 0.01) std::cout << "  (maximally entangled subsystem)\n";
    tutor_note("S=0 means the subsystem is separable. S=1 for a Bell pair. S=log2(dim_A) is maximum.");
    return true;
  }

  // --- RESET <j> ---
  if (cmd == "RESET") {
    int j = cli::get_arg(tokens, 1, "RESET");
    if (j == -1) return true;
    state->reset(j);
    std::cout << "RESET qubit " << j << " to |0>.\n";
    tutor_note("RESET non-unitarily collapses qubit j to |0> — needed for ancilla reuse.");
    state->display();
    return true;
  }

  // --- IF <cbit> <gate_tokens...> ---
  if (cmd == "IF") {
    if (tokens.size() < 3) {
      std::cerr << "Error: IF requires a classical bit index and a gate command.\n";
      return true;
    }
    int cbit = cli::get_arg(tokens, 1, "IF");
    if (cbit == -1) return true;
    int val = 0;
    try { val = state->get_cbit(static_cast<size_t>(cbit)); }
    catch (const std::out_of_range&) {
      std::cerr << "Error: classical bit " << cbit << " is out of range.\n";
      return true;
    }
    if (val == 1) {
      // Build sub-command tokens from tokens[2..]
      std::vector<std::string> sub_tokens(tokens.begin() + 2, tokens.end());
      tutor_note("IF: cbit=" + std::to_string(cbit) + " is 1 — applying conditional gate.");
      handle_command(sub_tokens);
    } else {
      std::cout << "IF: cbit " << cbit << " = 0 — gate skipped.\n";
    }
    return true;
  }

  // --- SWAP_TEST <ancilla> <j_start> <k_start> <n_qubits> ---
  if (cmd == "SWAP_TEST") {
    if (tokens.size() < 5) {
      std::cerr << "Error: SWAP_TEST requires ancilla, first_reg_start, second_reg_start, n.\n";
      std::cerr << "  e.g.  SWAP_TEST 0 1 2 1   (ancilla=0, reg_A=[1..1], reg_B=[2..2])\n";
      return true;
    }
    int anc     = cli::get_arg(tokens, 1, "SWAP_TEST");
    int a_start = cli::get_arg(tokens, 2, "SWAP_TEST");
    int b_start = cli::get_arg(tokens, 3, "SWAP_TEST");
    int n       = cli::get_arg(tokens, 4, "SWAP_TEST");
    if (anc == -1 || a_start == -1 || b_start == -1 || n == -1 || n < 1) return true;
    // H on ancilla
    state->h(anc);
    // CSWAP for each pair
    for (int i = 0; i < n; ++i)
      state->cswap(anc, a_start + i, b_start + i);
    // H on ancilla again
    state->h(anc);
    // Probability P(ancilla=0) = (1 + overlap) / 2  =>  overlap = 2*P(0) - 1
    double p0 = state->compute_probability_of_0(anc);
    double overlap = 2.0 * p0 - 1.0;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "SWAP_TEST: P(ancilla=0)=" << p0
              << "  => |<psi|phi>|^2 ~= " << std::max(0.0, overlap) << "\n";
    tutor_note("SWAP test: if |psi>=|phi> overlap=1 so P(0)=1; if orthogonal overlap=0 so P(0)=0.5.");
    state->display();
    return true;
  }

  return false;
}
