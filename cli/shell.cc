/*
 * This is a simple little shell that reads commands and manipulates
 * the underlying quantum state simulator. It's not bad, it's not great
 * but with it you can poke at the machine and see what it does.
 *
 * Beyond just basic quantum gates it's also able to call an implementation of
 * Grover's algorithm. Maybe I'll get ambitious and implement a few other
 * quantum algorithms.
 */
#include "state.hh"
#include "cli/commands.hh"
#include "cli/shell.hh"
#include "algorithms/api/grover_api.hh"
#include "algorithms/qrng.hh"
#include "demos/bernstein_vazirani_demo.hh"
#include "demos/deutsch_jozsa_demo.hh"
#include "demos/grover_demo.hh"
#include "demos/latin_demo.hh"
#include "demos/qubo_demo.hh"
#include "demos/shor_demo.hh"
#include "demos/vqa_demo.hh"
#include "demos/vqe_demo.hh"
#include "demos/anneal_demo.hh"
#include "logging.hh"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool require_initialized_state(const std::unique_ptr<State>& state)
{
  if (!state || !state->is_initialized()) {
    std::cerr << "Error: State must be initialized first. Use INIT <N> [C].\n";
    return false;
  }
  return true;
}

bool parse_square_matrix(const std::vector<std::string>& tokens,
                         size_t start_index,
                         int n,
                         const std::string& cmd,
                         std::vector<double>& out)
{
  if (n <= 0 || n >= 63) {
    std::cerr << "Error: " << cmd << " requires 1 <= n <= 62.\n";
    return false;
  }
  const size_t expected = static_cast<size_t>(n) * static_cast<size_t>(n);
  if (tokens.size() != start_index + expected) {
    std::cerr << "Error: " << cmd << " requires n*n matrix entries.\n";
    return false;
  }
  out.clear();
  out.reserve(expected);
  for (size_t i = 0; i < expected; ++i) {
    const double v = cli::get_double_arg(tokens, start_index + i, cmd);
    if (std::isnan(v)) {
      out.clear();
      return false;
    }
    out.push_back(v);
  }
  return true;
}

} // namespace


static std::string bits_to_string(const std::vector<int>& bits)
{
  if (bits.empty()) return "0";
  std::string out;
  out.reserve(bits.size());
  for (int i = static_cast<int>(bits.size()) - 1; i >= 0; --i) {
    out.push_back(bits[i] ? '1' : '0');
  }
  return out;
}

static uint64_t bits_to_u64(const std::vector<int>& bits)
{
  uint64_t value = 0;
  for (size_t i = 0; i < bits.size() && i < 63; ++i) {
    if (bits[i]) {
      value |= (1ULL << i);
    }
  }
  return value;
}

static int infer_qubits_from_targets(const std::vector<Bitstring>& targets)
{
  Bitstring max_target = 0ULL;
  for (size_t i = 0; i < targets.size(); ++i) {
    if (targets[i] > max_target) {
      max_target = targets[i];
    }
  }
  if (max_target == 0ULL) {
    return 1;
  }
  int n = 0;
  while (max_target != 0ULL) {
    ++n;
    max_target >>= 1;
  }
  return std::max(1, n);
}

static const char* log_level_name(qsim_log::Level level)
{
  switch (level) {
    case qsim_log::Level::Quiet:
      return "QUIET";
    case qsim_log::Level::Normal:
      return "NORMAL";
    case qsim_log::Level::Verbose:
      return "VERBOSE";
  }
  return "UNKNOWN";
}

void QuantumShell::handle_command(const std::vector<std::string>& tokens)
{
  if (tokens.empty()) return;

  const std::string cmd = tokens[0];

  // --- Initialization ---
  if (cmd == "INIT") {
    int N = cli::get_arg(tokens, 1, "INIT");
    int C = (tokens.size() > 2) ? cli::get_arg(tokens, 2, "INIT") : 0;

    if (N > 0 && N <= 64) { // Arbitrary limit for Bitstring (ULL)
      state.reset(new State(N, C));
      std::cout << "State initialized with " << N << " qubits and " << C << " classical register(s).\n";
    }
    return;
  }

  if (cmd == "QFTMODE") {
    if (tokens.size() < 2) {
      std::cerr << "Error: QFTMODE requires DIRECT or GATE.\n";
      return;
    }
    if (!require_initialized_state(state)) {
      return;
    }
    if (tokens[1] == "DIRECT") {
      state->set_qft_mode(State::QftMode::Direct);
      std::cout << "QFT mode set to DIRECT.\n";
    } else if (tokens[1] == "GATE") {
      state->set_qft_mode(State::QftMode::Gate);
      std::cout << "QFT mode set to GATE.\n";
    } else {
      std::cerr << "Error: QFTMODE must be DIRECT or GATE.\n";
    }
    return;
  }

  if (cmd == "QRNG") {
    int n = cli::get_arg(tokens, 1, "QRNG");
    int count = (tokens.size() > 2) ? cli::get_arg(tokens, 2, "QRNG") : 1;
    if (n <= 0) {
      std::cerr << "Error: QRNG requires n > 0.\n";
      return;
    }
    if (count <= 0) {
      std::cerr << "Error: QRNG count must be > 0.\n";
      return;
    }
    if (n > 63) {
      std::cout << "Note: QRNG displays integer values using only the lowest 63 bits.\n";
    }

    for (int i = 0; i < count; ++i) {
      std::vector<int> bits = qrng_bits(n);
      std::cout << "QRNG[" << i << "] bits=" << bits_to_string(bits)
                << " value=" << bits_to_u64(bits) << "\n";
    }
    return;
  }

  if (cmd == "VERBOSE" || cmd == "LOGLEVEL") {
    if (tokens.size() < 2) {
      std::cout << "Verbosity: " << log_level_name(qsim_log::get_level()) << "\n";
      return;
    }
    qsim_log::Level level;
    if (!qsim_log::parse_level(tokens[1], level)) {
      std::cerr << "Error: VERBOSE expects QUIET, NORMAL, VERBOSE, or 0/1/2.\n";
      return;
    }
    qsim_log::set_level(level);
    std::cout << "Verbosity set to " << log_level_name(level) << ".\n";
    return;
  }

  if (cmd == "DEUTSCH_JOZSA" || cmd == "DEUTSCH") {
    int n_inputs = cli::get_arg(tokens, 1, cmd);
    if (n_inputs == -1) return;
    if (tokens.size() < 3) {
      std::cerr << "Error: " << cmd << " requires an oracle (CONST0, CONST1, BALANCED_XOR0, BALANCED_PARITY).\n";
      return;
    }
    run_deutsch_jozsa_demo(n_inputs, tokens[2]);
    return;
  }

  if (cmd == "BV" || cmd == "BERNSTEIN_VAZIRANI") {
    int n_inputs = cli::get_arg(tokens, 1, cmd);
    int secret = cli::get_arg(tokens, 2, cmd);
    int bias = (tokens.size() > 3) ? cli::get_arg(tokens, 3, cmd) : 0;
    if (n_inputs == -1 || secret == -1 || bias == -1) return;
    run_bernstein_vazirani_demo(n_inputs,
                                static_cast<Bitstring>(secret),
                                bias);
    return;
  }

  if (cmd == "QUBO") {
    if (tokens.size() < 2) {
      std::cerr << "Error: QUBO requires a mode (DEMO, EXACT, GROVER).\n";
      return;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_qubo_demo();
      return;
    }
    if (mode == "EXACT") {
      int n = cli::get_arg(tokens, 2, "QUBO EXACT");
      if (n == -1) return;
      std::vector<double> matrix;
      if (!parse_square_matrix(tokens, 3, n, "QUBO EXACT", matrix)) return;
      run_qubo_exact_cli(n, matrix);
      return;
    }
    if (mode == "GROVER") {
      int n = cli::get_arg(tokens, 2, "QUBO GROVER");
      double threshold = cli::get_double_arg(tokens, 3, "QUBO GROVER");
      int iterations = cli::get_arg(tokens, 4, "QUBO GROVER");
      if (n == -1 || iterations == -1 || std::isnan(threshold)) return;
      std::vector<double> matrix;
      if (!parse_square_matrix(tokens, 5, n, "QUBO GROVER", matrix)) return;
      run_qubo_grover_cli(n, threshold, iterations, matrix);
      return;
    }
    std::cerr << "Error: QUBO mode must be DEMO, EXACT, or GROVER.\n";
    return;
  }

  if (cmd == "VQA") {
    if (tokens.size() < 2) {
      std::cerr << "Error: VQA requires a mode (DEMO, QAOA).\n";
      return;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_vqa_demo();
      return;
    }
    if (mode == "QAOA") {
      int n = cli::get_arg(tokens, 2, "VQA QAOA");
      int p_layers = cli::get_arg(tokens, 3, "VQA QAOA");
      int shots = cli::get_arg(tokens, 4, "VQA QAOA");
      int iters = cli::get_arg(tokens, 5, "VQA QAOA");
      double step = cli::get_double_arg(tokens, 6, "VQA QAOA");
      if (n == -1 || p_layers == -1 || shots == -1 || iters == -1 || std::isnan(step)) return;
      std::vector<double> matrix;
      if (!parse_square_matrix(tokens, 7, n, "VQA QAOA", matrix)) return;
      run_vqa_qaoa_cli(n, p_layers, shots, iters, step, matrix);
      return;
    }
    std::cerr << "Error: VQA mode must be DEMO or QAOA.\n";
    return;
  }

  if (cmd == "VQE") {
    if (tokens.size() < 2) {
      std::cerr << "Error: VQE requires a mode (DEMO, RUN).\n";
      return;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_vqe_demo();
      return;
    }
    if (mode == "RUN") {
      int n = cli::get_arg(tokens, 2, "VQE RUN");
      int layers = cli::get_arg(tokens, 3, "VQE RUN");
      int iters = cli::get_arg(tokens, 4, "VQE RUN");
      double step = cli::get_double_arg(tokens, 5, "VQE RUN");
      int shots = cli::get_arg(tokens, 6, "VQE RUN");
      int term_count = cli::get_arg(tokens, 7, "VQE RUN");
      if (n == -1 || layers == -1 || iters == -1 || shots == -1 || term_count == -1 || std::isnan(step)) {
        return;
      }
      if (term_count <= 0) {
        std::cerr << "Error: VQE RUN requires term_count > 0.\n";
        return;
      }

      VqeHamiltonian h;
      h.n_qubits = n;
      size_t idx = 8;
      for (int t = 0; t < term_count; ++t) {
        double coeff = cli::get_double_arg(tokens, idx, "VQE RUN");
        if (std::isnan(coeff)) return;
        ++idx;

        int op_count = cli::get_arg(tokens, idx, "VQE RUN");
        if (op_count == -1) return;
        if (op_count < 0) {
          std::cerr << "Error: VQE RUN term op_count must be >= 0.\n";
          return;
        }
        ++idx;

        VqePauliTerm term;
        term.coeff = coeff;
        for (int k = 0; k < op_count; ++k) {
          if (idx >= tokens.size()) {
            std::cerr << "Error: VQE RUN missing Pauli op token.\n";
            return;
          }
          const std::string op_token = tokens[idx++];
          if (op_token.size() != 1) {
            std::cerr << "Error: VQE RUN Pauli op must be one of X,Y,Z.\n";
            return;
          }
          int q = cli::get_arg(tokens, idx, "VQE RUN");
          if (q == -1) return;
          ++idx;
          term.ops.push_back(VqePauliOp{op_token[0], q});
        }
        h.terms.push_back(term);
      }

      if (idx != tokens.size()) {
        std::cerr << "Error: VQE RUN has extra trailing tokens.\n";
        return;
      }

      run_vqe_cli(h, layers, iters, step, shots);
      return;
    }
    std::cerr << "Error: VQE mode must be DEMO or RUN.\n";
    return;
  }

  if (cmd == "ANNEAL") {
    if (tokens.size() < 2) {
      std::cerr << "Error: ANNEAL requires a mode (DEMO, QUBO).\n";
      return;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_anneal_demo();
      return;
    }
    if (mode == "QUBO") {
      if (tokens.size() < 9) {
        std::cerr << "Error: ANNEAL QUBO requires method, n, steps, sweeps, beta_start, beta_end, replicas, matrix.\n";
        return;
      }
      const std::string method = tokens[2];
      int n = cli::get_arg(tokens, 3, "ANNEAL QUBO");
      int steps = cli::get_arg(tokens, 4, "ANNEAL QUBO");
      int sweeps = cli::get_arg(tokens, 5, "ANNEAL QUBO");
      double beta_start = cli::get_double_arg(tokens, 6, "ANNEAL QUBO");
      double beta_end = cli::get_double_arg(tokens, 7, "ANNEAL QUBO");
      int replicas = cli::get_arg(tokens, 8, "ANNEAL QUBO");
      if (n == -1 || steps == -1 || sweeps == -1 || replicas == -1 ||
          std::isnan(beta_start) || std::isnan(beta_end)) {
        return;
      }
      std::vector<double> matrix;
      if (!parse_square_matrix(tokens, 9, n, "ANNEAL QUBO", matrix)) return;
      run_anneal_qubo_cli(method, n, steps, sweeps, beta_start, beta_end, replicas, matrix);
      return;
    }
    std::cerr << "Error: ANNEAL mode must be DEMO or QUBO.\n";
    return;
  }

  if (cmd == "LATIN") {
    size_t idx = 1;
    std::string mode = "demo";
    if (idx < tokens.size()) {
      if (tokens[idx] == "DEMO" || tokens[idx] == "COUNT" || tokens[idx] == "PRINT-ALL") {
        mode = tokens[idx];
        ++idx;
      }
    }

    int iters = -1;
    std::vector<int> row0;

    if (idx < tokens.size()) {
      int maybe_iters = cli::get_arg(tokens, idx, "LATIN");
      if (maybe_iters == -1) return;
      if (mode == "DEMO") {
        iters = maybe_iters;
        ++idx;
      }
    }

    while (idx < tokens.size()) {
      int v = cli::get_arg(tokens, idx, "LATIN");
      if (v == -1) return;
      row0.push_back(v);
      ++idx;
    }

    int row_vals[3] = {0, 1, 2};
    if (!row0.empty()) {
      if (row0.size() != 3) {
        std::cerr << "LATIN row0 must have exactly 3 values.\n";
        return;
      }
      row_vals[0] = row0[0];
      row_vals[1] = row0[1];
      row_vals[2] = row0[2];
    }

    if (mode == "COUNT") {
      run_latin3_count_row0(row_vals);
      return;
    }
    if (mode == "PRINT-ALL") {
      run_latin3_print_all_row0(row_vals);
      return;
    }

    if (!row0.empty() || iters >= 0) {
      run_latin3_grover_demo_row0(row_vals, iters);
    } else {
      run_latin3_grover_demo(-1);
    }
    return;
  }
  if (cmd == "SHOR") {
    int N = cli::get_arg(tokens, 1, "SHOR");
    if (N == -1) return;
    run_shor_demo(static_cast<Bitstring>(N));
    return;
  }

  // --- Algorithms that self-initialize their own state ---
  if (cmd == "GROVER") {
    if (tokens.size() < 2) {
      std::cerr << "Error: GROVER requires at least one target.\n";
      return;
    }
    std::vector<Bitstring> targets;
    targets.reserve(tokens.size() - 1);
    for (size_t i = 1; i < tokens.size(); ++i) {
      int t = cli::get_arg(tokens, i, "GROVER");
      if (t == -1) return;
      if (t < 0) {
        std::cerr << "Error: GROVER targets must be >= 0.\n";
        return;
      }
      targets.push_back(static_cast<Bitstring>(t));
    }

    int n_qubits = infer_qubits_from_targets(targets);
    if (state && state->is_initialized()) {
      n_qubits = state->get_num_qubits();
    }

    GroverResult result = run_grover(n_qubits, targets);
    if (!result.ok) {
      std::cerr << "Grover error: " << result.error << "\n";
      return;
    }
    std::cout << "Grover iterations used: " << result.iterations
              << " (n=" << n_qubits << ")\n";
    return;
  }

  // Must be initialized before applying quantum operations
  if (!require_initialized_state(state)) {
    return;
  }

  // --- Gate Operations (Single Qubit) ---
  if (cmd == "H" || cmd == "X" || cmd == "Y" || cmd == "Z" || cmd == "S" || cmd == "T") {
    int j = cli::get_arg(tokens, 1, cmd);
    if (j == -1) return; 

    try {
      if (cmd == "H") state->h(j); // Implements flatMap + reduceByKey
      else if (cmd == "X") state->x(j); // Implements s.map(λb, a. (b¬j , a))
      else if (cmd == "Y") state->y(j); // Implements RZ(pi/2) X RZ(-pi/2)
      else if (cmd == "Z") state->z(j); // Implements RZ(pi)
      else if (cmd == "S") state->s(j); // Implements s.map(λb, a. (b, aibj ))
      else if (cmd == "T") state->t(j); // Implements s.map(λb, a. (b, a((1+i)/sqrt(2))^bj ))
      std::cout << cmd << "(" << j << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return;
  }

  if (cmd == "RX" || cmd == "RY" || cmd == "RZ") {
    int j = cli::get_arg(tokens, 1, cmd);
    double theta = cli::get_angle_arg_required(tokens, 2, cmd);
    if (j == -1 || std::isnan(theta)) return;

    try {
      if (cmd == "RX") state->rx(j, theta);
      else if (cmd == "RY") state->ry(j, theta);
      else state->rz(j, theta);
      std::cout << cmd << "(" << j << ", " << theta << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return;
  }

  if (cmd == "RU") {
    int j = cli::get_arg(tokens, 1, cmd);
    double theta = cli::get_angle_arg_required(tokens, 2, cmd);
    double phi = cli::get_angle_arg_required(tokens, 3, cmd);
    double lambda = cli::get_angle_arg_required(tokens, 4, cmd);
    if (j == -1 || std::isnan(theta) || std::isnan(phi) || std::isnan(lambda)) return;

    try {
      state->ru(j, theta, phi, lambda);
      std::cout << "RU(" << j << ", " << theta << ", " << phi << ", " << lambda << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return;
  }

  if (cmd == "CRZ" || cmd == "CRX" || cmd == "CRY") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    double theta = cli::get_angle_arg_required(tokens, 3, cmd);
    if (j == -1 || k == -1 || std::isnan(theta)) return;

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
    return;
  }

  if (cmd == "CU") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    double theta = cli::get_angle_arg_required(tokens, 3, cmd);
    double phi = cli::get_angle_arg_required(tokens, 4, cmd);
    double lambda = cli::get_angle_arg_required(tokens, 5, cmd);
    if (j == -1 || k == -1 || std::isnan(theta) || std::isnan(phi) || std::isnan(lambda)) return;

    try {
      state->cu(j, k, theta, phi, lambda);
      std::cout << "CU(" << j << ", " << k << ", " << theta << ", " << phi << ", " << lambda << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return;
  }

  // --- Gate Operations (Two Qubits) ---
  if (cmd == "CX" || cmd == "CZ" || cmd == "CY" || cmd == "CH" || cmd == "CNOT") {
    int j = cli::get_arg(tokens, 1, cmd);
    int k = cli::get_arg(tokens, 2, cmd);
    if (j == -1 || k == -1) return;

    if (cmd == "CX") {
      state->cx(j, k); // Implements s.map(λb, a. (ite(bj , b ¬k, b), a))
      std::cout << "CX(" << j << ", " << k << ") applied.\n";
    } else if (cmd == "CZ") {
      state->cz(j, k);
      std::cout << "CZ(" << j << ", " << k << ") applied.\n";
    } else if (cmd == "CY") {
      state->cy(j, k);
      std::cout << "CY(" << j << ", " << k << ") applied.\n";
    } else if (cmd == "CH") {
      state->ch(j, k);
      std::cout << "CH(" << j << ", " << k << ") applied.\n";
    } else {
      state->cnot(j, k);
      std::cout << "CNOT(" << j << ", " << k << ") applied.\n";
    }
    state->display();
    return;
  }

  if (cmd == "CCX" || cmd == "TOFFOLI") {
    int c1 = cli::get_arg(tokens, 1, cmd);
    int c2 = cli::get_arg(tokens, 2, cmd);
    int t = cli::get_arg(tokens, 3, cmd);
    if (c1 == -1 || c2 == -1 || t == -1) return;

    try {
      state->ccx(c1, c2, t);
      std::cout << cmd << "(" << c1 << ", " << c2 << ", " << t << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return;
  }
  if (cmd == "SWAP") {
    int j = cli::get_arg(tokens, 1, "SWAP");
    int k = cli::get_arg(tokens, 2, "SWAP");
    if (j == -1 || k == -1) return;

    state->swap(j, k); 
    std::cout << "SWAP(" << j << ", " << k << ") applied.\n";
    state->display();
    return;
  }

  // --- Measurement ---
  if (cmd == "MEASURE") {
    int j = cli::get_arg(tokens, 1, "MEASURE"); // Qubit index
    int c = cli::get_arg(tokens, 2, "MEASURE"); // Classical register index
    if (j == -1 || c == -1) return;

    state->measure(j, c); // Probability computation, collapse, and renormalization
    std::cout << "Qubit " << j << " measured.";    
    if (!state->get_cbits().empty())
      {
	std::cout << " Result stored in c[" << c << "]: " << state->get_cbits()[c] << "\n";
      }
    else
      {
	std::cout << " No Register to store in \n";
      }
    state->display();
    return;
  }

  // --- Display ---
  if (cmd == "DISPLAY") {
    state->display();
    return;
  }

  std::cout << "Unknown command: " << cmd << ". Use HELP for a list of commands.\n";
}

void QuantumShell::print_help()
{
  std::cout << "\n--- Quantum Simulator Commands ---\n";
  std::cout << "\n[Initialization]\n";
  std::cout << "INIT <N> [C]     : Initialize state with N qubits and C classical bits.\n";
  std::cout << "\n[Simple Gates]\n";
  std::cout << "H <j>            : Hadamard gate on qubit j (Superposition).\n";
  std::cout << "X <j>            : NOT gate on qubit j.\n";
  std::cout << "Y <j>            : Y gate on qubit j.\n";
  std::cout << "Z <j>            : Z gate on qubit j (Phase flip).\n";
  std::cout << "S <j>            : S gate on qubit j (Phase).\n";
  std::cout << "T <j>            : T gate on qubit j (Phase).\n";
  std::cout << "RX <j> <theta>   : Rotation around X by angle theta (radians by default).\n";
  std::cout << "RY <j> <theta>   : Rotation around Y by angle theta (radians by default).\n";
  std::cout << "RZ <j> <theta>   : Rotation around Z by angle theta (radians by default).\n";
  std::cout << "RU <j> <t> <p> <l>: General single-qubit gate U(t, p, l) (radians by default).\n";
  std::cout << "  Degrees: append DEG or add DEG token. Examples: RX 0 90 DEG, RY 1 45DEG\n";
  std::cout << "  Constants: PI, PI/2, PI/4, TAU (case-insensitive).\n";
  std::cout << "  Expressions: nPI, nPI/m (n,m numeric). Examples: RX 0 3PI/2, RZ 1 -PI/4\n";
  std::cout << "\n[Multi-Qubit Gates]\n";
  std::cout << "CX <j> <k>       : Controlled-X (CNOT) where j is control, k is target.\n";
  std::cout << "CNOT <j> <k>     : Controlled-X alias.\n";
  std::cout << "CZ <j> <k>       : Controlled-Z where j is control, k is target.\n";
  std::cout << "CY <j> <k>       : Controlled-Y where j is control, k is target.\n";
  std::cout << "CH <j> <k>       : Controlled-H where j is control, k is target.\n";
  std::cout << "CRZ <j> <k> <t>  : Controlled-RZ with angle t (radians by default).\n";
  std::cout << "CRX <j> <k> <t>  : Controlled-RX with angle t (radians by default).\n";
  std::cout << "CRY <j> <k> <t>  : Controlled-RY with angle t (radians by default).\n";
  std::cout << "CU <j> <k> <t> <p> <l>: Controlled-U with angles (radians by default).\n";
  std::cout << "CCX <c1> <c2> <t>: Toffoli (double-controlled X).\n";
  std::cout << "TOFFOLI <c1> <c2> <t>: Alias for CCX.\n";
  std::cout << "SWAP <j> <k>     : SWAP qubits j and, k maintaining amplitudes.\n";
  std::cout << "\n[Algorithms]\n";
  std::cout << "GROVER <t...>    : Run Grover search (self-initialized; no INIT required)\n";
  std::cout << "DEUTSCH_JOZSA <n> <oracle> : Run Deutsch-Jozsa demo (CONST0, CONST1, BALANCED_XOR0, BALANCED_PARITY)\n";
  std::cout << "BV <n> <secret> [bias]: Run Bernstein-Vazirani demo (bias defaults to 0)\n";
  std::cout << "SHOR <N>         : Run Shor's algorithm demo to factor N\n";
  std::cout << "QUBO DEMO        : Run built-in QUBO exact+Grover-threshold demo\n";
  std::cout << "QUBO EXACT <n> <n*n matrix entries> : Solve QUBO exactly by brute force\n";
  std::cout << "QUBO GROVER <n> <threshold> <iterations> <n*n matrix entries> : Grover threshold search\n";
  std::cout << "VQA DEMO         : Run built-in QAOA-on-QUBO demo\n";
  std::cout << "VQA QAOA <n> <p> <shots> <iters> <step> <n*n matrix entries> : QAOA optimize QUBO\n";
  std::cout << "VQE DEMO         : Run built-in VQE demo (H=-Z0)\n";
  std::cout << "VQE RUN <n> <layers> <iters> <step> <shots> <terms> <coeff op_count [op qubit]...>\n";
  std::cout << "ANNEAL DEMO      : Run built-in simulated quantum annealing demo\n";
  std::cout << "ANNEAL QUBO <SA|SQA> <n> <steps> <sweeps> <beta_start> <beta_end> <replicas> <n*n matrix entries>\n";
  std::cout << "LATIN [iters]               : Grover demo for 3x3 Latin squares (row0 fixed 0 1 2)\n";
  std::cout << "LATIN DEMO [iters] [r0 r1 r2]: Demo with custom row0 permutation\n";
  std::cout << "LATIN COUNT [r0 r1 r2]       : Count solutions for row0 permutation\n";
  std::cout << "LATIN PRINT-ALL [r0 r1 r2]   : Print all solutions for row0 permutation\n";
  std::cout << "  Examples:\n";
  std::cout << "    LATIN\n";
  std::cout << "    LATIN 8\n";
  std::cout << "    LATIN DEMO 6 0 2 1\n";
  std::cout << "    LATIN COUNT 2 0 1\n";
  std::cout << "    LATIN PRINT-ALL 1 2 0\n";
  std::cout << "\n[Utility]\n";
  std::cout << "MEASURE <j> <c>  : Measure qubit j, store result in classical register c.\n";
  std::cout << "DISPLAY          : Show the current quantum state.\n";
  std::cout << "HELP             : Show this help message.\n";
  std::cout << "QUIT             : Exit the simulator.\n";
  std::cout << "VERBOSE <level>  : Set verbosity (QUIET, NORMAL, VERBOSE or 0/1/2)\n";
  std::cout << "QFTMODE <mode>   : Set QFT mode (DIRECT or GATE, default DIRECT).\n";
  std::cout << "QRNG <n> [count] : Quantum random numbers from n qubits (count default 1)\n";
  std::cout << "----------------------------------\n";
}

void QuantumShell::run()
{
  State::set_default_log_stream(&std::cout);
  qsim_log::set_stream(&std::cout);
  print_help();
  std::string input_line;
  while (true) {
    std::cout << "QSIM> ";
    if (!std::getline(std::cin, input_line) || input_line == "QUIT" || input_line == "quit") {
      break;
    }
    if (input_line.empty()) continue;

    // Convert input to uppercase for case-insensitive matching
    std::transform(input_line.begin(), input_line.end(), input_line.begin(), ::toupper);
            
    if (input_line == "HELP") {
      print_help();
      continue;
    }
            
    handle_command(cli::parse_command(input_line));
  }
  std::cout << "Exiting simulator.\n";
}
