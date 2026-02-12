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
#include "logging.hh"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>


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
      if (state != nullptr) delete state; // Cleanup old state
      state = new State(N, C);
      std::cout << "State initialized with " << N << " qubits and " << C << " classical register(s).\n";
    }
    return;
  }

  if (cmd == "QFTMODE") {
    if (tokens.size() < 2) {
      std::cerr << "Error: QFTMODE requires DIRECT or GATE.\n";
      return;
    }
    if (state == nullptr) {
      std::cerr << "Error: State must be initialized first. Use INIT <N> [C].\n";
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
      if (n <= 0 || n >= 63) {
        std::cerr << "Error: QUBO EXACT requires 1 <= n <= 62.\n";
        return;
      }
      const size_t expected = static_cast<size_t>(n) * static_cast<size_t>(n);
      if (tokens.size() != 3 + expected) {
        std::cerr << "Error: QUBO EXACT requires n*n matrix entries.\n";
        return;
      }
      std::vector<double> matrix;
      matrix.reserve(expected);
      for (size_t i = 0; i < expected; ++i) {
        double v = cli::get_double_arg(tokens, 3 + i, "QUBO EXACT");
        if (std::isnan(v)) return;
        matrix.push_back(v);
      }
      run_qubo_exact_cli(n, matrix);
      return;
    }
    if (mode == "GROVER") {
      int n = cli::get_arg(tokens, 2, "QUBO GROVER");
      double threshold = cli::get_double_arg(tokens, 3, "QUBO GROVER");
      int iterations = cli::get_arg(tokens, 4, "QUBO GROVER");
      if (n == -1 || iterations == -1 || std::isnan(threshold)) return;
      if (n <= 0 || n >= 63) {
        std::cerr << "Error: QUBO GROVER requires 1 <= n <= 62.\n";
        return;
      }
      const size_t expected = static_cast<size_t>(n) * static_cast<size_t>(n);
      if (tokens.size() != 5 + expected) {
        std::cerr << "Error: QUBO GROVER requires n*n matrix entries.\n";
        return;
      }
      std::vector<double> matrix;
      matrix.reserve(expected);
      for (size_t i = 0; i < expected; ++i) {
        double v = cli::get_double_arg(tokens, 5 + i, "QUBO GROVER");
        if (std::isnan(v)) return;
        matrix.push_back(v);
      }
      run_qubo_grover_cli(n, threshold, iterations, matrix);
      return;
    }
    std::cerr << "Error: QUBO mode must be DEMO, EXACT, or GROVER.\n";
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

  // Must be initialized before applying quantum operations
  if (state == nullptr || !state->is_initialized()) {
    std::cerr << "Error: State must be initialized first. Use INIT <N> [C].\n";
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

  // --- Algorithims ---
  if (cmd == "GROVER") {
    if (tokens.size() < 2) {
      std::cerr << "Error: GROVER requires at least one target.\n";
      return;
    }
    std::vector<Bitstring> targets;
    for (size_t i = 1; i < tokens.size(); ++i) {
      int t = cli::get_arg(tokens, i, "GROVER");
      if (t == -1) return;
      targets.push_back(static_cast<Bitstring>(t));
    }
    GroverResult result = run_grover(*state, targets);
    if (!result.ok) {
      std::cerr << "Grover error: " << result.error << "\n";
      return;
    }
    std::cout << "Grover iterations used: " << result.iterations << "\n";
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

// Destructor to clean up the dynamically allocated State object
QuantumShell::~QuantumShell()
{
  if (state != nullptr) delete state;
}

void QuantumShell::print_help()
{
  std::cout << "\n--- Quantum Simulator Commands ---\n";
  std::cout << "INIT <N> [C]     : Initialize state with N qubits and C classical bits.\n";
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
  std::cout << "MEASURE <j> <c>  : Measure qubit j, store result in classical register c.\n";
  std::cout << "DISPLAY          : Show the current quantum state.\n";
  std::cout << "GROVER <t...>    : Run Grover's algorithm searching for one or more targets\n";
  std::cout << "VERBOSE <level>  : Set verbosity (QUIET, NORMAL, VERBOSE or 0/1/2)\n";
  std::cout << "DEUTSCH_JOZSA <n> <oracle> : Run Deutsch-Jozsa demo (CONST0, CONST1, BALANCED_XOR0, BALANCED_PARITY)\n";
  std::cout << "BV <n> <secret> [bias]: Run Bernstein-Vazirani demo (bias defaults to 0)\n";
  std::cout << "QUBO DEMO        : Run built-in QUBO exact+Grover-threshold demo\n";
  std::cout << "QUBO EXACT <n> <n*n matrix entries> : Solve QUBO exactly by brute force\n";
  std::cout << "QUBO GROVER <n> <threshold> <iterations> <n*n matrix entries> : Grover threshold search\n";
  std::cout << "QFTMODE <mode>   : Set QFT mode (DIRECT or GATE, default DIRECT).\n";
  std::cout << "QRNG <n> [count] : Quantum random numbers from n qubits (count default 1)\n";
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
  std::cout << "SHOR <N>         : Run Shor's algorithm demo to factor N\n";
  std::cout << "HELP             : Show this help message.\n";
  std::cout << "QUIT             : Exit the simulator.\n";
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
