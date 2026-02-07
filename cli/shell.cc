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
#include "cli/shell.hh"
#include "algorithms/api/grover_api.hh"
#include "algorithms/qrng.hh"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

extern void run_grover_search(State *s,Bitstring targer_w);
extern void run_latin3_grover_demo(int iterations);
extern void run_latin3_grover_demo_row0(const int row0[3], int iterations);
extern void run_latin3_count_row0(const int row0[3]);
extern void run_latin3_print_all_row0(const int row0[3]);
extern void run_shor_demo(Bitstring N);

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

std::vector<std::string> QuantumShell::parse_command(const std::string& line)
{
  std::vector<std::string> tokens;
  std::stringstream ss(line);
  std::string token;
  while (ss >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

// Helper to get integer arguments safely
int QuantumShell::get_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return -1;
  }
  try {
    return std::stoi(tokens[index]);
  } catch (const std::exception&) {
    std::cerr << "Error: Argument '" << tokens[index] << "' must be an integer.\n";
    return -1;
  }
}

double QuantumShell::get_double_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
  try {
    return std::stod(tokens[index]);
  } catch (const std::exception&) {
    std::cerr << "Error: Argument '" << tokens[index] << "' must be a number.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
}

double QuantumShell::get_angle_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::string token = tokens[index];
  bool is_deg = false;

  if (token.size() >= 3 && token.substr(token.size() - 3) == "DEG") {
    token = token.substr(0, token.size() - 3);
    is_deg = true;
  } else if (index + 1 < tokens.size() && tokens[index + 1] == "DEG") {
    is_deg = true;
  }

  std::string token_upper = token;
  std::transform(token_upper.begin(), token_upper.end(), token_upper.begin(), ::toupper);

  const double pi = std::acos(-1.0);
  double theta = std::numeric_limits<double>::quiet_NaN();

  if (token_upper == "PI") {
    theta = pi;
    is_deg = false;
  } else if (token_upper == "TAU") {
    theta = 2.0 * pi;
    is_deg = false;
  } else if (token_upper == "PI/2") {
    theta = pi / 2.0;
    is_deg = false;
  } else if (token_upper == "PI/4") {
    theta = pi / 4.0;
    is_deg = false;
  } else if (token_upper == "-PI") {
    theta = -pi;
    is_deg = false;
  } else if (token_upper == "-PI/2") {
    theta = -pi / 2.0;
    is_deg = false;
  } else if (token_upper == "-PI/4") {
    theta = -pi / 4.0;
    is_deg = false;
  } else {
    bool parsed_pi_expr = false;
    bool neg = false;
    std::string expr = token_upper;
    if (!expr.empty() && expr[0] == '-') {
      neg = true;
      expr = expr.substr(1);
    }

    size_t pi_pos = expr.find("PI");
    if (pi_pos != std::string::npos) {
      std::string coeff_str = expr.substr(0, pi_pos);
      if (!coeff_str.empty()) {
        try {
          double coeff = std::stod(coeff_str);
          theta = coeff * pi;
          parsed_pi_expr = true;
        } catch (const std::exception&) {
          parsed_pi_expr = false;
        }
      } else {
        theta = pi;
        parsed_pi_expr = true;
      }

      if (parsed_pi_expr) {
        size_t denom_pos = expr.find('/', pi_pos + 2);
        if (denom_pos != std::string::npos) {
          std::string denom_str = expr.substr(denom_pos + 1);
          try {
            double denom = std::stod(denom_str);
            if (denom == 0.0) {
              std::cerr << "Error: division by zero in '" << tokens[index] << "'.\n";
              return std::numeric_limits<double>::quiet_NaN();
            }
            theta /= denom;
          } catch (const std::exception&) {
            std::cerr << "Error: Invalid PI expression '" << tokens[index] << "'.\n";
            return std::numeric_limits<double>::quiet_NaN();
          }
        }
      }
    }

    if (parsed_pi_expr) {
      if (neg) theta = -theta;
    } else {
      try {
        theta = std::stod(token);
      } catch (const std::exception&) {
        std::cerr << "Error: Argument '" << tokens[index] << "' must be a number.\n";
        return std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  if (is_deg) {
    theta = theta * (pi / 180.0);
  }

  return theta;
}

double QuantumShell::get_angle_arg_required(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
  return get_angle_arg(tokens, index, cmd);
}

void QuantumShell::handle_command(const std::vector<std::string>& tokens)
{
  if (tokens.empty()) return;

  const std::string cmd = tokens[0];

  // --- Initialization ---
  if (cmd == "INIT") {
    int N = get_arg(tokens, 1, "INIT");
    int C = (tokens.size() > 2) ? get_arg(tokens, 2, "INIT") : 0;

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
    int n = get_arg(tokens, 1, "QRNG");
    int count = (tokens.size() > 2) ? get_arg(tokens, 2, "QRNG") : 1;
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

  // Must be initialized before applying quantum operations
  if (state == nullptr || !state->is_initialized()) {
    std::cerr << "Error: State must be initialized first. Use INIT <N> [C].\n";
    return;
  }

  // --- Gate Operations (Single Qubit) ---
  if (cmd == "H" || cmd == "X" || cmd == "Y" || cmd == "Z" || cmd == "S" || cmd == "T") {
    int j = get_arg(tokens, 1, cmd);
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
    int j = get_arg(tokens, 1, cmd);
    double theta = get_angle_arg_required(tokens, 2, cmd);
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
    int j = get_arg(tokens, 1, cmd);
    double theta = get_angle_arg_required(tokens, 2, cmd);
    double phi = get_angle_arg_required(tokens, 3, cmd);
    double lambda = get_angle_arg_required(tokens, 4, cmd);
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
    int j = get_arg(tokens, 1, cmd);
    int k = get_arg(tokens, 2, cmd);
    double theta = get_angle_arg_required(tokens, 3, cmd);
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
    int j = get_arg(tokens, 1, cmd);
    int k = get_arg(tokens, 2, cmd);
    double theta = get_angle_arg_required(tokens, 3, cmd);
    double phi = get_angle_arg_required(tokens, 4, cmd);
    double lambda = get_angle_arg_required(tokens, 5, cmd);
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
    int j = get_arg(tokens, 1, cmd);
    int k = get_arg(tokens, 2, cmd);
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

  if (cmd == "CCX") {
    int c1 = get_arg(tokens, 1, cmd);
    int c2 = get_arg(tokens, 2, cmd);
    int t = get_arg(tokens, 3, cmd);
    if (c1 == -1 || c2 == -1 || t == -1) return;

    try {
      state->ccx(c1, c2, t);
      std::cout << "CCX(" << c1 << ", " << c2 << ", " << t << ") applied.\n";
      state->display();
    } catch (const std::exception& e) {
      std::cerr << "Operation failed: " << e.what() << "\n";
    }
    return;
  }
  if (cmd == "SWAP") {
    int j = get_arg(tokens, 1, "SWAP");
    int k = get_arg(tokens, 2, "SWAP");
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
      int t = get_arg(tokens, i, "GROVER");
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
      int maybe_iters = get_arg(tokens, idx, "LATIN");
      if (maybe_iters == -1) return;
      if (mode == "DEMO") {
        iters = maybe_iters;
        ++idx;
      }
    }

    while (idx < tokens.size()) {
      int v = get_arg(tokens, idx, "LATIN");
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
    int N = get_arg(tokens, 1, "SHOR");
    if (N == -1) return;
    run_shor_demo(static_cast<Bitstring>(N));
    return;
  }

  // --- Measurement ---
  if (cmd == "MEASURE") {
    int j = get_arg(tokens, 1, "MEASURE"); // Qubit index
    int c = get_arg(tokens, 2, "MEASURE"); // Classical register index
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
  std::cout << "SWAP <j> <k>     : SWAP qubits j and, k maintaining amplitudes.\n";
  std::cout << "MEASURE <j> <c>  : Measure qubit j, store result in classical register c.\n";
  std::cout << "DISPLAY          : Show the current quantum state.\n";
  std::cout << "GROVER <t...>    : Run Grover's algorithm searching for one or more targets\n";
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
            
    handle_command(parse_command(input_line));
  }
  std::cout << "Exiting simulator.\n";
}
