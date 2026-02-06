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

  // Must be initialized before applying quantum operations
  if (state == nullptr || !state->is_initialized()) {
    std::cerr << "Error: State must be initialized first. Use INIT <N> [C].\n";
    return;
  }

  // --- Gate Operations (Single Qubit) ---
  if (cmd == "H" || cmd == "X" || cmd == "S" || cmd == "T") {
    int j = get_arg(tokens, 1, cmd);
    if (j == -1) return; 

    try {
      if (cmd == "H") state->h(j); // Implements flatMap + reduceByKey
      else if (cmd == "X") state->x(j); // Implements s.map(λb, a. (b¬j , a))
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
    double theta = get_double_arg(tokens, 2, cmd);
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

  // --- Gate Operations (Two Qubits) ---
  if (cmd == "CX") {
    int j = get_arg(tokens, 1, "CX");
    int k = get_arg(tokens, 2, "CX");
    if (j == -1 || k == -1) return;

    state->cx(j, k); // Implements s.map(λb, a. (ite(bj , b ¬k, b), a))
    std::cout << "CX(" << j << ", " << k << ") applied.\n";
    state->display();
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
  std::cout << "S <j>            : S gate on qubit j (Phase).\n";
  std::cout << "T <j>            : T gate on qubit j (Phase).\n";
  std::cout << "RX <j> <theta>   : Rotation around X by angle theta (radians).\n";
  std::cout << "RY <j> <theta>   : Rotation around Y by angle theta (radians).\n";
  std::cout << "RZ <j> <theta>   : Rotation around Z by angle theta (radians).\n";
  std::cout << "CX <j> <k>       : Controlled-X (CNOT) where j is control, k is target.\n";
  std::cout << "SWAP <j> <k>     : SWAP qubits j and, k maintaining amplitudes.\n";
  std::cout << "MEASURE <j> <c>  : Measure qubit j, store result in classical register c.\n";
  std::cout << "DISPLAY          : Show the current quantum state.\n";
  std::cout << "GROVER <t...>    : Run Grover's algorithm searching for one or more targets\n";
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
