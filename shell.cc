/*
 * This is a simple little shell that reads commands and manipulates
 * the underlying quantum state simulator. It's not bad, it's not great
 * but with it you can poke at the machine and see what it does.
 *
 * Beyond just basic quantum gates it's also able to call an implementation of
 * Grover's algorithim. Maybe I'll get ambitious and implement a few other
 * quantum algorithms.
 */
#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <sstream>
#include "state.hh"
#include "shell.hh"

extern void run_grover_search(State *s,Bitstring targer_w);

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

  // --- Algorithims ---
  if (cmd == "GROVER") {
    int target_w = get_arg(tokens, 1, "GROVER"); // Target Word
    run_grover_search(state,target_w);
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
  std::cout << "CX <j> <k>       : Controlled-X (CNOT) where j is control, k is target.\n";
  std::cout << "MEASURE <j> <c>  : Measure qubit j, store result in classical register c.\n";
  std::cout << "DISPLAY          : Show the current quantum state.\n";
  std::cout << "GROVER <t>       : Run Grover's algorithim searching for T\n";
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
