/*
 * This is a simple little shell that reads commands and manipulates
 * the underlying quantum state simulator. It's not bad, it's not great
 * but with it you can poke at the machine and see what it does.
 *
 * Beyond just basic quantum gates it's also able to call an implementation of
 * Grover's algorithm. Maybe I'll get ambitious and implement a few other
 * quantum algorithms.
 */
#include "cli/shell.hh"
#include "cli/commands.hh"
#include "cli/shell_detail.hh"
#include "logging.hh"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

void QuantumShell::handle_command(const std::vector<std::string>& tokens)
{
  if (tokens.empty()) return;
  const std::string cmd = tokens[0];
  if (handle_setup_commands(tokens, cmd)) return;
  if (handle_algorithm_commands(tokens, cmd)) return;
  if (!shell_detail::require_initialized_state(state)) return;
  if (handle_single_qubit_commands(tokens, cmd)) return;
  if (handle_parametric_gate_commands(tokens, cmd)) return;
  if (handle_multi_qubit_commands(tokens, cmd)) return;
  if (handle_measurement_and_display_commands(tokens, cmd)) return;
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
  std::cout << "QAOA DEMO        : Run built-in QAOA-on-QUBO demo\n";
  std::cout << "QAOA QUBO <n> <p> <shots> <iters> <step> <n*n matrix entries> : QAOA optimize QUBO\n";
  std::cout << "VQA DEMO         : Run built-in QAOA-on-QUBO demo\n";
  std::cout << "VQA QAOA <n> <p> <shots> <iters> <step> <n*n matrix entries> : QAOA optimize QUBO\n";
  std::cout << "VQE DEMO         : Run built-in VQE demo (H=-Z0)\n";
  std::cout << "VQE RUN <n> <layers> <iters> <step> <shots> <terms> <coeff op_count [op qubit]...>\n";
  std::cout << "ANNEAL DEMO      : Run built-in simulated quantum annealing demo\n";
  std::cout << "ANNEAL QUBO <SA|SQA> <n> <steps> <sweeps> <beta_start> <beta_end> <replicas> <n*n matrix entries>\n";
  std::cout << "TSP DEMO         : Run built-in 4-city TSP-to-QUBO exact demo\n";
  std::cout << "TSP EXACT <n> <penalty> <n*n distance entries> : Solve fixed-start TSP with exact QUBO search\n";
  std::cout << "QCOUNT DEMO      : Run built-in quantum counting demo\n";
  std::cout << "QCOUNT <n> <t1> [t2 ...] : Estimate number of marked states for targets in n-qubit space\n";
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
  std::cout << "QFTMODE <mode>   : Set QFT implementation mode (DIRECT matrix transform or GATE decomposition).\n";
  std::cout << "                   DIRECT is the default and is generally faster for this simulator.\n";
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

    std::transform(input_line.begin(), input_line.end(), input_line.begin(), ::toupper);

    if (input_line == "HELP") {
      print_help();
      continue;
    }

    handle_command(cli::parse_command(input_line));
  }
  std::cout << "Exiting simulator.\n";
}
