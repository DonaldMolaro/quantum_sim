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

void QuantumShell::tutor_note(const std::string& msg) const
{
  if (!tutor_mode) return;
  std::cout << "[TUTOR] " << msg << "\n";
}

void QuantumShell::handle_command(const std::vector<std::string>& tokens)
{
  if (tokens.empty()) return;
  const std::string cmd = tokens[0];
  // Record non-meta commands to history (skip SAVE/LOAD themselves).
  if (!loading_from_file_ && cmd != "SAVE" && cmd != "LOAD" && cmd != "HELP" && cmd != "QUIT") {
    // Reconstruct the raw token string for history.
    std::string raw;
    for (size_t i = 0; i < tokens.size(); ++i) {
      if (i > 0) raw += ' ';
      raw += tokens[i];
    }
    command_history_.push_back(raw);
  }
  if (handle_setup_commands(tokens, cmd)) return;
  if (handle_algorithm_commands(tokens, cmd)) return;
  if (!shell_detail::require_initialized_state(state)) return;
  if (handle_single_qubit_commands(tokens, cmd)) return;
  if (handle_parametric_gate_commands(tokens, cmd)) return;
  if (handle_multi_qubit_commands(tokens, cmd)) return;
  if (handle_measurement_and_display_commands(tokens, cmd)) return;
  std::cout << "Unknown command: " << cmd << ". Use HELP TOPICS or HELP ALL.\n";
}

void QuantumShell::print_help_summary()
{
  std::cout << "\n--- Quantum Simulator ---\n";
  std::cout << "Start here:\n";
  std::cout << "  INIT 2 2        create 2 qubits and 2 classical bits\n";
  std::cout << "  H 0             put q0 into superposition\n";
  std::cout << "  CX 0 1          entangle q0 and q1\n";
  std::cout << "  MEASURE 0 0     measure q0 into c0\n";
  std::cout << "  DISPLAY         show amplitudes and classical bits\n";
  std::cout << "\nCommon tasks:\n";
  std::cout << "  Gates      : H X Y Z S T RX RY RZ CX CZ SWAP CCX\n";
  std::cout << "  Algorithms : GROVER BV DEUTSCH_JOZSA SIMON SHOR QPE QCOUNT VQA VQE ANNEAL\n";
  std::cout << "  Utilities  : CHECK EXPECT BLOCH ENTROPY SHOTS LOAD SAVE SEED VERBOSE TUTOR\n";
  std::cout << "\nHelp topics:\n";
  std::cout << "  HELP TOPICS      list help sections\n";
  std::cout << "  HELP GATES       gate reference\n";
  std::cout << "  HELP ALGORITHMS  algorithm commands\n";
  std::cout << "  HELP UTILITY     measurement, checks, files, randomness\n";
  std::cout << "  HELP ALL         full command reference\n";
  std::cout << "\nTip: qubit indexing is little-endian, so q0 is the least-significant bit.\n";
  std::cout << "-------------------------\n";
}

void QuantumShell::print_help_topics()
{
  std::cout << "\nHelp topics:\n";
  std::cout << "  HELP             short getting-started guide\n";
  std::cout << "  HELP TOPICS      this topic list\n";
  std::cout << "  HELP GATES       gate syntax and angle notes\n";
  std::cout << "  HELP ALGORITHMS  demos and higher-level commands\n";
  std::cout << "  HELP UTILITY     measurement, display, files, noise, RNG\n";
  std::cout << "  HELP ALL         full reference\n";
}

void QuantumShell::print_help_gates()
{
  std::cout << "\n[Gates]\n";
  std::cout << "Single-qubit:\n";
  std::cout << "  H <j>   X <j>   Y <j>   Z <j>   S <j>   SDG <j>   T <j>   TDG <j>\n";
  std::cout << "  P <j> <phi>\n";
  std::cout << "  RX <j> <theta>   RY <j> <theta>   RZ <j> <theta>\n";
  std::cout << "  RU <j> <theta> <phi> <lambda>\n";
  std::cout << "Multi-qubit:\n";
  std::cout << "  CX/CNOT <c> <t>   CZ <c> <t>   CY <c> <t>   CH <c> <t>\n";
  std::cout << "  CRX/CRY/CRZ <c> <t> <theta>   CP <c> <t> <phi>   CU <c> <t> <th> <ph> <la>\n";
  std::cout << "  CCX/TOFFOLI <c1> <c2> <t>   MCX <c1> [c2 ...] <t>\n";
  std::cout << "  SWAP <j> <k>   ISWAP <j> <k>   CSWAP/FREDKIN <c> <a> <b>\n";
  std::cout << "  XX <j> <k> <theta>   YY <j> <k> <theta>   ZZ <j> <k> <theta>\n";
  std::cout << "Angles:\n";
  std::cout << "  Radians by default. Constants: PI PI/2 PI/4 TAU. Degrees: append DEG.\n";
  std::cout << "  Examples: RX 0 PI/2   RY 1 45 DEG   RZ 0 -PI/4\n";
}

void QuantumShell::print_help_algorithms()
{
  std::cout << "\n[Algorithms]\n";
  std::cout << "Search and oracles:\n";
  std::cout << "  GROVER <targets...>\n";
  std::cout << "  GROVER AUTO <n> <count_iters> <targets...>\n";
  std::cout << "  DEUTSCH_JOZSA <n> <CONST0|CONST1|BALANCED_XOR0|BALANCED_PARITY>\n";
  std::cout << "  BV <n> <secret> [bias]\n";
  std::cout << "  SIMON DEMO\n";
  std::cout << "  SIMON <n> <secret> [shots]\n";
  std::cout << "Factoring and phase estimation:\n";
  std::cout << "  SHOR <N>\n";
  std::cout << "  QPE DEMO\n";
  std::cout << "  QPE <m> <phase>\n";
  std::cout << "  QCOUNT DEMO\n";
  std::cout << "  QCOUNT RUN <n> <iters> <t1> [t2 ...]\n";
  std::cout << "Optimization and annealing:\n";
  std::cout << "  QUBO DEMO\n";
  std::cout << "  QUBO EXACT <n> <n*n matrix entries>\n";
  std::cout << "  QUBO GROVER <n> <threshold> <iterations> <n*n matrix entries>\n";
  std::cout << "  QAOA DEMO   |   VQA DEMO\n";
  std::cout << "  QAOA QUBO <n> <p> <shots> <iters> <step> <entries...>\n";
  std::cout << "  VQA QAOA <n> <p> <shots> <iters> <step> <entries...>\n";
  std::cout << "  VQE DEMO\n";
  std::cout << "  VQE RUN <n> <layers> <iters> <step> <shots> <terms> ...\n";
  std::cout << "  ANNEAL DEMO\n";
  std::cout << "  ANNEAL QUBO <SA|SQA> <n> <steps> <sweeps> <beta_start> <beta_end> <replicas> <entries...>\n";
  std::cout << "Other demos:\n";
  std::cout << "  TSP DEMO   |   TSP EXACT <n> <penalty> <distance entries...>\n";
  std::cout << "  QEC DEMO   |   QEC RUN <logical_bit> <error_qubit>\n";
  std::cout << "  LATIN [iters]   LATIN DEMO [iters] [r0 r1 r2]   LATIN COUNT [r0 r1 r2]\n";
}

void QuantumShell::print_help_utility()
{
  std::cout << "\n[Utility]\n";
  std::cout << "State and measurement:\n";
  std::cout << "  INIT <n> [c]   MEASURE <q> <c>   RESET <q>   IF <c> <gate...>\n";
  std::cout << "  DISPLAY   DISPLAY ALL   CHECK <mode> ...\n";
  std::cout << "Inspection:\n";
  std::cout << "  EXPECT <P> <q> [P q ...]   FIDELITY <idx>   BLOCH <q>   ENTROPY <j> [k]\n";
  std::cout << "Files and repeated runs:\n";
  std::cout << "  LOAD <file>   SAVE <file>   SHOTS <n> <file>\n";
  std::cout << "Runtime controls:\n";
  std::cout << "  VERBOSE <QUIET|NORMAL|VERBOSE|0|1|2>\n";
  std::cout << "  TUTOR <ON|OFF>   SEED <n>   NOISE <p>   QFTMODE <DIRECT|GATE>\n";
  std::cout << "Randomness:\n";
  std::cout << "  QRNG <n> [count]\n";
}

void QuantumShell::print_help_all()
{
  std::cout << "\n--- Quantum Simulator Commands ---\n";
  std::cout << "\n[Initialization]\n";
  std::cout << "INIT <N> [C]     : Initialize state with N qubits and C classical bits.\n";
  std::cout << "\n[Simple Gates]\n";
  std::cout << "H <j>            : Hadamard gate on qubit j (Superposition).\n";
  std::cout << "X <j>            : NOT gate on qubit j.\n";
  std::cout << "Y <j>            : Y gate on qubit j.\n";
  std::cout << "Z <j>            : Z gate on qubit j (Phase flip).\n";
  std::cout << "S <j>            : S gate on qubit j (Phase = +i on |1>).\n";
  std::cout << "SDG <j>          : S† gate on qubit j (Phase = -i on |1>, inverse of S).\n";
  std::cout << "T <j>            : T gate on qubit j (Phase = e^{i*pi/4} on |1>).\n";
  std::cout << "TDG <j>          : T† gate on qubit j (inverse of T).\n";
  std::cout << "P <j> <phi>      : Phase gate P(phi): applies e^{i*phi} to |1>, identity to |0>.\n";
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
  std::cout << "CP <j> <k> <phi> : Controlled-P: applies e^{i*phi} when both j and k are |1>.\n";
  std::cout << "CCX <c1> <c2> <t>: Toffoli (double-controlled X).\n";
  std::cout << "TOFFOLI <c1> <c2> <t>: Alias for CCX.\n";
  std::cout << "CSWAP <j> <k> <l>: Fredkin (controlled SWAP): swaps k,l when j=|1>.\n";
  std::cout << "FREDKIN <j> <k> <l>: Alias for CSWAP.\n";
  std::cout << "MCX <c1> [c2 ...] <t>: Multi-controlled X: flips t when ALL controls are |1>.\n";
  std::cout << "SWAP <j> <k>     : SWAP qubits j and k.\n";
  std::cout << "ISWAP <j> <k>    : iSWAP gate (native superconducting gate).\n";
  std::cout << "XX <j> <k> <t>   : Ising XX(theta) gate: exp(-i t/2 X⊗X).\n";
  std::cout << "YY <j> <k> <t>   : Ising YY(theta) gate: exp(-i t/2 Y⊗Y).\n";
  std::cout << "ZZ <j> <k> <t>   : Ising ZZ(theta) gate: exp(-i t/2 Z⊗Z).\n";
  std::cout << "\n[Algorithms]\n";
  std::cout << "GROVER <t...>    : Run Grover search (self-initialized; no INIT required)\n";
  std::cout << "GROVER AUTO <n> <count_iters> <t...> : Auto-tune Grover iterations using quantum counting\n";
  std::cout << "DEUTSCH_JOZSA <n> <oracle> : Run Deutsch-Jozsa demo (CONST0, CONST1, BALANCED_XOR0, BALANCED_PARITY)\n";
  std::cout << "BV <n> <secret> [bias]: Run Bernstein-Vazirani demo (bias defaults to 0)\n";
  std::cout << "SIMON DEMO       : Run built-in Simon's algorithm demo\n";
  std::cout << "SIMON <n> <secret> [shots]: Recover hidden xor-mask secret string\n";
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
  std::cout << "QEC DEMO         : Run built-in 3-qubit bit-flip quantum error correction demo\n";
  std::cout << "QEC RUN <log> <err>: Run QEC with logical qubit log (0|1) and error on qubit err (-1=none)\n";
  std::cout << "QPE DEMO         : Run built-in Quantum Phase Estimation demo\n";
  std::cout << "QPE <m> <phase>  : Run QPE with m precision qubits for P(phase) gate (radians)\n";
  std::cout << "QCOUNT DEMO      : Run built-in quantum counting demo\n";
  std::cout << "QCOUNT RUN <n> <iters> <t1> [t2 ...] : Estimate marked-state count with explicit fit iterations\n";
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
  std::cout << "RESET <j>        : Unconditionally reset qubit j to |0> (ancilla reuse).\n";
  std::cout << "IF <c> <gate...> : Apply gate only if classical bit c is 1.\n";
  std::cout << "DISPLAY          : Show the current quantum state.\n";
  std::cout << "CHECK <mode> ... : Run quick state sanity checks (NORMALIZED, BASIS, TARGETS, BELL)\n";
  std::cout << "EXPECT <P> <q> [P q ...]: Expectation value of Pauli product (I/X/Y/Z per qubit).\n";
  std::cout << "FIDELITY <idx>   : |<idx|psi>|^2 — fidelity with computational basis state.\n";
  std::cout << "BLOCH <j>        : Bloch sphere (x,y,z) for single-qubit reduced state of qubit j.\n";
  std::cout << "ENTROPY <j> [k]  : Von Neumann entropy (bits) of qubit subsystem j..k.\n";
  std::cout << "SWAP_TEST <anc> <a_start> <b_start> <n>: Estimate |<A|B>|^2 via CSWAP protocol.\n";
  std::cout << "HELP             : Show the short help message.\n";
  std::cout << "HELP ALL         : Show the full command reference.\n";
  std::cout << "QUIT             : Exit the simulator.\n";
  std::cout << "SHOTS <n> <file> : Run circuit file n times and print measurement histogram.\n";
  std::cout << "NOISE <p>        : Set per-gate depolarizing noise probability p in [0,1]. 0 = off.\n";
  std::cout << "LOAD <file>      : Load and execute circuit commands from a file.\n";
  std::cout << "SAVE <file>      : Save the current session's commands to a file.\n";
  std::cout << "VERBOSE <level>  : Set verbosity (QUIET, NORMAL, VERBOSE or 0/1/2)\n";
  std::cout << "TUTOR <ON|OFF>   : Enable or disable guided teaching narration\n";
  std::cout << "SEED <n>         : Set deterministic RNG seed for reproducible demos/runs\n";
  std::cout << "QFTMODE <mode>   : Set QFT implementation mode (DIRECT matrix transform or GATE decomposition).\n";
  std::cout << "                   DIRECT is the default and is generally faster for this simulator.\n";
  std::cout << "QRNG <n> [count] : Quantum random numbers from n qubits (count default 1)\n";
  std::cout << "Debug tips       : Qubit indexing is little-endian (q0 is LSB); use CHECK for fast validation.\n";
  std::cout << "----------------------------------\n";
}

void QuantumShell::print_help(const std::vector<std::string>& tokens)
{
  if (tokens.size() <= 1) {
    print_help_summary();
    return;
  }

  const std::string topic = tokens[1];
  if (topic == "TOPICS") {
    print_help_topics();
  } else if (topic == "GATES") {
    print_help_gates();
  } else if (topic == "ALGORITHMS") {
    print_help_algorithms();
  } else if (topic == "UTILITY") {
    print_help_utility();
  } else if (topic == "ALL") {
    print_help_all();
  } else {
    std::cout << "Unknown help topic: " << topic << ". Use HELP TOPICS.\n";
  }
}

void QuantumShell::run()
{
  State::set_default_log_stream(&std::cout);
  qsim_log::set_stream(&std::cout);
  print_help_summary();
  std::string input_line;
  while (true) {
    std::cout << "QSIM> ";
    if (!std::getline(std::cin, input_line) || input_line == "QUIT" || input_line == "quit") {
      break;
    }
    if (input_line.empty()) continue;

    std::transform(input_line.begin(), input_line.end(), input_line.begin(), ::toupper);

    std::vector<std::string> tokens = cli::parse_command(input_line);
    if (!tokens.empty() && tokens[0] == "HELP") {
      print_help(tokens);
      continue;
    }

    handle_command(tokens);
  }
  std::cout << "Exiting simulator.\n";
}
