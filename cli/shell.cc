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
#include "cli/help_catalog.hh"
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
  const std::string cmd = cli::upper_copy(tokens[0]);
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
  cli::help_catalog::print_section(std::cout, cli::help_catalog::Section::Gates);
}

void QuantumShell::print_help_algorithms()
{
  cli::help_catalog::print_section(std::cout, cli::help_catalog::Section::Algorithms);
}

void QuantumShell::print_help_utility()
{
  cli::help_catalog::print_section(std::cout, cli::help_catalog::Section::Utility);
}

void QuantumShell::print_help_all()
{
  cli::help_catalog::print_all(std::cout);
}

void QuantumShell::print_help(const std::vector<std::string>& tokens)
{
  if (tokens.size() <= 1) {
    print_help_summary();
    return;
  }

  const std::string topic = cli::upper_copy(tokens[1]);
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
    if (!std::getline(std::cin, input_line)) {
      break;
    }
    if (input_line.empty()) continue;

    std::vector<std::string> tokens = cli::parse_command(input_line);
    if (!tokens.empty() && cli::token_is(tokens[0], "QUIT")) {
      break;
    }
    if (!tokens.empty() && cli::token_is(tokens[0], "HELP")) {
      print_help(tokens);
      continue;
    }

    handle_command(tokens);
  }
  std::cout << "Exiting simulator.\n";
}
