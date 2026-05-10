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
#include "internal/limits.hh"
#include "logging.hh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

void QuantumShell::tutor_note(const std::string& msg) const
{
  if (!tutor_mode) return;
  std::cout << "[TUTOR] " << msg << "\n";
}

QuantumShell::TutorSnapshot QuantumShell::capture_tutor_snapshot() const
{
  TutorSnapshot snapshot;
  if (!tutor_mode || !state) {
    return snapshot;
  }

  snapshot.num_qubits = state->get_num_qubits();
  snapshot.cbits = state->get_cbits();
  for (const auto& entry : state->get_state()) {
    snapshot.amplitudes[entry.first] = entry.second;
  }
  return snapshot;
}

void QuantumShell::tutor_state_delta(const TutorSnapshot& before, int max_entries) const
{
  if (!tutor_mode || !state) return;

  TutorSnapshot after = capture_tutor_snapshot();
  struct Change {
    Bitstring basis = 0;
    ComplexNumber before_amp = 0.0;
    ComplexNumber after_amp = 0.0;
    double score = 0.0;
  };

  std::vector<Change> changes;
  std::map<Bitstring, bool> seen;
  for (const auto& entry : before.amplitudes) seen[entry.first] = true;
  for (const auto& entry : after.amplitudes) seen[entry.first] = true;

  for (const auto& entry : seen) {
    const Bitstring basis = entry.first;
    const std::map<Bitstring, ComplexNumber>::const_iterator before_it = before.amplitudes.find(basis);
    const std::map<Bitstring, ComplexNumber>::const_iterator after_it = after.amplitudes.find(basis);
    const ComplexNumber before_amp = (before_it != before.amplitudes.end()) ? before_it->second : ComplexNumber(0.0, 0.0);
    const ComplexNumber after_amp = (after_it != after.amplitudes.end()) ? after_it->second : ComplexNumber(0.0, 0.0);
    const double amp_delta = std::abs(after_amp - before_amp);
    const double prob_delta = std::abs(std::norm(after_amp) - std::norm(before_amp));
    const double score = std::max(amp_delta, prob_delta);
    if (score > qsim::limits::AMPLITUDE_EPSILON) {
      changes.push_back({basis, before_amp, after_amp, score});
    }
  }

  std::sort(changes.begin(), changes.end(), [](const Change& lhs, const Change& rhs) {
    if (std::abs(lhs.score - rhs.score) > qsim::limits::AMPLITUDE_EPSILON) {
      return lhs.score > rhs.score;
    }
    return lhs.basis < rhs.basis;
  });

  int cbit_changes = 0;
  const size_t max_cbits = std::max(before.cbits.size(), after.cbits.size());
  for (size_t i = 0; i < max_cbits; ++i) {
    const int before_bit = (i < before.cbits.size()) ? before.cbits[i] : -1;
    const int after_bit = (i < after.cbits.size()) ? after.cbits[i] : -1;
    if (before_bit != after_bit) {
      ++cbit_changes;
    }
  }

  if (changes.empty() && cbit_changes == 0) {
    tutor_note("No visible amplitude or classical-bit change.");
    return;
  }

  std::cout << "[TUTOR] What changed:\n";
  const int limit = std::max(0, max_entries);
  const int printed = std::min(limit, static_cast<int>(changes.size()));
  for (int i = 0; i < printed; ++i) {
    const Change& change = changes[static_cast<size_t>(i)];
    std::ostringstream line;
    line << "[TUTOR]   |" << bitstring_to_string(change.basis, after.num_qubits) << ">: amp "
         << complex_to_string(change.before_amp, 4) << " -> " << complex_to_string(change.after_amp, 4)
         << ", prob ";
    line.setf(std::ios::fixed);
    line.precision(3);
    line << std::norm(change.before_amp) << " -> " << std::norm(change.after_amp);
    std::cout << line.str() << "\n";
  }
  if (static_cast<int>(changes.size()) > printed) {
    std::cout << "[TUTOR]   ... " << (changes.size() - static_cast<size_t>(printed))
              << " more basis-state changes\n";
  }

  for (size_t i = 0; i < max_cbits; ++i) {
    const int before_bit = (i < before.cbits.size()) ? before.cbits[i] : -1;
    const int after_bit = (i < after.cbits.size()) ? after.cbits[i] : -1;
    if (before_bit != after_bit) {
      if (before_bit >= 0 && after_bit >= 0) {
        std::cout << "[TUTOR]   c[" << i << "]: " << before_bit << " -> " << after_bit << "\n";
      } else if (before_bit < 0) {
        std::cout << "[TUTOR]   c[" << i << "] initialized to " << after_bit << "\n";
      } else {
        std::cout << "[TUTOR]   c[" << i << "] cleared\n";
      }
    }
  }
}

void QuantumShell::print_circuit_history() const
{
  int num_qubits = state ? state->get_num_qubits() : 0;
  int num_cbits = state ? static_cast<int>(state->get_cbits().size()) : 0;
  for (const std::string& line : command_history_) {
    const std::vector<std::string> tokens = cli::parse_command(line);
    if (tokens.size() >= 3 && cli::token_is(tokens[0], "INIT")) {
      try {
        num_qubits = std::stoi(tokens[1]);
      } catch (const std::exception&) {
      }
      try {
        num_cbits = std::stoi(tokens[2]);
      } catch (const std::exception&) {
      }
    }
  }

  if (num_qubits <= 0) {
    std::cout << "No circuit to display. Initialize a state and run some commands first.\n";
    return;
  }

  struct Column {
    std::vector<std::string> qcells;
    std::vector<std::string> ccells;
  };

  const int cell_width = 7;
  const std::string wire = std::string(cell_width, '-');
  const std::string blank = std::string(cell_width, ' ');

  auto center_cell = [&](const std::string& text) {
    if (static_cast<int>(text.size()) >= cell_width) {
      return text.substr(0, static_cast<size_t>(cell_width));
    }
    const int left = (cell_width - static_cast<int>(text.size())) / 2;
    const int right = cell_width - static_cast<int>(text.size()) - left;
    return std::string(static_cast<size_t>(left), '-') + text + std::string(static_cast<size_t>(right), '-');
  };

  auto empty_column = [&]() {
    Column col;
    col.qcells.assign(static_cast<size_t>(num_qubits), wire);
    col.ccells.assign(static_cast<size_t>(num_cbits), blank);
    return col;
  };

  auto mark_vertical = [&](Column& col, int a, int b) {
    if (a > b) std::swap(a, b);
    for (int q = a + 1; q < b; ++q) {
      col.qcells[static_cast<size_t>(q)] = std::string(cell_width / 2, '-') + "|" + std::string(cell_width - (cell_width / 2) - 1, '-');
    }
  };

  auto set_qubit_label = [&](Column& col, int q, const std::string& label) {
    if (q < 0 || q >= num_qubits) return;
    col.qcells[static_cast<size_t>(q)] = center_cell(label);
  };

  auto set_cbit_label = [&](Column& col, int c, const std::string& label) {
    if (c < 0 || c >= num_cbits) return;
    col.ccells[static_cast<size_t>(c)] = center_cell(label);
  };

  std::vector<Column> columns;
  for (const std::string& line : command_history_) {
    const std::vector<std::string> tokens = cli::parse_command(line);
    if (tokens.empty()) continue;
    const std::string cmd = cli::upper_copy(tokens[0]);
    if (cmd == "INIT") continue;

    Column col = empty_column();
    bool recognized = true;
    if ((cmd == "H" || cmd == "X" || cmd == "Y" || cmd == "Z" || cmd == "S" || cmd == "T" ||
         cmd == "SDG" || cmd == "TDG" || cmd == "P" || cmd == "RX" || cmd == "RY" || cmd == "RZ" ||
         cmd == "RU" || cmd == "RESET") && tokens.size() >= 2) {
      try {
        set_qubit_label(col, std::stoi(tokens[1]), "[" + cmd.substr(0, std::min<size_t>(3, cmd.size())) + "]");
      } catch (const std::exception&) {
        recognized = false;
      }
    } else if ((cmd == "CX" || cmd == "CNOT" || cmd == "CZ" || cmd == "CY" || cmd == "CH" ||
                cmd == "CRX" || cmd == "CRY" || cmd == "CRZ" || cmd == "CP" || cmd == "CU") && tokens.size() >= 3) {
      try {
        const int c = std::stoi(tokens[1]);
        const int t = std::stoi(tokens[2]);
        set_qubit_label(col, c, "-o-");
        std::string target = "[X]";
        if (cmd == "CZ") target = "[Z]";
        else if (cmd == "CY") target = "[Y]";
        else if (cmd == "CH") target = "[H]";
        else if (cmd == "CRX") target = "[RX]";
        else if (cmd == "CRY") target = "[RY]";
        else if (cmd == "CRZ") target = "[RZ]";
        else if (cmd == "CP") target = "[P]";
        else if (cmd == "CU") target = "[U]";
        set_qubit_label(col, t, target);
        mark_vertical(col, c, t);
      } catch (const std::exception&) {
        recognized = false;
      }
    } else if ((cmd == "SWAP" || cmd == "ISWAP" || cmd == "XX" || cmd == "YY" || cmd == "ZZ") && tokens.size() >= 3) {
      try {
        const int a = std::stoi(tokens[1]);
        const int b = std::stoi(tokens[2]);
        const std::string label = (cmd == "SWAP") ? "[S]" : (cmd == "ISWAP" ? "[iS]" : "[" + cmd + "]");
        set_qubit_label(col, a, label);
        set_qubit_label(col, b, label);
        mark_vertical(col, a, b);
      } catch (const std::exception&) {
        recognized = false;
      }
    } else if ((cmd == "CCX" || cmd == "TOFFOLI") && tokens.size() >= 4) {
      try {
        const int c1 = std::stoi(tokens[1]);
        const int c2 = std::stoi(tokens[2]);
        const int t = std::stoi(tokens[3]);
        set_qubit_label(col, c1, "-o-");
        set_qubit_label(col, c2, "-o-");
        set_qubit_label(col, t, "[X]");
        mark_vertical(col, std::min(c1, c2), std::max(c1, c2));
        mark_vertical(col, std::min(c2, t), std::max(c2, t));
      } catch (const std::exception&) {
        recognized = false;
      }
    } else if (cmd == "MEASURE" && tokens.size() >= 3) {
      try {
        const int q = std::stoi(tokens[1]);
        const int c = std::stoi(tokens[2]);
        set_qubit_label(col, q, "[M]");
        set_cbit_label(col, c, "[c" + std::to_string(c) + "]");
      } catch (const std::exception&) {
        recognized = false;
      }
    } else if (cmd == "IF" && tokens.size() >= 3) {
      try {
        const int c = std::stoi(tokens[1]);
        set_cbit_label(col, c, "[?]");
      } catch (const std::exception&) {
        recognized = false;
      }
    } else {
      recognized = false;
    }

    if (!recognized) {
      set_qubit_label(col, 0, "[...]");
    }
    columns.push_back(col);
  }

  if (columns.empty()) {
    std::cout << "Circuit history is empty.\n";
    return;
  }

  std::cout << "Circuit view (" << num_qubits << " qubits";
  if (num_cbits > 0) {
    std::cout << ", " << num_cbits << " classical bits";
  }
  std::cout << "):\n";
  for (int q = num_qubits - 1; q >= 0; --q) {
    std::cout << "q" << q << ": ";
    for (const Column& col : columns) {
      std::cout << col.qcells[static_cast<size_t>(q)] << " ";
    }
    std::cout << "\n";
  }
  for (int c = num_cbits - 1; c >= 0; --c) {
    std::cout << "c" << c << ": ";
    for (const Column& col : columns) {
      std::cout << col.ccells[static_cast<size_t>(c)] << " ";
    }
    std::cout << "\n";
  }
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
  std::cout << "  Utilities  : CHECK EXPECT BLOCH ENTROPY SHOTS CIRCUIT LOAD SAVE SEED VERBOSE TUTOR\n";
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
