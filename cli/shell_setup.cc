#include "cli/shell.hh"
#include "cli/commands.hh"
#include "algorithms/qrng.hh"
#include "cli/shell_detail.hh"
#include "logging.hh"
#include <cmath>
#include <iostream>

bool QuantumShell::handle_setup_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  if (cmd == "INIT") {
    int N = cli::get_arg(tokens, 1, "INIT");
    int C = (tokens.size() > 2) ? cli::get_arg(tokens, 2, "INIT") : 0;
    if (N > 0 && N <= 64) {
      state.reset(new State(N, C));
      std::cout << "State initialized with " << N << " qubits and " << C << " classical register(s).\n";
    }
    return true;
  }

  if (cmd == "QFTMODE") {
    if (tokens.size() < 2) {
      std::cerr << "Error: QFTMODE requires DIRECT or GATE.\n";
      return true;
    }
    if (!shell_detail::require_initialized_state(state)) return true;
    if (tokens[1] == "DIRECT") {
      state->set_qft_mode(State::QftMode::Direct);
      std::cout << "QFT mode set to DIRECT.\n";
    } else if (tokens[1] == "GATE") {
      state->set_qft_mode(State::QftMode::Gate);
      std::cout << "QFT mode set to GATE.\n";
    } else {
      std::cerr << "Error: QFTMODE must be DIRECT or GATE.\n";
    }
    return true;
  }

  if (cmd == "QRNG") {
    int n = cli::get_arg(tokens, 1, "QRNG");
    int count = (tokens.size() > 2) ? cli::get_arg(tokens, 2, "QRNG") : 1;
    if (n <= 0) {
      std::cerr << "Error: QRNG requires n > 0.\n";
      return true;
    }
    if (count <= 0) {
      std::cerr << "Error: QRNG count must be > 0.\n";
      return true;
    }
    if (n > 63) {
      std::cout << "Note: QRNG displays integer values using only the lowest 63 bits.\n";
    }
    for (int i = 0; i < count; ++i) {
      std::vector<int> bits = qrng_bits(n);
      std::cout << "QRNG[" << i << "] bits=" << shell_detail::bits_to_string(bits)
                << " value=" << shell_detail::bits_to_u64(bits) << "\n";
    }
    return true;
  }

  if (cmd == "VERBOSE" || cmd == "LOGLEVEL") {
    if (tokens.size() < 2) {
      std::cout << "Verbosity: " << shell_detail::log_level_name(qsim_log::get_level()) << "\n";
      return true;
    }
    qsim_log::Level level;
    if (!qsim_log::parse_level(tokens[1], level)) {
      std::cerr << "Error: VERBOSE expects QUIET, NORMAL, VERBOSE, or 0/1/2.\n";
      return true;
    }
    qsim_log::set_level(level);
    std::cout << "Verbosity set to " << shell_detail::log_level_name(level) << ".\n";
    return true;
  }
  return false;
}
