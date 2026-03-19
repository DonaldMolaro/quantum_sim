#include "cli/shell.hh"
#include "cli/commands.hh"
#include "algorithms/qrng.hh"
#include "cli/shell_detail.hh"
#include "logging.hh"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

bool QuantumShell::handle_setup_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  if (cmd == "TUTOR") {
    if (tokens.size() < 2) {
      std::cout << "Tutor mode is " << (tutor_mode ? "ON" : "OFF") << ".\n";
      return true;
    }
    if (tokens[1] == "ON") {
      tutor_mode = true;
      std::cout << "Tutor mode enabled.\n";
      tutor_note("I will explain intent and expected outcomes for key commands.");
    } else if (tokens[1] == "OFF") {
      tutor_mode = false;
      std::cout << "Tutor mode disabled.\n";
    } else {
      std::cerr << "Error: TUTOR must be ON or OFF.\n";
    }
    return true;
  }

  if (cmd == "INIT") {
    int N = cli::get_arg(tokens, 1, "INIT");
    int C = (tokens.size() > 2) ? cli::get_arg(tokens, 2, "INIT") : 0;
    if (N > 0 && N <= 64) {
      state.reset(new State(N, C));
      std::cout << "State initialized with " << N << " qubits and " << C << " classical register(s).\n";
      tutor_note("INIT sets the Hilbert space size to 2^N and resets amplitudes to |0...0>.");
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
      tutor_note("DIRECT applies QFT as a transform; it is faster in this simulator.");
    } else if (tokens[1] == "GATE") {
      state->set_qft_mode(State::QftMode::Gate);
      std::cout << "QFT mode set to GATE.\n";
      tutor_note("GATE applies decomposed gates; use it to inspect algorithm structure.");
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
    tutor_note("QRNG uses Hadamard + measurement; each output bit is sampled from a superposition.");
    return true;
  }

  if (cmd == "SEED") {
    int seed = cli::get_arg(tokens, 1, "SEED");
    if (seed < 0) {
      std::cerr << "Error: SEED requires a non-negative integer.\n";
      return true;
    }
    std::srand(static_cast<unsigned int>(seed));
    std::string seed_s = std::to_string(seed);
    setenv("QSIM_RNG_SEED", seed_s.c_str(), 1);
    std::cout << "Random seed set to " << seed << ".\n";
    tutor_note("Using a fixed seed makes demos reproducible for teaching and grading.");
    return true;
  }

  if (cmd == "SHOTS") {
    // SHOTS <n> <circuit_file>  — run a circuit file n times, accumulate measurement histogram.
    if (tokens.size() < 3) {
      std::cerr << "Error: SHOTS requires <n> <circuit_file>.\n";
      return true;
    }
    int n = cli::get_arg(tokens, 1, "SHOTS");
    if (n <= 0) { std::cerr << "Error: SHOTS n must be > 0.\n"; return true; }
    const std::string& filename = tokens[2];

    // Read the circuit once.
    std::vector<std::string> circuit_lines;
    {
      std::ifstream in(filename);
      if (!in) { std::cerr << "Error: cannot open '" << filename << "'.\n"; return true; }
      std::string line;
      while (std::getline(in, line)) {
        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
        if (trimmed.empty() || trimmed[0] == '#') continue;
        std::transform(trimmed.begin(), trimmed.end(), trimmed.begin(), ::toupper);
        if (trimmed == "QUIT") break;
        circuit_lines.push_back(trimmed);
      }
    }
    if (circuit_lines.empty()) { std::cerr << "Error: circuit file is empty.\n"; return true; }

    // Suppress output during shots.
    std::streambuf* orig_cout = std::cout.rdbuf();
    std::streambuf* orig_cerr = std::cerr.rdbuf();
    std::ostringstream sink;
    std::cout.rdbuf(sink.rdbuf());
    std::cerr.rdbuf(sink.rdbuf());

    std::map<std::vector<int>, int> histogram;
    loading_from_file_ = true;
    for (int shot = 0; shot < n; ++shot) {
      for (const std::string& line : circuit_lines)
        handle_command(cli::parse_command(line));
      if (state) {
        const std::vector<int> cbits = state->get_cbits();
        histogram[cbits]++;
      }
    }
    loading_from_file_ = false;
    std::cout.rdbuf(orig_cout);
    std::cerr.rdbuf(orig_cerr);

    std::cout << "SHOTS " << n << " results from '" << filename << "':\n";
    // Sort by frequency descending.
    std::vector<std::pair<int, std::vector<int>>> sorted;
    for (const auto& p : histogram) sorted.push_back({p.second, p.first});
    std::sort(sorted.begin(), sorted.end(),
              [](const std::pair<int,std::vector<int>>& a, const std::pair<int,std::vector<int>>& b){ return a.first > b.first; });
    for (const auto& p : sorted) {
      std::cout << "  |";
      for (int i = static_cast<int>(p.second.size()) - 1; i >= 0; --i)
        std::cout << p.second[i];
      std::cout << "> : " << p.first << " (" << std::fixed << std::setprecision(1)
                << (100.0 * p.first / n) << "%)\n";
    }
    tutor_note("SHOTS repeats the circuit n times; measurement statistics reveal the true probability distribution.");
    return true;
  }

  if (cmd == "NOISE") {
    if (tokens.size() < 2) {
      double p = state ? state->get_noise_probability() : 0.0;
      std::cout << "Noise probability: " << p << (p > 0.0 ? " (active)" : " (off)") << "\n";
      return true;
    }
    double p = cli::get_double_arg(tokens, 1, "NOISE");
    if (std::isnan(p) || p < 0.0 || p > 1.0) {
      std::cerr << "Error: NOISE probability must be in [0.0, 1.0].\n";
      return true;
    }
    if (state) state->set_noise_probability(p);
    if (p == 0.0) {
      std::cout << "Noise disabled.\n";
    } else {
      std::cout << "Depolarizing noise enabled: p=" << p << " per primitive gate per qubit.\n";
      tutor_note("Each primitive gate now has a " + std::to_string(p * 100) + "% chance of applying a random Pauli error.");
    }
    return true;
  }

  if (cmd == "SAVE") {
    if (tokens.size() < 2) {
      std::cerr << "Error: SAVE requires a filename.\n";
      return true;
    }
    const std::string& filename = tokens[1];
    std::ofstream out(filename);
    if (!out) {
      std::cerr << "Error: cannot open '" << filename << "' for writing.\n";
      return true;
    }
    for (const std::string& line : command_history_) {
      out << line << "\n";
    }
    std::cout << "Circuit saved to '" << filename << "' (" << command_history_.size() << " commands).\n";
    return true;
  }

  if (cmd == "LOAD") {
    if (tokens.size() < 2) {
      std::cerr << "Error: LOAD requires a filename.\n";
      return true;
    }
    const std::string& filename = tokens[1];
    std::ifstream in(filename);
    if (!in) {
      std::cerr << "Error: cannot open '" << filename << "' for reading.\n";
      return true;
    }
    std::string line;
    int count = 0;
    loading_from_file_ = true;
    while (std::getline(in, line)) {
      // Skip blank lines and comments
      std::string trimmed = line;
      trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
      if (trimmed.empty() || trimmed[0] == '#') continue;
      std::transform(trimmed.begin(), trimmed.end(), trimmed.begin(), ::toupper);
      if (trimmed == "QUIT") break;
      std::cout << "LOAD> " << trimmed << "\n";
      handle_command(cli::parse_command(trimmed));
      ++count;
    }
    loading_from_file_ = false;
    std::cout << "Loaded " << count << " commands from '" << filename << "'.\n";
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
