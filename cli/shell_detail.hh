#pragma once

#include "state.hh"
#include "logging.hh"
#include "cli/commands.hh"
#include "internal/format_utils.hh"
#include "internal/limits.hh"
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace shell_detail {

inline bool require_initialized_state(const std::unique_ptr<State>& state)
{
  if (!state || !state->is_initialized()) {
    std::cerr << "Error: State must be initialized first. Use INIT <N> [C].\n";
    return false;
  }
  return true;
}

inline bool parse_square_matrix(const std::vector<std::string>& tokens,
                                size_t start_index,
                                int n,
                                const std::string& cmd,
                                std::vector<double>& out)
{
  if (!qsim::limits::valid_bitstring_qubit_count(n)) {
    std::cerr << "Error: " << cmd << " requires 1 <= n <= 62.\n";
    return false;
  }
  const size_t expected = static_cast<size_t>(n) * static_cast<size_t>(n);
  if (tokens.size() != start_index + expected) {
    std::cerr << "Error: " << cmd << " requires n*n matrix entries.\n";
    return false;
  }
  out.clear();
  out.reserve(expected);
  for (size_t i = 0; i < expected; ++i) {
    const double v = cli::get_double_arg(tokens, start_index + i, cmd);
    if (std::isnan(v)) {
      out.clear();
      return false;
    }
    out.push_back(v);
  }
  return true;
}

inline std::string bits_to_string(const std::vector<int>& bits)
{
  return qsim::format::bits_to_string_msb(bits);
}

inline uint64_t bits_to_u64(const std::vector<int>& bits)
{
  uint64_t value = 0;
  for (size_t i = 0; i < bits.size() && i < 63; ++i) {
    if (bits[i]) {
      value |= (1ULL << i);
    }
  }
  return value;
}

inline int infer_qubits_from_targets(const std::vector<Bitstring>& targets)
{
  return qsim::format::infer_qubits_from_targets(targets);
}

inline const char* log_level_name(qsim_log::Level level)
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

} // namespace shell_detail
