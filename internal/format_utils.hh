#pragma once

#include "state.hh"
#include <algorithm>
#include <string>
#include <vector>

namespace qsim {
namespace format {

inline std::string bits_to_string_msb(const std::vector<int>& bits)
{
  if (bits.empty()) return "0";
  std::string out;
  out.reserve(bits.size());
  for (int i = static_cast<int>(bits.size()) - 1; i >= 0; --i) {
    out.push_back(bits[i] ? '1' : '0');
  }
  return out;
}

inline int infer_qubits_from_targets(const std::vector<Bitstring>& targets)
{
  Bitstring max_target = 0ULL;
  for (size_t i = 0; i < targets.size(); ++i) {
    if (targets[i] > max_target) {
      max_target = targets[i];
    }
  }
  if (max_target == 0ULL) {
    return 1;
  }
  int n = 0;
  while (max_target != 0ULL) {
    ++n;
    max_target >>= 1;
  }
  return std::max(1, n);
}

} // namespace format
} // namespace qsim
