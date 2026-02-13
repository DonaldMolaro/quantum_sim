#pragma once

namespace qsim {
namespace limits {

// Maximum qubits where 1ULL << n is well-defined and fits in Bitstring.
static const int kMaxBitstringQubits = 62;

inline bool valid_bitstring_qubit_count(int n)
{
  return n > 0 && n <= kMaxBitstringQubits;
}

} // namespace limits
} // namespace qsim
