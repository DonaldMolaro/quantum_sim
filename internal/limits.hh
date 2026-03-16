#pragma once

namespace qsim {
namespace limits {

// Maximum qubits where 1ULL << n is well-defined and fits in Bitstring.
static const int kMaxBitstringQubits = 62;

// Maximum qubits supported by the Shor algorithm demo.
static const int kMaxShorQubits = 24;

// Threshold below which an amplitude is treated as zero.
constexpr double AMPLITUDE_EPSILON = 1e-9;

// Mathematical constant pi.
constexpr double PI = 3.14159265358979323846;

inline bool valid_bitstring_qubit_count(int n)
{
  return n > 0 && n <= kMaxBitstringQubits;
}

} // namespace limits
} // namespace qsim
