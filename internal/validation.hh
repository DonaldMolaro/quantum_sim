#pragma once
#include <stdexcept>
#include <string>

namespace qsim {

/// Throws std::out_of_range if j is not in [0, num_qubits).
inline void check_qubit(int j, int num_qubits, const char* gate_name)
{
  if (j < 0 || j >= num_qubits) {
    throw std::out_of_range(std::string(gate_name) + ": qubit index " +
                            std::to_string(j) + " out of range [0, " +
                            std::to_string(num_qubits - 1) + "]");
  }
}

} // namespace qsim
