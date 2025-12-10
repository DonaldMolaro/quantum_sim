#include "state.hh"

State& State::controlled_Rr(int control_j, int target_k, int r)
{
  // Phase applied only if both control and target bits are 1
  ComplexNumber phase = std::exp(2.0 * M_PI * IMAGINARY_UNIT_I / std::pow(2.0, r)); 

  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;
            
    // Check if control bit (j) is 1 AND target bit (k) is 1
    if (((current_b >> control_j) & 1) && ((current_b >> target_k) & 1)) {
      current_a *= phase;
    }
  }
  return *this;
}

/**
 * @brief Applies the controlled Inverse Phase Rotation (R_r_dag).
 * Inverse rotation phase is negative.
 */
State& State::controlled_Rr_dag(int control_j, int target_k, int r)
{
  // Inverse phase rotation: sign of angle is negated.
  ComplexNumber phase_dag = std::exp(-2.0 * M_PI * IMAGINARY_UNIT_I / std::pow(2.0, r)); 
        
  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;
            
    if (((current_b >> control_j) & 1) && ((current_b >> target_k) & 1)) {
      current_a *= phase_dag;
    }
  }
  return *this;
}

/* @brief Applies the Quantum Fourier Transform (QFT).
 * Complexity is O(n^2) gates, where n is the number of qubits.
 * Assumes qubit indices are 0 (LSB) to n-1 (MSB).
 * @param start_qubit Index of the least significant qubit in the register.
 * @param end_qubit Index of the most significant qubit in the register.
 */
State& State::qft(int start_qubit, int end_qubit)
{
  const int N = end_qubit - start_qubit + 1;

  // Apply Hadamards and Controlled Rotations
  for (int j = 0; j < N; ++j) {
    int target_q = start_qubit + j; // Target qubit index (LSB to MSB)
        
    // 1. Apply Hadamard H_j
    h(target_q); 

    // 2. Apply controlled phase rotations C-R_r
    // Controlled by preceding qubits (k < j)
    for (int k = 1; k < N - j; ++k) {
      int control_q = target_q + k;
      int rotation_r = k + 1; 

      // C-R_r is controlled by Q_{j+k} and acts on Q_j.
      controlled_Rr(control_q, target_q, rotation_r);
    }
  }
    
  // 3. Swap qubits (Reversing the output order)
  // The outputs are naturally produced in reverse order (MSB on Q_0, LSB on Q_{n-1}) 
  // so we swap q0 with q_{n-1}, q1 with q_{n-2}, etc.
  int mid = start_qubit + N / 2;
  for (int q = start_qubit; q < mid; ++q) {
    // Calculate the index corresponding to the swap partner
    int partner_q = end_qubit - (q - start_qubit); 
    swap(q, partner_q);
  }
  return *this;
}

/**
 * @brief Applies the Inverse Quantum Fourier Transform (IQFT).
 * Implemented by reversing the QFT circuit using inverse gates.
 * @param start_qubit Index of the least significant qubit in the register.
 * @param end_qubit Index of the most significant qubit in the register.
 */
State& State::iqft(int start_qubit, int end_qubit)
{
  const int N = end_qubit - start_qubit + 1;
    
  // 1. Swap qubits first (Undoing the QFT swap)
  int mid = start_qubit + N / 2;
  for (int q = start_qubit; q < mid; ++q) {
    int partner_q = end_qubit - (q - start_qubit); 
    swap(q, partner_q);
  }

  // Apply inverse controlled phase rotations and Hadamards in reverse order
  for (int j = N - 1; j >= 0; --j) {
    int target_q = start_qubit + j; 
    
    // 2. Apply inverse controlled phase rotations (reverse order of application)
    for (int k = N - 1 - j; k >= 1; --k) { 
      int control_q = target_q + k;
      int rotation_r = k + 1; 

      controlled_Rr_dag(control_q, target_q, rotation_r);
    }
        
    // 3. Apply Hadamard (H is its own inverse)
    h(target_q); 
  }
  return *this;
}
