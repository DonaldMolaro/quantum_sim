/*
 * This is the implemenation of the functions that Grover's search needs.
 * If there is something "interesting" here it's the fact that these methods
 * are unaware of the underlying non-linear algebra implementation.
 */

#include "state.hh"
#include <cstdlib>
#include <ostream>

static bool grover_verbose()
{
  static int initialized = 0;
  static bool enabled = false;
  if (!initialized) {
    const char* env = std::getenv("QSIM_GROVER_VERBOSE");
    enabled = (env && env[0] != '\0' && env[0] != '0');
    initialized = 1;
  }
  return enabled;
}

/**
 * @brief Applies the U_{0^perp} operator. 
 * This negates the amplitude of every state except |00...0>.
 * This is implemented as a simple s.map operation in the functional view.
 */
State& State::apply_U0_perp() {
        
  // Iterate over the sparse representation of the quantum state
  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;

    // If the bitstring is anything other than 0 (i.e., |00...0>), negate the amplitude.
    if (current_b != 0ULL) { 
      current_a *= -1.0; 
    }
  }
  // Note: Filtering zero amplitudes after this step is often useful, 
  // but unnecessary if the input state is already normalized and sparse.
  return *this;
}

/**
 * @brief Applies the Grover Diffusion Operator (U_psi_perp or U_s).
 * U_psi_perp = H^n * U_{0^perp} * H^n.
 * 
 * Assumes the existence of:
 * 1. h(j): Hadamard gate application to qubit j.
 */
State& State::grover_diffusion_Us()
{
  const int N = num_qubits_; 

  // 1. Apply H^n (Hadamard to all qubits)
  for (int j = 0; j < N; ++j) {
    // Note: Since qubit indexing uses little endian (right-to-left starting at 0),
    // we iterate from 0 to N-1.
    h(j); 
  }

  // 2. Apply U_{0^perp} (Reflection about |00...0>)
  apply_U0_perp(); 

  // 3. Apply H^n again
  for (int j = 0; j < N; ++j) {
    h(j); 
  }
        
  if (grover_verbose() && log_stream_) {
    (*log_stream_) << "Grover Diffusion Operator U_psi_perp applied.\n";
  }
        
  return *this;
}

/**
 * @brief Applies the Grover Phase Oracle (Uf) by negating the amplitude 
 *        of the single specified solution state |w>.
 *
 * @param solution_w The Bitstring index corresponding to the target solution |w>.
 * @return State& Reference to the modified State object.
 */
State& State::grover_oracle_Uf(Bitstring solution_w)
{
        
  // The oracle corresponds to the unitary Uf = I - 2|w><w|, 
  // achieved by phase kickback.
        
  // Iterate over the sparse representation of the quantum state
  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;

    // If the current basis state 'b' is the target solution 'w' (f(w)=1),
    // apply a -1 phase shift (negation).
    if (current_b == solution_w) {
      // Uf|w> = -|w>
      if (grover_verbose() && log_stream_) {
        (*log_stream_) << "Grover Oracle hit at " << current_b << " fliping  " << current_a << " to ";
      }
      current_a *= -1.0;
      if (grover_verbose() && log_stream_) {
        (*log_stream_) << current_a << "\n";
      }
                
      // If we assume a function f(x) that returns 1 only for x=w:
      // Uf|x> = (-1)^f(x) |x>. If f(x)=1, phase is -1.
    }
    // If current_b != solution_w (f(x)=0), amplitude remains unchanged (phase +1).
  }
  if (grover_verbose() && log_stream_) {
    (*log_stream_) << "Grover Oracle Uf applied. Solution state " << solution_w << " phase negated.\n";
  }
  return *this;
}

State& State::grover_oracle_Uf_multi(const std::unordered_set<Bitstring>& solutions)
{
  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;
    if (solutions.count(current_b)) {
      current_a *= -1.0;
    }
  }
  if (grover_verbose() && log_stream_) {
    (*log_stream_) << "Grover Oracle Uf applied to " << solutions.size() << " solution state(s).\n";
  }
  return *this;
}

State& State::grover_oracle_Uf_mask(const std::vector<uint8_t>& solution_mask)
{
  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;
    if (current_b < solution_mask.size() && solution_mask[current_b]) {
      current_a *= -1.0;
    }
  }
  if (grover_verbose() && log_stream_) {
    (*log_stream_) << "Grover Oracle Uf applied to " << solution_mask.size() << " mask entries.\n";
  }
  return *this;
}
