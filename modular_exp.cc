#include "modular_exp.hh"
#include "internal/limits.hh"
#include "math/bit_ops.hh"

#include <unordered_map>

using AmplitudeMap = std::unordered_map<Bitstring, ComplexNumber>;


// Classical Modular Expoentiation

Bitstring modular_exponentiation(Bitstring a, Bitstring power, Bitstring N)
{
  // Handle the case where the base is 0 (assuming N > 1, 0^0 is 1).
  if (a == 0 && N > 1) return (power == 0) ? 1ULL : 0ULL;
  if (a == 0 && N == 1) return 0ULL; // Modulo 1 always yields 0.

  if (power == 0) return 1ULL;

  Bitstring result = 1ULL;
  Bitstring base = a % N;

  while (power > 0) {
    // If power is odd, multiply result by base (mod N)
    if (power & 1) {
      result = (result * base) % N;
    }

    // Square the base (mod N) and halve the power
    base = (base * base) % N;
    power >>= 1;
  }
  return result;
}
/**
 * Complete and generalized implementation of Controlled Modular Exponentiation.
 * This performs the conditional unitary operation U_a^power |y> = |a^power * y mod N>.
 */
State& State::controlled_modular_exponentiation(
						int control_qubit, 
						int target_start, int target_end, 
						Bitstring a, Bitstring N, 
						Bitstring power)
{
        
  AmplitudeMap accum;
  accum.reserve(state_.size());
  Bitstring control_mask = 1ULL << control_qubit;

  // Step 1: Calculate the classical exponentiation result once.
  Bitstring exponent_factor = modular_exponentiation(a, power, N);

  // Step 2: Apply the transformation to every element in the superposition (map operation).
  for (auto const& pair : state_) {
    Bitstring b = pair.first;
    ComplexNumber amplitude = pair.second;

    if (std::abs(amplitude) < qsim::limits::AMPLITUDE_EPSILON) continue;

    Bitstring control_value = (b & control_mask) >> control_qubit;

    Bitstring b_new;
    if (control_value == 1) {
      // Control ON: Apply the function U^P to the target register.
      Bitstring current_target_value = extract_bits(b, target_start, target_end);
      // For unitarity: if y is outside [0, N-1], leave it unchanged.
      if (current_target_value >= N) {
        b_new = b;
      } else {
        Bitstring new_target_value = (exponent_factor * current_target_value) % N;
        b_new = replace_bits(b, target_start, target_end, new_target_value);
      }
    } else {
      // Control OFF: Apply Identity (keep state unchanged).
      b_new = b;
    }

    // Step 3: Aggregate the amplitude (O(1) amortized via unordered_map).
    accum[b_new] += amplitude;
  }

  // Convert accumulator back to QuantumState, filtering near-zero amplitudes.
  state_.clear();
  state_.reserve(accum.size());
  for (const auto& kv : accum) {
    if (std::abs(kv.second) > qsim::limits::AMPLITUDE_EPSILON) {
      state_.emplace_back(kv.first, kv.second);
    }
  }
  return *this;
}
