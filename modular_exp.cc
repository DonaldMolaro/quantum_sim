#include "modular_exp.hh"

#include <map>

// Assume ComplexNumber, Bitstring, and the sparse representation QuantumState
// (e.g., std::map<Bitstring, ComplexNumber>) are defined.

// --- ASSUMED HELPER FUNCTIONS (Internal to State Class Logic) ---
// Bitstring power_mod(Bitstring base, Bitstring exp, Bitstring mod) 
// Bitstring extract_bits(Bitstring b, int start, int end) 
// Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val) 
// void accumulate(Bitstring b, ComplexNumber a) 

using AmplitudeMap = std::unordered_map<Bitstring, ComplexNumber>;

const double EPSILON = 1e-9; 



/**
 * Creates a mask for a contiguous block of bits spanning [start, end].
 */
Bitstring create_mask(int start, int end) {
  if (start > end) return 0;
  int length = end - start + 1;
  // Creates a block of 'length' ones, then shifts it to the 'start' position.
  return ((1ULL << length) - 1) << start;
}

/**
 * Extracts the numerical value of the target register from a full bitstring.
 */
Bitstring get_target_value_general(Bitstring b, int start, int end) {
  Bitstring target_mask = create_mask(start, end);
  // Apply the mask and shift the result back to the 0 position
  return (b & target_mask) >> start;
}

/**
 * Replaces the numerical value of the target register within a full bitstring.
 */
Bitstring replace_target_value_general(Bitstring original, Bitstring new_target, int start, int end)
{
  Bitstring target_mask = create_mask(start, end);
  // 1. Clear the existing target bits
  Bitstring cleared = original & (~target_mask);
  // 2. Place the new target value (shifted into position)
  Bitstring shifted_new = new_target << start;
  // 3. Combine cleared and new values
  return cleared | shifted_new;
}


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
        
  QuantumState new_state_set;
  Bitstring control_mask = 1ULL << control_qubit;
        
  // Step 1: Calculate the classical exponentiation result once.
  Bitstring exponent_factor = modular_exponentiation(a, power, N);

  // Step 2: Apply the transformation to every element in the superposition (map operation).
  for (auto const& pair : state_) {
    Bitstring b = pair.first;
    ComplexNumber amplitude = pair.second;
            
    if (std::abs(amplitude) < EPSILON) continue;

    Bitstring control_value = (b & control_mask) >> control_qubit;
            
    Bitstring b_new;

    if (control_value == 1) {
      // Control ON: Apply the function U^P to the target register.
                
      // Extract the target register's numerical value (y).
      Bitstring current_target_value = get_target_value_general(b, target_start, target_end);

      // For unitarity: if y is outside [0, N-1], leave it unchanged.
      if (current_target_value >= N) {
        b_new = b;
      } else {
        // Calculate the new target value: y' = (F * y) mod N.
        Bitstring new_target_value = (exponent_factor * current_target_value) % N;

        // Construct the new full bitstring (b' = Control | Target_new).
        b_new = replace_target_value_general(b, new_target_value, target_start, target_end);
      }
    } else {
      // Control OFF: Apply Identity (keep state unchanged).
      b_new = b;
    }
            
    // Step 3: Aggregate the amplitude to the resulting state (reduceByKey operation).
    QubitAmplitudePair* entry = find_or_add(new_state_set, b_new);
    entry->second += amplitude;
  }

  // Update the state to the result of the linear transformation.
  state_ = new_state_set;
  return *this;
}
