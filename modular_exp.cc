#include "state.hh"

#include <map>

// Assume ComplexNumber, Bitstring, and the sparse representation QuantumState
// (e.g., std::map<Bitstring, ComplexNumber>) are defined.

// --- ASSUMED HELPER FUNCTIONS (Internal to State Class Logic) ---
// Bitstring power_mod(Bitstring base, Bitstring exp, Bitstring mod) 
// Bitstring extract_bits(Bitstring b, int start, int end) 
// Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val) 
// void accumulate(Bitstring b, ComplexNumber a) 

using AmplitudeMap = std::unordered_map<Bitstring, ComplexNumber>;

State& State::controlled_modular_exponentiation(
    int control_qubit, 
    int target_start, int target_end, 
    Bitstring a, Bitstring N, 
    Bitstring power) 
{
    // 1. Classical Pre-calculation
    // CME applies multiplication by R = a^power mod N
    // This value R is a classical constant throughout this unitary operation.
    // CME requires O((log N)^2 * polylog(log N)) classical time if decomposed.
    // NOTE: This relies on a highly efficient classical modular exponentiation function.
    Bitstring R = power_mod(a, power, N); 
    
    std::cout << "Applying Controlled Modular Exponentiation (a^power mod N = " 
              << R << ") controlled by qubit " << control_qubit << ".\n";

    // Prepare a temporary structure to hold the resulting state (b', a')
    AmplitudeMap new_state_map;

    // 2. Map Operation (Functional Transformation)
    // We iterate over the existing superposition states (b, a) and calculate 
    // the new state (b', a') based on the control qubit value.
    for (const auto& pair : state_) {
        Bitstring current_b = pair.first;
        ComplexNumber current_a = pair.second;

        // Check the control qubit (C) state
        bool C = (current_b >> control_qubit) & 1ULL;

        if (C) {
            // Control is |1>: Apply the modular multiplication U_R
            
            // a. Extract the current target register value (Y)
            Bitstring Y = extract_bits(current_b, target_start, target_end); 

            // b. Calculate the new target value Y' = (Y * R) mod N
            // The multiplication Y * R must preserve unitarity (i.e., be a permutation
            // within the computational subspace [0, N-1]).
            // Since gcd(a, N) = 1, the mapping Y -> Y' is indeed a permutation/reversible.
            
            Bitstring Y_prime;
            if (Y < N) {
                 // CME logic: Multiply and take modulus.
                Y_prime = (Y * R) % N; 
            } else {
                // If Y is outside the encoded space [0, N-1], it is usually defined 
                // to map to itself (identity) or handled by the reversible circuit design.
                Y_prime = Y; 
            }

            // c. Construct the new bitstring b' by replacing Y with Y'
            Bitstring new_b = replace_bits(current_b, target_start, target_end, Y_prime);
            
            // Map: (b, a) -> (b', a). Amplitudes are preserved (no phase kickback 
            // is required here as the phase estimation handles that separately).
            
            // d. Accumulate amplitudes for new_b (in case multiple inputs map to one output)
            new_state_map[new_b] += current_a; 
            
        } else {
            // Control is |0|: State is unchanged. Map: (b, a) -> (b, a)
            new_state_map[current_b] += current_a;
        }
    }
    
    // 3. Update the State
    QuantumState final_state_data;
    const double EPSILON = 1e-15; // Tolerance for floating point zero check
    for (const auto& entry : new_state_map) {
        // Only keep states with non-zero magnitude
        if (std::abs(entry.second) > EPSILON) { 
            final_state_data.push_back(entry); 
        }
    }

    // Replace the old internal state representation (assuming state_ is designed 
    // to be repopulated from this final set of pairs, e.g., if it is a vector)
    state_ = final_state_data; 
    return *this;
}
