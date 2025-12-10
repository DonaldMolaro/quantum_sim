#include "state.hh"    
/**
 * @brief Applies the SWAP gate between qubit j and qubit k.
 * This operation swaps the bit values (keys) in the quantum state representation 
 * but leaves the amplitudes unchanged, mirroring the core (b, a) -> (b', a) map.
 * 
 * @param j Index of the first qubit (0-indexed, LSB).
 * @param k Index of the second qubit.
 * @return State& Reference to the modified State object.
 */
State& State::swap(int j, int k)
{
  if (j == k) return *this;
        
  // Iterate over the sparse representation of the quantum state
  for (auto& pair : state_) {
    Bitstring& current_b = pair.first;
            
    // 1. Extract the bit values at positions j and k
    // We use 1ULL (unsigned long long 1) to define the masks.
    unsigned long long bit_j = (current_b >> j) & 1ULL;
    unsigned long long bit_k = (current_b >> k) & 1ULL;

    // 2. Only perform swaps if the bits differ (XOR mask technique)
    if (bit_j != bit_k) {
      // The swap requires flipping both bits.
      // Flip bit j: XOR with mask (1 shifted left by j)
      current_b ^= (1ULL << j);
                
      // Flip bit k: XOR with mask (1 shifted left by k)
      current_b ^= (1ULL << k);
    }
  }
        
  std::cout << "SWAP gate applied between qubit " << j << " and qubit " << k << ".\n";
  return *this;
}
