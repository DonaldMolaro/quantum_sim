#include "state.hh"
#include <iostream>

Bitstring power_mod(Bitstring base, Bitstring exp, Bitstring mod)
{
    Bitstring result = 1;
    // Ensure base is initially reduced modulo mod
    base %= mod;

    while (exp > 0) {
        // If exp is odd (LSB is 1), include this base term in the result
        if (exp & 1) {
            // Perform multiplication safely, ensuring intermediate results are reduced modulo mod
            result = (result * base) % mod;
        }

        // Square the base for the next iteration and reduce modulo mod
        base = (base * base) % mod;
        
        // Halve the exponent (bit shift right)
        exp >>= 1;
    }
    return result;
}


Bitstring extract_bits(Bitstring b, int start, int end)
{
    // Calculate the length of the segment to extract
    int length = end - start + 1;
    
    // 1. Create a mask containing 'length' number of ones.
    // (1ULL << length) creates a 1 followed by 'length' zeros. Subtracting 1 makes a bitmask of 'length' ones.
    Bitstring mask = (1ULL << length) - 1;

    // 2. Shift the original bitstring right by 'start' positions to align the desired section to index 0
    b >>= start;

    // 3. Apply the mask to keep only the 'length' bits
    return b & mask;
}

Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val) {
    int length = end - start + 1;

    // 1. Construct a mask to clear the bits in the target region of the original bitstring 'b'.
    // Mask of ones in the target position: ((1ULL << length) - 1) << start
    // Invert to get a mask of zeros in the target region: ~mask_ones
    Bitstring clear_mask = ~(((1ULL << length) - 1) << start);

    // 2. Clear the target section in the original bitstring 'b'
    b &= clear_mask;

    // 3. Shift the new value into the correct position
    Bitstring shifted_new_val = new_val << start;

    // 4. Use bitwise OR to insert the new value
    return b | shifted_new_val;
}

#include <complex>
#include <map>

// Assuming QuantumState is defined as std::map<Bitstring, ComplexNumber>
using ComplexNumber = std::complex<double>; 

// This function must be part of the State class or have direct access to its internal state map
void accumulate(std::map<Bitstring, ComplexNumber>& state_map, 
                Bitstring b, 
                ComplexNumber a) 
{
    // C++ map/unordered_map lookup performs summation of amplitudes (if key exists) 
    // or insertion (if key is new).
    // This is the functional equivalent of "reduceByKey(lambda x, y. x + y)" applied 
    // to collect amplitudes for identical keys (bitstrings).
    state_map[b] += a; 
}
