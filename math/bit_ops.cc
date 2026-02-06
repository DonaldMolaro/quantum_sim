#include "math/bit_ops.hh"

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
