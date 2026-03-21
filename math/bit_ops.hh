#pragma once
#include "state.hh"

/// Extract bits [start..end] (inclusive) from b, shifted to position 0.
/// Throws std::invalid_argument if the range is invalid.
Bitstring extract_bits(Bitstring b, int start, int end);

/// Replace bits [start..end] (inclusive) in b with new_val.
/// Throws std::invalid_argument if the range is invalid.
Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val);
