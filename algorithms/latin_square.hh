#pragma once
#include "state.hh"

// Assignment encodes row 1 and row 2, 6 cells total, 2 bits per cell.
// Row 0 is provided explicitly (must be a permutation of 0,1,2).
bool is_valid_latin3(Bitstring assignment, const int row0[3]);
bool is_valid_latin3_fixedrow(Bitstring assignment);
