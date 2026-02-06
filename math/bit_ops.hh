#pragma once
#include "state.hh"

Bitstring extract_bits(Bitstring b, int start, int end);
Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val);
