#include "math/bit_ops.hh"
#include "internal/limits.hh"
#include <stdexcept>
#include <string>

static void validate_bit_range(int start, int end, const char* func)
{
    if (start < 0 || end < 0 || start > end ||
        (end - start + 1) > qsim::limits::kMaxBitstringQubits) {
        throw std::invalid_argument(
            std::string(func) + ": invalid bit range [" +
            std::to_string(start) + ", " + std::to_string(end) + "]");
    }
}

Bitstring extract_bits(Bitstring b, int start, int end)
{
    validate_bit_range(start, end, "extract_bits");

    int length = end - start + 1;
    Bitstring mask = (1ULL << length) - 1;
    b >>= start;
    return b & mask;
}

Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val) {
    validate_bit_range(start, end, "replace_bits");

    int length = end - start + 1;
    Bitstring clear_mask = ~(((1ULL << length) - 1) << start);
    b &= clear_mask;
    Bitstring shifted_new_val = new_val << start;
    return b | shifted_new_val;
}
