#pragma once

#include <cstdint>
#include <vector>

// Generate random bits using H gates and measurement.
// Bits are returned LSB-first.
std::vector<int> qrng_bits(int n, const std::vector<double>* random_vals = nullptr);

// Generate up to 63 bits of randomness as a uint64_t (LSB-first).
uint64_t qrng_u64(int n, const std::vector<double>* random_vals = nullptr);
