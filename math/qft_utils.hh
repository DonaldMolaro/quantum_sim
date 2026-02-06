#pragma once
#include <cmath>

inline bool qft_range_valid(int n)
{
  return n > 0 && n < 63;
}

inline double qft_phase_exponent(unsigned long long j, unsigned long long k, unsigned long long N, int sign)
{
  const long double PI = std::acos(-1.0L);
  long double exponent = static_cast<long double>(sign) * 2.0L * PI
                         * static_cast<long double>(j) * static_cast<long double>(k)
                         / static_cast<long double>(N);
  return static_cast<double>(exponent);
}
