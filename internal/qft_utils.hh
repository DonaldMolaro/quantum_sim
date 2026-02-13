#pragma once
#include <cmath>

inline bool qft_range_valid(int n)
{
  return n > 0 && n < 63;
}

inline bool qft_dimension(int n, unsigned long long& out_N)
{
  if (!qft_range_valid(n)) {
    out_N = 0ULL;
    return false;
  }
  out_N = 1ULL << n;
  return true;
}

inline double qft_phase_exponent(unsigned long long j, unsigned long long k, unsigned long long N, int sign)
{
  const long double PI = std::acos(-1.0L);
  long double exponent = static_cast<long double>(sign) * 2.0L * PI
                         * static_cast<long double>(j) * static_cast<long double>(k)
                         / static_cast<long double>(N);
  return static_cast<double>(exponent);
}

inline double qft_rotation_angle(int r, int sign)
{
  if (!qft_range_valid(r)) {
    return 0.0;
  }
  unsigned long long N = 0ULL;
  qft_dimension(r, N);
  return qft_phase_exponent(1ULL, 1ULL, N, sign);
}
