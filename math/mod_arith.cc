#include "math/mod_arith.hh"

Bitstring gcd_bitstring(Bitstring a, Bitstring b)
{
  while (b != 0) {
    Bitstring t = a % b;
    a = b;
    b = t;
  }
  return a;
}

Bitstring mod_pow(Bitstring base, Bitstring exp, Bitstring mod)
{
  if (mod == 1ULL) return 0ULL;
  Bitstring result = 1ULL;
  Bitstring b = base % mod;
  Bitstring e = exp;
  while (e > 0) {
    if (e & 1ULL) {
      __int128 mul = static_cast<__int128>(result) * static_cast<__int128>(b);
      result = static_cast<Bitstring>(mul % mod);
    }
    __int128 sq = static_cast<__int128>(b) * static_cast<__int128>(b);
    b = static_cast<Bitstring>(sq % mod);
    e >>= 1ULL;
  }
  return result;
}
