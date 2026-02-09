#include "algorithms/shor_classical.hh"
#include "math/mod_arith.hh"
#include <vector>

static std::vector<Bitstring> continued_fraction(Bitstring numerator, Bitstring denominator)
{
  std::vector<Bitstring> terms;
  Bitstring n = numerator;
  Bitstring d = denominator;
  while (d != 0) {
    Bitstring a = n / d;
    terms.push_back(a);
    Bitstring r = n % d;
    n = d;
    d = r;
  }
  return terms;
}

static void convergents(const std::vector<Bitstring>& terms,
                        std::vector<Bitstring>& out_num,
                        std::vector<Bitstring>& out_den)
{
  out_num.clear();
  out_den.clear();
  Bitstring h1 = 1, h2 = 0;
  Bitstring k1 = 0, k2 = 1;
  for (size_t i = 0; i < terms.size(); ++i) {
    Bitstring a = terms[i];
    Bitstring h = a * h1 + h2;
    Bitstring k = a * k1 + k2;
    out_num.push_back(h);
    out_den.push_back(k);
    h2 = h1;
    h1 = h;
    k2 = k1;
    k1 = k;
  }
}

Bitstring estimate_order(Bitstring measured_x, int n_c, Bitstring a, Bitstring N)
{
  Bitstring q = 1ULL << n_c;
  std::vector<Bitstring> terms = continued_fraction(measured_x, q);

  std::vector<Bitstring> nums;
  std::vector<Bitstring> dens;
  convergents(terms, nums, dens);

  for (size_t i = 0; i < dens.size(); ++i) {
    Bitstring r = dens[i];
    if (r == 0) continue;
    if (r >= N) continue;
    if (mod_pow(a, r, N) == 1ULL) {
      return r;
    }
  }
  return 0ULL;
}
