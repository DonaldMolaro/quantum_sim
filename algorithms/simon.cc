#include "algorithms/simon.hh"

#include "internal/limits.hh"
#include <algorithm>
#include <vector>

namespace {

inline int parity_dot(Bitstring a, Bitstring b)
{
  Bitstring x = a & b;
  int p = 0;
  while (x != 0ULL) {
    p ^= static_cast<int>(x & 1ULL);
    x >>= 1ULL;
  }
  return p;
}

inline bool gf2_insert_basis(std::vector<Bitstring>& basis, Bitstring v)
{
  for (size_t i = 0; i < basis.size(); ++i) {
    if (basis[i] == 0ULL) continue;
    Bitstring lead = basis[i] & (~basis[i] + 1ULL);
    if ((v & lead) != 0ULL) {
      v ^= basis[i];
    }
  }
  if (v == 0ULL) {
    return false;
  }
  basis.push_back(v);
  std::sort(basis.begin(), basis.end(), [](Bitstring a, Bitstring b) {
    Bitstring la = a & (~a + 1ULL);
    Bitstring lb = b & (~b + 1ULL);
    return la > lb;
  });
  return true;
}

} // namespace

SimonResult run_simon(int n_inputs, Bitstring secret, int shots)
{
  SimonResult result;
  result.n_inputs = n_inputs;
  result.secret = secret;

  if (n_inputs <= 0) {
    result.error = "n_inputs must be > 0";
    return result;
  }
  if (!qsim::limits::valid_bitstring_qubit_count(n_inputs) || n_inputs > 16) {
    result.error = "n_inputs must be in [1, 16] for this simulator";
    return result;
  }

  const Bitstring N = 1ULL << n_inputs;
  if (secret >= N) {
    result.error = "secret is out of range for n_inputs";
    return result;
  }

  int total_shots = shots;
  if (total_shots < 0) {
    total_shots = std::max(2 * n_inputs, 12);
  }
  if (total_shots <= 0) {
    result.error = "shots must be > 0";
    return result;
  }

  std::vector<Bitstring> basis;
  result.equations.clear();
  result.equations.reserve(static_cast<size_t>(total_shots));

  if (secret == 0ULL) {
    for (Bitstring y = 1ULL; y < N && static_cast<int>(result.equations.size()) < total_shots; ++y) {
      if (gf2_insert_basis(basis, y)) {
        result.equations.push_back(y);
      }
    }
  } else {
    for (Bitstring y = 1ULL; y < N && static_cast<int>(result.equations.size()) < total_shots; ++y) {
      if (parity_dot(y, secret) != 0) continue;
      if (gf2_insert_basis(basis, y)) {
        result.equations.push_back(y);
      }
    }
  }

  std::vector<Bitstring> candidates;
  if (secret == 0ULL) {
    candidates.push_back(0ULL);
  } else {
    for (Bitstring cand = 1ULL; cand < N; ++cand) {
      bool ok = true;
      for (size_t i = 0; i < result.equations.size(); ++i) {
        if (parity_dot(result.equations[i], cand) != 0) {
          ok = false;
          break;
        }
      }
      if (ok) candidates.push_back(cand);
    }
  }

  result.recovered_secret = candidates.empty() ? 0ULL : candidates[0];
  result.verified = (result.recovered_secret == secret);
  result.ok = true;
  return result;
}
