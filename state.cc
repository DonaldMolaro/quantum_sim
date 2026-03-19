/*
 * Quantum state and gates simulation. The main point of this code
 * is to allow someone to poke away at set of simulated Qubits and see
 * what happens to the system.
 *
 * What is *different* here compared to the "normal" implementaion is
 * the lack of linear algebra the underlying principals of the quantum
 * machine are expressed in "map", "reduce_by_key" and "flatmap_and_reduce".
 *
 * Using those primatives it is possible to implement enough gates so that
 * we have a universal set of gates, which we can then simulate quantum
 * computation on.
 *
 * The most "interesting" thing here is the "measrure" method, which measures
 * the requested Qubit and puts it in a classical "bit", what is interesting here
 * is to watch the entangled state of the quantum machine gradually collapse
 * as the quantumness is taken out of it.
 */

#include "state.hh"
#include "internal/limits.hh"
#include <algorithm> // For std::remove_if
#include <cmath>
#include <complex>
#include <random>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

std::ostream* State::default_log_stream_ = nullptr;

// Single definitions for the global complex constants declared extern in state.hh.
const ComplexNumber ONE_COMPLEX(1.0, 0.0);
const ComplexNumber IMAGINARY_UNIT_I(0.0, 1.0);

// Internal RNG for non-deterministic measurements.
static std::mt19937 s_rng{std::random_device{}()};
static std::uniform_real_distribution<double> s_unit_dist(0.0, 1.0);

void State::seed_rng(unsigned int seed) {
  s_rng.seed(seed);
}
// --- Helper Functions (get_jth_bit, flip_jth_bit, set_jth_bit are inline in state.hh) ---
// --- s_map and s_flatMap_and_reduce are template methods defined in state.hh ---

/**
 * Implements .reduceByKey(λx, y. x+ y)
 * Sums amplitudes for identical bitstrings (keys), used after flatMap.
 */
QuantumState State::reduceByKey(const IntermediateState& intermediate_state) const
{
  std::unordered_map<Bitstring, ComplexNumber> reduced_map;
  reduced_map.reserve(intermediate_state.size());

  // Sum amplitudes
  for (const auto& pair : intermediate_state) {
    reduced_map[pair.first] += pair.second;
  }

  // Convert map back to QuantumState vector
  QuantumState final_state;
  final_state.reserve(reduced_map.size());
  for (const auto& pair : reduced_map) {
    // Filter out elements with zero amplitude
    if (std::abs(pair.second) > qsim::limits::AMPLITUDE_EPSILON) {
      final_state.push_back(pair);
    }
  }
  return final_state;
}


/** Constructor: Initializes the state to the ground state |00...0>. */
State::State(int N, int num_cbits)
  : num_qubits_(N), cbits_(num_cbits, 0), log_stream_(default_log_stream_)
{
  // Start state: bitstring 0, amplitude 1 (the ground state)
  state_.push_back({0ULL, ONE_COMPLEX}); 
}

// --- Public Gate Methods ---

/** X Gate (NOT): s.map(λb, a. (b¬j , a)) */
State& State::x(int j)
{
  s_map([j](const Bitstring& b, const ComplexNumber& a) {
    Bitstring b_prime = flip_jth_bit(b, j);
    return std::make_pair(b_prime, a);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

/** Z Gate (Phase flip): equivalent to RZ(pi). */
State& State::z(int j)
{
  constexpr double pi = qsim::limits::PI;
  return rz(j, pi);
}

/** Y Gate: composite via RZ(pi/2), X, RZ(-pi/2). */
State& State::y(int j)
{
  constexpr double pi = qsim::limits::PI;
  rz(j, pi / 2.0);
  x(j);
  rz(j, -pi / 2.0);
  return *this;
}

/** CZ Gate (Controlled Z): implemented via H on target, CX, H on target. */
State& State::cz(int j_control, int k_target)
{
  h(k_target);
  cx(j_control, k_target);
  h(k_target);
  return *this;
}

/** CNOT Gate (alias): implemented via H on target, CZ, H on target. */
State& State::cnot(int j_control, int k_target)
{
  return cx(j_control, k_target);
}

/** CY Gate (Controlled Y): RZ(-pi/2), CX, RZ(pi/2) on target. */
State& State::cy(int j_control, int k_target)
{
  constexpr double pi = qsim::limits::PI;
  rz(k_target, -pi / 2.0);
  cx(j_control, k_target);
  rz(k_target, pi / 2.0);
  return *this;
}

/** CH Gate (Controlled H): composite decomposition using RZ/RY and CX. */
State& State::ch(int j_control, int k_target)
{
  constexpr double pi = qsim::limits::PI;
  // Decomposition using A, B, C such that A X B X C = H and ABC = I.
  // For H = RZ(0) RY(pi/2) RZ(pi) (up to global phase):
  // A = RY(pi/4)
  // B = RY(-pi/4) RZ(-pi/2)
  // C = RZ(pi/2)
  ry(k_target, pi / 4.0);
  cx(j_control, k_target);
  ry(k_target, -pi / 4.0);
  rz(k_target, -pi / 2.0);
  cx(j_control, k_target);
  rz(k_target, pi / 2.0);
  // Correct relative phase to match H on the target when control is 1.
  cz(j_control, k_target);

  return *this;
}

/** CRZ Gate (Controlled RZ): composite using CX and RZ rotations. */
State& State::crz(int j_control, int k_target, double theta)
{
  rz(k_target, theta / 2.0);
  cx(j_control, k_target);
  rz(k_target, -theta / 2.0);
  cx(j_control, k_target);
  return *this;
}

/** CRX Gate (Controlled RX): composite via H, CRZ, H. */
State& State::crx(int j_control, int k_target, double theta)
{
  h(k_target);
  crz(j_control, k_target, theta);
  h(k_target);
  return *this;
}

/** CRY Gate (Controlled RY): composite using RY and CX. */
State& State::cry(int j_control, int k_target, double theta)
{
  ry(k_target, theta / 2.0);
  cx(j_control, k_target);
  ry(k_target, -theta / 2.0);
  cx(j_control, k_target);
  return *this;
}

/** CU Gate (Controlled U): composite using RZ/RY and CX. */
State& State::cu(int j_control, int k_target, double theta, double phi, double lambda)
{
  // Controlled-U3 decomposition.
  // RZ((phi+lambda)/2) RY(theta/2) CX RY(-theta/2) RZ(-(phi+lambda)/2) CX RZ((lambda-phi)/2)
  rz(k_target, (phi + lambda) / 2.0);
  ry(k_target, theta / 2.0);
  cx(j_control, k_target);
  ry(k_target, -theta / 2.0);
  rz(k_target, -(phi + lambda) / 2.0);
  cx(j_control, k_target);
  rz(k_target, (lambda - phi) / 2.0);
  return *this;
}

/** CCX Gate (Toffoli): composite using H, CNOT, and phase rotations. */
State& State::ccx(int c1, int c2, int target)
{
  constexpr double pi = qsim::limits::PI;
  h(target);
  cnot(c2, target);
  rz(target, -pi / 4.0); // T†
  cnot(c1, target);
  rz(target, pi / 4.0); // T
  cnot(c2, target);
  rz(target, -pi / 4.0); // T†
  cnot(c1, target);
  rz(c2, pi / 4.0); // T on control
  rz(target, pi / 4.0); // T on target
  h(target);
  cnot(c1, c2);
  rz(c1, pi / 4.0); // T on control
  rz(c2, -pi / 4.0); // T† on control
  cnot(c1, c2);
  return *this;
}
    
/** CX Gate (Controlled X): s.map(λb, a. (ite(b_j, b¬k, b), a)) */
State& State::cx(int j_control, int k_target)
{
  s_map([j_control, k_target](const Bitstring& b, const ComplexNumber& a) {
    if (get_jth_bit(b, j_control) == 1) {
      Bitstring b_prime = flip_jth_bit(b, k_target);
      return std::make_pair(b_prime, a);
    }
    return std::make_pair(b, a);
  });
  maybe_apply_depolarizing(j_control);
  maybe_apply_depolarizing(k_target);
  return *this;
}
    
/** S Gate (Phase): s.map(λb, a. (b, a * i^b_j)) */
State& State::s(int j)
{
  s_map([j](const Bitstring& b, const ComplexNumber& a) {
    ComplexNumber phase_factor = (get_jth_bit(b, j) == 1) ? IMAGINARY_UNIT_I : ONE_COMPLEX;
    return std::make_pair(b, a * phase_factor);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

/** T Gate (Phase): s.map(λb, a. (b, a * ((1+i)/sqrt(2))^b_j)) */
State& State::t(int j)
{
  const ComplexNumber T_GATE_CONSTANT(ONE_OVER_SQRT_TWO, ONE_OVER_SQRT_TWO);
  s_map([j, T_GATE_CONSTANT](const Bitstring& b, const ComplexNumber& a) {
    ComplexNumber phase_factor = (get_jth_bit(b, j) == 1) ? T_GATE_CONSTANT : ONE_COMPLEX;
    return std::make_pair(b, a * phase_factor);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

/** Hadamard Gate (Superposition): flatMap().reduceByKey() */
State& State::h(int j)
{
  s_flatMap_and_reduce([j](const Bitstring& b, const ComplexNumber& a, IntermediateState& out) {
    int bj = get_jth_bit(b, j);
    out.emplace_back(set_jth_bit(b, j, 0), a * ONE_OVER_SQRT_TWO);
    double coefficient = (1.0 - 2.0 * (double)bj) * ONE_OVER_SQRT_TWO;
    out.emplace_back(set_jth_bit(b, j, 1), a * coefficient);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

State& State::rx(int j, double theta)
{
  const double c = std::cos(theta / 2.0);
  const double s = std::sin(theta / 2.0);
  const ComplexNumber minus_i_s(0.0, -s);

  s_flatMap_and_reduce([j, c, minus_i_s](const Bitstring& b, const ComplexNumber& a, IntermediateState& out) {
    Bitstring b_flip = flip_jth_bit(b, j);
    out.emplace_back(b, a * c);
    out.emplace_back(b_flip, a * minus_i_s);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

State& State::ry(int j, double theta)
{
  const double c = std::cos(theta / 2.0);
  const double s = std::sin(theta / 2.0);

  s_flatMap_and_reduce([j, c, s](const Bitstring& b, const ComplexNumber& a, IntermediateState& out) {
    int bj = get_jth_bit(b, j);
    Bitstring b_flip = flip_jth_bit(b, j);
    if (bj == 0) {
      out.emplace_back(b, a * c);
      out.emplace_back(b_flip, a * s);
    } else {
      out.emplace_back(b_flip, a * (-s));
      out.emplace_back(b, a * c);
    }
  });
  maybe_apply_depolarizing(j);
  return *this;
}

State& State::rz(int j, double theta)
{
  const ComplexNumber phase0(std::cos(-theta / 2.0), std::sin(-theta / 2.0));
  const ComplexNumber phase1(std::cos(theta / 2.0), std::sin(theta / 2.0));

  s_map([j, phase0, phase1](const Bitstring& b, const ComplexNumber& a) {
    int bj = get_jth_bit(b, j);
    return std::make_pair(b, a * (bj ? phase1 : phase0));
  });
  maybe_apply_depolarizing(j);
  return *this;
}

State& State::ru(int j, double theta, double phi, double lambda)
{
  const double c = std::cos(theta / 2.0);
  const double s = std::sin(theta / 2.0);

  const ComplexNumber eiphi(std::cos(phi), std::sin(phi));
  const ComplexNumber eilambda(std::cos(lambda), std::sin(lambda));
  const ComplexNumber eiphilambda(std::cos(phi + lambda), std::sin(phi + lambda));

  s_flatMap_and_reduce([j, c, s, eiphi, eilambda, eiphilambda](const Bitstring& b,
                                                              const ComplexNumber& a, IntermediateState& out) {
    int bj = get_jth_bit(b, j);
    Bitstring b_flip = flip_jth_bit(b, j);
    if (bj == 0) {
      out.emplace_back(b, a * c);
      out.emplace_back(b_flip, a * (eiphi * s));
    } else {
      out.emplace_back(b_flip, a * (-eilambda * s));
      out.emplace_back(b, a * (eiphilambda * c));
    }
  });
  return *this;
}

State& State::phase_flip_if(const Oracle& predicate)
{
  s_map([&predicate](const Bitstring& b, const ComplexNumber& a) {
    if (predicate(b)) {
      return std::make_pair(b, a * -1.0);
    }
    return std::make_pair(b, a);
  });
  return *this;
}

/** Sdg Gate (S†): applies -i phase to |1>. */
State& State::sdg(int j)
{
  const ComplexNumber MINUS_I(0.0, -1.0);
  s_map([j, MINUS_I](const Bitstring& b, const ComplexNumber& a) {
    ComplexNumber phase_factor = (get_jth_bit(b, j) == 1) ? MINUS_I : ONE_COMPLEX;
    return std::make_pair(b, a * phase_factor);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

/** Tdg Gate (T†): applies (1-i)/sqrt(2) phase to |1>. */
State& State::tdg(int j)
{
  const ComplexNumber TDG_CONSTANT(ONE_OVER_SQRT_TWO, -ONE_OVER_SQRT_TWO);
  s_map([j, TDG_CONSTANT](const Bitstring& b, const ComplexNumber& a) {
    ComplexNumber phase_factor = (get_jth_bit(b, j) == 1) ? TDG_CONSTANT : ONE_COMPLEX;
    return std::make_pair(b, a * phase_factor);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

/** P Gate (arbitrary phase): P(φ)|0>=|0>, P(φ)|1>=e^{iφ}|1>. */
State& State::p(int j, double phi)
{
  const ComplexNumber phase(std::cos(phi), std::sin(phi));
  s_map([j, phase](const Bitstring& b, const ComplexNumber& a) {
    if (get_jth_bit(b, j) == 1) {
      return std::make_pair(b, a * phase);
    }
    return std::make_pair(b, a);
  });
  maybe_apply_depolarizing(j);
  return *this;
}

/** CP Gate (Controlled Phase): e^{iφ} when both control and target are |1>. */
State& State::cp(int j_control, int k_target, double phi)
{
  const ComplexNumber phase(std::cos(phi), std::sin(phi));
  s_map([j_control, k_target, phase](const Bitstring& b, const ComplexNumber& a) {
    if (get_jth_bit(b, j_control) == 1 && get_jth_bit(b, k_target) == 1) {
      return std::make_pair(b, a * phase);
    }
    return std::make_pair(b, a);
  });
  maybe_apply_depolarizing(j_control);
  maybe_apply_depolarizing(k_target);
  return *this;
}

/** CSWAP Gate (Fredkin): swaps k and l when control j is |1>. */
State& State::cswap(int j_control, int k, int l)
{
  s_map([j_control, k, l](const Bitstring& b, const ComplexNumber& a) {
    if (get_jth_bit(b, j_control) == 1) {
      int bk = get_jth_bit(b, k);
      int bl = get_jth_bit(b, l);
      Bitstring b_prime = set_jth_bit(set_jth_bit(b, k, bl), l, bk);
      return std::make_pair(b_prime, a);
    }
    return std::make_pair(b, a);
  });
  return *this;
}

/** MCX Gate (multi-controlled X): flips target when all controls are |1>. */
State& State::mcx(const std::vector<int>& controls, int target)
{
  s_map([&controls, target](const Bitstring& b, const ComplexNumber& a) {
    for (int c : controls) {
      if (get_jth_bit(b, c) == 0) return std::make_pair(b, a);
    }
    return std::make_pair(flip_jth_bit(b, target), a);
  });
  return *this;
}

/** RESET gate: collapses qubit j to |0> by discarding |1> components and renormalising. */
State& State::reset(int j)
{
  double norm_sq = 0.0;
  for (const auto& pair : state_) {
    if (get_jth_bit(pair.first, j) == 0) norm_sq += std::norm(pair.second);
  }
  if (norm_sq < qsim::limits::AMPLITUDE_EPSILON) {
    // All amplitude is on |1>_j — flip all those bits to 0.
    QuantumState next;
    next.reserve(state_.size());
    for (const auto& pair : state_)
      next.push_back({set_jth_bit(pair.first, j, 0), pair.second});
    state_ = reduceByKey(next);
  } else {
    const double scale = 1.0 / std::sqrt(norm_sq);
    QuantumState next;
    for (const auto& pair : state_)
      if (get_jth_bit(pair.first, j) == 0)
        next.push_back({pair.first, pair.second * scale});
    state_ = next;
  }
  return *this;
}

/** iSWAP: swaps j and k and multiplies off-diagonal swapped terms by i. */
State& State::iswap(int j, int k)
{
  s_flatMap_and_reduce([j, k](const Bitstring& b, const ComplexNumber& a, IntermediateState& out) {
    const int bj = get_jth_bit(b, j);
    const int bk = get_jth_bit(b, k);
    if (bj == bk) {
      out.emplace_back(b, a);
    } else {
      out.emplace_back(set_jth_bit(set_jth_bit(b, j, bk), k, bj), a * ComplexNumber(0.0, 1.0));
    }
  });
  maybe_apply_depolarizing(j);
  maybe_apply_depolarizing(k);
  return *this;
}

/** XX(θ): exp(-i θ/2 X⊗X) — equal superposition of b and flip-both(b). */
State& State::xx(int j, int k, double theta)
{
  const double c = std::cos(theta / 2.0);
  const double s = std::sin(theta / 2.0);
  s_flatMap_and_reduce([j, k, c, s](const Bitstring& b, const ComplexNumber& a, IntermediateState& out) {
    out.emplace_back(b,                                  a * c);
    out.emplace_back(flip_jth_bit(flip_jth_bit(b,j),k), a * ComplexNumber(0.0, -s));
  });
  maybe_apply_depolarizing(j);
  maybe_apply_depolarizing(k);
  return *this;
}

/** YY(θ): exp(-i θ/2 Y⊗Y) */
State& State::yy(int j, int k, double theta)
{
  const double c = std::cos(theta / 2.0);
  const double s = std::sin(theta / 2.0);
  s_flatMap_and_reduce([j, k, c, s](const Bitstring& b, const ComplexNumber& a, IntermediateState& out) {
    const int bj = get_jth_bit(b, j);
    const int bk = get_jth_bit(b, k);
    // Y⊗Y: |00>->-|11>, |01>->+|10>, |10>->+|01>, |11>->-|00>
    const double yy_sign = ((bj ^ bk) == 0) ? -1.0 : 1.0;
    out.emplace_back(b,                                   a * c);
    out.emplace_back(flip_jth_bit(flip_jth_bit(b,j),k),  a * ComplexNumber(0.0, -s * yy_sign));
  });
  maybe_apply_depolarizing(j);
  maybe_apply_depolarizing(k);
  return *this;
}

/** ZZ(θ): exp(-i θ/2 Z⊗Z) — diagonal, no state mixing. */
State& State::zz(int j, int k, double theta)
{
  const ComplexNumber phase_same(std::cos(theta / 2.0), -std::sin(theta / 2.0)); // same parity
  const ComplexNumber phase_diff(std::cos(theta / 2.0),  std::sin(theta / 2.0)); // diff parity
  s_map([j, k, phase_same, phase_diff](const Bitstring& b, const ComplexNumber& a) {
    const ComplexNumber ph = (get_jth_bit(b,j) == get_jth_bit(b,k)) ? phase_same : phase_diff;
    return std::make_pair(b, a * ph);
  });
  maybe_apply_depolarizing(j);
  maybe_apply_depolarizing(k);
  return *this;
}

// --- Analysis methods (non-mutating) -----------------------------------------

State::BlochVector State::bloch(int j) const
{
  std::unordered_map<Bitstring, ComplexNumber> amp0, amp1;
  for (const auto& pair : state_) {
    const Bitstring b_other = pair.first & ~(1ULL << j);
    if (get_jth_bit(pair.first, j) == 0) amp0[b_other] += pair.second;
    else                                  amp1[b_other] += pair.second;
  }
  double rho00 = 0.0, rho11 = 0.0;
  ComplexNumber rho01(0.0, 0.0);
  for (const auto& p : amp0) rho00 += std::norm(p.second);
  for (const auto& p : amp1) rho11 += std::norm(p.second);
  for (const auto& p0 : amp0) {
    auto it = amp1.find(p0.first);
    if (it != amp1.end()) rho01 += std::conj(p0.second) * it->second;
  }
  return {2.0 * rho01.real(), 2.0 * rho01.imag(), rho00 - rho11};
}

double State::expect_pauli(const std::vector<std::pair<char,int>>& ops) const
{
  ComplexNumber total(0.0, 0.0);
  for (const auto& kv : state_) {
    Bitstring b_out = kv.first;
    ComplexNumber phase(1.0, 0.0);
    for (const auto& op : ops) {
      const char p = op.first;
      const int  q = op.second;
      const bool bit1 = (b_out >> q) & 1ULL;
      if (p == 'X') {
        b_out ^= (1ULL << q);
      } else if (p == 'Y') {
        phase *= bit1 ? ComplexNumber(0.0, -1.0) : ComplexNumber(0.0, 1.0);
        b_out ^= (1ULL << q);
      } else if (p == 'Z') {
        if (bit1) phase *= -1.0;
      }
    }
    total += std::conj(kv.second) * phase * get_amplitude(b_out);
  }
  return total.real();
}

double State::entropy(int start_q, int end_q) const
{
  if (start_q > end_q) std::swap(start_q, end_q);
  const int n_A = end_q - start_q + 1;
  if (n_A < 1 || n_A > 10) return 0.0;
  const int dim_A = 1 << n_A;

  // Build mask for subsystem A bits.
  Bitstring mask_A = 0;
  for (int q = start_q; q <= end_q; ++q) mask_A |= (1ULL << q);

  // Group amplitude by complement bitstring → vector indexed by subsystem index.
  std::unordered_map<Bitstring, std::vector<ComplexNumber>> cols;
  for (const auto& kv : state_) {
    const int alpha = static_cast<int>((kv.first & mask_A) >> start_q);
    const Bitstring b_comp = kv.first & ~mask_A;
    auto& col = cols[b_comp];
    if (col.empty()) col.resize(dim_A, ComplexNumber(0.0, 0.0));
    col[alpha] += kv.second;
  }

  // Accumulate ρ_A[α][β] = Σ_{bB} conj(a[α,bB]) * a[β,bB]
  std::vector<std::vector<ComplexNumber>> rho(dim_A, std::vector<ComplexNumber>(dim_A, 0.0));
  for (const auto& entry : cols) {
    const auto& col = entry.second;
    for (int a2 = 0; a2 < dim_A; ++a2)
      for (int b2 = 0; b2 < dim_A; ++b2)
        rho[a2][b2] += std::conj(col[a2]) * col[b2];
  }

  // Eigenvalues via power deflation (ρ_A is Hermitian PSD, eigenvalues ≥ 0).
  auto rho_work = rho;
  double S = 0.0;
  for (int k = 0; k < dim_A; ++k) {
    std::vector<ComplexNumber> v(dim_A, ComplexNumber(0.0, 0.0));
    v[k % dim_A] = ComplexNumber(1.0, 0.0);
    double lam = 0.0;
    for (int iter = 0; iter < 500; ++iter) {
      std::vector<ComplexNumber> mv(dim_A, ComplexNumber(0.0, 0.0));
      for (int i = 0; i < dim_A; ++i)
        for (int j2 = 0; j2 < dim_A; ++j2)
          mv[i] += rho_work[i][j2] * v[j2];
      double norm2 = 0.0;
      ComplexNumber dot(0.0, 0.0);
      for (int i = 0; i < dim_A; ++i) { dot += std::conj(v[i]) * mv[i]; norm2 += std::norm(mv[i]); }
      lam = dot.real();
      if (norm2 < 1e-30) break;
      const double sc = 1.0 / std::sqrt(norm2);
      for (int i = 0; i < dim_A; ++i) v[i] = mv[i] * sc;
    }
    if (lam > 1e-12) S -= lam * std::log2(lam);
    for (int i = 0; i < dim_A; ++i)
      for (int j2 = 0; j2 < dim_A; ++j2)
        rho_work[i][j2] -= lam * v[i] * std::conj(v[j2]);
  }
  return std::max(S, 0.0);
}

/** Applies a stochastic single-qubit Pauli error to qubit j with probability noise_prob_. */
void State::maybe_apply_depolarizing(int j)
{
  if (noise_prob_ <= 0.0) return;
  const double r = s_unit_dist(s_rng);
  if (r >= noise_prob_) return;
  const double which = s_unit_dist(s_rng);
  if (which < 1.0 / 3.0) {
    // X error: flip bit j
    s_map([j](const Bitstring& b, const ComplexNumber& a) {
      return std::make_pair(flip_jth_bit(b, j), a);
    });
  } else if (which < 2.0 / 3.0) {
    // Y error: flip bit j and apply phase
    s_map([j](const Bitstring& b, const ComplexNumber& a) {
      ComplexNumber phase = (get_jth_bit(b, j) == 0)
                                ? ComplexNumber(0.0, 1.0)
                                : ComplexNumber(0.0, -1.0);
      return std::make_pair(flip_jth_bit(b, j), a * phase);
    });
  } else {
    // Z error: flip phase on |1>
    s_map([j](const Bitstring& b, const ComplexNumber& a) {
      ComplexNumber phase = (get_jth_bit(b, j) == 1) ? ComplexNumber(-1.0, 0.0) : ONE_COMPLEX;
      return std::make_pair(b, a * phase);
    });
  }
}


// --- Measurement Method ---

/**
 * Computes the probability of measuring 0 for qubit j: 
 * s.filter(λb, a. bj = 0) .map(λb, a. |a|^2) .sum()
 */
double State::compute_probability_of_0(int j) const
{
  double probability = 0.0;
        
  // 1. Filter and Map (in one loop for efficiency): 
  // Iterate through states and accumulate squared magnitudes (|a|^2) if b_j = 0.
  for (const auto& pair : state_) {
    if (get_jth_bit(pair.first, j) == 0) {
      // Sum: using std::norm to calculate |a|^2
      probability += std::norm(pair.second);
    }
  }
  return probability;
}

// Accessor for observation
const QuantumState& State::get_state() const { return state_; }

/**
 * @brief Measures the jth qubit, collapses the state, and stores the result.
 * 
 * @param j Qubit index to measure.
 * @param cbit_index Index of classical register to store the result.
 * @return Reference to the updated State object.
 */
State& State::measure_with_rng(int j, unsigned long cbit_index, double random_val) {
  // 1. Compute Pj,0: Probability of measuring 0
  double p0 = compute_probability_of_0(j);
  // 2. Determine the outcome (0 or 1)
  int outcome; 
  double p_outcome; // The probability corresponding to the chosen outcome
        
  if (random_val <= p0) {
    outcome = 0;
    p_outcome = p0;
  } else {
    outcome = 1;
    p_outcome = 1.0 - p0; // Pj,1
  }
        
  // Check probability stability (ensure p_outcome is not zero)
  if (p_outcome < qsim::limits::AMPLITUDE_EPSILON) {
    p_outcome = 1.0; 
  }

  // 3. Calculate the renormalization factor (1 / sqrt(p))
  ComplexNumber norm_factor(1.0 / std::sqrt(p_outcome), 0.0);
        
  // 4. Post-measurement State Update (Collapse and Renormalize)
  s_map([j, outcome, norm_factor](const Bitstring& b, const ComplexNumber& a) {
    int bj = get_jth_bit(b, j);
    if (bj == outcome) {
      return std::make_pair(b, a * norm_factor);
    } else {
      return std::make_pair(b, ComplexNumber(0.0, 0.0));
    }
  });

  // 5. Clean up the state vector by removing elements with zero amplitude
  state_.erase(std::remove_if(state_.begin(), state_.end(), 
			      [](const QubitAmplitudePair& pair) { 
				return std::abs(pair.second) < qsim::limits::AMPLITUDE_EPSILON; 
			      }), 
	       state_.end());

  // 6. Store the result in the classical register
  if (!cbits_.empty())
    {
      if (cbit_index < cbits_.size()) {
	cbits_[cbit_index] = outcome;
      }
    }
  else
    {
      // No classical register; measurement result is not stored.
    }
  return *this;
}

State& State::measure(int j, unsigned long cbit_index) {
  double random_val = s_unit_dist(s_rng);
  return measure_with_rng(j, cbit_index, random_val);
}

State& State::measure_all_with_rng(const std::vector<double>& random_vals, std::vector<int>& out) {
  out.clear();
  out.resize(num_qubits_, 0);

  for (int j = 0; j < num_qubits_; ++j) {
    double rv = (j < (int)random_vals.size()) ? random_vals[j] : s_unit_dist(s_rng);
    unsigned long cbit_index = cbits_.empty() ? 0 : static_cast<unsigned long>(j);
    measure_with_rng(j, cbit_index, rv);

    if (!cbits_.empty() && j < (int)cbits_.size()) {
      out[j] = cbits_[j];
    } else if (!state_.empty()) {
      out[j] = get_jth_bit(state_.front().first, j);
    }
  }

  return *this;
}

State& State::measure_all(std::vector<int>& out) {
  std::vector<double> empty;
  return measure_all_with_rng(empty, out);
}
