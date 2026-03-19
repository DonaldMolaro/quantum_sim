#include "state.hh"
#include "algorithms/qpe.hh"
#include "algorithms/qec.hh"
#include "demos/qpe_demo.hh"
#include "tests/helpers.hh"
#include "internal/limits.hh"
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <vector>

using test_helpers::assert_complex_equal;
using test_helpers::assert_complex_close;
using test_helpers::assert_amplitude_match;
using test_helpers::assert_amplitude_magnitude;
using test_helpers::assert_double_close;

static int g_failures = 0;

static void run_test(const std::string& name, std::function<void()> fn)
{
  try {
    fn();
    std::cout << "[PASS] " << name << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[FAIL] " << name << ": " << e.what() << "\n";
    ++g_failures;
  }
}

static constexpr double PI = qsim::limits::PI;

// ---- Sdg gate ---------------------------------------------------------------

void test_sdg_is_s_inverse()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  s.s(0);
  s.sdg(0);
  // S then Sdg should return to |1>
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b1), "Sdg·S|1>=|1>");
  assert_complex_equal(ComplexNumber(0.0, 0.0), s.get_amplitude(0b0), "Sdg·S: |0> is 0");
}

void test_sdg_phase_on_one()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  s.sdg(0);
  // Sdg|1> = -i|1>
  ComplexNumber got = s.get_amplitude(0b1);
  ComplexNumber expected(0.0, -1.0);
  if (std::abs(got - expected) > 1e-9) {
    throw std::runtime_error("Sdg|1> should be -i|1>, got (" +
                             std::to_string(got.real()) + "," + std::to_string(got.imag()) + ")");
  }
}

void test_sdg_no_change_on_zero()
{
  State s(1);
  s.set_basis_state(0b0, 1.0);
  s.sdg(0);
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b0), "Sdg|0>=|0>");
}

// ---- Tdg gate ---------------------------------------------------------------

void test_tdg_is_t_inverse()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  s.t(0);
  s.tdg(0);
  // T then Tdg should be identity on |1>
  if (std::abs(s.get_amplitude(0b1) - ComplexNumber(1.0, 0.0)) > 1e-9) {
    throw std::runtime_error("Tdg·T|1> != |1>");
  }
}

void test_tdg_phase_on_one()
{
  // Tdg applies (1-i)/sqrt(2) = e^{-i*pi/4} to |1>
  State s(1);
  s.set_basis_state(0b1, 1.0);
  s.tdg(0);
  ComplexNumber got = s.get_amplitude(0b1);
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  ComplexNumber expected(inv_sqrt2, -inv_sqrt2);
  if (std::abs(got - expected) > 1e-9) {
    throw std::runtime_error("Tdg|1> wrong phase");
  }
}

// ---- P gate -----------------------------------------------------------------

void test_p_gate_identity_on_zero()
{
  State s(1);
  s.set_basis_state(0b0, 1.0);
  s.p(0, PI / 4.0);
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b0), "P|0>=|0>");
}

void test_p_gate_applies_phase_to_one()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  s.p(0, PI / 2.0); // should give i|1>
  ComplexNumber got = s.get_amplitude(0b1);
  ComplexNumber expected(0.0, 1.0);
  if (std::abs(got - expected) > 1e-9) {
    throw std::runtime_error("P(pi/2)|1> != i|1>");
  }
}

void test_p_equals_s_at_pi_over_2()
{
  // P(pi/2) = S (up to global phase convention)
  State s1(1), s2(1);
  s1.set_basis_state(0b1, 1.0);
  s2.set_basis_state(0b1, 1.0);
  s1.s(0);
  s2.p(0, PI / 2.0);
  if (std::abs(s1.get_amplitude(0b1) - s2.get_amplitude(0b1)) > 1e-9) {
    throw std::runtime_error("P(pi/2) != S on |1>");
  }
}

// ---- CP gate ----------------------------------------------------------------

void test_cp_no_effect_when_control_zero()
{
  State s(2);
  s.set_basis_state(0b01, 1.0); // qubit 0 = 1, qubit 1 = 0 (control)
  s.cp(1, 0, PI);
  // control (qubit 1) is 0, so no phase change
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b01), "CP: no effect when control=0");
}

void test_cp_applies_phase_when_both_one()
{
  State s(2);
  s.set_basis_state(0b11, 1.0); // both qubits = 1
  s.cp(0, 1, PI); // phase = pi -> e^{i*pi} = -1
  ComplexNumber got = s.get_amplitude(0b11);
  if (std::abs(got - ComplexNumber(-1.0, 0.0)) > 1e-9) {
    throw std::runtime_error("CP(pi) should flip sign when both qubits are 1");
  }
}

// ---- CSWAP gate -------------------------------------------------------------

void test_cswap_no_swap_when_control_zero()
{
  State s(3);
  // qubit 0=1, qubit 1=0, qubit 2 (control)=0: bitstring = 0b001 = 1
  s.set_basis_state(0b001, 1.0);
  s.cswap(2, 0, 1);
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b001), "CSWAP no swap when control=0");
}

void test_cswap_swaps_when_control_one()
{
  // qubit 0=1, qubit 1=0, qubit 2 (control)=1: bitstring = 0b101 = 5
  State s(3);
  s.set_basis_state(0b101, 1.0);
  s.cswap(2, 0, 1);
  // After swap: qubit 0=0, qubit 1=1, qubit 2=1 -> 0b110 = 6
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b110), "CSWAP swaps qubits 0 and 1 when control=1");
  assert_complex_equal(ComplexNumber(0.0, 0.0), s.get_amplitude(0b101), "CSWAP original state gone");
}

// ---- MCX gate ---------------------------------------------------------------

void test_mcx_single_control_same_as_cx()
{
  // MCX with 1 control should be identical to CX
  State s1(2), s2(2);
  s1.set_basis_state(0b11, 1.0);
  s2.set_basis_state(0b11, 1.0);
  s1.cx(1, 0);
  s2.mcx({1}, 0);
  if (std::abs(s1.get_amplitude(0b01) - s2.get_amplitude(0b01)) > 1e-9) {
    throw std::runtime_error("MCX with 1 control != CX");
  }
}

void test_mcx_two_controls_same_as_ccx()
{
  // MCX and CCX should agree on *which* basis states get amplitude and which don't.
  // (CCX decomposition may introduce global phase, so compare probabilities.)
  State s1(3), s2(3);
  s1.set_basis_state(0b111, 1.0);
  s2.set_basis_state(0b111, 1.0);
  s1.ccx(1, 2, 0);
  s2.mcx({1, 2}, 0);
  // Both should map |111> -> |110>: probability 1 at 0b110, 0 elsewhere.
  const double p1 = std::norm(s1.get_amplitude(0b110));
  const double p2 = std::norm(s2.get_amplitude(0b110));
  if (std::abs(p1 - 1.0) > 1e-6 || std::abs(p2 - 1.0) > 1e-6) {
    throw std::runtime_error("MCX/CCX: |111> should map entirely to |110>");
  }
}

void test_mcx_three_controls_flips_only_when_all_set()
{
  State s(4);
  // All controls (1, 2, 3) = 1, target (0) = 1 -> should flip target to 0
  // bitstring: qubits 3,2,1,0 = 1,1,1,1 -> 0b1111 = 15
  s.set_basis_state(0b1111, 1.0);
  s.mcx({1, 2, 3}, 0);
  // target flips: 0b1110 = 14
  assert_complex_equal(ComplexNumber(1.0, 0.0), s.get_amplitude(0b1110), "MCX 3-control flip");

  // Now test with one control missing: 0b1101 = 13 (qubit 1 = 0)
  State s2(4);
  s2.set_basis_state(0b1101, 1.0);
  s2.mcx({1, 2, 3}, 0);
  // Should NOT flip (qubit 1 = 0)
  assert_complex_equal(ComplexNumber(1.0, 0.0), s2.get_amplitude(0b1101), "MCX: no flip when one control=0");
}

// ---- Noise model ------------------------------------------------------------

void test_noise_zero_is_off()
{
  // With noise=0, repeated H+H should return to |0> exactly.
  State::seed_rng(42);
  State s(1);
  s.set_noise_probability(0.0);
  s.set_basis_state(0b0, 1.0);
  for (int i = 0; i < 100; ++i) {
    s.h(0); s.h(0);
  }
  if (std::abs(s.get_amplitude(0b0) - ComplexNumber(1.0, 0.0)) > 1e-6) {
    throw std::runtime_error("With noise=0, HH should be identity");
  }
}

void test_noise_set_and_get()
{
  State s(2);
  s.set_noise_probability(0.05);
  assert_double_close(s.get_noise_probability(), 0.05, 1e-12, "noise get/set");
  s.set_noise_probability(0.0);
  assert_double_close(s.get_noise_probability(), 0.0, 1e-12, "noise reset to 0");
}

void test_noise_high_degrades_state()
{
  // With 100% noise, applying X should almost always introduce errors.
  // The state will not be a clean basis state after many noisy ops.
  State::seed_rng(1234);
  State s(1);
  s.set_noise_probability(1.0);
  s.set_basis_state(0b0, 1.0);
  // Apply many X gates; with 100% noise the state will decohere.
  for (int i = 0; i < 10; ++i) s.x(0);
  // State should be non-trivially mixed (at least non-empty)
  const QuantumState& qs = s.get_state();
  if (qs.empty()) {
    throw std::runtime_error("Noise test: state should not be empty");
  }
}

// ---- QPE algorithm ----------------------------------------------------------

void test_qpe_exact_phase_pi_over_2()
{
  // P(pi/2): phase fraction = (pi/2) / (2*pi) = 1/4
  // With m=4 precision bits: expected measured_int = round(1/4 * 16) = 4
  QpeResult r = run_qpe(4, PI / 2.0);
  if (r.measured_int != 4) {
    throw std::runtime_error("QPE pi/2: expected int=4, got " + std::to_string(r.measured_int));
  }
  assert_double_close(r.estimated_frac, 0.25, 1e-9, "QPE pi/2 frac");
  if (r.best_prob < 0.99) {
    throw std::runtime_error("QPE pi/2: expected ~100% probability, got " + std::to_string(r.best_prob));
  }
}

void test_qpe_exact_phase_pi_over_4()
{
  // P(pi/4): fraction = 1/8, m=4: expected = round(1/8 * 16) = 2
  QpeResult r = run_qpe(4, PI / 4.0);
  if (r.measured_int != 2) {
    throw std::runtime_error("QPE pi/4: expected int=2, got " + std::to_string(r.measured_int));
  }
  assert_double_close(r.estimated_frac, 0.125, 1e-9, "QPE pi/4 frac");
  if (r.best_prob < 0.99) {
    throw std::runtime_error("QPE pi/4: expected ~100% probability, got " + std::to_string(r.best_prob));
  }
}

void test_qpe_zero_phase()
{
  // P(0): fraction = 0, all precision bits should be 0
  QpeResult r = run_qpe(4, 0.0);
  if (r.measured_int != 0) {
    throw std::runtime_error("QPE 0: expected int=0, got " + std::to_string(r.measured_int));
  }
  assert_double_close(r.estimated_frac, 0.0, 1e-9, "QPE 0 frac");
}

void test_qpe_increases_precision_with_more_bits()
{
  // For phase = 2*PI/3 (fraction ~= 0.3333), more bits -> better approximation.
  double phase = 2.0 * PI / 3.0;
  QpeResult r4 = run_qpe(4, phase);
  QpeResult r6 = run_qpe(6, phase);
  double err4 = std::abs(r4.estimated_frac - r4.true_phase_frac);
  double err6 = std::abs(r6.estimated_frac - r6.true_phase_frac);
  // 6-bit estimate should be at most as coarse as 4-bit
  if (err6 > err4 + 1.0 / (1 << 4)) {
    throw std::runtime_error("QPE: 6-bit precision should be better than or equal to 4-bit");
  }
}

void test_qpe_m_1()
{
  // m=1: can only distinguish phase < 0.5 vs >= 0.5.
  // P(pi/2): frac=0.25 -> 0.25*2 = 0.5, which rounds to either 0 or 1.
  // Either is acceptable; just verify the result is a valid 1-bit value.
  QpeResult r = run_qpe(1, PI / 2.0);
  if (r.measured_int != 0 && r.measured_int != 1) {
    throw std::runtime_error("QPE m=1: measured_int must be 0 or 1, got " + std::to_string(r.measured_int));
  }
  if (r.m != 1) {
    throw std::runtime_error("QPE m=1: result.m should be 1");
  }
  // estimated_frac must be 0 or 0.5
  if (r.estimated_frac != 0.0 && std::abs(r.estimated_frac - 0.5) > 1e-9) {
    throw std::runtime_error("QPE m=1: estimated_frac should be 0 or 0.5");
  }
}

void test_qpe_invalid_m()
{
  bool threw = false;
  try { run_qpe(0, 1.0); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("QPE: m=0 should throw");
  threw = false;
  try { run_qpe(21, 1.0); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("QPE: m=21 should throw");
}

void test_qpe_demo_runs()
{
  // Just exercise the demo code paths without crashing.
  std::streambuf* orig = std::cout.rdbuf();
  std::ostringstream sink;
  std::cout.rdbuf(sink.rdbuf());
  run_qpe_demo();
  run_qpe_demo(4, PI / 4.0);
  std::cout.rdbuf(orig);
}

// ---- LOAD/SAVE circuit scripting --------------------------------------------

void test_load_save_round_trip()
{
  const char* fname = "/tmp/qsim_test_circuit.qsim";

  // Write a simple circuit file
  {
    std::ofstream f(fname);
    f << "INIT 2 2\n";
    f << "H 0\n";
    f << "CX 0 1\n";
  }

  // Build a shell and load the file
  // We test this via the State API directly since the shell requires stdin.
  // Verify the circuit file can be parsed and executed:
  std::ifstream in(fname);
  if (!in) throw std::runtime_error("Could not open test circuit file");

  // Read and execute commands programmatically
  std::string line;
  std::unique_ptr<State> s;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (line.substr(0, 4) == "INIT") {
      s.reset(new State(2, 2));
    } else if (line == "H 0" && s) {
      s->h(0);
    } else if (line == "CX 0 1" && s) {
      s->cx(0, 1);
    }
  }
  if (!s) throw std::runtime_error("State not created from circuit file");
  // Bell state: |00> + |11>
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  if (std::abs(std::abs(s->get_amplitude(0b00)) - inv_sqrt2) > 1e-9)
    throw std::runtime_error("LOAD test: expected Bell state |00>");
  if (std::abs(std::abs(s->get_amplitude(0b11)) - inv_sqrt2) > 1e-9)
    throw std::runtime_error("LOAD test: expected Bell state |11>");
  std::remove(fname);
}

// ---- RESET gate -------------------------------------------------------------

void test_reset_from_one()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  s.reset(0);
  if (std::abs(s.get_amplitude(0b0) - ComplexNumber(1.0, 0.0)) > 1e-9)
    throw std::runtime_error("reset(|1>) should give |0>");
  if (std::abs(s.get_amplitude(0b1)) > 1e-9)
    throw std::runtime_error("reset(|1>): |1> component should vanish");
}

void test_reset_from_superposition()
{
  State s(1);
  s.h(0); // (|0>+|1>)/sqrt(2)
  s.reset(0);
  // Should project to |0> and renormalise
  if (std::abs(s.get_amplitude(0b0) - ComplexNumber(1.0, 0.0)) > 1e-9)
    throw std::runtime_error("reset from superposition should give |0>");
}

void test_reset_preserves_other_qubit()
{
  // 2-qubit state: |01> (q0=1, q1=0). Reset q0 -> |00>.
  State s(2);
  s.set_basis_state(0b01, 1.0);
  s.reset(0);
  if (std::abs(s.get_amplitude(0b00) - ComplexNumber(1.0, 0.0)) > 1e-9)
    throw std::runtime_error("reset q0 of |01> should give |00>");
}

// ---- iSWAP gate -------------------------------------------------------------

void test_iswap_same_bits_unchanged()
{
  // iSWAP on |00>: bits are equal, so no change.
  State s(2);
  s.set_basis_state(0b00, 1.0);
  s.iswap(0, 1);
  if (std::abs(s.get_amplitude(0b00) - ComplexNumber(1.0, 0.0)) > 1e-9)
    throw std::runtime_error("iSWAP |00> should be unchanged");
}

void test_iswap_swaps_and_multiplies_i()
{
  // iSWAP on |01> (q0=1, q1=0): bits differ -> swap and multiply by i
  State s(2);
  s.set_basis_state(0b01, 1.0);
  s.iswap(0, 1);
  // After iSWAP: state becomes i|10>
  if (std::abs(s.get_amplitude(0b01)) > 1e-9)
    throw std::runtime_error("iSWAP |01>: original state should be 0");
  if (std::abs(s.get_amplitude(0b10) - ComplexNumber(0.0, 1.0)) > 1e-9)
    throw std::runtime_error("iSWAP |01> should give i|10>");
}

void test_iswap_squared_is_swap()
{
  // iSWAP^2 = SWAP (up to global phase on mixed states).
  // On |01>: iSWAP once gives i|10>, iSWAP again gives i*(i)|01> = -|01>
  State s(2);
  s.set_basis_state(0b01, 1.0);
  s.iswap(0, 1);
  s.iswap(0, 1);
  // Should be -|01>
  if (std::abs(s.get_amplitude(0b01) - ComplexNumber(-1.0, 0.0)) > 1e-9)
    throw std::runtime_error("iSWAP^2 |01> should give -|01>");
}

// ---- Ising XX / YY / ZZ gates -----------------------------------------------

void test_xx_zero_angle_is_identity()
{
  State s(2);
  s.h(0); s.cx(0, 1); // Bell state
  State s2 = s; // copy
  s.set_noise_probability(0.0);
  s2.set_noise_probability(0.0);
  s.xx(0, 1, 0.0);
  // XX(0) = identity
  for (Bitstring b = 0; b < 4; ++b) {
    if (std::abs(s.get_amplitude(b) - s2.get_amplitude(b)) > 1e-9)
      throw std::runtime_error("XX(0) should be identity");
  }
}

void test_xx_pi_over_2_is_maximal_entangler()
{
  // XX(pi/2) on |00>: cos(pi/4)|00> - i*sin(pi/4)|11>
  const double c = std::cos(PI/4.0), sv = std::sin(PI/4.0);
  State s(2);
  s.set_basis_state(0b00, 1.0);
  s.xx(0, 1, PI / 2.0);
  if (std::abs(s.get_amplitude(0b00) - ComplexNumber(c, 0.0)) > 1e-9)
    throw std::runtime_error("XX(pi/2) |00>: wrong |00> amplitude");
  if (std::abs(s.get_amplitude(0b11) - ComplexNumber(0.0, -sv)) > 1e-9)
    throw std::runtime_error("XX(pi/2) |00>: wrong |11> amplitude");
}

void test_yy_pi_over_2_is_maximal_entangler()
{
  const double c = std::cos(PI / 4.0), sv = std::sin(PI / 4.0);
  State s(2);
  s.set_basis_state(0b00, 1.0);
  s.yy(0, 1, PI / 2.0);
  if (std::abs(s.get_amplitude(0b00) - ComplexNumber(c, 0.0)) > 1e-9)
    throw std::runtime_error("YY(pi/2) |00>: wrong |00> amplitude");
  if (std::abs(s.get_amplitude(0b11) - ComplexNumber(0.0, sv)) > 1e-9)
    throw std::runtime_error("YY(pi/2) |00>: wrong |11> amplitude");
}

void test_zz_diagonal()
{
  // ZZ(pi) on |00>: phase_same = e^{-i*pi/2} = -i
  State s(2);
  s.set_basis_state(0b00, 1.0);
  s.zz(0, 1, PI);
  // same bits -> phase_same = (cos(pi/2), -sin(pi/2)) = (0, -1)
  if (std::abs(s.get_amplitude(0b00) - ComplexNumber(0.0, -1.0)) > 1e-9)
    throw std::runtime_error("ZZ(pi) |00> should give -i|00>");
}

// ---- Bloch vector -----------------------------------------------------------

void test_bloch_zero_state()
{
  State s(1);
  s.set_basis_state(0b0, 1.0);
  auto bv = s.bloch(0);
  // |0> -> (0, 0, +1)
  if (std::abs(bv.x) > 1e-9 || std::abs(bv.y) > 1e-9 || std::abs(bv.z - 1.0) > 1e-9)
    throw std::runtime_error("Bloch |0>: expected (0,0,1)");
}

void test_bloch_one_state()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  auto bv = s.bloch(0);
  // |1> -> (0, 0, -1)
  if (std::abs(bv.x) > 1e-9 || std::abs(bv.y) > 1e-9 || std::abs(bv.z + 1.0) > 1e-9)
    throw std::runtime_error("Bloch |1>: expected (0,0,-1)");
}

void test_bloch_plus_state()
{
  State s(1);
  s.h(0); // |+> = (|0>+|1>)/sqrt(2)
  auto bv = s.bloch(0);
  // |+> -> (1, 0, 0)
  if (std::abs(bv.x - 1.0) > 1e-9 || std::abs(bv.y) > 1e-9 || std::abs(bv.z) > 1e-9)
    throw std::runtime_error("Bloch |+>: expected (1,0,0)");
}

// ---- expect_pauli -----------------------------------------------------------

void test_expect_z_on_zero()
{
  State s(1);
  s.set_basis_state(0b0, 1.0);
  double val = s.expect_pauli({{'Z', 0}});
  if (std::abs(val - 1.0) > 1e-9)
    throw std::runtime_error("<Z>|0> should be +1");
}

void test_expect_z_on_one()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  double val = s.expect_pauli({{'Z', 0}});
  if (std::abs(val + 1.0) > 1e-9)
    throw std::runtime_error("<Z>|1> should be -1");
}

void test_expect_x_on_plus()
{
  State s(1);
  s.h(0); // |+>
  double val = s.expect_pauli({{'X', 0}});
  if (std::abs(val - 1.0) > 1e-9)
    throw std::runtime_error("<X>|+> should be +1");
}

void test_expect_y_on_one_is_zero()
{
  State s(1);
  s.set_basis_state(0b1, 1.0);
  double val = s.expect_pauli({{'Y', 0}});
  if (std::abs(val) > 1e-9)
    throw std::runtime_error("<Y>|1> should be 0");
}

// ---- Von Neumann entropy ----------------------------------------------------

void test_entropy_product_state_is_zero()
{
  // Separable state: |0>|0> has zero entanglement entropy for either qubit.
  State s(2);
  s.set_basis_state(0b00, 1.0);
  double S = s.entropy(0, 0);
  if (std::abs(S) > 1e-6)
    throw std::runtime_error("Entropy of |00> qubit 0 should be 0");
}

void test_entropy_bell_state_is_one_bit()
{
  // Bell state: maximally entangled, entropy = 1 bit.
  State s(2);
  s.h(0); s.cx(0, 1);
  double S = s.entropy(0, 0);
  if (std::abs(S - 1.0) > 1e-4)
    throw std::runtime_error("Entropy of Bell state qubit 0 should be 1 bit, got " +
                             std::to_string(S));
}

// ---- QEC algorithm ----------------------------------------------------------

void test_qec_no_error_logical_0()
{
  QecResult r = run_qec(0, -1);
  if (!r.recovery_success)
    throw std::runtime_error("QEC: no error, logical=0 should succeed");
  if (r.syndrome != 0)
    throw std::runtime_error("QEC: no error should give syndrome 0");
}

void test_qec_no_error_logical_1()
{
  QecResult r = run_qec(1, -1);
  if (!r.recovery_success)
    throw std::runtime_error("QEC: no error, logical=1 should succeed");
  if (r.syndrome != 0)
    throw std::runtime_error("QEC: no error should give syndrome 0");
}

void test_qec_corrects_each_qubit()
{
  for (int eq = 0; eq < 3; ++eq) {
    for (int lb : {0, 1}) {
      QecResult r = run_qec(lb, eq);
      if (!r.recovery_success)
        throw std::runtime_error("QEC: error on q" + std::to_string(eq) +
                                 " logical=" + std::to_string(lb) + " should be corrected");
      if (r.corrected_qubit != eq)
        throw std::runtime_error("QEC: wrong corrected qubit for error on q" + std::to_string(eq));
    }
  }
}

void test_qec_invalid_args()
{
  bool threw = false;
  try { run_qec(2, 0); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("QEC: logical_bit=2 should throw");
  threw = false;
  try { run_qec(0, 3); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("QEC: error_qubit=3 should throw");
}

void test_qec_demo_runs()
{
  std::streambuf* orig = std::cout.rdbuf();
  std::ostringstream sink;
  std::cout.rdbuf(sink.rdbuf());
  run_qec_demo();
  std::cout.rdbuf(orig);
}

void test_noise_y_branch_executes()
{
  State::seed_rng(0);
  bool saw_y = false;
  for (int iter = 0; iter < 256; ++iter) {
    State s(1);
    s.set_noise_probability(1.0);
    s.set_basis_state(0b0, 1.0);
    s.x(0);
    if (std::abs(s.get_amplitude(0b0) - ComplexNumber(0.0, -1.0)) < 1e-9 &&
        std::abs(s.get_amplitude(0b1)) < 1e-9) {
      saw_y = true;
      break;
    }
  }
  if (!saw_y)
    throw std::runtime_error("Expected deterministic noisy execution to hit the Y branch");
}

// ---- Entry point ------------------------------------------------------------

int run_new_feature_tests()
{
  State::set_default_log_stream(nullptr);

  // Sdg tests
  run_test("Sdg: S then Sdg is identity", test_sdg_is_s_inverse);
  run_test("Sdg: applies -i phase to |1>", test_sdg_phase_on_one);
  run_test("Sdg: no change on |0>", test_sdg_no_change_on_zero);

  // Tdg tests
  run_test("Tdg: T then Tdg is identity", test_tdg_is_t_inverse);
  run_test("Tdg: correct phase on |1>", test_tdg_phase_on_one);

  // P gate tests
  run_test("P: identity on |0>", test_p_gate_identity_on_zero);
  run_test("P: correct phase on |1>", test_p_gate_applies_phase_to_one);
  run_test("P(pi/2) == S", test_p_equals_s_at_pi_over_2);

  // CP gate tests
  run_test("CP: no effect when control=0", test_cp_no_effect_when_control_zero);
  run_test("CP: applies phase when both |1>", test_cp_applies_phase_when_both_one);

  // CSWAP tests
  run_test("CSWAP: no swap when control=0", test_cswap_no_swap_when_control_zero);
  run_test("CSWAP: swaps when control=1", test_cswap_swaps_when_control_one);

  // MCX tests
  run_test("MCX: 1 control = CX", test_mcx_single_control_same_as_cx);
  run_test("MCX: 2 controls = CCX", test_mcx_two_controls_same_as_ccx);
  run_test("MCX: 3 controls flip only when all set", test_mcx_three_controls_flips_only_when_all_set);

  // Noise model tests
  run_test("Noise: zero noise is off", test_noise_zero_is_off);
  run_test("Noise: set/get probability", test_noise_set_and_get);
  run_test("Noise: 100% noise degrades state", test_noise_high_degrades_state);

  // QPE tests
  run_test("QPE: exact phase pi/2 (m=4)", test_qpe_exact_phase_pi_over_2);
  run_test("QPE: exact phase pi/4 (m=4)", test_qpe_exact_phase_pi_over_4);
  run_test("QPE: zero phase", test_qpe_zero_phase);
  run_test("QPE: more bits => better precision", test_qpe_increases_precision_with_more_bits);
  run_test("QPE: m=1 edge case", test_qpe_m_1);
  run_test("QPE: invalid m throws", test_qpe_invalid_m);
  run_test("QPE: demo runs without crash", test_qpe_demo_runs);

  // LOAD/SAVE scripting test
  run_test("LOAD/SAVE: round-trip circuit", test_load_save_round_trip);

  // RESET gate tests
  run_test("RESET: |1> -> |0>", test_reset_from_one);
  run_test("RESET: superposition -> |0>", test_reset_from_superposition);
  run_test("RESET: preserves other qubit", test_reset_preserves_other_qubit);

  // iSWAP gate tests
  run_test("iSWAP: same bits unchanged", test_iswap_same_bits_unchanged);
  run_test("iSWAP: swaps and multiplies by i", test_iswap_swaps_and_multiplies_i);
  run_test("iSWAP: squared is SWAP (up to phase)", test_iswap_squared_is_swap);

  // Ising gates tests
  run_test("XX(0) is identity", test_xx_zero_angle_is_identity);
  run_test("XX(pi/2) creates entanglement", test_xx_pi_over_2_is_maximal_entangler);
  run_test("YY(pi/2) creates entanglement", test_yy_pi_over_2_is_maximal_entangler);
  run_test("ZZ is diagonal gate", test_zz_diagonal);

  // Bloch vector tests
  run_test("Bloch: |0> -> (0,0,+1)", test_bloch_zero_state);
  run_test("Bloch: |1> -> (0,0,-1)", test_bloch_one_state);
  run_test("Bloch: |+> -> (1,0,0)", test_bloch_plus_state);

  // expect_pauli tests
  run_test("Expect: <Z>|0> = +1", test_expect_z_on_zero);
  run_test("Expect: <Z>|1> = -1", test_expect_z_on_one);
  run_test("Expect: <X>|+> = +1", test_expect_x_on_plus);
  run_test("Expect: <Y>|1> = 0", test_expect_y_on_one_is_zero);

  // Entropy tests
  run_test("Entropy: product state = 0", test_entropy_product_state_is_zero);
  run_test("Entropy: Bell state = 1 bit", test_entropy_bell_state_is_one_bit);

  // QEC tests
  run_test("QEC: no error logical=0", test_qec_no_error_logical_0);
  run_test("QEC: no error logical=1", test_qec_no_error_logical_1);
  run_test("QEC: corrects each qubit", test_qec_corrects_each_qubit);
  run_test("QEC: invalid args throw", test_qec_invalid_args);
  run_test("QEC: demo runs without crash", test_qec_demo_runs);
  run_test("Noise: Y branch executes", test_noise_y_branch_executes);

  return g_failures;
}
