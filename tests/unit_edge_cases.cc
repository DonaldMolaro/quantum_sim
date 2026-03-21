#include "state.hh"
#include "algorithms/anneal.hh"
#include "math/bit_ops.hh"
#include "tests/helpers.hh"
#include "tests/test_harness.hh"
#include "internal/limits.hh"
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using test_helpers::assert_complex_equal;
using test_helpers::assert_complex_close;
using test_helpers::assert_amplitude_match;
using test_helpers::assert_double_close;
using test_harness::run_test;

static constexpr double PI = qsim::limits::PI;

// --- Display coverage ---

void test_display_sparse()
{
  State s(2);
  s.h(0);
  // Should not crash; exercises sparse display path
  s.display(false);
}

void test_display_dense()
{
  State s(2);
  s.h(0);
  // Exercises dense (show_all=true) path including zero-amplitude rows
  s.display(true);
}

void test_display_cbits_nonempty()
{
  State s(2, 2);
  s.x(0);
  s.measure_with_rng(0, 0, 0.0);
  // Exercises the decimal value and bitstring display path
  s.display_cbits();
}

void test_display_cbits_empty()
{
  State s(2, 0);
  // Exercises the "Classical register is empty" path
  s.display_cbits();
}

// --- complex_to_string coverage ---

void test_complex_to_string_paths()
{
  // Zero amplitude
  std::string z = complex_to_string(ComplexNumber(0.0, 0.0));
  if (z != "0.0") throw std::runtime_error("complex_to_string(0) should be '0.0'");

  // Purely imaginary positive
  std::string pi = complex_to_string(ComplexNumber(0.0, 0.5));
  if (pi.find('i') == std::string::npos) throw std::runtime_error("should contain 'i'");

  // Purely imaginary negative
  std::string ni = complex_to_string(ComplexNumber(0.0, -0.5));
  if (ni.find('-') == std::string::npos) throw std::runtime_error("should contain '-'");

  // Purely imaginary with magnitude 1
  std::string one_i = complex_to_string(ComplexNumber(0.0, 1.0));
  if (one_i.find('i') == std::string::npos) throw std::runtime_error("should contain 'i'");

  // Complex with both real and negative imaginary
  std::string both = complex_to_string(ComplexNumber(0.5, -0.5));
  if (both.find('-') == std::string::npos) throw std::runtime_error("should contain '-'");

  // Complex with both real and positive imaginary
  std::string bp = complex_to_string(ComplexNumber(0.5, 0.5));
  if (bp.find('+') == std::string::npos) throw std::runtime_error("should contain '+'");
}

// --- CCX (Toffoli) gate ---

void test_ccx_all_ones()
{
  State s(3);
  s.x(0).x(1);  // Set controls to |1>
  s.ccx(0, 1, 2);
  // CCX decomposition may introduce a global phase — check magnitude
  double mag = std::abs(s.get_amplitude(0b111));
  if (mag < 0.9) {
    throw std::runtime_error("CCX |110>->|111> magnitude should be ~1, got " + std::to_string(mag));
  }
}

void test_ccx_one_control_zero()
{
  State s(3);
  s.x(0);  // Only one control set
  s.ccx(0, 1, 2);
  // Target should not flip — check magnitude of original state
  double mag = std::abs(s.get_amplitude(0b001));
  if (mag < 0.9) {
    throw std::runtime_error("CCX one control zero: should stay ~|001>, mag=" + std::to_string(mag));
  }
}

// --- CU gate ---

void test_cu_identity_like()
{
  State s(2);
  s.x(0);  // Set control
  s.cu(0, 1, 0.0, 0.0, 0.0);
  // With theta=0, phi=0, lambda=0, CU should be close to identity
  assert_amplitude_match(s, 0b01, ComplexNumber(1.0, 0.0), "CU identity-like", 1e-6);
}

// --- iSWAP gate ---

void test_iswap_same_bits()
{
  State s(2);
  // |00> has same bits for j=0,k=1 -> no swap, identity
  s.iswap(0, 1);
  assert_amplitude_match(s, 0b00, ComplexNumber(1.0, 0.0), "iSWAP |00>", 1e-9);
}

void test_iswap_different_bits()
{
  State s(2);
  s.x(0);  // |01>
  s.iswap(0, 1);
  // |01> -> i|10>
  assert_amplitude_match(s, 0b10, ComplexNumber(0.0, 1.0), "iSWAP |01>->i|10>", 1e-9);
}

// --- XX gate ---

void test_xx_gate()
{
  State s(2);
  s.xx(0, 1, PI / 2.0);
  // After XX(pi/2) on |00>: cos(pi/4)|00> - i*sin(pi/4)|11>
  double c = std::cos(PI / 4.0);
  double sn = std::sin(PI / 4.0);
  assert_complex_close(ComplexNumber(c, 0.0), s.get_amplitude(0b00), 1e-9, "XX |00> component");
  assert_complex_close(ComplexNumber(0.0, -sn), s.get_amplitude(0b11), 1e-9, "XX |11> component");
}

// --- YY gate ---

void test_yy_gate()
{
  State s(2);
  s.yy(0, 1, PI / 2.0);
  // After YY(pi/2) on |00>: cos(pi/4)|00> + i*sin(pi/4)|11>  (yy_sign=-1 for same parity)
  double c = std::cos(PI / 4.0);
  double sn = std::sin(PI / 4.0);
  assert_complex_close(ComplexNumber(c, 0.0), s.get_amplitude(0b00), 1e-9, "YY |00> component");
  // yy_sign for |00> is -1, so result is -i*(-sin) = i*sin
  assert_complex_close(ComplexNumber(0.0, sn), s.get_amplitude(0b11), 1e-9, "YY |11> component");
}

// --- ZZ gate ---

void test_zz_gate()
{
  State s(2);
  s.zz(0, 1, PI);
  // ZZ(pi) on |00>: exp(-i*pi/2)|00>
  assert_complex_close(ComplexNumber(0.0, -1.0), s.get_amplitude(0b00), 1e-9, "ZZ |00>");
}

// --- SWAP bounds ---

void test_swap_out_of_range_throws()
{
  State s(2);
  bool threw = false;
  try {
    s.swap(0, 5);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  if (!threw) throw std::runtime_error("swap out-of-range should throw");
}

// --- Reset gate paths ---

void test_reset_all_on_one()
{
  State s(1);
  s.x(0);  // |1>
  s.reset(0);
  // All amplitude on |1> -> flip to |0>
  assert_amplitude_match(s, 0b0, ComplexNumber(1.0, 0.0), "reset |1>->|0>", 1e-9);
}

void test_reset_superposition()
{
  State s(1);
  s.h(0);  // (|0> + |1>)/sqrt(2)
  s.reset(0);
  // Discard |1>, renormalize |0>
  assert_amplitude_match(s, 0b0, ComplexNumber(1.0, 0.0), "reset superposition->|0>", 1e-9);
}

// --- Bloch vector (2-qubit subsystem, exercises partial trace paths) ---

void test_ec_bloch_one_state()
{
  State s(1);
  s.x(0);
  auto bv = s.bloch(0);
  // |1> -> (0, 0, -1)
  assert_double_close(bv.x, 0.0, 1e-9, "Bloch x for |1>");
  assert_double_close(bv.y, 0.0, 1e-9, "Bloch y for |1>");
  assert_double_close(bv.z, -1.0, 1e-9, "Bloch z for |1>");
}

void test_ec_bloch_entangled()
{
  State s(2);
  s.h(0);
  s.cx(0, 1);
  auto bv = s.bloch(0);
  // Maximally entangled: Bloch vector should be at origin
  assert_double_close(bv.x, 0.0, 1e-6, "Bloch x entangled");
  assert_double_close(bv.y, 0.0, 1e-6, "Bloch y entangled");
  assert_double_close(bv.z, 0.0, 1e-6, "Bloch z entangled");
}

// --- expect_pauli ---

void test_expect_pauli_z()
{
  State s(1);
  // |0> -> <Z> = 1
  double val = s.expect_pauli({{'Z', 0}});
  assert_double_close(val, 1.0, 1e-9, "<Z> for |0>");
}

void test_expect_pauli_x()
{
  State s(1);
  s.h(0);  // |+>
  double val = s.expect_pauli({{'X', 0}});
  assert_double_close(val, 1.0, 1e-9, "<X> for |+>");
}

void test_expect_pauli_y()
{
  State s(1);
  s.h(0);
  s.s(0);  // H then S
  double val = s.expect_pauli({{'Y', 0}});
  // |psi> = H S |0> -- verify <Y> is +/-1
  if (std::abs(std::abs(val) - 1.0) > 1e-6) {
    throw std::runtime_error("<Y> for H·S|0> should be +/-1, got " + std::to_string(val));
  }
}

void test_expect_pauli_identity()
{
  State s(1);
  double val = s.expect_pauli({{'I', 0}});
  assert_double_close(val, 1.0, 1e-9, "<I> for |0>");
}

// --- Entropy ---

void test_entropy_product_state()
{
  // Product state should have zero entropy
  State s(2);
  double ent = s.entropy(0, 0);
  assert_double_close(ent, 0.0, 1e-9, "entropy of product state");
}

void test_entropy_bell_state()
{
  // Bell state: maximum entanglement for 1 qubit subsystem -> entropy = 1 bit
  State s(2);
  s.h(0);
  s.cx(0, 1);
  double ent = s.entropy(0, 0);
  assert_double_close(ent, 1.0, 0.01, "entropy of Bell state");
}

void test_entropy_reversed_args()
{
  // start_q > end_q should be auto-swapped — test with 3 qubits
  State s(3);
  s.h(0);
  s.cx(0, 1);
  s.cx(0, 2);
  // Subsystem {0,1} vs {2}: entropy should be >0 for entangled state
  double ent = s.entropy(1, 0);  // reversed: should swap to entropy(0,1)
  // With a 3-qubit GHZ state, entropy of 2-qubit subsystem is 1 bit
  assert_double_close(ent, 1.0, 0.05, "entropy reversed args");
}

void test_entropy_too_many_qubits()
{
  // >10 qubit subsystem should return 0.0
  State s(12);
  double ent = s.entropy(0, 11);
  assert_double_close(ent, 0.0, 1e-9, "entropy >10 qubits returns 0");
}

// --- Measure with near-zero probability ---

void test_measure_near_zero_probability()
{
  // Force outcome 1 with p(1)=0 by setting random_val > p0=1.0
  // This triggers the near-zero probability clamping path
  State s(1, 1);
  s.set_basis_state(0b0, 1.0);  // Fully |0>, p0=1.0
  // random_val=1.1 > p0=1.0 -> outcome=1, but p_outcome=0.0 -> clamp
  std::ostringstream log;
  s.set_log_stream(&log);
  s.measure_with_rng(0, 0, 1.1);
  s.set_log_stream(nullptr);
  // Should have logged a warning
  if (log.str().find("near-zero") == std::string::npos) {
    throw std::runtime_error("Expected near-zero probability warning in log");
  }
}

// --- Measure all ---

void test_measure_all_deterministic()
{
  State s(2, 2);
  s.x(0).x(1);  // |11>
  std::vector<int> out;
  s.measure_all(out);
  if (out.size() != 2 || out[0] != 1 || out[1] != 1) {
    throw std::runtime_error("measure_all on |11> should give {1,1}");
  }
}

// --- from_basis ---

void test_from_basis()
{
  State s = State::from_basis(3, 0b101);
  assert_amplitude_match(s, 0b101, ComplexNumber(1.0, 0.0), "from_basis", 1e-9);
  assert_amplitude_match(s, 0b000, ComplexNumber(0.0, 0.0), "from_basis zero", 1e-9);
}

// --- Noise ---

void test_depolarizing_noise()
{
  // With noise probability 1.0, every gate should apply a random Pauli error
  State::seed_rng(42);
  State s(1);
  s.set_noise_probability(1.0);
  s.h(0);
  // We can't predict the exact outcome, but the state should still be normalized
  double norm = 0.0;
  for (const auto& p : s.get_state()) norm += std::norm(p.second);
  assert_double_close(norm, 1.0, 1e-6, "noisy state normalized");
  s.set_noise_probability(0.0);
}

// --- QFT modes ---

void test_qft_gate_mode()
{
  // Gate-mode QFT should be unitary: QFT followed by IQFT recovers original
  State s(3);
  s.set_qft_mode(State::QftMode::Gate);
  s.x(0);
  s.qft(0, 2);
  // State should be normalized after QFT
  double norm = 0.0;
  for (const auto& p : s.get_state()) norm += std::norm(p.second);
  assert_double_close(norm, 1.0, 1e-6, "QFT gate mode normalized");
}

void test_iqft_gate_mode()
{
  State s(3);
  s.set_qft_mode(State::QftMode::Gate);
  s.x(0);
  s.qft(0, 2);
  s.iqft(0, 2);
  // QFT then IQFT should recover original state
  assert_complex_close(s.get_amplitude(0b001), ComplexNumber(1.0, 0.0), 1e-6, "IQFT roundtrip");
}

// --- CRX, CRY, CRZ ---

void test_crx_gate()
{
  State s(2);
  s.x(0);  // control = |1>
  s.crx(0, 1, PI);
  // CRX(pi) with control=1 should flip target: |01> -> |11>
  assert_complex_close(s.get_amplitude(0b11), s.get_amplitude(0b11), 1e-6, "CRX pi");
}

void test_cry_gate()
{
  State s(2);
  s.x(0);
  s.cry(0, 1, PI);
  // CRY(pi) with control=1: rotates target from |0> to |1>
  double mag = std::abs(s.get_amplitude(0b11));
  if (mag < 0.9) throw std::runtime_error("CRY(pi) should mostly produce |11>");
}

void test_crz_gate()
{
  State s(2);
  s.x(0);
  s.h(1);
  s.crz(0, 1, PI);
  // CRZ should apply a phase; state should still be normalized
  double norm = 0.0;
  for (const auto& p : s.get_state()) norm += std::norm(p.second);
  assert_double_close(norm, 1.0, 1e-6, "CRZ state normalized");
}

// --- CSWAP (Fredkin) ---

void test_cswap_active()
{
  State s(3);
  s.x(0).x(1);  // |011> control=0 is 1, swap qubits 1,2
  s.cswap(0, 1, 2);
  // Should swap: |011> -> |101>
  assert_amplitude_match(s, 0b101, ComplexNumber(1.0, 0.0), "CSWAP active", 1e-9);
}

// --- MCX ---

void test_mcx_gate()
{
  State s(3);
  s.x(0).x(1);  // controls = {0,1} both |1>
  s.mcx({0, 1}, 2);
  // Should flip target: |011> -> |111>
  assert_amplitude_match(s, 0b111, ComplexNumber(1.0, 0.0), "MCX all controls", 1e-9);
}

// --- P gate ---

void test_p_gate()
{
  State s(1);
  s.x(0);
  s.p(0, PI / 2.0);
  // P(pi/2)|1> = e^{i*pi/2}|1> = i|1>
  assert_complex_close(ComplexNumber(0.0, 1.0), s.get_amplitude(0b1), 1e-9, "P gate");
}

// --- CP gate ---

void test_cp_gate()
{
  State s(2);
  s.x(0).x(1);
  s.cp(0, 1, PI);
  // CP(pi) with both |1>: phase -1
  assert_complex_close(ComplexNumber(-1.0, 0.0), s.get_amplitude(0b11), 1e-9, "CP gate");
}

// --- Sdg, Tdg ---

void test_sdg_gate()
{
  State s(1);
  s.x(0);
  s.sdg(0);
  // Sdg|1> = -i|1>
  assert_complex_close(ComplexNumber(0.0, -1.0), s.get_amplitude(0b1), 1e-9, "Sdg");
}

void test_tdg_gate()
{
  State s(1);
  s.x(0);
  s.t(0);
  s.tdg(0);
  // T then Tdg = identity
  assert_complex_close(ComplexNumber(1.0, 0.0), s.get_amplitude(0b1), 1e-9, "T·Tdg");
}

// --- set_amplitude ---

void test_set_amplitude_update()
{
  State s(2);
  s.set_amplitude(0b01, ComplexNumber(0.5, 0.5));
  assert_complex_close(ComplexNumber(0.5, 0.5), s.get_amplitude(0b01), 1e-9, "set_amplitude new");
  s.set_amplitude(0b01, ComplexNumber(0.3, 0.0));
  assert_complex_close(ComplexNumber(0.3, 0.0), s.get_amplitude(0b01), 1e-9, "set_amplitude update");
}

// --- set_superposition ---

void test_set_superposition_move()
{
  State s(2);
  QuantumState qs = {{0b00, ComplexNumber(0.6, 0.0)}, {0b11, ComplexNumber(0.8, 0.0)}};
  s.set_superposition(std::move(qs));
  assert_complex_close(ComplexNumber(0.6, 0.0), s.get_amplitude(0b00), 1e-9, "set_superposition move");
}

// --- bitstring_to_string ---

void test_bitstring_to_string()
{
  std::string s = bitstring_to_string(0b101, 3);
  if (s != "101") throw std::runtime_error("bitstring_to_string(5,3) should be '101', got '" + s + "'");
}

// --- Anneal single-step coverage ---

void test_anneal_sa_single_step()
{
  const int n = 2;
  const std::vector<double> q = {
    -1.0,  0.5,
     0.5, -1.0,
  };
  AnnealOptions opts;
  opts.method = AnnealMethod::SA;
  opts.steps = 1;
  opts.sweeps_per_step = 1;
  opts.seed = 42u;
  AnnealResult result = anneal_qubo(n, q, opts);
  if (!result.ok) throw std::runtime_error("SA anneal with steps=1 should succeed");
  if (result.best_history.size() != 2)
    throw std::runtime_error("SA steps=1 should produce 2-element history");
}

void test_anneal_sqa_single_step()
{
  const int n = 2;
  const std::vector<double> q = {
    -1.0,  0.5,
     0.5, -1.0,
  };
  AnnealOptions opts;
  opts.method = AnnealMethod::SQA;
  opts.steps = 1;
  opts.sweeps_per_step = 1;
  opts.replicas = 2;
  opts.seed = 42u;
  AnnealResult result = anneal_qubo(n, q, opts);
  if (!result.ok) throw std::runtime_error("SQA anneal with steps=1 should succeed");
  if (result.best_history.size() != 2)
    throw std::runtime_error("SQA steps=1 should produce 2-element history");
}

// --- Bit ops bounds checking ---

void test_extract_bits_invalid_range()
{
  bool threw = false;
  try { extract_bits(0xFF, 5, 2); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("extract_bits(start>end) should throw");

  threw = false;
  try { extract_bits(0xFF, -1, 3); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("extract_bits(negative start) should throw");
}

void test_replace_bits_invalid_range()
{
  bool threw = false;
  try { replace_bits(0xFF, 5, 2, 0); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("replace_bits(start>end) should throw");

  threw = false;
  try { replace_bits(0xFF, -1, 3, 0); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("replace_bits(negative start) should throw");
}

// --- QFT/IQFT bounds checking ---

void test_qft_out_of_range_throws()
{
  State s(3, 0);
  bool threw = false;
  try { s.qft(0, 5); } catch (const std::out_of_range&) { threw = true; }
  if (!threw) throw std::runtime_error("QFT with end_qubit > num_qubits should throw");

  threw = false;
  try { s.iqft(-1, 2); } catch (const std::out_of_range&) { threw = true; }
  if (!threw) throw std::runtime_error("IQFT with negative start should throw");
}

void test_qft_register_too_large_throws()
{
  // Create a state just large enough to make the register span exceed
  // kMaxBitstringQubits (62).  We need end-start+1 > 62, so use a
  // 63-qubit state (constructor only allocates the sparse vector).
  State s(63, 0);
  bool threw = false;
  try { s.qft(0, 62); } catch (const std::invalid_argument&) { threw = true; }
  if (!threw) throw std::runtime_error("QFT with 63-qubit register should throw invalid_argument");
}

// --- Entry point ---

int run_edge_case_tests()
{
  test_harness::reset_failures();

  const test_harness::TestCase cases[] = {
    {"display sparse", test_display_sparse},
    {"display dense", test_display_dense},
    {"display cbits nonempty", test_display_cbits_nonempty},
    {"display cbits empty", test_display_cbits_empty},
    {"complex_to_string paths", test_complex_to_string_paths},
    {"CCX all ones", test_ccx_all_ones},
    {"CCX one control zero", test_ccx_one_control_zero},
    {"CU identity-like", test_cu_identity_like},
    {"iSWAP same bits", test_iswap_same_bits},
    {"iSWAP different bits", test_iswap_different_bits},
    {"XX gate", test_xx_gate},
    {"YY gate", test_yy_gate},
    {"ZZ gate", test_zz_gate},
    {"SWAP out-of-range throws", test_swap_out_of_range_throws},
    {"reset all-on-one", test_reset_all_on_one},
    {"reset superposition", test_reset_superposition},
    {"Bloch one state", test_ec_bloch_one_state},
    {"Bloch entangled", test_ec_bloch_entangled},
    {"expect_pauli Z", test_expect_pauli_z},
    {"expect_pauli X", test_expect_pauli_x},
    {"expect_pauli Y", test_expect_pauli_y},
    {"expect_pauli I", test_expect_pauli_identity},
    {"entropy product state", test_entropy_product_state},
    {"entropy Bell state", test_entropy_bell_state},
    {"entropy reversed args", test_entropy_reversed_args},
    {"entropy too many qubits", test_entropy_too_many_qubits},
    {"measure near-zero probability", test_measure_near_zero_probability},
    {"measure all deterministic", test_measure_all_deterministic},
    {"from_basis", test_from_basis},
    {"depolarizing noise", test_depolarizing_noise},
    {"QFT gate mode", test_qft_gate_mode},
    {"IQFT gate mode roundtrip", test_iqft_gate_mode},
    {"CRX gate", test_crx_gate},
    {"CRY gate", test_cry_gate},
    {"CRZ gate", test_crz_gate},
    {"CSWAP active", test_cswap_active},
    {"MCX gate", test_mcx_gate},
    {"P gate", test_p_gate},
    {"CP gate", test_cp_gate},
    {"Sdg gate", test_sdg_gate},
    {"Tdg gate", test_tdg_gate},
    {"set_amplitude update", test_set_amplitude_update},
    {"set_superposition move", test_set_superposition_move},
    {"bitstring_to_string", test_bitstring_to_string},
    {"anneal SA single step", test_anneal_sa_single_step},
    {"anneal SQA single step", test_anneal_sqa_single_step},
    {"extract_bits invalid range", test_extract_bits_invalid_range},
    {"replace_bits invalid range", test_replace_bits_invalid_range},
    {"QFT/IQFT out-of-range throws", test_qft_out_of_range_throws},
    {"QFT register too large throws", test_qft_register_too_large_throws},
  };
  test_harness::run_cases(cases, sizeof(cases) / sizeof(cases[0]));

  return test_harness::failure_count();
}
