#include "state.hh"
#include "tests/helpers.hh"
#include "tests/unit_algorithms_and_cli.hh"
#include "tests/unit_qft_modexp.hh"
#include "algorithms/qrng.hh"
#include "algorithms/latin_square.hh"
#include "algorithms/deutsch_jozsa.hh"
#include "algorithms/bernstein_vazirani.hh"
#include "algorithms/qubo.hh"
#include "algorithms/vqa_qaoa.hh"
#include "algorithms/qaoa.hh"
#include "algorithms/vqe.hh"
#include "algorithms/anneal.hh"
#include "algorithms/api/grover_api.hh"
#include "algorithms/api/shor_api.hh"
#include "algorithms/shor_classical.hh"
#include "algorithms/shor_quantum.hh"
#include "cli/commands.hh"
#include "logging.hh"
#include "demos/grover_demo.hh"
#include "demos/latin_demo.hh"
#include "demos/qaoa_demo.hh"
#include "demos/shor_demo.hh"
#include "math/mod_arith.hh"
#include "math/bit_ops.hh"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <sys/wait.h>
#include <unistd.h>

const ComplexNumber I(0.0, 1.0); // Imaginary unit i

// --- Utility Functions ---

static int g_failures = 0;

void run_test(const std::string& name, std::function<void()> test_func) {
    try {
        test_func();
        std::cout << "[PASS] " << name << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] " << name << ": " << e.what() << std::endl;
        ++g_failures;
    } catch (...) {
        std::cerr << "[FAIL] " << name << ": Unknown error" << std::endl;
        ++g_failures;
    }
}

using test_helpers::assert_complex_equal;
using test_helpers::assert_complex_close;
using test_helpers::assert_equal;
using test_helpers::assert_amplitude_match;
using test_helpers::assert_amplitude_magnitude;

// --- Core Gate Tests ---

void test_x_gate_flips_basis() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    s.x(0);
    assert_complex_equal(0.0, s.get_amplitude(0b0), "X should clear |0> amplitude");
    assert_complex_equal(1.0, s.get_amplitude(0b1), "X should set |1> amplitude");
}

void test_h_gate_on_zero() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    s.h(0);
    ComplexNumber coeff = 1.0 / std::sqrt(2.0);
    assert_complex_equal(coeff, s.get_amplitude(0b0), "H|0> amplitude for |0>");
    assert_complex_equal(coeff, s.get_amplitude(0b1), "H|0> amplitude for |1>");
}

void test_h_gate_on_one() {
    State s(1, 0);
    s.set_basis_state(0b1, 1.0);
    s.h(0);
    ComplexNumber coeff = 1.0 / std::sqrt(2.0);
    assert_complex_equal(coeff, s.get_amplitude(0b0), "H|1> amplitude for |0>");
    assert_complex_equal(-coeff, s.get_amplitude(0b1), "H|1> amplitude for |1>");
}

void test_s_gate_phase() {
    State s(1, 0);
    s.set_basis_state(0b1, 1.0);
    s.s(0);
    assert_complex_equal(I, s.get_amplitude(0b1), "S should apply i phase to |1>");
}

void test_y_gate_on_zero() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    s.y(0);
    ComplexNumber a1 = s.get_amplitude(0b1);
    bool ok = (std::abs(a1 - ComplexNumber(0.0, 1.0)) < 1e-9) ||
              (std::abs(a1 - ComplexNumber(0.0, -1.0)) < 1e-9);
    if (!ok) {
        throw std::runtime_error("Y should map |0> to ±i|1> (global phase)");
    }
}

void test_z_gate_phase() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    s.h(0); // (|0> + |1>) / sqrt(2)
    s.z(0); // should flip phase of |1> relative to |0>

    ComplexNumber a0 = s.get_amplitude(0b0);
    ComplexNumber a1 = s.get_amplitude(0b1);
    assert_complex_close(a0, -a1, 1e-9, "Z should flip the relative phase between |0> and |1>");
}

void test_t_gate_phase() {
    State s(1, 0);
    s.set_basis_state(0b1, 1.0);
    s.t(0);
    ComplexNumber expected(ONE_OVER_SQRT_TWO, ONE_OVER_SQRT_TWO);
    assert_complex_equal(expected, s.get_amplitude(0b1), "T should apply (1+i)/sqrt(2) phase to |1>");
}

void test_rx_pi_on_zero() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    const double pi = std::acos(-1.0);
    s.rx(0, pi);
    assert_complex_close(0.0, s.get_amplitude(0b0), 1e-9, "RX(pi) should clear |0>");
    assert_complex_close(ComplexNumber(0.0, -1.0), s.get_amplitude(0b1), 1e-9, "RX(pi) should map |0> to -i|1>");
}

void test_ry_pi_on_zero() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    const double pi = std::acos(-1.0);
    s.ry(0, pi);
    assert_complex_close(0.0, s.get_amplitude(0b0), 1e-9, "RY(pi) should clear |0>");
    assert_complex_close(1.0, s.get_amplitude(0b1), 1e-9, "RY(pi) should map |0> to |1>");
}

void test_rz_pi_phases() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    const double pi = std::acos(-1.0);
    s.rz(0, pi);
    assert_complex_close(ComplexNumber(0.0, -1.0), s.get_amplitude(0b0), 1e-9, "RZ(pi) should map |0> to -i|0>");

    State s1(1, 0);
    s1.set_basis_state(0b1, 1.0);
    s1.rz(0, pi);
    assert_complex_close(ComplexNumber(0.0, 1.0), s1.get_amplitude(0b1), 1e-9, "RZ(pi) should map |1> to i|1>");
}

void test_qrng_deterministic_bits() {
    std::vector<double> random_vals = {0.1, 0.9, 0.1, 0.9};
    std::vector<int> bits = qrng_bits(4, &random_vals);
    assert_equal(bits.size(), static_cast<size_t>(4), "QRNG should return 4 bits");
    assert_equal(bits[0], 0, "QRNG bit0 should be 0");
    assert_equal(bits[1], 1, "QRNG bit1 should be 1");
    assert_equal(bits[2], 0, "QRNG bit2 should be 0");
    assert_equal(bits[3], 1, "QRNG bit3 should be 1");

    uint64_t value = qrng_u64(4, &random_vals);
    assert_equal(value, static_cast<uint64_t>(10), "QRNG value should match bits (1010b)");
}

void test_ru_matches_ry_pi() {
    State s(1, 0);
    s.set_basis_state(0b0, 1.0);
    const double pi = std::acos(-1.0);
    s.ru(0, pi, 0.0, 0.0); // U(pi,0,0) == RY(pi)
    assert_complex_close(0.0, s.get_amplitude(0b0), 1e-9, "RU(pi,0,0) should clear |0>");
    assert_complex_close(1.0, s.get_amplitude(0b1), 1e-9, "RU(pi,0,0) should map |0> to |1>");
}

void test_ru_zero_theta_phase() {
    const double pi = std::acos(-1.0);
    const double phi = pi / 3.0;
    const double lambda = -pi / 5.0;

    State s0(1, 0);
    s0.set_basis_state(0b0, 1.0);
    s0.ru(0, 0.0, phi, lambda);
    assert_complex_close(1.0, s0.get_amplitude(0b0), 1e-9, "RU(0,phi,lambda) should keep |0> amplitude");

    State s1(1, 0);
    s1.set_basis_state(0b1, 1.0);
    s1.ru(0, 0.0, phi, lambda);
    ComplexNumber expected(std::cos(phi + lambda), std::sin(phi + lambda));
    assert_complex_close(expected, s1.get_amplitude(0b1), 1e-9, "RU(0,phi,lambda) should apply e^{i(phi+lambda)} to |1>");
}

void test_cx_gate_control_on() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0); // control=1, target=0
    s.cx(1, 0);
    assert_complex_equal(0.0, s.get_amplitude(0b10), "CX should clear |10>");
    assert_complex_equal(1.0, s.get_amplitude(0b11), "CX should flip target when control=1");
}

void test_cx_gate_control_off() {
    State s(2, 0);
    s.set_basis_state(0b00, 1.0); // control=0
    s.cx(1, 0);
    assert_complex_equal(1.0, s.get_amplitude(0b00), "CX should preserve |00>");
    assert_complex_equal(0.0, s.get_amplitude(0b01), "CX should not flip target when control=0");
}

void test_cnot_alias_matches_cx() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0);
    s.cnot(1, 0);
    assert_complex_equal(0.0, s.get_amplitude(0b10), "CNOT should clear |10>");
    assert_complex_equal(1.0, s.get_amplitude(0b11), "CNOT should map |10> to |11>");
}

void test_cz_gate_phase_on_11() {
    State s(2, 0);
    s.set_basis_state(0b11, 1.0);
    s.cz(1, 0);
    assert_complex_equal(-1.0, s.get_amplitude(0b11), "CZ should apply -1 phase to |11>");
}

void test_cy_gate_on_10() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0); // control=1, target=0
    s.cy(1, 0);
    assert_complex_close(ComplexNumber(0.0, 1.0), s.get_amplitude(0b11), 1e-9, "CY should map |10> to i|11>");
}

void test_ch_gate_on_10() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0); // control=1, target=0
    s.ch(1, 0);
    ComplexNumber a10 = s.get_amplitude(0b10);
    ComplexNumber a11 = s.get_amplitude(0b11);
    double expected_mag = 1.0 / std::sqrt(2.0);
    assert_complex_close(expected_mag, std::abs(a10), 1e-9, "CH should give |10> amplitude magnitude 1/sqrt(2)");
    assert_complex_close(expected_mag, std::abs(a11), 1e-9, "CH should give |11> amplitude magnitude 1/sqrt(2)");
    assert_complex_close(a10, a11, 1e-9, "CH should keep |10> and |11> in-phase (up to global phase)");
}

void test_crz_gate_phase_on_control_on() {
    const double pi = std::acos(-1.0);
    const double theta = pi / 2.0;

    State s(2, 0);
    s.set_basis_state(0b11, 1.0); // control=1, target=1
    s.crz(1, 0, theta);
    ComplexNumber expected(std::cos(theta / 2.0), std::sin(theta / 2.0));
    assert_complex_close(expected, s.get_amplitude(0b11), 1e-9, "CRZ should apply phase e^{i theta/2} to |11>");
}

void test_crx_gate_on_control_on() {
    const double pi = std::acos(-1.0);
    State s(2, 0);
    s.set_basis_state(0b10, 1.0); // control=1, target=0
    s.crx(1, 0, pi);
    ComplexNumber a10 = s.get_amplitude(0b10);
    ComplexNumber a11 = s.get_amplitude(0b11);
    assert_complex_close(0.0, a10, 1e-9, "CRX(pi) should clear |10>");
    bool ok = (std::abs(a11 - ComplexNumber(0.0, 1.0)) < 1e-9) ||
              (std::abs(a11 - ComplexNumber(0.0, -1.0)) < 1e-9);
    if (!ok) {
        throw std::runtime_error("CRX(pi) should map |10> to ±i|11> (global phase)");
    }
}

void test_cry_gate_on_control_on() {
    const double pi = std::acos(-1.0);
    State s(2, 0);
    s.set_basis_state(0b10, 1.0); // control=1, target=0
    s.cry(1, 0, pi);
    ComplexNumber a10 = s.get_amplitude(0b10);
    ComplexNumber a11 = s.get_amplitude(0b11);
    assert_complex_close(0.0, a10, 1e-9, "CRY(pi) should clear |10>");
    assert_complex_close(1.0, std::abs(a11), 1e-9, "CRY(pi) should map |10> to |11> (up to phase)");
}

void test_cu_gate_matches_cnot_like() {
    const double pi = std::acos(-1.0);
    State s(2, 0);
    s.set_basis_state(0b10, 1.0); // control=1, target=0
    s.cu(1, 0, pi, 0.0, 0.0); // U(pi,0,0) == RY(pi)
    ComplexNumber a10 = s.get_amplitude(0b10);
    ComplexNumber a11 = s.get_amplitude(0b11);
    assert_complex_close(0.0, a10, 1e-9, "CU(pi,0,0) should clear |10>");
    assert_complex_close(1.0, std::abs(a11), 1e-9, "CU(pi,0,0) should map |10> to |11> (up to phase)");
}

void test_ccx_gate_on_110() {
    State s(3, 0);
    s.set_basis_state(0b110, 1.0); // c1=1, c2=1, t=0
    s.ccx(2, 1, 0);
    assert_complex_close(0.0, s.get_amplitude(0b110), 1e-9, "CCX should clear |110>");
    assert_complex_close(1.0, std::abs(s.get_amplitude(0b111)), 1e-9, "CCX should map |110> to |111> (up to phase)");
}

void test_swap_gate() {
    State s(2, 0);
    s.set_basis_state(0b01, 1.0); // |01>
    s.swap(1, 0); // swap qubits -> |10>
    assert_complex_equal(0.0, s.get_amplitude(0b01), "SWAP should clear |01>");
    assert_complex_equal(1.0, s.get_amplitude(0b10), "SWAP should move amplitude to |10>");
}

// --- Helper Function Tests ---

void test_extract_bits_basic() {
    Bitstring b = 0b10001111ULL; // 143
    Bitstring lower = extract_bits(b, 0, 3);
    assert_equal(lower, 15ULL, "extract_bits lower nibble");
    Bitstring bit7 = extract_bits(b, 7, 7);
    assert_equal(bit7, 1ULL, "extract_bits single bit");
}

void test_replace_bits_basic() {
    Bitstring b = 0b10000001ULL; // 129
    Bitstring replaced = replace_bits(b, 0, 3, 0b0100ULL); // 4
    assert_equal(replaced, 0b10000100ULL, "replace_bits preserves upper bits");
}

void test_set_amplitude_add_update() {
    State s(3, 0);
    s.set_amplitude(5ULL, 0.5);
    assert_complex_close(0.5, s.get_amplitude(5ULL), 1e-9, "set_amplitude add");
    size_t size_after_add = s.get_state().size();
    s.set_amplitude(5ULL, 0.25);
    assert_complex_close(0.25, s.get_amplitude(5ULL), 1e-9, "set_amplitude update");
    assert_equal(s.get_state().size(), static_cast<Bitstring>(size_after_add), "set_amplitude does not add duplicate");
}

void test_apply_U0_perp() {
    State s(2, 0);
    s.set_amplitude(0b00, 0.5);
    s.set_amplitude(0b01, 0.5);
    s.apply_U0_perp();
    assert_complex_close(0.5, s.get_amplitude(0b00), 1e-9, "U0_perp should keep |00>");
    assert_complex_close(-0.5, s.get_amplitude(0b01), 1e-9, "U0_perp should flip |01>");
}

void test_grover_oracle_mask() {
    State s(3, 0);
    s.set_amplitude(0b001, 0.5);
    s.set_amplitude(0b010, 0.5);
    std::vector<uint8_t> mask(1ULL << 3, 0);
    mask[0b010] = 1;
    s.grover_oracle_Uf_mask(mask);
    assert_complex_close(0.5, s.get_amplitude(0b001), 1e-9, "Mask oracle should keep non-target");
    assert_complex_close(-0.5, s.get_amplitude(0b010), 1e-9, "Mask oracle should flip target");
}

void test_phase_flip_if() {
    State s(2, 0);
    s.set_amplitude(0b01, 0.5);
    s.set_amplitude(0b10, 0.5);
    s.phase_flip_if([](Bitstring b) { return b == 0b10; });
    assert_complex_close(0.5, s.get_amplitude(0b01), 1e-9, "phase_flip_if should keep |01>");
    assert_complex_close(-0.5, s.get_amplitude(0b10), 1e-9, "phase_flip_if should flip |10>");
}

void test_measure_with_rng() {
    State s(1, 1);
    s.set_amplitude(0b0, 0.6);
    s.set_amplitude(0b1, 0.8); // norm=1
    s.measure_with_rng(0, 0, 0.2); // p0=0.36 -> outcome 0
    assert_equal(static_cast<Bitstring>(s.get_cbits()[0]), 0ULL, "measure_with_rng outcome 0");

    State s2(1, 1);
    s2.set_amplitude(0b0, 0.6);
    s2.set_amplitude(0b1, 0.8);
    s2.measure_with_rng(0, 0, 0.9); // outcome 1
    assert_equal(static_cast<Bitstring>(s2.get_cbits()[0]), 1ULL, "measure_with_rng outcome 1");
}

void test_measure_all_with_rng() {
    State s(2, 2);
    s.set_amplitude(0b00, 0.5);
    s.set_amplitude(0b01, 0.5);
    s.set_amplitude(0b10, 0.5);
    s.set_amplitude(0b11, 0.5);
    std::vector<int> out;
    s.measure_all_with_rng({0.1, 0.9}, out);
    assert_equal(static_cast<Bitstring>(out.size()), 2ULL, "measure_all_with_rng size");
}

void test_measure_all() {
    State s(2, 2);
    s.set_amplitude(0b00, 1.0);
    std::vector<int> out;
    s.measure_all(out);
    assert_equal(static_cast<Bitstring>(out.size()), 2ULL, "measure_all size");
}

void main_core_gate_tests() {
    std::cout << "Testing core gates (X, H, S, T, CX, SWAP):\n";
    std::cout << "------------------------------------------------------------------\n";

    run_test("X gate flips |0> to |1>", test_x_gate_flips_basis);
    run_test("H gate on |0>", test_h_gate_on_zero);
    run_test("H gate on |1>", test_h_gate_on_one);
    run_test("S gate phase on |1>", test_s_gate_phase);
    run_test("Y gate on |0>", test_y_gate_on_zero);
    run_test("Z gate phase on |0> and |1>", test_z_gate_phase);
    run_test("T gate phase on |1>", test_t_gate_phase);
    run_test("RX(pi) on |0> gives -i|1>", test_rx_pi_on_zero);
    run_test("RY(pi) on |0> gives |1>", test_ry_pi_on_zero);
    run_test("RZ(pi) phase on |0> and |1>", test_rz_pi_phases);
    run_test("RU(pi,0,0) matches RY(pi)", test_ru_matches_ry_pi);
    run_test("RU(0,phi,lambda) applies phase to |1>", test_ru_zero_theta_phase);
    run_test("QRNG deterministic bits", test_qrng_deterministic_bits);
    run_test("CX gate (control on)", test_cx_gate_control_on);
    run_test("CX gate (control off)", test_cx_gate_control_off);
    run_test("CNOT alias matches CX", test_cnot_alias_matches_cx);
    run_test("CZ gate (phase on |11>)", test_cz_gate_phase_on_11);
    run_test("CY gate (control on)", test_cy_gate_on_10);
    run_test("CH gate (control on)", test_ch_gate_on_10);
    run_test("CRZ gate (control on)", test_crz_gate_phase_on_control_on);
    run_test("CRX gate (control on)", test_crx_gate_on_control_on);
    run_test("CRY gate (control on)", test_cry_gate_on_control_on);
    run_test("CU gate (control on)", test_cu_gate_matches_cnot_like);
    run_test("CCX gate (control on)", test_ccx_gate_on_110);
    run_test("SWAP gate", test_swap_gate);
    run_test("extract_bits helper", test_extract_bits_basic);
    run_test("replace_bits helper", test_replace_bits_basic);
    run_test("set_amplitude add/update", test_set_amplitude_add_update);
    run_test("apply_U0_perp helper", test_apply_U0_perp);
    run_test("grover_oracle_Uf_mask helper", test_grover_oracle_mask);
    run_test("phase_flip_if helper", test_phase_flip_if);
    run_test("measure_with_rng helper", test_measure_with_rng);
    run_test("measure_all_with_rng helper", test_measure_all_with_rng);
    run_test("measure_all helper", test_measure_all);

    std::cout << "------------------------------------------------------------------\n";
}

// --- QFT / IQFT Tests ---

void test_qft_on_zero_uniform() {
    State s(2, 0);
    s.set_basis_state(0b00, 1.0);
    s.qft(0, 1);
    ComplexNumber coeff = 0.5; // 1/sqrt(4)
    assert_complex_close(coeff, s.get_amplitude(0b00), 1e-6, "QFT |00> -> uniform |00>");
    assert_complex_close(coeff, s.get_amplitude(0b01), 1e-6, "QFT |00> -> uniform |01>");
    assert_complex_close(coeff, s.get_amplitude(0b10), 1e-6, "QFT |00> -> uniform |10>");
    assert_complex_close(coeff, s.get_amplitude(0b11), 1e-6, "QFT |00> -> uniform |11>");
}

void test_qft_iqft_roundtrip() {
    State s(2, 0);
    s.set_amplitude(0b00, 0.5);
    s.set_amplitude(0b01, 0.5);
    s.set_amplitude(0b10, 0.5);
    s.set_amplitude(0b11, -0.5);

    s.qft(0, 1);
    s.iqft(0, 1);

    assert_complex_close(0.5, s.get_amplitude(0b00), 1e-6, "IQFT(QFT) |00>");
    assert_complex_close(0.5, s.get_amplitude(0b01), 1e-6, "IQFT(QFT) |01>");
    assert_complex_close(0.5, s.get_amplitude(0b10), 1e-6, "IQFT(QFT) |10>");
    assert_complex_close(-0.5, s.get_amplitude(0b11), 1e-6, "IQFT(QFT) |11>");
}

void test_qft_on_zero_uniform_gate() {
    State s(2, 0);
    s.set_qft_mode(State::QftMode::Gate);
    s.set_basis_state(0b00, 1.0);
    s.qft(0, 1);
    ComplexNumber coeff = 0.5; // 1/sqrt(4)
    assert_complex_close(coeff, s.get_amplitude(0b00), 1e-6, "Gate QFT |00> -> uniform |00>");
    assert_complex_close(coeff, s.get_amplitude(0b01), 1e-6, "Gate QFT |00> -> uniform |01>");
    assert_complex_close(coeff, s.get_amplitude(0b10), 1e-6, "Gate QFT |00> -> uniform |10>");
    assert_complex_close(coeff, s.get_amplitude(0b11), 1e-6, "Gate QFT |00> -> uniform |11>");
}

void test_qft_iqft_roundtrip_gate() {
    State s(2, 0);
    s.set_qft_mode(State::QftMode::Gate);
    s.set_amplitude(0b00, 0.5);
    s.set_amplitude(0b01, 0.5);
    s.set_amplitude(0b10, 0.5);
    s.set_amplitude(0b11, -0.5);

    s.qft(0, 1);
    s.iqft(0, 1);

    assert_complex_close(0.5, s.get_amplitude(0b00), 1e-6, "Gate IQFT(QFT) |00>");
    assert_complex_close(0.5, s.get_amplitude(0b01), 1e-6, "Gate IQFT(QFT) |01>");
    assert_complex_close(0.5, s.get_amplitude(0b10), 1e-6, "Gate IQFT(QFT) |10>");
    assert_complex_close(-0.5, s.get_amplitude(0b11), 1e-6, "Gate IQFT(QFT) |11>");
}

void main_qft_tests() {
    std::cout << "Testing QFT/IQFT:\n";
    std::cout << "------------------------------------------------------------------\n";
    run_test("QFT on |00> yields uniform superposition", test_qft_on_zero_uniform);
    run_test("QFT then IQFT round-trip", test_qft_iqft_roundtrip);
    run_test("Gate QFT on |00> yields uniform superposition", test_qft_on_zero_uniform_gate);
    run_test("Gate QFT then IQFT round-trip", test_qft_iqft_roundtrip_gate);
    std::cout << "------------------------------------------------------------------\n";
}

// --- Test Cases (r=2, Phase = I) ---

void test_shor_quantum_part_n15_a7_r4() {
  // Test Case: Factor N=15, choosing base a=7. Known period r=4.
  // Precision: m=5 qubits (2^5 = 32 states).
  // Expected output indices (x/2^m) = {0/4, 1/4, 2/4, 3/4} in the control register.
  // Corresponding control register values: {0, 8, 16, 24}.
  // Total qubits = 5 (control) + 4 (target) = 9

  const Bitstring N = 15;
  const Bitstring a = 7;
  const int R = 4;
  const double expected_mag = 1.0 / std::sqrt(static_cast<double>(R * R));

  State s(9, 0);
    
  s.display();
  s.run_shor_algorithm_quantum_part(N, a);

  std::cout << "Qubits 0-3 (target register) hold |a^k mod N> for k in [0, r-1]." << std::endl;
  std::cout << "Qubits 4-8 (control register) hold the period information." << std::endl;
  s.display();

  const Bitstring target_vals[R] = {1ULL, 7ULL, 4ULL, 13ULL};
  const Bitstring control_vals[R] = {0ULL, 8ULL, 16ULL, 24ULL};

  for (int i = 0; i < R; ++i) {
    for (int j = 0; j < R; ++j) {
      Bitstring b = target_vals[i] | (control_vals[j] << 4);
      assert_amplitude_magnitude(s, b, expected_mag, "expected Shor output basis state");
    }
  }

  // Test a state expected to have near zero amplitude (e.g., index 2, which is |00...02>)
  assert_amplitude_magnitude(s, 2ULL, 0.0, "state |00...02>");
    
  // Test another zero state (e.g., control register 1, target register 1)
  assert_amplitude_magnitude(s, 1ULL | (1ULL << 4), 0.0, "state x=1/2^m");
}

void assert_shor_grid(const State& s,
                      const Bitstring* target_vals,
                      const Bitstring* control_vals,
                      int r,
                      double expected_mag)
{
  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < r; ++j) {
      Bitstring b = target_vals[i] | (control_vals[j] << 4);
      assert_amplitude_magnitude(s, b, expected_mag, "expected Shor output basis state");
    }
  }
}

void test_shor_quantum_part_n15_a2_r4() {
  // Test Case: Factor N=15, choosing base a=2. Known period r=4.
  // Expected target values: 1,2,4,8.
  const Bitstring N = 15;
  const Bitstring a = 2;
  const int R = 4;
  const double expected_mag = 1.0 / std::sqrt(static_cast<double>(R * R));

  State s(9, 0);
  s.run_shor_algorithm_quantum_part(N, a);

  const Bitstring target_vals[R] = {1ULL, 2ULL, 4ULL, 8ULL};
  const Bitstring control_vals[R] = {0ULL, 8ULL, 16ULL, 24ULL};
  assert_shor_grid(s, target_vals, control_vals, R, expected_mag);

  // A non-expected state should be near zero.
  assert_amplitude_magnitude(s, 3ULL, 0.0, "state |00...03>");
}

void main_shor_tests() {
    std::cout << "Testing State::run_shor_algorithm_quantum_part (Order Finding):\n";
    std::cout << "------------------------------------------------------------------\n";
    
    run_test("Test 1: N=15, a=7, r=4 (Verification of QPE output concentration)", test_shor_quantum_part_n15_a7_r4);
    run_test("Test 2: N=15, a=2, r=4 (Verification of QPE output concentration)", test_shor_quantum_part_n15_a2_r4);
    
    std::cout << "------------------------------------------------------------------\n";
}

double target_probability(const State& s, int n_t, Bitstring target_val)
{
  double prob = 0.0;
  for (const auto& pair : s.get_state()) {
    Bitstring b = pair.first;
    Bitstring t = extract_bits(b, 0, n_t - 1);
    if (t == target_val) {
      prob += std::norm(pair.second);
    }
  }
  return prob;
}

void test_shor_small_semiprimes() {
  struct Case { Bitstring N; Bitstring a; int r; };
  const Case cases[] = {
    {21, 2, 6},  // 21 = 3*7
    {33, 10, 2}, // 33 = 3*11
    {35, 2, 12}, // 35 = 5*7
    {39, 2, 12}, // 39 = 3*13
  };

  for (const auto& c : cases) {
    int n_t = static_cast<int>(std::ceil(std::log2(static_cast<double>(c.N))));
    int n_c = 2 * n_t;
    int total_qubits = n_c + n_t;
    State s(total_qubits, 0);
    s.run_shor_algorithm_quantum_part(c.N, c.a);

    std::vector<Bitstring> expected_targets;
    expected_targets.reserve(c.r);
    for (int k = 0; k < c.r; ++k) {
      expected_targets.push_back(mod_pow(c.a, static_cast<Bitstring>(k), c.N));
    }

    double min_expected = 0.03;
    for (Bitstring t : expected_targets) {
      double prob = target_probability(s, n_t, t);
      if (prob < min_expected) {
        throw std::runtime_error("Expected target probability too low in Shor small semiprime test.");
      }
    }

    // Check non-targets are small.
    Bitstring max_target = (1ULL << n_t);
    for (Bitstring t = 0; t < max_target; ++t) {
      bool is_expected = false;
      for (Bitstring et : expected_targets) {
        if (et == t) {
          is_expected = true;
          break;
        }
      }
      if (is_expected) continue;
      double prob = target_probability(s, n_t, t);
      if (prob > 0.05) {
        throw std::runtime_error("Unexpected target probability too high in Shor small semiprime test.");
      }
    }
  }
}

int run_unit_tests()
{
  bool verbose = false;
  if (const char* env = std::getenv("QSIM_TEST_VERBOSE")) {
    verbose = (env[0] != '\0' && env[0] != '0');
  }
  bool slow = false;
  if (const char* env = std::getenv("QSIM_SLOW_TESTS")) {
    slow = (env[0] != '\0' && env[0] != '0');
  }
  bool demo = false;
  if (const char* env = std::getenv("QSIM_DEMO_TESTS")) {
    demo = (env[0] != '\0' && env[0] != '0');
  }
  if (const char* env = std::getenv("QSIM_TEST_SEED")) {
    unsigned int seed = static_cast<unsigned int>(std::strtoul(env, nullptr, 10));
    std::srand(seed);
  }
  State::set_default_log_stream(&std::cout);
  setenv("QSIM_GROVER_VERBOSE", "1", 1);
  run_test("Logging levels", test_logging_levels);
  qsim_log::set_level(qsim_log::Level::Verbose);
  main_test_controlled_Rr();
  main_test_controlled_Rr_dag();
  main_test_qft();
  main_mod_exp();
  main_all_cme_tests();
  main_core_gate_tests();
  main_qft_tests();
  run_test("Mod arithmetic (gcd)", test_mod_arith_gcd_basic);
  run_test("Mod arithmetic (mod_pow)", test_mod_arith_mod_pow_edges);
  run_test("QRNG edge cases", test_qrng_edges);
  run_test("Deutsch-Jozsa constant oracles", test_deutsch_jozsa_constant_oracles);
  run_test("Deutsch-Jozsa balanced oracles", test_deutsch_jozsa_balanced_oracles);
  run_test("Deutsch-Jozsa oracle helpers", test_deutsch_jozsa_oracle_helpers);
  run_test("Bernstein-Vazirani recovery", test_bernstein_vazirani_secret_recovery);
  run_test("Bernstein-Vazirani errors", test_bernstein_vazirani_errors);
  run_test("QUBO exact solver", test_qubo_exact_solver);
  run_test("QUBO Grover threshold", test_qubo_grover_threshold_solver);
  run_test("QUBO error paths", test_qubo_error_paths);
  run_test("TSP QUBO + exact solver", test_tsp_qubo_and_exact_solver);
  run_test("TSP demo paths", test_tsp_demo_paths);
  run_test("VQA QAOA good candidate", test_vqa_qaoa_finds_good_candidate);
  run_test("VQA QAOA improves objective", test_vqa_qaoa_improves_over_initial_state);
  run_test("VQA QAOA error paths", test_vqa_qaoa_error_paths);
  run_test("VQA QAOA shot mode path", test_vqa_qaoa_shot_mode_path);
  run_test("QAOA wrapper paths", test_qaoa_wrapper_paths);
  run_test("QAOA demo paths", test_qaoa_demo_paths);
  run_test("VQE single-qubit ground state", test_vqe_single_qubit_ground_state);
  run_test("VQE expectation Pauli terms", test_vqe_expectation_pauli_terms);
  run_test("VQE error and edge paths", test_vqe_error_and_edge_paths);
  run_test("Anneal SA finds good candidate", test_anneal_sa_finds_good_candidate);
  run_test("Anneal SQA finds good candidate", test_anneal_sqa_finds_good_candidate);
  run_test("Anneal option + parsing paths", test_anneal_options_and_parsing);
  run_test("Anneal SQA improvement branch", test_anneal_sqa_improvement_branch);
  run_test("Grover search helpers", test_grover_search_helpers);
  run_test("Grover API errors", test_grover_api_errors);
  run_test("Display output paths", test_display_output_paths);
  run_test("QFT invalid ranges", test_qft_invalid_ranges);
  run_test("QFT tiny amplitude continue", test_qft_tiny_amplitude_continue);
  run_test("Measurement branches", test_measure_branches);
  run_test("CLI command parsing", test_cli_command_parsing);
  run_test("Latin square validations", test_latin_square_validations);
  if (demo) {
    run_test("Latin square forced measurement", test_latin_square_demo_forced_measure);
  } else {
    std::cout << "Latin square demo tests skipped. Set QSIM_DEMO_TESTS=1 to run them.\n";
  }
  run_test("Latin square invalid row0", test_latin_square_invalid_row0);
  if (demo) {
    run_test("Latin square count/print", test_latin_square_count_print);
  } else {
    std::cout << "Latin square demo tests skipped. Set QSIM_DEMO_TESTS=1 to run them.\n";
  }
  run_test("Shor API paths", test_shor_api_paths);
  if (slow) {
    run_test("Shor classical estimate_order", test_shor_classical_estimate_order);
    run_test("Shor quantum free function", test_shor_quantum_free_function_paths);
  } else {
    std::cout << "Shor classical/quantum helper tests skipped. Set QSIM_SLOW_TESTS=1 to run them.\n";
  }
  if (demo) {
    run_test("Shor demo branches", test_shor_demo_branches);
    run_test("Algorithm demo wrappers", test_algorithm_demo_wrappers);
  } else {
    std::cout << "Shor demo branches skipped. Set QSIM_DEMO_TESTS=1 to run them.\n";
  }
  if (verbose && slow) {
    main_shor_tests();
    run_test("Small semiprimes order-finding (Shor)", test_shor_small_semiprimes);
  } else if (verbose) {
    std::cout << "Shor heavy tests skipped. Set QSIM_SLOW_TESTS=1 to run them.\n";
  } else {
    std::cout << "Shor tests skipped. Set QSIM_TEST_VERBOSE=1 to run them.\n";
  }

  return g_failures;
}
