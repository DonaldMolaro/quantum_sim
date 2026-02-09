#include "state.hh"
#include "tests/helpers.hh"
#include "algorithms/qrng.hh"
#include "algorithms/latin_square.hh"
#include "algorithms/api/grover_api.hh"
#include "algorithms/api/shor_api.hh"
#include "algorithms/shor_classical.hh"
#include "algorithms/shor_quantum.hh"
#include "cli/commands.hh"
#include "demos/grover_demo.hh"
#include "demos/latin_demo.hh"
#include "demos/shor_demo.hh"
#include "modular_exp.hh"
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

static void assert_double_close(double actual, double expected, double tol, const std::string& message)
{
    if (std::abs(actual - expected) > tol) {
        throw std::runtime_error("Assertion failed: " + message +
                                 " Expected: " + std::to_string(expected) +
                                 " Actual: " + std::to_string(actual));
    }
}


using test_helpers::assert_complex_equal;
using test_helpers::assert_complex_close;
using test_helpers::assert_equal;
using test_helpers::assert_amplitude_match;
using test_helpers::assert_amplitude_magnitude;

struct ScopedEnv {
    std::string key;
    std::string old_value;
    bool had_old;

    ScopedEnv(const std::string& k, const std::string& v) : key(k) {
        const char* old = std::getenv(key.c_str());
        had_old = (old != nullptr);
        if (had_old) {
            old_value = old;
        }
        setenv(key.c_str(), v.c_str(), 1);
    }

    ~ScopedEnv() {
        if (had_old) {
            setenv(key.c_str(), old_value.c_str(), 1);
        } else {
            unsetenv(key.c_str());
        }
    }
};

static void assert_exits_with_failure(const std::function<void()>& fn)
{
    pid_t pid = fork();
    if (pid == 0) {
        fn();
        std::exit(0);
    }
    if (pid < 0) {
        throw std::runtime_error("fork failed for exit test");
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        throw std::runtime_error("waitpid failed for exit test");
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) == 0) {
        throw std::runtime_error("expected child to exit with failure");
    }
}

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

const int CONTROL_Q = 1;
const int TARGET_Q = 0;
const int R_VALUE = 2; // Produces phase factor I

// Test Case 1: Control=0, Target=0 (State |00>, Bitstring 0)
void test_crr_00_no_change() {
    State s(2, 0);
    s.set_basis_state(0b00, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE); 
    assert_complex_equal(1.0, s.get_amplitude(0b00), "00 state should be unchanged.");
}

// Test Case 2: Control=0, Target=1 (State |01>, Bitstring 1)
void test_crr_01_no_change() {
    State s(2, 0);
    s.set_basis_state(0b01, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE);
    // Expect amplitude 1.0 (control bit 1 is 0)
    assert_complex_equal(1.0, s.get_amplitude(0b01), "01 state should be unchanged.");
}

// Test Case 3: Control=1, Target=0 (State |10>, Bitstring 2)
void test_crr_10_no_change() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE);
    // Expect amplitude 1.0 (target bit 0 is 0)
    assert_complex_equal(1.0, s.get_amplitude(0b10), "10 state should be unchanged.");
}

// Test Case 4: Control=1, Target=1 (State |11>, Bitstring 3)
void test_crr_11_applies_phase_I() {
    State s(2, 0);
    s.set_basis_state(0b11, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE); 
    // Expect amplitude I (phase rotation by pi/2 for r=2)
    assert_complex_equal(I, s.get_amplitude(0b11), "11 state should be multiplied by I.");
}

// Test Case 5: Superposition Test (Linearity)
void test_crr_superposition() {
    State s(2, 0);
    // State: 1/sqrt(2) |00> + 1/sqrt(2) |11>
    ComplexNumber coeff = 1.0 / std::sqrt(2.0);
    QuantumState initial_state = { {0b00, coeff}, {0b11, coeff} };
    s.set_superposition(initial_state);
    
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE);
    
    // |00> amplitude should remain unchanged (coeff)
    assert_complex_equal(coeff, s.get_amplitude(0b00), "|00> amplitude failure (linearity).");
    
    // |11> amplitude should be multiplied by I (I * coeff)
    assert_complex_equal(I * coeff, s.get_amplitude(0b11), "|11> amplitude failure (linearity).");
}


// --- Main Execution Block (Conceptual) ---

int main_test_controlled_Rr() {
    std::cout << "Testing State::controlled_Rr (Control=q1, Target=q0, r=2, Phase=I):\n";
    std::cout << "------------------------------------------------------------------\n";
    
    run_test("Test 1: |00> (C=0, T=0) -> No Change", test_crr_00_no_change);
    run_test("Test 2: |01> (C=0, T=1) -> No Change", test_crr_01_no_change);
    run_test("Test 3: |10> (C=1, T=0) -> No Change", test_crr_10_no_change);
    run_test("Test 4: |11> (C=1, T=1) -> Apply Phase I", test_crr_11_applies_phase_I);
    run_test("Test 5: Superposition Test (1/sqrt(2) |00> + 1/sqrt(2) |11>)", test_crr_superposition);
    
    std::cout << "------------------------------------------------------------------\n";

    return 0;
}

// Test Case 1: Control=0, Target=0 (State |00>, Bitstring 0)
void test_crrdag_00_no_change() {
    State s(2, 0);
    s.set_basis_state(0b00, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE); 
    assert_complex_equal(1.0, s.get_amplitude(0b00), "00 state should be unchanged.");
}

// Test Case 2: Control=0, Target=1 (State |01>, Bitstring 1)
void test_crrdag_01_no_change() {
    State s(2, 0);
    s.set_basis_state(0b01, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE);
    assert_complex_equal(1.0, s.get_amplitude(0b01), "01 state should be unchanged.");
}

// Test Case 3: Control=1, Target=0 (State |10>, Bitstring 2)
void test_crrdag_10_no_change() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE);
    assert_complex_equal(1.0, s.get_amplitude(0b10), "10 state should be unchanged.");
}

// Test Case 4: Control=1, Target=1 (State |11>, Bitstring 3) - Applies phase -I
void test_crrdag_11_applies_phase_minus_I() {
    State s(2, 0);
    s.set_basis_state(0b11, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE); 
    // Expected result: -I
    assert_complex_equal(-I, s.get_amplitude(0b11), "11 state should be multiplied by -I.");
}

// Test Case 5: Verification that R_r * R_r_dag = I
void test_crrdag_inversion() {
    State s(2, 0);
    s.set_basis_state(0b11, 1.0);
    
    // Step 1: Apply R_r (Phase I)
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE);
    // Check intermediate state: should be I
    assert_complex_equal(I, s.get_amplitude(0b11), "Intermediate state failed (should be I).");
    
    // Step 2: Apply R_r_dag (Phase -I)
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE);
    
    // Final state should be 1.0 + 0.0i (I * -I = -I^2 = 1)
    assert_complex_equal(1.0, s.get_amplitude(0b11), "R_r * R_r_dag failed to produce Identity.");
}

// Test Case 6: Superposition Test (Linearity)
void test_crrdag_superposition() {
    State s(2, 0);
    // State: 1/sqrt(2) |00> + 1/sqrt(2) |11>
    ComplexNumber coeff = 1.0 / std::sqrt(2.0);
    QuantumState initial_state = { {0b00, coeff}, {0b11, coeff} };
    s.set_superposition(initial_state);
    
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE);
    
    // |00> amplitude should remain unchanged (coeff)
    assert_complex_equal(coeff, s.get_amplitude(0b00), "|00> amplitude failure (linearity).");
    
    // |11> amplitude should be multiplied by -I (-I * coeff)
    assert_complex_equal(-I * coeff, s.get_amplitude(0b11), "|11> amplitude failure (linearity).");
}

int main_test_controlled_Rr_dag() {
    std::cout << "Testing State::controlled_Rr_dag (Control=q1, Target=q0, r=2, Phase=-I):\n";
    std::cout << "----------------------------------------------------------------------\n";
    
    // Test that the inverse phase is only applied when control=1 AND target=1
    run_test("Test 1: |00> (C=0, T=0) -> No Change", test_crrdag_00_no_change);
    run_test("Test 2: |01> (C=0, T=1) -> No Change", test_crrdag_01_no_change);
    run_test("Test 3: |10> (C=1, T=0) -> No Change", test_crrdag_10_no_change);
    
    // Test the specific condition where the inverse phase must be applied
    run_test("Test 4: |11> (C=1, T=1) -> Apply Phase -I", test_crrdag_11_applies_phase_minus_I);

    // Test for composition yielding identity
    run_test("Test 5: R_r * R_r_dag = I", test_crrdag_inversion);

    // Test that the operation correctly handles superposition (linearity)
    run_test("Test 6: Superposition Test (1/sqrt(2) |00> + 1/sqrt(2) |11>)", test_crrdag_superposition);
    
    std::cout << "----------------------------------------------------------------------\n";

    return 0;
}

// --- Test Cases (N=4, 2 Qubits: 0 and 1) ---
const int START = 0;
const int END = 1;

// Test 1: QFT on |00> (j=0) -> Should yield uniform superposition (amplitude 0.5 for all states).
void test_qft_zero_state() {
    State s(2, 0);
    s.set_basis_state(0b00, 1.0);
    s.qft(START, END);

    ComplexNumber expected_amp = 0.5;
    assert_complex_equal(expected_amp, s.get_amplitude(0b00), "|00> amplitude check");
    assert_complex_equal(expected_amp, s.get_amplitude(0b01), "|01> amplitude check");
    assert_complex_equal(expected_amp, s.get_amplitude(0b10), "|10> amplitude check");
    assert_complex_equal(expected_amp, s.get_amplitude(0b11), "|11> amplitude check");
}

// Test 2: QFT on |01> (j=1) -> Should yield characteristic phases (1/2, i/2, -1/2, -i/2).
void test_qft_phase_state_j1() {
    State s(2, 0);
    s.set_basis_state(0b01, 1.0);
    s.qft(START, END);
    
    // Output amplitudes scaled by 1/sqrt(4) = 0.5. Phase is e^(i*pi*k/2)
    // k=0: 0.5 * e^0 = 0.5
    // k=1: 0.5 * e^(i*pi/2) = 0.5 * I
    // k=2: 0.5 * e^(i*pi) = -0.5
    // k=3: 0.5 * e^(i*3pi/2) = -0.5 * I
    
    assert_complex_equal(0.5, s.get_amplitude(0b00), "|00> amplitude check");
    assert_complex_equal(0.5 * I, s.get_amplitude(0b01), "|01> amplitude check");
    assert_complex_equal(-0.5, s.get_amplitude(0b10), "|10> amplitude check");
    assert_complex_equal(-0.5 * I, s.get_amplitude(0b11), "|11> amplitude check");
}

// Test 3: QFT on |10> (j=2) -> Should yield alternating signs (1/2, -1/2, 1/2, -1/2).
void test_qft_phase_state_j2() {
    State s(2, 0);
    s.set_basis_state(0b10, 1.0);
    s.qft(START, END);
    
    // Phase is e^(i*pi*2*k/4) = e^(i*pi*k)
    // k=0: 0.5
    // k=1: -0.5
    // k=2: 0.5
    // k=3: -0.5
    
    assert_complex_equal(0.5, s.get_amplitude(0b00), "|00> amplitude check");
    assert_complex_equal(-0.5, s.get_amplitude(0b01), "|01> amplitude check");
    assert_complex_equal(0.5, s.get_amplitude(0b10), "|10> amplitude check");
    assert_complex_equal(-0.5, s.get_amplitude(0b11), "|11> amplitude check");
}


// Test 4: IQFT applied to the uniform superposition (result of QFT|00>) must return |00>.
void test_iqft_inversion_check() {
    State s(2, 0);
    s.set_basis_state(0b00, 1.0);
    
    // 1. Apply QFT -> creates uniform superposition
    s.qft(START, END); 
    
    // 2. Apply IQFT (using the special IQFT mock for superposition input)
    s.iqft(START, END); 
    
    // Expect: |00> amplitude 1.0, all others 0.0
    assert_complex_equal(1.0, s.get_amplitude(0b00), "Final |00> amplitude must be 1.0");
    assert_complex_equal(0.0, s.get_amplitude(0b11), "Final |11> amplitude must be 0.0");
}

// Test 5: Single qubit IQFT (N=2, indices 0-0) -> Reversal of Hadamard.
void test_iqft_single_qubit_h_inverse() {
    // Note: This requires redefining START/END for N=2 or mocking H->IQFT|0>.
    // Using the current N=4 architecture: Test IQFT on |00> (j=0) -> |00>
    State s(2, 0);
    s.set_basis_state(0b00, 1.0);
    s.iqft(START, END);

    // Phase factor e^(-2*pi*i*0*k / 4) = 1 for all k. Scaled by 1/2.
    ComplexNumber expected_amp = 0.5;
    assert_complex_equal(expected_amp, s.get_amplitude(0b00), "|00> amplitude check");
    assert_complex_equal(expected_amp, s.get_amplitude(0b01), "|01> amplitude check");
}


int main_test_qft() {
    std::cout << "Testing State::qft and State::iqft (N=4):\n";
    std::cout << "---------------------------------------\n";
    
    run_test("QFT |00> -> Uniform Superposition", test_qft_zero_state);
    run_test("QFT |01> -> Phase Shift N=4", test_qft_phase_state_j1);
    run_test("QFT |10> -> Alternating Sign N=4", test_qft_phase_state_j2);
    run_test("QFT followed by IQFT returns Identity", test_iqft_inversion_check);
    run_test("IQFT applied to |00> (IQFT definition check)", test_iqft_single_qubit_h_inverse);

    std::cout << "---------------------------------------\n";
    return 0;
}




void test_mod_exp_zero_power() {
    // Test case: a^0 mod N = 1
    Bitstring result = modular_exponentiation(42ULL, 0ULL, 100ULL);
    assert_equal(result, 1ULL, "Zero Power Test (42^0 mod 100)");
}

void test_mod_exp_mod_one() {
    // Test case: mod 1 should always yield 0
    Bitstring result = modular_exponentiation(5ULL, 123ULL, 1ULL);
    assert_equal(result, 0ULL, "Mod 1 Test (5^123 mod 1)");
}

void test_mod_exp_one_power() {
    // Test case: a^1 mod N = a mod N
    Bitstring result = modular_exponentiation(5ULL, 1ULL, 7ULL);
    assert_equal(result, 5ULL, "One Power Test (5^1 mod 7)");
    
    // Test case where a > N
    result = modular_exponentiation(10ULL, 1ULL, 7ULL);
    assert_equal(result, 3ULL, "One Power Test (10^1 mod 7)");
}

void test_mod_exp_standard() {
    // Test case: 3^3 mod 5 = 27 mod 5 = 2
    Bitstring result = modular_exponentiation(3ULL, 3ULL, 5ULL);
    assert_equal(result, 2ULL, "Standard Test (3^3 mod 5)");

    // Test case: 2^10 mod 100 = 1024 mod 100 = 24
    result = modular_exponentiation(2ULL, 10ULL, 100ULL);
    assert_equal(result, 24ULL, "Standard Test (2^10 mod 100)");
}

void test_mod_exp_shor_period() {
    // Shor's relevance: a=7, N=15. Period r=4.
    
    // Test 7^2 mod 15 = 49 mod 15 = 4
    Bitstring result = modular_exponentiation(7ULL, 2ULL, 15ULL);
    assert_equal(result, 4ULL, "Shor Period Test (7^2 mod 15)");

    // Test 7^4 mod 15 = 1 (End of period)
    result = modular_exponentiation(7ULL, 4ULL, 15ULL);
    assert_equal(result, 1ULL, "Shor Period Test (7^4 mod 15)");
    
    // Test 7^8 mod 15 = 1 (Used in CME controlled by q7)
    result = modular_exponentiation(7ULL, 8ULL, 15ULL);
    assert_equal(result, 1ULL, "Shor Period Test (7^8 mod 15)");

    // Test 7^16 mod 15 = 1 (Used in CME controlled by q8)
    result = modular_exponentiation(7ULL, 16ULL, 15ULL);
    assert_equal(result, 1ULL, "Shor Period Test (7^16 mod 15)");
}

/*
To test the **Controlled Modular Exponentiation (CME)** operation with parameters $a=7, P=1$, and $N=15$, we follow the **operational view** where the quantum state is a set of bitstring-amplitude pairs and the operation is a functional transformation.

### Test Parameters and Register Setup
*   **Total Register:** 7 Qubits (Index 0-6).
*   **Target Register:** Qubits $q_0$ through $q_3$ (to store values up to $N=15$).
*   **Control Qubit:** Qubit $q_6$.
*   **Modular Function:** $U_7^1 |y\rangle = |7^1 \cdot y \pmod{15}\rangle$.
*/

void test_cme_superposition_7_1_15() {
  //std::cout << " 1. Create initial state |0000000> " << std::endl;
  State s(7,1); 
  //s.display();
  //std::cout << " 2. Initialize Target Register to |1> (Bitstring 0000001) " << std::endl;
  //std::cout << "Multiplicative identity is required: y_new = (a^P * 1) mod N " << std::endl;
  s.x(0); 
  //s.display();
  //std::cout << " 3. Put Control Qubit (q6) into Superposition " << std::endl;
  //std::cout << "State becomes 1/sqrt(2) * (|0000001> + |1000001>)" << std::endl;
  s.h(6);
  //s.display();

  //std::cout << " 4. Execute CME: a=7, P=1, N=15, Control=q6, Target=" << std::endl;
  s.controlled_modular_exponentiation(6, 0, 3, 7, 15, 1);
  //s.display();

  // 5. Validation Logic
  // Expected Outcome 1 (Control 0): Identity applies. State remains |0000001> (Decimal 1)
  // Expected Outcome 2 (Control 1): U applies. 1 * 7 mod 15 = 7. 
  // Bitstring: |1(q6) 00(q5-q4) 0111(q3-q0)> = Binary 1000111 (Decimal 71)
    
  double prob1 = std::norm(s.get_amplitude(1));  // |0000001>
  double prob2 = std::norm(s.get_amplitude(71)); // |1000111>

    if (std::abs(prob1 - 0.5) < 1e-6 && std::abs(prob2 - 0.5) < 1e-6) {
        std::cout << "[PASS] CME Superposition Test" << std::endl;
    } else {
        std::cout << "[FAIL] CME Superposition Test" << std::endl;
        s.display(); // Detailed state dump in binary
    }
}

void test_mod_exp_large_power() {
    // Test 2^63 mod 100 = 8
    // Uses a power close to the maximum Bitstring value (if Bitstring is 64-bit).
    Bitstring result = modular_exponentiation(2ULL, 63ULL, 100ULL);
    assert_equal(result, 8ULL, "Large Power Test (2^63 mod 100)");
}

int main_mod_exp() {
     std::cout << "Testing modular_exponentiation utility:" << std::endl;
     run_test("Test 1: Zero Power", test_mod_exp_zero_power);
     run_test("Test 2: Modulo 1", test_mod_exp_mod_one);
     run_test("Test 2: One Power", test_mod_exp_one_power);
     run_test("Test 3: Standard Exponentiation", test_mod_exp_standard);
     run_test("Test 4: Shor Period Check", test_mod_exp_shor_period);
     run_test("Test 5: Large Exponent Check", test_mod_exp_large_power);
     return 0;
}


// --- Constants specific to the Shor test case ---
const int TARGET_START = 0;
const int TARGET_END = 3; 
//const int N_TARGET = 4; // Quanta in the target register (log2(15) rounded up)
const Bitstring N = 15; // Modulus
const Bitstring A = 7; // Base 'a'
const int CONTROL_QUBIT = 6; // A representative index outside the target range (e.g., q6)

// --- Unit Tests Helper ---

// --- Unit Tests Implementation ---

void test_cme_control_off() {
    // Test 1: Control qubit is 0. Operation should be Identity.
    // Start state: |0>_q6 |1>_target (Bitstring 1). Power=7.
    const Bitstring START_B = 1ULL;
    
    State s(7,4);
    s.x(0);
    
    // CME call must include TARGET_START and TARGET_END indices
    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 7);

    assert_amplitude_match(s, START_B, 1.0, "Control OFF: State should remain |1>");
    assert_amplitude_match(s, 7ULL, 0.0, "Control OFF: Other states must be zero");
}

void test_cme_control_on_p1() {
    // Test 2: Control qubit is 1, Power=1. 
    // Start state: |1>_q6 |1>_target (Bitstring 65). Expected result: 7^1 * 1 mod 15 = 7. 
    // Final state |1>_q6 |7>_target (Bitstring 71).
    
    const Bitstring START_B = 1ULL | (1ULL << CONTROL_QUBIT); // 65
    const Bitstring EXPECTED_TARGET_VAL = 7;
    const Bitstring EXPECTED_B = EXPECTED_TARGET_VAL | (1ULL << CONTROL_QUBIT); // 71

    //std::cout << "START_B (" << bitstring_to_string(START_B,15) << ")" << std::endl;
    State s(START_B,7);
    //s.display();
    
    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 1);
    
    //s.display();

    
    assert_amplitude_match(s, EXPECTED_B, 1.0, "Control ON P=1: State should shift to |7>");
    assert_amplitude_match(s, START_B, 0.0, "Control ON P=1: Initial state must be zero");
}

void test_cme_control_on_p2() {
    // Test 3: Control qubit is 1, Power=2. 
    // Start state: |1>_q6 |1>_target (Bitstring 65). Expected result: 7^2 * 1 mod 15 = 4. 
    // Final state |1>_q6 |4>_target (Bitstring 68).
    
    const Bitstring START_B = 1ULL | (1ULL << CONTROL_QUBIT); // 65
    const Bitstring EXPECTED_TARGET_VAL = 4;
    const Bitstring EXPECTED_B = EXPECTED_TARGET_VAL | (1ULL << CONTROL_QUBIT); // 68

    State s(START_B,7); 
    
    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 2);

    assert_amplitude_match(s, EXPECTED_B, 1.0, "Control ON P=2: State should shift to |4>");
    assert_amplitude_match(s, START_B, 0.0, "Control ON P=2: Initial state must be zero");
}

void test_cme_superposition_input() {
    // Test 4: Input 1/sqrt(2) (|0>|1> + |1>|1>). Power=1.
    // Expected: 1/sqrt(2) (|0>|1> + |1>|7>).
    
    const double AMP = 1.0 / std::sqrt(2.0);
    const Bitstring B_OFF = 1ULL; // |0>_q6 |1>_target
    const Bitstring B_ON_IN = 1ULL | (1ULL << CONTROL_QUBIT); // |1>_q6 |1>_target
    
    // Expected Output States:
    const Bitstring B_OFF_OUT = B_OFF; // Target remains |1> (7^0 mod 15 * 1 mod 15 = 1)
    const Bitstring B_ON_OUT = 7ULL | (1ULL << CONTROL_QUBIT); // Target shifts to |7> (7^1 mod 15 * 1 mod 15 = 7)

    State s(0, 0);
    s.set_amplitude(B_OFF, AMP);
    s.set_amplitude(B_ON_IN, AMP);
    
    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 1);

    assert_amplitude_match(s, B_OFF_OUT, AMP, "Superposition: Control OFF component (B=1) maintained magnitude");
    assert_amplitude_match(s, B_ON_OUT, AMP, "Superposition: Control ON component (B=71) maintained magnitude");
    assert_amplitude_match(s, B_ON_IN, 0.0, "Superposition: Initial ON input state must vanish");
}

void test_cme_out_of_range_target_identity() {
    // Test 5: If target value y >= N, the CME should act as identity (unitary embedding).
    // Use target value 15 (0b1111) with N=15 (out of range), control ON.
    const Bitstring START_B = 15ULL | (1ULL << CONTROL_QUBIT);
    State s(START_B, 7);

    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 1);

    assert_amplitude_match(s, START_B, 1.0, "Out-of-range target should remain unchanged");
}

void test_cme_zero_target_control_on() {
    // Test 7: If target value y=0, result should remain 0 when control is on.
    const Bitstring START_B = 0ULL | (1ULL << CONTROL_QUBIT);
    State s(START_B, 7);

    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 1);

    assert_amplitude_match(s, START_B, 1.0, "Target 0 should remain unchanged under CME");
}


// --- Main Test Harness ---

void run_test_cme(const std::string& name, void (*test_func)()) {
    try {
        test_func();
        std::cout << "[PASS] " << name << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "[FAIL] " << name << ": " << e.what() << std::endl;
    }
}

void main_all_cme_tests() {
    std::cout << "\nTesting State::controlled_modular_exponentiation (Modular Exponentiation):\n";
    std::cout << "------------------------------------------------------------------\n";
    
    run_test_cme("Test 1: Control OFF (Power=7)", test_cme_control_off);
    run_test_cme("Test 2: Control ON, Power=1 (7^1 mod 15 = 7)", test_cme_control_on_p1);
    run_test_cme("Test 3: Control ON, Power=2 (7^2 mod 15 = 4)", test_cme_control_on_p2);
    run_test_cme("Test 4: Superposition Input (Check entanglement)", test_cme_superposition_input);
    run_test_cme("Test 5: Superposition 7 1 15 (Check entanglement)", test_cme_superposition_7_1_15);
    run_test_cme("Test 6: Out-of-range target acts as identity", test_cme_out_of_range_target_identity);
    run_test_cme("Test 7: Zero target with control on", test_cme_zero_target_control_on);
    
    std::cout << "------------------------------------------------------------------\n";
}

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

void test_mod_arith_gcd_basic() {
    assert_equal(gcd_bitstring(0, 7), 7ULL, "gcd(0,7)=7");
    assert_equal(gcd_bitstring(54, 24), 6ULL, "gcd(54,24)=6");
}

void test_mod_arith_mod_pow_edges() {
    assert_equal(mod_pow(5, 3, 1), 0ULL, "mod_pow with mod=1 returns 0");
    assert_equal(mod_pow(2, 0, 15), 1ULL, "mod_pow exponent 0 returns 1");
    assert_equal(mod_pow(7, 4, 15), 1ULL, "mod_pow periodicity check");
}

void test_qrng_edges() {
    std::vector<int> empty = qrng_bits(0, nullptr);
    if (!empty.empty()) {
        throw std::runtime_error("qrng_bits(0) should return empty");
    }

    std::vector<int> bits_auto = qrng_bits(2, nullptr);
    if (bits_auto.size() != 2) {
        throw std::runtime_error("qrng_bits should return 2 bits");
    }

    std::vector<double> rv = {0.0, 1.0};
    std::vector<int> bits = qrng_bits(2, &rv);
    if (bits.size() != 2 || bits[0] != 0 || bits[1] != 1) {
        throw std::runtime_error("qrng_bits with deterministic RNG failed");
    }

    ScopedEnv cap("QSIM_QRNG_MAX_QUBITS", "8");
    std::vector<double> zeros(8, 0.0);
    uint64_t val = qrng_u64(70, &zeros);
    if (val != 0ULL) {
        throw std::runtime_error("qrng_u64 expected all-zero output with zero RNG");
    }

    if (qrng_u64(0, nullptr) != 0ULL) {
        throw std::runtime_error("qrng_u64(0) should return 0");
    }
}

void test_grover_search_helpers() {
    run_grover_search(nullptr, 1);
    State s(2, 0);
    run_grover_search(&s, 3);
    run_grover_search_multi(&s, std::vector<Bitstring>());
    run_grover_search_multi(&s, std::vector<Bitstring>{1, 2});
}

void test_grover_api_errors() {
    State s(2, 0);
    GroverResult empty = run_grover(s, std::vector<Bitstring>());
    if (empty.ok || empty.error.empty()) {
        throw std::runtime_error("Expected error on empty target list");
    }

    State big(63, 0);
    GroverResult too_big = run_grover(big, std::vector<Bitstring>{1});
    if (too_big.ok || too_big.error.empty()) {
        throw std::runtime_error("Expected error for too many qubits");
    }

    GroverResult out_of_range = run_grover(s, std::vector<Bitstring>{7});
    if (out_of_range.ok || out_of_range.error.empty()) {
        throw std::runtime_error("Expected out-of-range target error");
    }

    State s2(3, 0);
    GroverResult r0 = run_grover(s2, std::vector<Bitstring>{1}, 0);
    if (!r0.ok || r0.iterations != 0) {
        throw std::runtime_error("Expected R=0 Grover run to succeed");
    }

    {
        ScopedEnv env("QSIM_GROVER_FORCE_SET", "1");
        GroverResult forced = run_grover(s2, std::vector<Bitstring>{1}, 1);
        if (!forced.ok) {
            throw std::runtime_error("Expected forced set oracle Grover run to succeed");
        }
    }
}

void test_display_output_paths() {
    State s(2, 2);
    s.set_amplitude(0, ComplexNumber(1.0, 0.0));
    s.set_amplitude(1, ComplexNumber(0.0, 1.0));
    s.set_amplitude(2, ComplexNumber(0.0, -1.0));
    s.set_amplitude(3, ComplexNumber(2.0, 3.0));

    s.display(false);
    s.display(true);

    s.measure(0, 0);
    s.measure(1, 1);
    s.display_cbits();

    State empty_cbits(1, 0);
    empty_cbits.display_cbits();

    State sparse_only(2, 0);
    sparse_only.set_amplitude(0, ComplexNumber(1.0, 0.0));
    sparse_only.display(true);

    State one_bit(1, 1);
    one_bit.set_basis_state(1, ONE_COMPLEX);
    one_bit.measure_with_rng(0, 0, 1.0);
    one_bit.display_cbits();
}

void test_qft_invalid_ranges() {
    State s(2, 0);
    s.set_qft_mode(State::QftMode::Direct);
    s.qft(1, 0);
    s.iqft(1, 0);

    s.set_qft_mode(State::QftMode::Gate);
    s.qft(1, 0);
    s.iqft(1, 0);
}

void test_qft_tiny_amplitude_continue() {
    State s(2, 0);
    s.set_qft_mode(State::QftMode::Direct);
    s.set_amplitude(0, ComplexNumber(1e-12, 0.0));
    s.qft(0, 1);

    State s2(2, 0);
    s2.set_qft_mode(State::QftMode::Direct);
    s2.set_amplitude(0, ComplexNumber(1e-12, 0.0));
    s2.iqft(0, 1);
}

void test_measure_branches() {
    State s(1, 1);
    s.set_basis_state(1, ONE_COMPLEX);
    s.measure_with_rng(0, 0, 0.0);
    s.measure(0, 0);

    State no_cbits(1, 0);
    no_cbits.set_basis_state(1, ONE_COMPLEX);
    std::vector<int> out;
    std::vector<double> rv(1, 0.8);
    no_cbits.measure_all_with_rng(rv, out);

    State tiny(1, 1);
    tiny.set_amplitude(0, ComplexNumber(0.0, 0.0));
    tiny.measure_with_rng(0, 0, 0.0);
}

void test_cli_command_parsing()
{
    std::vector<std::string> tokens = cli::parse_command("  INIT  3   2 ");
    if (tokens.size() != 3 || tokens[0] != "INIT" || tokens[1] != "3" || tokens[2] != "2") {
        throw std::runtime_error("CLI parse_command failed to tokenize INIT");
    }

    std::vector<std::string> empty = cli::parse_command("   ");
    if (!empty.empty()) {
        throw std::runtime_error("CLI parse_command should return empty for whitespace");
    }

    std::vector<std::string> angle_tokens = {"RX", "0", "PI/2"};
    double theta = cli::get_angle_arg_required(angle_tokens, 2, "RX");
    assert_double_close(theta, std::acos(-1.0) / 2.0, 1e-9, "CLI PI/2 parsing");

    std::vector<std::string> deg_tokens = {"RY", "0", "90", "DEG"};
    double theta_deg = cli::get_angle_arg_required(deg_tokens, 2, "RY");
    assert_double_close(theta_deg, std::acos(-1.0) / 2.0, 1e-9, "CLI DEG parsing");

    std::vector<std::string> missing = {"X"};
    int arg = cli::get_arg(missing, 1, "X");
    if (arg != -1) {
        throw std::runtime_error("CLI get_arg should return -1 on missing argument");
    }
}

void test_latin_square_validations() {
    const int row0[3] = {0, 1, 2};
    Bitstring assignment = 0;
    assignment |= (1ULL << 0);  // cell 0 = 1
    assignment |= (2ULL << 2);  // cell 1 = 2
    assignment |= (0ULL << 4);  // cell 2 = 0
    assignment |= (2ULL << 6);  // cell 3 = 2
    assignment |= (0ULL << 8);  // cell 4 = 0
    assignment |= (1ULL << 10); // cell 5 = 1

    if (!is_valid_latin3(assignment, row0)) {
        throw std::runtime_error("Expected assignment to be valid latin square");
    }

    if (is_valid_latin3_fixedrow(0)) {
        throw std::runtime_error("Expected zero assignment to be invalid latin square");
    }
}

void test_latin_square_demo_forced_measure() {
    const int row0[3] = {0, 1, 2};
    Bitstring assignment = 0;
    assignment |= (1ULL << 0);
    assignment |= (2ULL << 2);
    assignment |= (0ULL << 4);
    assignment |= (2ULL << 6);
    assignment |= (0ULL << 8);
    assignment |= (1ULL << 10);

    {
        ScopedEnv env("QSIM_LATIN_FORCE_MEASURED", std::to_string(assignment));
        run_latin3_grover_demo_row0(row0, 0);
    }
    {
        ScopedEnv env("QSIM_LATIN_FORCE_MEASURED", "0");
        run_latin3_grover_demo_row0(row0, 0);
    }
    {
        ScopedEnv env("QSIM_LATIN_FORCE_M", "0");
        run_latin3_grover_demo_row0(row0, 0);
    }
}

void test_latin_square_invalid_row0() {
    const int bad_range[3] = {0, 1, 3};
    const int bad_perm[3] = {0, 0, 1};
    assert_exits_with_failure([&]() { run_latin3_count_row0(bad_range); });
    assert_exits_with_failure([&]() { run_latin3_count_row0(bad_perm); });
}

void test_latin_square_count_print() {
    const int row0[3] = {0, 1, 2};
    run_latin3_count_row0(row0);
    {
        ScopedEnv env("QSIM_LATIN_MAX_PRINT", "1");
        run_latin3_print_all_row0(row0);
    }
    run_latin3_grover_demo_row0(row0, 1);
    run_latin3_grover_demo(0);
}

void test_shor_api_paths() {
    ShorResult bad = run_shor_quantum_part(1, 2);
    if (bad.ok || bad.error.empty()) {
        throw std::runtime_error("Expected shor_api error for N<2");
    }

    Bitstring big = (1ULL << 32);
    ShorResult too_big = run_shor_quantum_part(big, 2);
    if (too_big.ok || too_big.error.empty()) {
        throw std::runtime_error("Expected shor_api error for large N");
    }

    ShorResult ok = run_shor_quantum_part(15, 7);
    if (!ok.ok) {
        throw std::runtime_error("Expected shor_api to succeed for N=15");
    }
}

void test_shor_classical_estimate_order()
{
    Bitstring r = estimate_order(64, 8, 2, 15);
    if (r != 4) {
        throw std::runtime_error("Expected estimate_order to return r=4");
    }

    Bitstring r_none = estimate_order(0, 8, 2, 15);
    if (r_none != 0) {
        throw std::runtime_error("Expected estimate_order to return 0 when no r is valid");
    }
}

void test_shor_quantum_free_function_paths()
{
    ShorQuantumResult too_big = run_shor_algorithm_quantum_part(1ULL << 32, 2);
    if (too_big.ok || too_big.error.empty()) {
        throw std::runtime_error("Expected control-register size error");
    }

    ShorQuantumResult max_qubits = run_shor_algorithm_quantum_part(257, 2);
    if (max_qubits.ok || max_qubits.error.empty()) {
        throw std::runtime_error("Expected max-qubits error");
    }

    ShorQuantumResult ok = run_shor_algorithm_quantum_part(15, 2);
    if (!ok.ok || ok.n_c == 0) {
        throw std::runtime_error("Expected shor_quantum free function to succeed");
    }
}

void test_shor_demo_branches() {
    ScopedEnv env_attempts("QSIM_SHOR_MAX_ATTEMPTS", "1");
    run_shor_demo(1);
    run_shor_demo(10);

    {
        ScopedEnv env_attempts2("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env("QSIM_SHOR_FORCE_A", "5");
        run_shor_demo(15);
    }

    run_shor_demo(1ULL << 20);
    {
        ScopedEnv env_attempts2b("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        run_shor_demo((1ULL << 20) + 1); // odd, triggers max-qubits guard
    }
    {
        ScopedEnv env_attempts2c("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        run_shor_demo((1ULL << 20) + 1); // env_x path to guard
    }

    {
        ScopedEnv env_attempts3("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts3b("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts3c("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts4("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "0");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts5("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "1");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "1");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts5b("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "1");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "2");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts6("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "14");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "2");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts7("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "6");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "6");
        run_shor_demo(9);
    }

    {
        ScopedEnv env_attempts8("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "4");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts9("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "4");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts10("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "0");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "0");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts11("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        run_shor_demo((1ULL << 32) + 1);
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

#ifndef ALL_TESTS
int main()
{
  return run_unit_tests() == 0 ? 0 : 1;
}
#endif
