#include <vector>
#include <complex>
#include <iostream>
#include <cmath>
#include <cassert>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include "state.hh"

const ComplexNumber I(0.0, 1.0); // Imaginary unit i
const double EPSILON = 1e-9; // Tolerance for complex number comparison

// --- Utility Functions ---

void run_test(const std::string& name, std::function<void()> test_func) {
    try {
        test_func();
        std::cout << "[PASS] " << name << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] " << name << ": " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "[FAIL] " << name << ": Unknown error" << std::endl;
    }
}

// Custom assertion for complex numbers
void assert_complex_equal(const ComplexNumber& expected, const ComplexNumber& actual, const std::string& message) {
    if (std::abs(expected - actual) > EPSILON) {
        throw std::runtime_error("Assertion failed: " + message + 
                                 " Expected: " + std::to_string(expected.real()) + 
                                 " Actual: " + std::to_string(actual.real()));
    }
}


// --- Test Cases (r=2, Phase = I) ---

const int CONTROL_Q = 1;
const int TARGET_Q = 0;
const int R_VALUE = 2; // Produces phase factor I

// Test Case 1: Control=0, Target=0 (State |00>, Bitstring 0)
void test_crr_00_no_change() {
    State s(2);
    s.set_basis_state(0b00, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE); 
    assert_complex_equal(1.0, s.get_amplitude(0b00), "00 state should be unchanged.");
}

// Test Case 2: Control=0, Target=1 (State |01>, Bitstring 1)
void test_crr_01_no_change() {
    State s(2);
    s.set_basis_state(0b01, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE);
    // Expect amplitude 1.0 (control bit 1 is 0)
    assert_complex_equal(1.0, s.get_amplitude(0b01), "01 state should be unchanged.");
}

// Test Case 3: Control=1, Target=0 (State |10>, Bitstring 2)
void test_crr_10_no_change() {
    State s(2);
    s.set_basis_state(0b10, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE);
    // Expect amplitude 1.0 (target bit 0 is 0)
    assert_complex_equal(1.0, s.get_amplitude(0b10), "10 state should be unchanged.");
}

// Test Case 4: Control=1, Target=1 (State |11>, Bitstring 3)
void test_crr_11_applies_phase_I() {
    State s(2);
    s.set_basis_state(0b11, 1.0);
    s.controlled_Rr(CONTROL_Q, TARGET_Q, R_VALUE); 
    // Expect amplitude I (phase rotation by pi/2 for r=2)
    assert_complex_equal(I, s.get_amplitude(0b11), "11 state should be multiplied by I.");
}

// Test Case 5: Superposition Test (Linearity)
void test_crr_superposition() {
    State s(2);
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
    State s(2);
    s.set_basis_state(0b00, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE); 
    assert_complex_equal(1.0, s.get_amplitude(0b00), "00 state should be unchanged.");
}

// Test Case 2: Control=0, Target=1 (State |01>, Bitstring 1)
void test_crrdag_01_no_change() {
    State s(2);
    s.set_basis_state(0b01, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE);
    assert_complex_equal(1.0, s.get_amplitude(0b01), "01 state should be unchanged.");
}

// Test Case 3: Control=1, Target=0 (State |10>, Bitstring 2)
void test_crrdag_10_no_change() {
    State s(2);
    s.set_basis_state(0b10, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE);
    assert_complex_equal(1.0, s.get_amplitude(0b10), "10 state should be unchanged.");
}

// Test Case 4: Control=1, Target=1 (State |11>, Bitstring 3) - Applies phase -I
void test_crrdag_11_applies_phase_minus_I() {
    State s(2);
    s.set_basis_state(0b11, 1.0);
    s.controlled_Rr_dag(CONTROL_Q, TARGET_Q, R_VALUE); 
    // Expected result: -I
    assert_complex_equal(-I, s.get_amplitude(0b11), "11 state should be multiplied by -I.");
}

// Test Case 5: Verification that R_r * R_r_dag = I
void test_crrdag_inversion() {
    State s(2);
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
    State s(2);
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
    State s(2);
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
    State s(2);
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
    State s(2);
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
    State s(2);
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
    State s(2);
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

void assert_amplitude_magnitude(const State& s, Bitstring b, double expected_mag, const std::string& msg) {
    double actual_mag = std::abs(s.get_amplitude(b));
    if (std::abs(actual_mag - expected_mag) > EPSILON) {
        throw std::runtime_error("Amplitude check failed for " + msg + " (Bitstring " + std::to_string(b) + 
                                 "). Expected magnitude: " + std::to_string(expected_mag) + 
                                 ", Actual: " + std::to_string(actual_mag));
    }
}

void test_shor_quantum_part_n15_a7_r4() {
    // Test Case: Factor N=15, choosing base a=7. Known period r=4.
    // Precision: m=5 qubits (2^5 = 32 states).
    // Expected output indices (x/2^m) = {0/4, 1/4, 2/4, 3/4} in the control register.
    // Corresponding control register values: {0, 8, 16, 24}.
    // Total qubits = 5 (control) + 4 (target) = 9
    State s(5,5);
    s.display();
    const Bitstring N = 15;
    const Bitstring a = 7;
    const int R = 4;
    const double expected_mag = 1.0 / std::sqrt(static_cast<double>(R));
    
    s.run_shor_algorithm_quantum_part(N, a);

    // Qubits 0-3 (target register) are fixed at |0001> (index 1).
    // Qubits 4-8 (control register) hold the period information.
    
    // Test expected states (x/2^m = s/r)
    // s=0 -> x=0 -> Bitstring index 1 (000000001)
    assert_amplitude_magnitude(s, 1ULL, expected_mag, "s=0/r=4 state (|00...01>)");
    
    // s=1 -> x=8 -> Bitstring index 1 + (8 << 4) = 129
    assert_amplitude_magnitude(s, 1ULL | (8ULL << 4), expected_mag, "s=1/r=4 state");
    
    // s=2 -> x=16 -> Bitstring index 1 + (16 << 4) = 257
    assert_amplitude_magnitude(s, 1ULL | (16ULL << 4), expected_mag, "s=2/r=4 state");
    
    // s=3 -> x=24 -> Bitstring index 1 + (24 << 4) = 385
    assert_amplitude_magnitude(s, 1ULL | (24ULL << 4), expected_mag, "s=3/r=4 state");

    // Test a state expected to have near zero amplitude (e.g., index 2, which is |00...02>)
    assert_amplitude_magnitude(s, 2ULL, 0.0, "state |00...02>");
    
    // Test another zero state (e.g., control register 1, target register 1)
    assert_amplitude_magnitude(s, 1ULL | (1ULL << 4), 0.0, "state x=1/2^m");
}

void main_shor_tests() {
    std::cout << "Testing State::run_shor_algorithm_quantum_part (Order Finding):\n";
    std::cout << "------------------------------------------------------------------\n";
    
    run_test("Test 1: N=15, a=7, r=4 (Verification of QPE output concentration)", test_shor_quantum_part_n15_a7_r4);
    
    std::cout << "------------------------------------------------------------------\n";
}


extern Bitstring modular_exponentiation(Bitstring a, Bitstring power, Bitstring N);

void assert_equal(Bitstring actual, Bitstring expected, const std::string& message) {
    if (actual != expected) {
        throw std::runtime_error("Test failed for " + message + 
                                 ". Expected: " + std::to_string(expected) + 
                                 ", Actual: " + std::to_string(actual));
    }
}


void test_mod_exp_zero_power() {
    // Test case: a^0 mod N = 1
    Bitstring result = modular_exponentiation(42ULL, 0ULL, 100ULL);
    assert_equal(result, 1ULL, "Zero Power Test (42^0 mod 100)");
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

void test_mod_exp_large_power() {
    // Test 2^63 mod 100 = 8
    // Uses a power close to the maximum Bitstring value (if Bitstring is 64-bit).
    Bitstring result = modular_exponentiation(2ULL, 63ULL, 100ULL);
    assert_equal(result, 8ULL, "Large Power Test (2^63 mod 100)");
}

int main_mod_exp() {
     std::cout << "Testing modular_exponentiation utility:" << std::endl;
     run_test("Test 1: Zero Power", test_mod_exp_zero_power);
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

void assert_amplitude_match(const State& s, Bitstring b, ComplexNumber expected_a, const std::string& msg) {
    ComplexNumber actual_a = s.get_amplitude(b);
    if (std::abs(actual_a - expected_a) > EPSILON) {
        throw std::runtime_error("Amplitude check failed for " + msg + " (Bitstring " + std::to_string(b) + 
                                 "). Expected: " + std::to_string(expected_a.real()) + 
                                 ", Actual: " + std::to_string(actual_a.real()));
    }
}

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

    State s(START_B);

    
    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 1);


    
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

    State s(START_B); 
    
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

    State s(0);
    s.set_amplitude(B_OFF, AMP);
    s.set_amplitude(B_ON_IN, AMP);
    
    s.controlled_modular_exponentiation(CONTROL_QUBIT, TARGET_START, TARGET_END, A, N, 1);

    assert_amplitude_match(s, B_OFF_OUT, AMP, "Superposition: Control OFF component (B=1) maintained magnitude");
    assert_amplitude_match(s, B_ON_OUT, AMP, "Superposition: Control ON component (B=71) maintained magnitude");
    assert_amplitude_match(s, B_ON_IN, 0.0, "Superposition: Initial ON input state must vanish");
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
    
    std::cout << "------------------------------------------------------------------\n";
}




int main()
{
  main_test_controlled_Rr();
  main_test_controlled_Rr_dag();
  main_test_qft();
  main_mod_exp();
  main_all_cme_tests();
  main_shor_tests();
}
