#include "state.hh"
#include "tests/helpers.hh"
#include "modular_exp.hh"
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void run_test(const std::string& name, std::function<void()> test_func);

using test_helpers::assert_complex_equal;
using test_helpers::assert_equal;
using test_helpers::assert_amplitude_match;
using test_helpers::assert_amplitude_magnitude;

const ComplexNumber I(0.0, 1.0);

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

