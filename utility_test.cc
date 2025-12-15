#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

using ComplexNumber = std::complex<double>; 
using Bitstring = unsigned long long; 
const double EPSILON = 1e-9; 

// --- Constants (matching the modular exponentiation context) ---
const int N_TARGET = 4; // Qubits 0-3 for target register (N=15)
const int CONTROL_Q6 = 6; // Example Control Qubit index (q6)
const int CONTROL_Q7 = 7; // Example Control Qubit index (q7)

// Each element of the state set is a (bitstring, amplitude) pair
using QubitAmplitudePair = std::pair<Bitstring, ComplexNumber>;
// The QuantumState is the set/vector of these pairs
using QuantumState = std::vector<QubitAmplitudePair>;

// --- Functions Under Test (Replicated from Correction) ---

// Helper function 1: Replaces the target register value in a full bitstring
Bitstring replace_target_value(Bitstring original, Bitstring new_target) {
    Bitstring target_mask = (1ULL << N_TARGET) - 1;
    return (original & ~target_mask) | new_target;
}

// Helper function 2: Finds an existing amplitude pair or adds a new one (amplitude 0.0)
QubitAmplitudePair* find_or_add(QuantumState& state, Bitstring b) {
    // Iterate through the vector to find the bitstring
    for (auto& pair : state) {
        if (pair.first == b) {
            return &pair;
        }
    }
    // Not found, add new entry with zero amplitude
    state.push_back({b, 0.0});
    return &state.back();
}


// --- Unit Test Helpers ---

void assert_equal(Bitstring actual, Bitstring expected, const std::string& message) {
    if (actual != expected) {
        throw std::runtime_error("Assertion failed for " + message + 
                                 ". Expected: " + std::to_string(expected) + 
                                 ", Actual: " + std::to_string(actual));
    }
}

void assert_amplitude_equal(ComplexNumber actual, ComplexNumber expected, const std::string& message) {
    if (std::abs(actual - expected) > EPSILON) {
        throw std::runtime_error("Amplitude assertion failed for " + message + 
                                 ". Expected: " + std::to_string(expected.real()) + 
                                 ", Actual: " + std::to_string(actual.real()));
    }
}

void run_test(const std::string& name, void (*test_func)()) {
    try {
        test_func();
        std::cout << "[PASS] " << name << std::endl;
    } catch (const std::runtime_error& e) {
        std::cerr << "[FAIL] " << name << ": " << e.what() << std::endl;
    }
}


// --- Unit Tests for replace_target_value ---

void test_replace_target_no_control() {
    // Initial state |0> (0) -> Target |15> (1111)
    Bitstring original = 0ULL;
    Bitstring expected = 15ULL; // 00001111
    Bitstring result = replace_target_value(original, 15ULL);
    assert_equal(result, expected, "Target replacement (0 -> 15)");
}

void test_replace_target_with_control() {
    // Initial state |1> + Control Q6 ON (65) -> Target |4> (0100)
    // 65 is 01000001. We want 01000100 (68).
    const Bitstring START_B = 1ULL | (1ULL << CONTROL_Q6); // 65
    const Bitstring EXPECTED_B = 4ULL | (1ULL << CONTROL_Q6); // 68
    
    Bitstring result = replace_target_value(START_B, 4ULL);
    assert_equal(result, EXPECTED_B, "Target replacement (1 -> 4) with Q6 control");
}

void test_replace_target_max_to_min() {
    // Initial state |15> + Control Q7 ON (143) -> Target |0> (0000)
    // 143 is 10001111. We want 10000000 (128).
    const Bitstring START_B = 15ULL | (1ULL << CONTROL_Q7); // 143
    const Bitstring EXPECTED_B = 0ULL | (1ULL << CONTROL_Q7); // 128
    
    Bitstring result = replace_target_value(START_B, 0ULL);
    assert_equal(result, EXPECTED_B, "Target replacement (15 -> 0) with Q7 control");
}


// --- Unit Tests for find_or_add ---

void test_find_or_add_new_bitstring() {
    QuantumState s = {{1ULL, 1.0}};
    size_t initial_size = s.size();

    // Search for a new bitstring (10ULL)
    QubitAmplitudePair* entry = find_or_add(s, 10ULL);

    // Assert size increased and initial amplitude is 0.0
    if (s.size() != initial_size + 1) {
        throw std::runtime_error("Size check failed for new bitstring.");
    }
    assert_equal(entry->first, 10ULL, "Check new bitstring key");
    assert_amplitude_equal(entry->second, 0.0, "Check initial amplitude (must be 0.0)");
}

void test_find_or_add_existing_bitstring() {
    ComplexNumber original_amp = {0.5, 0.0};
    QuantumState s = {{5ULL, original_amp}, {10ULL, 1.0}};
    size_t initial_size = s.size();

    // Search for an existing bitstring (5ULL)
    QubitAmplitudePair* entry = find_or_add(s, 5ULL);

    // Assert size remains the same and amplitude matches the original
    if (s.size() != initial_size) {
        throw std::runtime_error("Size check failed for existing bitstring.");
    }
    assert_amplitude_equal(entry->second, original_amp, "Check amplitude of existing bitstring (must match)");
}

void test_find_or_add_update_via_pointer() {
    QuantumState s = {{10ULL, 0.0}};
    
    // 1. Find the entry and update it
    QubitAmplitudePair* entry = find_or_add(s, 10ULL);
    entry->second = {0.5, 0.5};

    // 2. Find it again and check the updated value
    QubitAmplitudePair* updated_entry = find_or_add(s, 10ULL);

    assert_amplitude_equal(updated_entry->second, {0.5, 0.5}, "Check updated amplitude via pointer");
}


// --- Main Test Harness ---

int main() {
    std::cout << "Testing replace_target_value:\n";
    run_test("Test 1: Target (0 -> 15) no control", test_replace_target_no_control);
    run_test("Test 2: Target (1 -> 4) with Q6 control", test_replace_target_with_control);
    run_test("Test 3: Target (15 -> 0) with Q7 control", test_replace_target_max_to_min);

    std::cout << "\nTesting find_or_add (Vector State Management):\n";
    run_test("Test 4: Add new bitstring (check size and 0.0 init)", test_find_or_add_new_bitstring);
    run_test("Test 5: Find existing bitstring (check size and amplitude)", test_find_or_add_existing_bitstring);
    run_test("Test 6: Update amplitude via returned pointer", test_find_or_add_update_via_pointer);
    
    return 0;
}

