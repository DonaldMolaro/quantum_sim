#pragma once
#include "state.hh"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

namespace test_helpers {

constexpr double kEps = 1e-9;

inline void assert_complex_equal(const ComplexNumber& expected,
                                 const ComplexNumber& actual,
                                 const std::string& message,
                                 double eps = kEps) {
    if (std::abs(expected - actual) > eps) {
        throw std::runtime_error("Assertion failed: " + message +
                                 " Expected: " + std::to_string(expected.real()) +
                                 " Actual: " + std::to_string(actual.real()));
    }
}

inline void assert_complex_close(const ComplexNumber& expected,
                                 const ComplexNumber& actual,
                                 double tol,
                                 const std::string& message) {
    if (std::abs(expected - actual) > tol) {
        throw std::runtime_error("Assertion failed: " + message +
                                 " Expected: " + std::to_string(expected.real()) +
                                 " Actual: " + std::to_string(actual.real()));
    }
}

inline void assert_equal(Bitstring actual, Bitstring expected, const std::string& message) {
    if (actual != expected) {
        throw std::runtime_error("Test failed for " + message +
                                 ". Expected: " + std::to_string(expected) +
                                 ", Actual: " + std::to_string(actual));
    }
}

inline void assert_amplitude_match(const State& s,
                                  Bitstring b,
                                  ComplexNumber expected_a,
                                  const std::string& msg,
                                  double eps = kEps) {
    ComplexNumber actual_a = s.get_amplitude(b);
    if (std::abs(actual_a - expected_a) > eps) {
        throw std::runtime_error("Amplitude check failed for " + msg + " (Bitstring " + std::to_string(b) +
                                 "). Expected: " + std::to_string(expected_a.real()) +
                                 ", Actual: " + std::to_string(actual_a.real()));
    }
}

inline void assert_amplitude_magnitude(const State& s,
                                       Bitstring b,
                                       double expected_mag,
                                       const std::string& msg,
                                       double eps = kEps) {
    double actual_mag = std::abs(s.get_amplitude(b));
    if (std::abs(actual_mag - expected_mag) > eps) {
        throw std::runtime_error("Amplitude check failed for " + msg + " (Bitstring " + std::to_string(b) +
                                 "). Expected magnitude: " + std::to_string(expected_mag) +
                                 ", Actual: " + std::to_string(actual_mag));
    }
}

} // namespace test_helpers
