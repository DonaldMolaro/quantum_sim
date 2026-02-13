#pragma once
#include "state.hh"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <functional>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

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

inline void assert_double_close(double actual, double expected, double tol, const std::string& message)
{
    if (std::abs(actual - expected) > tol) {
        throw std::runtime_error("Assertion failed: " + message +
                                 " Expected: " + std::to_string(expected) +
                                 " Actual: " + std::to_string(actual));
    }
}

inline void assert_exits_with_failure(const std::function<void()>& fn)
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

} // namespace test_helpers
