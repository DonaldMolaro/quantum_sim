#pragma once

#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace test_harness {

typedef std::function<void()> TestFunc;

struct TestCase {
    const char* name;
    void (*fn)();
};

inline int& failure_count_storage()
{
    static int failures = 0;
    return failures;
}

inline void reset_failures()
{
    failure_count_storage() = 0;
}

inline int failure_count()
{
    return failure_count_storage();
}

inline int run_test(const std::string& name, const TestFunc& test_func)
{
    try {
        test_func();
        std::cout << "[PASS] " << name << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] " << name << ": " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "[FAIL] " << name << ": Unknown error" << std::endl;
    }
    ++failure_count_storage();
    return 1;
}

inline bool env_flag(const char* name)
{
    const char* env = std::getenv(name);
    return env && env[0] != '\0' && env[0] != '0';
}

inline bool env_unsigned(const char* name, unsigned int& out)
{
    const char* env = std::getenv(name);
    if (!env || env[0] == '\0') {
        return false;
    }
    out = static_cast<unsigned int>(std::strtoul(env, nullptr, 10));
    return true;
}

inline void print_skip(const std::string& label, const std::string& hint)
{
    std::cout << label << " skipped. " << hint << std::endl;
}

struct Suite {
    const char* name;
    int (*run)();
};

inline int run_suites(const Suite* suites, size_t count)
{
    int failures = 0;
    for (size_t i = 0; i < count; ++i) {
        failures += suites[i].run();
    }
    return failures;
}

inline void run_cases(const TestCase* cases, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        run_test(cases[i].name, cases[i].fn);
    }
}

} // namespace test_harness
