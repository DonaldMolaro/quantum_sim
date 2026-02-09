#define ALL_TESTS

#include "unit_tests.cc"
#include "grover_test.cc"
#include "grover_bench.cc"

int main()
{
  int failures = 0;
  failures += run_unit_tests();
  failures += run_grover_tests();

  if (const char* env = std::getenv("QSIM_GROVER_BENCH")) {
    if (env[0] != '\0' && env[0] != '0') {
      (void)run_grover_bench();
    }
  }

  if (failures != 0) {
    return 1;
  }
  return 0;
}
