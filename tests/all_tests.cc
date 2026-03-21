#include <cstdlib>
#include "tests/test_harness.hh"

int run_unit_tests();
int run_grover_tests();
int run_grover_bench();
int run_new_feature_tests();
int run_edge_case_tests();

int main()
{
  const test_harness::Suite suites[] = {
    {"unit", run_unit_tests},
    {"grover", run_grover_tests},
    {"new-features", run_new_feature_tests},
    {"edge-cases", run_edge_case_tests},
  };
  int failures = test_harness::run_suites(suites, sizeof(suites) / sizeof(suites[0]));

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
