cmake_minimum_required(VERSION 3.16)
project(quantum_sim LANGUAGES CXX)

include(GNUInstallDirs)

option(QSIM_ENABLE_COVERAGE "Enable gcov/coverage instrumentation" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  add_compile_options(-Wall -Wextra -O2)
endif()

if(QSIM_ENABLE_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-O0 -g --coverage)
    add_link_options(--coverage)
  else()
    message(FATAL_ERROR "QSIM_ENABLE_COVERAGE is only supported with Clang/GCC")
  endif()
endif()

set(LIB_SOURCES
  state.cc
  display.cc
  swap.cc
  qft.cc
  modular_exp.cc
  logging.cc
  math/bit_ops.cc
  math/mod_arith.cc
  algorithms/grover.cc
  algorithms/deutsch_jozsa.cc
  algorithms/bernstein_vazirani.cc
  algorithms/qubo.cc
  algorithms/vqa_qaoa.cc
  algorithms/qaoa.cc
  algorithms/vqe.cc
  algorithms/anneal.cc
  algorithms/shor_classical.cc
  algorithms/shor_quantum.cc
  algorithms/latin_square.cc
  algorithms/qrng.cc
  algorithms/tsp.cc
  algorithms/quantum_counting.cc
  algorithms/simon.cc
  algorithms/api/grover_api.cc
  algorithms/api/shor_api.cc
)

set(DEMO_SOURCES
  demos/latin_demo.cc
  demos/shor_demo.cc
  demos/grover_demo.cc
  demos/deutsch_jozsa_demo.cc
  demos/bernstein_vazirani_demo.cc
  demos/qubo_demo.cc
  demos/vqa_demo.cc
  demos/qaoa_demo.cc
  demos/vqe_demo.cc
  demos/anneal_demo.cc
  demos/tsp_demo.cc
  demos/quantum_counting_demo.cc
  demos/simon_demo.cc
)

set(CLI_SOURCES
  cli/commands.cc
  cli/shell.cc
  cli/shell_setup.cc
  cli/shell_algorithms.cc
  cli/shell_gates.cc
)

set(TEST_SOURCES
  tests/unit_tests.cc
  tests/unit_qft_modexp.cc
  tests/unit_algorithms_and_cli.cc
  tests/unit_cli_and_shor.cc
  tests/grover_test.cc
  tests/grover_bench.cc
)

add_library(quantum_sim_static STATIC ${LIB_SOURCES})
target_include_directories(quantum_sim_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(quantum_sim_static PROPERTIES OUTPUT_NAME quantum_sim)

add_library(quantum_sim_demos STATIC ${DEMO_SOURCES})
target_include_directories(quantum_sim_demos PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(quantum_sim_demos PUBLIC quantum_sim_static)

add_library(quantum_sim_cli STATIC ${CLI_SOURCES})
target_include_directories(quantum_sim_cli PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(quantum_sim_cli PUBLIC quantum_sim_static quantum_sim_demos)

add_executable(quantum_sim cli/main.cc)
target_link_libraries(quantum_sim PRIVATE quantum_sim_cli quantum_sim_demos quantum_sim_static)

add_executable(all_tests tests/all_tests.cc ${TEST_SOURCES} cli/commands.cc)
target_link_libraries(all_tests PRIVATE quantum_sim_demos quantum_sim_static)

enable_testing()
add_test(NAME test-fast COMMAND all_tests)
set_tests_properties(test-fast PROPERTIES ENVIRONMENT "QSIM_TEST_VERBOSE=1")

add_test(NAME test-slow COMMAND all_tests)
set_tests_properties(test-slow PROPERTIES ENVIRONMENT "QSIM_TEST_VERBOSE=1;QSIM_SLOW_TESTS=1")

add_test(NAME test-demo COMMAND all_tests)
set_tests_properties(test-demo PROPERTIES ENVIRONMENT "QSIM_TEST_VERBOSE=1;QSIM_SLOW_TESTS=1;QSIM_DEMO_TESTS=1")

add_custom_target(test-examples
  COMMAND ${CMAKE_COMMAND} -E env PATH=$ENV{PATH} ${CMAKE_SOURCE_DIR}/examples/run_examples.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS quantum_sim
  COMMENT "Running example smoke scripts")

add_custom_target(test-labs
  COMMAND ${CMAKE_COMMAND} -E env PATH=$ENV{PATH} ${CMAKE_SOURCE_DIR}/examples/check_labs.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS quantum_sim
  COMMENT "Running guided lab checks")

add_custom_target(coverage-check
  COMMAND ${CMAKE_COMMAND} -E env SOURCE_ROOT=${CMAKE_SOURCE_DIR} QSIM_TEST_VERBOSE=1 ./all_tests
  COMMAND ${CMAKE_COMMAND} -E env SOURCE_ROOT=${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/scripts/coverage_check.sh
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  DEPENDS all_tests
  COMMENT "Running coverage gate")

install(TARGETS quantum_sim_static
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES include/quantum_sim.hh include/quantum_sim_demos.hh state.hh logging.hh modular_exp.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES
  algorithms/latin_square.hh
  algorithms/deutsch_jozsa.hh
  algorithms/bernstein_vazirani.hh
  algorithms/qubo.hh
  algorithms/vqa_qaoa.hh
  algorithms/qaoa.hh
  algorithms/vqe.hh
  algorithms/anneal.hh
  algorithms/qrng.hh
  algorithms/tsp.hh
  algorithms/quantum_counting.hh
  algorithms/simon.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/algorithms)
install(FILES algorithms/api/grover_api.hh algorithms/api/shor_api.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/algorithms/api)
install(FILES
  demos/latin_demo.hh
  demos/shor_demo.hh
  demos/grover_demo.hh
  demos/deutsch_jozsa_demo.hh
  demos/bernstein_vazirani_demo.hh
  demos/qubo_demo.hh
  demos/vqa_demo.hh
  demos/qaoa_demo.hh
  demos/vqe_demo.hh
  demos/anneal_demo.hh
  demos/tsp_demo.hh
  demos/quantum_counting_demo.hh
  demos/simon_demo.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/demos)
install(FILES math/bit_ops.hh math/mod_arith.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/math)
