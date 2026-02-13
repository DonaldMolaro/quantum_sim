# Define variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -I.
LDFLAGS =
CXXDEPFLAGS = -MMD -MP

TARGETS = quantum_sim all_tests
LIB_NAME = libquantum_sim.a
LIB_SOURCES = state.cc display.cc swap.cc qft.cc modular_exp.cc \
	logging.cc \
	math/bit_ops.cc math/mod_arith.cc \
	algorithms/grover.cc algorithms/deutsch_jozsa.cc algorithms/bernstein_vazirani.cc \
	algorithms/qubo.cc algorithms/vqa_qaoa.cc algorithms/qaoa.cc algorithms/vqe.cc algorithms/anneal.cc \
	algorithms/shor_classical.cc algorithms/shor_quantum.cc \
	algorithms/latin_square.cc algorithms/qrng.cc algorithms/api/grover_api.cc \
	algorithms/api/shor_api.cc demos/latin_demo.cc demos/shor_demo.cc \
	demos/grover_demo.cc demos/deutsch_jozsa_demo.cc demos/bernstein_vazirani_demo.cc \
	demos/qubo_demo.cc demos/vqa_demo.cc demos/qaoa_demo.cc demos/vqe_demo.cc demos/anneal_demo.cc
LIB_OBJECTS = $(LIB_SOURCES:.cc=.o)
CLI_SOURCES = cli/commands.cc cli/shell.cc cli/main.cc
DRIVER_SOURCES = $(CLI_SOURCES)
DRIVER_OBJECTS = $(DRIVER_SOURCES:.cc=.o)
DEPS    = $(LIB_SOURCES:.cc=.d) $(DRIVER_SOURCES:.cc=.d) tests/all_tests.d

ALL_TESTS_SRC = tests/all_tests.cc
TEST_EXTRA_OBJECTS = cli/commands.o

# Default target: builds the executable
all: $(TARGETS)

# Rule to link the object files into the final executable

$(LIB_NAME): $(LIB_OBJECTS)
	ar rcs $(LIB_NAME) $(LIB_OBJECTS)

quantum_sim : $(LIB_NAME) $(DRIVER_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(DRIVER_OBJECTS) $(LIB_NAME) -o quantum_sim

all_tests : tests/all_tests.o $(LIB_NAME) $(TEST_EXTRA_OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) tests/all_tests.o $(LIB_NAME) $(TEST_EXTRA_OBJECTS) -o all_tests

# Rule to compile .cc files into .o files (using implicit rule)
.cc.o:
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -c $< -o $@

.cc.d:
	$(CXX) $(CXXDEPFLAGS) $< -o $@

# Cleanup target
clean:
	rm -f $(LIB_OBJECTS) $(DRIVER_OBJECTS) $(DEPS) $(TARGETS) $(LIB_NAME) tests/all_tests.o all_tests unit_tests grover_test grover_bench tests/unit_tests.o tests/grover_test.o tests/grover_bench.o
	rm -f *.o *.d
	find . -name "*.gcda" -delete
	find . -name "*.gcno" -delete
	find . -name "*.gcov" -delete

-include $(DEPS)

.PHONY: clean test

ifeq ($(COVERAGE),1)
CXXFLAGS += -O0 -g --coverage
LDFLAGS += --coverage
endif

test:
	$(MAKE) clean
	$(MAKE) COVERAGE=1 all_tests
	find . -name "*.gcda" -delete
	find . -name "*.gcov" -delete
	QSIM_TEST_VERBOSE=1 ./all_tests
	./scripts/coverage_check.sh

.PHONY: test-slow
test-slow:
	$(MAKE) clean
	$(MAKE) COVERAGE=1 all_tests
	find . -name "*.gcda" -delete
	find . -name "*.gcov" -delete
	QSIM_TEST_VERBOSE=1 QSIM_SLOW_TESTS=1 ./all_tests
	QSIM_SLOW_TESTS=1 ./scripts/coverage_check.sh

.PHONY: test-demo
test-demo:
	$(MAKE) clean
	$(MAKE) COVERAGE=1 all_tests
	find . -name "*.gcda" -delete
	find . -name "*.gcov" -delete
	QSIM_TEST_VERBOSE=1 QSIM_SLOW_TESTS=1 QSIM_DEMO_TESTS=1 ./all_tests
	QSIM_SLOW_TESTS=1 QSIM_DEMO_TESTS=1 ./scripts/coverage_check.sh

.PHONY: coverage
coverage:
	$(MAKE) COVERAGE=1 all_tests
	QSIM_TEST_VERBOSE=1 ./all_tests
	./scripts/coverage_check.sh

.PHONY: install
install: $(LIB_NAME)
	mkdir -p dist/lib dist/include/algorithms/api dist/include/math dist/include/demos
	cp $(LIB_NAME) dist/lib/
	cp include/quantum_sim.hh dist/include/
	cp state.hh dist/include/
	cp logging.hh dist/include/
	cp modular_exp.hh dist/include/
	cp algorithms/latin_square.hh dist/include/algorithms/
	cp algorithms/deutsch_jozsa.hh dist/include/algorithms/
	cp algorithms/bernstein_vazirani.hh dist/include/algorithms/
	cp algorithms/qubo.hh dist/include/algorithms/
	cp algorithms/vqa_qaoa.hh dist/include/algorithms/
	cp algorithms/qaoa.hh dist/include/algorithms/
	cp algorithms/vqe.hh dist/include/algorithms/
	cp algorithms/anneal.hh dist/include/algorithms/
	cp algorithms/qrng.hh dist/include/algorithms/
	cp algorithms/api/grover_api.hh dist/include/algorithms/api/
	cp algorithms/api/shor_api.hh dist/include/algorithms/api/
	cp demos/latin_demo.hh dist/include/demos/
	cp demos/shor_demo.hh dist/include/demos/
	cp demos/grover_demo.hh dist/include/demos/
	cp demos/deutsch_jozsa_demo.hh dist/include/demos/
	cp demos/bernstein_vazirani_demo.hh dist/include/demos/
	cp demos/qubo_demo.hh dist/include/demos/
	cp demos/vqa_demo.hh dist/include/demos/
	cp demos/qaoa_demo.hh dist/include/demos/
	cp demos/vqe_demo.hh dist/include/demos/
	cp demos/anneal_demo.hh dist/include/demos/
	cp math/bit_ops.hh dist/include/math/
	cp math/mod_arith.hh dist/include/math/
