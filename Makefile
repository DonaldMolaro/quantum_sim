# Define variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -I.
CXXDEPFLAGS = -MMD -MP

TARGETS = quantum_sim unit_tests grover_test
LIB_NAME = libquantum_sim.a
LIB_SOURCES = state.cc display.cc swap.cc qft.cc modular_exp.cc \
	math/bit_ops.cc math/mod_arith.cc \
	algorithms/grover.cc algorithms/grover_search.cc algorithms/shor.cc \
	algorithms/latin_square.cc algorithms/api/grover_api.cc \
	algorithms/api/shor_api.cc
LIB_OBJECTS = $(LIB_SOURCES:.cc=.o)
DRIVER_SOURCES = cli/shell.cc cli/main.cc
DRIVER_OBJECTS = $(DRIVER_SOURCES:.cc=.o)
DEPS    = $(LIB_SOURCES:.cc=.d) $(DRIVER_SOURCES:.cc=.d) tests/unit_tests.d tests/grover_test.d

UNIT_TEST_SRC = tests/unit_tests.cc
GROVER_TEST_SRC = tests/grover_test.cc

# Default target: builds the executable
all: $(TARGETS)

# Rule to link the object files into the final executable

$(LIB_NAME): $(LIB_OBJECTS)
	ar rcs $(LIB_NAME) $(LIB_OBJECTS)

quantum_sim : $(LIB_NAME) $(DRIVER_OBJECTS)
	$(CXX) $(CXXFLAGS) $(DRIVER_OBJECTS) $(LIB_NAME) -o quantum_sim

unit_tests : tests/unit_tests.o $(LIB_NAME)
	$(CXX) $(CXXFLAGS) tests/unit_tests.o $(LIB_NAME) -o unit_tests

grover_test : tests/grover_test.o $(LIB_NAME)
	$(CXX) $(CXXFLAGS) tests/grover_test.o $(LIB_NAME) -o grover_test

# Rule to compile .cc files into .o files (using implicit rule)
.cc.o:
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -c $< -o $@

.cc.d:
	$(CXX) $(CXXDEPFLAGS) $< -o $@

# Cleanup target
clean:
	rm -f $(LIB_OBJECTS) $(DRIVER_OBJECTS) $(DEPS) $(TARGETS) $(LIB_NAME) tests/unit_tests.o tests/grover_test.o

-include $(DEPS)

.PHONY: clean test

test: grover_test
	./grover_test

.PHONY: install
install: $(LIB_NAME)
	mkdir -p dist/lib dist/include/algorithms/api dist/include/math
	cp $(LIB_NAME) dist/lib/
	cp include/quantum_sim.hh dist/include/
	cp state.hh dist/include/
	cp algorithms/latin_square.hh dist/include/algorithms/
	cp algorithms/api/grover_api.hh dist/include/algorithms/api/
	cp algorithms/api/shor_api.hh dist/include/algorithms/api/
	cp math/bit_ops.hh dist/include/math/
	cp math/mod_arith.hh dist/include/math/
	cp math/qft_utils.hh dist/include/math/
	cp math/register_layout.hh dist/include/math/
