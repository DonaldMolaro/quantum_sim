# Define variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -I.
CXXDEPFLAGS = -MMD -MP

TARGETS = quantum_sim unit_tests grover_test
SOURCES = state.cc display.cc shell.cc swap.cc qft.cc modular_exp.cc \
	math/bit_ops.cc math/mod_arith.cc \
	algorithms/grover.cc algorithms/grover_search.cc algorithms/shor.cc
OBJECTS = $(SOURCES:.cc=.o)
DEPS    = $(SOURCES:.cc=.d) main.d tests/unit_tests.d tests/grover_test.d

UNIT_TEST_SRC = tests/unit_tests.cc
GROVER_TEST_SRC = tests/grover_test.cc

# Default target: builds the executable
all: $(TARGETS)

# Rule to link the object files into the final executable

quantum_sim : main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) main.o $(OBJECTS) -o quantum_sim

unit_tests : tests/unit_tests.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) tests/unit_tests.o $(OBJECTS) -o unit_tests

grover_test : tests/grover_test.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) tests/grover_test.o $(OBJECTS) -o grover_test

# Rule to compile .cc files into .o files (using implicit rule)
.cc.o:
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -c $< -o $@

.cc.d:
	$(CXX) $(CXXDEPFLAGS) $< -o $@

# Cleanup target
clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGETS) main.o tests/unit_tests.o tests/grover_test.o

-include $(DEPS)

.PHONY: clean test

test: grover_test
	./grover_test
