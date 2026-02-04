# Define variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
CXXDEPFLAGS = -MMD -MP

TARGETS = quantum_sim unit_tests grover_test
SOURCES = state.cc display.cc shell.cc grover.cc grover_search.cc \
	swap.cc qft.cc modular_exp.cc helper.cc shor.cc
OBJECTS = $(SOURCES:.cc=.o)
DEPS    = $(SOURCES:.cc=.d) main.d unit_tests.d grover_test.d

# Default target: builds the executable
all: $(TARGETS)

# Rule to link the object files into the final executable

quantum_sim : main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) main.o $(OBJECTS) -o quantum_sim

unit_tests : unit_tests.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) unit_tests.o $(OBJECTS) -o unit_tests

grover_test : grover_test.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) grover_test.o $(OBJECTS) -o grover_test

# Rule to compile .cc files into .o files (using implicit rule)
# This rule applies to both main.cc and state.cc
.cc.o:
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -c $< -o $@

.cc.d:
	$(CXX) $(CXXDEPFLAGS) $< -o $@

# Cleanup target
clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGETS) main.o unit_tests.o grover_test.o

-include $(DEPS)

.PHONY: clean test

test: grover_test
	./grover_test
