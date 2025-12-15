# Define variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGETS = quantum_sim unit_tests
SOURCES = state.cc display.cc shell.cc grover.cc grover_search.cc \
	swap.cc qft.cc modular_exp.cc helper.cc shor.cc
OBJECTS = $(SOURCES:.cc=.o)

# Default target: builds the executable
all: $(TARGETS)

# Rule to link the object files into the final executable

quantum_sim : main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) main.o $(OBJECTS) -o quantum_sim

unit_tests : unit_tests.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) unit_tests.o $(OBJECTS) -o unit_tests

# Rule to compile .cc files into .o files (using implicit rule)
# This rule applies to both main.cc and state.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Cleanup target
clean:
	rm -f $(OBJECTS) $(TARGET)
