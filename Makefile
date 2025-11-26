# Define variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGET = quantum_sim
SOURCES = main.cc state.cc display.cc shell.cc
OBJECTS = $(SOURCES:.cc=.o)

# Default target: builds the executable
all: $(TARGET)

# Rule to link the object files into the final executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

# Rule to compile .cc files into .o files (using implicit rule)
# This rule applies to both main.cc and state.cc
.cc.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Cleanup target
clean:
	rm -f $(OBJECTS) $(TARGET)
