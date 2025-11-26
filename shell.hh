#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <sstream>

class QuantumShell
{
private:
  // The core quantum state object we manipulate
  State* state = nullptr; 

  // Helper function to tokenize input line
  std::vector<std::string> parse_command(const std::string& line);
  // Helper to get integer arguments safely
  int get_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
  void handle_command(const std::vector<std::string>& tokens);
public:
  // Destructor to clean up the dynamically allocated State object
  ~QuantumShell();
  void print_help();
  void run();
};
