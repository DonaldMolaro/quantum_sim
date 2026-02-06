/*
 * For learning things which have some "weirdness" involved I like
 * command line interactive systems. This is a nice little shell wrapped
 * around the quantum simulator, with it you can poke and prod it and see
 * how the underlying machine works. In theroy you could write a program in
 * the shell.
 */
#pragma once
#include <complex>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class QuantumShell
{
private:
  // The core quantum state object we manipulate
  State* state = nullptr; 

  // Helper function to tokenize input line
  std::vector<std::string> parse_command(const std::string& line);
  // Helper to get integer arguments safely
  int get_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
  // Helper to get floating-point arguments safely
  double get_double_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
  // Helper to parse angles (radians by default, optional degrees)
  double get_angle_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
  // Helper to parse required angle arguments
  double get_angle_arg_required(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
  void handle_command(const std::vector<std::string>& tokens);
public:
  // Destructor to clean up the dynamically allocated State object
  ~QuantumShell();
  void print_help();
  void run();
};
