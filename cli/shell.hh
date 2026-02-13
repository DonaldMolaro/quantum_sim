/*
 * For learning things which have some "weirdness" involved I like
 * command line interactive systems. This is a nice little shell wrapped
 * around the quantum simulator, with it you can poke and prod it and see
 * how the underlying machine works. In theroy you could write a program in
 * the shell.
 */
#pragma once
#include "state.hh"
#include <complex>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

class QuantumShell
{
private:
  // The core quantum state object we manipulate.
  std::unique_ptr<State> state;
  void handle_command(const std::vector<std::string>& tokens);
  bool handle_setup_commands(const std::vector<std::string>& tokens, const std::string& cmd);
  bool handle_algorithm_commands(const std::vector<std::string>& tokens, const std::string& cmd);
  bool handle_single_qubit_commands(const std::vector<std::string>& tokens, const std::string& cmd);
  bool handle_parametric_gate_commands(const std::vector<std::string>& tokens, const std::string& cmd);
  bool handle_multi_qubit_commands(const std::vector<std::string>& tokens, const std::string& cmd);
  bool handle_measurement_and_display_commands(const std::vector<std::string>& tokens, const std::string& cmd);
public:
  void print_help();
  void run();
};
