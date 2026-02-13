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
public:
  void print_help();
  void run();
};
