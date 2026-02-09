#pragma once

#include "state.hh"
#include <string>
#include <vector>

enum class DeutschJozsaOracle {
  ConstantZero,
  ConstantOne,
  BalancedXor0,
  BalancedParity
};

struct DeutschJozsaResult {
  bool ok = false;
  bool is_constant = false;
  std::string error;
  std::vector<int> input_measurement;
};

DeutschJozsaResult run_deutsch_jozsa(int n_inputs, DeutschJozsaOracle oracle);

const char* deutsch_jozsa_oracle_name(DeutschJozsaOracle oracle);

bool parse_deutsch_jozsa_oracle(const std::string& token, DeutschJozsaOracle& out);
