#include "cli/commands.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace cli {

std::vector<std::string> parse_command(const std::string& line)
{
  std::vector<std::string> tokens;
  std::stringstream ss(line);
  std::string token;
  while (ss >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

int get_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return -1;
  }
  try {
    return std::stoi(tokens[index]);
  } catch (const std::exception&) {
    std::cerr << "Error: Argument '" << tokens[index] << "' must be an integer.\n";
    return -1;
  }
}

double get_double_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
  try {
    return std::stod(tokens[index]);
  } catch (const std::exception&) {
    std::cerr << "Error: Argument '" << tokens[index] << "' must be a number.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
}

double get_angle_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::string token = tokens[index];
  bool is_deg = false;

  if (token.size() >= 3 && token.substr(token.size() - 3) == "DEG") {
    token = token.substr(0, token.size() - 3);
    is_deg = true;
  } else if (index + 1 < tokens.size() && tokens[index + 1] == "DEG") {
    is_deg = true;
  }

  std::string token_upper = token;
  std::transform(token_upper.begin(), token_upper.end(), token_upper.begin(), ::toupper);

  const double pi = std::acos(-1.0);
  double theta = std::numeric_limits<double>::quiet_NaN();

  if (token_upper == "PI") {
    theta = pi;
    is_deg = false;
  } else if (token_upper == "TAU") {
    theta = 2.0 * pi;
    is_deg = false;
  } else if (token_upper == "PI/2") {
    theta = pi / 2.0;
    is_deg = false;
  } else if (token_upper == "PI/4") {
    theta = pi / 4.0;
    is_deg = false;
  } else if (token_upper == "-PI") {
    theta = -pi;
    is_deg = false;
  } else if (token_upper == "-PI/2") {
    theta = -pi / 2.0;
    is_deg = false;
  } else if (token_upper == "-PI/4") {
    theta = -pi / 4.0;
    is_deg = false;
  } else {
    bool parsed_pi_expr = false;
    bool neg = false;
    std::string expr = token_upper;
    if (!expr.empty() && expr[0] == '-') {
      neg = true;
      expr = expr.substr(1);
    }

    size_t pi_pos = expr.find("PI");
    if (pi_pos != std::string::npos) {
      std::string coeff_str = expr.substr(0, pi_pos);
      if (!coeff_str.empty()) {
        try {
          double coeff = std::stod(coeff_str);
          theta = coeff * pi;
          parsed_pi_expr = true;
        } catch (const std::exception&) {
          parsed_pi_expr = false;
        }
      } else {
        theta = pi;
        parsed_pi_expr = true;
      }

      if (parsed_pi_expr) {
        size_t denom_pos = expr.find('/', pi_pos + 2);
        if (denom_pos != std::string::npos) {
          std::string denom_str = expr.substr(denom_pos + 1);
          try {
            double denom = std::stod(denom_str);
            if (denom == 0.0) {
              std::cerr << "Error: division by zero in '" << tokens[index] << "'.\n";
              return std::numeric_limits<double>::quiet_NaN();
            }
            theta /= denom;
          } catch (const std::exception&) {
            std::cerr << "Error: Invalid PI expression '" << tokens[index] << "'.\n";
            return std::numeric_limits<double>::quiet_NaN();
          }
        }
      }
    }

    if (parsed_pi_expr) {
      if (neg) theta = -theta;
    } else {
      try {
        theta = std::stod(token);
      } catch (const std::exception&) {
        std::cerr << "Error: Argument '" << tokens[index] << "' must be a number.\n";
        return std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  if (is_deg) {
    theta = theta * (pi / 180.0);
  }

  return theta;
}

double get_angle_arg_required(const std::vector<std::string>& tokens, size_t index, const std::string& cmd)
{
  if (index >= tokens.size()) {
    std::cerr << "Error: " << cmd << " requires more arguments.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
  return get_angle_arg(tokens, index, cmd);
}

} // namespace cli
