#pragma once

#include <string>
#include <vector>

namespace cli {

std::vector<std::string> parse_command(const std::string& line);
std::string upper_copy(const std::string& token);
bool token_is(const std::string& token, const char* expected_uppercase);
int get_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
double get_double_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
double get_angle_arg(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);
double get_angle_arg_required(const std::vector<std::string>& tokens, size_t index, const std::string& cmd);

} // namespace cli
