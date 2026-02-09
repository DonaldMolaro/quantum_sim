#pragma once

#include <iosfwd>
#include <string>

namespace qsim_log {

enum class Level {
  Quiet = 0,
  Normal = 1,
  Verbose = 2
};

void set_level(Level level);
Level get_level();
bool parse_level(const std::string& token, Level& out);

void set_stream(std::ostream* out);
std::ostream* stream();

bool enabled(Level level);
void log(Level level, const std::string& message);

} // namespace qsim_log
