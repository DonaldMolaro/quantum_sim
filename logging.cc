#include "logging.hh"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>

namespace qsim_log {

namespace {
Level g_level = Level::Normal;
std::ostream* g_stream = &std::cout;
bool g_initialized = false;

Level level_from_env(const char* env)
{
  if (!env || env[0] == '\0') return g_level;
  std::string token(env);
  Level parsed;
  if (parse_level(token, parsed)) {
    return parsed;
  }
  return g_level;
}

void init_from_env()
{
  if (g_initialized) return;
  g_initialized = true;

  if (const char* env = std::getenv("QSIM_VERBOSE")) {
    g_level = level_from_env(env);
    return;
  }

  if (const char* env = std::getenv("QSIM_LOG_LEVEL")) {
    g_level = level_from_env(env);
    return;
  }

  if (const char* env = std::getenv("QSIM_GROVER_VERBOSE")) {
    if (env[0] != '\0' && env[0] != '0') {
      g_level = Level::Verbose;
    }
  }
}

std::string upper_copy(const std::string& token)
{
  std::string out = token;
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return out;
}
} // namespace

void set_level(Level level)
{
  g_level = level;
  g_initialized = true;
}

Level get_level()
{
  init_from_env();
  return g_level;
}

bool parse_level(const std::string& token, Level& out)
{
  std::string upper = upper_copy(token);
  if (upper == "0" || upper == "OFF" || upper == "QUIET") {
    out = Level::Quiet;
    return true;
  }
  if (upper == "1" || upper == "ON" || upper == "NORMAL") {
    out = Level::Normal;
    return true;
  }
  if (upper == "2" || upper == "VERBOSE") {
    out = Level::Verbose;
    return true;
  }
  return false;
}

void set_stream(std::ostream* out)
{
  if (out) {
    g_stream = out;
  }
}

std::ostream* stream()
{
  init_from_env();
  return g_stream;
}

bool enabled(Level level)
{
  return static_cast<int>(get_level()) >= static_cast<int>(level);
}

void log(Level level, const std::string& message)
{
  if (!enabled(level)) return;
  std::ostream* out = stream();
  if (!out) return;
  (*out) << message;
}

} // namespace qsim_log
