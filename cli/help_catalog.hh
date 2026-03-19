#pragma once

#include <iosfwd>

namespace cli {
namespace help_catalog {

enum class Section {
  Gates,
  Algorithms,
  Utility,
};

void print_section(std::ostream& out, Section section);
void print_all(std::ostream& out);

} // namespace help_catalog
} // namespace cli
