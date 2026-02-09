#include "algorithms/latin_square.hh"

static int decode_cell(Bitstring assignment, int cell_index)
{
  int shift = cell_index * 2;
  return static_cast<int>((assignment >> shift) & 0x3ULL);
}

static void decode_grid(Bitstring assignment, const int row0[3], int grid[3][3])
{
  grid[0][0] = row0[0];
  grid[0][1] = row0[1];
  grid[0][2] = row0[2];

  int idx = 0;
  for (int r = 1; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      grid[r][c] = decode_cell(assignment, idx);
      ++idx;
    }
  }
}

bool is_valid_latin3(Bitstring assignment, const int row0[3])
{
  int grid[3][3];
  decode_grid(assignment, row0, grid);

  // Validate values and row uniqueness.
  for (int r = 0; r < 3; ++r) {
    bool seen[3] = {false, false, false};
    for (int c = 0; c < 3; ++c) {
      int v = grid[r][c];
      if (v < 0 || v > 2) return false;
      if (seen[v]) return false;
      seen[v] = true;
    }
  }

  // Column uniqueness.
  for (int c = 0; c < 3; ++c) {
    bool seen[3] = {false, false, false};
    for (int r = 0; r < 3; ++r) {
      int v = grid[r][c];
      if (seen[v]) return false;
      seen[v] = true;
    }
  }

  return true;
}

bool is_valid_latin3_fixedrow(Bitstring assignment)
{
  const int row0[3] = {0, 1, 2};
  return is_valid_latin3(assignment, row0);
}
