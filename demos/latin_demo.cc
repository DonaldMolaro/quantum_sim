#include "demos/latin_demo.hh"
#include "algorithms/latin_square.hh"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

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

static void print_grid(Bitstring assignment, const int row0[3])
{
  int grid[3][3];
  decode_grid(assignment, row0, grid);
  std::cout << "Latin square (row0 fixed):\n";
  for (int r = 0; r < 3; ++r) {
    std::cout << "  ";
    for (int c = 0; c < 3; ++c) {
      std::cout << grid[r][c] << (c == 2 ? "" : " ");
    }
    std::cout << "\n";
  }
}

static int count_latin3_row0(const int row0[3], std::vector<Bitstring>* out_solutions)
{
  int count = 0;
  for (Bitstring b = 0; b < (1ULL << 12); ++b) {
    if (is_valid_latin3(b, row0)) {
      ++count;
      if (out_solutions) {
        out_solutions->push_back(b);
      }
    }
  }
  return count;
}

static void validate_row0_or_die(const int row0[3])
{
  bool seen[3] = {false, false, false};
  for (int i = 0; i < 3; ++i) {
    if (row0[i] < 0 || row0[i] > 2) {
      std::cerr << "Row0 values must be 0..2.\n";
      std::exit(1);
    }
    if (seen[row0[i]]) {
      std::cerr << "Row0 must be a permutation of 0,1,2.\n";
      std::exit(1);
    }
    seen[row0[i]] = true;
  }
}

void run_latin3_grover_demo_row0(const int row0[3], int iterations)
{
  const int n_qubits = 12; // 6 cells * 2 bits
  const double N = std::ldexp(1.0, n_qubits);
  const double PI = std::acos(-1.0);
  validate_row0_or_die(row0);

  int M = count_latin3_row0(row0, nullptr);
  if (const char* env = std::getenv("QSIM_LATIN_FORCE_M")) {
    M = std::atoi(env);
  }
  if (M == 0) {
    std::cout << "No solutions for given row0.\n";
    return;
  }

  const int R_default = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N / static_cast<double>(M))));
  const int R = (iterations >= 0) ? iterations : R_default;

  std::cout << "Grover Latin-3 demo (row0 fixed).\n";
  std::cout << "Row0: " << row0[0] << " " << row0[1] << " " << row0[2] << "\n";
  std::cout << "State size: 2^" << n_qubits << " = " << N << ", expected solutions: " << M << "\n";
  std::cout << "Iterations: " << R << "\n";

  State s(n_qubits, n_qubits);
  for (int j = 0; j < n_qubits; ++j) {
    s.h(j);
  }

  for (int k = 0; k < R; ++k) {
    s.phase_flip_if([&row0](Bitstring b) { return is_valid_latin3(b, row0); });
    s.grover_diffusion_Us();
  }

  // Find the best candidate among valid assignments.
  Bitstring best = 0;
  double best_prob = -1.0;
  for (const auto& pair : s.get_state()) {
    Bitstring b = pair.first;
    if (!is_valid_latin3(b, row0)) continue;
    double prob = std::norm(pair.second);
    if (prob > best_prob) {
      best_prob = prob;
      best = b;
    }
  }

  std::cout << "Best candidate probability: " << best_prob << "\n";
  print_grid(best, row0);

  // Sample one assignment by measuring all qubits.
  Bitstring measured = 0;
  if (const char* forced = std::getenv("QSIM_LATIN_FORCE_MEASURED")) {
    measured = static_cast<Bitstring>(std::strtoull(forced, nullptr, 10));
  } else {
    for (int q = 0; q < n_qubits; ++q) {
      s.measure(q, q);
    }
    for (int q = 0; q < n_qubits; ++q) {
      measured |= (static_cast<Bitstring>(s.get_cbit(q)) << q);
    }
  }
  std::cout << "Measured candidate:\n";
  if (is_valid_latin3(measured, row0)) {
    print_grid(measured, row0);
  } else {
    std::cout << "  (measurement not a valid Latin square; try more iterations)\n";
  }
}

void run_latin3_grover_demo(int iterations)
{
  const int row0[3] = {0, 1, 2};
  run_latin3_grover_demo_row0(row0, iterations);
}

void run_latin3_count_row0(const int row0[3])
{
  validate_row0_or_die(row0);
  int count = count_latin3_row0(row0, nullptr);
  std::cout << "Latin-3 solutions for row0: "
            << row0[0] << " " << row0[1] << " " << row0[2]
            << " => " << count << "\n";
}

void run_latin3_print_all_row0(const int row0[3])
{
  validate_row0_or_die(row0);
  std::vector<Bitstring> solutions;
  int count = count_latin3_row0(row0, &solutions);
  int max_print = -1;
  if (const char* env = std::getenv("QSIM_LATIN_MAX_PRINT")) {
    max_print = std::atoi(env);
  }
  std::cout << "Latin-3 solutions for row0: "
            << row0[0] << " " << row0[1] << " " << row0[2]
            << " => " << count << "\n";
  int printed = 0;
  for (Bitstring b : solutions) {
    if (max_print >= 0 && printed >= max_print) {
      break;
    }
    print_grid(b, row0);
    ++printed;
  }
}
