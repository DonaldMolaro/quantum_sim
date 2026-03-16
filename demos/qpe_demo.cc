#include "demos/qpe_demo.hh"
#include "algorithms/qpe.hh"
#include "internal/limits.hh"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

static void print_qpe_result(const QpeResult& r)
{
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  Precision qubits : " << r.m << " (resolution 1/" << (1 << r.m) << ")\n";
  std::cout << "  True phase (frac): " << r.true_phase_frac
            << "  (" << (r.true_phase_frac * 2.0 * qsim::limits::PI) << " rad)\n";
  std::cout << "  Measured integer : " << r.measured_int << " / " << (1 << r.m) << "\n";
  std::cout << "  Est. phase (frac): " << r.estimated_frac
            << "  (" << r.estimated_radians << " rad)\n";
  std::cout << "  Winning state prob: " << (r.best_prob * 100.0) << "%\n";
  const double err = std::abs(r.estimated_frac - r.true_phase_frac);
  std::cout << "  Absolute error   : " << err
            << (err < 1.0 / (1 << r.m) + 1e-9 ? "  [within 1 bit]\n" : "  [exceeds 1 bit]\n");
}

void run_qpe_demo(int m, double phase_radians)
{
  std::cout << "\n=== Quantum Phase Estimation ===\n";
  std::cout << "Estimating phase of U = P(" << std::fixed << std::setprecision(4)
            << phase_radians << " rad) using " << m << " precision qubits.\n";
  std::cout << "Eigenstate: |1>   (U|1> = e^{i*phi}|1>)\n\n";

  QpeResult r = run_qpe(m, phase_radians);
  print_qpe_result(r);
  std::cout << "================================\n";
}

void run_qpe_demo()
{
  const double PI = qsim::limits::PI;

  std::cout << "\n=== Quantum Phase Estimation Demo ===\n";
  std::cout << "We demonstrate QPE on several phases with increasing precision.\n\n";

  struct Case { const char* label; double phase; int m; };
  const std::vector<Case> cases = {
      {"phi = 1/4  (T^2 gate, exact at m=2)",      PI / 2.0,     4},
      {"phi = 1/8  (T gate, exact at m=3)",         PI / 4.0,     4},
      {"phi = 1/3  (approx, not a power-of-2 frac)", 2.0*PI/3.0,  6},
      {"phi = 3/8  (exact at m=3)",                 3.0*PI/4.0,   4},
  };

  for (const auto& c : cases) {
    std::cout << "--- " << c.label << " ---\n";
    QpeResult r = run_qpe(c.m, c.phase);
    print_qpe_result(r);
    std::cout << "\n";
  }

  std::cout << "Key insight: when phi is an exact m-bit binary fraction,\n"
            << "QPE succeeds with 100% probability. For other phases it\n"
            << "concentrates probability on the nearest integer, giving\n"
            << "an approximation within 1/(2^m).\n"
            << "=====================================\n";
}
