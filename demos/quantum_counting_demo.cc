#include "demos/quantum_counting_demo.hh"

#include "algorithms/quantum_counting.hh"
#include <iostream>
#include <vector>

void run_quantum_counting_cli(int n_qubits, const std::vector<Bitstring>& targets, int max_iterations)
{
  QuantumCountingResult result = run_quantum_counting(n_qubits, targets, max_iterations);
  if (!result.ok) {
    std::cerr << "QCOUNT error: " << result.error << "\n";
    return;
  }

  std::cout << "Quantum counting result:\n";
  std::cout << "  n_qubits=" << result.n_qubits << " (N=" << result.state_size << ")\n";
  std::cout << "  estimated_targets=" << result.estimated_targets
            << " (real=" << result.estimated_targets_real << ")\n";
  std::cout << "  estimated_theta=" << result.estimated_theta << "\n";
  std::cout << "  iterations_used=" << result.iterations_used << "\n";
}

void run_quantum_counting_demo()
{
  std::cout << "Quantum counting demo (n=3, targets={1,6})\n";
  run_quantum_counting_cli(3, std::vector<Bitstring>{1, 6}, -1);
}
