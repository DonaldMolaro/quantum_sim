#include "demos/grover_demo.hh"
#include "algorithms/api/grover_api.hh"
#include "internal/format_utils.hh"
#include <vector>

void run_grover_search(int n_qubits, Bitstring target_w)
{
  (void)run_grover(n_qubits, {target_w});
}

void run_grover_search_multi(int n_qubits, const std::vector<Bitstring>& targets)
{
  if (targets.empty()) return;
  (void)run_grover(n_qubits, targets);
}

void run_grover_search(State *s, Bitstring target_w)
{
  if (s) {
    (void)run_grover(*s, {target_w});
    return;
  }
  run_grover_search(1, target_w);
}

void run_grover_search_multi(State *s, const std::vector<Bitstring>& targets)
{
  if (targets.empty()) return;
  if (s) {
    (void)run_grover(*s, targets);
    return;
  }
  run_grover_search_multi(qsim::format::infer_qubits_from_targets(targets), targets);
}
