#include "demos/grover_demo.hh"
#include "algorithms/api/grover_api.hh"
#include <algorithm>
#include <cmath>
#include <vector>

static int infer_qubits_from_targets(const std::vector<Bitstring>& targets)
{
  Bitstring max_target = 0ULL;
  for (size_t i = 0; i < targets.size(); ++i) {
    if (targets[i] > max_target) {
      max_target = targets[i];
    }
  }
  if (max_target == 0ULL) {
    return 1;
  }
  int n = 0;
  while (max_target != 0ULL) {
    ++n;
    max_target >>= 1;
  }
  return std::max(1, n);
}

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
  run_grover_search_multi(infer_qubits_from_targets(targets), targets);
}
