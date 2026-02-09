#include "demos/grover_demo.hh"
#include "algorithms/api/grover_api.hh"
#include <vector>

void run_grover_search(State *s, Bitstring target_w)
{
  if (!s) return;
  run_grover(*s, {target_w});
}

void run_grover_search_multi(State *s, const std::vector<Bitstring>& targets)
{
  if (!s || targets.empty()) return;
  run_grover(*s, targets);
}
