#include "algorithms/qaoa.hh"

QaoaResult run_qaoa_qubo(int n,
                         const std::vector<double>& q,
                         const QaoaOptions& options)
{
  return run_vqa_qaoa_qubo(n, q, options);
}
