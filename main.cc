#include "state.hh"
#include "shell.hh"

int main()
{
  QuantumShell qs;
  qs.run();
  return 0;
};

/*

extern void run_grover_search(int n_qubits,Bitstring targer_w);

int main()
{
  const int N_QUBITS = 4;
  const Bitstring TARGET_W = 13;
  run_grover_search(N_QUBITS,TARGET_W);
  return 0;
};

*/



