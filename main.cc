/* Quantum Simulator *without* linear algebra.
 *
 * This is an implementation of an excellent paper by Aws Albarghouhi
 * Aws gives an implementation in Python, but I think in C++ hence this
 * implementation.
 *
 * If you actually want to *learn* something about quantum computation I
 * would strong suggest you give this paper a read.
 *
 * https://eprint.iacr.org/2025/1091.pdf
 *
 * Don - December 2025
 */
#include "state.hh"
#include "shell.hh"

int main()
{
  QuantumShell qs;
  qs.run();
  return 0;
};



