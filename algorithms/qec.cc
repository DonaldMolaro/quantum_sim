#include "algorithms/qec.hh"
#include "state.hh"
#include <iostream>
#include <stdexcept>
#include <string>

// Layout: physical qubits 0,1,2 = data; qubits 3,4 = syndrome ancilla.
// Classical registers: cbit 0 = s1, cbit 1 = s2.

QecResult run_qec(int logical_bit, int error_qubit)
{
  if (logical_bit != 0 && logical_bit != 1)
    throw std::invalid_argument("QEC: logical_bit must be 0 or 1");
  if (error_qubit < -1 || error_qubit > 2)
    throw std::invalid_argument("QEC: error_qubit must be -1 (none), 0, 1, or 2");

  // 5 physical qubits, 2 classical bits.
  State s(5, 2);

  // --- Encoding ---
  // Prepare logical qubit in q0.
  if (logical_bit == 1) s.x(0);
  // Spread to q1 and q2 via CX.
  s.cx(0, 1);
  s.cx(0, 2);
  // Now: |000> if logical_bit=0, |111> if logical_bit=1.

  // --- Inject error ---
  if (error_qubit >= 0) s.x(error_qubit);

  // --- Syndrome measurement ---
  // s1 = parity(q0, q1): ancilla q3
  s.cx(0, 3);
  s.cx(1, 3);
  s.measure(3, 0); // s1 into cbit 0
  // s2 = parity(q1, q2): ancilla q4
  s.cx(1, 4);
  s.cx(2, 4);
  s.measure(4, 1); // s2 into cbit 1

  const int s1 = s.get_cbit(0);
  const int s2 = s.get_cbit(1);
  const int syndrome = s1 | (s2 << 1);

  // --- Correction ---
  int corrected_qubit = -1;
  if      (syndrome == 0b01) { s.x(0); corrected_qubit = 0; } // s1=1,s2=0 → error on q0
  else if (syndrome == 0b11) { s.x(1); corrected_qubit = 1; } // s1=1,s2=1 → error on q1
  else if (syndrome == 0b10) { s.x(2); corrected_qubit = 2; } // s1=0,s2=1 → error on q2

  // --- Decode ---
  s.cx(0, 2);
  s.cx(0, 1);

  // --- Verify ---
  // Measure qubit 0: should match logical_bit.
  s.measure(0, 0);
  const bool success = (s.get_cbit(0) == logical_bit);

  return {logical_bit, error_qubit, syndrome, corrected_qubit, success};
}

void run_qec_demo()
{
  State::set_default_log_stream(nullptr);
  std::cout << "\n=== Quantum Error Correction: 3-Qubit Bit-Flip Code ===\n";
  std::cout << "Encoding: |0>_L = |000>, |1>_L = |111>\n";
  std::cout << "Syndrome: s1=parity(q0,q1), s2=parity(q1,q2)\n\n";

  for (int logical_bit : {0, 1}) {
    for (int err : {-1, 0, 1, 2}) {
      QecResult r = run_qec(logical_bit, err);
      std::cout << "  logical=" << logical_bit
                << "  error=" << (err == -1 ? "none" : ("q" + std::to_string(err)))
                << "  syndrome=" << r.syndrome
                << "  corrected=" << (r.corrected_qubit == -1 ? "none" : ("q" + std::to_string(r.corrected_qubit)))
                << "  result=" << (r.recovery_success ? "OK" : "FAIL") << "\n";
    }
  }

  std::cout << "\nKey insight: syndrome (s1,s2) uniquely identifies which physical qubit\n"
            << "flipped without revealing the logical state — no measurement on data qubits.\n"
            << "======================================================\n";
}
