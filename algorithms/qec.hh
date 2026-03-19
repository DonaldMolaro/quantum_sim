#pragma once
#include <string>

/**
 * Quantum Error Correction — 3-qubit bit-flip code.
 *
 * Encodes a single logical qubit (q0) using two ancilla qubits (q1, q2):
 *   |0>_L = |000>, |1>_L = |111>
 *
 * Syndrome measurement uses two additional ancilla qubits (q3, q4):
 *   s1 = parity(q0, q1), s2 = parity(q1, q2)
 *
 * Correction table:
 *   s1=0, s2=0 → no error
 *   s1=1, s2=0 → flip q0
 *   s1=0, s2=1 → flip q2
 *   s1=1, s2=1 → flip q1
 */
struct QecResult {
  int  logical_input;    // 0 or 1
  int  error_qubit;      // -1 = no error, else index of injected X error
  int  syndrome;         // 2-bit syndrome value (s1 | s2<<1)
  int  corrected_qubit;  // qubit the correction was applied to (-1 if none)
  bool recovery_success; // true if logical qubit was recovered correctly
};

/** Run the full 3-qubit bit-flip QEC demo. */
void run_qec_demo();

/**
 * Encode, optionally inject an error, measure syndrome, correct, and verify.
 * @param logical_bit  Input logical qubit value (0 or 1).
 * @param error_qubit  Physical qubit to flip (-1 = no error).
 */
QecResult run_qec(int logical_bit, int error_qubit);
