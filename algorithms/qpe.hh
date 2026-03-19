#pragma once

/**
 * Quantum Phase Estimation (QPE)
 *
 * Estimates the phase φ ∈ [0, 1) of a single-qubit phase gate P(2πφ)
 * using m precision qubits.  The unitary U = P(2πφ) has eigenstate |1>
 * with eigenvalue e^{2πiφ}, so QPE on |1> recovers φ to m bits.
 *
 * Circuit layout (qubit 0 = target/eigenstate, qubits 1..m = precision):
 *   1. X(0)                           — prepare |1>
 *   2. H(k) for k in 1..m            — superpose precision register
 *   3. CP(k, 0, 2πφ · 2^{k-1})       — controlled-U^{2^{k-1}}
 *   4. iQFT(1, m)                     — inverse QFT on precision register
 *   5. Read precision register        — most-probable bitstring → estimated φ
 */
struct QpeResult {
  int    m;                   // number of precision qubits
  double true_phase_frac;     // input φ as fraction of 2π (range 0..1)
  int    measured_int;        // raw integer readout from precision register
  double estimated_frac;      // measured_int / 2^m  (≈ true_phase_frac)
  double estimated_radians;   // estimated_frac * 2π
  double best_prob;           // probability of the winning bitstring
};

/**
 * Run QPE for a phase gate P(phase_radians) = P(2πφ).
 *
 * @param m              Number of precision qubits (1..20).
 * @param phase_radians  The angle φ·2π; the gate applies e^{i·phase_radians} to |1>.
 */
QpeResult run_qpe(int m, double phase_radians);
