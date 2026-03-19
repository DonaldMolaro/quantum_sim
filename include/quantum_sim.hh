#pragma once

#include "state.hh"
#include "logging.hh"
#include "algorithms/latin_square.hh"
#include "algorithms/deutsch_jozsa.hh"
#include "algorithms/bernstein_vazirani.hh"
#include "algorithms/qubo.hh"
#include "algorithms/vqa_qaoa.hh"
#include "algorithms/qaoa.hh"
#include "algorithms/vqe.hh"
#include "algorithms/anneal.hh"
#include "algorithms/qrng.hh"
#include "algorithms/tsp.hh"
#include "algorithms/quantum_counting.hh"
#include "algorithms/simon.hh"
#include "algorithms/api/grover_api.hh"
#include "algorithms/api/shor_api.hh"
#include "modular_exp.hh"
#include "math/bit_ops.hh"

namespace qsim {

// Release-facing convenience wrapper. This keeps common simulator operations
// on a smaller, more stable surface than exposing every State member directly.
class Simulator {
public:
  explicit Simulator(int qubits, int classical_bits = 0)
      : state_(qubits, classical_bits) {}

  Simulator& h(int qubit) { state_.h(qubit); return *this; }
  Simulator& x(int qubit) { state_.x(qubit); return *this; }
  Simulator& y(int qubit) { state_.y(qubit); return *this; }
  Simulator& z(int qubit) { state_.z(qubit); return *this; }
  Simulator& s(int qubit) { state_.s(qubit); return *this; }
  Simulator& t(int qubit) { state_.t(qubit); return *this; }
  Simulator& rx(int qubit, double theta) { state_.rx(qubit, theta); return *this; }
  Simulator& ry(int qubit, double theta) { state_.ry(qubit, theta); return *this; }
  Simulator& rz(int qubit, double theta) { state_.rz(qubit, theta); return *this; }
  Simulator& cx(int control, int target) { state_.cx(control, target); return *this; }
  Simulator& cz(int control, int target) { state_.cz(control, target); return *this; }
  Simulator& swap(int a, int b) { state_.swap(a, b); return *this; }
  Simulator& measure(int qubit, unsigned long classical_bit) { state_.measure(qubit, classical_bit); return *this; }
  Simulator& reset(int qubit) { state_.reset(qubit); return *this; }

  void display(bool show_all = false) const { state_.display(show_all); }
  std::vector<int> classical_bits() const { return state_.get_cbits(); }
  ComplexNumber amplitude(Bitstring basis_state) const { return state_.get_amplitude(basis_state); }

  State& raw_state() { return state_; }
  const State& raw_state() const { return state_; }

private:
  State state_;
};

} // namespace qsim
