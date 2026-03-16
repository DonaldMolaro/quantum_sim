#include "state.hh"
#include "internal/limits.hh"
#include "math/bit_ops.hh"
#include "internal/qft_utils.hh"
#include <iostream>
#include <unordered_map>

State& State::apply_controlled_Rr(int control_j, int target_k, int r, int sign)
{
  ComplexNumber phase = std::exp(IMAGINARY_UNIT_I * qft_rotation_angle(r, sign));
  for (auto& pair : state_) {
    if (((pair.first >> control_j) & 1) && ((pair.first >> target_k) & 1)) {
      pair.second *= phase;
    }
  }
  return *this;
}

State& State::controlled_Rr(int control_j, int target_k, int r)
{
  return apply_controlled_Rr(control_j, target_k, r, 1);
}

/** @brief Applies the controlled Inverse Phase Rotation (R_r_dag). */
State& State::controlled_Rr_dag(int control_j, int target_k, int r)
{
  return apply_controlled_Rr(control_j, target_k, r, -1);
}

using AmplitudeMap = std::unordered_map<Bitstring, ComplexNumber>;

static QuantumState qft_direct_transform(const QuantumState& input_state,
                                         int start_qubit,
                                         int end_qubit,
                                         Bitstring N,
                                         int sign)
{
  const ComplexNumber overall_scale = 1.0 / std::sqrt(static_cast<double>(N));
  const ComplexNumber I(0.0, 1.0);
  AmplitudeMap new_state_map;

  for (const auto& pair : input_state) {
    const Bitstring current_B = pair.first;
    const ComplexNumber A_j = pair.second;
    const Bitstring j = extract_bits(current_B, start_qubit, end_qubit);

    for (Bitstring k = 0; k < N; ++k) {
      const double exponent = qft_phase_exponent(j, k, N, sign);
      const ComplexNumber phase_factor = std::exp(I * exponent);
      const ComplexNumber A_jk = A_j * overall_scale * phase_factor;
      if (std::abs(A_jk) < qsim::limits::AMPLITUDE_EPSILON) {
        continue;
      }
      const Bitstring B_new = replace_bits(current_B, start_qubit, end_qubit, k);
      new_state_map[B_new] += A_jk;
    }
  }

  QuantumState output_state;
  output_state.reserve(new_state_map.size());
  for (const auto& entry : new_state_map) {
    if (std::abs(entry.second) > qsim::limits::AMPLITUDE_EPSILON) {
      output_state.push_back(entry);
    }
  }
  return output_state;
}

static State& qft_gate(State& s, int start_qubit, int end_qubit)
{
  int n = end_qubit - start_qubit + 1;
  if (!qft_range_valid(n)) {
    std::cerr << "Error: Invalid qubit range for QFT." << std::endl;
    return s;
  }

  for (int i = start_qubit; i <= end_qubit; ++i) {
    s.h(i);
    for (int j = i + 1; j <= end_qubit; ++j) {
      int r = j - i + 1;
      s.controlled_Rr(j, i, r);
    }
  }

  return s;
}

static State& iqft_gate(State& s, int start_qubit, int end_qubit)
{
  int n = end_qubit - start_qubit + 1;
  if (!qft_range_valid(n)) {
    std::cerr << "Error: Invalid qubit range for IQFT." << std::endl;
    return s;
  }

  for (int i = end_qubit; i >= start_qubit; --i) {
    s.h(i);
    for (int j = i - 1; j >= start_qubit; --j) {
      int r = i - j + 1;
      s.controlled_Rr_dag(i, j, r);
    }
  }

  return s;
}

State& State::qft(int start_qubit, int end_qubit) 
{
  if (qft_mode_ == QftMode::Gate) {
    return qft_gate(*this, start_qubit, end_qubit);
  }

  // 1. Determine register size and total dimension N
  int n = end_qubit - start_qubit + 1;
  if (!qft_range_valid(n)) {
      std::cerr << "Error: Invalid qubit range for QFT." << std::endl;
      return *this;
  }
  Bitstring N = 0ULL;
  qft_dimension(n, N);
  
  state_ = qft_direct_transform(state_, start_qubit, end_qubit, N, 1);
  return *this;
}

State& State::iqft(int start_qubit, int end_qubit) 
{
  if (qft_mode_ == QftMode::Gate) {
    return iqft_gate(*this, start_qubit, end_qubit);
  }

  // 1. Determine register size and total dimension N
  int n = end_qubit - start_qubit + 1;
  if (!qft_range_valid(n)) {
      std::cerr << "Error: Invalid qubit range for IQFT." << std::endl;
      return *this;
  }
  Bitstring N = 0ULL;
  qft_dimension(n, N);
  
  state_ = qft_direct_transform(state_, start_qubit, end_qubit, N, -1);
  return *this;
}
