#include "state.hh"
#include <iostream>

State& State::controlled_Rr(int control_j, int target_k, int r)
{
  // Phase applied only if both control and target bits are 1
  ComplexNumber phase = std::exp(2.0 * M_PI * IMAGINARY_UNIT_I / std::pow(2.0, r)); 

  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;
            
    // Check if control bit (j) is 1 AND target bit (k) is 1
    if (((current_b >> control_j) & 1) && ((current_b >> target_k) & 1)) {
      current_a *= phase;
    }
  }
  return *this;
}

/**
 * @brief Applies the controlled Inverse Phase Rotation (R_r_dag).
 * Inverse rotation phase is negative.
 */
State& State::controlled_Rr_dag(int control_j, int target_k, int r)
{
  // Inverse phase rotation: sign of angle is negated.
  ComplexNumber phase_dag = std::exp(-2.0 * M_PI * IMAGINARY_UNIT_I / std::pow(2.0, r)); 
        
  for (auto& pair : state_) {
    Bitstring current_b = pair.first;
    ComplexNumber& current_a = pair.second;
            
    if (((current_b >> control_j) & 1) && ((current_b >> target_k) & 1)) {
      current_a *= phase_dag;
    }
  }
  return *this;
}

using AmplitudeMap = std::unordered_map<Bitstring, ComplexNumber>;
const double EPSILON = 1e-9; 

State& State::qft(int start_qubit, int end_qubit) 
{
    // 1. Determine register size and total dimension N
    int n = end_qubit - start_qubit + 1;
    if (n <= 0 || n >= 63) {
        std::cerr << "Error: Invalid qubit range for QFT." << std::endl;
        return *this;
    }
    Bitstring N = 1ULL << n;
    
    // Calculate the overall scaling factor 1/sqrt(N)
    ComplexNumber overall_scale = 1.0 / std::sqrt((double)N);
    const ComplexNumber I(0.0, 1.0);
    const long double PI = std::acos(-1.0L);
    
    // 2. Prepare temporary map for accumulation (handling quantum interference)
    AmplitudeMap new_state_map;
    
    // 3. Iterate over every existing state in the superposition (Input j)
    for (const auto& pair : state_) {
        Bitstring current_B = pair.first;
        ComplexNumber A_j = pair.second;
        
        // Extract the current target register index j (input value)
        Bitstring j = extract_bits(current_B, start_qubit, end_qubit); 

        // 4. Calculate contribution to every possible output state k
        for (Bitstring k = 0; k < N; ++k) {
            
            // Calculate phase factor: exp(2 * pi * i * j * k / N)
            long double exponent = 2.0L * PI * static_cast<long double>(j) * static_cast<long double>(k)
                                   / static_cast<long double>(N);
            double exponent_d = static_cast<double>(exponent);
            ComplexNumber phase_factor = std::exp(I * exponent_d);
            
            // New amplitude contribution A_k = A_j * (1/sqrt(N)) * phase_factor
            ComplexNumber A_jk = A_j * overall_scale * phase_factor;
            
            if (std::abs(A_jk) < EPSILON) {
                continue;
            }

            // Construct the resulting full bitstring B_new by replacing j with k
            Bitstring B_new = replace_bits(current_B, start_qubit, end_qubit, k);
            
            // 5. Accumulate the amplitude (key step for interference/reduction)
            new_state_map[B_new] += A_jk;
        }
    }

    // 6. Update the internal state representation
    QuantumState final_state_data;
    for (const auto& entry : new_state_map) {
        if (std::abs(entry.second) > EPSILON) { 
            final_state_data.push_back(entry); 
        }
    }
    state_ = final_state_data; 
    
    return *this;
}

State& State::iqft(int start_qubit, int end_qubit) 
{
    // 1. Determine register size and total dimension N
    int n = end_qubit - start_qubit + 1;
    if (n <= 0 || n >= 63) {
        std::cerr << "Error: Invalid qubit range for IQFT." << std::endl;
        return *this;
    }
    Bitstring N = 1ULL << n;
    
    // Calculate the overall scaling factor 1/sqrt(N)
    ComplexNumber overall_scale = 1.0 / std::sqrt((double)N);
    const ComplexNumber I(0.0, 1.0); 
    const long double PI = std::acos(-1.0L);
    
    // 2. Prepare temporary map for accumulation
    AmplitudeMap new_state_map;
    
    // 3. Iterate over every existing state (Input j)
    for (const auto& pair : state_) {
        Bitstring current_B = pair.first;
        ComplexNumber A_j = pair.second;
        
        // Extract the current target register index j
        Bitstring j = extract_bits(current_B, start_qubit, end_qubit); 

        // 4. Calculate contribution to every possible output state k
        for (Bitstring k = 0; k < N; ++k) {
            
            // Calculate inverse phase factor: exp(-2 * pi * i * j * k / N)
            long double exponent = -2.0L * PI * static_cast<long double>(j) * static_cast<long double>(k)
                                   / static_cast<long double>(N);
            double exponent_d = static_cast<double>(exponent);
            ComplexNumber phase_factor = std::exp(I * exponent_d);
            
            // New amplitude contribution A_k = A_j * (1/sqrt(N)) * phase_factor
            ComplexNumber A_jk = A_j * overall_scale * phase_factor;
            
            if (std::abs(A_jk) < EPSILON) {
                continue;
            }

            // Construct the resulting full bitstring B_new by replacing j with k
            Bitstring B_new = replace_bits(current_B, start_qubit, end_qubit, k);
            
            // 5. Accumulate the amplitude
            new_state_map[B_new] += A_jk;
        }
    }

    // 6. Update the internal state representation
    QuantumState final_state_data;
    for (const auto& entry : new_state_map) {
        if (std::abs(entry.second) > EPSILON) { 
            final_state_data.push_back(entry); 
        }
    }
    state_ = final_state_data; 
    
    return *this;
}
