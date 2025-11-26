#include <vector>
#include <complex>
#include <utility>
#include <functional>
#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib> // For std::rand and RAND_MAX
#include <algorithm> // For std::remove_if
#include "state.hh"
// --- Type Definitions based on Source Material ---

// Amplitude is a complex number (using C++11 std::complex<double>)
using ComplexNumber = std::complex<double>; 
// Bitstring represents the combination of qubit values (the key)
using Bitstring = unsigned long long; 
// Each element of the state set is a (bitstring, amplitude) pair
using QubitAmplitudePair = std::pair<Bitstring, ComplexNumber>;
// The QuantumState is the set/vector of these pairs
using QuantumState = std::vector<QubitAmplitudePair>;
// Intermediate state used during flatMap, potentially containing duplicate keys
using IntermediateState = QuantumState; 

// Global Constants

// --- Helper Functions (used internally by gates) ---

/** Reads the jth bit (b_j) of the bitstring b. */
int get_jth_bit(const Bitstring& b, int j)
{
    return (b >> j) & 1; 
}

/** Flips the jth bit (b¬j). */
Bitstring flip_jth_bit(const Bitstring& b, int j) {
    return b ^ (1ULL << j);
}

/** Sets the jth bit of bitstring b to a specific value (0 or 1). */
Bitstring set_jth_bit(const Bitstring& b, int j, int value) {
    if (value == 0) {
        return b & ~(1ULL << j);
    } else {
        return b | (1ULL << j);
    }
}

/**
 * s.map(λb, a. (b', ca)) 
 * Applies a transformation independently to every element.
 */
void State::s_map(const std::function<QubitAmplitudePair(const Bitstring&, const ComplexNumber&)>& transformation_func)
{
  QuantumState new_state;
  // Use C++11 range-based for loop
  for (const auto& pair : state_) {
    new_state.push_back(transformation_func(pair.first, pair.second));
  }
  state_ = new_state;
}
    
/**
 * Implements .reduceByKey(λx, y. x+ y)
 * Sums amplitudes for identical bitstrings (keys), used after flatMap.
 */
QuantumState State::reduceByKey(const IntermediateState& intermediate_state) const
{
  std::map<Bitstring, ComplexNumber> reduced_map;
        
  // Sum amplitudes
  for (const auto& pair : intermediate_state) {
    reduced_map[pair.first] += pair.second;
  }

  // Convert map back to QuantumState vector
  QuantumState final_state;
  for (const auto& pair : reduced_map) {
    // Filter out elements with zero amplitude
    if (std::abs(pair.second) > 1e-9) { 
      final_state.push_back(pair);
    }
  }
  return final_state;
}

/**
 * s.flatMap(λb, a. {set_of_pairs}) 
 * Applies a transformation that can yield multiple resulting elements, 
 * followed by reduction (used specifically for Hadamard).
 */
void State::s_flatMap_and_reduce(const std::function<IntermediateState(const Bitstring&, const ComplexNumber&)>& transformation_func)
{
  IntermediateState intermediate_state;
        
  // Apply flatMap: take the union of all generated sets
  for (const auto& pair : state_) {
    IntermediateState generated_set = transformation_func(pair.first, pair.second);
    intermediate_state.insert(intermediate_state.end(), generated_set.begin(), generated_set.end());
  }
  // Apply reduceByKey to sum resulting amplitudes
  state_ = reduceByKey(intermediate_state);
}

/** Constructor: Initializes the state to the ground state |00...0>. */
State::State(int N, int num_cbits) : num_qubits_(N), cbits_(num_cbits, 0) 
{
  // Start state: bitstring 0, amplitude 1 (the ground state)
  state_.push_back({0ULL, ONE_COMPLEX}); 
}

// --- Public Gate Methods ---

/** X Gate (NOT): s.map(λb, a. (b¬j , a)) */
State& State::x(int j)
{
  s_map([j](const Bitstring& b, const ComplexNumber& a) {
    Bitstring b_prime = flip_jth_bit(b, j);
    return std::make_pair(b_prime, a);
  });
  return *this;
}
    
/** CX Gate (Controlled X): s.map(λb, a. (ite(b_j, b¬k, b), a)) */
State& State::cx(int j_control, int k_target)
{
  s_map([j_control, k_target](const Bitstring& b, const ComplexNumber& a) {
    // If control bit b_j is 1
    if (get_jth_bit(b, j_control) == 1) {
      Bitstring b_prime = flip_jth_bit(b, k_target);
      return std::make_pair(b_prime, a);
    }
    // Otherwise, b remains unchanged
    return std::make_pair(b, a);
  });
  return *this;
}
    
/** S Gate (Phase): s.map(λb, a. (b, a * i^b_j)) */
State& State::s(int j)
{
  s_map([j](const Bitstring& b, const ComplexNumber& a) {
    // Determine phase factor: i^b_j
    ComplexNumber phase_factor = (get_jth_bit(b, j) == 1) ? IMAGINARY_UNIT_I : ONE_COMPLEX;
    return std::make_pair(b, a * phase_factor);
  });
  return *this;
}
    
/** T Gate (Phase): s.map(λb, a. (b, a * ((1+i)/sqrt(2))^b_j)) */
State& State::t(int j)
{
  const ComplexNumber T_GATE_CONSTANT(ONE_OVER_SQRT_TWO, ONE_OVER_SQRT_TWO);
        
  s_map([j, T_GATE_CONSTANT](const Bitstring& b, const ComplexNumber& a) {
    // Determine phase factor
    ComplexNumber phase_factor = (get_jth_bit(b, j) == 1) ? T_GATE_CONSTANT : ONE_COMPLEX;
    return std::make_pair(b, a * phase_factor);
  });
  return *this;
}
    
/** Hadamard Gate (Superposition): flatMap().reduceByKey() */
State& State::h(int j)
{
  auto h_transformation = [j](const Bitstring& b, const ComplexNumber& a) -> IntermediateState {
    int bj = get_jth_bit(b, j);
    IntermediateState generated_set;

    // Pair 1: b_j=0, Amplitude: a / sqrt(2)
    generated_set.push_back(std::make_pair(set_jth_bit(b, j, 0), a * ONE_OVER_SQRT_TWO));
            
    // Pair 2: b_j=1, Amplitude: a * (1 - 2*b_j) / sqrt(2)
    double coefficient = (1.0 - 2.0 * (double)bj) * ONE_OVER_SQRT_TWO;
    generated_set.push_back(std::make_pair(set_jth_bit(b, j, 1), a * coefficient));
    return generated_set;
  };
  s_flatMap_and_reduce(h_transformation);
  return *this;
}
    
// --- Measurement Method ---

/**
 * Computes the probability of measuring 0 for qubit j: 
 * s.filter(λb, a. bj = 0) .map(λb, a. |a|^2) .sum()
 */
double State::compute_probability_of_0(int j) const
{
  double probability = 0.0;
        
  // 1. Filter and Map (in one loop for efficiency): 
  // Iterate through states and accumulate squared magnitudes (|a|^2) if b_j = 0.
  for (const auto& pair : state_) {
    if (get_jth_bit(pair.first, j) == 0) {
      // Sum: using std::norm to calculate |a|^2
      probability += std::norm(pair.second);
    }
  }
  return probability;
}

// Accessor for observation
const QuantumState& State::get_state() const { return state_; }

/*
void State::display() const
{
  std::cout << "\nQuantum State (N=" << num_qubits_ << "): \n";
        
  // Output Bitstrings (The Keys)
  std::cout << "{ ";
  for (const auto& pair : state_) {
    std::cout << bitstring_to_string(pair.first, num_qubits_) << " ";
  }
  std::cout << "}\n";

  // Output Amplitudes (The Values)
  std::cout << "{ ";
  for (const auto& pair : state_) {
    std::cout << complex_to_string(pair.second) << " ";
  }
  std::cout << "}\n";

  // Display Classical Registers (if applicable)
  if (!cbits_.empty()) {
    std::cout << "Classical Register: [";
    for (size_t i = 0; i < cbits_.size(); ++i) {
      std::cout << cbits_[i] << (i == cbits_.size() - 1 ? "" : ", ");
    }
    std::cout << "]\n";
  }
}
*/
/**
 * @brief Measures the jth qubit, collapses the state, and stores the result.
 * 
 * @param j Qubit index to measure.
 * @param cbit_index Index of classical register to store the result.
 * @return Reference to the updated State object.
 */
State& State::measure(int j, unsigned long cbit_index) {
  // 1. Compute Pj,0: Probability of measuring 0
  double p0 = compute_probability_of_0(j);
  // 2. Randomly determine the outcome (0 or 1)
  double random_val = (double)std::rand() / RAND_MAX;
  int outcome; 
  double p_outcome; // The probability corresponding to the chosen outcome
        
  if (random_val <= p0) {
    outcome = 0;
    p_outcome = p0;
  } else {
    outcome = 1;
    p_outcome = 1.0 - p0; // Pj,1
  }
        
  // Check probability stability (ensure p_outcome is not zero, though theoretically 
  // the state should already be clean if p_outcome=0)
  if (p_outcome < 1e-9) {
    // If probability is essentially zero, we cannot normalize. 
    // We set the factor to 1 and let the map zero out any remaining noise.
    p_outcome = 1.0; 
  }

  // 3. Calculate the renormalization factor (1 / sqrt(p))
  // The factor required is p^-0.5.
  ComplexNumber norm_factor(1.0 / std::sqrt(p_outcome), 0.0);
        
  // 4. Post-measurement State Update (Collapse and Renormalize)
  // This is implemented via s.map(λb, a. (b, new_amplitude))
  s_map([j, outcome, norm_factor](const Bitstring& b, const ComplexNumber& a) {
    int bj = get_jth_bit(b, j);
            
    if (bj == outcome) {
      // If the qubit matches the measured outcome, renormalize.
      // New amplitude = a * p^-0.5
      return std::make_pair(b, a * norm_factor);
    } else {
      // If the qubit does NOT match the measured outcome, zero out (collapse).
      // This mimics the (1 - bj) coefficient used in the measure_j,0 pseudocode.
      return std::make_pair(b, ComplexNumber(0.0, 0.0));
    }
  });

  // 5. Clean up the state vector by removing elements with zero amplitude
  // This step is vital for efficiency, though not explicitly in the pseudocode.
  state_.erase(std::remove_if(state_.begin(), state_.end(), 
			      [](const QubitAmplitudePair& pair) { 
				return std::abs(pair.second) < 1e-9; 
			      }), 
	       state_.end());

  // 6. Store the result in the classical register
  if (!cbits_.empty())
    {
      if (cbit_index >= 0 && cbit_index < cbits_.size()) {
	cbits_[cbit_index] = outcome;
      }
    }
  else
    {
      std::cout << "No Classical Registers \n";
    }
  return *this;
}


