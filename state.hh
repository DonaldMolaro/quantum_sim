/*
 * Almost all of the code for this quantum simulator was created by NotebookLM
 * The paper it is based on is public, and NotebookLM is public so in my opinion
 * this code should remain public.
 * 
 * This contains the "state" of the quantum machine, since this is a simulator
 * one can "peek" inside and see the state of the Qubits.
 */
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
const ComplexNumber ONE_COMPLEX(1.0, 0.0);
const ComplexNumber IMAGINARY_UNIT_I(0.0, 1.0);
const double ONE_OVER_SQRT_TWO = 1.0 / std::sqrt(2.0); 

// --- Helper Functions (used internally by gates) ---

/** Reads the jth bit (b_j) of the bitstring b. */
extern int get_jth_bit(const Bitstring& b, int j);

/** Flips the jth bit (b¬j). */
extern Bitstring flip_jth_bit(const Bitstring& b, int j);
/** Sets the jth bit of bitstring b to a specific value (0 or 1). */

extern Bitstring set_jth_bit(const Bitstring& b, int j, int value);

/** Converts the Bitstring (ULL) to a padded binary string based on N qubits. */
extern std::string bitstring_to_string(Bitstring b, int N);

/** Converts the ComplexNumber to a readable string (e.g., 0.707 + 0.707i). */
extern std::string complex_to_string(const ComplexNumber& a);
Bitstring power_mod(Bitstring base, Bitstring exp, Bitstring mod);
Bitstring extract_bits(Bitstring b, int start, int end);
Bitstring replace_bits(Bitstring b, int start, int end, Bitstring new_val);
void accumulate(Bitstring b, ComplexNumber a);

class State
{
private:
  QuantumState state_;
  int num_qubits_;
  std::vector<int> cbits_;

  // Accessor for the size of the classical register
  size_t cbits_size() const { return cbits_.size(); }

  // --- Core Functional Operations ---
  // These methods implement the set transformations defining quantum gates.
  /**
   * s.map(λb, a. (b', ca)) 
   * Applies a transformation independently to every element.
   */
  void s_map(const std::function<QubitAmplitudePair(const Bitstring&, const ComplexNumber&)>& transformation_func);
  /**
   * Implements .reduceByKey(λx, y. x+ y)
   * Sums amplitudes for identical bitstrings (keys), used after flatMap.
   */
  QuantumState reduceByKey(const IntermediateState& intermediate_state) const;
  /**
   * s.flatMap(λb, a. {set_of_pairs}) 
   * Applies a transformation that can yield multiple resulting elements, 
   * followed by reduction (used specifically for Hadamard).
   */
  void s_flatMap_and_reduce(const std::function<IntermediateState(const Bitstring&, const ComplexNumber&)>& transformation_func);
  // Helper to find or add a bitstring in the vector representation (O(N) operation)
  QubitAmplitudePair* find_or_add(QuantumState& state, Bitstring b) {
    // Find existing bitstring
    for (auto& pair : state) {
      if (pair.first == b) {
	return &pair;
      }
    }
    // Not found, add new entry with zero amplitude
    state.push_back({b, 0.0});
    return &state.back();
    }

public:
  /** Constructor: Initializes the state to the ground state |00...0>. */
  State(int N, int num_cbits);

  State(int n) : num_qubits_(n) {
    state_.push_back({0ULL, 1.0});
  }
  
  State(Bitstring initial_state_value = 0)
  {
    state_.push_back({initial_state_value, 1.0});
  }


  State& x(int j);
  /** CX Gate (Controlled X): s.map(λb, a. (ite(b_j, b¬k, b), a)) */
  State& cx(int j_control, int k_target);
  /** S Gate (Phase): s.map(λb, a. (b, a * i^b_j)) */
  State& s(int j);
  /** T Gate (Phase): s.map(λb, a. (b, a * ((1+i)/sqrt(2))^b_j)) */
  State& t(int j);
  /** Hadamard Gate (Superposition): flatMap().reduceByKey() */
  State& h(int j);
  // --- Measurement Method ---
  /**
   * Computes the probability of measuring 0 for qubit j: 
   * s.filter(λb, a. bj = 0) .map(λb, a. |a|^2) .sum()
   */
  /**
   * @brief Applies the SWAP gate between qubit j and qubit k.
   * This operation swaps the bit values (keys) in the quantum state representation 
   * but leaves the amplitudes unchanged, mirroring the core (b, a) -> (b', a) map.
   * 
   * @param j Index of the first qubit (0-indexed, LSB).
   * @param k Index of the second qubit.
   * @return State& Reference to the modified State object.
   */
  State& swap(int i, int k);
  /*
   * Methods needed for Shor's algorithm.
   */
  State& controlled_Rr(int control_j, int target_k, int r);
  State& controlled_Rr_dag(int control_j, int target_k, int r);
  State& qft(int start_qubit, int end_qubit);
  State& iqft(int start_qubit, int end_qubit);
  State& controlled_modular_exponentiation(int control_qubit, 
					   int target_start, int target_end, 
					   Bitstring a, Bitstring N, 
					   Bitstring power);
  double compute_probability_of_0(int j) const;
  const QuantumState& get_state() const;
  
  void display(bool sparse = false) const;
  void display_cbits() const;
  
  /**
   * @brief Measures the jth qubit, collapses the state, and stores the result.
   * 
   * @param j Qubit index to measure.
   * @param cbit_index Index of classical register to store the result.
   * @return Reference to the updated State object.
   */
  State& measure(int j, unsigned long cbit_index);
  bool is_initialized() const { return true; };
  int get_num_qubits() const { return num_qubits_; };
  std::vector<int> get_cbits() const { return cbits_; }
  /**
   * @brief Sets the state to a single basis vector. Clears previous entries.
   */
  void set_basis_state(Bitstring b, ComplexNumber a) {
    state_.clear(); // Clear the vector representation
    state_.push_back({b, a});
  }
    
  /**
   * @brief Sets the state from a new vector representation. 
   * Assumes the input vector already handles amplitude summation/uniqueness.
   */
  void set_superposition(const QuantumState& new_state) {
    state_ = new_state;
  }

  /**
   * @brief Retrieves the amplitude associated with a specific bitstring key.
   * Requires linear search through the vector.
   */
  ComplexNumber get_amplitude(Bitstring b) const {
    for (const auto& pair : state_) {
      if (pair.first == b) {
	return pair.second;
      }
    }
    // If the basis state is not explicitly listed, its amplitude is zero (sparse representation).
    return 0.0; 
  }
  /**
     * Sets the amplitude associated with a specific basis state (bitstring).
     * If the bitstring already exists, its amplitude is updated.
     * If it does not exist, a new entry is added.
     */
  void set_amplitude(Bitstring b, ComplexNumber a) {
        
    // Iterate through the vector to check if the bitstring already exists.
    for (auto& pair : state_) {
      if (pair.first == b) {
	// Found the existing basis state: update the amplitude (the value in the key-value pair)
	pair.second = a;
	return;
      }
    }
        
    // If the loop completes, the bitstring was not found, so add a new entry.
    // This corresponds to introducing a new basis state to the superposition.
    state_.push_back({b, a});
  }
  /**
   * @brief Retrieves the value (0 or 1) stored in the classical register 
   *        at the specified index.
   * 
   * This method accesses the result of a measurement stored in a 
   * classical register [1].
   * 
   * @param j The index of the classical bit (0-indexed).
   * @return int The value of the classical bit (0 or 1).
   * @throws std::out_of_range if the index j is invalid.
   */
  int get_cbit(size_t j) const {
    // Accesses the indexed classical register for retrieval, analogous to s.cbits[j] [2]
    // Using .at() ensures bounds checking.
    if (j < 0 || j >= cbits_.size()) {
      throw std::out_of_range("Classical bit index out of bounds.");
    }
    return cbits_.at(j);
  }
  /*
   */
  State& apply_U0_perp();
  State& grover_diffusion_Us();
  State& grover_oracle_Uf(Bitstring solution_w);
  State& run_shor_algorithm_quantum_part(Bitstring N, Bitstring a);  
};



