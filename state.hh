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
public:
  /** Constructor: Initializes the state to the ground state |00...0>. */
  State(int N, int num_cbits = 0);
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
  double compute_probability_of_0(int j) const;
  const QuantumState& get_state() const;
  
  void display() const;
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
  /*
   */
  State& apply_U0_perp();
  State& grover_diffusion_Us();
  State& grover_oracle_Uf(Bitstring solution_w);
};



