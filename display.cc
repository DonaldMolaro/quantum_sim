/*
 * So along with "measure" this is probably the most intesting part of the
 * simulator. This lets someone peek inside the state of a quantum machine
 * and see what is going on "inside". For Quantum obvious reasons one cannot
 * do that on a *real* quantum computer but if one is seeking to understand
 * what/how a quantum machine does what is does this is the intesting bit.
 *
 * Now the *actual* implementation of the Qubits is a sparce map, but for the
 * purposes of understanding the display method prints out everything, which
 * for me helps.
 */
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm> // Required for max/min if needed, or simply for utility
#include "state.hh"

// --- Helper Functions for Display (Refined for Alignment) ---

/** Converts the Bitstring (ULL) to a padded binary string based on N qubits. */
std::string bitstring_to_string(Bitstring b, int N) {
  std::string s;
  for (int i = N - 1; i >= 0; --i) {
    s += ((b >> i) & 1) ? '1' : '0';
  }
  return s;
}

/** Converts the ComplexNumber to a readable string (e.g., 0.707 + 0.707i). */
std::string complex_to_string(const ComplexNumber& a, int precision = 4) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision);
    
  const double real_part = a.real();
  const double imag_part = a.imag();
  const double TOLERANCE = 1e-9;
    
  if (std::abs(real_part) < TOLERANCE && std::abs(imag_part) < TOLERANCE) {
    return "0.0";
  }

  // Print real part if significant
  if (std::abs(real_part) > TOLERANCE) {
    ss << real_part;
  }
    
  // Print imaginary part if significant
  if (std::abs(imag_part) > TOLERANCE) {
    if (std::abs(real_part) > TOLERANCE) {
      // Include sign if real part is also present
      ss << (imag_part > 0 ? " + " : " - ");
    } else if (imag_part < 0) {
      // Handle negative sign for purely imaginary
      ss << "-";
    }
        
    double abs_imag = std::abs(imag_part);
    // Only print coefficient if it's not close to 1
    if (std::abs(abs_imag - 1.0) > TOLERANCE) {
      ss << abs_imag;
    }
    ss << "i";
  }
    
  // Fallback for purely imaginary parts where coefficient is 1
  if (ss.str().empty()) {
    if (std::abs(imag_part - 1.0) < TOLERANCE) return "i";
    if (std::abs(imag_part + 1.0) < TOLERANCE) return "-i";
  }
    
  return ss.str();
}

/** Calculates the probability, |a|^2. */
double amplitude_to_probability(const ComplexNumber& a) {
  // std::norm calculates |a|^2
  return std::norm(a); 
}

// --- State Class Definition (Partial, showing the new method) ---

void State::display(bool sparse) const
{
  const int N = get_num_qubits();
  const unsigned long long N_STATES = 1ULL << N;
  double total_probability = 0.0;
        
  // 1. Create a dense lookup table for quick access to amplitudes
  std::map<Bitstring, ComplexNumber> amplitude_map;
  for (const auto& pair : get_state()) {
    amplitude_map[pair.first] = pair.second;
  }

  // --- Formatting Setup ---
  // Bitstring column width N
  const int BIT_WIDTH = N + 2; 
  // Amplitude string formatting width
  const int AMP_WIDTH = 25; 
  // Probability column width
  const int PROB_WIDTH = 12;

  std::cout << "\n======================================================\n";
  std::cout << "Quantum State (N=" << N << " qubits, 2^N=" << N_STATES << " states):\n";
        
  // --- Print Header ---
  std::cout << std::left << std::setw(BIT_WIDTH) << "Bitstring";
  std::cout << std::left << std::setw(AMP_WIDTH) << "Amplitude (a + bi)";
  std::cout << std::right << std::setw(PROB_WIDTH) << "Probability (|a|^2)\n";
  std::cout << "------------------------------------------------------\n";

  // 2. Iterate through all 2^N possible bitstrings (j = 0 to 2^N - 1)
  for (Bitstring j = 0; j < N_STATES; ++j) {
    ComplexNumber amplitude;
            
    // Check if the amplitude exists in our sparse map; default to zero if not found
    if (amplitude_map.count(j)) {
      amplitude = amplitude_map.at(j);
    } else {
      amplitude = ComplexNumber(0.0, 0.0);
    }

    if ((amplitude != ComplexNumber(0.0, 0.0)) || (sparse == true))
      {
	double probability = amplitude_to_probability(amplitude);
	total_probability += probability;
            
	// Format Bitstring (Key)
	std::cout << std::left << std::setw(BIT_WIDTH) 
		  << ("|" + bitstring_to_string(j, N) + ">"); 
	
	// Format Amplitude (Value)
	std::cout << std::left << std::setw(AMP_WIDTH) 
		  << complex_to_string(amplitude, 6);
            
	// Format Probability (|a|^2)
	std::cout << std::right << std::setw(PROB_WIDTH) << std::fixed << std::setprecision(6) 
		  << probability << "\n";
      }
  }

  // --- Print Footer ---
  std::cout << "------------------------------------------------------\n";
  // Show normalization check: sum of probabilities must be 1
  std::cout << "Total Probability Sum: " << std::fixed << std::setprecision(6) << total_probability << "\n";
  // Show classical registers (if desired)
  display_cbits();
  std::cout << "======================================================\n";
}

/**
 * @brief Displays the classical register content as a bitstring (MSB -> LSB)
 *        and calculates its corresponding decimal value.
 *
 * Assumes cbits holds the LSB (little endian convention).
 */
void State::display_cbits() const
{
  if (cbits_.empty()) {
    std::cout << "Classical register is empty.\n";
    return;
  }

  // --- 1. Display as Bitstring (MSB -> LSB) ---

  std::string bitstring = "";
  // Iterate backward to print MSB first (index N-1 to 0)
  for (int i = cbits_.size() - 1; i >= 0; --i) {
    bitstring += std::to_string(cbits_[i]);
  }

  std::cout << "Classical Register Bitstring (MSB -> LSB): " << bitstring << "\n";

  // --- 2. Calculate and Display as Decimal Number ---

  // Use unsigned long long to accommodate large register sizes
  unsigned long long decimal_value = 0;
  unsigned long long power_of_two = 1;

  // Iterate forward (index 0 to N-1) to calculate decimal value:
  // cbits_[j] is the coefficient for 2^j.
  for (size_t j = 0; j < cbits_.size(); ++j) {
    if (cbits_[j] == 1) {
      // If the bit is set, add its corresponding power of 2
      decimal_value += power_of_two;
    }
    // Prepare for the next bit position (2^j+1)
    power_of_two *= 2;
  }
  std::cout << "Classical Register Decimal Number: " << decimal_value << "\n";
}
