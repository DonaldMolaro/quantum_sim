/*
 * Impementation of Grover's search algorithm. This is quite frankly
 * spooky when you figure out what it's doing. This implemenation is
 * ment to be verbose, printing out what it's doing at every step so
 * so using it and expecting it to be "fast" is not the point.
 *
 * It does "search" the solution space in O(root(N)) time, which for those of us
 * that are used to O(N) is well, spooky.
 *
 * Someday I'll extend the implementaion so that it can find multiple solutions
 * to the search space, but today is not that day.
 *
 * If you want to *understand* what/how Grover's algorithim works go watch this:
 *
 * from https://www.youtube.com/@3blue1brown
 *
 * https://www.youtube.com/watch?v=RQWpF2Gb-gU&t=7s
 *
 * Even if you *think* you understand Grover's algorithim, you should go watch
 * the video.
 *
 */
#include "state.hh"
#include <cmath>
#include <numeric> // For std::accumulate if needed

// Include Complex, Bitstring, State class definition (assumed from context)
// Assume State class methods grover_oracle_Uf, grover_diffusion_Us, 
// and the h(j) method are available.

void run_grover_search(State *s, Bitstring target_w)
{
    
  // 1. Initialization
  int n_qubits = s->get_num_qubits();
  if (n_qubits >= 63) {
    std::cerr << "Grover demo supports up to 62 qubits (got " << n_qubits << ").\n";
    return;
  }
  const double N = std::ldexp(1.0, n_qubits); // N = 2^n (safe for display and sqrt)

  std::cout << "Starting Grover search for solution index " << target_w 
	    << " in a database of size N=" << N << " (n=" << n_qubits << " qubits).\n";

  // Calculate optimal number of iterations R
  // R = floor( (PI/4) * sqrt(N) ) 
  const double PI = std::acos(-1.0);
  const int R = static_cast<int>(std::floor((PI / 4.0) * std::sqrt(N)));
  std::cout << "Optimal number of Grover iterations: R = " << R << "\n";
    
  // Check if R=0 is optimal (usually if N=1 or N=2)
  if (R == 0 && N > 1) {
    // For small N, we typically require R=1 iteration
    // R is kept at 1 or adjusted based on specific analysis. 
    // We will use the calculated R >= 1 if N > 1.
  }


  // 2. Apply H^n to initialize uniform superposition |s>
  // This is the starting state |ψ〉 = H|00 . . . 0〉.
  for (int j = 0; j < n_qubits; ++j) {
    s->h(j); 
  }

  std::cout << "State initialized to uniform superposition.\n";
  
  // 3. Apply the Grover Iterate G = U_psi_perp * Uf, R times
  // The Grover Iterate G is the operator that rotates the state towards the solution.
  for (int k = 0; k < R; ++k) {
    s->display();
    std::cout << "\n--- Iteration " << k + 1 << " / " << R << " ---\n";
        
    // A. Apply Oracle Uf (Phase shift)
    s->grover_oracle_Uf(target_w); 
        
    // B. Apply Diffusion Operator U_psi_perp (Inversion about the mean)
    s->grover_diffusion_Us();
  }

  std::cout << "\n--- Final State ---\n";
    
  // 4. Final Measurement / Display
  // The probability amplitude of |w〉 should now be boosted near 1.
  s->display(); // Use the display method to show amplitudes and probabilities
}


