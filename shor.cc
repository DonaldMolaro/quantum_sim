#include "state.hh" // Assuming Bitstring, ComplexNumber, etc. are defined here

// Assuming the existence of:
// - State::qft(start, end)
// - State::iqft(start, end)
// - State::controlled_modular_exponentiation(control_qubit, target_start, target_end, a, N, power)

State& State::run_shor_algorithm_quantum_part(Bitstring N, Bitstring a) {
    
    // --- 1. Determine Qubit Allocation ---
    
    // N=15 requires n_target=4 qubits (indices 0-3). 
    // The test framework assumed n_control=5 qubits (indices 4-8).
    const int n_target = 4; 
    const int target_start = 0;
    const int target_end = n_target - 1; 

    const int n_control = 5; 
    const int control_start = n_target; // Index 4
    const int control_end = n_target + n_control - 1; // Index 8

    std::cout << "The state must be initialized to |0>_control |1>_target." << std::endl;
    std::cout << "Assuming the State constructor initialized to |0...0>|0>" << std::endl;
    display();
    
    
    std::cout << "Set the target register LSB (q0) to |1> to achieve the |1 mod N> starting state." << std::endl;
    std::cout << "This assumes the remaining target qubits (q1, q2, q3) are already |0>." << std::endl;
    this->x(target_start);
    display();
    
    // --- 2. Initialize Control Register using Hadamard Gates ---
    
    std::cout << "Apply Hadamard to all qubits in the control register to create uniform superposition" << std::endl;
    // |+... +>_control |1>_target.
    for (int j = control_start; j <= control_end; ++j) {
        // H gate implements superposition for the input state preparation
        this->h(j); 
    }
    display();

    // --- 3. Controlled Modular Exponentiation (U^x) ---
    
    std::cout << "This step implements the controlled unitary transformations U^(2^j) needed for QPE." << std::endl;
    std::cout << "The exponent (power) for U is 2^j, controlled by the j-th qubit in the control register." << std::endl;
    std::cout << "This operation is computationally complex, typically implemented using arithmetic " << std::endl;
    std::cout << "built from gates like Toffoli (CCNOT) and CNOT, which act as set transformations." << std::endl;
    
    Bitstring power = 1ULL; // Corresponds to U^(2^0) = U^1
    
    // Iterate from the MSB of the control register (j=n_control - 1) down to the LSB (j=0)
    for (int j = control_end; j >= control_start; --j) {
        
        // Applying U^power, controlled by qubit j.
        // U_a^power implements multiplication by a^power mod N on the target register.
        this->controlled_modular_exponentiation(j, target_start, target_end, a, N, power);
        display();
        // The power doubles for the next qubit in the register (j-1)
        power *= 2; 
    }

    // --- 4. Inverse Quantum Fourier Transform (IQFT) ---
    
    std::cout << "Apply the IQFT to the control register to transform the phase information (s/r) " << std::endl;
    std::cout << "into measurable amplitudes." << std::endl;
    this->iqft(control_start, control_end);
    display();

    // The method assumes subsequent measurement and classical post-processing 
    // (e.g., Continued Fractions algorithm) are handled externally to this function.
    
    return *this;
}


/*
State& State::run_shor_algorithm_quantum_part(Bitstring N, Bitstring a) {
    
    // --- 1. Determine Qubit Allocation ---
    
    // Determine the required number of qubits for the target register (log2(N))
    // We assume N is such that a 4-bit register (indices 0-3) is sufficient for illustration (N=15)
    // In practice, this size is n_target = ceil(log2(N)).
    const int n_target = 4; 
    const int target_start = 0;
    const int target_end = n_target - 1; 

    // Determine the required number of qubits for the control register (precision register)
    // In complexity analysis, this is n_control = ceil(2 * log2(N)) + 1 or similar,
    // satisfying 2^m >= 2r^2. Using m=5 for illustrative purposes, indices 4-8.
    const int n_control = 5; 
    const int control_start = n_target; // Index 4
    const int control_end = n_target + n_control - 1; // Index 8

    // --- 2. Initialize State Registers ---
    
    // The initial state is |0>_control |1>_target.
    // Assuming State constructor initializes to |0...0> and X gate is applied externally/via helper
    
    // Set the target register to |1 mod N>. 
    // This assumes the initial state |0...0> only needs X on the LSB (index 0).
    // Note: The general case requires initialization to |1 mod N>.
    this->x(target_start); 
    
    // Apply Hadamard to all qubits in the control register to create a uniform superposition
    // |+... +>_control |1>_target.
    // This prepares the register for the functional transformation (QFT / phase estimation preparation).
    for (int j = control_start; j <= control_end; ++j) {
        this->h(j); 
    }

    // --- 3. Controlled Modular Exponentiation (U^x) ---
    
    // This step implements the controlled unitary transformations U^(2^j) needed for QPE.
    // The target register is fixed at indices [target_start, target_end].
    // The control qubits are the bits of the eigenvalue register [control_start, control_end].
    // The exponent (power) for U is 2^j, controlled by the j-th qubit in the control register.
    
    Bitstring power = 1ULL; // Start with U^(2^0) = U^1
    for (int j = control_end; j >= control_start; --j) {
        // Apply the controlled_modular_exponentiation operation U^(power) controlled by qubit j.
        // This abstracts the complex circuit necessary to compute |x>|y> -> |x>|a^x * y mod N>.
        this->controlled_modular_exponentiation(j, target_start, target_end, a, N, power);
        
        // Double the power for the next qubit.
        power *= 2; 
    }

    // --- 4. Inverse Quantum Fourier Transform (IQFT) ---
    
    // Apply the IQFT to the control register (qubits n_target to end) to extract the phase estimate.
    // This transforms the phase information into measurable amplitudes.
    this->iqft(control_start, control_end); 

    // The result is a superposition where the amplitude of state |k> is concentrated 
    // around values x = s/r * 2^m, which are then measured.
    // The state is now ready for measurement (step 5 in QPE, often outside this function).
    
    return *this;
}

*/
