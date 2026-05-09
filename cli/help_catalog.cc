#include "cli/help_catalog.hh"

#include <ostream>

namespace cli {
namespace help_catalog {
namespace {

struct CommandHelp {
  const char* syntax;
  const char* summary;
};

const CommandHelp kGateCommands[] = {
  {"H <j>", "Hadamard gate on qubit j."},
  {"X <j>", "NOT gate on qubit j."},
  {"Y <j>", "Pauli-Y gate on qubit j."},
  {"Z <j>", "Pauli-Z gate on qubit j."},
  {"S <j>   SDG <j>", "Quarter-turn phase gate and its inverse."},
  {"T <j>   TDG <j>", "Eighth-turn phase gate and its inverse."},
  {"P <j> <phi>", "Phase gate with angle phi."},
  {"RX/RY/RZ <j> <theta>", "Single-qubit rotations. Angles accept PI, TAU, and DEG forms."},
  {"RU <j> <theta> <phi> <lambda>", "General single-qubit U gate."},
  {"CX/CNOT <c> <t>", "Controlled-X gate."},
  {"CY/CH/CZ <c> <t>", "Controlled Pauli/Hadamard gates."},
  {"CRX/CRY/CRZ <c> <t> <theta>", "Controlled rotations."},
  {"CP <c> <t> <phi>", "Controlled phase gate."},
  {"CU <c> <t> <th> <ph> <la>", "Controlled general U gate."},
  {"CCX/TOFFOLI <c1> <c2> <t>", "Double-controlled X gate."},
  {"MCX <c1> [c2 ...] <t>", "Multi-controlled X gate."},
  {"SWAP <j> <k>   ISWAP <j> <k>", "Swap gates."},
  {"CSWAP/FREDKIN <c> <a> <b>", "Controlled swap gate."},
  {"XX/YY/ZZ <j> <k> <theta>", "Two-qubit Ising interaction gates."},
};

const CommandHelp kAlgorithmCommands[] = {
  {"GROVER <targets...>", "Run Grover search on the listed basis states."},
  {"GROVER AUTO <n> <count_iters> <targets...>", "Auto-tune Grover iterations using quantum counting."},
  {"DEUTSCH_JOZSA <n> <oracle>", "Run the Deutsch-Jozsa classification demo."},
  {"BV <n> <secret> [bias]", "Run Bernstein-Vazirani."},
  {"SIMON DEMO", "Run the built-in Simon demo."},
  {"SIMON <n> <secret> [shots]", "Recover Simon's hidden xor-mask."},
  {"SHOR <N>", "Factor N with the Shor demo wrapper."},
  {"QPE DEMO", "Run the built-in phase estimation demo."},
  {"QPE <m> <phase>", "Estimate a phase using m precision qubits."},
  {"QCOUNT DEMO", "Run the built-in quantum counting demo."},
  {"QCOUNT RUN <n> <iters> <t1> [t2 ...]", "Estimate marked-state count with explicit iterations."},
  {"QCOUNT <n> <t1> [t2 ...]", "Estimate marked-state count with default iterations."},
  {"QUBO DEMO", "Run the built-in QUBO demo."},
  {"QUBO EXACT <n> <entries...>", "Solve a QUBO exactly."},
  {"QUBO GROVER <n> <threshold> <iterations> <entries...>", "Threshold search for QUBO solutions with Grover."},
  {"QAOA DEMO", "Run the QAOA demo."},
  {"QAOA QUBO <n> <p> <shots> <iters> <step> <entries...>", "Optimize a QUBO with QAOA."},
  {"VQA DEMO", "Run the VQA demo."},
  {"VQA QAOA <n> <p> <shots> <iters> <step> <entries...>", "Optimize a QUBO with the VQA wrapper."},
  {"VQE DEMO", "Run the VQE demo."},
  {"VQE RUN <n> <layers> <iters> <step> <shots> <terms> ...", "Run VQE for an explicit Hamiltonian."},
  {"ANNEAL DEMO", "Run the annealing demo."},
  {"ANNEAL QUBO <SA|SQA> <n> <steps> <sweeps> <beta_start> <beta_end> <replicas> <entries...>", "Solve a QUBO with simulated annealing."},
  {"TSP DEMO", "Run the built-in TSP demo."},
  {"TSP EXACT <n> <penalty> <distance entries...>", "Solve the fixed-start TSP instance exactly."},
  {"QEC DEMO", "Run the 3-qubit bit-flip code demo."},
  {"QEC RUN <logical_bit> <error_qubit>", "Run QEC with explicit logical value and error location."},
  {"LATIN [iters]", "Grover demo for 3x3 Latin squares."},
  {"LATIN DEMO [iters] [r0 r1 r2]", "Latin-square demo with custom first row."},
  {"LATIN COUNT [r0 r1 r2]", "Count valid Latin squares for a first-row permutation."},
  {"LATIN PRINT-ALL [r0 r1 r2]", "Print all valid Latin squares for a first-row permutation."},
};

const CommandHelp kUtilityCommands[] = {
  {"INIT <n> [c]", "Initialize n qubits and optional classical bits."},
  {"MEASURE <q> <c>", "Measure qubit q into classical bit c."},
  {"RESET <q>", "Reset qubit q to |0>."},
  {"IF <c> <gate...>", "Run a gate only if classical bit c is 1."},
  {"DISPLAY   DISPLAY ALL   DISPLAY TOP <k>", "Show the current state, dense view, or largest amplitudes."},
  {"DISPLAY PROBS   DISPLAY PROBS ALL   DISPLAY PROBS TOP <k>", "Show probability-only views for faster scanning."},
  {"CHECK <mode> ...", "Run state sanity checks."},
  {"EXPECT <P> <q> [P q ...]", "Expectation value of a Pauli product."},
  {"FIDELITY <idx>", "Fidelity with computational basis state idx."},
  {"BLOCH <q>", "Bloch vector plus a compact single-qubit visual summary."},
  {"ENTROPY <j> [k]", "Von Neumann entropy for a subsystem."},
  {"SWAP_TEST <anc> <a_start> <b_start> <n>", "Estimate overlap with the SWAP test."},
  {"LOAD <file>", "Load and execute commands from a file."},
  {"SAVE <file>", "Save command history to a file."},
  {"SHOTS <n> <file>", "Run a circuit file n times and print a histogram with bars."},
  {"NOISE <p>", "Set per-gate depolarizing noise probability."},
  {"VERBOSE <level>", "Set verbosity to QUIET, NORMAL, VERBOSE, or 0/1/2."},
  {"TUTOR <ON|OFF>", "Enable or disable teaching narration."},
  {"SEED <n>", "Set deterministic RNG seed."},
  {"QFTMODE <DIRECT|GATE>", "Choose the QFT implementation mode."},
  {"QRNG <n> [count]", "Sample quantum random values."},
  {"HELP   HELP TOPICS   HELP ALL", "Open short, indexed, or full help."},
  {"QUIT", "Exit the simulator."},
};

const CommandHelp* commands_for(Section section, int& count, const char*& title)
{
  switch (section) {
    case Section::Gates:
      count = static_cast<int>(sizeof(kGateCommands) / sizeof(kGateCommands[0]));
      title = "[Gates]";
      return kGateCommands;
    case Section::Algorithms:
      count = static_cast<int>(sizeof(kAlgorithmCommands) / sizeof(kAlgorithmCommands[0]));
      title = "[Algorithms]";
      return kAlgorithmCommands;
    case Section::Utility:
      count = static_cast<int>(sizeof(kUtilityCommands) / sizeof(kUtilityCommands[0]));
      title = "[Utility]";
      return kUtilityCommands;
  }
  count = 0;
  title = "";
  return 0;
}

void print_commands(std::ostream& out, const CommandHelp* commands, int count)
{
  for (int i = 0; i < count; ++i) {
    out << "  " << commands[i].syntax << "\n";
    out << "    " << commands[i].summary << "\n";
  }
}

} // namespace

void print_section(std::ostream& out, Section section)
{
  int count = 0;
  const char* title = "";
  const CommandHelp* commands = commands_for(section, count, title);
  out << "\n" << title << "\n";
  print_commands(out, commands, count);
}

void print_all(std::ostream& out)
{
  out << "\n--- Quantum Simulator Commands ---\n";
  print_section(out, Section::Gates);
  print_section(out, Section::Algorithms);
  print_section(out, Section::Utility);
  out << "----------------------------------\n";
}

} // namespace help_catalog
} // namespace cli
