/*
 * JSON export support for the browser viewer (viewer/viewer.html).
 *
 * The shell records one "step" per circuit-relevant command: the parsed gate
 * plus a snapshot of the state right after it. EXPORT writes those steps out
 * so the viewer can draw the circuit and let a learner scrub through it.
 */
#pragma once

#include "state.hh"
#include <string>
#include <vector>

namespace viz {

/** Largest register for which per-step snapshots are recorded. */
constexpr int MAX_SNAPSHOT_QUBITS = 10;

/** Larger limit for algorithm recording (Shor on N=15 needs 12 qubits). */
constexpr int MAX_ALGORITHM_SNAPSHOT_QUBITS = 12;

/**
 * Amplitudes kept per snapshot (largest probabilities first). The viewer only lists the
 * top few dozen, and keeping every amplitude of a 12-qubit state makes exports multi-MB.
 * The omitted count and probability are recorded so totals stay exact.
 */
constexpr size_t MAX_SNAPSHOT_AMPS = 256;

/** Cap on recorded steps so a long LOAD cannot grow memory without bound. */
constexpr size_t MAX_STEPS = 1000;

/**
 * Parses a shell command into a JSON op object such as
 * {"name":"CX","qubits":[0,1],"params":[],"cbit":null,"raw":"CX 0 1"}.
 * Returns false (leaving `out` untouched) for commands that are not gates,
 * measurements or INIT.
 */
bool op_to_json(const std::vector<std::string>& tokens, std::string& out,
                std::vector<int>* qubits_out = nullptr);

/**
 * JSON op for an algorithm step (e.g. a Grover oracle or an inverse QFT), shown in the
 * viewer as one block across qubits first_qubit..last_qubit. `short_label` must fit a narrow box.
 */
std::string algorithm_op_json(const std::string& name, const std::string& short_label,
                              const std::string& detail, int first_qubit, int last_qubit);

/**
 * JSON snapshot: {"n":..,"total":nonzero count,"p_rest":omitted probability,
 * "amps":[[idx,re,im],...],"bloch":[[x,y,z],...],"cbits":[...]}.
 */
std::string snapshot_to_json(const State& state);

/** Assembles the full export document from recorded step objects. */
std::string build_document(int num_qubits, int num_cbits, const std::vector<std::string>& steps);

} // namespace viz
