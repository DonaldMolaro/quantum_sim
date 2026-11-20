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

/** JSON snapshot: {"n":..,"amps":[[idx,re,im],...],"bloch":[[x,y,z],...],"cbits":[...]}. */
std::string snapshot_to_json(const State& state);

/** Assembles the full export document from recorded step objects. */
std::string build_document(int num_qubits, int num_cbits, const std::vector<std::string>& steps);

} // namespace viz
