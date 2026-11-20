#include "cli/viz_export.hh"
#include "cli/commands.hh"
#include "internal/limits.hh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace viz {
namespace {

std::string json_escape(const std::string& s)
{
  std::string out;
  for (char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) out += ' ';
        else out += c;
    }
  }
  return out;
}

/** Number of leading qubit-index arguments for each supported gate. */
const std::map<std::string, int>& qubit_arity()
{
  static const std::map<std::string, int> table = {
    {"H", 1}, {"X", 1}, {"Y", 1}, {"Z", 1}, {"S", 1}, {"T", 1}, {"SDG", 1}, {"TDG", 1},
    {"P", 1}, {"RX", 1}, {"RY", 1}, {"RZ", 1}, {"RU", 1}, {"RESET", 1},
    {"CX", 2}, {"CNOT", 2}, {"CZ", 2}, {"CY", 2}, {"CH", 2}, {"CRX", 2}, {"CRY", 2},
    {"CRZ", 2}, {"CP", 2}, {"CU", 2}, {"SWAP", 2}, {"ISWAP", 2}, {"XX", 2}, {"YY", 2}, {"ZZ", 2},
    {"CCX", 3}, {"TOFFOLI", 3}, {"CSWAP", 3}, {"FREDKIN", 3},
  };
  return table;
}

bool parse_int(const std::string& s, int& v)
{
  try {
    size_t pos = 0;
    v = std::stoi(s, &pos);
    return pos == s.size();
  } catch (const std::exception&) {
    return false;
  }
}

std::string number(double v)
{
  if (std::abs(v) < 1e-12) return "0";
  std::ostringstream ss;
  ss << std::setprecision(9) << v;
  return ss.str();
}

} // namespace

bool op_to_json(const std::vector<std::string>& tokens, std::string& out, std::vector<int>* qubits_out)
{
  if (tokens.empty()) return false;
  const std::string cmd = cli::upper_copy(tokens[0]);
  std::vector<int> qubits;
  std::vector<std::string> params;
  int cbit = -1;

  if (cmd == "INIT") {
    int n = 0, c = 0;
    if (tokens.size() < 2 || !parse_int(tokens[1], n)) return false;
    if (tokens.size() >= 3 && !parse_int(tokens[2], c)) return false;
    std::ostringstream ss;
    ss << "{\"name\":\"INIT\",\"qubits\":[],\"params\":[" << n << "," << c
       << "],\"cbit\":null,\"raw\":\"" << json_escape(tokens[0]) << " " << n << " " << c << "\"}";
    out = ss.str();
    return true;
  }

  if (cmd == "MCX") {
    // MCX c1 ... cn t: every token is a qubit; the last one is the target.
    if (tokens.size() < 3) return false;
    for (size_t i = 1; i < tokens.size(); ++i) {
      int q = 0;
      if (!parse_int(tokens[i], q)) return false;
      qubits.push_back(q);
    }
  } else if (cmd == "MEASURE") {
    int q = 0;
    if (tokens.size() < 2 || !parse_int(tokens[1], q)) return false;
    qubits.push_back(q);
    if (tokens.size() >= 3 && !parse_int(tokens[2], cbit)) return false;
  } else {
    const auto it = qubit_arity().find(cmd);
    if (it == qubit_arity().end()) return false;
    const size_t arity = static_cast<size_t>(it->second);
    if (tokens.size() < 1 + arity) return false;
    for (size_t i = 0; i < arity; ++i) {
      int q = 0;
      if (!parse_int(tokens[1 + i], q)) return false;
      qubits.push_back(q);
    }
    for (size_t i = 1 + arity; i < tokens.size(); ++i) params.push_back(tokens[i]);
  }

  std::string raw;
  for (size_t i = 0; i < tokens.size(); ++i) {
    if (i > 0) raw += ' ';
    raw += tokens[i];
  }

  std::ostringstream ss;
  ss << "{\"name\":\"" << json_escape(cmd) << "\",\"qubits\":[";
  for (size_t i = 0; i < qubits.size(); ++i) ss << (i ? "," : "") << qubits[i];
  ss << "],\"params\":[";
  for (size_t i = 0; i < params.size(); ++i) ss << (i ? "," : "") << "\"" << json_escape(params[i]) << "\"";
  ss << "],\"cbit\":";
  if (cbit >= 0) ss << cbit; else ss << "null";
  ss << ",\"raw\":\"" << json_escape(raw) << "\"}";
  out = ss.str();
  if (qubits_out) *qubits_out = qubits;
  return true;
}

std::string algorithm_op_json(const std::string& name, const std::string& short_label,
                              const std::string& detail, int first_qubit, int last_qubit)
{
  std::ostringstream ss;
  ss << "{\"name\":\"" << json_escape(name) << "\",\"alg\":true,\"short\":\"" << json_escape(short_label)
     << "\",\"qubits\":[";
  for (int q = first_qubit; q <= last_qubit; ++q) ss << (q > first_qubit ? "," : "") << q;
  ss << "],\"params\":[],\"cbit\":null,\"raw\":\"" << json_escape(name + (detail.empty() ? "" : " - " + detail))
     << "\"}";
  return ss.str();
}

std::string snapshot_to_json(const State& state)
{
  const int n = state.get_num_qubits();
  std::vector<QubitAmplitudePair> entries(state.get_state().begin(), state.get_state().end());
  std::sort(entries.begin(), entries.end(),
            [](const QubitAmplitudePair& a, const QubitAmplitudePair& b) { return a.first < b.first; });

  std::ostringstream ss;
  ss << "{\"n\":" << n << ",\"amps\":[";
  bool first = true;
  for (const auto& e : entries) {
    if (std::norm(e.second) < qsim::limits::AMPLITUDE_EPSILON * qsim::limits::AMPLITUDE_EPSILON) continue;
    ss << (first ? "" : ",") << "[" << e.first << "," << number(e.second.real()) << "," << number(e.second.imag()) << "]";
    first = false;
  }
  ss << "],\"bloch\":[";
  for (int q = 0; q < n; ++q) {
    const State::BlochVector bv = state.bloch(q);
    ss << (q ? "," : "") << "[" << number(bv.x) << "," << number(bv.y) << "," << number(bv.z) << "]";
  }
  ss << "],\"cbits\":[";
  const auto& cbits = state.get_cbits();
  for (size_t i = 0; i < cbits.size(); ++i) ss << (i ? "," : "") << cbits[i];
  ss << "]}";
  return ss.str();
}

std::string build_document(int num_qubits, int num_cbits, const std::vector<std::string>& steps)
{
  std::ostringstream ss;
  ss << "{\"format\":\"qsim-viz\",\"version\":1,\"num_qubits\":" << num_qubits
     << ",\"num_cbits\":" << num_cbits << ",\"steps\":[";
  for (size_t i = 0; i < steps.size(); ++i) ss << (i ? ",\n" : "\n") << steps[i];
  ss << "\n]}\n";
  return ss.str();
}

} // namespace viz
