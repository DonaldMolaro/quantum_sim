#include "cli/shell.hh"
#include "cli/commands.hh"
#include "algorithms/api/grover_api.hh"
#include "demos/anneal_demo.hh"
#include "demos/bernstein_vazirani_demo.hh"
#include "demos/deutsch_jozsa_demo.hh"
#include "demos/grover_demo.hh"
#include "demos/latin_demo.hh"
#include "demos/qaoa_demo.hh"
#include "demos/qubo_demo.hh"
#include "demos/shor_demo.hh"
#include "demos/tsp_demo.hh"
#include "demos/vqa_demo.hh"
#include "demos/vqe_demo.hh"
#include "cli/shell_detail.hh"
#include <cmath>
#include <functional>
#include <iostream>
#include <unordered_map>

bool QuantumShell::handle_algorithm_commands(const std::vector<std::string>& tokens, const std::string& cmd)
{
  enum class AlgorithmCommand {
    DeutschJozsa,
    BernsteinVazirani,
    Qubo,
    Vqa,
    Qaoa,
    Vqe,
    Anneal,
    Latin,
    Shor,
    Grover,
    Tsp
  };

  static const std::unordered_map<std::string, AlgorithmCommand> kCommandMap = {
      {"DEUTSCH_JOZSA", AlgorithmCommand::DeutschJozsa},
      {"DEUTSCH", AlgorithmCommand::DeutschJozsa},
      {"BV", AlgorithmCommand::BernsteinVazirani},
      {"BERNSTEIN_VAZIRANI", AlgorithmCommand::BernsteinVazirani},
      {"QUBO", AlgorithmCommand::Qubo},
      {"VQA", AlgorithmCommand::Vqa},
      {"QAOA", AlgorithmCommand::Qaoa},
      {"VQE", AlgorithmCommand::Vqe},
      {"ANNEAL", AlgorithmCommand::Anneal},
      {"LATIN", AlgorithmCommand::Latin},
      {"SHOR", AlgorithmCommand::Shor},
      {"GROVER", AlgorithmCommand::Grover},
      {"TSP", AlgorithmCommand::Tsp},
  };

  std::unordered_map<std::string, AlgorithmCommand>::const_iterator it = kCommandMap.find(cmd);
  if (it == kCommandMap.end()) {
    return false;
  }

  switch (it->second) {
  case AlgorithmCommand::DeutschJozsa: {
    int n_inputs = cli::get_arg(tokens, 1, cmd);
    if (n_inputs == -1) return true;
    if (tokens.size() < 3) {
      std::cerr << "Error: " << cmd << " requires an oracle (CONST0, CONST1, BALANCED_XOR0, BALANCED_PARITY).\n";
      return true;
    }
    run_deutsch_jozsa_demo(n_inputs, tokens[2]);
    return true;
  }

  case AlgorithmCommand::BernsteinVazirani: {
    int n_inputs = cli::get_arg(tokens, 1, cmd);
    int secret = cli::get_arg(tokens, 2, cmd);
    int bias = (tokens.size() > 3) ? cli::get_arg(tokens, 3, cmd) : 0;
    if (n_inputs == -1 || secret == -1 || bias == -1) return true;
    run_bernstein_vazirani_demo(n_inputs, static_cast<Bitstring>(secret), bias);
    return true;
  }

  case AlgorithmCommand::Qubo: {
    if (tokens.size() < 2) {
      std::cerr << "Error: QUBO requires a mode (DEMO, EXACT, GROVER).\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_qubo_demo();
      return true;
    }
    if (mode == "EXACT") {
      int n = cli::get_arg(tokens, 2, "QUBO EXACT");
      if (n == -1) return true;
      std::vector<double> matrix;
      if (!shell_detail::parse_square_matrix(tokens, 3, n, "QUBO EXACT", matrix)) return true;
      run_qubo_exact_cli(n, matrix);
      return true;
    }
    if (mode == "GROVER") {
      int n = cli::get_arg(tokens, 2, "QUBO GROVER");
      double threshold = cli::get_double_arg(tokens, 3, "QUBO GROVER");
      int iterations = cli::get_arg(tokens, 4, "QUBO GROVER");
      if (n == -1 || iterations == -1 || std::isnan(threshold)) return true;
      std::vector<double> matrix;
      if (!shell_detail::parse_square_matrix(tokens, 5, n, "QUBO GROVER", matrix)) return true;
      run_qubo_grover_cli(n, threshold, iterations, matrix);
      return true;
    }
    std::cerr << "Error: QUBO mode must be DEMO, EXACT, or GROVER.\n";
    return true;
  }

  case AlgorithmCommand::Vqa: {
    if (tokens.size() < 2) {
      std::cerr << "Error: VQA requires a mode (DEMO, QAOA).\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_vqa_demo();
      return true;
    }
    if (mode == "QAOA") {
      int n = cli::get_arg(tokens, 2, "VQA QAOA");
      int p_layers = cli::get_arg(tokens, 3, "VQA QAOA");
      int shots = cli::get_arg(tokens, 4, "VQA QAOA");
      int iters = cli::get_arg(tokens, 5, "VQA QAOA");
      double step = cli::get_double_arg(tokens, 6, "VQA QAOA");
      if (n == -1 || p_layers == -1 || shots == -1 || iters == -1 || std::isnan(step)) return true;
      std::vector<double> matrix;
      if (!shell_detail::parse_square_matrix(tokens, 7, n, "VQA QAOA", matrix)) return true;
      run_vqa_qaoa_cli(n, p_layers, shots, iters, step, matrix);
      return true;
    }
    std::cerr << "Error: VQA mode must be DEMO or QAOA.\n";
    return true;
  }

  case AlgorithmCommand::Qaoa: {
    if (tokens.size() < 2) {
      std::cerr << "Error: QAOA requires a mode (DEMO, QUBO).\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_qaoa_demo();
      return true;
    }
    if (mode == "QUBO") {
      int n = cli::get_arg(tokens, 2, "QAOA QUBO");
      int p_layers = cli::get_arg(tokens, 3, "QAOA QUBO");
      int shots = cli::get_arg(tokens, 4, "QAOA QUBO");
      int iters = cli::get_arg(tokens, 5, "QAOA QUBO");
      double step = cli::get_double_arg(tokens, 6, "QAOA QUBO");
      if (n == -1 || p_layers == -1 || shots == -1 || iters == -1 || std::isnan(step)) return true;
      std::vector<double> matrix;
      if (!shell_detail::parse_square_matrix(tokens, 7, n, "QAOA QUBO", matrix)) return true;
      run_qaoa_qubo_cli(n, p_layers, shots, iters, step, matrix);
      return true;
    }
    std::cerr << "Error: QAOA mode must be DEMO or QUBO.\n";
    return true;
  }

  case AlgorithmCommand::Vqe: {
    if (tokens.size() < 2) {
      std::cerr << "Error: VQE requires a mode (DEMO, RUN).\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_vqe_demo();
      return true;
    }
    if (mode == "RUN") {
      int n = cli::get_arg(tokens, 2, "VQE RUN");
      int layers = cli::get_arg(tokens, 3, "VQE RUN");
      int iters = cli::get_arg(tokens, 4, "VQE RUN");
      double step = cli::get_double_arg(tokens, 5, "VQE RUN");
      int shots = cli::get_arg(tokens, 6, "VQE RUN");
      int term_count = cli::get_arg(tokens, 7, "VQE RUN");
      if (n == -1 || layers == -1 || iters == -1 || shots == -1 || term_count == -1 || std::isnan(step)) {
        return true;
      }
      if (term_count <= 0) {
        std::cerr << "Error: VQE RUN requires term_count > 0.\n";
        return true;
      }

      VqeHamiltonian h;
      h.n_qubits = n;
      size_t idx = 8;
      for (int t = 0; t < term_count; ++t) {
        double coeff = cli::get_double_arg(tokens, idx, "VQE RUN");
        if (std::isnan(coeff)) return true;
        ++idx;
        int op_count = cli::get_arg(tokens, idx, "VQE RUN");
        if (op_count == -1) return true;
        if (op_count < 0) {
          std::cerr << "Error: VQE RUN term op_count must be >= 0.\n";
          return true;
        }
        ++idx;
        VqePauliTerm term;
        term.coeff = coeff;
        for (int k = 0; k < op_count; ++k) {
          if (idx >= tokens.size()) {
            std::cerr << "Error: VQE RUN missing Pauli op token.\n";
            return true;
          }
          const std::string op_token = tokens[idx++];
          if (op_token.size() != 1) {
            std::cerr << "Error: VQE RUN Pauli op must be one of X,Y,Z.\n";
            return true;
          }
          int q = cli::get_arg(tokens, idx, "VQE RUN");
          if (q == -1) return true;
          ++idx;
          term.ops.push_back(VqePauliOp{op_token[0], q});
        }
        h.terms.push_back(term);
      }

      if (idx != tokens.size()) {
        std::cerr << "Error: VQE RUN has extra trailing tokens.\n";
        return true;
      }
      run_vqe_cli(h, layers, iters, step, shots);
      return true;
    }
    std::cerr << "Error: VQE mode must be DEMO or RUN.\n";
    return true;
  }

  case AlgorithmCommand::Anneal: {
    if (tokens.size() < 2) {
      std::cerr << "Error: ANNEAL requires a mode (DEMO, QUBO).\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_anneal_demo();
      return true;
    }
    if (mode == "QUBO") {
      if (tokens.size() < 9) {
        std::cerr << "Error: ANNEAL QUBO requires method, n, steps, sweeps, beta_start, beta_end, replicas, matrix.\n";
        return true;
      }
      const std::string method = tokens[2];
      int n = cli::get_arg(tokens, 3, "ANNEAL QUBO");
      int steps = cli::get_arg(tokens, 4, "ANNEAL QUBO");
      int sweeps = cli::get_arg(tokens, 5, "ANNEAL QUBO");
      double beta_start = cli::get_double_arg(tokens, 6, "ANNEAL QUBO");
      double beta_end = cli::get_double_arg(tokens, 7, "ANNEAL QUBO");
      int replicas = cli::get_arg(tokens, 8, "ANNEAL QUBO");
      if (n == -1 || steps == -1 || sweeps == -1 || replicas == -1 || std::isnan(beta_start) || std::isnan(beta_end)) {
        return true;
      }
      std::vector<double> matrix;
      if (!shell_detail::parse_square_matrix(tokens, 9, n, "ANNEAL QUBO", matrix)) return true;
      run_anneal_qubo_cli(method, n, steps, sweeps, beta_start, beta_end, replicas, matrix);
      return true;
    }
    std::cerr << "Error: ANNEAL mode must be DEMO or QUBO.\n";
    return true;
  }

  case AlgorithmCommand::Latin: {
    size_t idx = 1;
    std::string mode = "demo";
    if (idx < tokens.size() && (tokens[idx] == "DEMO" || tokens[idx] == "COUNT" || tokens[idx] == "PRINT-ALL")) {
      mode = tokens[idx];
      ++idx;
    }
    int iters = -1;
    std::vector<int> row0;
    if (idx < tokens.size()) {
      int maybe_iters = cli::get_arg(tokens, idx, "LATIN");
      if (maybe_iters == -1) return true;
      if (mode == "DEMO") {
        iters = maybe_iters;
        ++idx;
      }
    }
    while (idx < tokens.size()) {
      int v = cli::get_arg(tokens, idx, "LATIN");
      if (v == -1) return true;
      row0.push_back(v);
      ++idx;
    }
    int row_vals[3] = {0, 1, 2};
    if (!row0.empty()) {
      if (row0.size() != 3) {
        std::cerr << "LATIN row0 must have exactly 3 values.\n";
        return true;
      }
      row_vals[0] = row0[0];
      row_vals[1] = row0[1];
      row_vals[2] = row0[2];
    }
    if (mode == "COUNT") {
      run_latin3_count_row0(row_vals);
      return true;
    }
    if (mode == "PRINT-ALL") {
      run_latin3_print_all_row0(row_vals);
      return true;
    }
    if (!row0.empty() || iters >= 0) run_latin3_grover_demo_row0(row_vals, iters);
    else run_latin3_grover_demo(-1);
    return true;
  }

  case AlgorithmCommand::Shor: {
    int N = cli::get_arg(tokens, 1, "SHOR");
    if (N == -1) return true;
    run_shor_demo(static_cast<Bitstring>(N));
    return true;
  }

  case AlgorithmCommand::Grover: {
    if (tokens.size() < 2) {
      std::cerr << "Error: GROVER requires at least one target.\n";
      return true;
    }
    std::vector<Bitstring> targets;
    targets.reserve(tokens.size() - 1);
    for (size_t i = 1; i < tokens.size(); ++i) {
      int t = cli::get_arg(tokens, i, "GROVER");
      if (t == -1) return true;
      if (t < 0) {
        std::cerr << "Error: GROVER targets must be >= 0.\n";
        return true;
      }
      targets.push_back(static_cast<Bitstring>(t));
    }
    int n_qubits = shell_detail::infer_qubits_from_targets(targets);
    if (state && state->is_initialized()) {
      n_qubits = state->get_num_qubits();
    }
    GroverResult result = run_grover(n_qubits, targets);
    if (!result.ok) {
      std::cerr << "Grover error: " << result.error << "\n";
      return true;
    }
    std::cout << "Grover iterations used: " << result.iterations << " (n=" << n_qubits << ")\n";
    return true;
  }

  case AlgorithmCommand::Tsp: {
    if (tokens.size() < 2) {
      std::cerr << "Error: TSP requires a mode (DEMO, EXACT).\n";
      return true;
    }
    const std::string mode = tokens[1];
    if (mode == "DEMO") {
      run_tsp_demo();
      return true;
    }
    if (mode == "EXACT") {
      int n_cities = cli::get_arg(tokens, 2, "TSP EXACT");
      double penalty = cli::get_double_arg(tokens, 3, "TSP EXACT");
      if (n_cities == -1 || std::isnan(penalty)) return true;
      std::vector<double> distance;
      if (!shell_detail::parse_square_matrix(tokens, 4, n_cities, "TSP EXACT", distance)) return true;
      run_tsp_exact_cli(n_cities, penalty, distance);
      return true;
    }
    std::cerr << "Error: TSP mode must be DEMO or EXACT.\n";
    return true;
  }
  }
  return false;
}
