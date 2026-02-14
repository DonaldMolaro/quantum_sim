#include "state.hh"
#include "tests/helpers.hh"
#include "algorithms/latin_square.hh"
#include "algorithms/api/shor_api.hh"
#include "algorithms/shor_classical.hh"
#include "algorithms/shor_quantum.hh"
#include "cli/commands.hh"
#include "demos/anneal_demo.hh"
#include "demos/bernstein_vazirani_demo.hh"
#include "demos/deutsch_jozsa_demo.hh"
#include "demos/grover_demo.hh"
#include "demos/latin_demo.hh"
#include "demos/qaoa_demo.hh"
#include "demos/quantum_counting_demo.hh"
#include "demos/qubo_demo.hh"
#include "demos/shor_demo.hh"
#include "demos/simon_demo.hh"
#include "demos/vqa_demo.hh"
#include "demos/vqe_demo.hh"
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

using test_helpers::assert_equal;
using test_helpers::ScopedEnv;
using test_helpers::assert_double_close;
using test_helpers::assert_exits_with_failure;

void test_display_output_paths() {
    State s(2, 2);
    s.set_amplitude(0, ComplexNumber(1.0, 0.0));
    s.set_amplitude(1, ComplexNumber(0.0, 1.0));
    s.set_amplitude(2, ComplexNumber(0.0, -1.0));
    s.set_amplitude(3, ComplexNumber(2.0, 3.0));

    s.display(false);
    s.display(true);

    s.measure(0, 0);
    s.measure(1, 1);
    s.display_cbits();

    State empty_cbits(1, 0);
    empty_cbits.display_cbits();

    State sparse_only(2, 0);
    sparse_only.set_amplitude(0, ComplexNumber(1.0, 0.0));
    sparse_only.display(true);

    State one_bit(1, 1);
    one_bit.set_basis_state(1, ONE_COMPLEX);
    one_bit.measure_with_rng(0, 0, 1.0);
    one_bit.display_cbits();
}

void test_qft_invalid_ranges() {
    State s(2, 0);
    s.set_qft_mode(State::QftMode::Direct);
    s.qft(1, 0);
    s.iqft(1, 0);

    s.set_qft_mode(State::QftMode::Gate);
    s.qft(1, 0);
    s.iqft(1, 0);
}

void test_qft_tiny_amplitude_continue() {
    State s(2, 0);
    s.set_qft_mode(State::QftMode::Direct);
    s.set_amplitude(0, ComplexNumber(1e-12, 0.0));
    s.qft(0, 1);

    State s2(2, 0);
    s2.set_qft_mode(State::QftMode::Direct);
    s2.set_amplitude(0, ComplexNumber(1e-12, 0.0));
    s2.iqft(0, 1);
}

void test_measure_branches() {
    State s(1, 1);
    s.set_basis_state(1, ONE_COMPLEX);
    s.measure_with_rng(0, 0, 0.0);
    s.measure(0, 0);

    State no_cbits(1, 0);
    no_cbits.set_basis_state(1, ONE_COMPLEX);
    std::vector<int> out;
    std::vector<double> rv(1, 0.8);
    no_cbits.measure_all_with_rng(rv, out);

    State tiny(1, 1);
    tiny.set_amplitude(0, ComplexNumber(0.0, 0.0));
    tiny.measure_with_rng(0, 0, 0.0);
}

void test_cli_command_parsing()
{
    std::vector<std::string> tokens = cli::parse_command("  INIT  3   2 ");
    if (tokens.size() != 3 || tokens[0] != "INIT" || tokens[1] != "3" || tokens[2] != "2") {
        throw std::runtime_error("CLI parse_command failed to tokenize INIT");
    }

    std::vector<std::string> empty = cli::parse_command("   ");
    if (!empty.empty()) {
        throw std::runtime_error("CLI parse_command should return empty for whitespace");
    }

    std::vector<std::string> angle_tokens = {"RX", "0", "PI/2"};
    double theta = cli::get_angle_arg_required(angle_tokens, 2, "RX");
    assert_double_close(theta, std::acos(-1.0) / 2.0, 1e-9, "CLI PI/2 parsing");

    std::vector<std::string> deg_tokens = {"RY", "0", "90", "DEG"};
    double theta_deg = cli::get_angle_arg_required(deg_tokens, 2, "RY");
    assert_double_close(theta_deg, std::acos(-1.0) / 2.0, 1e-9, "CLI DEG parsing");

    std::vector<std::string> missing = {"X"};
    int arg = cli::get_arg(missing, 1, "X");
    if (arg != -1) {
        throw std::runtime_error("CLI get_arg should return -1 on missing argument");
    }

    std::vector<std::string> qcount_run = cli::parse_command("QCOUNT RUN 3 6 1 6");
    if (qcount_run.size() != 6 || qcount_run[0] != "QCOUNT" || qcount_run[1] != "RUN") {
        throw std::runtime_error("CLI parse_command failed on QCOUNT RUN");
    }

    std::vector<std::string> grover_auto = cli::parse_command("GROVER AUTO 3 6 1 6");
    if (grover_auto.size() != 6 || grover_auto[1] != "AUTO") {
        throw std::runtime_error("CLI parse_command failed on GROVER AUTO");
    }

    std::vector<std::string> simon = cli::parse_command("SIMON 4 10 8");
    if (simon.size() != 4 || simon[0] != "SIMON") {
        throw std::runtime_error("CLI parse_command failed on SIMON");
    }

    std::vector<std::string> vqa = cli::parse_command("VQA QAOA 3 1 0 40 0.25 -2 0 2 0 1 0 2 0 -3");
    if (vqa.size() != 16 || vqa[0] != "VQA" || vqa[1] != "QAOA") {
        throw std::runtime_error("CLI parse_command failed on VQA QAOA transcript");
    }

    std::vector<std::string> vqe = cli::parse_command("VQE RUN 2 2 50 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1");
    if (vqe.size() != 22 || vqe[0] != "VQE" || vqe[1] != "RUN") {
        throw std::runtime_error("CLI parse_command failed on VQE transcript");
    }

    std::vector<std::string> qubo = cli::parse_command("QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3");
    if (qubo.size() != 12 || qubo[0] != "QUBO" || qubo[1] != "EXACT") {
        throw std::runtime_error("CLI parse_command failed on QUBO EXACT transcript");
    }

    std::vector<std::string> seed = cli::parse_command("SEED 1234");
    if (seed.size() != 2 || seed[0] != "SEED" || seed[1] != "1234") {
        throw std::runtime_error("CLI parse_command failed on SEED command");
    }

    std::vector<std::string> check = cli::parse_command("CHECK TARGETS 3 5");
    if (check.size() != 4 || check[0] != "CHECK" || check[1] != "TARGETS") {
        throw std::runtime_error("CLI parse_command failed on CHECK command");
    }
}

void test_latin_square_validations() {
    const int row0[3] = {0, 1, 2};
    Bitstring assignment = 0;
    assignment |= (1ULL << 0);  // cell 0 = 1
    assignment |= (2ULL << 2);  // cell 1 = 2
    assignment |= (0ULL << 4);  // cell 2 = 0
    assignment |= (2ULL << 6);  // cell 3 = 2
    assignment |= (0ULL << 8);  // cell 4 = 0
    assignment |= (1ULL << 10); // cell 5 = 1

    if (!is_valid_latin3(assignment, row0)) {
        throw std::runtime_error("Expected assignment to be valid latin square");
    }

    if (is_valid_latin3_fixedrow(0)) {
        throw std::runtime_error("Expected zero assignment to be invalid latin square");
    }
}

void test_latin_square_demo_forced_measure() {
    const int row0[3] = {0, 1, 2};
    Bitstring assignment = 0;
    assignment |= (1ULL << 0);
    assignment |= (2ULL << 2);
    assignment |= (0ULL << 4);
    assignment |= (2ULL << 6);
    assignment |= (0ULL << 8);
    assignment |= (1ULL << 10);

    {
        ScopedEnv env("QSIM_LATIN_FORCE_MEASURED", std::to_string(assignment));
        run_latin3_grover_demo_row0(row0, 0);
    }
    {
        ScopedEnv env("QSIM_LATIN_FORCE_MEASURED", "0");
        run_latin3_grover_demo_row0(row0, 0);
    }
    {
        ScopedEnv env("QSIM_LATIN_FORCE_M", "0");
        run_latin3_grover_demo_row0(row0, 0);
    }
}

void test_latin_square_invalid_row0() {
    const int bad_range[3] = {0, 1, 3};
    const int bad_perm[3] = {0, 0, 1};
    assert_exits_with_failure([&]() { run_latin3_count_row0(bad_range); });
    assert_exits_with_failure([&]() { run_latin3_count_row0(bad_perm); });
}

void test_latin_square_count_print() {
    const int row0[3] = {0, 1, 2};
    run_latin3_count_row0(row0);
    {
        ScopedEnv env("QSIM_LATIN_MAX_PRINT", "1");
        run_latin3_print_all_row0(row0);
    }
    run_latin3_grover_demo_row0(row0, 1);
    run_latin3_grover_demo(0);
}

void test_shor_api_paths() {
    ShorResult bad = run_shor_quantum_part(1, 2);
    if (bad.ok || bad.error.empty()) {
        throw std::runtime_error("Expected shor_api error for N<2");
    }

    Bitstring big = (1ULL << 32);
    ShorResult too_big = run_shor_quantum_part(big, 2);
    if (too_big.ok || too_big.error.empty()) {
        throw std::runtime_error("Expected shor_api error for large N");
    }

    ShorResult ok = run_shor_quantum_part(15, 7);
    if (!ok.ok) {
        throw std::runtime_error("Expected shor_api to succeed for N=15");
    }
}

void test_shor_classical_estimate_order()
{
    Bitstring r = estimate_order(64, 8, 2, 15);
    if (r != 4) {
        throw std::runtime_error("Expected estimate_order to return r=4");
    }

    Bitstring r_none = estimate_order(0, 8, 2, 15);
    if (r_none != 0) {
        throw std::runtime_error("Expected estimate_order to return 0 when no r is valid");
    }
}

void test_shor_quantum_free_function_paths()
{
    ShorQuantumResult too_big = run_shor_algorithm_quantum_part(1ULL << 32, 2);
    if (too_big.ok || too_big.error.empty()) {
        throw std::runtime_error("Expected control-register size error");
    }

    ShorQuantumResult max_qubits = run_shor_algorithm_quantum_part(257, 2);
    if (max_qubits.ok || max_qubits.error.empty()) {
        throw std::runtime_error("Expected max-qubits error");
    }

    ShorQuantumResult ok = run_shor_algorithm_quantum_part(15, 2);
    if (!ok.ok || ok.n_c == 0) {
        throw std::runtime_error("Expected shor_quantum free function to succeed");
    }
}

void test_shor_demo_branches() {
    ScopedEnv env_attempts("QSIM_SHOR_MAX_ATTEMPTS", "1");
    run_shor_demo(1);
    run_shor_demo(10);

    {
        ScopedEnv env_attempts2("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env("QSIM_SHOR_FORCE_A", "5");
        run_shor_demo(15);
    }

    run_shor_demo(1ULL << 20);
    {
        ScopedEnv env_attempts2b("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        run_shor_demo((1ULL << 20) + 1); // odd, triggers max-qubits guard
    }
    {
        ScopedEnv env_attempts2c("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        run_shor_demo((1ULL << 20) + 1); // env_x path to guard
    }

    {
        ScopedEnv env_attempts3("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts3b("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts3c("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts4("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "0");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts5("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "1");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "1");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts5b("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "1");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "2");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts6("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "14");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "2");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts7("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "6");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "6");
        run_shor_demo(9);
    }

    {
        ScopedEnv env_attempts8("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "1");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        ScopedEnv env_r("QSIM_SHOR_FORCE_R", "4");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts9("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "4");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "4");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts10("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "7");
        ScopedEnv env_x("QSIM_SHOR_FORCE_X", "0");
        ScopedEnv env_nc("QSIM_SHOR_FORCE_NC", "0");
        run_shor_demo(15);
    }

    {
        ScopedEnv env_attempts11("QSIM_SHOR_MAX_ATTEMPTS", "1");
        ScopedEnv env_a("QSIM_SHOR_FORCE_A", "2");
        run_shor_demo((1ULL << 32) + 1);
    }
}

void test_algorithm_demo_wrappers() {
    const std::vector<double> q3 = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    run_deutsch_jozsa_demo(3, "CONST0");
    run_deutsch_jozsa_demo(3, "BALANCED_PARITY");
    run_deutsch_jozsa_demo(0, "CONST0");
    run_deutsch_jozsa_demo(3, "NOT_AN_ORACLE");

    run_bernstein_vazirani_demo(3, 0b101, 0);
    run_bernstein_vazirani_demo(0, 0, 0);

    run_grover_search_multi(3, std::vector<Bitstring>());
    run_grover_search_multi(3, std::vector<Bitstring>{3});
    run_grover_search_multi(static_cast<State*>(nullptr), std::vector<Bitstring>{0});
    run_grover_search_multi(static_cast<State*>(nullptr), std::vector<Bitstring>{5});

    run_qubo_exact_cli(3, q3);
    run_qubo_exact_cli(3, std::vector<double>{1.0, 0.0});
    run_qubo_grover_cli(3, -3.0, -1, q3);
    run_qubo_grover_cli(3, -100.0, -1, q3);
    run_qubo_grover_cli(3, -3.0, -1, std::vector<double>{1.0, 0.0});
    run_qubo_demo();

    run_vqa_qaoa_cli(3, 1, 0, 5, 0.2, q3);
    run_vqa_qaoa_cli(3, -1, 0, 5, 0.2, q3);
    run_vqa_demo();

    run_qaoa_qubo_cli(3, 1, 0, 5, 0.2, q3);
    run_qaoa_qubo_cli(3, -1, 0, 5, 0.2, q3);
    run_qaoa_demo();

    run_anneal_qubo_cli("SQA", 3, 8, 4, 0.1, 2.0, 4, q3);
    run_anneal_qubo_cli("BAD", 3, 8, 4, 0.1, 2.0, 4, q3);
    run_anneal_qubo_cli("SQA", 3, 0, 4, 0.1, 2.0, 4, q3);
    run_anneal_demo();

    VqeHamiltonian h;
    h.n_qubits = 1;
    VqePauliTerm term;
    term.coeff = -1.0;
    term.ops.push_back(VqePauliOp{'Z', 0});
    h.terms.push_back(term);
    run_vqe_cli(h, 1, 3, 0.25, 0);
    run_vqe_cli(VqeHamiltonian(), 1, 3, 0.25, 0);
    run_vqe_demo();
    run_quantum_counting_cli(3, std::vector<Bitstring>{1, 6}, 4);
    run_quantum_counting_demo();
    run_simon_cli(4, 0b1010, 8);
    run_simon_demo();

    pid_t pid = fork();
    if (pid < 0) {
        throw std::runtime_error("fork failed in shor seeded branch test");
    }
    if (pid == 0) {
        setenv("QSIM_RNG_SEED", "1234", 1);
        setenv("QSIM_SHOR_MAX_ATTEMPTS", "1", 1);
        setenv("QSIM_SHOR_FORCE_A", "7", 1);
        setenv("QSIM_SHOR_FORCE_X", "1", 1);
        setenv("QSIM_SHOR_FORCE_NC", "4", 1);
        setenv("QSIM_SHOR_FORCE_R", "4", 1);
        run_shor_demo(15);
        std::exit(0);
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        throw std::runtime_error("waitpid failed in shor seeded branch test");
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        throw std::runtime_error("child failed in shor seeded branch test");
    }
}
