#include "state.hh"
#include "tests/helpers.hh"
#include "algorithms/qrng.hh"
#include "algorithms/latin_square.hh"
#include "algorithms/deutsch_jozsa.hh"
#include "algorithms/bernstein_vazirani.hh"
#include "algorithms/qubo.hh"
#include "algorithms/vqa_qaoa.hh"
#include "algorithms/qaoa.hh"
#include "algorithms/vqe.hh"
#include "algorithms/anneal.hh"
#include "algorithms/api/grover_api.hh"
#include "algorithms/api/shor_api.hh"
#include "algorithms/shor_classical.hh"
#include "algorithms/shor_quantum.hh"
#include "cli/commands.hh"
#include "logging.hh"
#include "demos/grover_demo.hh"
#include "demos/latin_demo.hh"
#include "demos/qaoa_demo.hh"
#include "demos/shor_demo.hh"
#include "modular_exp.hh"
#include "math/mod_arith.hh"
#include "math/bit_ops.hh"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

using test_helpers::assert_complex_equal;
using test_helpers::assert_complex_close;
using test_helpers::assert_equal;
using test_helpers::assert_amplitude_match;
using test_helpers::assert_amplitude_magnitude;
using test_helpers::ScopedEnv;

void test_mod_arith_gcd_basic() {
    assert_equal(gcd_bitstring(0, 7), 7ULL, "gcd(0,7)=7");
    assert_equal(gcd_bitstring(54, 24), 6ULL, "gcd(54,24)=6");
}

void test_mod_arith_mod_pow_edges() {
    assert_equal(mod_pow(5, 3, 1), 0ULL, "mod_pow with mod=1 returns 0");
    assert_equal(mod_pow(2, 0, 15), 1ULL, "mod_pow exponent 0 returns 1");
    assert_equal(mod_pow(7, 4, 15), 1ULL, "mod_pow periodicity check");
}

void test_qrng_edges() {
    std::vector<int> empty = qrng_bits(0, nullptr);
    if (!empty.empty()) {
        throw std::runtime_error("qrng_bits(0) should return empty");
    }

    std::vector<int> bits_auto = qrng_bits(2, nullptr);
    if (bits_auto.size() != 2) {
        throw std::runtime_error("qrng_bits should return 2 bits");
    }

    std::vector<double> rv = {0.0, 1.0};
    std::vector<int> bits = qrng_bits(2, &rv);
    if (bits.size() != 2 || bits[0] != 0 || bits[1] != 1) {
        throw std::runtime_error("qrng_bits with deterministic RNG failed");
    }

    ScopedEnv cap("QSIM_QRNG_MAX_QUBITS", "8");
    std::vector<double> zeros(8, 0.0);
    uint64_t val = qrng_u64(70, &zeros);
    if (val != 0ULL) {
        throw std::runtime_error("qrng_u64 expected all-zero output with zero RNG");
    }

    if (qrng_u64(0, nullptr) != 0ULL) {
        throw std::runtime_error("qrng_u64(0) should return 0");
    }
}

void test_deutsch_jozsa_constant_oracles()
{
    DeutschJozsaResult r0 = run_deutsch_jozsa(3, DeutschJozsaOracle::ConstantZero);
    if (!r0.ok || !r0.is_constant) {
        throw std::runtime_error("Expected ConstantZero to be classified as constant");
    }
    for (int bit : r0.input_measurement) {
        if (bit != 0) {
            throw std::runtime_error("Expected ConstantZero to measure all zeros");
        }
    }

    DeutschJozsaResult r1 = run_deutsch_jozsa(3, DeutschJozsaOracle::ConstantOne);
    if (!r1.ok || !r1.is_constant) {
        throw std::runtime_error("Expected ConstantOne to be classified as constant");
    }
    for (int bit : r1.input_measurement) {
        if (bit != 0) {
            throw std::runtime_error("Expected ConstantOne to measure all zeros");
        }
    }
}

void test_deutsch_jozsa_balanced_oracles()
{
    DeutschJozsaResult r0 = run_deutsch_jozsa(3, DeutschJozsaOracle::BalancedXor0);
    if (!r0.ok || r0.is_constant) {
        throw std::runtime_error("Expected BalancedXor0 to be classified as balanced");
    }
    bool any_one = false;
    for (int bit : r0.input_measurement) {
        if (bit != 0) {
            any_one = true;
            break;
        }
    }
    if (!any_one) {
        throw std::runtime_error("Expected BalancedXor0 to measure a non-zero input register");
    }

    DeutschJozsaResult r1 = run_deutsch_jozsa(3, DeutschJozsaOracle::BalancedParity);
    if (!r1.ok || r1.is_constant) {
        throw std::runtime_error("Expected BalancedParity to be classified as balanced");
    }
}

void test_deutsch_jozsa_oracle_helpers()
{
    DeutschJozsaOracle oracle;
    if (!parse_deutsch_jozsa_oracle("CONST0", oracle)) {
        throw std::runtime_error("Expected CONST0 to parse");
    }
    if (std::string(deutsch_jozsa_oracle_name(oracle)).empty()) {
        throw std::runtime_error("Expected oracle name to be non-empty");
    }
    if (!parse_deutsch_jozsa_oracle("CONST1", oracle)) {
        throw std::runtime_error("Expected CONST1 to parse");
    }
    if (std::string(deutsch_jozsa_oracle_name(oracle)).empty()) {
        throw std::runtime_error("Expected CONST1 oracle name to be non-empty");
    }
    if (!parse_deutsch_jozsa_oracle("balanced_xor0", oracle)) {
        throw std::runtime_error("Expected BALANCED_XOR0 to parse case-insensitive");
    }
    if (std::string(deutsch_jozsa_oracle_name(oracle)).empty()) {
        throw std::runtime_error("Expected BALANCED_XOR0 oracle name to be non-empty");
    }
    if (!parse_deutsch_jozsa_oracle("balanced_parity", oracle)) {
        throw std::runtime_error("Expected BALANCED_PARITY to parse case-insensitive");
    }
    if (std::string(deutsch_jozsa_oracle_name(oracle)).empty()) {
        throw std::runtime_error("Expected BALANCED_PARITY oracle name to be non-empty");
    }
    if (parse_deutsch_jozsa_oracle("unknown", oracle)) {
        throw std::runtime_error("Expected unknown oracle to fail parsing");
    }

    DeutschJozsaResult bad = run_deutsch_jozsa(0, DeutschJozsaOracle::ConstantZero);
    if (bad.ok || bad.error.empty()) {
        throw std::runtime_error("Expected invalid n_inputs to return error");
    }

    DeutschJozsaResult too_large = run_deutsch_jozsa(63, DeutschJozsaOracle::ConstantZero);
    if (too_large.ok || too_large.error.empty()) {
        throw std::runtime_error("Expected too-large n_inputs to return error");
    }

    const char* unknown = deutsch_jozsa_oracle_name(static_cast<DeutschJozsaOracle>(-1));
    if (std::string(unknown).empty()) {
        throw std::runtime_error("Expected unknown oracle name fallback");
    }
}

void test_bernstein_vazirani_secret_recovery()
{
    BernsteinVaziraniResult r0 = run_bernstein_vazirani(5, 0b10110, 0);
    if (!r0.ok) {
        throw std::runtime_error("Expected BV run to succeed");
    }
    if (r0.measured_secret != 0b10110ULL) {
        throw std::runtime_error("Expected BV to recover hidden string");
    }

    BernsteinVaziraniResult r1 = run_bernstein_vazirani(5, 0b10110, 1);
    if (!r1.ok) {
        throw std::runtime_error("Expected BV run with bias=1 to succeed");
    }
    if (r1.measured_secret != 0b10110ULL) {
        throw std::runtime_error("Expected BV to recover hidden string regardless of bias");
    }
}

void test_bernstein_vazirani_errors()
{
    BernsteinVaziraniResult n_bad = run_bernstein_vazirani(0, 0, 0);
    if (n_bad.ok || n_bad.error.empty()) {
        throw std::runtime_error("Expected BV error for n_inputs <= 0");
    }

    BernsteinVaziraniResult secret_bad = run_bernstein_vazirani(3, 8, 0);
    if (secret_bad.ok || secret_bad.error.empty()) {
        throw std::runtime_error("Expected BV error for out-of-range secret");
    }

    BernsteinVaziraniResult bias_bad = run_bernstein_vazirani(3, 1, 2);
    if (bias_bad.ok || bias_bad.error.empty()) {
        throw std::runtime_error("Expected BV error for invalid bias");
    }

    BernsteinVaziraniResult too_large = run_bernstein_vazirani(70, 1, 0);
    if (too_large.ok || too_large.error.empty()) {
        throw std::runtime_error("Expected BV error for too-large n_inputs");
    }
}

void test_qubo_exact_solver()
{
    const int n = 3;
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    QuboExactResult result = qubo_solve_exact(n, q);
    if (!result.ok) {
        throw std::runtime_error("Expected QUBO exact solver to succeed");
    }
    if (result.argmin != 4ULL) {
        throw std::runtime_error("Expected QUBO exact argmin to be 4 (0b100)");
    }
    if (std::abs(result.min_value - (-3.0)) > 1e-12) {
        throw std::runtime_error("Expected QUBO exact minimum value to be -3");
    }
    if (std::abs(qubo_evaluate(4ULL, n, q) - (-3.0)) > 1e-12) {
        throw std::runtime_error("Expected qubo_evaluate(4) to be -3");
    }

    const std::vector<double> zero_q = {
        0.0, 0.0,
        0.0, 0.0
    };
    QuboExactResult tie = qubo_solve_exact(2, zero_q);
    if (!tie.ok) {
        throw std::runtime_error("Expected zero QUBO exact solve to succeed");
    }
    if (tie.minimizers.size() != 4) {
        throw std::runtime_error("Expected all four assignments to tie at value 0");
    }
}

void test_qubo_grover_threshold_solver()
{
    const int n = 3;
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    QuboGroverResult none = qubo_solve_grover_threshold(n, q, -100.0, -1);
    if (!none.ok || none.found) {
        throw std::runtime_error("Expected no candidate for very low threshold");
    }

    QuboGroverResult result = qubo_solve_grover_threshold(n, q, -3.0, -1);
    if (!result.ok || !result.found) {
        throw std::runtime_error("Expected QUBO Grover threshold to find a candidate");
    }
    if (result.candidate != 4ULL) {
        throw std::runtime_error("Expected QUBO Grover candidate to be 4 (0b100)");
    }
    if (std::abs(result.candidate_value - (-3.0)) > 1e-12) {
        throw std::runtime_error("Expected QUBO Grover candidate value to be -3");
    }
}

void test_qubo_error_paths()
{
    std::string error;
    if (qubo_matrix_valid(0, std::vector<double>(), error)) {
        throw std::runtime_error("Expected qubo_matrix_valid to reject n=0");
    }

    std::vector<double> bad_matrix = {1.0, 2.0, 3.0};
    if (qubo_matrix_valid(2, bad_matrix, error)) {
        throw std::runtime_error("Expected qubo_matrix_valid to reject invalid matrix size");
    }

    QuboExactResult exact_bad = qubo_solve_exact(2, bad_matrix);
    if (exact_bad.ok || exact_bad.error.empty()) {
        throw std::runtime_error("Expected qubo_solve_exact to fail on invalid matrix");
    }

    QuboGroverResult grover_bad = qubo_solve_grover_threshold(70, std::vector<double>(), 0.0, -1);
    if (grover_bad.ok || grover_bad.error.empty()) {
        throw std::runtime_error("Expected qubo_solve_grover_threshold to fail for invalid n");
    }
}

void test_vqa_qaoa_finds_good_candidate()
{
    const int n = 3;
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    VqaQaoaOptions opts;
    opts.p_layers = 1;
    opts.max_iters = 50;
    opts.shots = 0;
    opts.step_size = 0.3;
    opts.seed = 17u;

    VqaQaoaResult result = run_vqa_qaoa_qubo(n, q, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected VQA QAOA solve to succeed");
    }
    if (result.energy_history.empty()) {
        throw std::runtime_error("Expected VQA QAOA to record optimization history");
    }
    if (result.evaluations <= 0) {
        throw std::runtime_error("Expected VQA QAOA evaluations to be positive");
    }

    QuboExactResult exact = qubo_solve_exact(n, q);
    if (!exact.ok) {
        throw std::runtime_error("Expected exact QUBO solve to succeed");
    }
    if (result.best_energy > exact.min_value + 1.00) {
        throw std::runtime_error("Expected VQA QAOA best energy to be close to exact minimum");
    }
}

void test_vqa_qaoa_improves_over_initial_state()
{
    const int n = 2;
    const std::vector<double> q = {
        -1.0,  0.0,
         0.0, -2.0
    };

    VqaQaoaOptions opts;
    opts.p_layers = 1;
    opts.max_iters = 30;
    opts.shots = 0;
    opts.step_size = 0.25;
    opts.seed = 11u;

    VqaQaoaResult result = run_vqa_qaoa_qubo(n, q, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected VQA QAOA solve to succeed");
    }
    if (result.energy_history.size() < 2) {
        throw std::runtime_error("Expected at least two energy history points");
    }
    if (result.best_energy > result.energy_history.front() + 1e-12) {
        throw std::runtime_error("Expected VQA QAOA to not regress from initial energy");
    }
}

void test_vqa_qaoa_error_paths()
{
    const std::vector<double> q2 = {
        1.0, 0.0,
        0.0, 1.0
    };

    VqaQaoaOptions bad_p;
    bad_p.p_layers = -1;
    VqaQaoaResult p_result = run_vqa_qaoa_qubo(2, q2, bad_p);
    if (p_result.ok || p_result.error.empty()) {
        throw std::runtime_error("Expected VQA QAOA to reject negative p_layers");
    }

    VqaQaoaOptions bad_step;
    bad_step.step_size = 0.0;
    VqaQaoaResult step_result = run_vqa_qaoa_qubo(2, q2, bad_step);
    if (step_result.ok || step_result.error.empty()) {
        throw std::runtime_error("Expected VQA QAOA to reject non-positive step_size");
    }

    VqaQaoaOptions bad_shots;
    bad_shots.shots = -4;
    VqaQaoaResult shot_result = run_vqa_qaoa_qubo(2, q2, bad_shots);
    if (shot_result.ok || shot_result.error.empty()) {
        throw std::runtime_error("Expected VQA QAOA to reject negative shots");
    }

    VqaQaoaOptions bad_iters;
    bad_iters.max_iters = -2;
    VqaQaoaResult iter_result = run_vqa_qaoa_qubo(2, q2, bad_iters);
    if (iter_result.ok || iter_result.error.empty()) {
        throw std::runtime_error("Expected VQA QAOA to reject negative max_iters");
    }

    std::vector<double> bad_matrix = {1.0, 0.0, 0.0};
    VqaQaoaResult matrix_result = run_vqa_qaoa_qubo(2, bad_matrix, VqaQaoaOptions());
    if (matrix_result.ok || matrix_result.error.empty()) {
        throw std::runtime_error("Expected VQA QAOA to reject invalid matrix size");
    }
}

void test_vqa_qaoa_shot_mode_path()
{
    const int n = 2;
    const std::vector<double> q = {
        -1.0, 0.0,
         0.0, -2.0
    };

    VqaQaoaOptions opts;
    opts.p_layers = 1;
    opts.max_iters = 3;
    opts.shots = 64;
    opts.step_size = 0.2;
    opts.seed = 123u;

    VqaQaoaResult result = run_vqa_qaoa_qubo(n, q, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected VQA QAOA shot-mode solve to succeed");
    }
    if (result.energy_history.size() != 4) {
        throw std::runtime_error("Expected shot-mode history to include init + max_iters entries");
    }
}

void test_qaoa_wrapper_paths()
{
    const int n = 3;
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    QaoaOptions opts;
    opts.p_layers = 1;
    opts.max_iters = 20;
    opts.shots = 0;
    opts.step_size = 0.25;

    QaoaResult result = run_qaoa_qubo(n, q, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected QAOA wrapper run to succeed");
    }
    if (result.energy_history.empty()) {
        throw std::runtime_error("Expected QAOA wrapper to produce history");
    }

    QaoaOptions bad_step = opts;
    bad_step.step_size = 0.0;
    QaoaResult bad = run_qaoa_qubo(n, q, bad_step);
    if (bad.ok || bad.error.empty()) {
        throw std::runtime_error("Expected QAOA wrapper error propagation for bad options");
    }
}

void test_qaoa_demo_paths()
{
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };
    run_qaoa_qubo_cli(3, 1, 0, 3, 0.2, q);
    run_qaoa_demo();
}

void test_vqe_single_qubit_ground_state()
{
    VqeHamiltonian h;
    h.n_qubits = 1;
    VqePauliTerm z_term;
    z_term.coeff = 1.0;
    z_term.ops.push_back(VqePauliOp{'Z', 0});
    h.terms.push_back(z_term);

    VqeOptions opts;
    opts.layers = 1;
    opts.max_iters = 30;
    opts.step_size = 0.3;
    opts.shots = 0;

    VqeResult result = run_vqe(h, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected VQE run to succeed");
    }
    if (result.best_energy > -0.95) {
        throw std::runtime_error("Expected VQE to approach Z ground-state energy -1");
    }
}

void test_vqe_expectation_pauli_terms()
{
    VqeHamiltonian hx;
    hx.n_qubits = 1;
    VqePauliTerm x_term;
    x_term.coeff = 1.0;
    x_term.ops.push_back(VqePauliOp{'X', 0});
    hx.terms.push_back(x_term);

    State plus(1, 0);
    plus.h(0);
    if (std::abs(vqe_expectation_exact(plus, hx) - 1.0) > 1e-9) {
        throw std::runtime_error("Expected <+|X|+> = 1");
    }

    VqeHamiltonian hy;
    hy.n_qubits = 1;
    VqePauliTerm y_term;
    y_term.coeff = 1.0;
    y_term.ops.push_back(VqePauliOp{'Y', 0});
    hy.terms.push_back(y_term);

    State zero(1, 0);
    if (std::abs(vqe_expectation_exact(zero, hy)) > 1e-9) {
        throw std::runtime_error("Expected <0|Y|0> = 0");
    }

    VqeHamiltonian hz;
    hz.n_qubits = 1;
    VqePauliTerm z_term;
    z_term.coeff = 1.0;
    z_term.ops.push_back(VqePauliOp{'Z', 0});
    hz.terms.push_back(z_term);

    State one(1, 0);
    one.x(0);
    if (std::abs(vqe_expectation_exact(one, hz) - (-1.0)) > 1e-9) {
        throw std::runtime_error("Expected <1|Z|1> = -1");
    }
}

void test_vqe_error_and_edge_paths()
{
    VqeHamiltonian bad_h;
    bad_h.n_qubits = 0;
    std::string error;
    if (vqe_hamiltonian_valid(bad_h, error)) {
        throw std::runtime_error("Expected VQE Hamiltonian validation to reject n_qubits=0");
    }

    VqeHamiltonian empty_terms;
    empty_terms.n_qubits = 1;
    if (vqe_hamiltonian_valid(empty_terms, error)) {
        throw std::runtime_error("Expected VQE Hamiltonian validation to reject empty term list");
    }

    VqeResult invalid_h_run = run_vqe(empty_terms, VqeOptions());
    if (invalid_h_run.ok || invalid_h_run.error.empty()) {
        throw std::runtime_error("Expected run_vqe to fail on invalid Hamiltonian");
    }

    VqeHamiltonian bad_op;
    bad_op.n_qubits = 1;
    VqePauliTerm bad_term;
    bad_term.coeff = 1.0;
    bad_term.ops.push_back(VqePauliOp{'Q', 0});
    bad_op.terms.push_back(bad_term);
    if (vqe_hamiltonian_valid(bad_op, error)) {
        throw std::runtime_error("Expected VQE Hamiltonian validation to reject invalid Pauli op");
    }

    VqeHamiltonian out_of_range;
    out_of_range.n_qubits = 1;
    VqePauliTerm oob_term;
    oob_term.coeff = 1.0;
    oob_term.ops.push_back(VqePauliOp{'X', 3});
    out_of_range.terms.push_back(oob_term);
    if (vqe_hamiltonian_valid(out_of_range, error)) {
        throw std::runtime_error("Expected VQE Hamiltonian validation to reject out-of-range qubit");
    }

    VqeHamiltonian simple;
    simple.n_qubits = 1;
    VqePauliTerm z_term;
    z_term.coeff = 1.0;
    z_term.ops.push_back(VqePauliOp{'Z', 0});
    simple.terms.push_back(z_term);

    VqeOptions bad_layers;
    bad_layers.layers = -1;
    VqeResult r_layers = run_vqe(simple, bad_layers);
    if (r_layers.ok || r_layers.error.empty()) {
        throw std::runtime_error("Expected VQE to reject negative layers");
    }

    VqeOptions bad_iters;
    bad_iters.max_iters = -1;
    VqeResult r_iters = run_vqe(simple, bad_iters);
    if (r_iters.ok || r_iters.error.empty()) {
        throw std::runtime_error("Expected VQE to reject negative max_iters");
    }

    VqeOptions bad_step;
    bad_step.step_size = 0.0;
    VqeResult r_step = run_vqe(simple, bad_step);
    if (r_step.ok || r_step.error.empty()) {
        throw std::runtime_error("Expected VQE to reject non-positive step_size");
    }

    VqeOptions bad_shots;
    bad_shots.shots = 16;
    VqeResult r_shots = run_vqe(simple, bad_shots);
    if (r_shots.ok || r_shots.error.empty()) {
        throw std::runtime_error("Expected VQE to reject shots!=0 for now");
    }

    VqeOptions zero_layer_opts;
    zero_layer_opts.layers = 0;
    zero_layer_opts.max_iters = 2;
    zero_layer_opts.step_size = 0.2;
    VqeResult zero_layer = run_vqe(simple, zero_layer_opts);
    if (!zero_layer.ok) {
        throw std::runtime_error("Expected VQE with zero layers to succeed");
    }
    if (zero_layer.energy_history.size() != 3) {
        throw std::runtime_error("Expected VQE zero-layer history to include init + max_iters");
    }

    // Exercise 2-qubit ansatz entangler path.
    VqeHamiltonian two_qubit;
    two_qubit.n_qubits = 2;
    VqePauliTerm t2;
    t2.coeff = 1.0;
    t2.ops.push_back(VqePauliOp('Z', 0));
    two_qubit.terms.push_back(t2);
    VqeOptions two_opts;
    two_opts.layers = 1;
    two_opts.max_iters = 1;
    two_opts.step_size = 0.2;
    VqeResult two_result = run_vqe(two_qubit, two_opts);
    if (!two_result.ok) {
        throw std::runtime_error("Expected 2-qubit VQE run to succeed");
    }
}

void test_anneal_sa_finds_good_candidate()
{
    const int n = 3;
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    AnnealOptions opts;
    opts.method = AnnealMethod::SA;
    opts.steps = 60;
    opts.sweeps_per_step = 20;
    opts.beta_start = 0.1;
    opts.beta_end = 5.0;
    opts.seed = 19u;

    AnnealResult result = anneal_qubo(n, q, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected SA anneal to succeed");
    }
    if (result.best_history.size() != static_cast<size_t>(opts.steps + 1)) {
        throw std::runtime_error("Expected SA best history to match steps+1");
    }
    if (result.attempted_moves <= 0 || result.accepted_moves <= 0) {
        throw std::runtime_error("Expected SA to attempt and accept moves");
    }
    if (result.best_value > -2.0) {
        throw std::runtime_error("Expected SA to find a low-energy assignment");
    }
}

void test_anneal_sqa_finds_good_candidate()
{
    const int n = 3;
    const std::vector<double> q = {
        -2.0,  0.0,  2.0,
         0.0,  1.0,  0.0,
         2.0,  0.0, -3.0
    };

    AnnealOptions opts;
    opts.method = AnnealMethod::SQA;
    opts.steps = 60;
    opts.sweeps_per_step = 20;
    opts.beta_start = 0.1;
    opts.beta_end = 5.0;
    opts.replicas = 8;
    opts.seed = 19u;

    AnnealResult result = anneal_qubo(n, q, opts);
    if (!result.ok) {
        throw std::runtime_error("Expected SQA anneal to succeed");
    }
    if (result.best_history.size() != static_cast<size_t>(opts.steps + 1)) {
        throw std::runtime_error("Expected SQA best history to match steps+1");
    }
    if (result.attempted_moves <= 0 || result.accepted_moves <= 0) {
        throw std::runtime_error("Expected SQA to attempt and accept moves");
    }
    if (result.best_value > -2.5) {
        throw std::runtime_error("Expected SQA to find a low-energy assignment");
    }
}

void test_anneal_options_and_parsing()
{
    const std::vector<double> q2 = {
        1.0, 0.0,
        0.0, 1.0
    };
    std::string error;

    AnnealOptions bad_steps;
    bad_steps.steps = 0;
    if (anneal_options_valid(2, q2, bad_steps, error)) {
        throw std::runtime_error("Expected anneal_options_valid to reject steps <= 0");
    }

    AnnealOptions bad_beta;
    bad_beta.beta_end = 0.0;
    if (anneal_options_valid(2, q2, bad_beta, error)) {
        throw std::runtime_error("Expected anneal_options_valid to reject non-positive beta");
    }

    AnnealOptions bad_schedule;
    bad_schedule.beta_start = 2.0;
    bad_schedule.beta_end = 1.0;
    if (anneal_options_valid(2, q2, bad_schedule, error)) {
        throw std::runtime_error("Expected anneal_options_valid to reject beta_end < beta_start");
    }

    AnnealOptions bad_sqa;
    bad_sqa.method = AnnealMethod::SQA;
    bad_sqa.replicas = 1;
    if (anneal_options_valid(2, q2, bad_sqa, error)) {
        throw std::runtime_error("Expected anneal_options_valid to require replicas >= 2 for SQA");
    }

    AnnealOptions bad_sweeps;
    bad_sweeps.sweeps_per_step = 0;
    if (anneal_options_valid(2, q2, bad_sweeps, error)) {
        throw std::runtime_error("Expected anneal_options_valid to reject sweeps_per_step <= 0");
    }

    AnnealOptions bad_replicas_sa;
    bad_replicas_sa.method = AnnealMethod::SA;
    bad_replicas_sa.replicas = 0;
    if (anneal_options_valid(2, q2, bad_replicas_sa, error)) {
        throw std::runtime_error("Expected anneal_options_valid to reject non-positive replicas in SA mode");
    }

    std::vector<double> bad_q = {1.0, 0.0, 0.0};
    if (anneal_options_valid(2, bad_q, AnnealOptions(), error)) {
        throw std::runtime_error("Expected anneal_options_valid to reject invalid matrix shape");
    }

    AnnealResult invalid_result = anneal_qubo(2, bad_q, AnnealOptions());
    if (invalid_result.ok || invalid_result.error.empty()) {
        throw std::runtime_error("Expected anneal_qubo to fail when options/matrix are invalid");
    }

    AnnealMethod method = AnnealMethod::SA;
    if (!parse_anneal_method("SA", method) || method != AnnealMethod::SA) {
        throw std::runtime_error("Expected parse_anneal_method to parse SA");
    }
    if (!parse_anneal_method("quantum", method) || method != AnnealMethod::SQA) {
        throw std::runtime_error("Expected parse_anneal_method to parse SQA aliases");
    }
    if (parse_anneal_method("???", method)) {
        throw std::runtime_error("Expected parse_anneal_method to reject unknown method");
    }
    if (std::string(anneal_method_name(AnnealMethod::SA)) != "SA") {
        throw std::runtime_error("Expected anneal_method_name(SA) == SA");
    }
    if (std::string(anneal_method_name(static_cast<AnnealMethod>(99))) != "UNKNOWN") {
        throw std::runtime_error("Expected anneal_method_name unknown enum to return UNKNOWN");
    }
}

void test_anneal_sqa_improvement_branch()
{
    const int n = 1;
    const std::vector<double> q = {-1.0};

    bool observed_improvement = false;
    for (unsigned seed = 1; seed <= 24; ++seed) {
        AnnealOptions opts;
        opts.method = AnnealMethod::SQA;
        opts.steps = 4;
        opts.sweeps_per_step = 4;
        opts.beta_start = 0.1;
        opts.beta_end = 2.0;
        opts.replicas = 4;
        opts.seed = seed;

        AnnealResult result = anneal_qubo(n, q, opts);
        if (!result.ok) {
            throw std::runtime_error("Expected SQA anneal run to succeed in improvement branch test");
        }
        if (!result.best_history.empty() && result.best_value < result.best_history.front() - 1e-12) {
            observed_improvement = true;
            break;
        }
    }
    if (!observed_improvement) {
        throw std::runtime_error("Expected at least one SQA run to exercise best-value improvement path");
    }
}

void test_logging_levels()
{
    auto run_child = [](const char* key, const char* value, qsim_log::Level expected) {
        pid_t pid = fork();
        if (pid < 0) {
            throw std::runtime_error("fork failed in logging test");
        }
        if (pid == 0) {
            if (key) {
                setenv(key, value, 1);
            }
            qsim_log::Level got = qsim_log::get_level();
            std::exit(got == expected ? 0 : 1);
        }
        int status = 0;
        if (waitpid(pid, &status, 0) == -1) {
            throw std::runtime_error("waitpid failed in logging test");
        }
        if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
            throw std::runtime_error("logging env test failed");
        }
    };

    run_child("QSIM_VERBOSE", "VERBOSE", qsim_log::Level::Verbose);
    run_child("QSIM_VERBOSE", "INVALID", qsim_log::Level::Normal);
    run_child("QSIM_LOG_LEVEL", "QUIET", qsim_log::Level::Quiet);
    run_child("QSIM_GROVER_VERBOSE", "1", qsim_log::Level::Verbose);

    qsim_log::Level parsed = qsim_log::Level::Quiet;
    if (!qsim_log::parse_level("NORMAL", parsed) || parsed != qsim_log::Level::Normal) {
        throw std::runtime_error("parse_level failed for NORMAL");
    }
    if (qsim_log::parse_level("???", parsed)) {
        throw std::runtime_error("parse_level should fail for unknown token");
    }

    std::ostringstream out;
    qsim_log::set_stream(&out);
    qsim_log::set_level(qsim_log::Level::Quiet);
    qsim_log::log(qsim_log::Level::Normal, "hidden\n");
    if (!out.str().empty()) {
        throw std::runtime_error("log should be quiet at QUIET level");
    }
    qsim_log::set_level(qsim_log::Level::Normal);
    qsim_log::log(qsim_log::Level::Normal, "visible\n");
    if (out.str().find("visible") == std::string::npos) {
        throw std::runtime_error("log should emit at NORMAL level");
    }
    qsim_log::set_stream(&std::cout);
}

void test_grover_search_helpers() {
    run_grover_search(nullptr, 1);
    State s(2, 0);
    run_grover_search(&s, 3);
    run_grover_search_multi(&s, std::vector<Bitstring>());
    run_grover_search_multi(&s, std::vector<Bitstring>{1, 2});
}

void test_grover_api_errors() {
    GroverResult bad_n = run_grover(0, std::vector<Bitstring>{1});
    if (bad_n.ok || bad_n.error.empty()) {
        throw std::runtime_error("Expected error for invalid auto-created qubit count");
    }

    State s(2, 0);
    GroverResult empty = run_grover(s, std::vector<Bitstring>());
    if (empty.ok || empty.error.empty()) {
        throw std::runtime_error("Expected error on empty target list");
    }

    State big(63, 0);
    GroverResult too_big = run_grover(big, std::vector<Bitstring>{1});
    if (too_big.ok || too_big.error.empty()) {
        throw std::runtime_error("Expected error for too many qubits");
    }

    GroverResult out_of_range = run_grover(s, std::vector<Bitstring>{7});
    if (out_of_range.ok || out_of_range.error.empty()) {
        throw std::runtime_error("Expected out-of-range target error");
    }

    State s2(3, 0);
    GroverResult r0 = run_grover(s2, std::vector<Bitstring>{1}, 0);
    if (!r0.ok || r0.iterations != 0) {
        throw std::runtime_error("Expected R=0 Grover run to succeed");
    }

    {
        ScopedEnv env("QSIM_GROVER_FORCE_SET", "1");
        GroverResult forced = run_grover(s2, std::vector<Bitstring>{1}, 1);
        if (!forced.ok) {
            throw std::runtime_error("Expected forced set oracle Grover run to succeed");
        }
    }
}
