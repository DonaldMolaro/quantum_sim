# Quantum Computing for the Classical Engineer
## A Self-Study Course Using the `quantum_sim` Simulator

---

**Who this is for:** Software or hardware engineers with no quantum background.
You should be comfortable with boolean logic, basic linear algebra (vectors,
matrices), and the idea of a Fourier transform. You do not need physics.

**What you will use:** The `quantum_sim` CLI (command-line simulator) built in
this repository. Every experiment in this course runs directly inside it.

**Approach:** Each module explains *why* something works using classical
intuition first, then shows the quantum version, then lets you experiment.
The emphasis is on *understanding* rather than memorizing formulas.

**How long:** Plan for roughly 10–15 hours of focused study. Each module is
designed to stand alone, so you can pause and return.

---

## How to Use This Course

Start the simulator from the project root:

```
./quantum_sim
```

To see algorithm internals as they run:

```
VERBOSE VERBOSE
```

To enable step-by-step teaching narration:

```
TUTOR ON
```

For reproducible results (same random outcomes each run), set a seed at the
start of any session:

```
SEED 42
```

Quick sanity checks at any time:

```
CHECK NORMALIZED
CHECK TARGETS 3 5
```

---

## One Critical Fact About Bit Ordering

Before anything else: in this simulator, **`q0` is the least-significant bit
(LSB)**. This matches the convention used in many textbooks on quantum circuits.

When you `DISPLAY` a state like `|101>`, the bits read right-to-left by qubit
index: `q2=1`, `q1=0`, `q0=1`. This is the opposite of how you probably write
binary numbers. Keep this in mind whenever you inspect a bitstring.

---

---

# MODULE 1 — The Qubit

## 1.1 Classical Intuition: a Weighted Coin

A classical bit is a light switch: it is either 0 or 1, and that's the whole
story. A qubit before measurement is more like a coin *while it is still
spinning in the air*: it has a well-defined mathematical description of its
two-outcome tendency, but it has not yet "decided" which face will land.

More precisely, a qubit state is described by two complex numbers — called
**amplitudes** — one for `|0>` and one for `|1>`:

```
|ψ> = α|0> + β|1>
```

The **probability** of measuring 0 is `|α|²`, and the probability of measuring
1 is `|β|²`. These must sum to 1 (the coin must land on *some* face).

The fact that α and β are *complex* — not just real-valued probabilities — is
what makes quantum computing qualitatively different. The phase (angle) of
these complex numbers has no classical analogue, and it enables interference
(Module 2).

## 1.2 The Most Important Single-Qubit Gate: H

The **Hadamard gate** (H) puts a qubit into equal superposition. Starting from
`|0>`:

```
H|0> = (|0> + |1>) / √2
```

Both outcomes are equally likely — 50/50. But the amplitudes are *not* a
probability distribution; they are complex weights that can cancel or reinforce
each other when operations are combined. This is the key insight you will
return to in every module.

Applied to `|1>`:

```
H|1> = (|0> - |1>) / √2
```

Same probabilities, but the amplitude for `|1>` is negative. This difference
(the **phase**) is invisible to a single measurement but is critical for
interference.

## 1.3 Simulator Experiment: Single-Qubit Superposition

```
INIT 1
H 0
DISPLAY
```

You should see two entries with equal amplitude (≈0.7071) and equal probability
(0.5). Now measure and watch the state collapse:

```
MEASURE 0 0
DISPLAY
```

The state is now either `|0>` or `|1>` — the superposition is gone. Run this
several times (using `INIT 1` then `H 0` each time) to see that the outcome is
random but averaged 50/50.

Compare with what H does to `|1>`:

```
INIT 1
X 0
H 0
DISPLAY
```

Same probabilities, but you should see a negative amplitude on `|1>`. This
negative sign is invisible to measurement here — but *not* when we compose
gates.

## 1.4 Key Insight

> A qubit in superposition is not "secretly 0 or 1, we just don't know which."
> It genuinely represents *both possibilities simultaneously*, with complex
> amplitudes that can interfere. The complex phase of those amplitudes is the
> resource that quantum algorithms exploit.

## 1.5 Self-Check

1. After `H 0` on `|0>`, what is the probability of measuring 1?
2. After `X 0` then `H 0`, is the probability of measuring 1 different from
   after just `H 0`?
3. If α = 0.6 and β = 0.8 for a qubit, do the probabilities sum to 1?
   (Hint: compute `|0.6|² + |0.8|²`.)

---

---

# MODULE 2 — Phase: The Hidden Degree of Freedom

## 2.1 Classical Intuition: Signal Phase

In signal processing, two sine waves at the same frequency can be in phase (and
add) or out of phase (and cancel). The *power* of each individual wave is the
same either way — but their combined effect is completely different. Quantum
phase works the same way.

Phase is the *angle* component of a complex amplitude. When you look at just
one amplitude in isolation, |e^{iφ}| = 1 regardless of φ. Phase is therefore
invisible in a single measurement. But when two paths combine — as they do in
any multi-gate circuit — phase determines whether they interfere constructively
(amplitudes add) or destructively (amplitudes cancel).

## 2.2 Phase Gates

The **Z gate** flips the phase of the `|1>` component without changing probabilities:

```
Z(α|0> + β|1>) = α|0> - β|1>
```

The **S gate** applies a 90° (i) phase to `|1>`.

The **T gate** applies a 45° (π/4) phase to `|1>`.

None of these change measurement probabilities on their own. They only matter
in combination with H gates, which convert phase differences into amplitude
differences (probability differences).

## 2.3 Simulator Experiment: Phase Becomes Visible

Start with `|+> = H|0> = (|0> + |1>)/√2`. Apply Z to add a phase, then apply
H again to convert the phase difference into a probability difference:

```
INIT 1
H 0
DISPLAY
```

Both |0> and |1> have probability 0.5. Now add a phase and decode:

```
Z 0
H 0
DISPLAY
```

You should now see the state is `|1>` with probability 1.0. The Z gate flipped
the sign of `|1>`, and the second H made that sign difference *visible* as a
certain measurement outcome.

Contrast with *not* applying Z:

```
INIT 1
H 0
H 0
DISPLAY
```

H applied twice is the identity — you return to `|0>` with probability 1.0.

This H → (optional phase gate) → H pattern is the basic mechanism behind
**phase kickback** and nearly every quantum algorithm.

## 2.4 Key Insight

> Phase is invisible to measurement in isolation, but when two amplitude paths
> combine (via another H or rotation), phase determines whether they add or
> cancel. This is *interference*, and it is the engine of quantum speedup.

## 2.5 Self-Check

1. Does `Z 0` change the probability of measuring 1? Verify with `DISPLAY`.
2. What does `H 0` → `Z 0` → `H 0` do to `|0>`? Predict first, then run.
3. What does `H 0` → `S 0` → `H 0` do to `|0>`? (S is a 90° phase.)

---

---

# MODULE 3 — Two Qubits, Entanglement, and CX

## 3.1 Classical Intuition: State Space Explosion

With 1 bit you have 2 states (0 or 1). With 2 bits you have 4 states (00, 01,
10, 11). With n bits you have 2^n states — an exponential explosion in the
number of possible configurations.

A classical computer with n bits holds *one* of those 2^n configurations at any
time. An n-qubit quantum state holds a *superposition of all 2^n configurations
simultaneously*, each with its own complex amplitude. This is the origin of
quantum parallelism — but as we will see, reading it out collapses it back to
one classical value. The trick is to arrange interference so the *right* answer
has high amplitude when you measure.

## 3.2 The CX (CNOT) Gate

The **CX gate** (Controlled-X, also called CNOT) is the fundamental two-qubit
gate. It flips the *target* qubit if and only if the *control* qubit is `|1>`:

```
|00> → |00>
|01> → |01>
|10> → |11>
|11> → |10>
```

As a classical gate, this is XOR: target = target ⊕ control. In the quantum
setting it applies to superpositions, creating **entanglement**.

## 3.3 Entanglement

Prepare a Bell state:

```
INIT 2
H 0
CX 0 1
DISPLAY
```

You should see two entries: `|00>` with amplitude ≈0.707, and `|11>` with
amplitude ≈0.707. There is no entry for `|01>` or `|10>`.

This state — `(|00> + |11>) / √2` — is **entangled**. It cannot be written as
`(α₀|0> + β₀|1>) ⊗ (α₁|0> + β₁|1>)` for any choice of α and β. The two
qubits are not independent. Measuring one immediately determines the other.

**Classical correlation vs entanglement:** A classical bit pair could be (0,0)
or (1,1) with equal probability and still show correlated measurements. The
quantum difference is that:
- The superposition is a single coherent mathematical object, not a mixture.
- The amplitudes can interfere.
- Measurements in *other bases* (not just 0/1) reveal correlations that no
  classical probability distribution can reproduce (Bell inequalities — see
  `chsh_bell_test.md`).

## 3.4 Simulator Experiment: All Four Bell States

The four Bell states are the canonical maximally-entangled two-qubit states:

```
INIT 2
H 0
CX 0 1
DISPLAY
```

Expected: |Φ+> = `(|00> + |11>) / √2`

```
INIT 2
H 0
CX 0 1
Z 0
DISPLAY
```

Expected: |Φ-> = `(|00> - |11>) / √2`

```
INIT 2
H 0
CX 0 1
X 1
DISPLAY
```

Expected: |Ψ+> = `(|01> + |10>) / √2`

```
INIT 2
H 0
CX 0 1
X 1
Z 0
DISPLAY
```

Expected: |Ψ-> = `(|01> - |10>) / √2`

Verify correlations: measure both qubits on |Φ+>:

```
INIT 2
H 0
CX 0 1
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

Repeat several times. The outcomes are always 00 or 11 — never 01 or 10.

## 3.5 Key Insight

> An n-qubit state lives in a 2^n-dimensional complex vector space. A CX gate
> applied to a superposition creates entanglement — joint states that carry
> information in their *correlations*, not just in individual qubits. Measuring
> one entangled qubit instantly collapses the joint state.

## 3.6 Self-Check

1. For |Ψ+>, repeat the `MEASURE 0 0` then `MEASURE 1 1` experiment several
   times. Are the outcomes always different?
2. Why can't a 10-qubit state be described by 10 independent qubits? How many
   complex amplitudes are needed in general?
3. Starting from |Φ->, apply `H 0` and `H 1`, then measure. Are the outcomes
   correlated or anti-correlated? Why? (See `bell_states.md` for the answer.)

---

---

# MODULE 4 — Measurement and the Classical Register

## 4.1 Classical Intuition: Sampling

Measurement is how quantum information becomes classical information. It is
*destructive* and *probabilistic*: once you measure a qubit, the superposition
collapses to whichever outcome was sampled, weighted by the squared amplitudes.

Think of it like sampling from a probability distribution — but the distribution
is not given to you in advance; it is *defined by* the quantum state. And once
you sample, the state *becomes* the sample.

## 4.2 The Born Rule

The probability of measuring basis state `|x>` is:

```
P(x) = |amplitude of |x> |²
```

This is called the **Born rule**. It is the only way quantum amplitudes become
observable probabilities.

Note that the probabilities must sum to 1 over all basis states — this is why
valid quantum states always have unit-norm amplitude vectors.

## 4.3 Partial Measurement and Entanglement Collapse

When you measure *one qubit* of an entangled pair, the joint state collapses
*consistently*. If you measure qubit 0 and get 1, the remaining state (qubit 1)
is instantly in the state consistent with that outcome — regardless of how far
apart the qubits are.

Verify this:

```
INIT 2
H 0
CX 0 1
MEASURE 0 0
DISPLAY
```

After measuring qubit 0, `DISPLAY` shows the remaining state has only *one*
basis state with nonzero amplitude — the one consistent with your measurement.

## 4.4 Classical Registers

The simulator maintains a classical register alongside the quantum state. Each
`MEASURE j c` operation:
1. Samples qubit `j` according to Born-rule probabilities.
2. Collapses qubit `j` to the measured value.
3. Stores the result in classical bit `c`.

You read the classical register with `DISPLAY` (it appears at the bottom of the
output) or access individual bits in C++ via `get_cbit(c)`.

## 4.5 Simulator Experiment: Measurement Statistics

Use the QRNG (quantum random number generator) to see raw measurement statistics:

```
QRNG 8
```

This measures 8 independent qubits (each prepared in superposition), giving 8
truly random classical bits. Re-run multiple times to see the randomness.

For a controlled experiment with a biased state:

```
INIT 1
RY 0 1.047
DISPLAY
```

`RY 0 1.047` is a rotation by ≈ 60°, giving P(1) ≈ 0.25. Measure:

```
MEASURE 0 0
DISPLAY
```

Repeat 10 times and observe roughly 25% ones.

## 4.6 Key Insight

> Measurement is the *only* interface between the quantum and classical worlds.
> It is irreversible, probabilistic, and collapses the state. All quantum
> algorithms must be designed so that the *useful* answer has high probability
> amplitude at the moment of measurement.

## 4.7 Self-Check

1. A state has amplitudes `(0.6, 0.8i)`. What are the measurement probabilities?
2. If you measure qubit 0 of |Ψ+> and get 0, what state is qubit 1 in?
3. Why can you not use quantum measurement to transmit information faster than
   light? (Hint: what controls *which* outcome you get?)

---

---

# MODULE 5 — Interference: Why Quantum Algorithms Work

## 5.1 The Core Mechanism

Quantum speedups do not come from "trying all inputs simultaneously." If they
did, measurement — which collapses the state to one random outcome — would
erase the advantage immediately.

The speedup comes from **interference**: carefully arranging amplitudes so that
wrong answers cancel out (destructive interference) and correct answers
reinforce (constructive interference) before measurement.

Every quantum algorithm follows this template:
1. **Prepare** a superposition of all possible inputs.
2. **Apply** an operation that encodes the problem (the oracle).
3. **Interfere**: use gates to route amplitude toward the right answer.
4. **Measure**: read out the amplified answer with high probability.

## 5.2 Deutsch-Jozsa: The Simplest Demonstration

The Deutsch-Jozsa problem is deliberately toy-like, but it is the clearest
example of quantum advantage from interference.

**Problem:** You have a black-box function f: {0,1}^n → {0,1} that is *promised*
to be either:
- **Constant**: f(x) = 0 for all x, or f(x) = 1 for all x.
- **Balanced**: f(x) = 0 for exactly half of all inputs, and 1 for the other half.

**Classical cost:** In the worst case (deterministically), you need 2^(n-1)+1
queries to be certain. Even with randomness, you need O(1) queries but cannot
guarantee correctness.

**Quantum cost:** 1 query. Always.

### Why 1 Query is Enough

The algorithm prepares *all inputs simultaneously* in superposition, queries the
oracle *once* (which applies the function to all inputs at once via linearity),
and then uses Hadamard gates to interfere the amplitudes:

- If f is **constant**: all phase shifts align the same way. After interference,
  the state returns to `|00...0>`. Measuring all zeros tells you "constant."

- If f is **balanced**: the phase shifts are +1 for half the inputs, -1 for the
  other half. These cancel at `|00...0>`. Measuring reveals at least one 1 (it
  did not return to all-zeros). This tells you "balanced."

The oracle query is being used not to check one input, but to *simultaneously*
adjust the phases of all 2^n inputs. Interference then reads out the global
property (constant vs balanced) in one shot.

## 5.3 Simulator Experiment: Constant vs Balanced

```
DEUTSCH_JOZSA 3 CONST0
```

Expected: measured register is `000` (all zeros = constant).

```
DEUTSCH_JOZSA 3 CONST1
```

Expected: again `000`.

```
DEUTSCH_JOZSA 3 BALANCED_XOR0
```

Expected: measured register has at least one `1` (= balanced).

```
DEUTSCH_JOZSA 3 BALANCED_PARITY
```

Expected: at least one `1`.

### Build it by hand (n=2)

This manual version shows exactly what is happening at each step:

```
INIT 3
X 2
H 0
H 1
H 2
DISPLAY
```

After this: uniform superposition over all input combinations, ancilla in (|0>-|1>)/√2
(the phase-kickback trick that encodes the oracle output as amplitude phase).

Apply the oracle (BALANCED_XOR0 means: flip ancilla when q0=1):

```
CX 0 2
DISPLAY
```

Notice the phase pattern has changed — different inputs have different phases.
Now interfere (Hadamard on input qubits):

```
H 0
H 1
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

At least one measured bit should be 1.

## 5.4 Phase Kickback (The Oracle Mechanism)

The ancilla qubit prepared in `(|0>-|1>)/√2` is the standard trick for turning
an oracle's output into a *phase* on the input state. When the oracle flips the
ancilla (f(x)=1), the effect is to multiply the input amplitude by -1. This is
**phase kickback**: the classical output is "kicked back" into the phase of the
input.

This is used in *every* oracle-based quantum algorithm. The ancilla absorbs the
function output and returns unchanged; only the input amplitudes carry the mark.

## 5.5 Key Insight

> Quantum speedup is not about parallelism — it is about interference. The oracle
> marks a *global property* (constant/balanced) as a phase pattern. Hadamard
> gates then convert that pattern into a *measurable outcome*. This requires
> only 1 oracle call regardless of n.

## 5.6 Self-Check

1. For CONST0, why does measuring the input register always give `00...0`?
2. For BALANCED_XOR0 with n=3, which input values have f(x)=1?
3. Why would a classical randomized algorithm need more than 1 query to be
   *certain* of the answer?

---

---

# MODULE 6 — Oracle Algorithms: Query Complexity

## 6.1 The Pattern

Oracle algorithms exploit quantum superposition and phase kickback to extract
global properties of a function using exponentially fewer queries than any
classical algorithm requires.

The general pattern:
1. Prepare an input register in uniform superposition.
2. Query the oracle once (all inputs in parallel).
3. Use interference to convert the global property into a measurable bit pattern.

## 6.2 Bernstein-Vazirani: Hidden Linear Structure

**Problem:** A secret bitstring `s` is hidden inside a linear oracle:

```
f(x) = s · x  (mod 2)    [dot product of binary vectors]
```

**Classical cost:** n queries — one per bit of s (query e₀, e₁, ..., eₙ₋₁).

**Quantum cost:** 1 query.

**Why:** The query encodes all n bits of s simultaneously into the phases of
the superposed input states. One Hadamard decode reads all n bits at once.

Run the demo:

```
BV 5 22
```

The binary representation of 22 is `10110`. The recovered secret should be `22`.
This used one oracle query to determine 5 bits simultaneously.

### By hand (n=4, secret=1011 → decimal 11)

```
INIT 5 5
X 4
H 0
H 1
H 2
H 3
H 4
CX 0 4
CX 1 4
CX 3 4
H 0
H 1
H 2
H 3
MEASURE 0 0
MEASURE 1 1
MEASURE 2 2
MEASURE 3 3
DISPLAY
```

Expected classical bits (LSB first): `1`, `1`, `0`, `1` → secret = `1011` = 11.

## 6.3 Simon's Algorithm: Hidden XOR Structure

**Problem:** An oracle f has a hidden *xor mask* s such that:

```
f(x) = f(x ⊕ s)    for all x
```

If s=0, f is one-to-one. Otherwise f is two-to-one (pairs differing by s give
the same output).

**Classical cost:** O(√(2^n)) random queries before a collision is expected
(birthday paradox).

**Quantum cost:** O(n) queries. Exponentially fewer.

**Why:** Each quantum query generates a random vector y satisfying `y · s = 0 (mod 2)`.
After n such vectors, Gaussian elimination over GF(2) reveals s. Classically,
you would need to find a collision; quantumly, you extract structural equations.

```
SIMON 4 10 8
```

Secret `10` in decimal is `1010` in binary. The algorithm recovers it from
8 sampled equations over GF(2).

## 6.4 Classical vs Quantum Query Complexity (Summary)

| Problem         | Classical (worst) | Quantum   | Mechanism        |
|-----------------|-------------------|-----------|------------------|
| Deutsch-Jozsa   | 2^(n-1)+1         | 1         | Interference     |
| Bernstein-Vazirani | n              | 1         | Phase + Hadamard |
| Simon's          | O(√(2^n))        | O(n)      | Equation sampling|
| Grover (search)  | O(N)             | O(√N)     | Amplitude amp.   |
| Shor (factoring) | Subexponential   | O(n³)     | QFT + period finding |

## 6.5 Key Insight

> Oracle query complexity separates the number of times you *call* the function
> from the amount of *computation* inside each call. Quantum computers do not
> speed up arbitrary computation — they exploit *structure* (linearity, periodicity,
> amplitude concentration) to need exponentially or polynomially fewer function
> calls.

## 6.6 Self-Check

1. In BV, why does using n=10 instead of n=5 still require only 1 oracle query?
2. What kind of structure does Simon's algorithm exploit that Grover does not?
3. Deutsch-Jozsa provides an exponential speedup over classical deterministic
   algorithms. Is this a "practical" speedup? Why or why not?

---

---

# MODULE 7 — Amplitude Amplification (Grover's Algorithm)

## 7.1 Classical Intuition: Probabilistic Search

Suppose you have an unsorted list of N items and you want to find the one item
that satisfies a predicate f(x) = 1. Classically, random search requires O(N)
calls to f on average. There is no way to do better than O(N) without structure.

Grover's algorithm finds the marked item in O(√N) oracle calls. This is an
unconditionally proven *quadratic speedup* for unstructured search.

For N = 10^12, classical search needs ~ 10^12 calls; Grover needs ~ 10^6.

## 7.2 How Grover Works: Amplitude as a Weight

Think of each basis state `|x>` as having a "weight" — its amplitude. Initially
all weights are equal (uniform superposition). The Grover iteration does two
things in each round:

**Step 1 — Oracle phase flip:** The oracle negates the amplitude of the marked
state(s). Visually: the correct answer's weight flips sign.

**Step 2 — Diffusion (reflection about the mean):** Every amplitude is
reflected about the *average* amplitude. This is the key step. Because the
marked state was flipped negative, it pulls the average slightly down. After
reflection, the marked state's amplitude is now larger than before (it bounced
above the average), and all others are slightly smaller.

After √(N/M) iterations (where M is the number of marked states), the marked
state has near-100% amplitude.

**What goes wrong with too many iterations:** The amplitude oscillates. Past the
optimal point, the marked state's amplitude starts decreasing again. This is why
iteration count matters.

## 7.3 Simulator Experiment: Grover Search

Basic search (3 qubits, N=8, one target):

```
GROVER 5
```

This searches for a single target (5) in 8 possibilities. Run with verbose mode
to watch the probability grow iteration by iteration.

```
VERBOSE VERBOSE
GROVER 5
```

Multi-target search:

```
GROVER 3 5
```

Auto-tuned search (count targets first, then set optimal iterations):

```
GROVER AUTO 3 6 3 5
```

### Manual Grover (2 qubits, target = |11>)

Step through the algorithm by hand to see each piece:

```
INIT 2
H 0
H 1
DISPLAY
```

All 4 states have equal amplitude (0.5). Apply oracle (phase-flip |11> via CZ):

```
CZ 0 1
DISPLAY
```

|11> now has amplitude -0.5; others still +0.5. Apply diffusion:

```
H 0
H 1
X 0
X 1
CZ 0 1
X 0
X 1
H 0
H 1
DISPLAY
```

|11> should now have amplitude ≈ 1.0 (one iteration is exact for N=4, M=1).

## 7.4 Key Insight

> Grover's algorithm achieves a quadratic speedup over classical unstructured
> search — provably the best possible for this problem class. It works by
> iteratively rotating the amplitude vector toward the target(s). The oracle
> marks winners; diffusion amplifies them. The O(√N) bound is tight.

## 7.5 Practical Notes

- Grover is provably optimal for unstructured search (no quantum algorithm can
  do better than O(√N) for this problem).
- If you know how many marked items there are (M), the optimal number of
  iterations is approximately π√(N/M)/4.
- If M is unknown, use **Quantum Counting** (next module) first.
- Grover's speedup is quadratic, not exponential. It is useful but does not
  break NP-hard problems in general.

## 7.6 Self-Check

1. For N=1024 (10 qubits, 1 target), approximately how many Grover iterations
   are needed? (Hint: ≈ π/4 × √N)
2. What happens if you run 100 iterations when only 10 are optimal? Try it
   with GROVER on a small example and watch the probability.
3. Why does Grover not solve NP-complete problems in polynomial time?

---

---

# MODULE 8 — Quantum Counting and the Grover + Counting Workflow

## 8.1 Problem: Unknown Target Count

Grover's algorithm requires knowing how many marked states (M) exist to set the
optimal iteration count. Quantum Counting estimates M using a different technique
(Quantum Phase Estimation applied to the Grover operator), which is itself
interesting.

## 8.2 Quantum Counting in Practice

```
QCOUNT RUN 3 6 1 5
```

This counts marked states `{1, 5}` in a 3-qubit space using 6 phase estimation
qubits. The output is an estimate of M (should be ≈2).

The workflow for unknown-M search:

```
QCOUNT RUN 3 6 3 5
GROVER 3 5
```

Or the all-in-one:

```
GROVER AUTO 3 6 3 5
```

## 8.3 Lab: Counting and Grover Together

```
VERBOSE VERBOSE
QCOUNT RUN 4 8 0 3 5 12
GROVER 0 3 5 12
```

Observe:
- Counting phase: QCOUNT estimates M=4 (four targets in a 16-state space).
- Grover phase: algorithm runs the appropriate number of iterations.

## 8.4 Key Insight

> Quantum Counting is an application of Quantum Phase Estimation — one of the
> most widely used primitives in quantum algorithms. It gives you the iteration
> budget for Grover when you do not know M, converting a parameter-sensitive
> algorithm into a parameter-free one.

---

---

# MODULE 9 — The Quantum Fourier Transform

## 9.1 Classical Intuition: FFT

The Fast Fourier Transform (FFT) converts a signal from the time domain to the
frequency domain in O(n 2^n) operations. It is used everywhere in signal
processing, compression, and numerical computation.

The **Quantum Fourier Transform** (QFT) does the same mathematical operation on
the *amplitude vector* of a quantum state. In the quantum setting, it requires
only O(n²) gates (where n is the number of qubits, not 2^n amplitudes). This
is an exponential reduction in circuit size — but with a critical catch: you
cannot directly *read out* the transformed amplitudes without measurement
destroying them.

The QFT is not useful on its own. Its power is as a subroutine: it allows
algorithms to detect *periodicity* in quantum states, which is the key to
Shor's algorithm.

## 9.2 What the QFT Does

Applied to a state that has periodic amplitude structure with period r, the QFT
concentrates the amplitude at states that are integer multiples of N/r (where
N=2^n). Measuring then gives a clue about r.

## 9.3 Simulator Experiment: QFT on Uniform Superposition

```
INIT 4
QFT 0 3
DISPLAY
```

Starting from `|0000>` (which is a uniform superposition in the "frequency domain"),
the QFT produces a single peak at `|0000>` in the "time domain." This is the
quantum analog of "DC component only."

Try a basis state input:

```
INIT 4
X 0
QFT 0 3
DISPLAY
```

The QFT of a single basis state creates a uniform superposition with a phase
gradient — like the DFT of a single spike.

Verify QFT is its own inverse (up to ordering):

```
INIT 4
H 0
H 1
H 2
H 3
QFT 0 3
IQFT 0 3
DISPLAY
```

Should return to the uniform superposition.

## 9.4 Key Insight

> The QFT is exponentially faster than the FFT (O(n²) gates vs O(n·2^n)
> classical operations) but you cannot read out all 2^n Fourier coefficients
> directly. Its value is as a subroutine for detecting periodicity in quantum
> amplitude patterns — the foundation of Shor's algorithm.

## 9.5 Self-Check

1. How many QFT operations would it take to transform an n-qubit state
   classically (via FFT)? How many gates does the quantum QFT use?
2. Why can you not use the QFT to compute FFTs faster on classical data?
   (Hint: what is the cost of loading classical data into a quantum state?)
3. The QFT and IQFT are inverses. Verify this with a 3-qubit state of your
   choice.

---

---

# MODULE 10 — Period Finding and Shor's Algorithm

## 10.1 Classical Intuition: Factoring via Period Finding

The RSA cryptosystem's security rests on the hardness of factoring large
integers. The best classical algorithms (general number field sieve) run in
sub-exponential but super-polynomial time.

Shor's key insight: **factoring reduces to period finding**, and period finding
can be done exponentially faster with a quantum computer.

The reduction: to factor N, pick a random `a` coprime to N and find the period r
of the function `f(x) = a^x mod N`. Once you have r:
- If r is even and `a^(r/2) ≢ -1 (mod N)`, then `gcd(a^(r/2) ± 1, N)` gives
  a non-trivial factor.

Finding r classically requires trying many values of x — essentially as hard as
factoring itself. The QFT makes period finding efficient.

## 10.2 The Quantum Circuit for Period Finding

1. **Initialize** a control register in uniform superposition and a target
   register in `|1>`.
2. **Apply** controlled modular exponentiation: `|x>|y> → |x>|a^x mod N>`.
   This entangles the control and target registers.
3. **QFT** on the control register. The entanglement ensures the QFT output
   peaks at multiples of `2^n / r`.
4. **Measure** the control register. Get a random multiple of `2^n / r`.
5. **Classical post-processing**: use continued fractions to extract r from the
   measured value.

## 10.3 Simulator Experiment

```
SHOR 15
```

Expected: the demo finds that 15 = 3 × 5 or 15 = 5 × 3.

```
VERBOSE VERBOSE
SHOR 15
```

With verbose mode, you can watch the quantum order-finding and the classical
factor extraction steps.

Try other small composites:

```
SHOR 21
SHOR 33
SHOR 35
```

Note: the simulator is educational scale. Real Shor's algorithm would require
thousands of (error-corrected) qubits to factor cryptographically relevant
numbers.

## 10.4 Why This Matters

Shor's algorithm runs in polynomial time O(n³) in the number of bits n=log₂(N).
The best known classical algorithm is O(exp(n^(1/3))) — subexponential but
still faster than polynomial growth at large n. For 2048-bit RSA keys:
- Classical: ~10^14 years (at current hardware speeds).
- Shor (on a large fault-tolerant quantum computer): hours.

This is why post-quantum cryptography is an active research area.

## 10.5 Key Insight

> Shor's algorithm combines two quantum primitives: superposition to evaluate
> a^x mod N for all x simultaneously, and the QFT to extract the period r from
> the resulting amplitude pattern. The period determines the factors via simple
> classical number theory. The algorithm runs in polynomial time — an
> exponential speedup over the best known classical approach.

## 10.6 Self-Check

1. For N=15, a=7: what is the period r of 7^x mod 15? (Check: 7^1=7, 7^2=49≡4,
   7^4≡1. So r=4.) Does gcd(7^2-1, 15) give a non-trivial factor?
2. Why does Shor's algorithm not work if r is odd?
3. The quantum part of Shor's algorithm finds a multiple of 2^n/r. How does
   classical continued-fraction approximation extract r from this?

---

---

# MODULE 11 — Quantum Communication Protocols

## 11.1 Superdense Coding

Entanglement can be used to transmit 2 classical bits by sending only 1 qubit
(plus a pre-shared Bell pair). This protocol — superdense coding — shows that
entanglement is a resource that can be "consumed" to increase channel capacity.

Alice and Bob share |Φ+>. Alice encodes 2 classical bits by applying a local
gate to her qubit, then sends that single qubit to Bob. Bob decodes with CNOT
and H.

```
INIT 2
H 0
CNOT 0 1
Z 0
CNOT 0 1
H 0
MEASURE 0 0
MEASURE 1 1
DISPLAY
```

Expected: classical bits = `10` (the message Alice encoded with `Z 0`).

See `superdense_coding.md` for all four message cases.

## 11.2 Quantum Teleportation

Teleportation transmits an *unknown quantum state* from Alice to Bob using a
Bell pair and 2 classical bits. Unlike superdense coding, the payload is quantum
(the state being teleported could be any qubit state, including one you don't
know).

Key points:
- The original state is *destroyed* by Alice's measurement (no-cloning theorem).
- No information travels faster than light: Alice must send her 2 classical bits
  to Bob before he can complete the correction.
- Teleportation does not transmit *matter* — only quantum information.

See `quantum_teleportation.md` for the full CLI walkthrough.

## 11.3 Key Insight

> Entanglement is a physical resource. Superdense coding uses it to double
> classical channel capacity; teleportation uses it to transmit quantum states
> using only classical channels. Both protocols are foundational to quantum
> networking and quantum cryptography research.

---

---

# MODULE 12 — Variational Quantum Algorithms (QAOA and VQE)

## 12.1 Classical Intuition: Parameterized Optimization

Variational algorithms are the current NISQ-era (Noisy Intermediate-Scale
Quantum) approach to quantum utility. The idea is:

1. Define a **parameterized quantum circuit** (called the *ansatz*) with tunable
   rotation angles.
2. Evaluate the **expected energy** of the current parameters by running the
   circuit and measuring many times.
3. **Optimize** the parameters using a classical optimizer (coordinate descent,
   BFGS, COBYLA, etc.).
4. Repeat until convergence.

From a classical engineer's perspective: this is derivative-free optimization
of an expensive black-box function. The quantum circuit is the function
evaluator; the classical optimizer is the outer loop.

## 12.2 QAOA: Quantum Approximate Optimization Algorithm

QAOA targets **combinatorial optimization** problems encoded as QUBO (Quadratic
Unconstrained Binary Optimization):

```
minimize: x^T Q x    where x ∈ {0,1}^n
```

QUBO captures many important problems: Max-Cut, portfolio optimization, Sudoku,
traveling salesman (see `max_cut.md`, `tsp_usage.md`).

QAOA alternates between two parameterized operations (p layers):
- **Cost layer**: phase shift proportional to the objective value.
- **Mixing layer**: rotation across all bit values (typically RX gates).

The `γ` parameters control cost emphasis; `β` parameters control mixing.

### Running QAOA

```
QUBO EXACT 3 -2 0 2 0 1 0 2 0 -3
```

Find the exact minimum first (use as baseline):

```
VQA QAOA 3 1 0 40 0.25 -2 0 2 0 1 0 2 0 -3
```

- `n=3` variables, `p=1` layer, exact expectation (`shots=0`), 40 optimizer
  iterations, step size 0.25.

Compare with increasing depth:

```
VQA QAOA 3 2 0 60 0.20 -2 0 2 0 1 0 2 0 -3
```

More layers generally get closer to the exact minimum — but also require more
parameters to tune and more circuit depth (more noise on real hardware).

## 12.3 VQE: Variational Quantum Eigensolver

VQE targets **eigenvalue problems** — finding the ground state energy of a
quantum Hamiltonian. This is important for quantum chemistry (finding molecular
ground states) and condensed matter physics.

The Hamiltonian is expressed as a sum of Pauli operators (X, Y, Z). The VQE
optimizer minimizes the expected energy by tuning the ansatz circuit parameters.

### Single-qubit ground state

The ground state energy of H = Z₀ is -1 (eigenvalue of Z for `|1>`):

```
VQE RUN 1 1 30 0.3 0 1 1.0 1 Z 0
```

Expected: `best_energy ≈ -1`.

### Two-qubit Hamiltonian

```
VQE RUN 2 2 50 0.2 0 3 1.0 1 Z 0 1.0 1 Z 1 0.5 2 X 0 X 1
```

This minimizes the energy of `H = Z₀ + Z₁ + 0.5 X₀X₁`. Increase `layers` and
`iters` to approach the true ground state.

## 12.4 Simulated Annealing: Classical Baseline

The simulator also includes simulated annealing (SA) and simulated quantum
annealing (SQA) for solving QUBO classically. These are useful baselines:

```
ANNEAL SA 3 -2 0 2 0 1 0 2 0 -3
ANNEAL SQA 3 -2 0 2 0 1 0 2 0 -3
```

Compare results and runtimes with `VQA QAOA`. On small instances like this,
the classical annealer will typically find the same minimum — QAOA's advantage
(if any) appears at larger scale or with specific circuit structures.

## 12.5 Key Insight

> Variational algorithms are the current best strategy for NISQ-era quantum
> computers because they are noise-tolerant (short circuits, classical outer
> loop absorbs errors). Their advantage over classical methods is unproven
> at scale — this is an active research frontier. Treat them as quantum-assisted
> heuristics, not guaranteed speedups.

## 12.6 Self-Check

1. What is the difference between QAOA and VQE in terms of the problem they
   target?
2. Why does increasing the number of QAOA layers (p) generally improve solution
   quality but have diminishing returns?
3. Run both `ANNEAL SA` and `VQA QAOA` on the same 3-variable QUBO. Do they
   find the same minimum? Which runs faster?

---

---

# MODULE 13 — What Quantum Can and Cannot Do

## 13.1 Proven Quantum Advantages

These are cases where quantum computers are provably faster than *any*
classical algorithm, not just known ones:

| Problem | Classical | Quantum | Proven? |
|---------|-----------|---------|---------|
| Unstructured search (Grover) | O(N) | O(√N) | Yes (tight) |
| Integer factoring (Shor) | Subexponential | Polynomial | Quantum yes; classical lower bound unknown |
| Query complexity (DJ, BV, Simon) | Exponential | Polynomial | Yes (in query model) |
| Simulation of quantum systems | Exponential | Polynomial | Yes (this is the original motivation) |

## 13.2 Not Proven Quantum Advantages

- **NP-complete problems**: Grover gives a quadratic speedup but NP-complete
  problems remain exponential. Quantum computers do not solve SAT, TSP, or
  graph coloring in polynomial time (unless BQP = NP, which is not believed).
- **QAOA/VQE**: Advantages over classical methods are not proven at scale.
  These are heuristics.
- **Machine learning**: Claims of "quantum ML" speedup are often conditional
  on input assumptions (data in quantum form) that are hard to satisfy.

## 13.3 Key Limitations for Engineers

**Decoherence and noise:** Real qubits lose their quantum state quickly.
Current hardware has error rates of 0.1–1% per gate. Useful Shor's algorithm
requires fault-tolerant qubits with error rates below ~10⁻¹⁵ — requiring
thousands of physical qubits per logical qubit. We are not there yet.

**Input/output bottleneck:** Loading n classical bits into a quantum state
takes O(2^n) operations in the worst case. Many proposed quantum speedups
assume you already have data "in quantum form," which is often unrealistic.

**No quantum RAM yet:** The QRAM hypothesis (efficient quantum random access
memory) underlies many quantum speedup claims. QRAM does not exist at useful
scale.

**Measurement collapses state:** You get one sample per circuit run. To
estimate probabilities, you must run thousands of times (shots). This overhead
is often ignored in theoretical analyses.

## 13.4 Where Quantum Will Matter First

Based on current trajectory, quantum advantage is most credible for:

1. **Quantum simulation**: Simulating molecules and materials — the original
   motivation (Feynman, 1982). Drug discovery, materials science.
2. **Cryptography**: Shor's algorithm will eventually break RSA/ECC. Post-quantum
   cryptography (CRYSTALS, NTRU, etc.) is being standardized now.
3. **Structured optimization**: Problems with exploitable periodicity or
   entanglement structure (not generic NP-hard instances).

## 13.5 Key Insight

> Quantum computers are not magic. They offer provable advantages for specific
> problem structures: unstructured search (quadratic), period-finding
> (exponential), and quantum simulation (exponential). They do not solve NP-hard
> problems efficiently and are still limited by hardware noise. The right question
> is not "is my problem hard?" but "does my problem have the structure that
> quantum algorithms exploit?"

---

---

# MODULE 14 — Practical Workflow and Next Steps

## 14.1 Decision Flowchart: Should You Use a Quantum Algorithm?

```
Does my problem have UNSTRUCTURED SEARCH as its core?
  Yes → Grover (quadratic speedup)
  No →
    Does it have PERIODICITY or HIDDEN LINEAR STRUCTURE?
      Yes → QFT-based (Shor, Simon, BV)
      No →
        Is it COMBINATORIAL OPTIMIZATION?
          Yes → QAOA/VQE (heuristic, unproven advantage)
          No →
            Is it QUANTUM SYSTEM SIMULATION?
              Yes → Quantum simulation (proven advantage)
              No → Use classical algorithms
```

See `algorithm_comparison.md` for a more detailed selection guide.

## 14.2 The C++ API

For programmatic use, the simulator exposes a clean C++ API:

```cpp
#include "include/quantum_sim.hh"

State s(3, 3);        // 3 qubits, 3 classical bits
s.h(0).cx(0, 1).cx(0, 2);  // GHZ state
s.measure(0, 0).measure(1, 1).measure(2, 2);
s.display();
```

See `api_quickstart.md` for build instructions and Grover/Shor/QAOA API examples.

## 14.3 Recommended Learning Path After This Course

1. **Nielsen & Chuang** — "Quantum Computation and Quantum Information" (the
   standard textbook). Chapters 1–6 cover everything in this course with full
   mathematical treatment.
2. **IBM Quantum Learning** — hands-on with real hardware (Qiskit).
3. **Preskill's lecture notes** (available free online) — more theoretical depth.
4. **arXiv quant-ph** — follow active research. Variational algorithms, error
   correction, and quantum advantage are the hot topics.

## 14.4 Labs in This Repository

Work through these in order for structured practice:

| Lab | File | Skills |
|-----|------|--------|
| Bell States | `lab_bell_states.md` | Superposition, entanglement, measurement |
| Grover + Counting | `lab_grover_counting.md` | Amplitude amplification, QCOUNT |
| Shor (small) | `lab_shor.md` | QFT, period finding, classical post-processing |
| QAOA | `lab_qaoa.md` | Variational optimization, QUBO encoding |
| VQE | `lab_vqe.md` | Hamiltonian minimization, ansatz tuning |

---

---

# Appendix A — Quick Reference: CLI Commands

| Command | Description |
|---------|-------------|
| `INIT n` | Initialize n-qubit state to `\|00...0>` |
| `INIT n m` | n qubits, m classical bits |
| `H j` | Hadamard on qubit j |
| `X j` | Pauli-X (NOT) on qubit j |
| `Y j` / `Z j` | Pauli-Y / Pauli-Z on qubit j |
| `S j` / `T j` | Phase gates (90° / 45°) on qubit j |
| `RX j θ` / `RY j θ` / `RZ j θ` | Rotation gates (radians) |
| `CX j k` / `CNOT j k` | Controlled-X: control j, target k |
| `CZ j k` | Controlled-Z |
| `CCX c1 c2 t` | Toffoli gate |
| `SWAP j k` | Swap qubits j and k |
| `MEASURE j c` | Measure qubit j into classical bit c |
| `DISPLAY` | Show full state vector |
| `DISPLAY ALL` | Show all 2^n basis states (including zeros) |
| `QFT j k` | Quantum Fourier Transform on qubits j..k |
| `IQFT j k` | Inverse QFT on qubits j..k |
| `GROVER targets...` | Run Grover search |
| `GROVER AUTO n count targets...` | Auto-tuned Grover |
| `QCOUNT RUN n pebits targets...` | Quantum counting |
| `DEUTSCH_JOZSA n oracle` | Deutsch-Jozsa demo |
| `BV n secret` | Bernstein-Vazirani |
| `SIMON n secret` | Simon's algorithm |
| `SHOR N` | Shor's factoring demo |
| `QUBO EXACT n entries...` | Brute-force QUBO solve |
| `VQA QAOA n p shots iters step entries...` | QAOA optimizer |
| `VQE RUN n layers iters step shots terms...` | VQE optimizer |
| `ANNEAL SA n entries...` | Simulated annealing |
| `ANNEAL SQA n entries...` | Simulated quantum annealing |
| `CHECK NORMALIZED` | Verify state has unit norm |
| `CHECK TARGETS t...` | Show probability on listed targets |
| `SEED n` | Set RNG seed for reproducibility |
| `VERBOSE VERBOSE` | Enable detailed algorithm logs |
| `TUTOR ON` | Enable step-by-step teaching narration |

---

# Appendix B — Glossary (Classical-First)

| Term | Classical Analogue | Quantum Meaning |
|------|--------------------|-----------------|
| Qubit | Bit | A two-state system described by complex amplitudes |
| Superposition | Parallel processing (rough analogy) | A coherent combination of multiple basis states |
| Amplitude | Square root of probability (signed/complex) | Complex weight on each basis state |
| Phase | Signal phase angle | Imaginary rotation of amplitude; invisible to measurement alone |
| Interference | Signal add/cancel | Amplitude paths reinforcing or cancelling |
| Entanglement | Correlated random variables (but stronger) | Joint state that cannot be factored into independent parts |
| Measurement | Sampling | Collapses quantum state to one classical outcome |
| Oracle | Black-box function call | A reversible unitary implementing f(x) |
| QFT | FFT | Fourier transform on amplitude vector; O(n²) gates |
| Ansatz | Parameterized model family | Parameterized circuit used in variational algorithms |
| QUBO | Integer programming (binary) | Quadratic binary objective: x^T Q x |
| Hamiltonian | Cost function / energy | Operator defining energy in physics-inspired algorithms |
| Decoherence | Noise / environmental coupling | Loss of quantum coherence due to interaction with environment |
| NISQ | Prototype hardware | Noisy Intermediate-Scale Quantum: current ~50–1000 qubit, high-error devices |

---

# Appendix C — Mathematics Cheat Sheet

**State vector:** A unit vector in ℂ^(2^n). For n=2, it has 4 complex entries.

**Inner product:** `<ψ|φ> = Σ conj(ψᵢ) φᵢ`. Probability = `|<x|ψ>|²`.

**Unitary gate:** A matrix U with `U†U = I`. All quantum gates are unitary.
This guarantees norm preservation (probabilities always sum to 1) and
reversibility (every gate has an inverse).

**Tensor product:** `|a>⊗|b>` is the joint state of two independent registers.
For `|0>⊗|1>` in 2-qubit notation: `|01>`.

**Hadamard matrix:**
```
H = (1/√2) [[1,  1],
             [1, -1]]
```

**Pauli gates:**
```
X = [[0,1],[1,0]]    Y = [[0,-i],[i,0]]    Z = [[1,0],[0,-1]]
```

**Born rule:** P(outcome x) = `|<x|ψ>|²` = squared magnitude of amplitude.

**No-cloning theorem:** You cannot copy an unknown quantum state. This has
deep implications for quantum cryptography and error correction.

---
