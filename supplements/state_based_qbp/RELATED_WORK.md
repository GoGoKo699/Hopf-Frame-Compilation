# Literature context for state-based QBP

[State-based task theorem](STATE_BASED_QBP_THEOREM.md) · [Core related work](../../docs/RELATED_WORK.md) · [Source catalogue](../../docs/reference/SOURCE_CATALOGUE.md)

This supplement has its own state, coherent-reference, and decoder contracts.
Its fine state preparation does not implement the prescribed complete frame.
The following literature comparisons retain their original section numbers
and contribution boundaries; they are not extra claims of the core A–D package.

## Reading index

- [State preparation and the changed gradient decoder](#13-state-only-preparation-and-a-separate-gradient-decoder)
- [Gauge-fixed complex extension](#14-gauge-fixed-complex-state-and-gradient-extension)
- [Certified bounded inputs and classical comparison](#15-bounded-input-construction-and-the-classical-comparator)

## 13. State-only preparation and a separate gradient decoder

The [state-only construction](STATE_ONLY_COMPILER.md) uses the retained
borrowed compiler for a coarse circuit, the full-operator rotation primitive
for a single correction table, and standard amplitude amplification.
[Brassard–Høyer–Mosca–Tapp, Section 2, Eq. (8)](https://arxiv.org/pdf/quant-ph/0005055)
gives the one-step success law at accepted amplitude one half. The local
proof supplies the actual coarse residual, two-flag word, coherent branch
extension, and literal dirty-space and error ledgers. At $`L=N`$ its linear
state-preparation T count agrees with the order of
[GKW's unrestricted-ancilla benchmark, Theorem 1.1](https://quantum-journal.org/papers/q-2026-07-22-2168/pdf/);
the additional local statement is the two-clean allocation and arbitrary
dirty/reference return. No complete-frame conclusion follows.

The [reference-state QBP protocol](REFERENCE_STATE_QBP.md) changes the
measurement and classical scores. Logarithmic wave-function derivatives
already appear in variational Monte Carlo, for example
[Toulouse–Umrigar, Eq. (45)](https://arxiv.org/pdf/physics/0701039).
Hadamard-test interference is also standard; see
[Mitarai–Fujii, Fig. 1](https://arxiv.org/pdf/1901.00015).
The general overlap identity is not a novelty claim. The local analysis
combines Hopf's one-node-per-depth leaf support with a positive derivative
envelope and the charged state-only compiler. A fixed lower bound on local
branch probabilities yields constant score overhead; arbitrary angles can
require an n-dependent factor in this sufficient sampling bound. This is a
specific alternative for raw gradients, not a generic training speedup or a
priority claim over all importance-sampling protocols.

The subsequent [coarse-frame decoder](COARSE_FRAME_QBP.md) combines the
same standard X/Y interference identity with a coarse full-frame word and
the original Walsh spreading step. Its local contribution is the
constant depth-record bound from coarse marker agreement, the coherent
reference matched to the actual native word, and the charged quantum and
classical reconstruction. Exact correction in the decoder removes coarse
bias without synthesizing a fine inverse frame. No claim of priority for
basis-dependent overlap estimation or classical adjoint accumulation is
made. The result covers all real angle tuples, including singular ones;
the complex phase-gradient stream remains separate.

## 14. Gauge-fixed complex state and gradient extension

The prefix phase cascade and its residual arithmetic-mean global phase
are inherited from [Möttönen–Vartiainen–Bergholm–Salomaa, Section III,
Eqs. (4), (5), and (7)](https://arxiv.org/pdf/quant-ph/0407010v1).
The [complex coarse compiler](COMPLEX_COARSE_COMPILER.md) fixes the sign
and half-angle convention, cancels the common phase for a state-based
task, and applies the retained literal reflection interpreter to each
determinant-one table. Its local result is the actual logical coarse C
with exact dirty return, linear T count, recorded possibly nondiagonal
rows, and the unchanged two-flag fine state-preparation reservation.

The [two-stream decoder](COMPLEX_COARSE_QBP.md) combines that interface
with the existing X/Y coarse-frame reconstruction and direct leaf-phase
records. It supplies the consistent gauge correction, arbitrary-angle
gradient identities, and complete quantum/classical precision ledger.
Neither the standard phase cascade nor global-phase invariance is a new
claim. The result does not compile the prescribed literal common phase of
the complete frame, prove universal gradient-cost improvement, or emit a
general native complex state-preparation package.

## 15. Bounded-input construction and the classical comparator

The [residual-table procedure](RESIDUAL_TABLE_PREPROCESSING.md) uses
ordinary half-angle algebra and the standard three-rotation SU(2)
factorization. Its local statement is a certified coefficient construction
for the particular residual completion: one shared half-phase preserves
the literal branch phase, finite rational decisions cover zero and boundary
cases, and direct sine/cosine programming retains the existing native error
and workspace contracts. No new general Euler decomposition is claimed.

The [bounded-input audit](BOUNDED_INPUT_QBP.md) uses GKW Lemma 2.3 only for
short-word existence. Exhaustive search at coarse precision gives a
conservative polynomial-in-N construction. This must be distinguished
from efficient fine synthesis: [Ross–Selinger](https://arxiv.org/abs/1403.2975v3)
give optimal synthesis with a factoring oracle and prove the efficient
expected runtime without that oracle under a number-theoretic hypothesis.
Neither is assumed in the bounded-input construction. Its residual fine
precision instead comes from explicit source masks.

For explicit Pauli inputs, applying signed permutations and propagating
derivatives backward are standard classical operations. Reverse
differentiation has a much broader established theory; see
[Baur–Strassen, *The complexity of partial derivatives*](https://www.sciencedirect.com/science/article/pii/030439758390110X),
*Theoretical Computer Science* **22**(3), 317–330 (1983).
The local audit provides the division-free Hopf recurrence, certified
dyadic rounding, duplicate-term accounting, and a term-sampling comparator.
These baselines limit the interpretation of the quantum T-count result;
they are not a claim to have invented reverse differentiation or classical
importance sampling. Unknown controlled observables retain their distinct
access model.
