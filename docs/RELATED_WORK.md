# Related work and contribution boundary

[← Source map](SOURCE_MAP.md) · [Complete narrative](../REVIEW.md) · [Landing page →](../README.md)

The selected paper studies one prescribed Hopf differential-frame operator
in exact arbitrary-one-qubit+CNOT and coherent Clifford+T models. Its
contribution is the complete-operator construction and the resulting
precision/workspace guarantees. The ingredients retain their established
attribution; no priority claim follows from a bounded literature search.

This page summarizes the lineage of Results A–D. The
[source catalogue](reference/SOURCE_CATALOGUE.md#extended-lineage-for-the-core-and-its-supporting-constructions)
retains the detailed comparisons, precise locators, and bibliography.
[Research comparisons](../research/RELATED_WORK.md) and the
[state-based QBP supplement](../supplements/state_based_qbp/RELATED_WORK.md)
have separate reading routes and claim boundaries.

## 1. Exact state preparation and the target frontier

P. Yuan and S. Zhang, *Quantum* **7**, 956 (2023), give the active exact
compiler toolkit and all-workspace state-preparation frontier:

```math
S=\Theta(2^n),\qquad
D=\Theta\!\left(n+\frac{2^n}{n+m}\right).
```

The published `arXiv:2202.11302v2` statements were also checked in v3.
Their Theorem 2, Lemma 5, Lemma 6, and Lemma 9 supply the comparison
frontier, ancilla-free multi-controlled X, complete UCG synthesis, and
coherent copy–uncopy. X. Sun, G. Tian, S. Yang, P. Yuan, and S. Zhang's
earlier time–space result is the historical predecessor.

These state-preparation tools can be adapted after the Hopf completion
is factored into addressed complete-operator layers. Preserving designated
marker columns is the additional obligation; preparing the correct state
alone does not meet it.

## 2. Uniformly controlled gates and the Möttönen route

Möttönen and coauthors developed uniformly controlled rotations, and
Bergholm and coauthors the general UCG language. The inherited one-clean
Hopf compiler computes a common zero-suffix flag and applies a prefix/flag
multiplexor. At zero workspace, padding one full-width UCG with identity
blocks leaves a generally dense Walsh/Gray-code angle table. The local
borrowed-suffix factorization instead retains two narrow half-angle UCGs
per addressed depth.

## 3. Borrowed workspace and controlled-unitary roots

Barenco and coauthors' controlled-unitary roots, later borrowed-wire
constructions, and conditional-clean/toggle-detection methods supply the
established cancellation context. Khattar–Gidney organize conditional work
and dirty-helper predicates. The local strict-zero echo uses
$C_p^2=R_y(\theta_p)$ and $XC_pX=C_p^{-1}$.

### Narrow strict-zero contribution

One original logical suffix bit carries a predicate temporarily. The target
operation depends on its original value, and the complete word restores it
on arbitrary entangled inputs. Two total-width-$`(d+2)`$ UCGs and linear
predicate toggles realize one addressed depth. Summing this construction
proves the strict-zero complete-frame frontier; borrowed wires, UCGs, and
conjugation identities themselves are inherited.

## 4. Sparse and restricted multiplexors

Xu and coauthors study restricted UCG structures and savings from limited
control participation. The Hopf specialization has one long zero-suffix
predicate shared by all prefix-selected rotations. Factoring that predicate
reduces physical UCG width; logical sparsity alone does not prune the
full-width rotation table. The [catalogue comparison](reference/SOURCE_CATALOGUE.md#4-sparse-and-restricted-multiplexors)
retains the primary reference and its distinct scope.

## 5. Coherent routing and workspace parallelism

The Hopf tail is an exact direct sum of subtree frames. Standard CNOT
fanout, Fredkin routing, controlled operations, and uncomputation realize
that sum. The local contributions are the tree-cut identity, register
layout, route–operate–unroute schedule, and common workspace envelope
for prefix and tail. The [explicit router](../compiler_robust_hopf/router.py)
is checked on complex prefix–suffix-entangled inputs.

## 6. Hopf coordinates and inverse-frame gradients

The first Hopf work supplies the balanced chart, inverse map, metric, and
coordinate directions. Hopf-QBP supplies the addressed state-and-marker
frame, shared magnitude record, separate phase stream, and checkpoint
interface. The compiler preserves these prescribed columns:

```math
W_{\mathbb R}|0^n\rangle=|\psi\rangle,\qquad
W_{\mathbb R}|\lambda(j)\rangle=|e_j\rangle.
```

The oriented incoming amplitude obeys
$\partial_{\theta_j}|\psi\rangle=a_j|e_j\rangle$ and $g_{j,j}=a_j^2$.
Canonical amplitudes are nonnegative; singular raw derivatives vanish
while marker continuations remain selected by the full parameter tuple.
The complex magnitude frame is $D_{\mathrm{ph}}W_{\mathbb R}$, with a
separate leaf-phase gradient stream.

## 7. Comparison by mathematical role

| Result | Established ingredients | Additional construction here |
|---|---|---|
| A: exact frame | State-preparation UCGs, borrowed-wire cancellation, reversible routing, Hopf geometry | Prescribed columns, strict-zero echo, tree cut and routed all-workspace schedule |
| B: sufficient-clean T count | Geometric bit sources, dirty lookup, compression, LCU/amplification, counting and diagonal lower bounds | Complete residual representation and shared-source composition, uniformly in precision |
| C: one-clean compilation | Full-space loaders, Clifford overlap algebra, dirty queries, standard amplification | Conjugated scalar-source word, native head-and-two-tail source, grouped frame and literal-phase corollaries |
| D: count and T-depth | Returned dirty indicators, controlled Clifford circuits, phase synthesis, arithmetic, phase-reference reuse | Charged unary-source reuse, complete-input query words, rectangular width allocation and uniform frame composition |

## 8. Contribution at theorem level

The [publication scope](../manuscript/PUBLICATION_SCOPE.md) gives the four
principal results and their full reservations. A matches exact size and
depth at every clean budget. B matches the precision-uniform T-count
frontier under sufficient clean work. C reduces the clean allocation to
one flag with its stated grouped upper bound. D gives simultaneous count
and T-depth bounds for real frames and matches both on explicit ranges.

The [frame-safety proof](FRAME_SAFE_COMPILATION.md) and two-qubit
counterexample explain why the completion matters. The
[QBP consequence](QBP_CONSEQUENCE.md) transfers the compiler result under
matched program, observable-access, and raw-gradient output conventions.
It neither proves a general gradient-query optimum nor differentiates a
discretely synthesized gate word.

## 9. Fault-tolerant lineage and the precision register

Gosset–Kothari–Wu (GKW), Theorems 1.1–1.2 and 4.1–4.2, supply optimal
unrestricted state/diagonal T-count benchmarks and ancilla-independent
lower bounds. Appendix B treats complete one-qubit multiplexors and
geometric error allocation. The local frame reduction uses the diagonal
subfamily; state preparation by itself leaves the completion unspecified.

Low–Kliuchnikov–Schaeffer (LKS), Section 2 and Appendix C, supply exact
dirty-bank lookup, reversible selection, and dirty indicators; Section 5
supplies finite-width circuit counting. Their bank identity extends to
reference-entangled inputs by linearity. The dirty banks do not turn an
arbitrary dirty output word into an initialized instruction register.

Bausch's Eqs. (4), (6), and Section 2.3.3 supply geometric precision weights,
an address-and-digit oracle, and a capped source. Its exact preparation
uses unary temporary work before binary conversion. The compact local
Gray construction proves the simultaneous linear T-count/logarithmic
peak-clean-width guarantee with ordinary binary labels and returned work.
Neither geometric sampling nor general dirty lookup is claimed as new.

## 10. What the full-frame T theorem adds

The [sufficient-clean construction](FAULT_TOLERANT_COMPILER.md) retains a
shared source across complete residual corrections and applies one final
amplification. Failure-history compression is inherited from Low–Wiebe,
Lemma 13, in the form explained by Fang–Lin–Tong, Lemma 3 and Appendix D.
LCU and normalization-two amplification follow Berry and coauthors.
The new specialization is the residual representation and charged frame
schedule, with complete initialized-isometry error.

This gives $\Theta(\sqrt{NL}+L+NL/(n+a+b))$ T count under its sufficient
clean reservation, with exact dirty return and $O(NL)$ Clifford cost.
The simpler direct sampler already matches this in some polynomial-
precision regimes. Shared-source composition removes its repeated $nL$
precision charge uniformly; those easier regimes are not a separate gain
attributed to the composition.

## 11. Operator sources and the one-clean refinement

Kerenidis–Prakash's full-space anticommuting loaders, the overlap algebra
also used by Chee and coauthors, and Bravyi's paired-Majorana representation
supply the source algebra. LKS supplies dirty lookup; Berry and coauthors
supply the oblivious amplification setting. The
[two-clean baseline](OPERATOR_SOURCE_COMPILER.md) gives $O(N+nL)$ T count
at $b\ge L+n+7$ by programming full-space operator identities on a dirty
core. Lookup banks and selectors return exactly; core return is approximate
within the joint isometry norm.

### One initialized flag

The [one-clean word](ONE_CLEAN_COMPILER.md) conjugates scalar-source blocks
and routes rejected components through different logical Pauli operators.
Its head-and-two-tail source packs the scalar coefficients into the stated
dirty reservation. Five-call amplification at $\sin(\pi/10)$ specializes
the standard Brassard–Høyer–Mosca–Tapp success law, with a local proof of
literal phase, full-output error, and core return.

[Conditional-suffix grouping](CONDITIONAL_SUFFIX_COMPILER.md) uses logical
suffix work only in its selected sector and streams coarse programs. With
the one-clean refinement it gives $`O(N+L\ell_*(n))`$ T gates at
$b\ge L+n+7$, where $`\ell_*(n)=1+\log_2^*(n+2)`$.
Its final conditional-work return error is included in the norm.
The high-precision full-frame gap remains $\Omega(N)$ to
$`O(N\ell_*(n))`$.

## 12. Contemporary comparisons and the broader compiler contribution

The retained comparisons distinguish initialization, literal phase,
complete-input error, and precision uniformity. Tan and Yuan–Zhang–Zi treat
generic unitary synthesis with initialized ancillary allocations; the
complete initialized-isometry criterion already occurs there. GKW's
instruction-based diagonal construction and Yamazaki–Akibue's precision
improvements retain initialized precision-length instruction work in their
displayed dirty-assisted constructions. General uncomputation/compression
results keep their original query or dilation hypotheses. These are
different contracts, not missing attribution.

The one-clean construction covers every literal diagonal and complete
one-target U(2) multiplexor, with certified approximate Euler coordinates,
its explicit dirty reservations, and worst-case matching banked T counts.
It preserves literal block phases and joint core-return error. It asserts
neither better leading constants nor efficient general classical coordinate
search. [Detailed comparisons and exact thresholds](reference/SOURCE_CATALOGUE.md#12-contemporary-comparisons-and-the-broader-compiler-contribution)
remain in the catalogue.

For D, Khattar–Gidney's two-dirty logarithmic-depth MCX provides predicates;
LKS-based dirty echoes supply lookup. Established controlled-Clifford,
commuting-selection, phase-synthesis, arithmetic, and phase-reference ideas
are credited individually in F29–F37. The local result composes charged
native source boundaries and complete-input queries with a rectangular
workspace allocation. Its error and resource constants are uniform in
precision. It concerns real frames and T-depth with Clifford interlayers;
standalone source obstructions and radial filters supply no additional
premise for the current theorem. The [dated development record](../research/RELATED_WORK.md)
retains the narrower intermediate comparisons and later optional studies.

## References highlighted here

- Exact toolkit: P. Yuan and S. Zhang; historical predecessor: X. Sun and coauthors.
- Controlled operations and UCGs: Barenco; Möttönen; Bergholm; Khattar–Gidney.
- Fault-tolerant count and lookup: Bausch; LKS; GKW.
- Coherent compression and amplification: Low–Wiebe; Fang–Lin–Tong; Berry and coauthors; Brassard and coauthors.
- Full bibliographic entries and precise imported propositions: [source catalogue](reference/SOURCE_CATALOGUE.md#references-highlighted-here).
