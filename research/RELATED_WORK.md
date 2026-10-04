# Literature context for research studies

[Research questions](../docs/OPEN_PROBLEM.md) · [Core related work](../docs/RELATED_WORK.md) · [Source catalogue](../docs/reference/SOURCE_CATALOGUE.md)

This record preserves the dated literature comparisons that accompanied
the depth and endpoint investigations. These sections describe the scope
of individual constructions and obstructions as they were developed;
they are not a second statement of the selected Results A–D. Several
constructions became supporting proof dependencies, while other source,
filter, scheduling, and cache studies remain optional research results.
The [core source map](../docs/SOURCE_MAP.md) identifies the current dependencies.
Restricted source lower bounds do not become unrestricted frame lower bounds.

## Reading index

| Topic | Retained sections |
|---|---|
| Source depth, filtering, and literal phase | [16](#16-precision-depth-and-workspace-assumptions-2-october-2026), [17](#17-fixed-point-filtering-and-literal-phase-2-october-2026) |
| Grouped programs and the depth constructions supporting D | [18](#18-grouped-programs-and-chunked-dirty-queries-2-october-2026), [19](#19-unary-phase-source-reuse-and-linear-t-depth-3-october-2026), [20](#20-blocked-bilinear-queries-and-the-dirty-width-tradeoff-3-october-2026), [21](#21-uniform-precision-and-rectangular-query-allocation-3-october-2026) |
| Limits of external depth comparisons | [22](#22-scope-of-recent-depth-lower-bounds-3-october-2026) |
| Nonuniform chunks, protected sources, and predicate caches | [23](#23-nonuniform-chunks-in-a-returned-dirty-indicator-3-october-2026), [24](#24-a-protected-source-across-early-groups-3-october-2026), [25](#25-cached-activity-predicates-across-group-windows-3-october-2026), [26](#26-exact-shared-prefix-cache-capacity-3-october-2026) |
| Retained-source factors and invariant sectors | [27](#27-retained-source-factors-and-invariant-completion-sectors-3-october-2026) |

## 16. Precision depth and workspace assumptions (2 October 2026)

The [source-depth audit](depth/SOURCE_T_DEPTH.md) separates a restriction of the
geometric source implementation from the unrestricted
[Hopf T-depth problem](../docs/T_DEPTH_COMPILER.md#4-lower-bounds-and-the-remaining-depth-gap).
The relevant denominator technique is already present in
[Casas et al., *Matchgate synthesis via Clifford matchgates and T gates*,
arXiv:2602.05425v1, Section III.2.2, Eq. (29)](https://arxiv.org/html/2602.05425v1#S3.SS2.SSS2).
In their Majorana representation, Clifford matchgates permute signed modes,
and one layer of disjoint non-Clifford mode rotations increases the least
square-root-of-two denominator exponent by at most one. Their resulting
depth lower bound concerns exact synthesis within this restricted gate
set. Applying that mechanism to the geometric source supplies a local
obstruction to parallelizing that source in the same representation; it
does not introduce a new general lower-bound method. Arbitrary Clifford
interlayers need not preserve the Majorana span, and the frame compiler
need not implement this source exactly, or use it at all.

There is also an established alternative when precision-sized **clean**
workspace is available. [Vasconcelos, *Depth-Optimal Quantum Compilation*,
arXiv:2609.34659v1, Theorem 8](https://arxiv.org/html/2609.34659v1#S3.SS4.SSS1)
constructs a fully unitary single-qubit rotation approximation with
$`O(L)`$ initialized ancillas, $`O(L)`$ gates, and $`O(\log L)`$
elementary depth, for $`L=\Theta(\log(1/\varepsilon))`$. Its error
includes ancilla return and a stated common phase. Fixed exact Toffoli
decompositions preserve these orders over Clifford+T. This demonstrates
that linear precision depth is not intrinsic to rotation approximation.
It does not give the same construction with two clean qubits and arbitrary
dirty work. Its bounded-arity elementary-depth lower bound also does not
bound T-depth when unrestricted Clifford circuits between T layers are free.

[Kim, *Catalytic z-rotations in constant T-depth*, arXiv:2506.15147v3,
Section 3](https://arxiv.org/pdf/2506.15147), published in *Quantum*
**10**, 2191 (2026), explicitly leaves constant-T-depth rotation using
only clean or dirty ancillas open. The depth-three construction assumes
a prepared nonstabilizer catalyst; the supplied-catalyst resource cannot
be replaced by arbitrary borrowed qubits. Charging catalyst preparation
and its initialized workspace is necessary before composition with the
present compiler. Its final note records subsequent depth-two and
measurement-assisted depth-one refinements. The universal-catalyst
comparison in Section 19 below includes Kim–Laakkonen's later construction.

These comparisons identify usable techniques and their workspace
conditions. The 2 October source-restriction pass yielded no asymptotic
improvement to the unrestricted complete-frame depth bound, and no
matching lower bound. Section 19 records the later unary-source compiler.
The exact source restriction must therefore remain separate from both
optimal approximate Hopf T-depth and the constant-clean T-count endpoint.

The [two-layer obstruction](depth/SHALLOW_SOURCE_OBSTRUCTION.md) uses a
different, full-input argument. [Aaronson–Gottesman,
arXiv:quant-ph/0406196v5, Section III](https://arxiv.org/pdf/quant-ph/0406196v5)
gives the discrete magnitudes of stabilizer overlaps.
[Zhang–Zhang, arXiv:2409.13809v2, Theorem III.1,
Eqs. (10)–(11)](https://arxiv.org/html/2409.13809v2#S3.SS1)
explicitly uses the fact that a T layer conjugates Paulis to Hermitian
Cliffords. Splitting a two-layer transfer coefficient at the middle
Clifford reduces it to a normalized Clifford trace. The local proof
then derives a transfer alphabet and constant approximation gaps for
the geometric sources with any dirty width. These standard ingredients
are attributed; no growing unrestricted depth lower bound is inferred.

The [conditional geometric source](../docs/CONDITIONAL_GEOMETRIC_SOURCE.md)
supplies a different positive interface. It prepares the geometric
one-hot state using reversible prefix ORs and native controlled
Hadamards, clears noninvariant tree scratch before that batch, and
uses conditional suffix workspace rather than an externally initialized
precision register. Its preparation, selective reflection, actual
inverse, inactive-sector identity, and return error are charged locally.
This is an explicit construction from standard reversible and
block-encoding ingredients; it does not import a catalyst or claim
generic fast rotation synthesis with arbitrary dirty helpers. The
complete-frame query and suffix-predicate bottlenecks require a separate
composition, supplied in Section 18 below.

## 17. Fixed-point filtering and literal phase (2 October 2026)

[Grover, *Fixed-point quantum search*, PRL **95**, 150501 (2005)](https://arxiv.org/abs/quant-ph/0503205),
Eq. (1) and Section 3, gives the three-call pi-over-three phase sequence
that cubes the failure probability. The
[radial-filter chapter](depth/HOPF_RADIAL_FILTER.md) applies this established
sequence to the actual amplified Hopf source. It retains the complex
accepted scalar: cubic suppression of rejected amplitude leaves a
quadratic full operator error because an accepted phase remains.

The selective phase is not a free exact Clifford+T gate. Its
determinant-one representative factors into three one-qubit Pauli-Z
rotations on the existing two flags, and a literal Clifford global
correction calibrates the complete word. The standard phase-sensitive
one-qubit approximation theorem, [GKW Lemma 2.3](https://arxiv.org/html/2411.04790v3#S2),
supplies the six ancilla-free native phase words at their charged
precision. The source-width cap can then use a half-logarithmic
coefficient, while the added calls and phase words leave the established
asymptotic T-count and T-depth unchanged. The fine-word search distinction
discussed above still applies.

The [short-echo audit](depth/HOPF_FLAG_ECHO.md) is separate: its exact error
formulas rule out uniform radial cancellation for four specific Pauli
interpositions, and include a full-space equal-mask exception. They are
not a no-go theorem for fixed-point amplification or general composites.

## 18. Grouped programs and chunked dirty queries (2 October 2026)

The [grouped-program theorem](../docs/GROUPED_PROGRAM_PREFETCH.md#8-complete-frame-theorem-at-fixed-accuracy)
combines the conditional source with exact prefetch into an active zero
suffix. The program stays read-only through source leakage, allowing its
actual inverse to erase it exactly. Private predicate trees and literal
enable-dependent amplification phases complete the local interface.
The [chunked indicator](../docs/CHUNKED_DIRTY_INDICATOR.md) interpolates between
the existing routed and read-only counter constructions. A dirty linear
tree conjugates a root flip into a selected-path translation, and the
outer echo removes every unknown dirty offset.

Conditional workspace, reversible conjunctions, linear conjugation,
bilinear lookup, and block-encoding amplification are established tools.
The local contribution is their charged complete-frame composition:
an adaptive early prefetch schedule and a distinct late chunk allocation
give fixed-accuracy $`O(n\log\log(n+2))`$ T-depth with optimal-order
$`O(\sqrt N)`$ T-count, $`O(N)`$ Clifford count, two external clean
flags, and sufficient $`C_\eta\sqrt N`$ dirty workspace. No generic
priority claim is inferred. The variable-precision theorem keeps its
previous bounds. This geometric-source schedule remains a valid
predecessor to the linear fixed-accuracy bound below; unrestricted
large-width depth optimality and the constant-clean high-precision
endpoint remain open.

## 19. Unary phase-source reuse and linear T-depth (3 October 2026)

The [unary phase-source compiler](../docs/UNARY_PHASE_GRADIENT.md) uses established
phase kickback. [Jones et al., arXiv:1204.0567, Section 2.1,
Eqs. (2)–(4), and Section 4.1, Fig. 15](https://arxiv.org/pdf/1204.0567)
describe Fourier eigenstates of modular shifts, programmable phases,
and repeated reference reuse. The local source uses a unary encoding,
coherently selects its cyclic shift from a loaded one-hot program, and
prepares/unprepares it inside the active zero suffix. Neither phase
kickback nor shared phase-reference preparation is claimed as new.

The cyclic convolution uses the standard three-product Karatsuba
identity; see [Iggy van Hoof, arXiv:1910.02849v2, Section 4.2](https://arxiv.org/html/1910.02849v2).
Recursive scalar products give the bilinear rank bound, and reduction
modulo $`X^q-1`$ is linear. The constant T-depth does not come from
van Hoof's space-efficient reversible multiplication schedule. It comes
from the local guarded trilinear-phase construction with private
conditional work and the retained native Toffoli word. Parallel
parity-phase synthesis with initialized ancillas is already explicit in
[Selinger, arXiv:1210.0974v2, Section 2, Eqs. (5)–(6), and
Theorem 4.1](https://arxiv.org/pdf/1210.0974). The local proof must
additionally establish identity on arbitrary inactive inputs and exact
temporary return on every active source state.

[Kim–Laakkonen, arXiv:2512.24982v1, Theorems 3, 5 and 6,
and Section 5.1](https://arxiv.org/html/2512.24982v1) already give
constant-depth controlled CNOT/Clifford circuits and a universal
logarithmic-size catalyst for rotations. Their catalytic rotation
chooses an angle-dependent CNOT matrix classically; its depth-one
implementation uses measurement-assisted uncomputation. Their charged
preparation uses measured phase estimation, expected repetitions, and
dynamic compilation after selecting the catalyst eigenvalue. This does
not directly supply the unary compiler's coherently loaded angle table,
unitary source boundary pair, or conditional-work return contract.
The present bound does not settle Kim's generic clean/dirty-only
constant-T-depth rotation question: source preparation is charged and
the complete frame has linear, rather than constant, T-depth.

A recent comparison is [Wu et al., *Shared Phase Arithmetic for Parallel
Quantum Rotations*, arXiv:2609.36574v1, Sections II.C, III.B and
IV.C–E](https://arxiv.org/html/2609.36574v1), submitted 29 September 2026
and checked here on 3 October. Their Theorem 1 expresses reversible
phase-function evaluation, one binary addition, and decoding; Proposition 1
gives Clifford-only encoding for disjoint binary supports. The displayed
ripple adder permits measurement/feedforward and has depth linear in
the phase-register width. Preparation is charged separately, and parallel
batches require separate resources. These are useful shared-arithmetic
precedents, not the constant-depth unary-program interface used here.

The additional result is the complete charged composition: one unitary
source boundary pair per group, coherent one-hot shift selection,
inactive-sector identity, the Hopf angle-stability bound, and an
early/late workspace and query allocation. For every fixed accuracy
$`\eta`$, it gives one complete real-frame circuit with

```math
D_T=O_\eta(n),\qquad T=O_\eta(\sqrt N),\qquad G=O_\eta(N),
```

two external clean flags, and sufficient $`C_\eta\sqrt N`$ dirty
workspace. The initialized-isometry error includes all returned work
and arbitrary dirty-reference entanglement. It uses no supplied phase
state or intermediate measurement. This improves the fixed-accuracy
depth upper bound; it proves neither an unrestricted matching depth
lower bound nor the constant-clean high-precision endpoint. The sources
above identify inherited ingredients and interface distinctions, not
priority for the composite construction. This square-root-width result
is also a corollary of the broader fixed-accuracy tradeoff below.

## 20. Blocked bilinear queries and the dirty-width tradeoff (3 October 2026)

The [blocked bilinear lookup](../docs/BLOCKED_BILINEAR_LOOKUP.md) combines the
existing [two-pass dirty traversal](../docs/AMORTIZED_DIRTY_LOOKUP.md#2-selecting-a-shear-with-dirty-unary-traversal)
with the [bilinear indicator echo](../docs/PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction).
Their lineage remains [LKS, Appendix C](https://arxiv.org/html/1812.00954v2),
[Khattar–Gidney, Sections 4 and 7](https://arxiv.org/html/2407.17966v2),
and the controlled-linear and commuting-basis constructions of
Kim–Laakkonen and Boyd discussed in Section 9. Standard Boolean phase
polarization supplies the controlled leaf: for a bilinear phase P,
apply $`(-1)^{dP}`$, toggle the dirty bit d by hz, apply the phase
again, and undo the toggle. This full four-step word leaves exactly
$`(-1)^{hzP}`$ and returns d.
Neither that algebra nor generic dirty cancellation is a novelty claim.

The local proof establishes a rank-sensitive native leaf with preserved
controls and exact helper return, a selected-block traversal, and a chunk
allocation whose depth sums across the final Hopf layers. Together with
the unary early groups, it gives, for each fixed $`0\lt\eta\le1/64`$,
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$ and
$`b\ge17(L+n+7)`$, one complete real-frame circuit with

```math
T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad G=O_\eta(N),
\qquad D_T=O_\eta\!\left(\frac N{b^2}+n\right).
```

Two external clean flags suffice. Full initialized-isometry error,
literal phase and arbitrary dirty-reference return retain their existing
contracts. Sources, actual inverses and all queries are charged.

When the interval is nonempty, the inherited count lower bound and the
physical width $`n+2+b=\Theta(b)`$ give simultaneous worst-case matches:

```math
17(L+n+7)\le b\le\sqrt{N/n},\qquad
T^\star=\Theta_\eta(N/b),\quad D_T^\star=\Theta_\eta(N/b^2).
```

This extends the fixed-accuracy matching window; it introduces no new
external synthesis premise or depth lower-bound method. The preceding
variable-accuracy bounds remain valid; Section 21 gives a uniform
precision refinement. The unary
square-root-width schedule remains a valid predecessor and corollary;
unrestricted large-width depth optimality and the constant-clean
high-precision endpoint remain open. No generic lookup priority is claimed.

## 21. Uniform precision and rectangular query allocation (3 October 2026)

The [uniform precision theorem](../docs/UNIFORM_PRECISION_DEPTH.md) uses the same
native selected-block query, dirty traversal and bilinear echo as Section
20. Rectangular block dimensions balance word precision against indicator
cost. The proof makes the unary-source cutoff uniform at low precision
and uses the existing hybrid when its source term absorbs the logarithmic
depth contribution. These are allocation and composition results; the
external ingredients in Sections 19–20 are unchanged, and no new native
query or synthesis premise is assumed.

With $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$, two clean flags
and $`b\ge17(L+n+7)`$, one complete real-frame circuit has absolute,
precision-independent constants in

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+nL\right).
```

For $`6\le L\le\log_2(n+2)/16`$, the sharper bounds are

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}\right),\qquad G=O(NL),
\qquad D_T=O\!\left(\frac{NL}{b^2}+n\right).
```

The existing count lower bounds divided by physical width give
simultaneous worst-case $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$ through $`b\le\sqrt{N/n}`$
generally, and through $`b\le\sqrt{NL/n}`$ in the low-precision
regime, above the literal threshold and when the intervals are nonempty.
The low-precision count is optimal at every eligible width; the general
count retains nL and is asserted optimal outside its matching interval
only under a sufficient condition such as $`L\le N/n^2`$.
All preparation, queries, actual inverses and work return remain charged.
The selected high-precision endpoint and unrestricted large-width depth
optimality remain open. No priority claim follows from this composition.

## 22. Scope of recent depth lower bounds (3 October 2026)

These comparisons guide further research; neither is an imported compiler
premise. [Parham, arXiv:2504.19966v1](https://arxiv.org/html/2504.19966v1),
Proposition 1.8, relates T-depth to alternations of unrestricted Clifford
and shallow circuits. Theorems 1.14–1.15 connect sufficiently strong
explicit-state and Boolean-function lower bounds to classical threshold
circuit lower bounds, with polynomial clean workspace in the model.
Section 6 explicitly says that no analogous reduction is known for
general prescribed-unitary implementation. This is therefore not a
blanket complexity barrier to complete-frame operator lower bounds.

[Al-Ghattas–Gamarnik–Kiani, arXiv:2610.02166v1](https://arxiv.org/html/2610.02166v1),
submitted 1 October 2026, treats arbitrary Clifford blocks. Corollary
1.7(iii) and Section 4.2 cover every fixed number of shallow blocks at
total width $`M=O(n)`$; arbitrary-width extensions concern specific
one-round classes. Lemma 4.4 retains an $`O(kM^2)`$ entropy term, so
this theorem does not cover the current $`b\asymp\sqrt N`$ allocation.
No applicable growing depth bound was identified in this comparison;
that audit outcome does not prove that such a bound is impossible.

## 23. Nonuniform chunks in a returned dirty indicator (3 October 2026)

The [nonuniform indicator](depth/NONUNIFORM_DIRTY_INDICATOR.md) refines the
existing [chunk-boundary tree echo](../docs/CHUNKED_DIRTY_INDICATOR.md) using the
same [read-only dirty conjunction](../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences).
Its lineage remains [LKS, Appendix C, Theorems 1–2](https://arxiv.org/html/1812.00954v2)
for dirty indicators and count-efficient XOR queries, and
[Khattar–Gidney, Sections 4, 5.5 and 7](https://arxiv.org/html/2407.17966v2)
for toggle detection, conjunctions and unary tree traversal. LKS's
displayed s-bit indicator has $`O(2^s)`$ T-count and $`O(s^2)`$
depth without extra ancillas; its depth includes Clifford operations.
Parallel phase synthesis retains the Selinger attribution in Section 19.

The local refinement consumes large chunks while few prefixes are live,
then reduces the remaining address length to
$`\lceil4\log_2(r+2)\rceil`$ at each step, finishing at bounded
length. The weighted tree ledger gives

```math
T,G,w=O(2^s),\qquad D_T=O(\log_2(s+2)).
```

The indicator is exact on all inputs, including arbitrary outputs and
dirty-reference correlations. This is a T-depth bound with charged
Clifford count, not a bound on total elementary depth. For $`Q=2^r`$
rows and word width m, the existing rectangular bilinear construction,
under its stated additional dirty-width reservation B, consequently has

```math
T=O\!\left(\sqrt{Qm}+m+\frac{Qm}{B}\right),\qquad G=O(Qm),
\qquad D_T=O\!\left(m\left[1+\frac{Q}{B^2}\right]+\log_2(r+2)\right).
```

In the existing [low-precision unary regime](../docs/UNIFORM_PRECISION_DEPTH.md) at sufficient
$`b=\Theta(\sqrt{NL})`$, the late tail becomes $`o(n)`$.
The early groups retain their $`O(n)`$ contribution, so this does not
change the proved complete-frame depth order. The added result is the
nonuniform allocation and charged query composition; no generic
priority or new external synthesis premise is claimed.

## 24. A protected source across early groups (3 October 2026)

The [protected-source refinement](depth/PROTECTED_UNARY_SOURCE.md) extends the
charged unary reference from one group to the entire early segment.
Reusable phase references remain the Jones et al. precedent in Section
19; the Kim and Kim–Laakkonen catalyst comparisons retain their different
preparation and workspace contracts. No new external premise is needed.

A fixed terminal logical bank B stays outside every early target set.
The initially zero flag H records $`[B=0]`$ before the unconditional
native preparation U. The second initialized flag h records
$`H[Z=0]`$ for each group's outer suffix Z with B removed. Each group
returns h and its temporary work exactly for arbitrary source-core states.
On initial nonzero-B inputs, the entire middle word is identity, so
U cancels with its actual inverse. On the initial zero sector, the ideal
Fourier source survives every group. One two-boundary comparison therefore
charges $`2\delta`$ for the entire early segment, including final H
erasure. Final source-bank and H return are approximate within that bound;
the inactive cancellation requires the stated initialized flags.

With a disjoint bank and the revised sufficient cutoff, preparation,
unpreparation and the initial-zero predicate pair contribute
$`O(L+\log(n+2))`$ T-depth, uniformly in the stated low-precision
regime. Source preparation is fully charged and no source is supplied.
The target stages/selectors, group predicates, and program query pairs
retain separate $`O(n)`$ early contributions. Thus this is a source-use
and error-accounting refinement, with no new complete-frame depth order,
matching interval, high-precision endpoint, or generic reuse-priority claim.

## 25. Cached activity predicates across group windows (3 October 2026)

The [windowed-predicate schedule](depth/WINDOWED_GROUP_PREDICATES.md) combines
the existing borrowed multi-control X construction attributed to
Khattar–Gidney with the
[consume-before-change selector schedule](../docs/GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables)
and the protected source of Section 24. A window shares one predicate
for its unchanged outer suffix. Short block predicates and a suffix-product
chain supply each group's activity; each cached predicate is erased while
its logical controls still have their original values.

For $`J=\lceil\log_2(n+2)\rceil`$ groups per full window, the
protected logical bank gains $`2J+1`$ cache bits. They start at zero
on H equal to one and may be arbitrary on H equal to zero. The exact
schedule restores the cache in both sectors, with h initially zero, and reuses the
same two returned dirty predicate helpers. No additional external clean
flag or dirty-helper reservation is introduced. The protected source's
global $`2\delta`$ error bound, including final bank and H return,
is unchanged.

In $`6\le L\le\log_2(n+2)/16`$, the total early activity depth is

```math
O\!\left(\frac{n\log\log(n+2)}{\log(n+2)}+\log(n+2)\right).
```

All cache preparation, consumed cleanup, and native flag toggles are
charged. Logical stages/selectors and program query pairs retain separate
$`O(n)`$ early allowances. This is a predicate-allocation refinement,
with no new external premise, generic priority claim, complete-frame
depth order, matching interval, or high-precision endpoint result.

## 26. Exact shared-prefix cache capacity (3 October 2026)

The [shared-prefix audit](depth/SHARED_PREFIX_QUERY_AUDIT.md) applies
[Nielsen–Chuang's programming theorem](https://arxiv.org/pdf/quant-ph/9703032),
pp. 1–2, Eq. (3) and the Result/Eqs. (6)–(10). Encoding preserves read-only
prefix x and acts reversibly on b arbitrary dirty bits and c initialized
cache bits, giving each program support rank $`2^b`$. A fixed exact
processor without access to x, implementing R distinct XOR-query
signatures, requires pairwise orthogonal supports. Consequently
$`R2^b\le2^{b+c}`$, or $`R\le2^c`$. All other prefix-bearing
bits must be counted in the program; conditional logical zeros count
among c. Phase-query signatures instead require distinctness modulo
global phase. This restricts the stated interface, without a T-depth
lower bound or priority claim.

The [charged continuation](depth/SHARED_PREFIX_QUERY_AUDIT.md#6-a-charged-two-group-program-refresh)
retains prefix access and so lies outside that capacity restriction. Its
exact two-group identity reuses a program bank and charges an intervening
table-difference query. Fixing one legal row to a constant reduces an
ordinary fresh query to that correction by a Clifford output XOR. A cheap
affine-factor example uses the existing dirty echo and literal Toffoli
words; low rank without cheap factors does not suffice. These are local
algebraic consequences of the existing query interface, with no new
external synthesis premise or unrestricted depth lower bound.

## 27. Retained-source factors and invariant completion sectors (3 October 2026)

The [retained-source continuation](depth/RETAINED_SOURCE_FUSION.md#9-retained-source-fusion-without-a-larger-modulus)
uses the already attributed Vatan–Williams magic basis (F21), with
Laurent entries on one common cyclic source. Opposite determinant
monomials keep the original source modulus; they are retained operator
phases, not discarded scalars. The Bell identity separates a factor that
fixes the shared boundary from its complementary transport.

The [completion schedule](depth/RETAINED_SOURCE_FUSION.md#10-two-source-shifts-for-the-stabilizer-completion)
then uses a property of the Hopf word: these factors have mutually
orthogonal last-nonzero-pair sectors. Cached sector labels survive both
paired-axis rounds, so the existing guarded cyclic-shift primitive serves
the whole completion twice. Parallel conjunctions, conditional work,
phase kickback, and bilinear multiplication retain their existing source
attribution. The local contribution is this exact factorization and
invariant-cache schedule, including full inactive cancellation and charged
derived program rows. No new external synthesis premise is imported.
The remaining transport still has a sequential height-linear schedule;
the result does not improve the complete-frame depth order or establish
an unrestricted lower bound.
