<a id="two-resource-questions-for-complete-hopf-frames"></a>
<a id="endpoint-theorem-and-remaining-depth-question"></a>
<a id="two-clean-count-theorems-and-remaining-depth-question"></a>
<a id="one-clean-count-frontier-and-remaining-depth-question"></a>
<a id="one-clean-count-frontier-and-sublinear-large-width-depth"></a>

# Resource frontiers and the remaining depth question

[Selected claims](../manuscript/PUBLICATION_SCOPE.md) · [Proof map](README.md) ·
[Research catalogue](../research/README.md) · [Continuing workspace](../WORKSPACE.md)

This page states the current resource frontier and its scope. Detailed proofs,
component constructions and historical outcomes have their own homes in the
research catalogue; they are not additional premises of the general theorem.

Write $`N=2^n`$, with integer $`n\ge1`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. The coherent Clifford+T
model permits arbitrary Clifford interlayers, a clean bank of size a and
arbitrary dirty work of size b. If $`J_a`$ appends clean zero inputs,
the complete-frame contract is

```math
\|VJ_a-J_a(W\otimes I_b)\|_{\rm op}\le\eta.
```

It includes all logical columns, dirty inputs entangled with any reference,
literal signs and phases, singular continuations, leakage and return.
T-depth is distinct from total depth; elementary Clifford count G is charged
separately. Supplied angles and leaf phases admit terminating certified
evaluation, without an efficient-preprocessing claim.

## Current status

| Question | Strongest established statement or remaining obligation |
|---|---|
| Selected count frontier | Results B–C are settled at their stated clean/dirty allocations; the exact-base and banked ranges below remain distinct |
| General depth upper | R96 at all precisions. R99–R102 and [R104](../research/depth/EXACT_KEY_COVER_PRECISION.md) improve the precision term at the same budgets. R104 gives sublinear additive depth through $`L=o(N)`$ in its high-precision range. |
| Precision-dominated depth | [R103](../research/depth/PRECISION_COUNT_DEPTH_TRADEOFF.md) matches R96 at $`L\ge N`$ with one clean input and $`T=O(L)`$. [R105](../research/depth/CLEAN_SPACE_PRECISION_DEPTH.md#structured-coarse-prefix) reaches polylogarithmic depth using nearly linear clean space, with the same count and dirty width. It comes within a logarithmic factor of R103 up to $`t=L/P(n)`$ for a fixed polynomial. |
| Selected Result D | R79 remains the selected all-precision count/depth theorem |
| Optimal large-width depth | At fixed accuracy and square-root dirty width, logarithmic general upper depth versus a five-band unrestricted lower bound remains open |
| Main constructive question | A sublogarithmic complete indicator or coherent signed-gap source at fixed accuracy. The [work checkpoint](../WORKSPACE.md) gives a weaker per-row target. |
| General lower question | A growing native T-depth bound on the prescribed Boolean-table frame family, allowing table-specific compilation, arbitrary Clifford interlayers and the admitted dirty width |

The [working checkpoint](../WORKSPACE.md#current-checkpoint-native-source-scope-and-complete-readout)
identifies the current question and first unpaid operation. The [route decisions](../research/ROUTE_DECISIONS.md)
record stopped mechanisms and the conditions that would justify reopening them.

<a id="precision-uniform-two-clean-count"></a>
<a id="precision-uniform-one-clean-count"></a>
<a id="optimal-banked-count"></a>

## Result C: the one-clean count frontier

For every integer $`n\ge1,L\ge6`$, the
[one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
has full error at most $`2^{-L}`$ and universal constants. Its base
dirty reservation is $`B=L+n+7`$ for complete real frames, or
$`B=L+n+8`$ for phase-dressed complex magnitude frames.

| Dirty allocation | T-count guarantee |
|---|---|
| Every $`b\ge B`$ | $`T=O(N+L)`$ |
| Exactly $`b=B`$, with $`L\ge\max\{6,n\}`$ | Worst-case $`T^\star=\Theta(N+L)`$ |
| Every $`b\ge2B`$ | Worst-case $`T^\star=\Theta(\sqrt{NL}+L+NL/b)`$ |

All upper constructions have $`G=O(NL)`$ and charge sources, queries,
actual inverses and work return. At base width with $`L\lt n`$, or
in the transition $`B\lt b\lt2B`$, retain the upper bound without
a general matching claim. These count results add no depth guarantee.
[Publication scope](../manuscript/PUBLICATION_SCOPE.md#result-c-the-one-clean-count-theorem)
is the selected claim ledger, including the distinct sufficient-clean
Result B and diagonal/multiplexor corollaries.

## Constant-clean high-precision count

At the separate literal allocation

```math
a=2,\qquad \eta=2^{-N},\qquad L=N,\qquad b=N+n+7,\qquad n\ge3,
```

both complete frame families have worst-case
$`2N-1\le T^\star\le O(N)`$, hence $`T^\star=\Theta(N)`$.
The [growing-source construction and literal phase composition](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
give the upper bound; the
[literal-phase root witness](FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness)
gives the lower bound even at unrestricted width. This two-clean allocation
saves one dirty qubit for the complex family relative to Result C.

## Large-width T-depth

**All-precision upper R96.** For every integer $`n\ge1,L\ge6`$,
the [all-precision theorem](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md)
uses one external clean qubit and retains the literal thresholds

| Complete target | Sufficient dirty allocation |
|---|---|
| Real Hopf frame | $`b\ge17(L+n+7)`$ |
| Prescribed phase-dressed complex magnitude frame | $`b\ge17(L+n+8)`$ |

It gives the following bounds on one native circuit with full error at most $`2^{-L}`$.

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\ell),\qquad \ell=\log_2(n+2).
```

The T-count is optimal at eligible allocations. At fixed accuracy and square-root dirty width, the general depth gap is constant lower versus logarithmic upper.
**Precision lower bound R103.** With k clean inputs and unrestricted dirty width, a uniform deepest rotation forces $`D_T T\ge (L-\log_2 25)^2/((k+1)\ln 4)`$ in the complete initialized-isometry norm. At $`L\ge N`$ and the R96 width thresholds, $`D_T=\Theta(L)`$ is worst-case optimal at the simultaneous one-clean $`T=O(L)`$, $`G=O(NL)`$ budgets. [Proof](../research/depth/PRECISION_COUNT_DEPTH_TRADEOFF.md)

[R99](../research/depth/PROGRAMMED_ROTATION_DEPTH.md) gives the same budgets
with depth $`O(NL/b^2+\ell+\log_2^2(L+2))`$ through $`L\le cn`$.
[R100](../research/depth/GEOMETRIC_CHUNK_DEPTH.md) retains them for
$`n^{1/4}\le L\le2^{c'n}`$ with $`a=\log_2(L+2)`$ and depth

```math
D_T=O\!\left(NL/b^2+\min\{(1+L/n)a^2,\ (1+L/n)a+L^{1/8}\ell\}\right).
```

At $`b=\Theta(\sqrt{NL})`$, R100 gives $`o(L)`$ for $`\log_2 L=o(n)`$. [R102](../research/depth/COMPACT_PROGRAM_PRECISION.md) gives $`O(NL/b^2+(1+L/n)\ell)`$ for $`6\le L\le2^{c''n}`$.
The bounds match worst-case depth where $`NL/b^2`$ dominates.

**Selected Result D: [R79](../research/depth/UNIFORM_SUBTREE_DEPTH.md).**
Use the same simultaneous guarantees and write the width threshold as $`17B_D`$. Define

```math
F(n,L)=\min\{nL,(L+\ell)\log_2(n+L+2)\}.
```

R79 gives $`D_T=O(NL/b^2+F(n,L))`$. On every nonempty interval
$`17B_D\le b\le\sqrt{NL/F(n,L)}`$, count and depth match
simultaneously: $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$. Its high-precision branch has
$`D_T=O(NL/b^2+L)`$ when $`L\gt N/n^2`$.

The [depth catalogue](../research/README.md#large-workspace-t-depth)
indexes the earlier theorem and component proofs.

### Complete real frames with shared rotations

For $`q=n^{o(1)}`$, allow arbitrary rotations on the upper q levels. Inside
each resulting prefix subtree, require a shared actual rotation at each lower
depth. At fixed precision and square-root dirty width, the complete frame
has $`D_T=O_L(\log(q+\log(n+4)+2))`$, hence $`o(\log n)`$, with
one clean qubit, $`T=O_L(\sqrt N)`$, $`G=O_L(N)`$ and full return. Polylogarithmic q
gives doubly logarithmic depth. The general family permits independent
deepest-level entries. [Proof and related families](../research/depth/PREFIX_ROTATION_PROFILES.md)

## Complete query transfer and constructive target

A [prescribed scalar sign](../research/depth/PRESCRIBED_SCALAR_TABLES.md)
also suffices under simultaneous exact query bounds across admissible widths.
At fixed L the resulting depth is

```math
D_T=O_L\!\left(N/b^2+\log\log(n+4)+(1+\log^*(n+2))d(n+O_L(1))\right).
```

Here d is the additive source depth; routing and the bounded-height terminal
are paid. A general doubly logarithmic source remains open.

For a sufficiently accurate fixed-error scalar source, exact preservation of
the original computational address sectors on every physical input permits
[pair-reflection amplification](../research/depth/PRESCRIBED_SCALAR_TABLES.md#pair-reflection-scalar-amplification).
Writing $`h=1+\log^*(n+2)`$, the resulting conditional depth is

```math
D_T=O_L\!\left(N/b^2+\log\log(n+4)
 +h\log(h+2)d(n+O_L(1))\right).
```

The same allocation and count budgets hold, with actual inverses and no
initialized vote bank or controlled source.

Two weaker algebraic targets are a [zero-clean full-space signed gap](../research/depth/PRESCRIBED_SCALAR_TABLES.md#full-space-signed-gap-source)
and a [one-clean coherent decision](../research/depth/PRESCRIBED_SCALAR_TABLES.md#one-clean-decision-source).
For the latter, J initializes only c and $`M=J^\dagger V^\dagger Z_cVJ`$
has sign $`D_f\otimes I`$ and a uniform positive gap on the entire
initialized code. V may move the address and leave garbage; the rejected
clean sector needs no prescribed sign. Amplify on the same signal before
paid XOR promotion. For either source, a first improvement at
$`b=\lceil K_L\sqrt N\rceil`$ needs only T-count and dirty width
$`O(2^{a/2})`$, Clifford count $`O(2^a)`$ and depth d(a). The
[square-root-width corollary](../research/depth/PRESCRIBED_SCALAR_TABLES.md#square-root-width-signed-gap-transfer)
gives $`T=O_L(\sqrt N)`$, $`G=O_L(N)`$ and

```math
D_T=O_L(\log\log(n+4)+h\log(h+2)d(n+C_L)).
```

Thus $`d(a)=O(\log a/\log\log a)`$ already suffices for
$`o(\log n)`$ frame depth. This is conditional, does not interpolate
to every width, and does not establish optimal depth even for a
constant-depth source. The [two-band gap screen](../research/depth/FOUR_LAYER_SCALAR_GAP.md#two-band-signed-gap-sources)
rules out a constant-gap source with at most two raw T bands for all
tables; it is not a general growing-depth lower bound.

The [raw-decision count bound](../research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md#raw-decision-count)
forces worst-case $`\Omega(2^{a/2})`$ T-count with arbitrary width, garbage and arbitrarily small strict bias at a fixed Pauli readout.
It is tight: [controlled Gauss decisions](../research/depth/GAUSS_SUM_QUERY_ENCODING.md#constant-depth-raw-decisions)
attain square-root T-count/width and constant depth for every table's strict basis bias.
Their compressed sign fails on fixed product superpositions. A growing-depth invariant must
use more, such as a uniform bias, coherent sign or return; single-address predicates also admit
[constant-depth fixed-error reflections](../research/depth/EXPANDER_ZERO_REFLECTION.md).

The [four-call shear reduction](../research/depth/PRESCRIBED_SCALAR_TABLES.md#self-inverse-shear-tables) makes scalar entry signs
for known self-inverse binary shears sufficient, with actual inverses,
fixed phases, borrowed address bits and error paid. A cheap physical shear
does not provide its coherent binary-index entry readout.

[Conditional R98](../research/depth/QUERY_DEPTH_FRAME_TRANSFER.md)
accepts a zero-clean arbitrary-table XOR query whose completed address
and helpers return. Intermediate address motion is allowed. For every
a-bit address, $`Q=2^a`$ rows, q output bits and every eligible dirty
allocation $`B\ge C_Q(a+2)`$, its hypothesis requires simultaneous

```math
T=O(\sqrt{Qq}+Qq/B),\quad G=O(Qq),\quad
D_T=O(d(a)+Qq/B^2),\quad \text{extra width}\le B.
```

For fixed integer $`L\ge6`$, that hypothesis gives both frame families,
at the original one-clean thresholds, optimal-order T-count,
$`G=O_L(N)`$ and

```math
D_T=O_L\!\left(N/b^2+d(2n+10)+\log_2\log_2(n+4)\right).
```

The [approximate extension](../research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md)
needs only constant source accuracy for each fixed L, with uniform
simultaneous budgets at that requested accuracy. It supplies no source.
R91 has logarithmic d; a smaller d is not obtained merely by relaxing
the address convention. No uniform growing-L improvement is asserted.

An alternative sufficient target, for $`m=2^p`$, is the complete indicator

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle,\qquad |Y|=m,
```

at every fixed tolerance tau, with zero clean work,
$`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$ and
$`D_T=o(\log(p+2))`$. Accuracy is in full operator norm on arbitrary
coherent address, output, dirty and reference inputs.

[Shared-layer fusion](../research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#shared-layer-indicator-fusion)
allows fixed singly exponential source costs on a constant fraction
of the address. One source and its inverse serve every parent in a
layer, and the full error is at most eight times the sum of the
stage-source tolerances. A common cost exponent suffices over fixed
tolerances. Exact address preservation permits one fixed source error.

## Paid source interfaces and their limits

The missing operation is a complete shallow native source. The reductions
below charge its readout and work return; none supplies that source or
changes the general bounds above. Use each proof home's full error and
resource ledger at the requested source precision.

| Interface | Authoritative contract and recovery |
|---|---|
| Prescribed scalar entry | [Self-inverse shear tables](../research/depth/PRESCRIBED_SCALAR_TABLES.md#self-inverse-shear-tables) suffice through four calls with no additional initialized work. Their scalar readout, rather than applying or inverting the bank shear, is the remaining source problem. |
| Complete Pauli derivative | [Observable transport](../research/depth/OBSERVABLE_TRANSPORT_SOURCE.md) and the [phase/indicator equivalence](../research/depth/BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) include arbitrary input and helper action, actual adjoints and complete return. The canonical phase is an equivalent target, not a cheaper primitive. |
| Encoded dirty copy | [Joint carrier/copy interface](../research/depth/ENCODED_DIRTY_COPY_INTERFACE.md): a complete echo cancels an abstract decoder on returned dirty work. The [exact grading criterion](../research/depth/DIRTY_COPY_GRADING.md) permits noncommuting families, while [native copy screens](../research/depth/NATIVE_COPY_INTERACTION.md) restrict specified interactions. The encoded transition and physical output interaction still need jointly paid native words. |
| Accepted projector Gram | A normalized source must supply $`K^\dagger K\approx\Pi_{\rm eq}/2`$ on every data/dirty input. [One-clean recovery](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout) and [joint two-clean recovery](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout) permit arbitrary accepted private fibers and intermediate address motion; their ledgers include leakage, actual inverses and dirty promotion. |
| Alternative complete readouts | The [single HRT residual probe](../research/depth/DIRTY_CHOI_OBSERVABLE_SOURCE.md#one-bank-polarization), [partial central character basis](../research/depth/COHERENT_CENTRAL_TRANSPORT.md#partial-character-diagonalizer) and [balanced two-unitary acceptance](../research/depth/COUPLED_ORDER_GRAM_SOURCE.md#balanced-two-unitary-gram-readout) retain distinct sufficient promises and paid extraction words. |
| Clean-source promotion | [Serial promotion](../research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) and [spectral promotion](../research/depth/MSS_SPECTRAL_PROMOTION.md) retain their own clean-count, native-sharing, precision and leakage hypotheses. A completed source boundary alone does not establish native address sharing. |

The [route decisions](../research/ROUTE_DECISIONS.md) own the assessment of
each mechanism, its stopping reason and the condition for reopening it.
Quantitative exclusions and component ledgers remain in their proof homes.
The [component catalogue](../research/README.md#constructive-components-and-source-candidates)
and [scoped-limit catalogue](../research/README.md#source-specific-limits-and-failed-compositions)
locate the detailed constructions and failed words. A source-class obstruction,
logical-query count or component improvement does not establish unrestricted
native-depth optimality. The [working checkpoint](../WORKSPACE.md) selects the
next complete-operator question.

For the zero-clean Pauli-source interface, the
[three-band source theorem](../research/depth/THREE_BAND_LOOKUP_SOURCE_OBSTRUCTION.md)
excludes native source depth at most three for addresses of at least seven
bits at complete observable error below
$`\sin(\pi/8)/(2\sqrt2)`$. Clifford interlayers, intermediate address
motion and arbitrary dirty width are unrestricted. This is a constant
bound on a source whose Pauli derivative is lookup; it is not a growing
frame-depth lower bound or a theorem about initialized-clean sources.

## Constant-depth frontier

### Zero-clean error floors and fixed-clean thresholds

For the complete frame $`W_n=I_{2^{n-1}}\otimes R_y(\theta)`$,
$`\theta=\tfrac12\arccos(1/3)`$, the unrestricted zero-clean
[R73 theorem](../research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md) has
the sharp three-layer floor
$`c_0=\sqrt{2-(2+\sqrt2)/\sqrt3}\approx0.169714`$.
[R75](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md) gives the
four-layer floor $`\rho/2\approx0.135299`$, where
$`\rho=\sin(\pi/8)/\sqrt2`$, at arbitrary dirty width.
The [rare two-level cycle](../research/depth/ROBUST_DECISION_BIAS.md#zero-clean-six-band-floor) raises the zero-clean floor to six bands at $`\eta\lt\sin(\pi/24)/2`$ for sufficiently large n.
For typical [complete deepest-corner frames](../research/depth/ROBUST_DECISION_BIAS.md#complete-return-four-band-floor),
four bands have the sharp limiting worst-case error $`\sqrt{2-\sqrt2}`$
at any fixed clean count and arbitrary dirty width. A [two-table witness](../research/depth/ROBUST_DECISION_BIAS.md#complete-return-two-table-frame) allows any fixed polynomial clean budget. With one clean input,
five bands require [large axis intersections](../research/depth/ROBUST_DECISION_BIAS.md#complete-return-five-band-intersection) at both ends.

### One-clean midpoint restrictions

[R76](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md) considers
$`V=V_2V_1`$, each half of depth at most two, and the actual
encoded clean observable $`S=V_1Z_CV_1^\dagger`$.
A Pauli S gives floor $`\rho/2`$; a Hermitian-Clifford S gives
an explicit positive width-independent floor. Terminal return alone
implies neither premise. An exactly returned occupied block of finite
projective order separately has floor $`1/196608`$; generic
Clifford+T circuits need not have finite projective order.

### One-clean inner spread

[R77](../research/depth/ONE_CLEAN_INNER_SPREAD.md) gives a positive
floor for every fixed number t of first-band axes anticommuting with the
incoming clean Pauli, independently of dirty width and crossing rank.
Minimizing over both two-band encoders and every native representation,
error $`2^{-L}`$ requires $`t_{\min}=\Omega(\log\log L)`$.
For the rotation witness, the [one-clean interface argument](../research/depth/ONE_CLEAN_NATIVE_INTERFACES.md)
sharpens this to a linear-in-L count of clean-active axes at both ends.

### Outer-axis overlap

[R97](../research/depth/OUTER_OVERLAP_RIGIDITY.md) gives a positive
floor depending only on $`r=\min\{\mu_A,\mu_B\}`$. For outer native
axis spaces, $`I=L_A\cap L_B`$ and $`\mu_A`$ is
$`\dim I`$ minus the number of individual A axes in I; define
$`\mu_B`$ similarly. Dirty width, inner spreads, crossing rank and
$`\dim I`$ may grow. Arbitrarily accurate four-layer sequences
must have both misalignments diverge in every admissible representation.
Exact return and large representation parameters alone give no obstruction;
[alternative native splits](../research/depth/FOUR_LAYER_REPRESENTATION_LIMITS.md)
must remain allowed.

### Middle crossing rank and the remaining gap

For the R74 rotation witness with fixed clean count a and middle crossing
rank r, the [precision estimate](../research/depth/FOUR_LAYER_LINEAR_PRECISION.md)
is $`2^{-O_a(r+1)}`$. Achieving error $`2^{-L}`$ requires $`r_{\min}=\Omega_a(L)`$
as L grows. With one clean qubit, the [adjacent crossings and common middle
directions](../research/depth/ONE_CLEAN_NATIVE_INTERFACES.md) also grow linearly
in L. These rank requirements allow large ranks within four layers.

At zero clean count the [scalar gap](../research/depth/FOUR_LAYER_SCALAR_GAP.md)
settles the four-layer Hilbert–Schmidt question and rules out an exactly
silent single-address response at seven or more address bits. An unrestricted
one-clean five-layer lower bound for the rotation witness remains open.

Two further restrictions apply at fixed precision. Some targets require
$`\Omega(\sqrt N)`$ T gates [outside the two widest layers](../research/depth/NATIVE_LAYER_COUNT_DISTRIBUTION.md)
under an $`O(\sqrt N)`$ T budget. An [explicit quartic family](../research/depth/QUARTIC_FRAME_MIXING.md)
requires linearly growing basis-mixing resources and admits constant native
depth. Both restrictions are compatible with constant depth at large width.

## Boolean-query lower bounds and the focused scout

| Actual native model | Established degree statement |
|---|---|
| At most one clean qubit, two bands | [Boolean XOR degree at most two](../research/depth/TWO_LAYER_DIAGONAL_RIGIDITY.md); higher degree has full-isometry floor $`\sin(\pi/8)/(2\sqrt2)`$ |
| Zero clean qubits, three bands | [Boolean XOR degree at most four](../research/depth/THREE_LAYER_NATIVE_STRUCTURE.md), with the same floor for higher degree |
| At most one clean qubit, three bands | [Exhaustive original-pair/no-pair closure](../research/depth/BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md#9-complete-native-three-corollary-with-the-no-pair-constant-retained), sharpened by [central-parity reduction](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md): sign degree at most $`D_{\rm np}=D_3+576=480\,2^{43}+903`$, nonconstant XOR degree at most $`D_{\rm np}-1`$, with rounding radius strictly below $`2^{-530690}`$ |

These statements allow unrestricted intermediate Clifford/address mixing
and arbitrary dirty width. The tiny three-layer radius does not exclude
every fixed tolerance. The paired branch separately has sign/XOR bounds
six/five, or five/four at common active rank at least four.
[Two-call frame corner reductions](../research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md)
transfer suitable unrestricted query lower bounds; source-specific
grammars do not become unrestricted bounds by that transfer.

The [original physical-pair scout](../research/depth/FOUR_BAND_ORIGINAL_PAIR_SCOUT.md)
retains actual independent signed lists and the physical clean code.
Fixed-endpoint/orbit closure, literal pair removal and
[center transport](../research/depth/FOUR_BAND_CENTER_TRANSPORT.md)
are scoped components. The affected-silent reduction leaves at most
$`\lfloor4k/3\rfloor`$ original axes per band at pairing nullity k;
it adds no depth or workspace and is entirely Clifford for $`k\le2`$.
The [coupled six-wire trial](../research/depth/FOUR_BAND_CENTER_TRANSPORT.md#coupled-six-wire-endpoint-obstruction)
now has a structural endpoint obstruction: its ten-T, two-band middle
transports a probe to four Paulis of magnitude one half, forcing complete
clean-isometry error at least $`1/4`$ for every one-band output endpoint,
even with spectator mixing and arbitrary logical action, for all $`2^{10}`$
T/T-dagger orientations of the fixed middle axes. This candidate
is paused before Boolean return. A changed actual middle or a paid
non-Clifford endpoint is required to reopen it. The coefficient argument
can inspect longer words, but extra non-Clifford layers can repair its
witness. No general four-band theorem or fixed-accuracy growing-depth
bound follows; the scalar zero-test construction rules out conjunction
alone as such a witness.

## Earlier endpoint constructions

<a id="proved-linear-count-families"></a>
<a id="bounded-diagonal-factorization"></a>
<a id="a-regular-fixed-tree-cayley-family"></a>
<a id="source-program-resynthesis"></a>

The [research catalogue](../research/README.md#constant-clean-high-precision-count)
retains promised families, fixed-tree/Cayley/scattering reductions,
source programs, diagonal/reflection architectures and their evidence.
Their method-specific missing steps are not gaps in Result C.
Proofs, negative results and historical outcomes remain at their original
paths; [route decisions](../research/ROUTE_DECISIONS.md) distinguish
current questions from superseded instructions.

## Full-circuit requirement

Charge source construction, tables, basis changes, queries, actual inverses,
amplification and workspace return on the same circuit, at the actual
requested precision. No QRAM, supplied precision state, measurement,
reset, postselection or uncharged coherent evaluator is permitted.
Certified preprocessing has no uniform efficiency claim. A lower bound
must apply to this full model; failure of one source, query interface,
schedule or factorization establishes only its stated limitation.
