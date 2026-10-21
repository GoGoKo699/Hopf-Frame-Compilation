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
| Strongest general depth upper | R96: $`O(NL/b^2+L+\log_2(n+2))`$ on the same circuit as optimal-order T-count and $`G=O(NL)`$ |
| Selected Result D | R79 remains the selected all-precision count/depth theorem; R96 is the stronger research result |
| Optimal large-width depth | At fixed accuracy and square-root dirty width, logarithmic general upper depth versus constant unrestricted lower depth remains open |
| Main constructive question | A complete approximate dirty indicator with sublogarithmic address depth and simultaneous resource budgets |
| Focused lower question | Actual original physical-pair transport beyond the proved four-band eligibility classes; no general four-band degree ceiling is known |

The [current readout assessment](../research/depth/READOUT_AND_NATIVE_SOURCE_REVISION_RESULT.md)
separates general theorems, paid components, source-specific limitations and
the first unpaid operation. The [route decisions](../research/ROUTE_DECISIONS.md)
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

**Strongest general upper: R96.** For every integer $`n\ge1,L\ge6`$,
the [all-precision theorem](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md)
uses one external clean qubit and retains the literal thresholds

| Complete target | Sufficient dirty allocation |
|---|---|
| Real Hopf frame | $`b\ge17(L+n+7)`$ |
| Prescribed phase-dressed complex magnitude frame | $`b\ge17(L+n+8)`$ |

On one native circuit, with absolute constants and full error at most
$`2^{-L}`$, it gives

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\ell),\qquad \ell=\log_2(n+2).
```

The T-count is worst-case optimal at eligible allocations. The depth
formula is an upper bound, not an unrestricted optimality theorem.
At sufficient $`b=\Theta(\sqrt{NL}+L+n)`$ it gives
$`D_T=O(L+\ell)`$; at fixed accuracy and square-root dirty width,
the unresolved interval is $`\Omega(1)`$ versus $`O(\log n)`$.
This concerns the prescribed complex magnitude completion, not arbitrary
complex-unitary synthesis. No practical crossover or scalable native
emitter is claimed.

**Selected Result D: R79.** The
[uniform subtree-depth theorem](../research/depth/UNIFORM_SUBTREE_DEPTH.md)
retains the same count, Clifford, error and one-clean width guarantees.
With $`B_D=L+n+7`$ for real frames or $`B_D=L+n+8`$ for the
complex completion, define

```math
F(n,L)=\min\{nL,(L+\ell)\log_2(n+L+2)\}.
```

R79 gives $`D_T=O(NL/b^2+F(n,L))`$. On every nonempty interval
$`17B_D\le b\le\sqrt{NL/F(n,L)}`$, count and depth match
simultaneously: $`T^\star=\Theta(NL/b)`$ and
$`D_T^\star=\Theta(NL/b^2)`$. Its high-precision branch has
$`D_T=O(NL/b^2+L)`$ when $`L\gt N/n^2`$.

R78's subtree decoder, R80–R85's promised-family/terminal results,
R86–R90's fixed-accuracy improvements and R92–R95's precision refinements
retain independent component proofs. Their statements and dependencies
belong in the [depth catalogue](../research/README.md#large-workspace-t-depth),
not in competing current-frontier narratives.

## Complete query transfer and constructive target

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

The main sufficient target, for $`m=2^p`$, is

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle,\qquad |Y|=m,
```

at every fixed tolerance tau, with zero clean work,
$`T,b_{\rm dirty}=O_\tau(m)`$, $`G=O_\tau(m^2)`$ and
$`D_T=o(\log(p+2))`$. Accuracy is in full operator norm on arbitrary
coherent address, output, dirty and reference inputs.

[Polynomial-overhead fusion](../research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
allows discovery at fixed polynomial overhead on
$`p-O_\tau(\log(p+2))`$ address bits. Eight source/actual-inverse
calls incur at most eight times source error and add
$`O_\tau(\log\log(p+4))`$ depth. A complete sublogarithmic source,
not count normalization or its controlled use, remains unpaid.

## Paid source interfaces and their limits

| Interface or component | What is paid; what is not |
|---|---|
| [Clean-source promotion](../research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) | Bounded clean banks and constant-amplitude complete blocks can be removed with charged actual inverses, leakage and arbitrary off-code cancellation; this supplies no shallow source |
| [Spectral promotion](../research/depth/MSS_SPECTRAL_PROMOTION.md) | At fixed tolerance, depth $`O_\tau(d_V+k+1)`$ and polynomial overhead for $`k=o(\log(p+2))`$; native address-sharing and uniform precision budgets are premises. Exact sector invariance avoids the general leakage multiplier |
| [Scalar zero reflection](../research/depth/EXPANDER_ZERO_REFLECTION.md) | With $`L_\tau=1+\log_2(1/\tau)`$, zero clean work, $`T,G=O(pL_\tau+L_\tau^2)`$, $`b=O(p)`$, $`D_T=O(L_\tau^2)`$. Native address motion and whole-bank error prevent an automatic indicator |
| [Observable transport](../research/depth/OBSERVABLE_TRANSPORT_SOURCE.md) | Two calls extract a complete addressed observable with charged return; nonlinear observable synthesis is still missing |
| [Paired projector-Gram readout](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#paired-projector-readout) | One clean signal with accepted Gram $`K^\dagger K=\Pi_{\rm eq}/2`$ suffices; accepted unitary fibers and intermediate address motion are unrestricted. Normalized shallow U is not supplied |

For the last interface, let $`n_c=pm`$ denote one carrier-bank size.
Six U/actual-adjoint calls followed by serial dirty promotion give the
complete zero-clean ledger

```math
\begin{aligned}
T_{\rm out}&\le48T_U+48m+8,&
G_{\rm out}&\le48G_U+16n_c+98m+152,\\
D_T({\rm out})&\le48d_U+32,&
b_{\rm out}&\le2n_c+b_U+2.
\end{aligned}
```

Gram error delta gives complete error at most
$`32\sqrt{24\delta}`$; request $`\delta\le\tau^2/24576`$.
The stronger block promise
$`\|K-B_x/\sqrt2\|\le\eta\le1/8`$ gives error
$`32(3\sqrt2\eta+2\eta^2)`$ instead. All source costs are evaluated
at the requested precision. Recovery and full return are paid; obtaining
the constant-amplitude source is the first unpaid native operation.

The relevant exclusions have different hypotheses:

| Source class | Established exclusion and boundary |
|---|---|
| General two-band unitary compression | [Rank–norm theorem](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen): exact $`K=\alpha I_p`$, $`\alpha>0`$, needs $`k\ge(2/3)\log_2p-(1/3)\log_2\log_2p-O(1)`$. Fixed-error growing-clean sources remain outside its shrinking approximate floor |
| Native read-only two-band compression | [Complete-bank Fourier bound](../research/depth/READ_ONLY_TWO_BAND_SOURCE_GAP.md): amplitude $`1/\sqrt2`$ and relative error at most $`1/8`$ require $`k\ge\lceil(p-3)/2\rceil`$ when actual Clifford boundaries fix the address. No flatness, exact return or dirty-width bound is assumed |
| One-clean two-band projector-Gram source | [Principal-angle screen](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-band-projector-gram-gap): for $`p\ge3`$, Gram error at most $`1/256`$ is impossible even with intermediate address motion. More source bands or another clean allocation are outside this screen |
| Wide Clifford filters | [Sharp lookup correlation](../research/depth/CLIFFORD_LOOKUP_FIXED_SPACE.md) is $`(p+1)/2^p`$ for $`p\ge2`$; normalized Clifford averages cannot separate the two full lookup eigenspaces by a constant gap. Coefficient and call bounds are not native-depth bounds |

Completed address return does not supply native read-only sharing.
The direct read-only two-band/current-spectral-wrapper combination fails;
this does not lower-bound every promotion scheme or unrestricted depth.
The next source question is a normalized higher-band Gram mechanism, or
another complete wide non-Clifford readout with a paid native schedule.

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
Thus integer $`L\ge3`$ requires at least five T layers with zero
clean work. Sharpness is not claimed for the four-layer floor.
R73 also gives a positive three-layer threshold at every fixed clean
budget; a separate explicit one-clean witness needs three layers for
$`L\ge8`$. Zero leaf phases include these witnesses in the complex family.

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
This excludes bounded spread, not growing spread or general four layers.

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

[R74](../research/depth/FOUR_LAYER_CROSSING_RANK.md) gives a positive
four-layer gap for every fixed clean budget a and bounded middle crossing
rank, with explicit rigidity gap $`2^{-O_a(r^2)}`$.
For $`2^{-L}\lt c_0`$, its minimum over decompositions obeys
$`r_{\min}=\Omega_a(\sqrt L)`$. Together with R77, all four
layer sizes must grow in a hypothetical arbitrarily accurate one-clean
sequence. Growing rank is not excluded. The one-clean unrestricted
five-layer lower bound and the separate four-layer Hilbert–Schmidt gap
remain open. These constant-layer results do not match R96 asymptotically.

## Boolean-query lower bounds and the focused scout

| Actual native model | Established degree statement |
|---|---|
| At most one clean qubit, two bands | [Boolean XOR degree at most two](../research/depth/TWO_LAYER_DIAGONAL_RIGIDITY.md); higher degree has full-isometry floor $`\sin(\pi/8)/(2\sqrt2)`$ |
| Zero clean qubits, three bands | [Boolean XOR degree at most four](../research/depth/THREE_LAYER_NATIVE_STRUCTURE.md), with the same floor for higher degree |
| At most one clean qubit, three bands | [Exhaustive original-pair/no-pair closure](../research/depth/BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md#9-complete-native-three-corollary-with-the-no-pair-constant-retained), sharpened by [central-parity reduction](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md): sign degree at most $`D_{\rm np}=D_3+576`$, nonconstant XOR degree at most $`D_{\rm np}-1`$, with conservative rounding floor $`2^{-530690}`$ |

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
Genuinely coupled transport, an actual output endpoint and full Boolean
return remain the focused test. No general four-band theorem follows.

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
