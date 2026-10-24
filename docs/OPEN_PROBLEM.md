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
| [Scalar zero reflection](../research/depth/EXPANDER_ZERO_REFLECTION.md) | Zero clean work, $`T,G=O(pL_\tau+L_\tau^2)`$, $`b=O(p)`$, $`D_T=O(L_\tau^2)`$, with $`L_\tau=1+\log_2(1/\tau)`$. [Coupled quarter-phase sources](../research/depth/COUPLED_QUARTER_PHASE_SOURCE.md) and [full matching stages](../research/depth/QUARTER_PHASE_MATCHING_LIMITS.md) are paid components; four-path acceptance retains coherent error, while the changed expander gap covers one row only. The complete source remains unpaid |
| [Complete fermionic lift](../research/depth/BANK_MOTION_WALSH_INDICATOR.md#expander-fermionic-lift) | Native lifting and arbitrary auxiliary-mode decoding are paid. The full spectral-phase contract and sparse-column obstruction force growing filter length for the literal fixed-auxiliary expander word; displayed coherent centering also has linear address depth |
| [Observable transport](../research/depth/OBSERVABLE_TRANSPORT_SOURCE.md) | Two calls extract a complete addressed observable with charged return; nonlinear observable synthesis is still missing |
| [MOD3 joint Pauli transport](../research/depth/MOD3_ZETA_TRANSPORT.md#three-band-transport-exactness) | At most three source T layers give either exact complete transport or error at least $`\sin(\pi/8)/(2\sqrt2)`$, with arbitrary dirty width and register motion. Two layers are excluded for p at least three; unrestricted exact three-layer existence remains open. [Joint two-coordinate transport](../research/depth/WEAK_PAULI_ORBIT_COMPRESSION.md#two-coordinate-lift) is paid at twelve T gates, three layers and no extra work, but does not cover the remaining address coordinates |
| [Complete coefficient assembly](../research/depth/WEAK_PAULI_ORBIT_COMPRESSION.md#complete-coefficient-shear) | The four-band nested source remains a partial component. A separate two-band shear completes all p=3 coefficients with 28 T gates and five source layers, no helpers and exact complete return. Its uniform complete source has $`T=4m-4`$, $`G=pm+8m+36`$ and emitted depth $`2p-1`$. This linear-depth recurrence is paused; the unpaid operation is a sublogarithmic joint schedule, not coefficient coverage |
| [Paired projector-Gram readout](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout) | One clean signal with accepted Gram $`K^\dagger K=\Pi_{\rm eq}/2`$ suffices; accepted unitary fibers and intermediate address motion are unrestricted. A shallow source selecting this lookup sector is not supplied |
| [Joint two-clean Gram readout](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-clean-joint-gram-readout) | Joint-code reflection is Clifford. Four-call recovery and full dirty promotion give error $`128\sqrt\delta`$; $`\delta\le\tau^2/2^{20}`$ supplies tolerance tau/8 before fusion. The source is still unpaid |
| [Squared-product interface](../research/depth/COUPLED_ORDER_GRAM_SOURCE.md#squared-product-interface) | Successful coupled-order interference already gives a direct zero-clean reflection. More generally, a [balanced two-unitary source](../research/depth/COUPLED_ORDER_GRAM_SOURCE.md#balanced-two-unitary-gram-readout) with Gram error delta supplies a squared relative word with error at most $`4\delta+4\sqrt\delta`$. These reductions supply no qualifying source |

For the one-clean interface, let $`n_c=pm`$ denote one carrier-bank size.
[Four U/actual-adjoint calls](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout) followed by serial dirty promotion give the
complete zero-clean ledger

```math
\begin{aligned}
T_{\rm out}&\le32T_U+48m+8,&
G_{\rm out}&\le32G_U+16n_c+98m+136,\\
D_T({\rm out})&\le32d_U+32,&
b_{\rm out}&\le2n_c+b_U+2.
\end{aligned}
```

Gram error delta gives complete error at most
$`64\sqrt\delta`$; request $`\delta\le\tau^2/4096`$.
The six-call alternative retains its linear-error and 48-call promoted ledger.
All costs use the requested precision. The normalized lookup sector remains unpaid.

The relevant exclusions have different hypotheses:

| Source class | Established exclusion and boundary |
|---|---|
| General two-band unitary compression | [Rank–norm theorem](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen): exact $`K=\alpha I_p`$, $`\alpha>0`$, needs $`k\ge(2/3)\log_2p-(1/3)\log_2\log_2p-O(1)`$. Fixed-error growing-clean sources remain outside its shrinking approximate floor |
| Native read-only two-band compression | [Complete-bank Fourier bound](../research/depth/READ_ONLY_TWO_BAND_SOURCE_GAP.md): amplitude $`1/\sqrt2`$ and relative error at most $`1/8`$ require $`k\ge\lceil(p-3)/2\rceil`$ when actual Clifford boundaries fix the address. No flatness, exact return or dirty-width bound is assumed |
| One-clean two-band projector-Gram source | [Principal-angle screen](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#two-band-projector-gram-gap): for $`p\ge3`$, Gram error at most $`1/256`$ is impossible even with intermediate address motion. More source bands or another clean allocation are outside this screen |
| Bounded-clean two-band projector-Gram source | [Direct Gram-defect screen](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#bounded-clean-two-band-gram-screen): tiny positive asymptotic error floor for each fixed clean count k, using $`M_k=2^{2k+1}(2^k-1)`$ Clifford summands; $`M_2=96`$. Not an exclusion at every moderate tolerance or for growing/accuracy-dependent k |
| Actual Clifford signal image | [Compression screen](../research/depth/CROSSED_CUBIC_KERNEL_SOURCE.md#clifford-signal-gram-gap): if $`U^\dagger Z_aU`$ is Clifford, Gram error is at least $`1/(2\sqrt2)`$ for $`p\ge1`$, including arbitrary Clifford address motion. General higher-band words need not satisfy this premise |
| Three-layer source with bank-independent final block | [Bank-fiber theorem](../research/depth/COUPLED_FINAL_ENCODER.md#bank-independent-final-boundary): for p at least three, a completed final one-band block independent of the nonconstant coefficient bank cannot supply exact transport or error below the retained three-layer gap. Arbitrary middle Clifford boundaries, dirty width and a fixed output Clifford are covered. Bank-dependent final blocks and more source layers remain outside |
| Specified full-bank X-motion final encoder | [Rotated-bank centralizer theorem](../research/depth/COUPLED_FINAL_ENCODER.md#full-bank-x-motion-boundary): for p at least three, the paid syndrome-CNOT/transversal-X-rotation encoder has no one-band middle completion, even with arbitrary Clifford core and identity helpers. Its affine bank-X response code gives a nonisotropic witness. This excludes that bank-moving word, not general bank-dependent final blocks or four-or-more-layer sources |
| Syndrome reflections around zero or fixed-row acceptance | [Thin-source bounds](../research/depth/COUPLED_FINAL_ENCODER.md#syndrome-star-pulse-boundaries) and [normalized broad trial](../research/depth/COUPLED_FINAL_ENCODER.md#broad-baseline-syndrome-source): the paid broad trial has exact ideal Gram distance one half for p at least three. A baseline-relative rank bound excludes fewer than $`(m-1)/4`$ thin controlled-reflection insertions into a word whose deletion baseline is the fixed-row source, with arbitrary interlayers. A different lookup-dependent baseline or jointly synthesized growing-call word remains outside; call count is not native depth |
| Specified shifted-syndrome commutator | [Four-band native trial](../research/depth/COUPLED_FINAL_ENCODER.md#shifted-syndrome-wide-commutator): 4m T gates and no helpers, but complete Pauli-transport error at least square root two for p at least two, even at arbitrary rotation angles on the same axes; the native p=1 gap is at least one quarter. The missing row response excludes this commutator family without a thin-source assumption. A changed middle action, axes, full word or accepted Gram remains outside |
| Wide Clifford filters | [Sharp lookup correlation](../research/depth/CLIFFORD_LOOKUP_FIXED_SPACE.md) is $`(p+1)/2^p`$ for $`p\ge2`$; normalized Clifford averages cannot separate the two full lookup eigenspaces by a constant gap. Coefficient and call bounds are not native-depth bounds |

Completed address return does not establish native read-only sharing.
The [route decisions](../research/ROUTE_DECISIONS.md) record the precise
limitations and reopening conditions for normalized kernel, signal-interleaved,
scalar/Bell, orbit and fermionic sources. The
[working checkpoint](../WORKSPACE.md#main-unfinished-route-and-exact-missing-operation)
identifies the next complete-operator question. None of the recorded source
trials supplies the missing complete shallow source or changes the general frontier.

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
even with spectator mixing and arbitrary logical action. This candidate
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
