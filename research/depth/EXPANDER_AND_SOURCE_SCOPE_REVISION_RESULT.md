# Revision: spectral promotion, zero tests and complete sources

This continuation starts from main commit
`49d58b07edb99f4ad4cabfee5993ecbe3df1e9dd`, tree
`452d9a8ac1d57d818c8b3237f2894ed5605ed382`, merged as PR #215 on
20 October 2026. The [preceding revision](PARALLEL_PROMOTION_AND_KERNEL_REVISION_RESULT.md)
paid parallel promotion and separated complete source identities from
their missing native readout. This revision pays a larger conditional
promotion window and a useful zero-test primitive, and strengthens
several scope tests. It does not claim a
new complete indicator or an improved unrestricted frame frontier.

## 1. Destination, strongest results, components and limitations

| Role | Current assessment |
| --- | --- |
| General destination | Optimal native T-depth for complete frames under simultaneous T-count, Clifford count, dirty width, precision and return budgets |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md) remains unchanged; [R79](UNIFORM_SUBTREE_DEPTH.md) remains selected Result D |
| Concrete constructive target | For p address bits and m equal to two to the p, a complete zero-clean indicator with fixed-accuracy linear T-count/dirty width, quadratic Clifford count, and depth sublogarithmic in p |
| New paid positive component | A single approximate zero reflection has zero clean work, linear count/width and constant T-depth at every fixed tolerance; its full uniform-precision ledger is below |
| New conditional source window | [MSS spectral promotion](MSS_SPECTRAL_PROMOTION.md) pays extraction for clean count k sublogarithmic in p, with constant requested source accuracy under exact sector invariance and the actual native-sharing premise |
| Retained alternative | [Parallel exact packing](PARALLEL_CLEAN_PROMOTION.md) retains its smaller sufficient window and its two-times-source-error bound |
| New arithmetic components | Arbitrary-function addition, odd power-of-two multiplication and interval rotation have paid exact dirty-work constructions; the desired full in-place nonzero-radix packing is still a separate operation |
| Current strongest paired lower result | Actual one-clean three-band paired Boolean sources have sign/XOR degree at most six/five, or five/four at high active rank. The general three-band ceiling also retains its separate no-pair theorem |
| New lower-scout component | Actual two-band middle transport requires overlap on the probe-active native coordinates before noncommuting unmatched orbits can survive. This is eligibility, not a general four-band degree or depth theorem |
| Representation limits | Joint relabeling invariants, bounded operator-Schmidt rank, graph-projector normalization, fixed-fiber packing and Majorana boundary theorems keep their stated premises |

All complete contracts include arbitrary coherent address, output,
dirty-work and reference inputs, literal phase and actual inverse.
Approximate return is charged in the same operator norm as the target
error. A prepared code, trace average or favorable column does not
replace that contract.

The general R96 upper remains, on the same complete circuit,

```math
\begin{aligned}
D_T&=O(NL/b^2+L+\log_2(n+2)),\\
T&=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
\end{aligned}
```

with one clean qubit and the retained real/complex dirty thresholds.
At fixed accuracy and square-root dirty width, the unrestricted
upper/lower interval remains logarithmic versus constant. None of the
new source screens closes it.

## 2. A complete zero reflection is now paid

The [expander construction](EXPANDER_ZERO_REFLECTION.md) imports the
explicit fourteen-generator CNOT group expander in its actual
permutation representation. A local Z test removes its unwanted
uniform-nonzero invariant vector. The proof gives a positive constant
gap on the entire nonzero subspace, rather than inferring one from
Hilbert--Schmidt mixing in another representation.

Four clean selector bits give a Hermitian compression A that fixes
zero and has norm at most a constant c less than one elsewhere. The
same selectors support the literal general-block cubing recurrence

```math
S_{j+1}=S_j R S_j^\dagger R^\dagger S_j,
\qquad PS_jP=A^{3^j}.
```

R is a determinant-one selector phase. Its fifteen Pauli rotations,
parity exposures, synthesis precision and actual inverses are paid.
Conjugating an exact four-selector XOR by this word produces a
complete four-clean zero-test query. If N equals three to the j and
each R occurrence has full-norm error at most nu, its complete error is

```math
\epsilon_V\le2c^N+2(N-1)\nu.
```

Bounded-clean promotion then gives a zero-clean query with error at
most 128 times epsilon_V. A two-call arbitrary-bit phase echo gives
the zero reflection with error at most 256 times epsilon_V. At most
four dirty-padding conjugates cover every input length and cost at
most 1024 times epsilon_V. Every previously clean selector and padding
bit is arbitrary in the final circuit, and its return is included.

Choosing N proportional to L and rotation precision with logarithm
proportional to L, where L is one plus log of inverse tolerance, gives
the simultaneous final ledger

```math
\begin{aligned}
T,G&=O(pL+L^2),\\
b_{\rm dirty}&=O(p),\qquad D_T=O(L^2),\\
\|\widetilde R_p-R_p\otimes I\|&\le\tau,\qquad
R_p=I-2|0^p\rangle\langle0^p|.
\end{aligned}
```

All constants include the small supplied expander gap. Thus fixed
tolerance gives constant T-depth, but no practical crossover or
efficient general emitter is claimed. The uniform precision depth is
quadratic, not linear, in L.

The first unpaid operation has moved from this scalar predicate to
its **complete-bank use**. The native generators move the tested
address. Independent row copies therefore do not meet the existing
shared-address schedule, even though the completed ideal reflections
commute. Serial constant-error rows do not give a uniform full-bank
error bound. Substituting the primitive into isolated guards also
leaves the other full word-query terms in R91/R96/R98 intact.

## 3. Complete-bank revisions must leave the tested representations

The same proof home gives three distinct tests of proposed readout.

First, completed joint permutations of the address and bank labels,
together with arbitrary address-independent quantum bank/helper
operations, preserve the uniform-address subspace. The desired lookup
phase sends a balanced computational bank input orthogonally out of
that subspace. Hence the compression error is at least the requested
amplitude, and complete initialized-unitary error is at least square
root of two. Arbitrary bank-changing operations alone cannot fix this
invariant; an address-sensitive operation is needed.

Second, bounded-choice joint SELECTs have bounded operator-Schmidt
rank across the address versus everything-else cut. The target has
m equal Schmidt coefficients. Arbitrary local address and bank
operations do not increase rank, so a bounded-rank product needs a
number of factors proportional to p at fixed relative error. This is
not native T-depth: p parallel cross-cut CNOTs have rank m in one
Clifford boundary. Wide coupling is an available escape, not a free
solution of its complete readout.

Third, a paid forward quadratic query and a paid zero reflection give
a reflection about the graph of the nilpotent pencil. Ordinary graph
readout still has a small normalization and a large prepared-code
requirement. The [graph audit](EXPANDER_ZERO_REFLECTION.md#11-a-paid-graph-reflection-still-needs-complete-readout)
states the principal-angle and thin-rank bounds with their precise
interfaces. These do not exclude arbitrary wide native coupling or
an alternative word whose erased baseline already implements the
desired action. The first unpaid operation is a full nonvanishing
compression or addressed observable, not its forward graph predicate.

The [latent-center refinement](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#latent-center-readout)
gives a constructive alternative interface. The forward pencil pays
actual Clifford generators of a Pauli subgroup whose center contains
the addressed observable. Both central sectors carry nontrivial
irreducible actions; ordinary invariant-vector filtering selects
neither sector. Conjugation averaging from one fixed bounded probe
has a shrinking addressed coefficient. Neither statement excludes
an alternative full-space central-sign readout. That readout, with
constant amplitude and returned arbitrary work, is the first unpaid
operation in this representation.

These results select a question rather than retire the destination:
can a genuinely wide address-bank operation supply a complete source
outside these invariants, with a constant-fraction coherent readout
and returned arbitrary work?

## 4. The amplitude-one-over-square-root-two recurrence is retained

Let U use one clean signal and let its entire accepted block satisfy

```math
\|J^\dagger UJ-I_p/\sqrt2\|\le\epsilon\le1/8.
```

Here I_p denotes the complete XOR indicator, not an identity matrix.
The [paid phase-S recovery and promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#6-a-one-signal-block-is-another-sufficient-source-interface)
give, for q specified query outputs and common source/inverse costs,

```math
\begin{aligned}
T_{\rm prom}&\le24T_U+48q+8,\\
G_{\rm prom}&\le24G_U+82q+104,\\
D_T({\rm prom})&\le24d_U+32,\\
b_{\rm prom}&\le q+b_U+2,\\
\epsilon_{\rm prom}
 &\le16(3\sqrt2\epsilon+2\epsilon^2).
\end{aligned}
```

The actual inverse, all off-code action and every return error are
included. Under separately proved exact flatness the recovered clean
sector is exactly invariant; the refined promoted error is instead
at most two square root of two times epsilon, with exact signal and
toggle return. For general leakage, the refined bound is unavailable.
The growing-clean version and its native parallel schedule retain
their own complete ledgers in the linked proof homes.

No qualifying new U is supplied here. Two-band sources with uniformly
bounded clean count remain obstructed at some fixed positive relative
accuracy. A fresh [Clifford rank--norm bound](FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen)
strengthens the exact-source threshold to

```math
k\ge\frac23\log_2p-\frac13\log_2\log_2p-O(1).
```

It therefore excludes exact two-band sources throughout the new
sublogarithmic-clean promotion window. This application permits
arbitrary defect coefficients and uses the smaller direct defect,
rather than the old fixed-coefficient gap. Its elementary elimination
bound alone excludes exact one-clean sources for p≥18, and sources
with relative error at most $`2^{-10}`$ for p≥20. The constants are
independent of positive amplitude, dirty width and gate count.

The proof also gives an explicit approximate inequality for each k.
Its excluded tolerance still shrinks as k grows, so fixed-error
growing-clean sources remain open. The old coefficient-gap proof,
higher-band mechanisms and unrestricted sublogarithmic depth retain
their separate scopes. The next source question is genuinely
approximate two-band compression with growing clean allocation, or
a different mechanism, rather than an exact source in this window.

## 5. Spectral paving changes the extraction mechanism

The [new spectral wrapper](MSS_SPECTRAL_PROMOTION.md) avoids the
unpaid full in-place radix conversion. Take h equal to two to the k
private source copies and let K count their source-clean zero keys.
Retain only sectors with K equal to one or two. The controlled global
selective echo is exactly identity outside those sectors; under exact
source-sector preservation its error is at most four times the source
error, independently of h. This is proved for the whole native word,
including the shared predicate and all arbitrary off-code actions.

The resulting Boolean phase has positive support density greater than
one quarter. Its bad projector therefore has density below three
quarters. A global Hadamard makes that projector's diagonal flat.
Marcus, Spielman and Srivastava's partition theorem gives a 64-color
paving; Fourier signs of those colors supply 64 global frames whose
averaged bad projector has norm at most 63/64. Thus the positive
Hamiltonian in the retained filter has a whole-space gap at least one.

These frames are actual non-Clifford native operations. An exact dirty
indicator and a bank-parity echo implement each Boolean sign, returning
the entire arbitrary truth-table bank. Their count and dirty width are
exponential in the packing-register size, and their depth is logarithmic
in that size. The paving can be found by finite enumeration with exact
rational positive-semidefinite certification; no efficient classical
search is claimed. The result does not assume the earlier Clifford-frame
gap or contradict its representation-specific limits.

The proof substitutes these costs into the retained two-clean
Hamiltonian simulation, fixed-point filter, clean promotion and phase
echo. All precisions, actual inverses and return errors are charged.
At each fixed requested final tolerance, the resulting depth is

```math
D_T=O_\tau(d_V+k+1),
```

with source-copy cost proportional to $`2^k`$ and additional
count/width bounded by $`2^{O(k2^k)}`$.
For

```math
k=o(\log(p+2)),
```

these resources fit the polynomial-overhead complete-source envelope.
Exact sector invariance retains a source tolerance depending only on
tau. Generic initialized-isometry error retains its factor of two to
the k and therefore needs the separately stated accuracy envelope.
An exactly flat amplitude-one-over-square-root-two block inherits the
constant-accuracy case after the paid recovery and its extra clean flag.

The [seven-selector refinement](MSS_SPECTRAL_PROMOTION.md#8-a-lazy-seven-selector-filter-improves-uniform-precision)
improves the uniform precision recurrence. Half of its SELECT branches
are identity and half are the paid framed phases. The compressed block
is identity when the desired bit is zero and a positive contraction of
norm at most 63/64 when it is one. General-block cubing uses the same
seven selectors at every level. A returned zero test and final
bounded-clean promotion recover the complete query.

For $`L=1+\log_2(1/\tau)`$, put $`h=2^k`$ and
$`\Lambda=2^{O(h(k+1))}\mathrm{poly}(h,k)`$.
The simultaneous pre-fusion bounds are

```math
\begin{aligned}
d_{\rm out}&=O\bigl(L(d_V+k+1)+L^2\bigr),\\
T_{\rm out}&=O\bigl(L\{h(T_V+q)+\Lambda\}+L^2\bigr),\\
G_{\rm out}&=O\bigl(L\{h(G_V+q)+\Lambda\}+L^2\bigr),\\
w_{\rm out}&=O\bigl(h(b_V+q+k)+\Lambda\bigr).
\end{aligned}
```

The source terms are evaluated at accuracy proportional to tau/L
under exact sector preservation, or tau divided by hL with leakage.
This removes the Hamiltonian-simulation inverse-tolerance factor.
It does not supply a source whose own budgets meet these conditions;
the retained root fusion still uses its fixed-tolerance hypotheses.

This is a genuine improvement of the conditional promotion frontier.
It is not a qualifying source. The first unpaid operation is now that
source's complete logical transition, with native address sharing and
all simultaneous budgets, rather than radix packing or an unspecified
global-frame gap.

## 6. Arithmetic and dirty mixing: useful components, unchanged costs

The [structured packing revision](PARALLEL_CLEAN_PROMOTION.md#7-structured-packing-concentration-is-paid-radix-compression-is-not)
pays concentration of the marked weights. The remaining operation is
a total in-place key permutation taking the all-nonzero digit set
into one computational halfspace. Elementary interval shuffles are
paid, but the current multi-bit radix schedule has linearly many
dependent stages. It also lacks the refined one-active-fiber error
cancellation of the earlier packing.

The [virtual-function arithmetic proof](VIRTUAL_FUNCTION_ARITHMETIC.md)
pays exact addition of an arbitrary function of stable inputs without
first preparing that function in a clean register. Recursive carry
predicates, odd power-of-two multiplication, and invariant-predicate
interval echoes retain all arbitrary inputs and actual inverses.
Their Boolean-query count can be exponential in the control dimension.
That cost can fit the small packing register; applying it separately
to every full-address row can exceed the complete source budget.
Neither an out-of-place query nor paired modular multiplication
automatically supplies the required single-bank radix permutation.
The spectral wrapper bypasses that operation; these arithmetic
components retain their separate uses and full-space contracts.

The [one-pass filter audit](COHERENT_COUNTER_SOURCE.md#9-one-pass-row-filters-clean-dimension-error-and-dirty-mixing)
separates fixed finite-dimensional fibers from arbitrary dirty mixing.
The latter has sharp component clean-dimension tradeoff

```math
r=\Theta(\min\{p+1,\rho^{-2}\}).
```

Its reflected row error is at most twice rho. The stronger fixed-error
logarithmic dimension bound needs invariant finite-dimensional fibers.
A fully dirty Majorana commutator gives an exact complete-bank identity
and exact return, but its present native depth is at most 16 times
p plus one, worse than the retained depth-four-p construction. The
separate cycle lower bound concerns actual Majorana-preserving,
read-only-address boundaries, not a resynthesized whole echo.

## 7. The focused original physical-pair scout

The [native transport proof](FOUR_BAND_CENTER_TRANSPORT.md) keeps the
actual independent lists and actual Clifford endpoint probes. A sole
isolated anticommuting orbit forces an actual axis merger in the other
band. Its two T rotations become a paid Pauli Clifford or cancel,
with literal scalar and full physical return preserved. The first
multiple-factor matching exception is a commuting cubic component
with a cheaper known two-band realization. It changes the local
representation, not the frontier.

The stronger new criterion is expressed by the cross-commutation
matrix C of the first middle band and its probe images. Let A_S span
the actual axes in nonzero C rows, and H be the intersection of the
two actual middle spans. Then

```math
\mathrm{rank}\,C\le\dim(H\cap A_S).
```

In particular, if their only overlap is the center of a genuinely
matched physical pair, every such cross-commutator vanishes. This
holds for arbitrary support of that center in the other native list;
the earlier dependent-code support calculation remains a useful
component, but its residual light-coupling escape is now closed.
Commuting unmatched images are not actual-axis matching, and no heavy
phase is deleted merely because its span is invariant.

A separate [paid removal](FOUR_BAND_CENTER_TRANSPORT.md#i4-a-nonsingular-silent-block-has-a-paid-clifford-removal)
applies when the affected first-band axes silent on the second band
have a nonsingular restricted pairing. Their two actual subbands form
a full-space Clifford. Replacing it removes twice the subspace
dimension in T gates, adds quadratic Clifford synthesis cost and
preserves the literal scalar, actual inverse and every input column.
Singular blocks, including the commuting cubic guard, remain outside.

The [Pauli-domain fiber inequality](HIGHER_SUPPORT_PAIR_BOOLEANITY.md)
also rejects the stated common quantum-spectator extension of the
higher-support trial. Its two factors have a two-plane domain deficit
on an actual fiber, which prevents any full logical Pauli Lagrangian.
This concerns that attachment, not all additional quantum coupling.

The [minimal singleton-tail trial](FOUR_BAND_CENTER_TRANSPORT.md#singleton-sector-assembly)
shows why agreeing sector symplectic maps is insufficient. Its exact
Pauli corrections have quadratic labels, so the assembled actual
endpoint is non-Clifford. The remaining doubly controlled Pauli has
an explicit four-T, two-band receipt with no helper and full literal
return. Paying that residual restores the old cost; the mechanism is
therefore set aside without turning its failure into a general bound.

The next noncommuting scout needs overlap on C-active coordinates,
or a change of the original center's transport. It must supply an
actual output endpoint and complete Boolean action on every dirty and
reference input before any degree statement is inferred. The
unrestricted higher-layer and dimension-growing depth problems remain
open.

## 8. Research decision

Continue the complete source through wide address-bank readout and
the enlarged sublogarithmic-clean source window, testing the first
native operation whose cost or return has not been paid. Retain the
structured arithmetic components without making their radix conversion
a required next step. Keep the focused physical-pair scout
on actual active overlap and Boolean endpoint realizability. Stop
repackaging a mechanism when its completed bill retains the old cost.
The new zero test, exact arithmetic and eligibility theorems are
useful progress; a complete shallow source with simultaneous budgets
is still the frontier-changing missing result.
