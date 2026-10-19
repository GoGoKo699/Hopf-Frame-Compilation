# Revision: where the remaining native cubic defect can live

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Middle-defect closure](MIDDLE_DEFECT_SUPPORT_CLOSURE.md) ·
[Approximate indicator limits](APPROXIMATE_INDICATOR_FILTER_LIMITS.md) ·
[Previous checkpoint](CENTRAL_COCYCLE_SESSION_RESULT.md)

This continuation began at verified main `6876094` (PR #193) on
19 October 2026. It revisited the nonempty-middle gap before asking for
another decoder. The result is a new degree closure for bounded
center–middle overlap, and a sharper description of the remaining
unbounded case. On the upper side, two approximate predicate mechanisms
are assessed against the complete dirty indicator. Neither supplies a
new complete-frame circuit, and the general resource frontier is unchanged.

## 1. Destination, theorems, components, and mechanism limits

| Role | Current conclusion |
| --- | --- |
| Destination | Optimal general depth for the prescribed complete frame, with every count, width and full-operator requirement met on the same circuit |
| Strongest general upper theorem | R96 gives $`D_T=O(NL/b^2+L+\log_2(n+2))`$, optimal $`T=O(\sqrt{NL}+L+NL/b)`$, and $`G=O(NL)`$ |
| Literal domain | $`n\ge1,L\ge6`$; one clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full initialized-isometry error at most $`2^{-L}`$ |
| Selected result | R79 remains Result D; the selected count frontier remains settled at its stated allocations |
| Unrestricted native components | At most one clean and two layers give XOR degree at most two; zero clean and three layers give sign/XOR bounds five/four |
| Retained one-clean components | Both fixed endpoint clean-active ranks close, including arbitrary middle pairs; fixed quadratic-control polar rank closes in the actual no-pair reduction |
| Stronger retained subcase | Constant translations on the entire common isotropic support give sign/XOR bounds three/two; an empty middle band closes uniformly |
| New component | Constant control-phase commutators give sign/XOR bounds five/four even with a nonempty middle; bounded exceptional center–middle overlap closes uniformly in all other ranks |
| Method limits | A fixed Clifford decoder, additive polynomial phases, and short hash summaries have specific obstructions; these do not exclude unrestricted native synthesis |
| Remaining questions | Growing-rank one-clean three-layer degree, faster complete queries satisfying every budget, and matching general depth bounds |

At fixed accuracy and square-root dirty width, the complete-frame depth
interval remains a constant lower bound versus a logarithmic upper bound.
The new closure reduces the possible lower-side counterexamples; it does
not close that interval. The retained exact-rounding floor and all
initialized-isometry hypotheses remain unchanged.

## 2. Peel only what the control needs

The actual no-pair reduction has three independent native bands and
one retained insertion $`F=H^{b(z)}`$, where R is a genuine shared
computational center and b is quadratic. For every computational probe,
global matching gives a central affine parity p_u and the literal identity

```math
K_u=D_u(-1)^{b(z)p_u(z)},
\qquad D_u=E_BJ_uE_B^\dagger J_u^\dagger,
```

with K_u Clifford. The previous proposed decoder tried to normalize the
whole common support, including the middle band. That was more than the
degree argument needs.

If all p_u are constant, only the endpoint translations of
$`S_F=R+\langle h\rangle`$ must be constant. Native independence
identifies actual diagonal endpoint sublists active for this support.
Peeling them and matching the residual Pauli transports makes F a
diagonal C3 input factor. The remaining full-space diagonal word has
three actual independent native bands, so the retained full complex
three-layer theorem gives ordinary phase degree at most five.
Restoring the diagonal C3 factors gives sign degree at most five.

This proof retains the native middle band even if its endpoint images
still have center-dependent translations. Identity-mode Hadamards are
used only through an equality for the whole residual word; they are not
undone separately on its factors. Actual affine input changes and literal
scalars remain accounted for. The earlier cubic theorem still gives the
sharper result when its stronger whole-support hypothesis holds.

## 3. Localize the unpaid cubic cancellation

Remove only actual middle-axis orbits that are permuted by every J_u.
Their signed native phase defects are Clifford. The remaining actual
coordinate sublist, denoted B_exc, is invariant, and every cubic middle
defect is supported in its span. Define

```math
W=R\cap\mathrm{span}(B_{\mathrm{exc}}),
\qquad \beta=b_2,\qquad \ell_u=(p_u)_1.
```

An exact diagonal Clifford square root of the middle defect squared
leaves a Boolean C3 residual. Comparing its cubic symbol in the matching
identity proves

```math
\beta\wedge\ell_u\in\Lambda^3W
\quad\hbox{for every computational probe }u.
```

This uses Boolean cubic symbols, not an unsupported coordinate invariance
for raw weighted modulo-eight phase coefficients. The removed factors
are actual axis orbits; no new orthogonal basis is declared native.

At polar rank at least four, any nonzero ell_u must lie in W, and beta
must lie in $`\Lambda^2W+\ell_u\wedge R`$. Hence control rank
greater than $`\dim W+2`$ forces every p_u to be constant and gives
the degree-five closure. If $`\dim W=d`$ is fixed, the complementary
low-rank controls have insertion support at most d+3 and use the retained
finite-insertion theorem. Thus fixed d closes uniformly in middle width,
endpoint ranks and control dimension.

The support statement also closes some unbounded-W cases at control rank
at least four: a nonzero quadratic component on R/W, or mixed W-to-R/W
rank at least two, forces all p_u constant. A genuinely remaining
high-rank candidate must place its quadratic control wholly inside W,
or have only one mixed direction and one possible nonconstant parity
line. This is a necessary condition in the actual normal form, not a
construction or a general degree theorem.

## 4. Approximation must be joint across the dirty bank

The retained approximate transfer shows that fixed frame precision needs
only constant-in-address accuracy for a completed source. That permits a
change of mechanism, but it does not permit replacing a complete query
by unrelated bounded-error equality tests.

The first tested mechanism uses a degree-below-p real polynomial of
Hamming distance as an additive phase coefficient for each dirty output
bit. Its complete operator error is at least $`\sqrt2`$. The bound is
sharp even though every individual coefficient can have error only
$`2^{-p}`$: arbitrary dirty subsets expose their accumulated phase.
This is a phase-program obstruction, not an unrestricted approximate
query lower bound.

A changed candidate uses coherent random linear hashes. A short hash
can approximate one equality-controlled phase on a uniformly initialized
seed. The first unpaid native operation is that initialization: an
arbitrary dirty seed can select a bad hash and give operator error two.
There is also a separate whole-bank failure even with a clean coherent
seed. A rank-k hash with k<p selects an even-size bucket, so the all-one
phase-basis output gives phase +1 while the complete indicator requires
phase -1. Nonlinear postprocessing of the same insufficient parity
summaries does not restore the missing information.

These failures leave general interference and richer direct access to the
dirty output outside their scope. No independently accurate predicate is
credited as a jointly accurate, work-returning indicator.

## 5. Research decision and evidence

The next lower question is whether the original clean-derived active and
silent lists allow the remaining exceptional overlap to carry the
required simultaneous cubic cancellations. An arbitrary isotropic model
or a complex diagonal completion does not settle that actual-source
question. Bounding a raw center dimension before removing permuted axes
would also miss the parameter that the new theorem actually uses.
An explicit independent native-list fixture permits a nonconstant
displacement inside the exceptional overlap and satisfies the clean
Gram identity and literal total-product relation, but fails the Boolean
cubic-support identity. Thus native independence, Gram matching and
total products alone do not finish the argument;
the next step must retain the literal clean-transition phase as well.

The next upper question is a jointly accurate native indicator mechanism
with linear T-count and dirty width, quadratic Clifford count, zero
initialized work, and depth below the retained logarithmic address
construction. Its actual inverse and full dirty/reference return must
be supplied at the same budgets. The hash and polynomial limits tell us
which interfaces have lost information; they do not justify abandoning
the full query objective.

The bounded-overlap closure and approximation obstructions have analytic
proofs. Bounded exact tests check exterior-support edge cases, native
orbit identities, and complete small phase/seed witnesses. They are
checks of those mechanisms, not empirical proofs of unrestricted lower
bounds. Final manuscript drafting remains separate.
