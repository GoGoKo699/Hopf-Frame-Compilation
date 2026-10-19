# Revision: constant translations and approximate query transfer

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Constant-translation closure](CONSTANT_TRANSLATION_NATIVE_CLOSURE.md) ·
[Center saturation limits](NATIVE_CENTER_SATURATION_LIMITS.md) ·
[Actual native witness](ACTUAL_NATIVE_CENTER_WITNESS.md) ·
[Approximate query transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Previous checkpoint](SHARED_CENTER_SESSION_RESULT.md)

This continuation began at verified main `871941a` (PR #192) on
19 October 2026. It tested the two unpaid steps from the previous
revision: decoding the common monomial algebra, and compiling its upper
query counterpart at the accuracy actually needed by the frame. The
result is a new restricted degree theorem and an approximate conditional
transfer, together with exact limits on the proposed decoding mechanism.
The complete-frame frontier is unchanged.

## 1. Destination and strongest results

The destination remains optimal general depth for the prescribed complete
Hopf differential frame, with all resource and operator requirements met
on the same circuit.

| Role | Retained conclusion |
| --- | --- |
| Strongest general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$ |
| Literal domain | $`n\ge1,L\ge6`$; one clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full initialized-isometry error at most $`2^{-L}`$ |
| Selected theorem | R79 remains Result D; the selected count frontier is settled at its stated allocations |
| Unrestricted native degree results | At most one clean and two layers give XOR degree at most two; zero clean and three layers give sign/XOR bounds five/four |
| Retained components | Both fixed endpoint clean-active ranks close, including arbitrary middle pairs; fixed quadratic-control polar rank closes in the no-pair reduction |
| New lower component | Constant translations in the common monomial endpoint family force a full diagonal C3 operator, hence Boolean sign/XOR bounds three/two |
| New upper component | Constant-in-address readout accuracy suffices at fixed frame precision, under a uniform simultaneous resource hypothesis for the approximate source |
| Still open | Growing-rank one-clean three-layer degree, a faster complete query meeting every budget, and matching general depth bounds |

At fixed accuracy and square-root dirty width, the general depth interval
remains a constant lower bound versus a logarithmic upper bound. None of
the representation-specific obstructions below closes that interval.

## 2. A positive closure before asking for a decoder

The previous global matching theorem gives a common invariant isotropic
space S and, after paid normalization of the genuine shared center,
commuting monomial Clifford endpoint images. Their computational action
has translations of the form

```math
(z,y)\longmapsto(z,y+a_s+L_sz).
```

The new theorem covers the case $`L_s=0`$ for every s, without bounding
the ranks or dimension of S. It also permits an arbitrary actual commuting
C3 middle operator supported on S. Native independence makes every
endpoint axis activated by S computational-Z in the corresponding
endpoint coordinates. Those actual sublists peel as diagonal C3 factors.
The remaining Pauli transports admit a literal diagonal Clifford match;
normalizing their shared S makes the middle phase computational diagonal.
The remaining full-space two-native-band diagonal word is C3 by the
retained theorem. Restoring only diagonal C3 factors proves the claim.

This is exact degree accounting, not a new implementation at the original
depth. It uses neither a free quartic correction nor a family of separately
chosen fiber encoders. For the actual clean-derived word, restoring the
previous central cubic and quadratic signs keeps the sign bound three.
The XOR bound applies after checking these hypotheses in the actual
output-Hadamard sign representation.

## 3. What enlarging the center does and does not buy

Three genuine center constructions are now explicit: displacement
directions, endpoint active-axis sums, and whole-band axis products. The
active-axis-sum argument retains the possible fourth-root scalar in
$`Q\overline Q`$; the whole-word Galois argument retains the literal
native signs. Each uses the original global lists. Repeating displacement
normalization for the same S and encoders creates no new direction.

A general exact criterion also settles the proposed fixed Clifford
decoder. For an independent native band E and isotropic Pauli space S,
$`E^\dagger\mathrm{Alg}(S)E`$ admits a fixed Clifford decoder
to a Pauli algebra if and only if every E axis already commutes S.
Nonzero coefficients on the original and activated Pauli labels prove
this. The paid native decoder remains available with its old T cost.
This criterion concerns that encoded algebra and that decoder mechanism;
it is not an obstruction to a different compiler or to a degree theorem.

The new four-wire source realizes literal CCZ on the occupied clean
sector using three independent native bands. Both interior clean images
have active rank two, there are no middle pairs, and there is no shared
Pauli clean flip at all. Its unused clean sector has a nontrivial affine
permutation. On the occupied logical register, all three center
constructions are already saturated but a work translation still depends
on $`r+s`$. Thus the failure of constant-translation reduction occurs in
an actual native source, not merely an abstract monomial representation.
Its target has degree three and belongs to the already closed fixed-rank
case. It is not a high-degree or growing-rank counterexample.

Safe-axis derivative pruning gives a sharper quartic-plus-quotient
description using only endpoint axes that fail to commute the actual
middle operator. In the genuine quadratic normal form these are precisely
the axes failing to centralize S. The actual source retains a nonzero
residual parity, and an entire endpoint Clifford image need not remain
Clifford on passage through the middle, even though the final derivative
is just CZ. This isolates a missing cancellation between distinct
endpoints; it does not rule out another paid derivative implementation.

## 4. The first lower inference still missing

The genuinely open monomial case has some $`L_s\ne0`$. After the
displacement center is normalized, the quadratic-controlled Pauli
contributes a central cubic defect $`b(z)\ell_u(z)`$ to each native
matching equation. Its cancellation must come from the actual middle
band's phase defect under the induced action. Removing a central phase
from that band does not change the defect, because the probes fix the
center pointwise.

The next lower question is therefore whether these native phase-defect
identities, simultaneously for all computational probes and both
endpoints, bound the completed Boolean phase independently of ranks and
width. Support invariance, center saturation and a formal Clifford
involution do not supply that inference. The constant-translation theorem
shows exactly one regime where actual sublist peeling pays it.

## 5. The upper accuracy requirement is now paid through the transfer

A zero-clean approximate Walsh dual readout with full operator error tau,
including arbitrary dirty inputs and its actual inverse, gives a complete
indicator with error at most $`6\tau`$. The R91 word-query assembly
uses at most four such indicators, so its error is at most $`24\tau`$,
independently of the table length and batching. A diagonal table-phase
gauge cancels in the wrapper.

The conditional frame construction needs only a constant number K of
completed long-query occurrences. All short tables and routers remain
exact. Choosing

```math
\tau=\frac{2^{-L}}{96K}
```

adds at most $`2^{-L}/4`$ to the full error and fits the retained
precision slack. Thus for each fixed L, the source needs constant accuracy
as the address grows. A source compiler must meet its count, width,
Clifford and depth budgets uniformly at that accuracy and on every shorter
address. Any precision cost inside its native implementation remains its
own obligation. No faster source is constructed here.

The exact standalone readout count obstruction remains valid for arbitrary
Clifford interlayers and returned helpers, but does not extend to ordinary
approximation. A separate robust theorem excludes the familiar strategy
of a separately synthesized Majorana-preserving Gaussian conjugator:
small constant error still forces logarithmic conjugator depth and
superlinear T-count. That hypothesis concerns the native decomposition,
not merely a Gaussian endpoint. It does not exclude an unrestricted
approximate readout or a jointly synthesized complete indicator.

## 6. Decision and evidence

Continue with the simultaneous phase-defect question on the lower side,
and unrestricted approximate readout or joint complete-query synthesis on
the upper side. Further center normalization by itself retains the old
cost and does not justify a new frontier claim. Requiring vanishing source
error with address size would also impose an unnecessary condition on the
upper route.

The witness fixtures check bounded exact native words, all columns,
literal phases, actual inverses and their ledgers. The approximate-source
fixture checks a bounded native perturbation with a full-operator error
certificate and the pair-plane geometry of the scoped obstruction.
Arbitrary-rank closure, decoder
rigidity, shared-center identities, resource splices and the robust
asymptotic bound are analytic proofs. Every claim retains its stated
scope; manuscript drafting remains separate.
