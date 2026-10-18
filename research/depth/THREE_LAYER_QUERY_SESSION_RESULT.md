# Three-layer queries: execution outcome

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Authoritative proof and native fixtures](READ_ONLY_THREE_LAYER_RIGIDITY.md) ·
[Previous round](RETURNED_QUERY_SESSION_RESULT.md)

This continuation began from verified main
`3913263b862a01a56ec38405dc90abab8be9e92a`, the merge of PR #166,
on `research/three-layer-query-round`. It resolves the previous round's
selected three-layer question within the preserved-address interface.
The complete-frame upper bounds and selected manuscript claims are unchanged.

## 1. Compact revision and the missing invariant

The count package is settled under its stated hypotheses. R96 remains
the strongest unconditional general research upper law, for every
$`n\ge1,L\ge6`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It uses one clean qubit, dirty width $`b\ge17(L+n+7)`$ for real
frames or $`b\ge17(L+n+8)`$ for phase-dressed complex magnitude
frames, and initialized-isometry error $`2^{-L}`$. Every logical
column, literal sign and phase, singular and zero coefficient, arbitrary
dirty/reference input, actual inverse, leakage and complete workspace
return remains part of the contract. R79 remains selected Result D.
R97 gives a scoped four-layer obstruction; unrestricted one-clean
four-layer rigidity and dimension-growing depth lower bounds remain open.

Previously open general-magnitude recovery, recursive terminal overhead
and precision-window restrictions have been resolved. At fixed accuracy
and square-root width, optimal general depth still lies between a constant
lower bound and an $`O_L(\log(n+2))`$ upper bound. Conditional R98
would improve the complete circuit through a faster exact **zero-clean,
multi-output** dirty query with R91's full resource envelope. Its current
logarithmic address-indicator depth, rather than its already constant
eight-layer bilinear echo, is the remaining cost on that route.

This round addressed the first missing native invariant beyond the
two-layer proof: can three physical T layers produce unbounded Boolean
query degree with one clean qubit and arbitrarily many dirty helpers?
Completed Clifford interlayers must preserve the address and T gates
act only on work. The intervening dressed Pauli families can leave the
initialized subspace, so logical compression or iteration of the old
three-band argument would be invalid. Unlike the stopped Gaussian
route, the present proof allows non-Clifford work interference.

## 2. Supported results and what changed

The single proof home is
[three-layer rigidity](READ_ONLY_THREE_LAYER_RIGIDITY.md). At full
operator or initialized-isometry error
$`\eta\lt\rho=\sqrt2/6`$, it proves:

| Simultaneous initialized work | XOR-query conclusion | Literal Boolean scalar-phase conclusion |
|---|---|---|
| Zero clean; arbitrary finite dirty width | Every output is affine | Degree at most two |
| One clean; arbitrary finite dirty width | Degree at most four | Degree at most five |

These are necessary conditions independent of T-count and elementary
Clifford count. They do not assert attainability for all tables below
the degree bound. In particular, five-bit AND cannot be implemented or
approximated within the stated radius by a one-clean three-layer query
in this interface. The previously selected unbounded-degree escape is
therefore closed.

The new identity factors the four address-control bands into endpoint
Paulis and two exact projective Pauli representations. An actual-inverse
comparison produces a Pauli-label cocycle with error at most $`6\eta`$.
With zero clean qubits, the distinct physical T axes make that cocycle
additive. A generic affine symplectic action would not suffice: a small
nonlinear CNOT cocycle is retained as a falsifier of that weaker argument.

For one clean qubit, the retained three-layer gap first freezes only
**terminal** clean return. The two dressed representations transport the
same signed Pauli orbit of intermediate codes. Their product commutes
with that whole orbit. Anticommuting orbit members extend the occupied
Pauli approximation to full work; commuting members allow comparisons
on nonzero code intersections. Opposite sectors are excluded explicitly.
This proves degree at most two on each endpoint-syndrome coset, and the
two endpoint syndrome bits give degree at most four globally. No
intermediate return, bounded dirty rank or free clean-removal step is assumed.

## 3. Constructive tests and the phase cost that survives

Exact native Toffoli queries test the invariant at its first nonlinear
boundary. In both rows the address a,b is read-only and the target t is
arbitrary; the dirty column counts only additional helper d.

| Literal query | Clean | Dirty helpers | T-count | T-depth | Elementary Cliffords | Error/leakage |
|---|---:|---:|---:|---:|---:|---:|
| $`t\mapsto t\oplus ab`$ | 0 | 1 | 8 | 4 | 18 | 0 |
| $`t\mapsto t\oplus ab`$ | 1 | 1 | 8 | 3 | 28 | 0 |

The proof home gives the eight-term integer phase identity, every actual
CNOT basis change and its inverse, the complete complementary-sector
action, and the work-return proof. All phases are literal. The four-layer
zero-clean construction is depth-optimal **within this interface**, by
the new theorem, including the stated positive error radius.

A different native mechanism gives a two-layer candidate using one clean
qubit, no additional dirty helper, four T gates and 16 elementary
Cliffords. It returns its clean qubit exactly but implements

```math
(-i)^{ab}X_t^{ab},
```

not the required literal query. Its missing correction is the
non-Clifford address phase $`i^{ab}`$, a controlled-S. Thus the apparent
saving retains a concrete phase cost; the phase cannot be hidden in a
Clifford band or discarded as a common global phase. The exact full-space
fixture records this useful interference mechanism without promoting it
to a literal two-layer construction.

The phase diagnosis led to a second, exact-only result in the
[two-layer proof home](READ_ONLY_TWO_LAYER_RIGIDITY.md#5-exact-one-clean-queries-are-affine):
every literal one-clean two-layer XOR query in the same interface is
affine. The old robust degree/polar-rank bound reduces a hypothetical
nonlinear query to two endpoint-syndrome cases. In the independent case
the occupied block is scalar. In the common-syndrome case a genuine
one-T-conjugated Pauli with distinct logical sector labels forces those
labels to anticommute. The remaining occupied block is then
anti-Hermitian up to a **real** sign, and cannot equal the literal X
target. Keeping the address Clifford's mixed phase real is essential;
discarding it would miss precisely the controlled-S debt above.

This makes the one-clean three-layer AND fixture **exactly** optimal
within the interface. No positive error radius is asserted for this
stronger two-layer affine statement; the earlier robust degree-two and
polar-rank-two theorem remains the approximate result. The two-bit
restriction substitutes linear address controls on the same work
register and supplies no extra initialized address helpers.

The retained unrestricted two-layer dirty Toffoli is also a decisive
scope test: it applies T gates to transformed address data and therefore
lies outside the new theorem. Neither choosing another native
representation nor returning the address only at the end preserves the
intermediate-boundary hypothesis automatically.

## 4. Precise remaining step and next question

**Subsequent resolution:** the
[four-layer follow-up](FOUR_LAYER_QUERY_SESSION_RESULT.md) proves
zero-clean XOR degree at most two at error below $`\sqrt2/12`$.
The four-layer test below is therefore resolved. The uniform
degree-growth question remains current, with five layers its first
unclosed test; the following passage records this round's handoff.

The most useful next question is a **degree-growth law for exact
zero-clean read-only queries as the number of physical T layers grows**:
is there an absolute C such that every such d-layer XOR query satisfies
$`\deg f\le2^{Cd}`$, independent of dirty width, or is there a literal
native family violating every width-independent degree bound at some
fixed d? This is a question, not a consequence of the present proof.
Intermediate growth rates are also possible; failure of the exponential
bound alone would not imply the stronger fixed-depth escape.

The concrete first case is four layers. Its five-band normal form has
three dressed Pauli representations. After the outer Clifford family
is removed, the product of the remaining two is generally **not** a
projective representation. Consequently the present exact cocycle and
shared Pauli-orbit identity do not iterate. The missing step is to
control that product's commutator defect, perhaps through a second
finite-difference invariant, without replacing dependent transported
axes by independent physical axes. Start with zero clean qubits to
isolate this issue from the separate unresolved four-layer clean-return
gap. A counterexample must give a literal circuit, actual inverses,
full arbitrary-work return and all native resource charges.

An exponential degree bound would force $`\Omega(\log a)`$ depth
for an a-address AND query and hence for worst-case tables in R98's
existing interface. It would explain why a sublogarithmic complete-frame
improvement must change that interface or use a different construction.
A fixed-depth high-degree mechanism would instead identify interference
to test against R91's multi-output count, Clifford and simultaneous-width
envelope; high degree by itself does not compile arbitrary tables.

For a general frame lower bound, two further obligations remain: a
dimension-growing bound robust at a fixed positive error, and an
invariant valid under unrestricted address mixing. The
[two-call query/frame reduction](BOOLEAN_QUERY_FRAME_REDUCTION.md)
requires that broader native premise. It cannot import the current
read-only obstruction. R98 remains conditional, R96 remains unchanged,
and optimal general depth remains the continuing destination.

## 5. Evidence and integration

Seventeen bounded tests cover the literal four-band identity and actual
inverse with non-Pauli zero-address action, independent-axis additivity
and its generic symplectic falsifier, stabilizer-intersection gaps,
opposite sectors, the anticommuting-orbit extension, all 65,536
four-variable truth tables for the derivative lemma, the scalar-to-query
conversion, and complete native matrices and resource ledgers for the
three constructions above. Four of these checks target the exact
two-layer strengthening: physical flat-support commutation, an abstract
three-band literal AND outside the one-T restriction, imaginary phase
retention and cancellation of the actual complementary-sector unitary.
These finite tests support the stated
identities and boundaries; the width-independent theorems are analytic.
Independent reviews checked both clean allocations and all phase,
leakage, intersection and coset steps.

The repository's retained validation and presentation gates are recorded
with the merge. No selected manuscript claim, frame compiler or immutable
scientific receipt is changed.
