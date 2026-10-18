# Unbounded paired cores and controlled-Pauli transport

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Isolated cores](ISOLATED_NATIVE_CORE_STRUCTURE.md) ·
[Controlled-Pauli transport](CONTROLLED_PAULI_TRANSPORT.md) ·
[Constructive boundary](NONCOMMUTATIVE_CYCLE_CHARACTER.md) ·
[Previous checkpoint](TWO_OVERLAP_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`ad67843633e6a3129e915b330035077bf43a1857` (PR #178).
Its question is whether growing control rank is an actual transport
obstruction, and whether block inversion can remove the cycle source's
sequential native cost. This is a research checkpoint; final manuscript
drafting remains separate.

## 1. Compact revision

The count package is settled under its stated hypotheses. For integers
$`n\ge1,L\ge6`$ and $`N=2^n`$, R96 remains the strongest unconditional
complete-frame research upper law:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It uses one clean qubit and dirty width $`b\ge17(L+n+7)`$ for real
frames or $`b\ge17(L+n+8)`$ for the prescribed complex magnitude
completion, with full initialized-isometry error at most $`2^{-L}`$.
Every logical column, literal phase and sign, singular and zero
coefficient, arbitrary dirty input and external reference, leakage,
actual inverse and complete workspace return remains included. Counts,
T-depth, precision and simultaneous workspace are separate budgets;
every source, query, control, recovery and inverse is charged.
R79 remains selected Result D; R97 is a scoped four-layer obstruction;
R98 is conditional on a faster exact multi-output native read-only query.

The exact general question is the optimal worst-case T-depth of the
complete prescribed frame under those simultaneous budgets. At fixed
accuracy and square-root width, does some family require growing depth,
or can the additive logarithmic-in-n upper bound be reduced?
Unrestricted one-clean four-layer rigidity and dimension-growing
unrestricted lower bounds remain open. This round does not improve R96.

Already resolved questions include the one-clean count frontier, the
all-precision general upper law, zero-clean read-only five-layer degree
14, six-layer Pauli-comparison degree 51, diagonal-C3 transport 78/26,
and exact native source closure through overlap two. The bounded-core
hybrid theorem also separates wide diagonal coupling from bounded quantum
complexity. Old next-step passages are historical evidence, not reasons
to reopen these questions.

## 2. Why this round changed the source question

The overlap-two proof used the small number of controls when bounding
sign functions. Merely increasing that finite cutoff would not supply a
native depth invariant: overlap can grow in a shallow physical source.
The [isolated-core proof](ISOLATED_NATIVE_CORE_STRUCTURE.md) instead
identifies a full-space invariant that works with any number of controls.

After native pair extraction, each active target has a fixed Pauli axis.
The Clifford images of all control-X operators commute and are
involutions. Their symplectic actions force every control derivative of
that target-axis sign to be constant. Thus the sign is affine, without
assuming a bound on the number of controls. A complete commuting Pauli
frame has Pauli images, giving diagonal C3 coordinates. The exact paired
phase formula separately proves that the extracted source is C3 before
remaining native inner axes are restored.

The resulting theorem covers arbitrarily many isolated paired cores.
A checkable sufficient condition is that zero-character B axes are good
and every C axis anticommutes with at most two B axes. In particular,
no zero characters and nonzero activation classes of size at most three
suffice. The total first/middle overlap is unrestricted. The existing
six-layer degree bounds 78/26 therefore apply below full error
$`\sqrt2/8`$ in these source families.

A further full-space exclusion handles a single nontrivial C orbit when
it contains all actual B/C noncommutation. Its incidence row must have
weight two, irrespective of total overlap. If the row had even weight
k at least four, the isolated control image would force a target Pauli,
while pulling it back through the middle band leaves a nonzero
rank-$`k-2`$ Clifford action unless it is the fixed axis. That remaining
axis gives an impossible odd relation in a nonzero activation class.
Good zero-character axes are allowed; additional bad zero axes are not.
Thus a single wider orbit is also closed, while coupled orbits remain.

A literal three-core fixture already has overlap three and phases coupling
a target to two other controls. Its fair simplified baseline remains one
T layer. The theorem is a source certificate, not a claim that the
displayed three-layer source is depth optimal.

## 3. A failed diagonalization extension leads to direct transport

Allowing targets with constant binary Clifford action invalidates the
broader proposed diagonalization certificate. The known seven-qubit
Gottesman–Mochon source, recorded in Beigi–Shor, Section 1.1, has exact
coordinates

```math
F_z=X^{f(z)}Z^{g(z)},\qquad
f=(z_1z_2,z_1z_3,z_1z_4),\qquad
g=(z_3z_4,z_2z_4,z_2z_3).
```

Its Pauli-preserving domain is only the four control-Z directions; it
is C3 but not semi-Clifford, and its inverse is not C3. The source and
its actual inverse have a paid schedule with 42 T calls, 66 elementary
Cliffords, T-depth at most 30, seven arbitrary data wires and no helpers.
This is an attributed scope fixture, not a new native three-band source
or a returned query. Its literal square is a quartic control phase.

The successful reformulation is [direct controlled-Pauli transport](CONTROLLED_PAULI_TRANSPORT.md).
Assume the specified source has controlled-Pauli blocks, is exactly C3,
and has the theorem's canonical exact Pauli domain, or the stated
explicit compatibility with the actual inner Clifford family. The target
and control widths may both grow without bound.

C3 forces quadratic block labels. It also forces the raw scalar phase
into a modulo-eight form with linear, twice-quadratic, four-times-cubic
and four-times-quartic terms. The quartic term is necessary even after
a fixed target Hadamard change of coordinates. These terms still give
bounded binary phase degrees. Exact C3 endpoint rounding and at most
four fixed control blocks per output coordinate then yield

```math
\deg S_K\le58,\qquad \deg\kappa_K(P)\le174,
\qquad \eta\lt\sqrt2/8.
```

The fixed-probe sign includes all modulo-four carries from the Hermitian
Pauli convention. The same bounds give 58 for visible outer target shear
and 174 for general target degree. They are independent of source width
and Pauli deficiency. Tensor powers of the known source satisfy the
hypotheses, so diagonalization and bounded residual size are not necessary
for this transport argument.

The Pauli-domain or explicit compatibility premise remains essential to
the proof. It is not inferred from complete return or arbitrary C3
membership. No coordinate transformation is supplied as a free physical
gate, and no inverse-C3 assumption is used.

## 4. Constructive block inversion: the exact missing return

The [cycle proof](NONCOMMUTATIVE_CYCLE_CHARACTER.md) now gives a block
version of its commuting pivots. For one contiguous middle block A,
write $`R=(I+A)^{-1}`$. Pivoting through it sends its incident blocks
and cross block by

```math
(A,B,C,D)\longmapsto(R-I,BR,RC,D+BRC).
```

This is an exact all-input involution. A fully charged constant-depth
product XOR into a separate arbitrary array follows by a four-call
dirty echo. The unresolved operation is the in-place incident-block
update, not that product XOR.

For an invertible linear block operator M, three alternating shears and
a swap implement $`(B,E)\mapsto(MB,M^{-1}E)`$. The helper counter-action
is literal and cannot be discarded for arbitrary E. The inverse shear
also needs the actual original coefficient blocks; restoring and
recomputing them adds two child-inversion positions. Including their
initial inversion gives three sequential positions before helper repair,
so this schedule does not justify the proposed logarithmic recurrence.

A scoped determinant invariant sharpens the diagnosis. With frozen
coefficients, every uniform polynomial whole-array shear and swap has
determinant one over the faithful ring generated by M. A fully returned
action $`\mathrm{diag}(M,I,\ldots,I)`$ has determinant M, excluding it
when M is nonidentity. This remains true even if each permitted shear is
granted without cost. Nonuniform coordinate operations, coefficient
changes and Hadamards lie outside the invariant; it is not an unrestricted
native lower bound. The retained cubic-count, linear-depth cycle schedule
and the faster conjunction baselines remain unchanged.

## 5. Most useful resumption question

**Does full C3 membership force overlapping orbit-incidence components of
an actual three-native-band source into a class with width-independent
signed transport, or can such components realize a native high-degree
transport family?**

Within a nonzero repeated activation class, the even C-orbit incidence
rows form a full-row-rank self-orthogonal binary code. The first coupled
test has rows 111100 and 001111. They are legal independent isotropic
B/C axis data; activation matching alone does not make their complete
source C3. Weighted overlaps seen by an isolated A probe need not vanish,
so the isolated-target argument cannot simply be repeated.

A useful next proof must use the complete source's Clifford involutions,
including their signs, on those coupled components. The new controlled-Pauli
result offers a transport route even if diagonalization fails. A true
countermechanism must retain a native realization, the actual inner
family and full return; an abstract high-degree C3 action is insufficient.
General zero-character interactions remain a separate source issue.

A separate possible extension is to remove the controlled-Pauli theorem's
canonical-domain premise. The missing lemma is a fixed Clifford
normalization of the entire signed Pauli domain that preserves scalar-
Pauli control blocks on both sides. Straightening selected domain
directions alone is insufficient; no such general normalization is
claimed in this round.

Success would remove a structural restriction in the native
degree-versus-depth program. It would still leave non-C3 prefixes,
one-clean work, address mixing and a robust recurrence at growing depth
before an unrestricted lower bound transfers to complete frames.
On the constructive side, resume with a returned in-place block map that
leaves the frozen-coefficient whole-array interface, and charge all
sources, controls, inverses and overlapping workspace before comparing
against the logarithmic-depth conjunction baseline or the R98 interface.

## 6. Evidence and integration

The authoritative homes contain analytic arbitrary-width proofs.
Independent audits checked the commuting-control sign invariant,
C3-derived phase budget, fixed-block recovery and constructive determinant
scope. Bounded checks retain literal phases, actual inverse words,
all-input identities and returned private helpers. They are mechanism
checks, not substitutes for the proofs. Validation and rendering receipts
belong to the merging pull request. No manuscript claim selection changes.
