# Zero-character source certificates and an intrinsic transport test

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Native source proof](NATIVE_BAND_SOURCE_STRUCTURE.md) ·
[Paired-orbit proof](ISOLATED_NATIVE_CORE_STRUCTURE.md) ·
[Intrinsic criterion](CONTROLLED_PAULI_TRANSPORT.md) ·
[Query cost audit](GAUSS_SUM_QUERY_ENCODING.md) ·
[Previous checkpoint](COUPLED_INCIDENCE_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`5e53931827efef4b16ba395d430dd561fe70bc1a` (PR #180).
It addresses the bad zero-character source interactions left by that
checkpoint, and the native read-only conversion cost left by the complete
logarithmic-depth inverse query. Proofs and bounded evidence live in the
linked authoritative homes. This is a research checkpoint; final manuscript
drafting remains separate.

## 1. Compact revision

The count package is settled under its stated hypotheses. R96 remains
the strongest unconditional complete-frame research upper law. For every
integer $`n\ge1,L\ge6`$, with $`N=2^n`$,

```math
\begin{aligned}
T&=O(\sqrt{NL}+L+NL/b),&G&=O(NL),\\
D_T&=O(NL/b^2+L+\log_2(n+2)).
\end{aligned}
```

The same circuit uses one clean qubit, dirty width
$`b\ge17(L+n+7)`$ for real frames or $`b\ge17(L+n+8)`$ for the
prescribed complex magnitude frames, and full initialized-isometry error
at most $`2^{-L}`$. Every logical column, literal sign and phase,
singular and zero coefficient, arbitrary dirty input and external
reference, actual inverse, leakage and complete workspace return remains
included. T-count, elementary Clifford count, T-depth, precision and
simultaneous workspace are separate budgets. Every source, query,
control, recovery and inverse is charged.

The exact general question remains the optimal worst-case T-depth of
this prescribed complete frame under those simultaneous budgets. At
fixed accuracy and square-root width, can the additive logarithmic-in-n
upper term be improved, or must some family have growing depth? R79
remains selected Result D; R97 is scoped; R98 remains conditional.
Unrestricted one-clean four-layer rigidity and dimension-growing
unrestricted lower bounds remain open. No complete-frame upper law is
improved in this round.

Previously resolved questions include one-clean optimal count, the
all-precision upper law, zero-clean read-only five-layer degree 14,
six-layer Pauli-comparison degree 51, diagonal-C3 transport 78/26,
source closure through overlap two, arbitrary coupled nonzero-character
incidence when zero-character axes are good, signed normalization for
scalar-Pauli sources, and the complete logarithmic-depth structured
inverse query. Historical next-step passages are not current work orders.

The useful next source question is now an explicit **coordinate-independent
tableau identity**, stated in section 5. It is weaker than demanding
diagonalization and sufficient for the retained 174/58 transport theorem.
The constructive audit in section 4 explains why the new inverse query
does not yet improve R98, despite an existing arbitrary-table representation.

## 2. Two complementary source certificates

The source is the specified exact full-space C3 product
$`F=E_CE_BE_A`$ of three independent isotropic physical axis lists,
with their actual signed T/T-dagger phases. C3 membership is an extra
premise, never inferred from return.

**First certificate:** if every middle B axis commutes with all of A
or with all of C, F is diagonal C3 in fixed Clifford coordinates.
This includes arbitrarily many bad zero-character B axes, which lie in
the A span, and arbitrary good B axes outside it.

Partition B into bad axes and those commuting with C. The two resulting
commuting phase sources are C3 even when their merged lists have
relations. An isolating A probe certifies that each actual A axis
interacting with the good B subset has a Pauli image under F. A signed
stripping lemma then removes precisely those A rotations while preserving
C3. The remaining B/A source commutes, so the existing native left
extension applies; native right extension restores the removed rotations.
No non-Clifford factor is cancelled without its Pauli-image certificate,
and no inverse-C3 premise is used.

The same derivative gives a general preprocessing reduction without
the commutation-cover hypothesis. Remove every A axis interacting with
a good middle axis by certified stripping; the reduced source then
fixes those good middle axes, which can also be stripped. This yields
an exact C3 source $`W=E_CE_{B_0}E_{A_0}`$ with all remaining middle
axes bad, and the literal factorization $`F=E_{B_g}WE_{A_1}`$.
F is diagonal C3 if and only if W is. The companion domain-axis phase
closure in the controlled-Pauli home proves the analogous equivalence
for scalar-Pauli form. All removed rotations remain charged in the
original circuit. Activation must be recomputed for the reduced lists.

A bounded five-wire fixture has four composite bad zero-character axes,
none duplicating an A axis, plus a good B axis outside the A span and an
A axis requiring certified stripping. Thus duplicate-only cancellation
would miss the mechanism. Its literal inverse and complete source
conjugations are checked on all columns.

**Second certificate:** if every singleton C activation orbit commutes
with all B axes, then every zero-character B axis is good. The previous
arbitrary-overlap theorem therefore applies, again giving diagonal C3.
This excludes zero-axis incidence confined to paired C orbits.

The full incidence matrix now has extra zero-character columns. Its
old diagonal Gram bit can vanish, so the former nonvanishing argument
cannot simply be copied. Restricted activation characters instead prove
that every active B-side label is nonzero. The full involution defect
still has twice their span rank. Native A commutation with each zero B
axis supplies the additional equation that kills every paired zero
incidence after the nonzero classes localize.

Each certificate gives the existing six-layer 78/26 degree bounds below
full error $`\sqrt2/8`$. Neither changes the physical circuit ledger.
Together they imply that an uncovered specified source must contain both
a middle axis interacting with both outer bands and a bad singleton C
axis. They do not prove that these features belong to one independently
factorable component; native A coupling and overlapping radicals cannot
be discarded.

## 3. Scalar-Pauli blocks have an intrinsic test

For an exact C3 source F, let E be its exact Pauli-preserving input
domain and $`E'`$ its output label image. For every Pauli label p, let
$`S_p`$ be the symplectic action of the Clifford
$`F\tau(p)F^\dagger`$, and put $`N_p=S_p-I`$ over the binary field.
The [controlled-Pauli proof](CONTROLLED_PAULI_TRANSPORT.md) establishes
the equivalence

```math
\begin{aligned}
&F\text{ has scalar-Pauli blocks in fixed left/right Clifford coordinates}\\
&\quad\Longleftrightarrow\quad
\mathrm{im}(N_p)\subseteq\mathrm{rad}(E')
       \quad\text{for every }p\in E^\perp\\
&\quad\Longleftrightarrow\quad
N_pN_q=0\quad\text{for every }p\in E^\perp\text{ and every }q.
\end{aligned}
```

The proof first gives a Pauli-support coset criterion for arbitrary
unitaries. Signed symplectic coordinates match the exact input/output
domains and separate identity spectators and control Z axes. The support
condition then forces every target block to have one common Clifford
action, leaving only scalar Paulis after removing that fixed action.
For C3 sources, Clifford support has direction $`\mathrm{im}(S_p-I)`$,
and the common fixed space of all $`S_q`$ is exactly $`E'`$.
These two facts give the displayed test. Signs are retained in the
normalization; the binary test is not substituted for the literal operator.

Multiplying a scalar-Pauli source by commuting phase gadgets on its
exact domain axes preserves this class, in both directions. After
signed normalization those axes act only on identity spectators and
control Z labels. Diagonalizing their commuting spectator components
and promoting those coordinates to controls proves the claim with all
phases retained. This is the closure used by the general good-middle
reduction; it requires no extra wire or physical synthesis claim.

For an exact C3 source there is an equivalent inverse test: its actual
inverse must send every output Pauli in $`E'^\perp`$ to a Clifford.
Only a basis of $`E'^\perp/\mathrm{rad}(E')`$ needs checking.
The inverse need not be globally C3. Any inverse Pauli image already
known to be Clifford has deviation contained in the input-domain
radical, by an exact commutator argument. This proves sufficiency via
the support criterion for the inverse; necessity follows directly in
scalar-Pauli coordinates. Known singleton-C inverse images are Clifford,
but they have not been shown to span the required residual quotient.

The known He–Robitaille–Tan four-generator exterior-algebra permutation
fails this invariant. A symbolic Pauli image sends a deviation outside
the radical of its exact output domain. Therefore no alternative fixed
Clifford coordinates put that source in scalar-Pauli form. This is a
scope check using an attributed C3 source, not a native-three-band
counterexample or a returned high-degree query. The smaller three-generator
case is Clifford-equivalent to the retained Gottesman–Mochon example.

## 4. The constructive saving does not survive generic conversion

An arbitrary Q-row scalar table already has an affine unitriangular
inverse-entry representation of order $`O(\sqrt Q)`$ in the retained
prefix-tree construction. Representation existence is therefore resolved.
Physical coefficient access and returned extraction within R98's budgets
are the missing steps.

Substituting the generic completed inverse-word circuit into that
representation costs $`\Theta(Q^{3/2})`$ T gates and work; its stated
depth guarantee becomes $`O(\log Q)`$. Neither matches the existing
general reader. A counting audit shows that a universal affine matrix
embedding needs $`\Omega(Q/\log Q)`$ coefficient positions. More
generally, a fixed number of completed generic inverse queries with
computational-basis affine wiring has a feature-span obstruction to the
square-root T budget. These are charged-substitution and representation
boundaries, not lower bounds for arbitrary native synthesis, quantum
interference or sparse specialization.

The simple dirty loading $`E\mapsto E+U`$ leaves mixed terms even for
a three-by-three inverse. Graded finite-field interpolation does cancel
the common arbitrary mask exactly, using U only as native Clifford
controls and returning every work wire. But it needs linearly many
translated source evaluations at maximum entry separation. Its fully
charged implementation has $`O(d^4\log^2d)`$ T count,
$`O(d^4\log^3d)`$ Clifford count, $`O(d\log d)`$ depth and
$`O(d^3\log^2d)`$ dirty work, with zero clean work and zero error.
The direct read-only path-propagation baseline already has cubic counts,
linear depth and quadratic dirty work. Interpolation is therefore stopped
as an improvement candidate, while its identity and exact interface are
retained in the existing query proof home.

Independent arbitrary masks cannot simply parallelize their cancellation;
each uncoupled factor would have to be mask-independent itself. Likewise,
conjugating a translation by the new garbage inversion gives a nonlinear
fractional transformation of the dirty mask, not the required inverse
translation. Coupled-mask processing, nonlinear work and sparse native
specialization remain outside these failed mechanisms and the counting
boundaries.

## 5. Precise resumption question

**Does an exact three-native-band C3 source necessarily satisfy
$`N_pN_q=0`$ for $`p\in E^\perp`$ and arbitrary q in the cases
left by the two new certificates?**

The immediate missing step is the bad singleton C factor. Its restriction
to an orbit-control sector can remain a native T phase and need not
commute with B. Moving it past B or cancelling its conjugation would not
preserve the Clifford membership used by the half-angle argument. The
present certificates do not justify that move. A proof may instead
certify additional Pauli images for stripping, or control the intrinsic
tableau product directly. Equivalently, it may prove Clifford inverse
images on the missing residual output generators after good-middle
elimination. A counterexample must give the actual native
bands and literal full operator; changing Clifford coordinates cannot
repair a nonzero intrinsic product.

Success would extend 174/58 transport to this remaining C3 source class,
advancing a native degree-versus-depth recurrence without assuming
diagonalization. General non-C3 prefixes, robustness at growing depth,
one-clean work and address mixing would still be open. On the constructive
side, the useful resumption question is a returned native read-only
extractor using coupled masks or sparse implicit coefficients within
R98's simultaneous budgets, rather than another unpriced inverse identity.
Neither stopped candidate retires the optimal-general-depth objective.
