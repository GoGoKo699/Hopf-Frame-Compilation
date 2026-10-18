# Overlap-two sources, bounded-core transport, and commuting inverse pivots

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Source proof](TWO_OVERLAP_SOURCE_STRUCTURE.md) ·
[Hybrid transport](HYBRID_DIAGONAL_CORE_TRANSPORT.md) ·
[Cycle source](NONCOMMUTATIVE_CYCLE_CHARACTER.md) ·
[Previous checkpoint](ONE_OVERLAP_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`d2ecebfd36a965eb126db057a80aefade76954c1` (PR #177).
It addresses the preceding checkpoint's overlap-two source question and
the complete native cost of unitriangular inversion. This is a research
checkpoint; final manuscript drafting remains separate.

## 1. Compact revision

The count package remains settled under its stated hypotheses. R96 is
the strongest unconditional research upper law, for integers
$`n\ge1,L\ge6`$, $`N=2^n`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It uses one clean qubit, with dirty width $`b\ge17(L+n+7)`$ for real
frames or $`b\ge17(L+n+8)`$ for the prescribed complex magnitude
completion, and full initialized-isometry error at most $`2^{-L}`$.
Every logical column, literal sign and phase, singular and zero
coefficient, arbitrary dirty input and external reference, actual inverse,
leakage and complete work return is included. T-count, elementary Clifford
count, T-depth, precision and simultaneous workspace are distinct budgets;
every source, query, control, recovery and inverse is charged.
R79 remains selected Result D, R97 remains a scoped four-layer
obstruction, and R98 remains conditional on a faster exact multi-output
native read-only query.

The exact general question is the optimal worst-case T-depth of that
complete frame under these simultaneous budgets. At fixed accuracy and
square-root width, can the logarithmic-in-n additive depth be reduced,
or must some legal frames require growing depth? Unrestricted one-clean
four-layer rigidity and dimension-growing unrestricted lower bounds are
still open. None of this round's results changes R96.

Formerly open issues already resolved include the one-clean count
frontier, the all-precision general upper law, zero-clean read-only
five-layer degree at most 14, six-layer Pauli-comparison degree at most
51, diagonal-C3 transport degree 78/26, and the overlap-one source
certificate. The bounded-deficiency transport theorem also remains
valid, but its parameter can grow within one native T band.

## 2. The source question addressed

The previous theorem bounded middle/outer noncommutation by at most 2t
middle axes, where $`t=\dim(\mathrm{span}(A)\cap\mathrm{span}(B))`$
for the specified source $`F=E_CE_BE_A`$. That count did not control
coupling to a wide diagonal source or prove a small residual core.
This round examines the two-dimensional relation space analytically,
retaining the full source C3 hypothesis and each band's independent
physical T axes. No intermediate C3 premise or inverse-C3 premise is
inserted.

The [authoritative source home](TWO_OVERLAP_SOURCE_STRUCTURE.md) proves:
**every exact C3 three-native-band source with first/middle span overlap
at most two is diagonal C3 in fixed Clifford coordinates.** Neither
adjacent partial source needs to be C3. The existing six-layer query
degree bounds **78**, or **26** for visible outer shear, therefore apply
below full error $`\sqrt2/8`$ throughout this larger source class.

The activation-code classification separates repeated nonzero characters
from zero characters. The former give paired controlled cores; the
latter can have arbitrarily large raw outer-axis spread through exact
cancellation. Raw exceptional-axis counting is therefore not a valid
replacement for the full source argument.

With no zero activation character, at most two isolated paired cores
remain. A rank-four obstruction prevents one isolated inner axis from
flipping both controls. After extracting the inner pairs, the remaining
source is itself C3; its control involutions force the fixed target-axis
signs to be affine. A complete commuting Pauli frame is preserved,
giving common diagonal coordinates before restoring all remaining
inner axes by native extension.

Two zero characters reduce to duplicate absorption or a three-qubit
core. The nonduplicate core violates the full-space symplectic
involution condition. The mixed zero-plus-pair case first localizes the
inner axes with paired probes, rules out a pair of isolated outer axes
by independence, then excludes its last core by another full-space
involution. No radical eigensector is fixed in these contradictions.
The surviving duplicate case reduces to overlap one with its literal
right Clifford and scalar retained.

This resolves the preceding checkpoint's exact overlap-two question.
It does not infer C3 membership from return, give an approximate source
classification, or bound the overlap of an arbitrary shallow source.

## 3. A transport theorem that does not require diagonal classification

The [hybrid transport theorem](HYBRID_DIAGONAL_CORE_TRANSPORT.md)
handles either $`D(U_A\otimes I)`$ or $`(U_A\otimes I)D`$, up to fixed
input/output Cliffords. Here D is an arbitrarily wide diagonal C3
source, U is an arbitrary fixed unitary on a qubits, and the full source
is C3. Put

```math
N(a)=3\,2^{5a}+15a+2,\qquad
B(a)=24N(a)+8a^2+12a.
```

For positive a and full error $`\eta\lt\sqrt2/8`$, the actual native
six-layer comparison has label degree at most $`4B(a)+2`$ and literal
fixed-probe sign degree at most $`12B(a)+6`$. The same bounds apply to
visible outer shear and general target degree, respectively. The
constants depend only on the core size, even when the exact
Pauli-preserving domain has unbounded codimension. A zero- or one-qubit
core has the sharper existing bounds 26/78.

The proof strips exact spectator directions, matches the entire signed
Pauli domain, and retains an explicit bounded list of phase features for
each remaining block. At most four fixed control blocks recover each output label coordinate.
The fourth-root scalar and modulo-four Hermitian-section carry then
recover the literal sign; otherwise a wide sum of signs would reopen
the degree gap. The radius is uniform in the source and dirty width.

A reusable application handles a diagonal-C3 source followed by a
commuting T-gadget list with relation dimension r, when the full product
is C3. Duplicate axes are collected with their actual scalar and Clifford
factors before applying paired half-angle rigidity. At most
$`4(2^r-1)`$ surviving axes have non-Pauli images. A monomial Clifford
localizes their isotropic span while preserving the wide diagonal phase,
giving a hybrid core of at most that many qubits. This yields a transport
bound depending only on r. It neither makes dependent gadgets one
physical T band nor turns r into an unrestricted native-depth invariant.

All coordinate changes here are analytical. They supply no free physical
Clifford, helper, source, inverse or control, and they do not infer C3
membership from approximate workspace return.

## 4. Constructive mechanism and its stopping diagnosis

The [cycle-source proof](NONCOMMUTATIVE_CYCLE_CHARACTER.md) already
identified its character with $`\mathrm{tr}((I+U)^{-1}L)`$ over the
binary field. This round tests an exact reversible inversion mechanism.
For every internal index k, one whole pivot adds
$`U_{ik}U_{kj}`$ to every entry $`U_{ij}`$ with $`i\lt k\lt j`$.
Each pivot is an involution. The whole pivots commute, and their product
maps every strictly upper triangular U to $`(I+U)^{-1}-I`$.

The path expansion proves this on every input. An arbitrary subset of
pivots inserts exactly the paths whose internal vertices lie in that
subset; composing pivot subsets takes their symmetric difference.
The native implementation retains the private-helper parallel Toffoli
schedule and charges

```math
T\le8\binom d3,\qquad G\le20\binom d3,\qquad
D_T\le8\max(d-2,0).
```

For $`d\ge2`$, at most $`\lfloor(d-1)^2/4\rfloor`$ simultaneous
arbitrary helpers suffice in addition to the $`d(d-1)/2`$ data bits;
no clean work is needed. The map and its actual inverse return every
helper exactly, with zero error and literal phase one.

Commutativity of the whole pivots does not parallelize their overlapping
native gates. Individual Toffolis cross control and target roles;
cancellation occurs only after the whole pivot families are included.
Applying all quadratic updates from the original entries gives only
$`U+U^2`$ and loses longer paths. The charged construction still has
linear depth and changes U internally, so it supplies neither a faster
native read-only inverse query nor R98's arbitrary-table reader. It
does not beat the fairly simplified conjunction baseline. The useful
identity is retained; its unsupported depth saving is stopped.

## 5. Most useful resumption question

**Does the actual six-layer inner family have a width-independent signed
transport-degree bound through an arbitrary three-native-band C3 source,
without a bound on the first/middle overlap?**

The hybrid theorem resolves the wide-diagonal coupling step once a bounded
core factorization is available. The remaining challenge is to control
native relation components as their rank grows, or to produce a genuinely
native high-degree transport family that escapes such control. A useful
first test is a connected component with at least three independent
relation directions. A bound that merely grows with displayed overlap
does not settle the question. Alternative native representations and
fair Clifford simplifications must be checked before treating any source
parameter as a lower-bound invariant.

Success would advance the native degree-versus-depth program by removing
an additional source restriction, rather than only increasing a finite
cutoff. Even that would leave non-C3 prefixes, one-clean work, address
mixing, and a recurrence with a robust radius for growing depth. The
separate query/frame reduction requires an unrestricted query lower
bound before it yields a general frame lower bound.

On the constructive side, the precise unresolved cost is a returned
native read-only implementation of the whole inverse-word interface
faster than the paid pivot schedule, together with a route from that
structured inverse problem to the arbitrary-table query in R98. No free
inverse, initialized cache or control wrapper may be inserted. Exhausting
this pivot candidate does not retire the general objective.

## 6. Evidence and integration

The proof homes contain the arbitrary-width arguments. Focused bounded
checks exercise signed source identities, symplectic involution
obstructions, phase carries, monomial localization and the all-input
pivot identity with its actual inverse and returned helpers. These checks
support the explicit mechanisms; they do not replace the proofs.
Validation and rendering receipts are recorded in the merging pull
request. No manuscript claim selection changes in this round.
