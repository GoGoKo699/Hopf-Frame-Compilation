# Coupled incidence, signed normalization and a complete inverse query

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Native source proof](ISOLATED_NATIVE_CORE_STRUCTURE.md) ·
[Controlled-Pauli transport](CONTROLLED_PAULI_TRANSPORT.md) ·
[Inverse-query construction](NONCOMMUTATIVE_CYCLE_CHARACTER.md) ·
[Previous checkpoint](UNBOUNDED_CONTROL_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`b7646d3bf846791b01eefba149f9ecefac4070be` (PR #179).
It addresses the recorded coupled-incidence gap, the extra domain
premise in controlled-Pauli transport, and the cost of repeated child
inversions in the constructive block schedule. Each result is proved
in its existing authoritative home linked above. This checkpoint records
the decisions and remaining questions; final manuscript drafting is separate.

## 1. Compact revision

The count package is settled under its stated hypotheses. R96 remains
the strongest unconditional complete-frame research upper law. For all
integers $`n\ge1,L\ge6`$, with $`N=2^n`$,

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It uses one clean qubit, dirty width $`b\ge17(L+n+7)`$ for real
frames or $`b\ge17(L+n+8)`$ for the prescribed complex magnitude
completion, and full initialized-isometry error at most $`2^{-L}`$.
Every logical column, literal sign and phase, singular and zero
coefficient, arbitrary dirty input and external reference, actual inverse,
leakage and complete workspace return remains included. T-count,
elementary Clifford count, T-depth, precision and simultaneous workspace
are separate budgets. Every source, query, control, recovery and inverse
is charged. R79 remains selected Result D, R97 is a scoped four-layer
obstruction, and R98 remains conditional on its exact multi-output native
read-only query interface.

The exact general question is the optimal worst-case T-depth of the
prescribed complete frame under those simultaneous resource budgets.
At fixed accuracy and square-root width, must some family have growing
depth, or can the additive logarithmic-in-n upper term be reduced?
Unrestricted one-clean four-layer rigidity and dimension-growing
unrestricted lower bounds remain open. This round does not improve R96.

Formerly open questions already resolved before this round include the
one-clean count frontier, the all-precision general upper law,
zero-clean read-only five-layer degree 14, six-layer Pauli-comparison
degree 51, diagonal-C3 transport 78/26, source closure through overlap
two, arbitrarily many isolated paired cores, and transport through the
previously scoped controlled-Pauli class. Old next-step passages are
historical evidence. The three restrictions removed in this round are
described below; none is a claim of unrestricted optimal depth.

## 2. Coupled incidence: the missing native restriction

The first candidate tried to bound the rank of the entire weighted
incidence Gram matrix from a single source involution. That inference
was too strong: the involution controls only the active rows in that
probe. Even control-flip weight also permits binary involution corrections,
so an odd-order-only argument would miss genuine cases. A generic
wide-row derivative need not be nonzero when several relations coexist.

The successful reformulation retains the exact relative source block.
The half-angle lemma gives a Pauli correction without assuming the
middle-band restriction is Clifford. The full involution defect has rank
twice the rank of the active B-side labels. A single Pauli correction
forces those labels to coincide and localizes the correction in one
plane, including for even control flips and with all spectator actions.

Native A isotropy then supplies the missing restriction. B axes in the
same nonzero activation class differ by an A-span label, so every actual
A axis commutes with them uniformly. The localized correction's complete
B-neighborhood must therefore be the whole class indicator. The independent
control relations require their A-coordinate columns to span all control
directions. Together these facts force one actual C orbit per interacting
class. Only after that reduction is the single-row derivative rank
$`k-2`$ applicable, forcing class size two.

Thus, for an exact full-space three-native-band C3 source, **if every
zero-character B axis commutes with C, all B/C interactions are isolated
paired cores**, at arbitrary width and overlap. The retained isolated-core
proof gives diagonal C3 in fixed Clifford coordinates and six-layer
transport degree 78/26 below full error $`\sqrt2/8`$.
All native phases and remaining axes stay in the exact operator proof.
This does not infer C3 from return, bound optimal source depth, or exclude
bad zero-character interactions.

The overlapping rows 111100 and 001111 provide a useful bounded check:
all 64 selections satisfy the rank identity. Among the 48 nonzero
activations and 256 correction labels per activation, 104 are binary
involutory corrections, but none survives the necessary native
same-class condition after excluding the pure-B correction. The positive
binary cases are retained as evidence against the rejected stronger rank
and parity shortcuts, rather than deleted as inconvenient examples.

## 3. The controlled-Pauli domain premise is derived

The [signed normalization](CONTROLLED_PAULI_TRANSPORT.md) works for
any scalar-Pauli computational-block source, before C3 is assumed.
Ordered controlled-Pauli Cliffords straighten its mixed input translations
and match their output translations. Their target axes may anticommute;
the order-dependent scalar is kept. A diagonal control Clifford matches
the remaining output shear and every sign. The resulting identity factors
are extracted on all sectors, without initialized wires.

A second signed correction makes the target Pauli domain commute with
the source. A target Clifford separates its symplectic part from its
radical; the former is an identity spectator and the latter becomes
additional computational controls. This produces fixed left and right
Cliffords A,B such that $`H=AFB`$ still has scalar-Pauli blocks and
its exact domain is $`\mathcal P_S\oplus Z_R`$. No domain or
fiber-compatibility premise remains. Both endpoint coordinate changes,
all ordered phases and actual inverses are explicit.

If the original source is C3, H is C3. Re-deriving the quadratic labels
and weighted scalar phases in these coordinates permits the existing
four-block transport argument: label/visible degree at most **58** and
literal sign/general degree at most **174**, below $`\sqrt2/8`$,
independently of control and target width. The normalization needs no
extra physical wire and is an analytical coordinate change, not a free
native synthesis. Arbitrary C3 sources and arbitrary controlled-Clifford
blocks are not classified by this result.

A five-wire fixture checks every column, its full exact domain, mixed
translations with noncommuting section targets, essential affine signs,
radical promotion and both endpoint directions. It complements the
retained attributed Gottesman–Mochon source and phase-carry checks.

## 4. A complete logarithmic-depth structured query

The earlier three-shear proposal demanded a returned helper inside each
in-place block inversion. It left an inverse counter-action on the helper
and required three sequential child-inversion positions. Its
frozen-coefficient determinant obstruction remains correct for that
interface. The useful reformulation is to require only the principal
inverse coordinates during computation and return the helper at the
end of the completed query.

For unit upper blocks A,C, put $`N(E)=AEC`$. With arbitrary B,E,
perform chronologically: $`B\mathrel{+}=N(E)`$; recursively invert
A,C once in parallel; $`E\mathrel{+}=N^{-1}(B)`$ using their live
inverse coefficients; swap B,E. Exactly,

```math
B_{\rm out}=A^{-1}B_{\rm in}C^{-1},\qquad
E_{\rm out}=B_{\rm in}+AE_{\rm in}C.
```

The principal inverse is independent of all helper inputs. The helper
counter-action remains an invertible permutation and is kept live.
The paid two-sided product XOR returns its temporary work at each call,
even when that work contains correlated child garbage. There is one
recursive child position, giving logarithmic T-depth.

For a d-by-d matrix, write $`q=\binom d2`$,
$`a=\lfloor d/2\rfloor`$, $`c=\lceil d/2\rceil`$, and
$`h_d=\max\{\lceil\log_2d\rceil-1,0\}`$. Compute the principal
inverse, CNOT it into the arbitrary output word, flip the diagonal
output bits, and execute the actual inverse computation. This complete
query returns the coefficient inputs and every helper, on all columns
and with arbitrary external references. Its exact ledger is

```math
\begin{aligned}
T_Q&=64\binom d3,&G_Q&\le160\binom d3+23q+4d,\\
D_T(Q)&\le128h_d,&W_Q&=ac(c+3)/2.
\end{aligned}
```

W counts extra arbitrary work beyond q coefficient inputs and
$`d(d+1)/2`$ arbitrary output bits; d equal to one has W zero.
There is no clean work and the error is zero. The d equal to two
computation retains Clifford-only helper garbage, so W is two there.
The actual reversed query replaces the last 4d Clifford term by 8d
in the fixed H,S,CNOT alphabet. No ideal self-inverse substitution hides
that charge.

A direct controlled readout borrows q current garbage bits as private
Toffoli helpers, returning their current state before uncomputation.
It has $`T\le64\binom d3+8q`$,
$`G\le160\binom d3+42q+d`$,
$`D_T\le128h_d+8`$, the same W and the same actual-inverse ledger.
The phase source uses the analogous diagonal readout and has the same
T/depth/work orders, with its own literal Clifford ledger in the proof.
No precision state, query, helper initialization or inverse is free.

This improves the complete structured source/query from linear to
logarithmic T-depth by increasing work from quadratic to cubic. The
coefficient inputs are modified internally, so the result is **not**
native read-only in those inputs. The separate transport in the cycle
conjunction is unchanged and still loses to the fairly simplified
conjunction baseline. Arbitrary-table encoding, R98's read-only interface
and its simultaneous resource budgets are additional missing steps.

## 5. Most useful next question

**Do bad zero-character B/C interactions in an exact three-native-band
C3 source admit scalar-Pauli block coordinates, or another invariant
that bounds signed transport independently of width?**

Nonzero-character incidence is now closed under the stated good-zero
hypothesis. The broader controlled-Pauli theorem means the next proof
need not force diagonalization: a signed scalar-Pauli block reduction
would already give 174/58 transport. Conversely, a countermechanism must
be an actual native three-band source with complete signs and inverses,
and must be tested under alternative native representations before it is
used against a lower-bound parameter.

Success would complete another source class needed for a native
degree-versus-depth recurrence. General non-C3 prefixes, growing-depth
error control, one-clean work and address mixing would still separate
that result from an unrestricted frame lower bound. On the constructive
side, the concrete resumption question is whether the completed inverse
query can be given a native read-only coefficient interface and useful
arbitrary-table encoding within R98's simultaneous budgets, or whether a
different representation avoids that conversion cost. The source's
logarithmic depth alone does not answer this.

Exhausting either mechanism leaves the general optimal-depth objective
active. The proof homes state the exact unresolved interfaces; this
checkpoint retires only the disproved stronger rank shortcut and the
unnecessary requirement of returning every helper before readout.
