# Joint encodings: unrestricted one-clean two-layer Boolean rigidity

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Authoritative proof](TWO_LAYER_DIAGONAL_RIGIDITY.md) ·
[Previous checkpoint](JOINT_TRANSFER_SESSION_RESULT.md) ·
[Conditional upper transfer](QUERY_DEPTH_FRAME_TRANSFER.md)

This round began on 18 October 2026 at verified main
`3c0fc6284c2e9cbf4e231b20c5244e441a053adc` (PR #183). It addresses
the joint-encoding gap left by the six-bit cubic bent query, after
separating the general destination from the present representation.
The supported result is an unrestricted one-clean two-layer degree
obstruction, including a width-independent error floor. An explicit
clean-removal construction and a native scope counterexample identify
why the proof needs both the complete operator and literal Boolean phase.
The complete-frame resource frontier is unchanged.

## 1. Revision before choosing the question

The destination is optimal worst-case general T-depth for the prescribed
complete Hopf differential frame, under simultaneous count, precision
and workspace budgets. R96 remains the strongest unconditional general
upper law. For $`N=2^n`$, all integers $`n\ge1,L\ge6`$, it gives

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

There is one clean qubit and $`b\ge17(L+n+7)`$ dirty qubits for
real frames, or $`b\ge17(L+n+8)`$ for complex magnitude frames.
The full initialized-isometry error is at most $`2^{-L}`$.
Every logical column, literal sign and phase, singular or zero
coefficient, arbitrary dirty input, external reference, actual inverse,
leakage and complete work return remains covered. All source, query,
control and recovery costs are charged at simultaneous allocations.
R79 remains selected Result D. R97 retains its outer-overlap hypothesis.

The count package is settled under its stated hypotheses. The exact
large-width question is whether the additive depth can be improved
while keeping these budgets, or a matching unrestricted lower bound
can be proved. At fixed accuracy and square-root dirty width, must
depth grow with n? Count divided by physical layer capacity does not
settle this regime. Unrestricted one-clean four-layer rigidity remains
open as well.

Useful completed components are distinct from this destination:
R98 accepts fully restored addresses, including internal changes;
the arbitrary-table inverse-entry representation and structured
logarithmic-depth inverse query already exist; the unrestricted
two-layer occupied transfer alphabet and AND3 obstruction are proved.
Read-only degree bounds, source certificates and outer-overlap bounds
retain their hypotheses. Historical next-step passages do not reopen
resolved count or representation questions.

## 2. Constructive check and first unpaid cost

The target left by the preceding round is
$`f(x,y)=x\cdot y+x_1x_2x_3`$, whose complete scalar transfer
matrix passes the one-clean two-layer alphabet. Put
$`A=\mathrm{Toffoli}(x_2,x_3;y_1)`$ and
$`B_i=\mathrm{Toffoli}(x_i,y_i;z)`$. The exact identity

```math
Q_f=A(B_1B_2B_3)A
```

restores every address and output column with literal phase one.
The two A occurrences are paid nonlinear operations. Using the
existing two-layer dirty-helper Toffoli for each occurrence gives a
transparent **40 T, 90 elementary Clifford, ten T-layer** baseline,
with no additional clean or dirty workspace: other address wires
outside each triple are arbitrary returned helpers. The actual reversed
word has the same ledger and exact error zero. This is a finite
benchmark, not an optimality claim. Sharing z does not make the three
B decompositions a legal simultaneous physical band.

For the stated Clifford count, the first four-wire parity basis is
$`(a+d,b+d,c+d,d)`$, computed by three CNOTs, with T signs
$`(-,-,-,+)`$. The second is
$`(b+c+d,a+c+d,a+b+c+d,a+b+d)`$, computed chronologically by
$`a\to b,b\to d,d\to c,c\to a,a\to b`$, with signs
$`(+,+,-,+)`$. Compute, phase and actual uncompute each basis;
the resulting eight-parity word is CCZ independent of d. Two output
Hadamards give Toffoli, costing 16 CNOT and two H. The
[exact baseline fixtures](../../tests/test_bent_query_baseline.py)
check all logical columns, helpers, inverses and this ledger.

Output and address Hadamards leave the quartic phase
$`(-1)^{zx_1x_2x_3}`$ alongside quadratic shears. Multiplying by
$`Q_{x\cdot y}`$ is non-Clifford, and setting y to zero consumes
three new clean promises. Neither is a free reduction to AND3. More
generally Clifford pre/post multiplication only signed-permutes transfer
rows and columns: this bent query has magnitudes
$`0,1,1/2,1/4,1/8`$, whereas AND3 has an entry of magnitude 3/4.
Dirty identities do not remove that mismatch. No faster arbitrary-table
query within R98's budgets results from these identities.

This diagnosis motivated a joint lower-bound question: what does
preserving the entire diagonal Pauli algebra impose on the two native
encoders, beyond their separate scalar overlaps?

## 3. Supported theorem and the necessary scope boundary

The [proof home](TWO_LAYER_DIAGONAL_RIGIDITY.md) proves that an exact
one-clean native circuit of T-depth at most two, whose occupied action
is diagonal and real up to one common phase, has diagonal C3 occupied
action. Arbitrary dirty width and Clifford/address mixing are allowed.
Conjugating an XOR output by H gives the sign phase
$`(-1)^{z f(x)}`$. Therefore **every exact qualifying Boolean XOR has
algebraic degree at most two**. The six-bit cubic bent query is now
excluded in the unrestricted one-clean model, despite passing the entire
scalar transfer alphabet.

With the retained rounding results, every nonquadratic Boolean XOR has
full initialized-isometry error at least

```math
\tau=\frac{\sin(\pi/8)}{2\sqrt2}\approx0.135299
```

at two native T layers. Thus precision $`2^{-L}`$ for integer
$`L\ge3`$ requires at least three layers, independently of arity,
dirty width, T-count and elementary Clifford count. This replaces the
shrinking AND-arity floors for this depth. Leakage, arbitrary references,
literal phases and complete work return are included.

The proof first handles zero clean work. Preserving the entire diagonal
Pauli algebra forces the two physical bands to share one commuting axis
family, using intrinsic active spans and their affine Pauli-support
cosets. The actual full-space operation is consequently diagonal C3.
For a non-Pauli common clean midpoint, a shared anticommuting axis gives
two explicit controlled-Pauli Clifford dressings. They extend the occupied
action to its tensor product with identity on the clean wire, at unchanged
physical width, T-count and T-depth. The forward Clifford overhead is
at most $`10(q-1)+6`$; the transformed inverse has this overhead relative
to its **separately charged** inverse word, using compiled exact dressing
inverses. No equality of a native word's forward and inverse Clifford
counts is assumed.

The blanket one-clean diagonal extension is false without the real-phase
condition. A fully native circuit on five data wires and one clean wire
has occupied action $`(iZ_t)^{q(x)}`$, where
$`q=\sum_i x_i+\sum_{i\lt j}x_ix_j`$ on four x bits. It uses
12 T gates, two T layers and 44 elementary Cliffords, returns the clean
qubit exactly, and has a nonzero quartic phase coefficient. Its raw
inverse costs 52 Cliffords; adjacent-S reduction gives 44. The unused
sector is nondiagonal. Exact full-column checks retain every sign.
After output Hadamards this becomes $`i^qX_t^q`$; dropping the phase
would change the requested operator.

## 4. How the remaining native phase was resolved

For a Pauli clean midpoint, restriction to its stabilizer code creates
at most one binary relation in each physical band's projected axis list.
The relation can hide one direction from the intrinsic active span.
The exceptional case reduces to a conditional Pauli with quadratic
control q. Two such hidden directions can leave the factor
$`i^{\kappa q}`$. An arbitrary diagonal C3 correction can cancel its
imaginary parity while retaining a quartic carry. Replacing the remaining
native factors by an unrestricted C3 operation would therefore invalidate
the obstruction.

The successful reduction retains those factors literally. Silent phases
are first moved to the input side of their own encoder, where they are
diagonal; the active transition remains diagonal. Fixed and affine
Clifford factors preserve the remaining physical axes. When the two
quadratic controls are independent, the hidden directions split into
commuting C3 factors. When their nonzero quadratic parts agree, the only
possible higher-degree phase is the shared $`i^{\kappa q}`$.
A signed simultaneous Clifford extraction then identifies the two full
Pauli algebras. It reinserts both silent lists in common commuting
coordinates, with each native list still independent modulo the
quadratic-control support.

The resulting row-Gram identities force the quadratic part of any
available imaginary-phase correction to vanish. It cannot cancel q
when kappa is odd. This argument uses the original independent native
rows, rather than a freely supplied phase correction. A coordinate issue
is resolved by squaring the full residual operator P: projective reality
makes $`P^2`$ a scalar Pauli. Because P belongs to the identified commuting
Pauli algebra, that Pauli belongs to the same algebra. Its phase modulo
four is therefore affine in those intrinsic coordinates. No claim that
the entire commuting algebra is a computational Z subspace is needed.
Kappa must be even, and all surviving occupied factors are diagonal C3.

The robust statement does not assume exact reality from approximation.
Below tau, R75 forces the same physical circuit to return its clean
qubit exactly. The occupied transfer alphabet then rounds all required
Z-preserving overlaps to one, making the occupied action diagonal.
The retained unit-circle lemma makes its diagonal entries eighth roots
of unity. Distinct ratios to the Boolean sign target have optimized
common-phase distance at least $`2\sin(\pi/16)\gt\tau`$; hence an
approximation below tau must equal that target up to one common phase.
The exact theorem then applies.

## 5. Current resumption question and its relation to general depth

The six-bit bent query and the unrestricted one-clean two-layer Boolean
degree question are resolved. The first unresolved lower question is
now: **is Boolean degree bounded at three native T layers with one clean
qubit and arbitrarily many returned dirty wires, or can a complete native
family realize unbounded degree at that depth?** A construction must
retain literal phases and all occupied columns; an obstruction must
survive arbitrary intermediate Clifford mixing. The present proof's
one-band midpoint and single-relation structure cannot simply be assigned
to a two-band encoder.

This is the next useful test of whether the new joint invariant can
support dimension-growing depth lower bounds. A uniform degree bound
would extend the unrestricted base case; a fully charged counterfamily
would identify a mechanism the general lower-bound program must admit.
Neither outcome by itself determines optimal general depth.

The current result leaves the complete-frame frontier unchanged. The
existing two-call query/frame reduction turns a three-layer query lower
bound into only a two-layer frame lower bound, weaker than the retained
complete-frame obstruction. R96, selected Result D R79, and the scoped
R97 theorem are unchanged. On the constructive side the continuing
question is an exact sublogarithmic restored-address arbitrary-table
query within R98's simultaneous T-count, Clifford, inverse and workspace
budgets. The paid 40-T, ten-layer bent-query benchmark above is an exact
finite baseline, not that missing arbitrary-table compiler. Final
manuscript drafting remains separate.
