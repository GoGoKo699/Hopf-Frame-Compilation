# Joint bank routing and the complete dirty indicator

[Complete approximate transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Retained indicator](NONUNIFORM_DIRTY_INDICATOR.md) ·
[Shared-control native batch](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)

Routing a selected bank qubit is a legitimate joint mechanism: the
compiler may temporarily change the address and use noncommuting changes
of basis. A complete selected swap is equivalent, up to constant source
calls and linear overhead, to the complete dirty indicator. Thus a faster
selected-swap compiler would be useful, but naming that operation does
not supply a weaker interface. Two concrete attempts below expose their
remaining costs. Neither changes R96 or gives an unrestricted depth
lower bound.

Throughout, let p be positive, m equal to two to the power p, x be a
p-bit address, Y an arbitrary m-bit bank, and t an arbitrary bit. Write

```math
I_p:(x,Y)\longmapsto(x,Y+e_x),\qquad
S_p:(x,t,Y)\longmapsto(x,Y_x,Y^{(x\leftarrow t)}).
```

Addition is over the binary field. S swaps t with the selected physical
bank coordinate. Every helper is arbitrary; all assertions concern the
complete operator, including external references. Common forward and
actual-inverse ledgers use the positive-S elementary Clifford alphabet.

## 1. A selected swap supplies the complete indicator

The literal identity is

```math
I_p\otimes I_t=S_p^\dagger X_t S_p.                 \qquad(1)
```

The chronological word is S, X on t, then the actual S adjoint. The
original value of t returns, even when it is entangled with x, Y or a
reference. There is no condition on intermediate address motion in the
native S word.

Suppose V approximates S tensor identity on b_S arbitrary helpers in
full operator norm by epsilon. Replacing S and its adjoint in (1) by V
and its actual adjoint gives indicator error at most twice epsilon. With
common source bounds, the emitted word has

```math
\begin{aligned}
T&\le2T_S,&D_T&\le2d_S,\\
G&\le2G_S+4,&b&\le b_S+1.
\end{aligned}                                                   \qquad(2)
```

X is charged as HSSH. The extra arbitrary bit is t; no clean bit is
introduced. Leakage and address/helper return errors are included in the
two source norm estimates, and actual reversal has the same common
ledger.

## 2. Six complete indicators supply a selected swap

Reserve an arbitrary m-bit bank A and the b_I private helpers of an
indicator whose output is A. Two exact bilinear shears are

```math
B_G:t\longmapsto t+\sum_i A_iY_i,\qquad
B_S:Y\longmapsto Y+tA.                                \qquad(3)
```

Both are conjugates of the same diagonal operation

```math
C=\prod_{i=0}^{m-1}\mathrm{CCZ}(t,A_i,Y_i).
```

For B_G conjugate C by H on t. For B_S conjugate it by Hadamards on all
Y. These different choices are retained in the actual word; they cannot
be dropped merely because the abstract shears are classical.

The shared-control phase identity gives four native layers for C. The
six nonconstant parities per pair contribute six m T gates. The common
term is T to the power m on t. Since m is even, that common term is a
Clifford and costs at most three positive S gates in either orientation.
The parity compute/uncompute circuits use ten m CNOTs. Consequently

| Word | T-count | T-depth | Clifford count | Extra dirty |
| --- | ---: | ---: | ---: | ---: |
| $`B_G`$ | $`6m`$ | 4 | $`10m+5`$ | 0 |
| $`B_S`$ | $`6m`$ | 4 | $`12m+3`$ | 0 |

These are common forward/actual-inverse upper bounds. In particular,
neither shear needs an initialized or copied t. The shared t participates
as a phase target only in the collected common term.

For B equal to either shear, execute chronologically

```math
B,\quad I_p^{(A)},\quad B^\dagger,\quad
\bigl(I_p^{(A)}\bigr)^\dagger.                         \qquad(5)
```

The two uses of B cancel every arbitrary initial A value. For B_G the
remaining map is G: t goes to t plus Y_x. For B_S it is the scatter map
R: Y goes to Y plus t e_x. The complete word G,R,G is precisely the
selected swap. On the selected pair its chronological action is

```math
(t,Y_x)\longmapsto(t+Y_x,Y_x)
\longmapsto(t+Y_x,t)\longmapsto(Y_x,t).
```

All other Y coordinates, A, x and the indicator's private work return.
The same A and private pool are reused only after completed source
calls. No parallel interleaving of sources that might move x is used.

If V_I is within epsilon of the complete indicator on its whole space,
each of the six source or actual-adjoint appearances has error at most
epsilon. Unit-norm product telescoping therefore gives selected-swap
error at most six epsilon. Its simultaneous common ledger is

```math
\begin{aligned}
T_S&\le6T_I+36m,&D_T(S)&\le6d_I+24,\\
G_S&\le6G_I+64m+26,&b_S&\le m+b_I.
\end{aligned}                                                   \qquad(6)
```

The Clifford constant in (6) includes four B_G and two B_S appearances.
There is no extra bilinear helper pool. The six ideal calls preserve the
arbitrary work, so the hybrid remains valid if actual approximants leak
or fail to return exactly between calls. Every such error is charged.

Equations (2) and (6) show equivalence of the requested fixed-accuracy
envelopes: linear T-count and dirty width, quadratic Clifford count and
sublogarithmic address T-depth. Accuracy is rescaled by a constant, not
by m. This conclusion concerns the complete selected swap. A source V
promising only that V adjoint times X on t times V approximates the
indicator has a weaker operator promise than a full swap, but its
two-call extraction is already (1). Such a source has the same
asymptotic native status when its two calls fit (2); the weaker promise
is a potential synthesis freedom, not an independently cheaper query.

## 3. A direct pruned router has linear count but linear address depth

There is a simple exact routing word with no private work. At level j,
starting with j zero, partition the bank into consecutive blocks of
size two to the power j+1. In each block conditionally swap its first
coordinate and its coordinate with offset two to the power j, using
address bit x_j. There are

```math
R_j=2^{p-j-1}
```

disjoint pairs at that level. After level j the first coordinate of
each block contains the bank entry selected by the low j+1 address
bits. The complete word V_x therefore puts the original Y_x in
coordinate zero. The unselected values undergo a reversible permutation;
none is assumed zero or discarded. Thus

```math
V_x^\dagger X_0V_x=I_p.                               \qquad(7)
```

The shared-control Fredkin schedule uses at most four T layers,
six R_j plus its parity in T-count, and fourteen R_j plus three
Cliffords, with no helpers. Only the last R_j is odd. Summing all
levels and charging both router orientations in (7) gives

```math
\begin{aligned}
T(I_p)&\le12(m-1)+2,&D_T(I_p)&\le8p,\\
G(I_p)&\le28(m-1)+6p+4,&b&=0.
\end{aligned}                                                   \qquad(8)
```

Equation (8) is a complete exact baseline, including literal phases.
It has a useful count and no helper, but its critical path still visits
p successive decisions. The ordinary full bank translation, which
swaps all m/2 pairs at every address bit, also retains those p positions
and instead has order p m T-count. Global simultaneous commutativity of
the ideal full translations does not supply a simultaneous native
schedule for their overlapping physical targets.

This is the first unpaid step in the routing proposal: replacing the
successive coherent decisions by a joint native implementation. Grouping
several address bits into one selected router merely moves that task
into the grouped primitive unless its whole native ledger is supplied.
The retained indicator already improves (8) to logarithmic address
depth using arbitrary dirty work; (8) is not offered as a frontier
improvement.

## 4. Fixed Clifford re-encoding cannot diagonalize the bank translations

One proposed flattening uses a fixed Clifford to simultaneously convert
the full XOR bank translations into diagonal Cliffords. This fails
already at the symplectic action, without a Gaussian restriction.

Let P_i be the m by m binary permutation matrix for translation of bank
labels by the ith address coordinate. The physical bank permutation has
symplectic action

```math
S_i=\begin{pmatrix}P_i&0\\0&P_i\end{pmatrix},\qquad
N_i=S_i-I.                                             \qquad(9)
```

The two blocks agree because P_i is an involutive permutation. For any
set J of distinct address coordinates, the product of P_i minus identity
over i in J sends the zero-coordinate basis vector to the sum of the
distinct coordinate vectors labeled by the span of those directions.
It is nonzero. In particular,

```math
\prod_{i\in J}N_i\ne0
\quad\text{for every nonempty }J\subseteq\{1,\ldots,p\}.
                                                               \qquad(10)
```

An arbitrary computational-diagonal Clifford has a symplectic matrix
of the form

```math
D_B=\begin{pmatrix}I&0\\B&I\end{pmatrix}.
```

The product of the differences from identity of any two such matrices
is zero. A single fixed Clifford change of coordinates conjugates all
S_i by the same invertible symplectic matrix and preserves whether the
products in (10) vanish. It therefore cannot turn every bank translation
into a diagonal Clifford once p is at least two.

Adding arbitrary helpers that return identically to each completed
translation tensors its symplectic matrix with identity. The nonzero
product persists, even if the proposed fixed Clifford mixes bank and
helper coordinates. Literal phases cannot repair a contradiction in
the symplectic action. Conversely, this statement does not constrain a
non-Clifford change of coordinates, an address-dependent re-encoding,
a joint native word that never supplies these translations separately,
or a source whose helper action is retained for cancellation elsewhere.

This is a limitation of simultaneous diagonalization, not a lower bound
on the unrestricted complete indicator or selected swap. Its nonzero
products count neither native T layers nor T gates.

## 5. Classical round boundaries retain a logarithmic degree requirement

Temporary address/bank mixing is permitted by the target contract, but
it does not by itself accelerate a recurrence whose completed rounds
remain classical coordinate maps. This distinction can be stated without
requiring any register to stay read-only.

Suppose a candidate is a composition of r monomial reversible rounds,
with arbitrary affine changes of all coordinates between them. Each
round may have arbitrary phases, but its underlying Boolean permutation
has coordinate degree at most a fixed integer k, at least two. All
address, output and private dirty coordinates may be mixed. Assume the
whole native word has full-operator distance less than square root two
from the complete indicator tensor identity on its helpers. Then

```math
k^r\ge p,\qquad r\ge\lceil\log_k p\rceil.             \qquad(11)
```

To prove this, a monomial unitary sends each computational basis vector
to one basis vector with a unit-modulus coefficient. If that vector has
the wrong destination, its distance from the required output is exactly
square root two, regardless of phase. The error hypothesis therefore
forces the entire final basis permutation to equal the indicator's,
including helper and address return. Polynomial substitution bounds the
degree of the composed coordinate maps by k to the power r; affine
changes do not increase it. Fix every helper to zero. The coordinate
corresponding to Y at label zero is

```math
Y_0+\prod_{j=0}^{p-1}(1+x_j),
```

whose degree is p. Restriction cannot increase degree, proving (11).
No output or dirty bit was assumed initialized in the source contract;
restriction to one allowed input subspace is only used for this necessary
condition.

For quadratic reversible rounds, (11) retains logarithmically many
nonlinear rounds even after arbitrary affine address/output mixing. It
applies to proposed recurrences with those classically closed round
boundaries, including exact circuits and monomial approximants. It does
not convert r into a lower bound on unrestricted native T-depth. A
Clifford between native T layers may be nonmonomial, and a final
monomial action does not imply monomial intermediate boundaries. A joint
word that keeps Hadamard interference across those proposed cuts is
outside the premise. Growing-degree rounds are also allowed by the
unrestricted model, but their native implementation must be paid rather
than counted as one elementary decision.

## 6. What a new routing mechanism must change

A complete selected swap is a sufficient source with the same asymptotic
status as the indicator, and it permits intermediate address motion.
Its native construction remains open at sublogarithmic depth. The
pruned router pays every switch and retains a linear address critical
path; a fixed Clifford simultaneous diagonalization cannot eliminate
that path by making all full translations diagonal. Allowing arbitrary
affine address/output mixing still leaves logarithmically many quadratic
classical rounds in the precise setting of (11).

The remaining constructive opportunity is a whole word that transports
the selected Pauli, or performs the swap, without separately synthesizing
the full translation family or retaining the router's sequential
decisions. Nonmonomial interference across the classical round boundaries
is one remaining freedom; a fixed re-encoding or a relabeled selected
swap does not supply it. Such a word may use noncommuting basis changes
and temporary address/bank mixing. It must supply full-space error, actual reversal,
and the simultaneous native count, width and depth ledger before the
paid R98 transfer applies. The results here impose no unrestricted
obstruction on that opportunity.
