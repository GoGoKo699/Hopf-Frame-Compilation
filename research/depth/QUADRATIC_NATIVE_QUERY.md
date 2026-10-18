# Exact two-layer quadratic queries with one returned dirty qubit

[Unrestricted two-layer degree bound](TWO_LAYER_DIAGONAL_RIGIDITY.md) ·
[Existing full-input one-layer floor](../../docs/PARALLEL_DIRTY_LOOKUP.md#11-a-full-input-lower-bound) ·
[Exact native fixtures](../../tests/test_quadratic_native_query.py)

Every Boolean quadratic XOR query admits two native T layers with zero
clean work and one arbitrary returned dirty qubit. This supplies the
constructive side of the existing degree obstruction. The earlier
full-input one-layer floor makes depth two optimal for every nonaffine
quadratic query at this allocation. The count and Clifford bounds below
are constructive upper bounds; no count optimality or external priority
is asserted.

This is a one-output query component. It does not directly change the
complete-frame upper law or provide a dimension-growing depth lower
bound. In particular, several queries cannot share an address or a dirty
helper inside these bands without a new simultaneous schedule.

## 1. Statement and native ledger

Let f be a Boolean polynomial of degree at most two on m address bits.
Its polar matrix B has zero diagonal and off-diagonal entries equal to
its quadratic coefficients. Write r for its even binary rank. For r>0,
there is an exact Clifford+T circuit

```math
|x,z,d\rangle\longmapsto|x,z\mathbin\oplus f(x),d\rangle
```

using zero initialized qubits and exactly one arbitrary dirty helper d,
with

```math
T=2r+4,\qquad D_T=2,\qquad G=O(m^2+r).
```

All m address bits, both output values and both dirty values are
included in the literal equality. The construction has zero error and
returns every physical wire other than the prescribed output. It
therefore preserves arbitrary entanglement with external references.
No additional simultaneous workspace is hidden in the coordinate
changes or the phase evaluation.

More precisely, let L_A be the CNOT count of the address coordinate
change from section 3, let ell be the residual linear coefficient vector,
and let beta=f(0). Define

```math
s_r=\begin{cases}0&r=0,2\pmod8,\\3&r=4,6\pmod8,\end{cases}
\qquad
s'_r=\begin{cases}0&r=0,2\pmod8,\\5&r=4,6\pmod8.\end{cases}
```

The forward word and its actual reversed native inverse have respective
primitive Clifford counts

```math
G=2L_A+8r+8+s_r+|\ell|+4\beta,\qquad
G_{\rm inv}=2L_A+8r+8+s'_r+|\ell|+4\beta.       \qquad(1)
```

Here inverse S gates are expanded into three S gates and only adjacent
S powers are reduced modulo four. X is explicitly HSSH, costing four
primitive Cliffords. Gaussian elimination gives
$`L_A\le m(m-1)+3(m-1)`$ for m positive: at most one three-CNOT row
swap per nonfinal pivot and at most m-1 row additions per column.
Neither a linear map nor its inverse is supplied for free.

If r=0, f is affine and its query uses only at most m CNOTs and one
HSSH word. Its optimal T-depth is zero. If r>0, f is nonaffine and the
existing full-input theorem excludes T-depth zero and one, even with
arbitrarily many returned dirty helpers. Thus

```math
D_T^\star(f;\hbox{zero clean, one dirty})=2
\quad\hbox{for every nonaffine quadratic }f.              \qquad(2)
```

The lower bound in (2) is a zero-clean statement. Merely allowing one
initialized qubit does not license importing that full-space argument
onto its occupied sector.

## 2. A complete native word for an even-arity symmetric quadratic

First let n be positive and even, put
$`e_2(x)=\sum_{i\lt j}x_ix_j\pmod2`$, and let
$`p=\bigoplus_i x_i`$. Every parenthesized sum in the parity expressions
below is binary. As an integer identity modulo eight,

```math
4z e_2(x)=(n-1)z+\sum_i x_i-\sum_i(z\oplus x_i)-p+(z\oplus p).
                                                               \qquad(3)
```

Indeed, $`(\sum_i x_i-p)/2\pmod2=e_2(x)`$; substituting
$`z\oplus u=z+u-2zu`$ proves (3). Replace every displayed parity l
on the right by $`d\oplus l`$ and append $`-(n-1)d`$. The sum of
the original coefficients is n-1. For d=0 the result is (3); for d=1
it is the negative of (3), which agrees modulo eight because (3) is
either zero or four. The complete phase is therefore independent of d,
including its coherent off-diagonal matrix elements.

The resulting 2n+4 odd parities form two physical bases on the n+2
wires. In physical order $`x_1,\ldots,x_n,z,d`$ they are

```math
\begin{array}{c|c|c}
 &\hbox{parities in physical slots}&\hbox{integer coefficients}\\ \hline
A&(d\oplus x_i)_i,\ d\oplus z,\ d
 &(+1)_i,\ n-1,\ 1-n\\
B&(d\oplus z\oplus x_i)_i,\ d\oplus p,\ d\oplus z\oplus p
 &(-1)_i,\ -1,\ +1.
\end{array}                                                    \qquad(4)
```

Their invertibility is certified by the following CNOT words, written
chronologically. A uses n+1 CNOTs:

```text
CX(d,x_i) for every i; CX(d,z).
```

B uses 3n+2 CNOTs:

```text
CX(z,d);
CX(d,x_i) for every i;
CX(x_i,z) for every i;
CX(d,z);
CX(x_i,d) for every i.
```

Direct propagation gives exactly the rows in (4), using n even. Each
word is an invertible map on all input wires, and its reverse is its
actual inverse. No parity is copied into a promised-zero wire.

For each basis, compute its CNOT word, apply the specified odd phase
coefficient on each wire, and reverse the CNOT word. Express an odd
coefficient as one T or T-dagger plus an even part implemented by S
powers. All S gates can precede the simultaneous odd band, since they
act on the same physical phase coordinates. Every band has n+2 distinct
physical T targets. The ordinary coefficients +1 and -1 need no S.
The two exceptional coefficients n-1 and 1-n need no S if n is zero or
two modulo eight, and a total of three S gates otherwise. Their
reversed native inverses need a total of five S gates after adjacent
power reduction.

Hadamards on z before and after the two diagonal blocks give the
literal XOR of e2. The total ledger, including both coordinate inverses
and both output Hadamards, is

```math
T=2n+4,\qquad D_T=2,\qquad
G=2(n+1)+2(3n+2)+2+s_n=8n+8+s_n.                 \qquad(5)
```

The actual reversed native inverse has the same T-count and two bands,
with $`8n+8+s'_n`$ Cliffords. The forward logical operator is also its
own inverse on the entire physical Hilbert space, since it is exactly
a Boolean XOR tensor identity on d. This stronger full-space fact,
unlike an occupied-sector identity alone, permits reusing the forward
word when that is the cheaper implementation of the actual inverse.

## 3. Reduction of any quadratic by a paid address coordinate change

The polar form of e2 on an even number r of variables is
$`J_r+I_r`$, where J is the all-ones matrix. It is nonsingular: a vector
in its kernel must equal its total parity times the all-ones vector;
for r even that vector has zero total parity and hence vanishes.

Every alternating binary form of rank r has a symplectic-pair basis
followed by a basis of its radical. Here is a constructive elimination.
Choose vectors a,b with $`B(a,b)=1`$. For every remaining vector w,
replace it by

```math
w+B(w,b)a+B(w,a)b.                                  \qquad(6)
```

This removes its pairing with a and b while preserving a basis. Recur
on the remaining subspace. If no such pair remains, that entire
subspace is the radical. Applying the same elimination to
$`J_r+I_r`$ gives an invertible r-dimensional change to the identical
standard paired form. Compose these two coordinate changes, extend by
identity on the radical, and call the resulting address map A. Then
f(x) and $`e_2((Ax)_1,\ldots,(Ax)_r)`$ have the same polar form.
Their difference is affine:

```math
f(x)=e_2((Ax)_1,\ldots,(Ax)_r)+\ell\cdot x+\beta.
                                                               \qquad(7)
```

The constant is beta=f(0); each linear coefficient is obtained by
evaluating the difference at its corresponding unit vector and
subtracting beta. These are finite classical synthesis calculations,
not supplied quantum queries or runtime state preparations.

Synthesize A by CNOT Gaussian elimination on the m physical address
wires. Chronologically apply A, the section 2 query on its first r
coordinates with the original output z and dirty d, and the actual
reverse CNOT word for A. Finally apply $`\mathrm{CX}(x_i,z)`$
for every nonzero ell_i and the literal HSSH word if beta=1. This
proves (7) as an operator equality and gives the ledger (1). The actual
inverse reverses and adjoints every emitted native gate; all phases,
affine corrections, address restoration and dirty return remain paid.

## 4. A fairly simplified symmetric cubic baseline

This component diagnoses the extra cost of using the three-layer
[e4 phase family](THREE_LAYER_NATIVE_STRUCTURE.md#84-growing-arity-does-not-remove-the-remaining-address-phase)
twice to remove its residual address phase. For an odd number m of
addresses, m at least three, set n=m-1. A CNOT change of coordinates
keeps y=(x_1,...,x_n) and replaces the last address by
$`p=\bigoplus_{i=1}^m x_i`$. The Boolean identities
$`e_3=e_1e_2`$ and $`e_1^2=e_1`$ give

```math
e_3(x)=p e_2(y).                                      \qquad(8)
```

Compute the section 2 quadratic XOR into one clean c, apply the exact
[eight-T, two-layer dirty-helper Toffoli](DIAGONAL_ENCODING_SESSION_RESULT.md#2-constructive-check-and-first-unpaid-cost)
from p,c into the original output z, and uncompute the quadratic.
All three calls reuse the same arbitrary dirty d, because each returns
it before the next call. Finally reverse the parity coordinate change.
The quadratic's full-space self-inverse property justifies using its
forward word for the actual uncomputation. Even for arbitrary c, the
whole physical word is the phase-free permutation
$`z\mapsto z\oplus p(c\oplus e_2(y))`$, leaving c, d and every
address unchanged. Thus the complete physical operator is its own
inverse; no unused sector has been discarded.

The certified ledger is

```math
T=4m+12,\qquad D_T\le6,\qquad
G=18m+16+2s_{m-1},
```

with one clean c and one arbitrary dirty d. Its actual reversed native
word instead has $`18m+16+2s'_{m-1}`$ Cliffords. The error is zero.
At m=7 this is 40 T gates and 148 Cliffords in six bands; the literal
reversed native word has 152 Cliffords. Six is an upper bound, not an
optimality claim. At m=3 the retained three-layer cubic construction
has smaller depth, and the retained relative-Toffoli baseline has
smaller counts than this generic construction.

The fixed-degree e3 family does not resolve the general depth question.
A construction in three, four or five native bands uniformly in m
would improve this displayed component frontier; the six-band e4 echo
does not do so. Neither construction supplies a constant-depth
arbitrary-table query.

## 5. Evidence and limits of reuse

The exact fixtures verify the phase identity and actual inverse on
every basis input for n=2,4,6,8,10, all coefficient residue cases and
native band/resource ledgers, every quadratic function through four
addresses at the congruence level, and literal full native matrices for
all functions through two addresses and representative higher ranks.
The higher-rank matrices include both dirty inputs and the actual
native adjoint. A seven-address cubic fixture verifies all 512 occupied
columns, including both output values and both dirty values. These
finite checks support the displayed general algebra; they are not a
substitute for it.

During the two bands, this construction freely mixes and applies T to
transformed address wires. Consequently it lies outside the read-only
source architectures. In particular, applying it separately to many
bilinear rows does not make those rows simultaneous: their original
address wires can be targets in several different parity changes, and
a single dirty helper cannot occupy several native T targets at once.
The shared-control and output batching charges in
[R91](PARALLEL_BILINEAR_WORD_QUERY.md) therefore remain unless a complete
new simultaneous schedule and allocation are supplied.
