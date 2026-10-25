# Native fusion of complete coefficient shears

[Complete shear source](WEAK_PAULI_ORBIT_COMPRESSION.md#complete-coefficient-shear) ·
[Dirty cubic-phase compiler](DIRTY_CUBIC_PHASE_COMPILER.md) ·
[Complete derivative reduction](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence)

Completed coordinate shears can be reordered with an exact phase correction.
Their physical half-bands cannot simply be regrouped: a zero-address-calibrated
word has complete probe error two. A different bipartite grouping has a paid
constant-depth native compiler, but a paired invariant subspace excludes every
word of its uniform-angle generators. These are statements about specified
source grammars. They do not bound unrestricted native T-depth, and they do
not change R96 or supply the missing complete indicator.

## 1. Complete coordinate exchanges and physical half-bands

Let $`p\ge2`$, $`m=2^p`$, and label the arbitrary coefficient-bank bits
by subsets $`S\subseteq\{0,\ldots,p-1\}`$. Write $`a_S`$ for a bit,
$`n_S=(I-Z_{A_S})/2`$ for its number projector, and

```math
f_x(a)=\sum_S a_S\prod_{i\in S}x_i\pmod2,\qquad
Q_x|a\rangle=(-1)^{f_x(a)}|a\rangle.
```

The fixed invertible Möbius bank wall turns this complete coefficient
phase into ordinary lookup. Every bank, address and helper input is arbitrary;
all equalities retain scalar phases and external references.

For each coordinate i, define

```math
\begin{aligned}
H_i&=\sum_{S:i\notin S}X_{A_S}n_{S\cup\{i\}},\\
L_i&=\exp(i\pi H_i/4),&
R_i(x_i)&=\exp(-i\pi Z_{x_i}H_i/4),\\
\Gamma_i(x_i)&=R_i(x_i)L_i.
\end{aligned}                                                     \qquad(1)
```

Each L or R has the literal one-band schedule obtained by splitting
the retained two-band shear, with its own Hadamard and CNOT walls.
In particular, $`R_i(0)=L_i^\dagger`$ and $`R_i(1)=L_i`$.
The complete source factor has action

```math
\Gamma_i|x,a\rangle
=i^{\,x_i\sum_{S:i\notin S}a_{S\cup\{i\}}}
 |x,M_i(x_i)a\rangle,\qquad M_i(x_i)=I+x_iN_i,
```

where $`(N_i a)_S=a_{S\cup\{i\}}`$ when i is absent from S and is zero
otherwise. The scalar sum is an integer sum. The maps satisfy
$`N_i^2=0`$ and $`N_iN_j=N_jN_i`$. Consequently any completed
coordinate ordering has the same binary action and transports
$`Z_{A_\varnothing}`$ to $`Q_x`$. This does not discard its phases:
they cancel in that ordering's actual source/adjoint echo.

For distinct i,j the full native unitaries obey

```math
\begin{aligned}
\Gamma_j\Gamma_i&=\Phi_{ij}\Gamma_i\Gamma_j,\\
\Phi_{ij}|x,a\rangle
&=(-1)^{x_ix_j
 \sum_{S:S\cap\{i,j\}=\varnothing}
 a_{S\cup\{i,j\}}\left(a_{S\cup\{i\}}+a_{S\cup\{j\}}\right)}
 |x,a\rangle.
\end{aligned}                                                     \qquad(2)
```

The exponent in the sign is binary. To verify (2), partition the bank
into four-cells with labels $`S,Si,Sj,Sij`$. In each cell, compare the
two integer phase sums using $`u\mathbin\oplus v=u+v-2uv`$.
The linear contributions from the top bit cancel; the remaining phase
is exactly the displayed sign. Its value is unchanged by the common
completed permutation, so it can be written on the left as in (2).
The cells have disjoint bank supports, giving the all-size identity.

The completed all-bank source retains its serial ledger
$`T=2pm`$, $`D_T\le2p`$, and $`4pm`$ elementary Clifford gates,
with no helpers. Counting T gates as well gives $`6pm`$ total native gates.
Equation (2) exchanges completed factors; it is not an interchange
rule for their individual physical T rotations.

### A paid physical half-band exchange

An actual interchange has a nontrivial correction. Set

```math
K_{ij}(x_i)=R_i(x_i)^\dagger L_j^\dagger R_i(x_i)L_j,
\qquad L_jR_i=R_iL_jK_{ij}^\dagger.                            \qquad(3)
```

This correction has a literal four-band implementation with
$`T=3m`$, no helpers, and at most $`13m+16`$ elementary native gates,
including T gates. These are common bounds for the source and actual inverse.

On a four-cell $`a=A_S,b=A_{Si},c=A_{Sj},d=A_{Sij}`$,
the fixed halves have the four rotation axes

```math
L_i:\ X_a,\ X_aZ_b,\ X_c,\ X_cZ_d;\qquad
L_j:\ X_a,\ X_aZ_c,\ X_b,\ X_bZ_d.
```

For $`x_i=1`$, write $`K_0=L_i^\dagger L_j^\dagger L_iL_j`$.
The four rotations about $`X_a`$ in this word commute with every
other displayed axis and cancel literally, including their scalar phases.
Delete the corresponding physical T or T-dagger on wire a in each
of the four native bands, retaining every Clifford parity wall.
Each band now has three distinct targets. This gives twelve T gates
and thirty-two elementary Clifford gates per cell.

Put $`P=X_bZ_c`$. Direct conjugation gives
$`PL_iP=e^{i\pi X_a/4}L_i^\dagger`$ and $`PL_jP=L_j`$.
The extra $`X_a`$ rotation is central, so
$`K_{ij}(0)=PK_0P`$. Across $`m/4`$ disjoint cells, apply their
product P conditional on $`x_i=0`$ before and after K_0.
Each such dressing uses one X on the address before and after,
plus one CNOT and one CZ per cell. With
$`X=HS^2H`$ and a three-gate CZ, both dressings cost at most
$`2m+16`$ Clifford gates. The middle costs $`11m`$ total native
gates, proving the stated simultaneous ledger. Every temporary address
flip is restored.

The correction is generally non-diagonal. The paid identity therefore
does not authorize dropping it as a phase or commutant gauge.
Its four bands also do not supply a sublogarithmic schedule for all
crossings in a growing-coordinate source.

<a id="palindromic-half-band-gap"></a>

## 2. Calibrating the zero address does not repair half-band regrouping

Consider the regrouped source

```math
W_x=R_0(x_0)\cdots R_{p-1}(x_{p-1})
             L_{p-1}\cdots L_0.                                \qquad(4)
```

Products act from right to left. Unlike an arbitrary regrouping, this
word satisfies $`W_0=I`$ exactly. Nevertheless,

```math
\left\|W Z_{A_\varnothing}W^\dagger-Q\right\|=2
\qquad\text{for every }p\ge2.                                  \qquad(5)
```

For the two-coordinate face, write
$`a=A_\varnothing,b=A_{\{0\}},c=A_{\{1\}},d=A_{\{0,1\}}`$.
Choose address $`(x_0,x_1)=(0,1)`$ and the input sector $`d=1`$.
Then $`W=L_0^\dagger\Gamma_1L_0`$. On the $`X_a=\lambda`$
fiber, with $`\lambda=\pm1`$,

```math
L_{0,\lambda}
 =e^{i\lambda\pi n_b/4}e^{i\pi X_c/4},\qquad
\Gamma_{1,\lambda}=iX_b e^{i\lambda\pi n_c/2}.
```

Substitution in $`W_\lambda=L_{0,\lambda}^\dagger\Gamma_{1,\lambda}L_{0,\lambda}`$ gives
$`W_+W_-^\dagger=-iZ_bY_c`$, and therefore

```math
WZ_aW^\dagger=-Y_aZ_bY_c,\qquad Q_x=Z_aZ_c.                    \qquad(6)
```

On the normalized state

```math
|\psi\rangle=
\frac{|00\rangle_{a,c}-|11\rangle_{a,c}}{\sqrt2}
\otimes|0\rangle_b|1\rangle_d,
```

the actual observable has eigenvalue minus one and the target has
eigenvalue plus one. Their distance is at least two, and two is the
maximum distance between unitary observables.

For larger p, choose all other address bits and all coefficients whose
labels leave this face to be zero in this witness. That subspace is
invariant: every remaining half-band has a zero upper control and acts
identically. These are allowed inputs used to disprove the complete
promise, not clean resources granted to the circuit.

Right multiplication by any J commuting with $`Z_{A_\varnothing}`$
leaves the wrong probe image unchanged. A left correction commuting with
Q also preserves its norm defect. Even an arbitrary computational-monomial
left correction retains an off-diagonal unit-modulus column from (6);
its distance from the diagonal target is at least square root two.
This last statement allows extra arbitrary helpers and address permutations.
General nonmonomial corrections remain outside these conclusions.

<a id="common-clifford-basis-gap"></a>

## 3. Commuting completed permutations need a non-Clifford joint eigenbasis

One alternative is to discard the relative phases and use the commuting,
phase-free CNOT layers $`C_i|a\rangle=|M_i(1)a\rangle`$.
Their spectral square roots commute abstractly. A common Clifford
diagonalizer, however, does not exist for any two distinct coordinates,
even after tensoring arbitrary identity helpers.

In X/Z label order, the binary symplectic matrix of C_i is

```math
F_i=\mathrm{diag}(I+N_i,I+N_i^{\mathsf T}).
```

For distinct i,j,

```math
(F_i-I)(F_j-I)
=\mathrm{diag}(N_iN_j,N_i^{\mathsf T}N_j^{\mathsf T})
\ne0.                                                        \qquad(7)
```

For example $`N_iN_j e_{\{i,j\}}=e_\varnothing`$.
Every computational-diagonal Clifford has symplectic form

```math
\begin{pmatrix}I&0\\K&I\end{pmatrix}.
```

The product of the differences from identity of any two such matrices
is zero. Simultaneous Clifford conjugation preserves whether (7) is zero,
which proves the contradiction. Identity helpers add zero blocks to
these differences and cannot change it.

This identifies the missing basis change in the proposed shortcut.
It does not prohibit a paid non-Clifford common basis, an
address-dependent construction, or another complete native source.

<a id="uniform-bipartite-mixing-gap"></a>

## 4. A paid bipartite grammar and its complete-sector limitation

For arbitrary real angles, define

```math
\begin{aligned}
E(\alpha)&=\prod_{\substack{|S|\ {\rm even}\\i\notin S}}
 e^{i\alpha X_{A_S}n_{S\cup\{i\}}},&
O(\beta)&=\prod_{\substack{|S|\ {\rm odd}\\i\notin S}}
 e^{i\beta X_{A_S}n_{S\cup\{i\}}},\\
E_x(\alpha)&=\prod_{\substack{|S|\ {\rm even}\\i\notin S}}
 e^{-i\alpha Z_{x_i}X_{A_S}n_{S\cup\{i\}}},&
O_x(\beta)&=\prod_{\substack{|S|\ {\rm odd}\\i\notin S}}
 e^{-i\beta Z_{x_i}X_{A_S}n_{S\cup\{i\}}}.
\end{aligned}                                                   \qquad(8)
```

All axes within each product commute: its X supports lie in one
parity class and its Z supports in the other. This grouping changes the
source grammar, rather than asserting that overlapping coordinate
half-bands commute.

At quarter angles these products have a paid constant-depth compiler.
Hadamards on the target parity class turn one unconditional edge into
the phase $`\omega^{b-2ab}`$, and a conditional edge into
$`\omega^{-b+2ab+2xb-4xab}`$, where
$`\omega=e^{i\pi/4}`$ and a,b are the two bank bits.
The [dirty cubic-phase compiler](DIRTY_CUBIC_PHASE_COMPILER.md)
pays every shared-control term, arbitrary helper return and actual inverse.
There are $`r=pm/4`$ edges per parity. Thus every quarter-angle
generator in (8) has constant T-depth and $`O(pm)`$ count, dirty width
and elementary native gate count. No native cost is asserted here for
arbitrary real angles.

Despite this paid schedule, every finite word V in the generators (8)
and their actual inverses obeys

```math
\left\|V Z_{A_\varnothing}V^\dagger-Q\right\|
\ge\sqrt2.                                                    \qquad(9)
```

Each occurrence may use its own angle, but that angle is uniform across
the edges in that occurrence. There are no extra bank interlayers or
address-changing operations in this theorem.

To prove (9), compare two invariant two-coordinate faces, retaining the
notation a,b,c,d of Section 2:

| Face | Address | Top coefficient | Middle coefficients | Target on a |
| --- | --- | --- | --- | --- |
| A | $`(1,0)`$ | $`d=0`$ | $`|10\rangle_{b,c}`$ | $`-Z_a`$ |
| B | $`(1,1)`$ | $`d=1`$ | $`(|10\rangle-|01\rangle)_{b,c}/\sqrt2`$ | $`+Z_a`$ |

On both faces, every $`E(\alpha)`$ and $`E_x(\alpha)`$
acts as $`e^{i\alpha X_a}`$. Every O or O_x acts as identity:
on face A the top control is zero; on face B,
$`(X_b+X_c)|\psi^-\rangle=0`$ and $`n_b+n_c=1`$.
The conditional signs agree because both address bits in face B are one,
whereas the unselected coefficient in face A is zero.
The faces are invariant for every generator and inverse.

Consequently the whole word has the same restriction
$`e^{i\theta X_a}`$ on both faces, for the signed sum theta of its
even-layer angles. Its two errors are
$`2|\cos\theta|`$ and $`2|\sin\theta|`$, whose maximum is at least
square root two. The face embedding from Section 2 gives the same
statement for every p. For the two calibrated quarter-angle words

```math
V_1=E_x O_x O E,\qquad V_2=O_x E_x E O,                         \qquad(10)
```

theta is pi/2, so the complete error is exactly two.
Both words are identity at the zero address. Neither that calibration
nor arbitrary repetitions and shared-angle tuning repair (9).

Coordinate-dependent weights, new bank interlayers and address motion
can break this paired-face premise. They require a new complete-action
calculation; they are not excluded by (9).

If an actual implementation is within epsilon of one of these ideal
words on the complete space, its transported probe is within
$`2\epsilon`$ of the ideal probe. Thus (9) becomes
$`\sqrt2-2\epsilon`$, and the calibrated words retain the lower
bound $`2-2\epsilon`$. Every helper and return error is included
in that complete-space comparison.

## 5. What parallel dirty snapshots do and do not change

Dirty replication can parallelize controlled operations on disjoint
targets, but independent snapshots do not supply the mixed coordinate
terms of the complete shear product. For example, the successive maps
$`b\gets b+xa`$ and $`d\gets d+yb`$ produce
$`d'=d+yb+xya`$. Independent evaluation on the old data misses xya.

More generally, retain a held address, fixed address-independent Clifford
encoding/decoding, and parallel branches containing at most h completed
affine-address shear queries each. Every entry of a branch's symplectic
matrix has address degree at most h. Its actual inverse has the same
bound by reversing its h queries. Direct sums and fixed Clifford walls
preserve that bound, regardless of the number or width of replicas.
The coefficient of $`A_{\{0,\ldots,p-1\}}`$ in the target Q has
degree p. Therefore h must be at least p for exact complete transport.
If h is smaller, some address has a different Pauli label, giving
operator distance at least square root two. In the narrower
computational-linear grammar the two signs are diagonal and the gap is two.

This is a bound on completed affine-query branches, not on their jointly
resynthesized native depth. Non-Clifford interlayers, address motion and
nonlinear interaction between replicas are outside it. Independent
diagonal-echo replicas likewise fall under the existing
[returned-echo theorem](RETURNED_DIAGONAL_ECHO_DEGREE.md), by tensor-grouping
their source calls into parallel slots; that observation is not a new
degree theorem.

## 6. Consequence for the complete source question

The native phase-exchange law, paid half-band correction and dirty cubic compiler are reusable
components. The regrouped-coordinate and uniform bipartite words do
not provide a shallow lookup source, and the common-Clifford-basis and
independent-replica shortcuts retain the stated limitations.

The next unpaid operation is a complete native transition that
distinguishes the two invariant lookup sectors, for example through a
specified coordinate-dependent interaction or changed address motion.
Its complete Pauli derivative, actual inverse, uniform requested accuracy,
count, width and return must be proved together. Right commutant freedom
does not correct an already wrong derivative.

The [two-call reduction](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence)
and [polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
remain available for a qualifying source. The [R98 transfer](QUERY_DEPTH_FRAME_TRANSFER.md)
still has its separate log-log term and fixed-precision scope.
No unconditional general frame bound changes here.

Exact [coordinate-half fixtures](../../tests/test_native_shear_fusion.py),
[bipartite fixtures](../../tests/test_bipartite_shear_mixing.py) and
[dirty-phase fixtures](../../tests/test_dirty_cubic_phase_compiler.py)
check literal native phases, actual adjoints, arbitrary dirty columns and
the invariant-sector witnesses. The all-size identities and scopes are
established by the arguments above.
