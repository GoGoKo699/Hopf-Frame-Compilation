# Coordinate-weighted bipartite pulses

[Uniform-angle obstruction](COEFFICIENT_SHEAR_NATIVE_FUSION.md#uniform-bipartite-mixing-gap) ·
[Dirty cubic-phase compiler](DIRTY_CUBIC_PHASE_COMPILER.md) ·
[Complete derivative contract](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence)

Coordinate-dependent quarter angles retain a paid thirteen-band compiler.
They do not remove the need to distinguish coordinates across an entire
pulse word: two equal angle signatures retain the paired-sector obstruction,
and sufficiently close signatures give a quantitative version. This bounds
the number of specified pulses, not native T-depth. In particular, useful
coordinate-dependent Clifford pulses cost no T-depth and remain available
as bank interlayers. No general frame bound changes.

## 1. Pulse grammar and a paid weighted compiler

Let $`p\ge2`$, $`m=2^p`$, and label arbitrary coefficient-bank bits by
subsets S of the p address coordinates. Put $`n_S=(I-Z_{A_S})/2`$.
For parity $`\rho\in\{0,1\}`$ and a real vector theta, define

```math
\begin{aligned}
P_{\rho}(\theta)&=
 \prod_{\substack{|S|\equiv\rho\pmod2\\i\notin S}}
 e^{i\theta_i X_{A_S}n_{S\cup\{i\}}},\\
P_{\rho,x}(\theta)&=
 \prod_{\substack{|S|\equiv\rho\pmod2\\i\notin S}}
 e^{-i\theta_i Z_{x_i}X_{A_S}n_{S\cup\{i\}}}.
\end{aligned}                                                     \qquad(1)
```

The angle depends on the coordinate and pulse occurrence, but not on the
tail S. Every factor within a pulse commutes. There are
$`r=pm/4`$ edges of each parity. Each angle has exact period $`2\pi`$:
the edge generator has eigenvalues zero, plus one and minus one. This
periodicity retains the scalar phase on every logical column.

Every pulse whose angles are integer multiples of $`\pi/4`$ has common
forward and actual-inverse bounds

```math
\begin{aligned}
T&\le m+32r,&D_T&\le13,&w_{\rm dirty}&\le5r,\\
G_{\rm all}&\le6m+71r,&a_{\rm clean}&=0,&\epsilon&=0.
\end{aligned}                                                     \qquad(2)
```

Here $`G_{\rm all}`$ counts every elementary H, S, CX, T and T-dagger gate.
The p address bits and m bank bits are excluded from the helper width.
The address is preserved and every arbitrary helper returns, including
coherent inputs correlated with any external reference.

To prove (2), reduce each angle to $`k_i\pi/4`$ with
$`k_i\in\{0,\ldots,7\}`$, and write $`k_i=2\ell_i+e_i`$,
where $`0\le\ell_i\le3`$ and $`e_i\in\{0,1\}`$.
For an edge with tail a and head b, its even-angle part is Clifford:

```math
\begin{aligned}
e^{i\pi X_a n_b/2}
 &=S_b\,\mathrm{CX}_{b\to a},\\
e^{-i\pi Z_{x_i}X_a n_b/2}
 &=S_b^\dagger\,\mathrm{CX}_{b\to a}\,\mathrm{CZ}_{x_i,b}.
\end{aligned}                                                     \qquad(3)
```

For b zero both sides act identically; for b one they act as
$`iX_a`$ or $`-i(-1)^{x_i}X_a`$, respectively. Raise (3) to
the power $`\ell_i`$. Reduce the adjacent S block modulo four and the CX/CZ
blocks modulo two. At most three S gates, one CX and one three-gate CZ
are needed per edge, hence at most $`7r`$ Clifford gates.
The same bound holds for the literal inverse after reducing its S block.

The remaining odd-angle subset is a subset of the unit quarter-angle
edges. The [cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md#2-a-paid-parity-class-of-coefficient-cube-edges)
uses at most m active linear variables, $`Q\le2r`$ quadratic occurrences
and $`C\le r`$ cubic occurrences. Its bounds are
$`T\le m+32r`$, thirteen bands, $`w\le5r`$ and
$`G_{\rm all}\le6m+64r`$. Removing edges does not change the shared-control
scheduling proof. All even and odd edge factors commute within this parity,
so prepend the Clifford factors and add their gate charge. This proves (2)
without a clean assumption or an unpriced phase.

For a word of t such pulses the same helper pool can be reused after each
pulse returns it. Thus T and $`G_{\rm all}`$ multiply by t, $`D_T\le13t`$, and
the helper width remains at most $`5r`$. These are upper bounds for a
specified emitter; they are not lower bounds on another native synthesis.
No native synthesis cost is asserted for arbitrary real angles.

## 2. Equal and nearby coordinate signatures

Let $`V=P_t\cdots P_1`$ be a word of t pulses of the four types in
(1), with no additional bank interlayers or address motion. Each occurrence
may choose its parity, signed or unsigned type, and coordinate angles
independently. Actual inverses are included by negating their angles.
The complete target is

```math
Q|x,a\rangle=(-1)^{\sum_S a_S\prod_{i\in S}x_i}|x,a\rangle,
\qquad
\mathcal E(V)=\|VZ_{A_\varnothing}V^\dagger-Q\|.
```

Define the coordinate signature and circular distance by

```math
\sigma_i=(\theta_{1,i},\ldots,\theta_{t,i})
 \in(\mathbb R/2\pi\mathbb Z)^t,\qquad
d_1(\sigma_i,\sigma_j)=
 \sum_{k=1}^t\min_{z\in\mathbb Z}
       |\theta_{k,i}-\theta_{k,j}+2\pi z|.
```

For every distinct pair i,j,

```math
\mathcal E(V)\ge\sqrt2-2d_1(\sigma_i,\sigma_j).
                                                               \qquad(4)
```

In particular, any equal pair of signatures forces error at least
square root two, regardless of the number, order or angles of the pulses.
Distinct signatures are necessary for error below that threshold; they
are not sufficient for correct lookup.

To prove (4), restrict to the invariant four-coefficient face
$`a=A_\varnothing,b=A_{\{i\}},c=A_{\{j\}},d=A_{\{i,j\}}`$.
Set every outside address and every coefficient whose label leaves this
face to zero in the witness subspace. An exterior edge has zero head
control and acts identically. This is a choice of allowed inputs, not
a grant of clean qubits to the source.

First suppose the two signatures agree. On this face compare the same
two invariant sectors used by the uniform-angle theorem:

| Sector | Address | d | State of b,c | Target on a |
| --- | --- | --- | --- | --- |
| A | $`(x_i,x_j)=(1,0)`$ | zero | $`\lvert10\rangle`$ | $`-Z_a`$ |
| B | $`(x_i,x_j)=(1,1)`$ | one | $`(\lvert10\rangle-\lvert01\rangle)/\sqrt2`$ | $`+Z_a`$ |

An even pulse, signed or unsigned, restricts on both sectors to
$`e^{i\theta_{k,i}X_a}`$. An odd pulse is identity on both: its top
control vanishes in A, while $`(X_b+X_c)`$ annihilates the singlet
in B. Thus V restricts to the same $`e^{i\Theta X_a}`$ in both.
Its two target errors are $`2|\cos\Theta|`$ and
$`2|\sin\Theta|`$, whose maximum is at least square root two.

For unequal signatures, form $`V'`$ by replacing coordinate j's angle
by coordinate i's angle in every pulse. On the invariant four-coefficient
face this changes exactly one commuting edge factor per pulse: a from c
in an even pulse, or b from d in an odd pulse. With the shortest circular
angle difference $`\delta_k`$, its operator difference obeys

```math
\|P_k-P'_k\|_{\rm face}
 \le\|e^{i\delta_k X n}-I\|
 =2|\sin(\delta_k/2)|\le|\delta_k|.
                                                               \qquad(5)
```

The addressed sign leaves the same bound. Telescoping the product gives
$`\|V-V'\|_{\rm face}\le d_1(\sigma_i,\sigma_j)`$.
Conjugation of a norm-one probe changes by at most twice this distance.
Apply the equal-signature obstruction to $`V'`$ on the same face to obtain
(4). No claim about the distance between V and $`V'`$ on the entire larger
bank is needed.

If a native implementation $`\widetilde V`$ obeys the complete operator
bound $`\|\widetilde V-V\otimes I_{\rm work}\|\le\varepsilon`$,
then its target error is at least
$`\sqrt2-2d_1(\sigma_i,\sigma_j)-2\varepsilon`$.
This comparison includes all helper return errors and arbitrary references.

## 3. Pulse-count consequences

Suppose every occurrence draws its coordinate angles from an alphabet of
at most K values modulo $`2\pi`$. At most $`K^t`$ signatures exist, so

```math
\mathcal E(V)\lt\sqrt2\quad\Longrightarrow\quad p\le K^t.
                                                               \qquad(6)
```

For quarter-angle words the alphabet has eight values, including all
inverse angles. Consequently $`t\ge\lceil\log_8 p\rceil`$ is necessary.
The weighted compiler in Section 1 pays these pulses, but its thirteen-band
upper bound does not turn (6) into a native-depth lower bound.

For arbitrary real angles and fixed $`0\le\tau\lt\sqrt2`$, put
$`\Delta=\sqrt2-\tau`$. For $`t\ge1`$ let
$`q=\lfloor4\pi t/\Delta\rfloor+1`$. Partition the angle circle into
q half-open intervals of length $`2\pi/q`$. If $`p>q^t`$, two signatures
fall in the same product cell, giving
$`d_1\lt2\pi t/q\lt\Delta/2`$ and contradicting (4). Therefore

```math
\mathcal E(V)\le\tau\quad\Longrightarrow\quad
p\le\left(\left\lfloor\frac{4\pi t}{\sqrt2-\tau}\right\rfloor+1\right)^t,
\qquad
t=\Omega_\tau\!\left(\frac{\log p}{\log\log p}\right)
\quad\text{as }p\to\infty.                                   \qquad(7)
```

For $`t=0`$, V is identity and its complete target error is two: choose
an address with one selected coordinate and set that singleton coefficient
to one. A fixed number of real-angle pulses is therefore insufficient for
arbitrarily many coordinates at any target error below square root two.
The real-angle estimate does not exclude every sublogarithmic pulse count.

## 4. The free-Clifford escape must remain open

Equation (3) prevents charging every pulse against native T-depth. On the
two-coordinate face, an odd pulse selecting only coordinate i with angle
$`\pi/2`$ is identity in sector A and acts as $`iX_c`$ in sector B.
It sends the singlet to
$`i(\lvert11\rangle-\lvert00\rangle)/\sqrt2`$ and hence breaks the
common invariant-sector calculation. This operation is a paid Clifford,
with no helper and zero T-depth.

Accordingly, a growing number of Clifford pulses interleaved with a few
non-Clifford pulses is not excluded by (6) or (7). The total-pulse signature
can grow while the T-depth charge does not. No normal form eliminating
those Clifford interlayers has been proved here. Arbitrary edge-dependent
weights, which can depend on S as well as i, also lie outside (1), as do
additional bank transformations and address motion.

The useful next source question is whether a specified sequence of paid
Clifford bank transformations and a few nonlinear pulses can supply the
complete lookup derivative. Its all-size action, actual inverse, accuracy,
count, width and work return remain necessary. Signature separation alone
is not a construction. The complete-indicator and frame frontiers remain
unchanged.

The [exact fixtures](../../tests/test_coordinate_weighted_pulses.py) check
weighted native phases, actual inverses and arbitrary helper return, the
embedded signature collision, and the Clifford escape from a chosen paired
sector. The all-size bounds follow from the identities and packing argument
above, rather than from the finite fixtures.
