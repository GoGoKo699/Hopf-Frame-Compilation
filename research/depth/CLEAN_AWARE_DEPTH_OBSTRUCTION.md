# Constant T-depth obstructions with initialized clean work

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Full-space two-layer lemma](SHALLOW_SOURCE_OBSTRUCTION.md) · [Frame convention](../../docs/HOPF_INTERFACE.md)

Arbitrarily many dirty qubits do not make three T layers sufficient for
arbitrarily accurate complete-frame compilation when the number of
initialized clean qubits is fixed. This statement allows arbitrary
Clifford interlayers and does not restrict the compiler to a source,
lookup, or frame-factorization architecture.

**Theorem.** For each fixed integer $`a\ge0`$, there is a constant
$`\varepsilon_a>0`$ such that, for every $`n\ge1`$ and dirty width
$`b\ge0`$, the complete real frame

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3),\qquad R_y(\theta)=e^{-i\theta Y},
```

has initialized-isometry distance at least $`\varepsilon_a`$ from
every coherent Clifford+T circuit of T-depth at most three using at
most a initialized clean qubits. Accuracy strictly below
$`\varepsilon_a`$ therefore requires at least four T layers,
independently of b. Zero leaf phases give the same witness in the
phase-dressed complex magnitude family.

For positive clean budgets, Section 4 establishes the constant by
compactness. The elementary [initialized-sector proof in R74,
Section 5.1](FOUR_LAYER_CROSSING_RANK.md#51-quantitative-rigidity-on-the-initialized-sector)
also gives an explicit, conservative threshold without the imported
polynomial-limit premise. With no clean qubits, Section 4.2 sharpens
the same witness to the exact error floor

```math
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}\approx0.169714,
```

so integer precision $`L\ge3`$ requires at least four T layers.
The separate [four-layer operator theorem](FOUR_LAYER_OPERATOR_RIGIDITY.md),
R75, strengthens this zero-clean conclusion to at least five layers at
the same precision, with arbitrary dirty width.
A separate explicit one-clean witness with $`\cos(2\theta)=13/16`$
has distance at least

```math
c_* = \frac{3-2\sqrt2}{32}>2^{-8}
```

from every circuit of T-depth at most two. Thus one-clean compilation
at integer precision $`L\ge8`$ requires at least three T layers in
the worst case. These are constant lower bounds; they do not grow with
n or L. Optimal depth remains an open problem.

## 1. Initialized work and a robust transfer witness

Let S be the n-qubit logical register, D the b-qubit arbitrary dirty
register, and C the a-qubit initialized register. Put $`d=n+b`$,
$`q=d+a`$, and

```math
J|\psi\rangle_{SD}=|\psi\rangle_{SD}|0^a\rangle_C.
```

For logical Hermitian Paulis P,Q and an actual circuit V on all q wires,
define

```math
x_V(P,Q)=2^{-d}\mathrm{Tr}_{SD}\!\left[
(Q\otimes I_D)J^\dagger V^\dagger
(P\otimes I_{DC})VJ\right].
```

This is a logical-channel transfer coefficient with the dirty input
averaged over its maximally mixed state. The average is only a proof
witness: the compiler must still work for every dirty input and reference.
Expanding the clean projector gives

```math
|0^a\rangle\langle0^a|=2^{-a}
\sum_{R\in\{I,Z\}^{\otimes a}}R,
```

```math
x_V(P,Q)=\sum_{R\in\{I,Z\}^{\otimes a}}
2^{-q}\mathrm{Tr}\!\left[
(P\otimes I_{DC})V(Q\otimes I_D\otimes R)V^\dagger\right].
```

There are $`2^a`$ full-space transfer entries and no additional
prefactor: the clean-projector normalization is already in $`2^{-q}`$.

Suppose $`\|VJ-J(U\otimes I_D)\|\le\eta`$. Set
$`F=VJ`$, $`G=J(U\otimes I_D)`$, and $`O=P\otimes I_{DC}`$.
Both F and G are isometries, so

```math
\|F^\dagger OF-G^\dagger OG\|
\le\|(F-G)^\dagger OF\|+\|G^\dagger O(F-G)\|
\le2\eta.
```

Normalized trace against $`Q\otimes I_D`$ then gives

```math
\left|x_V(P,Q)-2^{-n}\mathrm{Tr}(PUQU^\dagger)\right|
\le2\eta.
```

Only initialized columns have been used; no full-unitary approximation
on occupied C is assumed. The norm includes clean leakage and dirty
return, and is unchanged after tensoring a reference identity.
Optimizing a common global phase cannot evade the transfer condition.
Fewer clean qubits are covered by appending untouched initialized
spectators up to the stated budget.

## 2. An explicit one-clean two-layer gap

The [full-space lemma](SHALLOW_SOURCE_OBSTRUCTION.md#1-a-discrete-pauli-transfer-alphabet-at-t-depth-at-most-two)
places every transfer entry of a circuit with at most two T layers in

```math
\mathcal A=\{0\}\cup\{\pm2^{-j/2}:j\ge0\}.
```

It uses Pauli conjugation by a T layer into a Hermitian Clifford and
stabilizer-overlap quantization, with arbitrary Clifford interlayers.
Section 1 therefore gives $`x_V(P,Q)\in\mathcal A+\mathcal A`$
for one clean qubit. The exact distance

```math
\mathrm{dist}(13/16,\mathcal A+\mathcal A)
=\frac{3-2\sqrt2}{16}
```

can be checked without an infinite enumeration. Put $`s=\sqrt2`$
and order the two summands $`u\ge v`$. If $`u\le s/4`$, their
sum is at most $`s/2`$, too far below $`13/16`$. The other possible
largest summands and their nearest candidates are

| u | Desired v | Adjacent alphabet points |
|---|---|---|
| $`1/2`$ | $`5/16`$ | $`1/4,\ s/4`$ |
| $`s/2`$ | $`(13-8s)/16`$ | $`s/16,\ 1/8`$ |
| $`1`$ | $`-3/16`$ | $`-1/4,\ -s/8`$ |

The first two rows have distances strictly greater than
$`(3-2s)/16`$; the relevant comparisons include
$`10>7s`$ and $`10s>14`$. The last row attains equality at
$`u+v=1-s/8`$. Negative and zero summands are included.

For every n, set all frame angles above depth n−1 to zero and every
angle at depth n−1 to $`\theta=\tfrac12\arccos(13/16)`$.
The earlier layers are literal identities. The last layer has an empty
suffix and the same rotation at every prefix, hence its complete operator is

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta).
```

This holds on every column, including at the allowed chart-boundary
angle tuples, and covers n=1. With $`P=Q=I_{2^{n-1}}\otimes Z`$,
its transfer entry is $`\cos(2\theta)=13/16`$. Section 1 gives
$`\eta\ge c_*`$, and $`2^{-8}\lt c_*`$ proves the explicit precision
statement. The witness is not claimed to optimize all possible gaps.

## 3. Three layers as a stabilizer overlap with one T layer

Let $`V=C_3D_3C_2D_2C_1D_1C_0`$, with arbitrary Clifford blocks
and physical T/T-dagger layers D. Absorb the exterior Cliffords into
the witness Paulis. Cyclicity of trace puts a full-space transfer entry
in the form

```math
2^{-q}\mathrm{Tr}(A D_2 B D_2^\dagger),
```

where

```math
A=C_2^\dagger(D_3^\dagger P'D_3)C_2,\qquad
B=C_1(D_1Q'D_1^\dagger)C_1^\dagger.
```

A and B are Hermitian Clifford unitaries, by the same one-layer
conjugation fact used for two layers. For normalized vectorization

```math
|K\rangle\!\rangle=2^{-q/2}\sum_{x,y}K_{xy}|x,y\rangle,
```

the entry is

```math
\langle\!\langle A|
(D_2\otimes\overline{D_2})|B\rangle\!\rangle.
```

Both boundary vectors are stabilizer states on 2q qubits. The middle
operator is one parallel layer of T and T-dagger gates.

We use the affine-support quadratic-phase description of stabilizer
states from [Dehaene–De Moor, arXiv:quant-ph/0304125v1,
Section V](https://arxiv.org/pdf/quant-ph/0304125). Each boundary
state is flat on an affine binary support, with relative phases given
by a linear phase modulo four and a quadratic sign. Its common phase
is retained: here it is an eighth root of unity. Indeed the vectorized
matrix is Hermitian. If one nonzero entry has phase $`\alpha i^r`$,
its transposed entry has phase $`\alpha i^t`$ and equals its complex
conjugate. Hence $`\alpha^2=i^{-r-t}`$ and $`\alpha^8=1`$.

If the supports are disjoint the entry is zero. Otherwise let their
dimensions be $`r_A,r_B`$ and their intersection dimension k.
Parameterize the intersection by $`u\in\mathbb F_2^k`$. The overlap is

```math
x=2^{-s/2}\mathbb E_u\,\zeta_8^{p(u)},\qquad
s=r_A+r_B-2k\ge0,\qquad \zeta_8=e^{i\pi/4}.
```

Here p is a multilinear polynomial modulo eight of degree at most three,
including the constant relative phase. To see this, expand the affine
parities in the T phases modulo eight:

```math
\bigoplus_{j\in S}u_j
\equiv\sum_{j\in S}u_j
-2\sum_{i\lt j\in S}u_i u_j
+4\sum_{i\lt j\lt\ell\in S}u_i u_j u_\ell\pmod8.
```

Complementing a parity changes only constants and signs. Stabilizer
linear phases have coefficient two and quadratic signs coefficient four;
their affine restrictions also have degree at most three modulo eight.
Crucially, the prefactor is at most one, independently of width.

## 4. Polynomial limits and the closed coefficient set

The additional imported premise is [Khodabandeh–Shinkar,
*On Sampling Lower Bounds for Polynomials*, arXiv:2605.00995v1,
Appendix A, Lemma A.6](https://arxiv.org/html/2605.00995v1#A1):
if output probabilities of Boolean polynomials with a fixed degree bound
converge, their limits are dyadic. The number of input bits may vary.
We use this limit statement, not only their gap theorem for a single
specified probability. No quantitative separation constant is imported.

Every event $`[p(u)=r]`$, $`r\in\{0,\ldots,7\}`$, in Section 3
is a Boolean polynomial of degree at most 21. For a direct sufficient
bound, represent each coefficient modulo eight by an integer from zero
to seven and repeat its monomial that many times. This writes p modulo
eight as an integer sum of Boolean tokens $`t_1,\ldots,t_M`$ of
degree at most three. Its three binary digits are, over $`\mathbb F_2`$,

```math
b_0=e_1(t),\qquad b_1=e_2(t),\qquad b_2=e_4(t),
```

where $`e_j`$ is the jth elementary symmetric polynomial. Comparing
coefficients of $`z^{2^j}`$ in

```math
\prod_\ell(1+t_\ell z)=(1+z)^{\sum_\ell t_\ell}
```

over $`\mathbb F_2`$ proves these identities pointwise. The digit
degrees are at most 3,6,12, so the product of their three required-value
indicators has degree at most 21. This conservative bound suffices.

Let $`\mathcal S_3`$ contain all full-space transfer entries of circuits
with at most three T layers, for every finite width and every Hermitian
Pauli pair, and put $`\mathcal K_3=\overline{\mathcal S_3}`$.
It is compact in $`[-1,1]`$. We claim

```math
\mathcal K_3\subseteq\mathbb Z[1/\sqrt2].
```

Take a convergent sequence of such coefficients with limit x. If x=0
there is nothing to prove. Otherwise the bound $`|x_j|\le2^{-s_j/2}`$
forces $`s_j`$ to be bounded along a tail. Pass to a subsequence
with fixed s and then to one on which all eight probabilities

```math
\pi_r=\Pr_u[p(u)=r]
```

converge. Each limiting probability is dyadic by Lemma A.6 applied to
its degree-at-most-21 event. Consequently

```math
x=2^{-s/2}\sum_{r=0}^7\pi_r^{\rm lim}\zeta_8^r
\in\mathbb Z[1/\sqrt2,i].
```

Transfer entries are real, so x lies in the real part of this ring,
$`\mathbb Z[1/\sqrt2]`$. This restricts limits across unbounded
widths. Membership of each finite circuit in a dyadic ring alone would
not imply it.

For fixed a, define the compact finite Minkowski sum

```math
\mathcal B_{3,a}=
\underbrace{\mathcal K_3+\cdots+\mathcal K_3}_{2^a\ \text{terms}}.
```

It remains a subset of $`\mathbb Z[1/\sqrt2]`$. Every element of
that ring is $`u+v\sqrt2`$ with dyadic u,v; its rational elements
are therefore dyadic. Thus $`1/3\notin\mathcal B_{3,a}`$, and
compactness gives

```math
\varepsilon_a=\tfrac12\mathrm{dist}(1/3,\mathcal B_{3,a})>0.
```

Section 1 places every initialized-channel coefficient in this sum.
Use Section 2's complete-frame tuple with
$`\theta=\tfrac12\arccos(1/3)`$. Its coefficient is 1/3, so its
initialized-isometry error is at least $`\varepsilon_a`$.
This proves the theorem uniformly over n and b. The constant may depend
on a; no uniform bound is asserted when a grows.

### 4.1 A sharp gap below a preserved Pauli

The phase polynomial in Section 3 has the stronger weighted form

```math
p(u)=c+\sum_i\ell_i u_i
+2\sum_{i\lt j}q_{ij}u_i u_j
+4\sum_{i\lt j\lt k}t_{ijk}u_i u_j u_k\pmod8.
```

The parity expansion, the factor two in stabilizer linear phases,
and the factor four in quadratic signs establish these divisibilities.
They survive affine restriction. Put

```math
\delta=\frac{1-1/\sqrt2}{2}.
```

Every full-space transfer entry x of a circuit with at most three
T layers obeys the sharp alternative

```math
x=1\quad\text{or}\quad
x\le1-\delta=\frac{1+1/\sqrt2}{2}.
```

Zero entries satisfy this alternative directly. For the normal form
$`x=2^{-s/2}\mathbb E\zeta_8^{p}`$, the case
$`s\ge1`$ gives the stronger bound $`x\le1/\sqrt2`$. For s=0,
use $`1-x=\mathbb E[1-\cos(\pi p/4)]`$. If p modulo two is
nonzero, it is affine and is one on at least half the inputs. Odd
residues contribute loss at least $`1-1/\sqrt2`$, proving the bound.
Otherwise all multilinear coefficients are even. If p/2 modulo two
is nonzero, it has degree at most two and is one with probability
at least 1/4; those inputs have $`p\equiv2\pmod4`$ and contribute
loss one. Otherwise p/4 is a Boolean polynomial of degree at most
three; if nonzero it is one with probability at least 1/8, contributing
loss two. Both latter cases give loss at least 1/4. If all three
reductions vanish, p is identically zero and x=1.

The elementary support fact used here is that a nonzero Boolean
polynomial of degree d is nonzero on at least a $`2^{-d}`$ fraction
of the cube. Choose a maximal-degree monomial S: its iterated difference
along S is the constant one, so every S-subcube contains a nonzero
value. Multilinear coefficients are unique on the Boolean cube, which
justifies division by two when the preceding parity vanishes. The
native family in Section 5 at m=1 attains $`x=1-\delta`$.
This gap does not use the polynomial-limit theorem.

Suppose a target full-space unitary U preserves a Hermitian Pauli P
and $`\|V-U\|\le\eta`$. With the normalized Hilbert--Schmidt norm,

```math
1-x_V(P,P)
=\tfrac12\|V^\dagger PV-P\|_{2,\mathrm{tr}}^2
\le2\eta^2.
```

Consequently, for

```math
\eta\lt\rho,\qquad
\rho=\sqrt{\delta/2}
=\frac{\sin(\pi/8)}{\sqrt2}\approx0.270598,
```

the circuit preserves P exactly. This applies simultaneously to any
number of target-preserved Paulis, without a dimension or union-bound
factor. A common global phase does not change the conclusion.

### 4.2 Exact dirty return and the explicit zero-clean floor

Take the theorem's target with $`\theta=\tfrac12\arccos(1/3)`$
and no initialized clean work. Its initialized-isometry condition is
then a full-space operator condition. It preserves Y on the final
logical wire and X,Z on every spectator logical and dirty wire.
At error below $`\rho`$, Section 4.1 forces V to preserve all these
operators exactly. Commutation with the spectator Pauli algebra gives

```math
V=e^{i\varphi}R_y(\alpha)\otimes I_{\mathrm{spectators}}
```

with the tensor factors ordered as in the target. The two transfer
entries $`c=\cos(2\alpha)`$, $`s=\sin(2\alpha)`$ lie in
$`\mathbb Z[1/\sqrt2]`$ and satisfy $`c^2+s^2=1`$. These constraints
force $`c+is`$ to be an eighth root of unity, hence
$`\alpha\in(\pi/8)\mathbb Z`$ modulo a common phase.

Here is an elementary proof of this unit-circle classification. Write

```math
c+is=\frac{(a+b\sqrt2)+i(e+f\sqrt2)}{2^{k/2}}
```

with integer coefficients. Unit norm gives

```math
a^2+e^2+2(b^2+f^2)=2^k,\qquad ab+ef=0.
```

The integer vectors $`(a,e)`$ and $`(b,f)`$ are perpendicular.
Choose a primitive integer direction $`(p,q)`$ so these vectors
are $`g(p,q)`$ and $`h(-q,p)`$, allowing g or h to vanish. Then

```math
(p^2+q^2)(g^2+2h^2)=2^k.
```

Primitivity and parity force $`p^2+q^2`$ to be one or two. Also,
$`g^2+2h^2=2^j`$ forces g or h to vanish: for $`j\ge2`$ both
are even, so divide by four until the immediate cases j=0 or j=1.
Normalizing the remaining axial or diagonal vector gives exactly the
eight roots. Using transfer entries keeps this argument independent
of the circuit's common phase.

The nearest allowed angle to $`\theta`$ is $`\pi/4`$. Even after
optimizing a common phase, its operator distance is

```math
c_0=2\sin\!\left(\frac{\pi/4-\theta}{2}\right)
=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}},\qquad
\frac18\lt c_0\lt\rho.
```

If a circuit had error below $`c_0`$, the preceding rigidity argument
would quantize its angle and contradict this nearest-angle distance.
Conversely the literal Clifford $`R_y(\pi/4)=HZ`$, tensored with
identity on all other wires, attains $`c_0`$ on every column.
Thus $`c_0`$ is the exact best error at at most three T layers in the
zero-clean, arbitrary-dirty model. Every integer $`L\ge3`$ requires
at least four layers for this complete-frame witness. This stronger
explicit statement does not extend to a positive clean budget by
ignoring its initialized-sector restriction.

### 4.3 One clean qubit: exact preservation of the terminal sectors

For one initialized clean qubit, write the circuit in its two clean-bit
blocks as

```math
V=\begin{pmatrix}A&B\\C&D\end{pmatrix}.
```

The initialized-isometry error bound gives $`\|C\|\le\eta`$.
Unitarity implies $`\|B\|=\|C\|`$: the matrices
$`BB^\dagger=I-AA^\dagger`$ and $`C^\dagger C=I-A^\dagger A`$
have the same eigenvalues. Hence

```math
\|[V,Z_C]\|\le2\eta,\qquad
1-x_V(Z_C,Z_C)\le2\eta^2.
```

At at most three T layers and $`\eta\lt\rho`$, Section 4.1 forces
$`[V,Z_C]=0`$ exactly. Both terminal clean-bit sectors are therefore
invariant, even though approximation was required only on the zero
sector. This does not assert intermediate clean return or control the
logical unitary on the occupied sector. For example a CNOT from the
clean wire to a logical wire is identity on the initialized sector,
preserves $`Z_C`$, and does not preserve the logical Y globally.
A controlled Pauli is Clifford rather than Pauli; applying Section 4.1
to it would require an additional argument.

## 5. Why the two-layer alphabet does not extend

Three-layer transfer coefficients can have nonzero accumulation points.
The following exact native family also rules out replacing their set by
$`\mathcal A+\mathcal A`$. It is compatible with Section 4's stronger
restriction on limits.

Use arbitrary-input wires $`s,a,b_1,\ldots,b_m`$, with $`m\ge1`$, and put

```math
X_B=\prod_{j=1}^m X_{b_j},\qquad
F_m=|0\rangle\langle0|_a\otimes I_B+
|1\rangle\langle1|_a\otimes X_B
=\frac{I+Z_a+X_B-Z_aX_B}{2}.
```

The controlled fanout $`F_m`$ is Clifford. The three Pauli axes

```math
R_1=Y_sZ_a,\qquad R_2=Y_sX_B,\qquad R_3=Y_sZ_aX_B
```

commute and are independent: the exponents of $`Y_s,Z_a,X_B`$ in a
product force all three axis exponents to vanish. Each anticommutes
with $`P=X_s`$. Define

```math
S_m=e^{-i\pi R_1/8}e^{-i\pi R_2/8}e^{+i\pi R_3/8}.
```

A Clifford maps the three independent commuting axes to Z on distinct
wires. Thus $`S_m`$ is one physical T layer, containing two T gates and
one T-dagger, conjugated by Cliffords. Its scalar phase cancels in the
conjugations below. No initialized wire is required. Direct expansion gives

```math
(I+iR_1)(I+iR_2)(I-iR_3)
=(I+iY_s)(I+Z_a+X_B-Z_aX_B),
```

```math
S_mPS_m^\dagger
=\frac{P(I+iR_1)(I+iR_2)(I-iR_3)}{\sqrt8}
=\frac{X_s-Z_s}{\sqrt2}F_m=:A_m.
```

Let $`D_m=T^{\otimes m}`$ act on B and set $`V_m=S_m^\dagger D_mS_m`$.
This is a native circuit with at most three T layers and $`m+6`$ T gates.
Cyclicity of trace and cancellation of the s reflection give

```math
\widehat V_{m;P,P}
=2^{-m-2}\mathrm{Tr}(A_mD_mA_mD_m^\dagger)
=\frac{1+2^{-m/2}}2.
```

Indeed the $`a=0`$ block is identity; the $`a=1`$ block factors into m
copies of $`XTXT^\dagger`$, each with normalized trace $`1/\sqrt2`$.
These distinct coefficients converge to $`1/2`$. Moreover, disjoint
tensor products preserve T-depth and multiply transfer coefficients, so
$`V_4\otimes V_4`$ has coefficient $`(5/8)^2=25/64`$.

This rational number is outside $`\mathcal A+\mathcal A`$. A nonzero
rational sum of two alphabet entries must be a signed sum of two dyadic
powers: irrational terms either remain irrational or cancel to zero.
After scaling by 64, two unequal signed powers representing the odd
integer 25 would have their smaller power equal to one; the other would
then be 24 or 26, neither a power of two. Equal powers also cannot give 25.

All wires here have arbitrary inputs. The construction makes no claim
that a designated dirty subset returns to identity while a prescribed
logical gate is performed. It establishes a limitation of simple
transfer alphabets, not a shallow complete-frame compiler or interval
density of three-layer coefficients.

## 6. Scope and the next obstruction

The two-layer bound is explicit at one clean qubit. The three-layer
bound covers every fixed clean budget, including the one- and two-clean
frame theorems. Its zero-clean threshold is explicit and sharp for the
stated witness; R74 now supplies an explicit elementary threshold for
positive clean budgets as well.
Neither assumes a bound on T-count or Clifford count. Both permit all
intermediate mixing with dirty work and cover the reference-stable
return contract.

The proof does not iterate automatically. Peeling the outer T layers
from a four-layer circuit leaves two middle layers separated by a
Clifford block. A path expansion of that block can introduce a
normalization factor larger than one, depending on width. Small
polynomial biases multiplied by growing factors are not controlled by
Section 4's limiting-probability theorem. Bounding that factor, or
replacing the scalar invariant, is an additional problem. R74 controls
the middle crossing rank and gives a quantitative initialized-sector
gap. R75 replaces the scalar invariant by operator norm, proving the
unrestricted zero-clean five-layer obstruction. The corresponding
positive-clean five-layer question remains open.

A fixed lower bound of four implies no growth with dimension or precision,
nor an additive cost per frame layer or source call. Measurements,
postselection, and supplied nonstabilizer catalysts have different
contracts, recorded in the [literature comparison](../RELATED_WORK.md#29-constant-depth-obstructions-with-initialized-work).

The [bounded checks](../../tests/test_clean_aware_depth_obstruction.py)
audit the initialized-projector normalization, witness arithmetic,
complete-frame embedding, binary carry identities, the weighted phase
gap, the literal zero-clean attainer, terminal clean-sector rigidity,
and the native three-layer counterexample. The uniform
three-layer conclusion is analytic and uses the cited polynomial-limit
theorem; finite circuits do not certify compactness or its error constant.
