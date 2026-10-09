# Bounded diagonal factorization as a linear-T compiler

[Packed diagonal compiler](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames) · [Complete-frame contract](../../docs/HOPF_INTERFACE.md) · [Structural restrictions](STRUCTURAL_COMPILATION_LIMITS.md)

A uniformly bounded number of computational-basis diagonals, interleaved
with explicit Clifford circuits, would suffice for the complete-frame
linear-T endpoint. The packed diagonal compiler makes every fixed factor
count affordable at the prescribed width. Exact factor coverage alone
would suffice: a terminating certified search supplies approximate factors
without a computable or continuous phase-selection map.

## 1. Conditional compilation theorem

Let $`n\ge3`$, $`N=2^n`$, $`\eta=2^{-N}`$, and let
$`W(\theta)`$ be the prescribed complete real Hopf frame for any real
angle tuple. Each angle admits terminating certified sine and cosine
evaluation. No uniform running-time bound on those evaluators is assumed.

**Coverage hypothesis.** There is an absolute integer $`K\ge1`$ and,
for each n, an effectively specified finite menu of Clifford skeletons
with at most K diagonal slots such that every $`W(\theta)`$ has an
exact literal factorization in that menu. A skeleton has the form

```math
F=C_KD_KC_{K-1}\cdots C_1D_1C_0,
```

where every $`D_j`$ is computational-basis diagonal and unitary, and
each $`C_j`$ is supplied as an explicit Clifford gate word on the n
logical qubits. Identity slots can pad shorter skeletons. The factors
may have non-Clifford baselines that cancel at the identity; no factor
is required to approach identity as $`\theta\to0`$.

**Theorem.** Under this hypothesis, a deterministic finite compiler
produces a circuit over $`\{H,S,\mathrm{CNOT},T,T^\dagger\}`$ with

```math
a=2,\qquad b=N+n+7,\qquad
\|VJ_2-J_2(W(\theta)\otimes I_b)\|\le\eta,\qquad T\le C_KN,
```

where $`C_K`$ is independent of n and the angle tuple. The physical
width is $`N+2n+9`$. The norm includes all logical and arbitrary dirty
inputs, their references, clean leakage and approximate dirty return.
All inverses are actual circuit inverses, and no intermediate
initialization is used.

The hypothesis is a sufficient algebraic route, not an established
factorization of the Hopf family. Local tangent coverage or finite-size
fits do not imply it. Restrictions on particular three-slot or anchored
factorizations have the narrower scope recorded in the linked structural
chapter. Its
[global Haar obstruction](STRUCTURAL_COMPILATION_LIMITS.md#global-obstruction-for-unrestricted-haar-orientations)
excludes every uniformly bounded balanced-Haar word in either orientation,
without any regularity assumption on its masks. Dense Clifford mixers
remain compatible with this conditional theorem.

## 2. Native precision and the finite-dimensional fallback

The [packed diagonal theorem](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames)
compiles one literal diagonal to error $`2^{-\ell}`$, for
$`\ell\ge6`$, with

```math
q=\left\lceil\frac{\ell+4}{2}\right\rceil,\qquad
b_{\rm diag}=q+n+1,\qquad
T\le1344N+36q-18,\qquad G=O(Nq).
```

Here q counts dirty core qubits. Two clean flags and the same
$`n+1`$ dirty selectors suffice; the source and masks require no other
work. Set

```math
h=\lceil\log_2K\rceil,\qquad
\ell=N+h+2,\qquad
q=\left\lceil\frac{N+h+6}{2}\right\rceil.
```

If $`h\le N+6`$, then $`q+n+1\le N+n+7`$ and

```math
K2^{-\ell}\le\eta/4.
```

Thus the native diagonal cost is $`O_K(N)`$ T gates and
$`O_K(N^2)`$ Clifford gates. Add the finite, explicitly emitted
Clifford interlayers to the latter count. Their cost is not assumed
free; the endpoint imposes no bound on the total finite Clifford count.
Every $`K\le16384`$ satisfies this width inequality for all
$`n\ge3`$, since $`N+6\ge14`$.

Only finitely many allowed n can violate $`h\le N+6`$. For those,
use the [layerwise complete-frame compiler](../../docs/OPERATOR_SOURCE_COMPILER.md#6-precision-workspace-and-full-frame-composition)
at $`L=N`$. It already has two clean qubits,
$`b=N+n+7`$, full-isometry error at most $`\eta`$, and
$`T=O(N+nN)`$. In this finite range, n is bounded by a constant
depending only on K, so this is also $`C_KN`$. This fallback needs
neither extra work nor an unrestricted unitary-synthesis oracle.

The diagonal stages reuse one initialized embedding $`J=J_2`$.
For actual unitary stages $`A_j`$ and ideal stages $`U_j`$ (including
the exact Clifford interlayers),

```math
A_s\cdots A_1J-JU_s\cdots U_1
=\sum_{j=1}^s A_s\cdots A_{j+1}
 (A_jJ-JU_j)U_{j-1}\cdots U_1.
```

Taking norms bounds the total by the sum of stage errors. Actual
earlier leakage propagates under the unitary suffix; flags and dirty
work are not reset between stages. Tensoring with an arbitrary
reference leaves the bound unchanged.

## 3. Certified search replaces factor extraction

For $`r=1,2,\ldots`$, define the finite rational unit-circle grid

```math
G_r=\{-1\}\cup\left\{
\frac{1-t^2+2it}{1+t^2}:
t=\frac{a}{2^r},\quad a\in\mathbb Z,\quad |a|\le2^{2r}
\right\}.
```

These grids are nested, their union is dense in the unit circle, and
every entry has exact rational real and imaginary parts and modulus one.
At stage r, enumerate every skeleton in the finite menu and every
assignment of $`G_r`$ to its diagonal entries. All candidate products
F retain the literal scalar phases of their native Clifford words.

The supplied angle evaluators and finite interval matrix multiplication
give certified enclosures for every W entry, including at zero or
singular angles. Obtain complex entry disks of radius at most
$`\eta/(32N)`$ for W, and evaluate each candidate F to the same
entry radius. The two Frobenius uncertainties sum to at most
$`\eta/16`$. Accept a candidate only when a rational upper bound on
its squared Frobenius residual is strictly below $`\eta^2/4`$:

```math
\|F-W\|_F\lt\eta/2.
```

The upper bound can be computed with finite rational arithmetic from
the disks; rational square-root upper bounds may be refined by a fixed
additional margin, for example $`\eta/32`$ in norm. A candidate
failing this finite enclosure test is skipped, including a possible
boundary case. No real equality test is required. Clifford entries are
computable from their finite native words, so every candidate test
terminates.

To prove that some candidate is accepted, choose any exact factorization
guaranteed by coverage. At a sufficiently fine grid, approximate each
of its diagonal entries within $`\eta/(8K\sqrt N)`$. Telescoping
products of unitary factors gives

```math
\|F-W\|_F
\le\sqrt N\,\|F-W\|
\le\sqrt N\sum_{j=1}^K\|D_j-\widehat D_j\|
\lt\eta/8.
```

Such a candidate passes the finite enclosure test with strict room to
spare and is eventually enumerated. The only continuity used is that of
a finite matrix product in its entries. The argument requires no
continuous factor map, no nonsingular inverse chart, and no effective
selection among exact factorizations.

Compile the accepted rational diagonals directly from their real and
imaginary entries using Section 2's precision. Since operator norm is
bounded by Frobenius norm, the total error is

```math
\|VJ_2-J_2(W\otimes I_b)\|
\lt\eta/2+\eta/4=3\eta/4.
```

The search and input evaluation may be extremely expensive, but they
are finite classical preprocessing. The output is a finite explicit
native gate list with the stated width and T count.

The same argument allows target-dependent Clifford interlayers when
coverage is quantified over the finite native n-qubit Clifford group:
enumerate that finite group with explicit words and all literal scalar
representatives. Quotienting by global phase would not meet the literal
error contract. An unbounded factor count $`K(n)`$ still incurs
$`O(K(n)N)`$ through this construction; extra precision does not
amortize repeated source calls.

## 4. Complete coverage through nested projectors

For every dyadic interval $`I=[a,a+s)`$, write $`\Pi_I`$ for its
coordinate projection. Let $`\psi_I`$ be the normalized root column of
the local subtree frame, using only the angles inside I, and define

```math
E_I=\Pi_I-|a\rangle\langle a|,\qquad
P_I=\Pi_I-|\psi_I\rangle\langle\psi_I|.
```

These definitions use no division by the global root-state mass of I.
They therefore apply to every admitted real angle tuple, including a
subtree with zero global mass.

**Exact criterion.** The prescribed complete frame satisfies

```math
W E_I W^\dagger=P_I\qquad\text{for every dyadic }I.
```

Conversely, a unitary F satisfying all these identities is $`F=WD`$
for a computational unitary diagonal D. Thus one actual right phase
mask corrects every column phase.

The inputs in $`I\setminus\{a\}`$ are precisely the subtree's internal
markers. Ancestor rotations fix those inputs; their final columns are
the orthonormal local details spanning $`\psi_I^\perp`$ inside I.
This proves the identities. If I has children L,R and right-child left
endpoint b, then

```math
E_I-E_L-E_R=|b\rangle\langle b|,\qquad
P_I-P_L-P_R=W|b\rangle\langle b|W^\dagger.
```

Also $`I-E_{[0,N)}=|0\rangle\langle0|`$. These differences recover
each individual column projector. Equal unit-vector projectors differ
by a unit scalar, proving the converse without changing the prescribed
completion.

**Operator-norm stability.** Suppose a logical unitary F obeys

```math
\|F E_I F^\dagger-P_I\|\le\delta\quad\text{for every }I,
\qquad \kappa=(2n-1)\delta\lt1.
```

There is a computational unitary diagonal D with

```math
\|FD-W\|\le\kappa+1-\sqrt{1-\kappa^2}
\le\kappa+\kappa^2.
```

To prove it, set $`U=W^\dagger F`$. Unitary invariance gives
$`\|[U,E_I]\|\le\delta`$. Split the whole space into
$`|0\rangle`$ and $`E_{[0,N)}`$. Its off-diagonal part has norm at
most delta, since the two off-diagonal blocks have orthogonal domains
and ranges. Recursively split each internal marker space by

```math
E_I=E_L+E_R+|b\rangle\langle b|.
```

Every off-diagonal block in this three-part split has norm at most
delta, by the commutator with $`E_L`$ or $`E_R`$. Comparing component
norms with the scalar matrix
$`\delta(\boldsymbol1\boldsymbol1^{\mathsf T}-I_3)`$ bounds this part
by $`2\delta`$. At one tree depth the parent spaces are mutually
orthogonal, so their direct sum has the same bound. There are $`n-1`$
nontrivial depths, down to intervals of size four; size-two marker
spaces are already one-dimensional. The root split and these depth
parts partition every off-diagonal entry once. Consequently

```math
\|U-\mathrm{diag}(U)\|\le\delta+2(n-1)\delta=\kappa.
```

Each unit column then gives $`|U_{jj}|\ge\sqrt{1-\kappa^2}`$.
Choose $`D_{jj}=\overline{U_{jj}}/|U_{jj}|`$. The triangle inequality
proves the stated bound. If F admits certified matrix-entry evaluation,
as does a concrete candidate word, these overlaps are bounded away from
zero and their phases admit terminating certified evaluation without an
equality oracle. D is an actual factor, not a phase quotient of the error
contract.

For example, $`\delta=\eta/[8(2n-1)]`$ gives
$`\|FD-W\|\le9\eta/64\lt\eta/4`$. Any implementation errors of F
and D then compose through the full-output identity in Section 2. The
criterion is a complete logical coverage test; separately preparing the
local states $`\psi_I`$ does not implement its simultaneous identities.

### Hierarchical rank-one cuts

The same marker support also gives, at every dyadic interval I,

```math
\mathrm{rank}((I-\Pi_I)W\Pi_I)\le1,\qquad
\mathrm{rank}(\Pi_IW(I-\Pi_I))\le1,\qquad
\mathrm{rank}([\Pi_I,W])\le2.
```

Only the anchor input a can propagate outside I; all other inputs there
are internal markers. Conversely, ancestors inject outside inputs into
I only through a. Descendant rotations send all such contributions to
scalar multiples of the same local root state $`\psi_I`$. This proves
both rank-one bounds, including zero branches. The two off-diagonal
commutator blocks have orthogonal domains and ranges. These identities
describe the coupled unitary hierarchy; they do not factor it into a
bounded number of global native gates.

## 5. Direct weighted-Haar linear normal form

There is a different all-column identity on the chart
$`0\lt\theta_{d,p}\lt\pi/2`$.
Set $`\psi=We_0`$ and, for each node v with interval
$`I_v=L_v\cup R_v`$, define
$`\mu_v=\sum_{x\in I_v}\psi_x^2`$ and child masses
$`\mu_L,\mu_R`$. Its marker column is

```math
W_{x,\lambda(v)}=
\begin{cases}
-\psi_x\sqrt{\mu_R/(\mu_L\mu_v)},&x\in L_v,\\
\phantom{-}\psi_x\sqrt{\mu_L/(\mu_R\mu_v)},&x\in R_v,\\
0,&x\notin I_v.
\end{cases}
```

To verify the formula, the node rotation contributes its negative sine
or cosine to the marker column, and the descendant rotations supply
the normalized restriction of psi to that child. Earlier layers act
trivially on the marker input. This proves every prescribed column.

Let $`u=N^{-1/2}(1,\ldots,1)^{\mathsf T}`$ and let $`S_n`$
have root column u and marker column
$`|I_v|^{-1/2}1_{I_v}`$. Define diagonal tables by

```math
\begin{aligned}
a_0&=\sqrt N,& b_0&=0,\\
a_v&=\frac{\sqrt{|I_v|}\sqrt{\mu_v}}{2\sqrt{\mu_L\mu_R}},&
b_v&=\frac{\sqrt{|I_v|}(\mu_L-\mu_R)}
 {2\sqrt{\mu_L\mu_R\mu_v}}.
\end{aligned}
```

Combining the constant and signed columns gives the literal identity

```math
W=\mathrm{diag}(\psi)
 \bigl[Q_n\mathrm{diag}(a)+S_n\mathrm{diag}(b)\bigr].
```

This is a linear normal form, not a unitary factorization: $`S_n`$
is not a free mixer. Separately block encoding the three factors of the
first term has no uniform product-normalization bound, even at fixed N.
Choose a positive psi with almost
all norm at leaf zero and every other amplitude epsilon. A deepest
pair disjoint from zero has $`a_v=1/\varepsilon`$, while
$`\|\mathrm{diag}(\psi)\|`$ approaches one. The product of these
separate factor normalizations is therefore at least order
$`1/\varepsilon`$.
The separate-factor normalization therefore does not supply the
constant-cost unitary coverage required in Section 1.
