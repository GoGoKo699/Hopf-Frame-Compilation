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
fits do not imply it. The structural chapter proves the restrictions on
[three-mask](STRUCTURAL_COMPILATION_LIMITS.md#three-masks-with-two-flat-cliffords-a-width-independent-gap)
and [four-mask flat-Clifford words](STRUCTURAL_COMPILATION_LIMITS.md#four-masks-with-three-flat-cliffords-a-width-independent-gap),
[bounded balanced-Haar words](STRUCTURAL_COMPILATION_LIMITS.md#global-obstruction-for-unrestricted-haar-orientations),
and [transport into one terminal Clifford gauge](STRUCTURAL_COMPILATION_LIMITS.md#transport-into-one-terminal-clifford-gauge-is-rigid).
These restrictions leave the general coverage hypothesis, including
larger constant factor counts and target-dependent Clifford skeletons,
undecided.

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

### An explicit nonflat factorization with n diagonal slots

Let $`B_d`$ act as SH on logical qubit d and identity elsewhere. Define
$`D_d`$ to have entry $`e^{-i(-1)^z\theta_{d,p}}`$ on prefix p,
target bit z, and zero lower suffix; on every other suffix its entry is
one. Since $`B_dZ_dB_d^\dagger=Y_d`$ and $`B_d`$ leaves the prefix
and suffix unchanged, the exact depth layer is

```math
L_d=B_dD_dB_d^\dagger,\qquad
W=B_{n-1}D_{n-1}B_{n-1}^\dagger\cdots B_0D_0B_0^\dagger.
```

This retains literal phases and applies to zero, negative, and unequal
angles. Each Clifford basis change is explicit and nonflat on the whole
logical register. With

```math
h=\lceil\log_2 n\rceil,\qquad
\ell=N+h+2,\qquad q=\left\lceil\frac{N+h+6}{2}\right\rceil,
```

the n packed diagonal calls have total initialized-isometry error at most
$`n2^{-\ell}\le\eta/4`$ and the explicit ledger

```math
a=2,\qquad b=q+n+1\le N+n+7,\qquad
T\le n(1344N+36q-18)=O(nN).
```

The same flags and dirty wires are reused, with Section 2's telescoping
bound accounting for leakage and return. The Clifford count is
$`O(nNq)`$, including the explicit basis changes. For a final-depth-only
frame there is one slot: take $`\ell=N+2`$ and
$`q=\lceil(N+6)/2\rceil`$ to obtain error at most $`\eta/4`$ and
T-count $`1344N+36q-18=O(N)`$ within the same workspace allowance.
This includes the witness for the
[four-flat-mask obstruction](STRUCTURAL_COMPILATION_LIMITS.md#four-masks-with-three-flat-cliffords-a-width-independent-gap).

The full-frame word has $`K(n)=n`$, so it does not meet the bounded-slot
hypothesis. Although the $`B_d`$ commute with one another, $`B_d`$
generally fails to commute with $`D_e`$ when $`e\ne d`$: for
$`e\lt d`$ the diagonal uses qubit d in its lower-suffix projector,
and for $`e\gt d`$ it can depend on prefix bit d through the unequal
angle table. Collecting these basis changes at the endpoints therefore
does not preserve the complete operator.

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

<a id="positive-conditional-root-stability"></a>

## 6. Stability from every conditional root state

Let $`W_\mu,W_\nu`$ be two prescribed positive real frames on
the same ordered binary tree, including endpoint angles in
$`[0,\pi/2]`$. For each interval I, let $`p_I,q_I`$ be their
normalized local root vectors. Their definition uses only subtree
angles. Put

```math
c=\min_I\langle p_I,q_I\rangle,\qquad
\delta=\sqrt{1-c^2}
       =\max_I\|p_Ip_I^{\mathsf T}-q_Iq_I^{\mathsf T}\|.
```

Then the literal frame distance satisfies

```math
\delta\le\|W_\mu-W_\nu\|\le\min\{2,4\delta+\delta^2\}.
```

For the lower bound, the two conjugations of the common subtree
marker projection differ by the corresponding root-projector
difference. The norm of a difference of two orthogonal
projections is the larger of the two cross-range norms.
Each cross-range norm is at most $`\|W_\mu-W_\nu\|`$ by
inserting the other frame.

For the upper bound first assume positive branch masses, with
global leaf probability measures $`\mu,\nu`$. Write
$`h_i^\mu,h_i^\nu`$ for the corresponding marker vectors
at internal interval i, and set $`a_i=\langle h_i^\nu,p_i\rangle`$.
The restriction of $`p_I`$ to $`i\subseteq I`$ is
$`\sqrt{\mu(i)/\mu(I)}p_i`$. Subtree Parseval, including its
root vector, gives

```math
\sum_{i\subseteq I}\mu(i)|a_i|^2
  =\mu(I)\bigl(1-|\langle p_I,q_I\rangle|^2\bigr)
  \le\delta^2\mu(I).
```

In $`U=W_\nu^{\mathsf T}W_\mu`$, disjoint marker intervals
have zero overlap. Treat the global root column as an ancestor
of every marker. For the part with row i a strict descendant of
column j, the ancestor vector restricts on i to a scalar multiple
of $`p_i`$. All other columns have zero inner product with $`p_i`$.
Thus this part acts by
$`(U_{\rm desc}z)_i=a_i\langle p_i,W_\mu z\rangle`$.

For a physical vector x, put $`f(j)=x_j/\sqrt{\mu(j)}`$.
Then $`\langle p_i,x\rangle=\sqrt{\mu(i)}\,\mathbb E_{\mu|i}f`$.
At a level t, cover the intervals with conditional average
larger than t in modulus by their maximal disjoint members.
Their total weight $`\mu(i)|a_i|^2`$ is at most $`\delta^2`$
times the measure of their union. Integrating the levels gives

```math
\sum_i \mu(i)|a_i|^2|\mathbb E_{\mu|i}f|^2
 \le \delta^2\|Mf\|_{L^2(\mu)}^2
 \le 4\delta^2\|f\|_{L^2(\mu)}^2 ,
```

where $`Mf`$ is the dyadic maximal conditional average of
$`|f|`$. For the last inequality, maximal disjoint intervals
give $`t\mu\{Mf\gt t\}\le\int_{\{Mf\gt t\}}|f|\,d\mu`$.
Integration and Cauchy–Schwarz yield $`\|Mf\|_2\le2\|f\|_2`$.
The two strict triangular parts of U each have norm at most
$`2\delta`$, by exchanging the frames for the other part.

For a node i with children L and R put
$`r_J=\langle p_J,q_J\rangle`$ and
$`b_J=\sqrt{\mu(J)\nu(J)/(\mu(i)\nu(i))}`$. The prescribed
positive rotations give
$`r_i=b_Lr_L+b_Rr_R`$ and
$`\langle h_i^\mu,h_i^\nu\rangle=b_Rr_L+b_Lr_R`$.
Since $`b_L+b_R\ge r_i\ge c`$ and both child overlaps are
at least c, the marker overlap is at least $`c^2`$. The root
overlap is at least c. Hence $`\|\mathrm{diag}(U)-I\|\le\delta^2`$.
Endpoint angles follow by continuity of every subtree vector in
its own angles. This includes subtrees of zero incoming mass.

### A weighted angle certificate

Let $`\theta_v,\phi_v`$ denote the two positive angle tuples
and $`\mu(v|I),\nu(v|I)`$ the conditional branch masses
computed from I. The root affinities obey

```math
1-\langle p_I,q_I\rangle
 =\sum_{v\subseteq I}\sqrt{\mu(v|I)\nu(v|I)}
     \bigl(1-\cos(\theta_v-\phi_v)\bigr),
```

where v ranges over internal nodes including I. At one split the
two child overlap coefficients sum to the cosine of the local angle
difference. Repeated substitution with terminal affinities one proves
the identity. In particular, a uniform bound

```math
\sup_I\sum_{v\subseteq I}
   \frac{\mu(v|I)+\nu(v|I)}{2}|\theta_v-\phi_v|^2\le A
```

gives complete error at most $`4\sqrt A+A`$. The supremum includes
rare subtrees regardless of their incoming mass. Any paid reference
frame within this certificate can be used at its own native costs,
with its implementation error added in operator norm.

For example, fix $`L\ge6`$ and use the positive balanced
reference. Allow deviations of magnitude at most $`2^{-L}/64`$
on the leftmost internal spines of mutually disjoint subtrees,
with arbitrary spine lengths. Each path enters at most one such
subtree. Once on a spine, both frames continue along it with
probability at most 3/4 per split, so the expected number of
perturbed nodes from any starting interval is at most four.
The exact balanced reference then has complete error less than
$`2^{-(L+1)}`$, even if the angle deviations accumulate to a
large sum on some individual path. Compiling balance at precision
L+1 by the
[prefix-shared construction](../depth/PREFIX_ROTATION_PROFILES.md)
adds at most $`2^{-(L+1)}`$ and supplies its one-clean,
square-root-dirty, $`O_L(\log\log(n+4))`$ depth bound for
this promised family. The comparison includes every column
and tensors with the identity on work and references.

### Arbitrary orthogonal candidates can incur height loss

The two-frame bound uses the exact conditional marker structure.
For the general projector test of Section 4, take the prescribed
comparison frame to be the identity. On $`m=n-1\ge1`$ address
bits let $`S_j`$ exchange the normalized uniform vectors of the
two children of every dyadic interval at depth j, with zero action
on the orthogonal complement. Put $`K=\sum_{j=0}^{m-1}S_j`$.
Then $`K|+^m\rangle=m|+^m\rangle`$. Each crossing block of K
for an address interval has norm at most $`\sqrt2`$: only
ancestor levels contribute and the outputs occupy disjoint sibling
intervals, with squared coefficients bounded by $`1,1/2,1/4,\ldots`$.
The same estimate holds for each single column of K.

For $`0\lt\xi\le\pi/2`$ use the real orthogonal candidate

```math
F=\exp[(\xi/m)K\otimes XZ].
```

A dyadic projection of size at least two is an address-interval
projection tensored with identity on the last bit. Its commutator
with the generator is bounded by the crossing estimate. Removing
the anchor coordinate costs at most the single-column bound.
Integrating the commutator gives
$`\max_I\|FE_IF^{\mathsf T}-E_I\|\le2\sqrt2\,\xi/m`$.
For singleton intervals $`E_I=0`$. Yet
$`F|+^m\rangle|0\rangle`$ has last-bit-one component of norm
$`\sin\xi`$. A computational diagonal preserves that last-bit
sector on the same input, so the distance from every such
diagonal is at least $`\sin\xi`$.
