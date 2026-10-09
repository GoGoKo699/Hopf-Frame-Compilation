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
[global Haar obstruction](STRUCTURAL_COMPILATION_LIMITS.md#global-obstruction-for-alternating-balanced-haar-masks)
also excludes every uniformly bounded strictly alternating balanced-Haar
architecture, without any regularity assumption on its masks.

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

## 4. Nonalternating Haar candidate

This subsection preserves an unresolved coverage candidate outside the
[strict-alternation obstruction](STRUCTURAL_COMPILATION_LIMITS.md#global-obstruction-for-alternating-balanced-haar-masks).
Let $`Q_n`$ be the fixed balanced Hopf frame and $`V_n=Q_n^2`$.
The proposed five-mask equation is

```math
W=D_0V_nD_1V_n^\dagger D_2V_nD_3V_n^\dagger D_4.
```

It contains eight Haar calls, including consecutive equal orientations.
Coverage by this word is unproved. The exact algebraic problem is

```math
V_n^\dagger D_0^\dagger W D_4^\dagger V_n
=D_1(V_n^\dagger D_2V_n)D_3.
```

Since $`V_n`$ is real, the matrix inside parentheses is complex
symmetric. The free exterior masks must first make the left side have
symmetric entrywise absolute values, and then permit left/right phase
removal into the fixed algebra $`V_n^\dagger\mathcal ZV_n`$, where
$`\mathcal Z`$ is the complex diagonal algebra.
The first condition alone does not imply the second. A useful next proof
would establish both conditions under tree growth with these same five
slots; parameter counting, local coverage, and finite fits do not do so.

**An exact failed simplification.** The exterior choice
$`D_0=D_4=I`$ cannot be imposed as a gauge. At $`n=3`$, let W
have root angle $`\pi/4`$ and all other angles zero, so
$`W=R_{0,4}(\pi/4)`$. With $`V=Q_3^2`$, exact arithmetic gives

```math
Z=V^\dagger WV,\qquad
Z_{0,1}=-\frac{\sqrt2}{16},\qquad Z_{1,0}=-\frac18,\qquad
|Z_{0,1}|^2-|Z_{1,0}|^2=-\frac1{128}.
```

This violates the necessary modulus symmetry for every choice of the
three internal masks. The stdlib exact certificate is in
[`test_global_haar.py`](../../tests/test_global_haar.py). It refutes
the anchored simplification, not the free five-mask equation.

**Conditional native accounting.** If the free five-mask equation were
proved for all frames, certified search from Section 3 would apply to
its computable fixed mixers. At $`\ell=N+5`$, each packed diagonal
has error $`\eta/32`$ and needs

```math
b_D=n+1+\left\lceil\frac{N+9}{2}\right\rceil\le N+n+7.
```

The five calls cost $`O(N)`$. The eight Haar calls and their actual
inverses cost $`O(n^3)=O(N)`$ by the
[fixed-mixer construction](STRUCTURAL_COMPILATION_LIMITS.md#balanced-haar-boundary-identities);
they are not free Clifford interlayers. The same two clean flags and
dirty bank suffice. A searched product within $`\eta/2`$, followed by
the five compiled masks, has initialized-isometry error at most
$`\eta/2+5\eta/32=21\eta/32`$. The full-output telescoping
identity in Section 2 supplies work return and reference correctness
without resetting intermediate leakage. The missing premise is global
coverage, not factor extraction or the resource ledger.

### Direct weighted-Haar linear normal form

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
A useful construction would need a coupled unitary dilation that
preserves the cancellation between local row and column weights, with
priced work return; neither that dilation nor an all-chart extension
is provided by this identity.
