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
chapter.

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
