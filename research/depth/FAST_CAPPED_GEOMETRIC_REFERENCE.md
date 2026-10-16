# Fast capped geometric references by coherent truncation

[Exact capped source](../../docs/FAULT_TOLERANT_COMPILER.md#3-an-exact-geometric-source-of-logarithmic-width) ·
[Complete compact-coin flags](CAUSAL_COMPACT_COIN_DEPTH.md#7-whole-coarse-flags-reference-work-and-actual-inverses) ·
[Small-workspace early stage](RELAXED_COUNT_LOCAL_DEPTH.md#2-an-early-stage-with-a-logarithmic-precision-register)

The reference part of
the small-workspace early stage admits depth linear in its requested
precision, even when the original geometric coin is much longer. Keep
the exact original coin and causal matching in the mathematical embedding;
approximate only the native unitary preparing its reference column.
This chapter does not by itself prove a global complete-frame bound.

## 1. Two exact capped columns and their distance

For a power of two $`M\ge2`$, write

```math
g_M=\sum_{j=0}^{M-2}2^{-(j+1)/2}|j\rangle
       +2^{-(M-1)/2}|M-1\rangle.
```

The retained recursive Gray-path circuit prepares this column exactly
with $`O(M)`$ T and elementary Clifford gates and
$`O(\log_2M)`$ initialized width, including its exactly returned
private scratch. Its actual reverse is the inverse native word.

Let $`2\le R\lt M`$ also be a power of two. Regard $`g_R`$ as
a column on the original $`\log_2 M`$ coin bits, with its high bits
zero. The first $`R-1`$ amplitudes agree. At the final occupied
position the truncated cap has amplitude $`2^{-(R-1)/2}`$, while
the original column has amplitude $`2^{-R/2}`$. Hence

```math
\langle g_M,g_R\rangle
=1-2^{-(R-1)}+\sqrt2\,2^{-R}
=1-(2-\sqrt2)2^{-R},                                  \qquad(1)
```

```math
\|g_M-g_R\|_2^2=(4-2\sqrt2)2^{-R}\lt2^{1-R}.          \qquad(2)
```

For a tolerance $`\epsilon=2^{-s}`$, with $`s\ge1`$, choose

```math
R_*=2^{\lceil\log_2(2s+2)\rceil},\qquad
R=\min\{M,R_*\}.                                     \qquad(3)
```

If $`R=M`$ the preparation is exact. Otherwise (2) gives an error
strictly below $`\epsilon`$. Physically apply the exact $`P_R`$
word to the low $`\log_2R`$ coin bits, leave the higher coin bits
untouched, and reserve its own logarithmic initialized private scratch.
This is a full native unitary on every occupied input, although its
specified column is the padded $`g_R`$.

## 2. A full-space comparison for transported zero tests

Let C be the full original coin register and A the private scratch of
the truncated preparation. Include any uniform-column Hadamards in its
actual unitary U. The initialized column u is the uniform column tensored
with padded $`g_R`$ and zero A. The target column v uses the SAME
uniform column, the original $`g_M`$, and zero A. The source has
literal positive amplitudes, so $`c=\langle u,v\rangle`$ is real
and positive. Equations (1)–(2) give $`\|u-v\|\le\epsilon`$.

There is an explicit analysis-only unitary S mapping u to v while
obeying $`\|S-I\|=\|u-v\|`$. If u equals v, take S to be
identity. Otherwise put

```math
w=\frac{v-cu}{\sqrt{1-c^2}},\qquad
Su=cu+\sqrt{1-c^2}\,w,\qquad
Sw=-\sqrt{1-c^2}\,u+cw,
```

and let S be identity on the orthogonal complement of this plane.
Its two nontrivial eigenvalues are the conjugate rotations with cosine
c, giving the asserted norm. Both u and v have A zero, so this plane
is within the same initialized-scratch subspace. No gate for S is
implemented.

Set $`\widehat U=SU`$. Then, on the ENTIRE physical space,

```math
\|\widehat U-U\|\le\epsilon,\qquad
\widehat U|0^{CA}\rangle=v,
```

```math
\widehat U|0^{CA}\rangle\langle0^{CA}|\widehat U^\dagger
=|v\rangle\langle v|.                                \qquad(4)
```

The actual word needs no dirty register. Tensoring U, S and
$`\widehat U`$ with identity on arbitrary dirty registers and
external references preserves all these equalities and bounds. The
comparison is therefore stronger than a single successful-root overlap;
it supplies the full-space unitary estimate needed to transport a source
anchor. It does not identify U's other columns with a chosen ideal
completion.

## 3. Keep the original causal matching and all flag columns

Use the padded exact reference geometry from
[the local compiler](RELAXED_COUNT_LOCAL_DEPTH.md#3-exact-reference-coins-without-a-precision-bank).
For each bin i, its original exponent is
$`s_i'=2^{\lceil\log_2(q_i\tau_0+1)\rceil}-1`$ and
$`M_i=s_i'+1`$. The mathematical reference remains

```math
\gamma_i=|+\rangle^{\otimes q_i}\otimes g_{M_i}.
```

Every alias acceptance probability, positive atom, weight class, hash
catalogue and causal matching is built for this ORIGINAL reference.
Nothing is rounded or truncated in these data. In particular, every
original geometric weight class retains its original positive mass and
at least $`2^{q_i}`$ atoms. Its inverse-polynomial matching density
and rank-window bound are unchanged. The smaller physical source is
used only when reading the already matched coarse flag.

For p bins, let $`F=\bigotimes_i U_i`$ be the actual parallel
reference word and $`\widehat F=\bigotimes_i\widehat U_i`$
its comparison word from (4). Unitary telescoping gives

```math
\|F-\widehat F\|\le\sum_{i=1}^p\epsilon_i.            \qquad(5)
```

The computational nonzero test for each bin includes its entire original
coin C_i and EVERY private source bit A_i. All other fixed encoding
padding is zero on the ideal matched image and remains untouched.
Equation (4) then gives the retained complete-band identity: the root
maps to all zero; band i maps into the sector whose last nonzero
augmented block is i. Orthogonality of a marker to $`\gamma_i`$
is sufficient for this statement, independently of the other columns of
$`\widehat U_i`$. It includes all marker columns and their coherent
superpositions, not only the prepared root.

Let $`L_H`$ be the retained H-guarded label XOR, with its original
mapping of the last nonzero block into the old coarse label. In ordinary
operator order, the actual sandwich is $`F L_H F^\dagger`$.
Its full-operator distance from the comparison sandwich obeys

```math
\|F L_H F^\dagger-
  \widehat F L_H\widehat F^\dagger\|
\le2\sum_i\epsilon_i.                                \qquad(6)
```

If the preceding actual compact embedding has initialized-isometry
error $`\delta_E`$, the complete embedding, reference inverse,
label XOR, reference forward and actual embedding inverse has error
at most

```math
2\delta_E+2\sum_i\epsilon_i.                          \qquad(7)
```

This is a norm bound on every initialized logical column and arbitrary
tag, dirty and reference input. Earlier embedding leakage propagates
through the actual unitaries; it is never reset or treated as fresh zero
work.

The reference preparations run unconditionally on disjoint blocks and
private scratch. They do not share an enabled-H source control. On H
zero the central XOR is exactly identity, so
$`F I F^\dagger=I`$ on all occupied coin and private work. Thus
the complete sandwich has literal inactive identity without assuming
that an individual reference preparation is inactive there. Its actual
inverse is always used.

## 4. Native resources and the scope of the improvement

For stage precision l, take
$`s=l+\lceil\log_2(p+1)\rceil+c`$ and use (3) for each bin,
where the absolute c is fixed by the surrounding error budget.
Then $`\sum_i\epsilon_i\le2^{-(l+c)}`$ and
$`R_*=O(l+\log_2(p+1))`$. The complete reference word has

```math
T_{\rm ref}=O\!\left(p[l+\log_2(p+1)]\right),
```

```math
G_{\rm ref}=O\!\left(\sum_i q_i+
                  p[l+\log_2(p+1)]\right),\qquad
D_{T,\rm ref}=O(l+\log_2(p+1)).                        \qquad(8)
```

The extra initialized private scratch is bounded by
$`\sum_i O(\log_2M_i)`$. At the retained coarse geometric
precision $`\tau_0=O(\log_2(m+2))`$, one has
$`p=O(\log_2^2(m+2))`$ and $`M_i=O(m)`$. Thus this scratch
is $`O(\log_2^3(m+2))`$, with no additional precision-sized
initialized bank. The original coin-coordinate extension already has
the same bound; it is not reduced by dropping the original high coin
bits. All these bits remain in the augmented nonzero predicates.

There is no reference truth table, coherent evaluator, approximate
single-qubit synthesis, supplied state or measurement in this construction.
The only extra classical calculation is the explicit integer cutoff (3).
All private reference pools are disjoint, and the native inverses have
the same cost. Therefore (8) replaces the earlier reference-depth
$`O(m)`$ allowance by $`O(l+\log_2\log_2(m+4))`$ while
preserving the same initialized-work asymptotics.

Using this component in a globally selected early stage still requires
the charged address, query, router, reflection and simultaneous dirty
allocation ledger for every OTHER component. In particular a generic
standalone bound $`2^{C(m+\chi)}`$ cannot be absorbed at m=n
without being sharpened. The reference construction removes that
possible exponential standalone cost only for the reference sources
themselves. It does not assert the complete-frame depth bound from this
component alone.
