# Transport through a literal permutation source

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Six-layer comparison](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[Diagonal-source transport](NATIVE_CLIFFORD_TRANSPORT.md)

This independently audited resumption result supplies an alternative scoped
transport theorem. It requires a literal permutation source of bounded forward
degree and a monomial native inner family. No actual three-band source outside
the retained diagonal-source population has yet been supplied. It therefore
does not establish a new family of actual six-layer queries, an unrestricted
six-layer theorem, a native depth lower bound, or a resource improvement.

## Hypotheses and native input

Let $`p`$ be a fixed permutation of the work bit strings, with every forward
coordinate having Boolean degree at most $`h\ge1`$. Let $`F`$ be its literal
permutation matrix, up to one irrelevant fixed global phase. No assumption is
made on the degree of $`p^{-1}`$ or on the Clifford-hierarchy membership of
$`F`$.

Let $`L_x`$ be the actual inner Clifford family from the six-layer reduction,
or any Clifford family with the same native estimates: for $`L_x`$ and its
actual adjoint, every fixed Hermitian Pauli probe has output label degree at
most two and literal sign degree at most six. Assume additionally that every
$`L_x`$ is monomial in the computational basis.

Let $`K_x`$ be the actual specified Clifford-valued comparison, including its
actual scalar, and assume

```math
\|F L_xF^\dagger-K_x\|_{\mathrm{op}}\le2\eta
\quad\text{for every }x.
```

Projective bounds suffice in the rounding lemma below. These are statements
about operators on the entire arbitrary work Hilbert space, without an
initialized-subspace promise.

**Permutation-source transport theorem.** If $`\eta\lt\sin(\pi/8)`$,
then the endpoints agree projectively and their signed tableaus obey

```math
\deg_x S_K\le6h+8,\qquad
\deg_x\kappa_K(q)\le24h+16
\quad\text{for every fixed Hermitian Pauli probe }q.
```

Consequently the retained actual outer extractor would give XOR degree at
most $`6h+8`$ for visible outer shear, or $`24h+16`$ in general. For
$`h=2`$ these are 20 and 64. These consequences require all the hypotheses
above on the specified physical circuit; the source and its inverse remain
charged. The degree constants are conservative; no sharpness is claimed.

## 1. A rounding gap independent of third-level membership

Write

```math
\delta=\sqrt{2-\sqrt2}=2\sin(\pi/8),\qquad
 d_{\mathrm{pr}}(U,K)=\inf_{|\alpha|=1}\|U-\alpha K\|_{\mathrm{op}}.
```

For any monomial Clifford $`L`$ and literal permutation $`F`$, the operator
$`U=FLF^\dagger`$ is monomial, and ratios of any two of its nonzero
computational entries are fourth roots of unity. Suppose
$`d_{\mathrm{pr}}(U,K)\lt\delta`$ for a Clifford $`K`$.

Every Clifford computational column is a stabilizer state and has equal
nonzero magnitudes $`2^{-r/2}`$ for some nonnegative integer $`r`$.
If $`K`$ is not monomial, at least one column has $`r\ge1`$. The overlap
of this normalized column with the one basis vector supporting the
corresponding column of $`U`$ has magnitude at most $`1/\sqrt2`$. Thus
their vector distance, after any common phase choice, is at least $`\delta`$.
Hence $`K`$ must be monomial.

If the monomial permutations differ, some two corresponding columns have
orthogonal supports, so operator distance is at least $`\sqrt2`$,
independently of global phase. Thus $`K`$ and $`U`$ have the same
permutation.

The diagonal ratios $`UK^\dagger`$ have, up to one common scalar, values
among the fourth roots of unity. If two such values differ, their angular
separation is at least $`\pi/2`$; the least possible maximum distance from
a common unit scalar is at least $`2\sin(\pi/8)=\delta`$. Therefore
$`d_{\mathrm{pr}}(U,K)\lt\delta`$ implies that $`U`$ is exactly a scalar
multiple of $`K`$. This is a full-operator gap and does not make any
intermediate band Clifford.

Apply this with $`d_{\mathrm{pr}}(U,K)\le2\eta\lt\delta`$. No
third-level condition was used.

## 2. Computational permutation coefficients

Since $`L_x`$ is monomial Clifford, write

```math
L_x|z\rangle=\alpha_x\omega^{\theta_x(z)}|A_xz+b_x\rangle,
\qquad\omega=e^{i\pi/4},
```

where $`\alpha_x`$ is arbitrary, $`\theta_x(0)=0`$, and $`\theta_x`$
is even modulo eight with weighted Clifford-polynomial dependence on
$`z`$. The matrix $`A_x`$ is invertible. Its entries and those of
$`A_x^{-1}`$ have degree at most two by the native label bound. The bit
$`(b_x)_j`$ is the literal sign in $`L_x^\dagger Z_jL_x`$, so
$`\deg_x b_x\le6`$. The displacement of $`L_x^{-1}`$ also has degree
at most six by applying the same statement to $`L_x^\dagger`$.

By the rounding lemma, the monomial permutation of $`K_x`$ is

```math
k_x(y)=p\bigl(A_xp^{-1}(y)+b_x\bigr).
```

For any fixed $`y`$, the string $`p^{-1}(y)`$ is fixed, so
$`\deg_x k_x(y)\le6h`$. Because $`K_x`$ is Clifford, $`k_x`$ is affine
in $`y`$. Its translation and each matrix column are extracted by
evaluation at $`y=0`$ and $`y=e_j`$, and therefore have degree at most
$`6h`$ in $`x`$. Apply the same argument to the actual inverse
$`L_x^\dagger`$ to see that the affine inverse matrix and inverse
displacement $`d_x=k_x^{-1}(0)`$ also have degree at most $`6h`$.

The degree of $`p^{-1}`$ is never used: it is evaluated only at fixed work
strings in this argument.

## 3. Relative computational phase at each fixed input

Fix a work string $`z`$. Write the native fixed-probe image as

```math
L_xX^zL_x^\dagger
 =(-1)^{\sigma_x(z)}\tau\bigl(a_x(z),c_x(z)\bigr),
\qquad
\tau(a,c)=i^{a\cdot c}X^aZ^c,
```

with an integer dot product in the Hermitian section. Both labels have degree
at most two, and $`\sigma_x(z)`$ has degree at most six. Applying this
equality to $`|b_x\rangle`$ gives the literal modulo-eight identity

```math
\theta_x(z)
 =4\sigma_x(z)+2\bigl[a_x(z)\cdot c_x(z)\bigr]
                +4\bigl[c_x(z)\cdot b_x\bigr]\pmod8.
```

The first nontrivial bit, at weight two, has degree at most four. The
weight-four bit has degree at most eight: weight-two sums of degree-four
bits have a second symmetric carry of degree at most eight; the last term
has degree at most $`2+6=8`$; and $`\sigma_x(z)`$ has degree at most six.
Both bounds are width independent.

The relative computational phase of $`K_x`$ at a fixed input $`y`$ is

```math
\phi_x(y)=\theta_x\bigl(p^{-1}(y)\bigr)
             -\theta_x\bigl(p^{-1}(0)\bigr)\pmod8.
```

Its weight-two bit has degree at most four and its weight-four bit has
degree at most eight. Subtraction preserves these bounds, including its
carry.

Since $`K_x`$ is monomial Clifford and $`\phi_x(0)=0`$, write

```math
\phi_x(y)
 =2\sum_i\ell_i(x)y_i
  +4\sum_{i\lt j}m_{ij}(x)y_iy_j\pmod8.
```

The low bit of $`\ell_i`$ has degree at most four, and its high bit has
degree at most eight, by $`\phi_x(e_i)/2`$. The polar coefficient
$`m_{ij}`$ is extracted from the second phase difference at
$`0,e_i,e_j,e_i+e_j`$, divided by four; its Boolean degree is at most
eight, with all carries retained. Consequently the symmetric phase shear
$`B_x`$, with diagonal $`\ell_i\bmod2`$ and off-diagonal entries
$`m_{ij}`$, has degree at most eight.

## 4. Symplectic and literal sign bounds

Let $`A'_x`$ be the affine linear matrix of $`K_x`$ and $`d_x`$ its
inverse displacement. Conjugation of a fixed $`\tau(a,b)`$ gives output
label

```math
a'_x=A'_xa,\qquad
c'_x=(A'_x)^{-\mathsf T}(b+B_xa).
```

Thus $`\deg_x a'_x\le6h`$ and $`\deg_x c'_x\le6h+8`$, proving the
symplectic bound.

Evaluate this conjugation at $`|0\rangle`$ to obtain its literal sign
$`\kappa_x(a,b)`$:

```math
\begin{aligned}
4\kappa_x(a,b)
 ={}&\phi_x(d_x+a)-\phi_x(d_x)
     +2(a\cdot b)+4(b\cdot d_x)\\
    &-2(a'_x\cdot c'_x)\pmod8.
\end{aligned}
```

The left side is indeed divisible by four because the transported operator
is a Hermitian Pauli; no unsigned replacement is made.

For fixed $`a`$, the phase difference has the weighted expansion

```math
2\sum_i\ell_i a_i(1-2d_i)
 +4\sum_{i\lt j}m_{ij}(a_id_j+a_jd_i+a_ia_j)\pmod8,
```

where the dependence of $`\ell_i,m_{ij},d_i`$ on $`x`$ is suppressed.
Its weight-four bit has degree at most $`6h+8`$, and its weight-two bit has
degree at most four. The term $`4(b\cdot d_x)`$ has degree at most
$`6h`$. Each Boolean product in $`a'_x\cdot c'_x`$ has degree at most
$`12h+8`$. Multiplication of their sum by two gives a weight-four bit of
degree at most $`2(12h+8)=24h+16`$, including the second symmetric carry.
Combining these terms does not exceed this bound: possible carries between
weight-two bits have degree at most $`2(12h+8)`$. Dividing the known
multiple of four yields $`\deg_x\kappa_x(a,b)\le24h+16`$.

This proves the claimed bounds using the retained native input ledger.

## 5. Scope and possible extension

This argument allows $`F`$ outside the third level, but requires a literal
permutation source and monomial inner family. The rounded endpoint $`K_x`$
is necessarily monomial. The monomial condition on $`L_x`$ remains an
additional hypothesis; full workspace return alone does not supply it.
A fixed-Clifford input change is permissible if the conjugated native
$`L_x`$ is monomial. A fixed-Clifford output change leaves the label bound
unchanged and contributes a sign of degree at most $`3(6h+8)`$. The resulting
sign-degree bound is therefore

```math
\max(24h+16,18h+24),
```

still 64 for $`h=2`$.

All coordinate permutations are analytical descriptions of the given source.
They do not create a free oracle. In the six-layer application $`F`$ must
be the specified charged three-band prefix. In particular, the quadratic
HRT map with high-degree inverse is not known here to admit a three-band
native source; invoking it abstractly would not enlarge the settled
six-layer population. A quadratic permutation outside the third level,
realized in three physical bands, would make this source class concretely
distinct from the retained diagonal third-level class, provided a compatible
monomial $`L_x`$ and actual returned query exist.

## 6. A four-band candidate does not supply the missing three-band source

Consider the chronological pair of Toffolis with controls $`(a,b)`$, target
$`c`$, then controls $`(a,c)`$, target $`b`$. Its literal permutation is

```math
p(a,b,c)=(a,\ b+ac+ab,\ c+ab).
```

It is quadratic. In the $`a=1`$ sector it sends $`(b,c)`$ to $`(c,b+c)`$,
of order three; in the $`a=0`$ sector it is identity. Its inverse is also
quadratic, but the source is not third-level: conjugating $`X_a`$ by the
permutation source gives

```math
(a,b,c)\longmapsto(a+1,\ c+ab,\ b+c+ac),
```

which is nonaffine and therefore not Clifford. This is an explicit source
outside the diagonal third-level class, but it has not been placed in three
native bands here.

The retained [two-band dirty Toffoli](../../docs/PARALLEL_DIRTY_LOOKUP.md#12-the-two-bit-indicator-has-optimal-t-depth-two)
implements each gate on $`a,b,c`$ and one arbitrary helper $`d`$. Reuse the
same returned helper in the second gate. Concatenating the retained Gaussian-
elimination parity-basis words gives the fully charged conservative ledger

```math
T=16,\qquad G=92,\qquad D_T=4.
```

Here $`G`$ counts 88 CNOTs and four Hadamards; there are zero clean qubits,
three logical work wires and one additional arbitrary helper. Its actual
reversed inverse has the same ledger. Both words act on all sixteen columns
with literal permutation phases and return $`d`$, including its external
reference correlations. No source call or actual inverse is free.

The middle two bands in this particular schedule do not directly merge.
In common physical coordinates, the second band of the first Toffoli has
axis $`Z_dZ_aZ_b`$, while the first band of the second has axis $`Z_dX_b`$.
These axes anticommute. A physical T band in any fixed Clifford coordinates
has commuting axes, so these two displayed sets cannot simply be combined
as one band. This failure does not rule out a different three-band
resynthesis of the complete permutation. Nor does the four-band source
constitute the specified three-band prefix of the six-layer comparison.

The [focused checks](../../tests/test_permutation_source_transport.py)
verify the sharp rounding-gap examples, phase and Hermitian-section identities,
the candidate's complete native operator and actual inverse, its non-C3
witness, its emitted ledger, and the incompatible middle axes. The analytic
argument establishes the arbitrary-width theorem; these finite checks do not
populate the missing actual six-layer case.
