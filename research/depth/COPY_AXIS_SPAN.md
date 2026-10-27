# A non-Clifford copy needs extensive output-axis span

[Encoded dirty copy](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Field-shear source](FIELD_SHEAR_INTERTWINER.md)

A non-Clifford copy can evade the fixed-Clifford copy obstruction only
by supplying enough independent output directions through its rotation
axes. At fixed error below square root two, a bounded-degree Clifford
query requires an extensive span, even for an arbitrary seed and
arbitrary dirty helpers. The result bounds copy T-count, not T-depth.
It is compatible with the linear count budget and does not improve R96.

## 1. Full-space copy and absorbed axes

Let x have p bits, let $`m=2^p`$, and let R be any unitary on a
carrier and its arbitrary helpers. Let $`Q_x`$ be Clifford unitaries
whose binary symplectic matrices have Boolean degree at most d, and set
$`R_x=Q_xRQ_x^\dagger`$. Scalar phases of Q do not matter. The
physical output Y has m arbitrary qubits. A fixed copy C, independent
of x, is compared on the complete space by

```math
\eta=\max_x\|C(R_x\otimes I_Y)C^\dagger
                    -R_x\otimes X_{Y_x}\|.
```

All private wires mixed by C must be included, with identity in the
ideal comparison where appropriate. No initialized-sector promise is used.
Absorb fixed Clifford interlayers in a chosen native word to write

```math
C=U F,\qquad U=\prod_{j=1}^t\exp(i\theta_j P_j),
```

where F is Clifford and each $`P_j`$ is a signed Hermitian Pauli.
Literal scalar phases can be retained without changing conjugation.
The theorem allows arbitrary angles. For a Clifford+T word, absorb
all Clifford gates into F and conjugate each T/T-dagger axis accordingly;
there are t axes when the copy has t T-type gates.

Let $`L_Y\subseteq\mathbb F_2^{2m}`$ be the span of the restrictions
of these axes to the physical output. Let $`\mathsf X_Y`$ denote the
m-dimensional subspace of X-only output labels and put

```math
r_X=\dim(L_Y\cap\mathsf X_Y),\qquad
M_d(p)=\sum_{j=0}^{\min(d,p)}\binom pj.
```

The intersection is important: a large span of output-Z axes alone
does not erase any independent one-hot output-X direction.

## 2. Robust span bound

Every such representation obeys

```math
\eta\ge\sqrt{2-2\sqrt{\min\{1,(M_d(p)+r_X)/m\}}}.     \qquad(1)
```

In particular, for $`\eta\lt\sqrt2`$,

```math
r_X\ge m(1-\eta^2/2)^2-M_d(p),\qquad
 t\ge r_X.                                           \qquad(2)
```

For exact copying, $`r_X\ge m-M_d(p)`$. Integer ranks may be
rounded up. These conclusions hold for every absorbed-axis description,
so choosing another Clifford factorization cannot evade the lower bound.

### Conjugation preserves quotient-sector weight

Partition the full Pauli operator basis according to its output label
modulo $`L_Y`$. These are orthogonal subspaces in normalized
Hilbert--Schmidt inner product. Conjugation by a Pauli rotation maps a
Pauli either to itself or to a linear combination of its label and that
label plus the rotation axis. It therefore preserves each output coset.
Its inverse does too. Conjugation by U is unitary on each of these
operator subspaces and preserves its total squared coefficient norm.

Individual coefficients can interfere and cancel within a coset. The
proof uses only its total norm, which cannot flow to a different coset.
It does not assume that a non-Clifford copy permutes Pauli coefficients.

### Count matches after taking the quotient

Expand R in the Hermitian Pauli basis as $`R=\sum_v r_v\tau(v)`$.
Unitarity gives $`\sum_v|r_v|^2=1`$; R need not be Hermitian.
After conjugation by Q and F, the output label of the v term is a
vector-valued degree-d function $`F_v(x)`$. The values of this function
lie in the span $`A_v`$ of its at most $`M_d(p)`$ ANF coefficient
vectors, including the constant coefficient.

Write $`e_x`$ for the physical one-hot X label. The set

```math
S_v=\{x:F_v(x)+L_Y=e_x+L_Y\}
```

has at most $`M_d(p)+r_X`$ elements. To see this, the labels
$`\{e_x:x\in S_v\}`$ span a space E of dimension $`|S_v|`$.
Its image in the quotient by $`L_Y`$ lies in the image of $`A_v`$,
which has dimension at most $`M_d(p)`$. The kernel on E is
$`E\cap L_Y\subseteq\mathsf X_Y\cap L_Y`$. Rank-nullity gives
$`|S_v|\le M_d(p)+r_X`$, without treating quotient labels as
independent or distinct.

Before U acts, the squared norm in the target output coset is

```math
w_x=\sum_v|r_v|^2[\,x\in S_v\,].
```

Distinct seed labels remain distinct full labels under Q and F, even
when their output cosets coincide. Thus this expression has no cross
terms. The preceding sector argument shows that the same w is present
after U. Its address average is at most
$`\min\{1,(M_d(p)+r_X)/m\}`$.

The ideal $`R_x\otimes X_{Y_x}`$ lies entirely in that coset and
has normalized Hilbert--Schmidt norm one. Cauchy--Schwarz bounds its
real overlap with the actual conjugated operator by $`\sqrt{w_x}`$.
The full operator error implies

```math
1-\eta^2/2\le\frac1m\sum_x\sqrt{w_x}
 \le\sqrt{\frac1m\sum_x w_x}
 \le\sqrt{\min\{1,(M_d(p)+r_X)/m\}},
```

which proves (1). Squaring only when its left side is positive proves
(2). Operator norms already include arbitrary external references.

## 3. Consequence and limits

For the one-shear field orbit, $`d=1`$ and $`M_d(p)=p+1`$.
Exact copying requires at least $`2^p-p-1`$ independent X directions
in the projected rotation-axis span, and at least that many T gates.
At every fixed error below square root two the requirement is
$`\Omega(2^p)`$. A localized non-Clifford correction of a Clifford
copy, with $`r_X=o(2^p)`$, has asymptotic error at least square root two.

Arbitrarily broad seeds and arbitrary Clifford mixing are already
allowed. Arbitrary numbers of identity helpers do not change the
argument. A native band can contain an extensive independent axis
family, so this result gives no growing depth lower bound. The
required count is compatible with the discovery budget. A constructive
next candidate must distribute its non-Clifford copy action across the
output directions while still paying the seed and complete echo.

## 4. Exact evidence

[The fixtures](../../tests/test_copy_axis_span.py) exhaust the two-bit
X-label affine maps against all zero- or one-generator X subspaces,
verify the sharper intersection dimension with mixed X/Z axes, and
check exact full-matrix quotient-sector norm preservation under
noncommuting native Pauli rotations, including interference within a
coset. They test the ingredients; the all-size bounds follow from the
proof above, not from finite enumeration.
