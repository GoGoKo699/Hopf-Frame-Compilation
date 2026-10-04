# Ideal-angle stability of complete Hopf frames

[Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Capped precision](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision) · [Current depth frontier](OPEN_PROBLEM.md)

Exact angle perturbations admit a dimension-independent root-sum-square
bound for the complete real Hopf frame. This is the ideal-angle estimate
used by the unary groups in Result D. It does not replace a full
initialized-isometry certificate for approximate source circuits; those
[leakage diagnostics](../research/depth/SOURCE_LEAKAGE_DIAGNOSTICS.md)
are preserved separately.

## 1. A sharp angle-error bound on the full frame

Use the repository convention $`R_y(\theta)=e^{-i\theta Y}`$ and
$`W(\theta)=L_{n-1}(\theta)\cdots L_0(\theta)`$. For two real
angle tuples, choose their actual real differences and put

```math
e_d=\max_p|\Delta\theta_{d,p}|,\qquad
\Lambda_0=0,\qquad
\Lambda_{d+1}=\frac{\Lambda_d+
 \sqrt{\Lambda_d^2+4e_d^2}}2.
```

Then the complete operator, including every marker column, satisfies

```math
\|W(\theta+\Delta\theta)-W(\theta)\|
\le\min\{2,\Lambda_n\}
\le\min\!\left\{2,
 \sqrt{e_0^2+2\sum_{d=1}^{n-1}e_d^2}\right\}.
```

In particular $`\sqrt{2\sum_d e_d^2}`$ is a uniform bound. There
is no factor one-half in e: the Hopf rotation convention differs from
the common half-angle gate convention. Tensoring with untouched work or
a reference preserves the same operator norm.

### Nested support proves the recurrence

Interpolate $`\theta(t)=\theta+t\Delta\theta`$, for
$`0\le t\le1`$, and set $`V_d=L_{d-1}\cdots L_0`$,
with $`V_0=I`$. Define the old and new logical subspaces

```math
\mathcal S_d=\mathrm{span}\{|p\rangle|0^{n-d}\rangle\},\qquad
\mathcal M_d=\mathrm{span}\{|p\rangle|1\rangle|0^{n-d-1}\rangle\}.
```

Both have dimension $`2^d`$, and
$`\mathcal S_{d+1}=\mathcal S_d\oplus\mathcal M_d`$.
Earlier layers act inside $`\mathcal S_d`$ and are identity outside
it. Their logarithmic derivative $`K_d=V_d^\dagger\dot V_d`$
is therefore supported there. Since rotations in one layer have disjoint
supports, its new generator is the diagonal angle-difference map between
these old and new spaces. In the displayed decomposition,

```math
K_{d+1}=
\begin{pmatrix}K_d&-B_d^\dagger\\B_d&0\end{pmatrix},\qquad
B_d=D_d(V_d|_{\mathcal S_d}),\qquad
D_d=\mathrm{diag}(\Delta\theta_{d,p}).
```

Thus $`\|B_d\|=e_d`$. Applying the Rayleigh bound to the
Hermitian matrix $`iK_{d+1}`$ reduces its norm to the largest
eigenvalue of

```math
\begin{pmatrix}\|K_d\|&e_d\\e_d&0\end{pmatrix}.
```

This proves $`\|K_{d+1}\|\le\Lambda_{d+1}`$.
The recurrence also gives
$`\Lambda_{d+1}^2\le\Lambda_d^2+2e_d^2`$, and
$`\Lambda_1=e_0`$. These statements hold at every t. Integrating
$`\|\dot W(t)\|=\|K_n(t)\|`$ proves the finite bound.

### The square-root factor cannot be removed

If all angle differences at depth d have common magnitude e_d, then
$`B_d^\dagger B_d=e_d^2I`$. A unitary change of basis on the new
subspace turns $`iK_{d+1}`$ into

```math
\begin{pmatrix}iK_d&e_d I\\e_d I&0\end{pmatrix}.
```

Each old eigenvalue lambda gives the two new eigenvalues
$`(\lambda\pm\sqrt{\lambda^2+4e_d^2})/2`$.
The real skew-symmetric generators have symmetric Hermitian spectra, so
their norm attains the recurrence exactly. This holds at every base
angle tuple, with either sign allowed at each node.

For common magnitudes one, write the exact norm as $`\lambda_j`$.
Then

```math
\lambda_j^2-\lambda_{j-1}^2=2-\lambda_j^{-2},\qquad
2n-H_n\le\lambda_n^2\le2n-1,
\qquad H_n=\sum_{j=1}^n j^{-1}.
```

The lower inequality uses $`\lambda_j^2\ge j`$. Hence
$`\lambda_n/\sqrt n\to\sqrt2`$. This is also a sharp finite-error
limit: at the regular base tuple with every angle $`\pi/4`$, choose
every change $`\delta_n=n^{-2}`$. The product's second derivative
has norm at most $`n^2`$, so the Taylor remainder is at most
$`n^2\delta_n^2/2`$. Dividing by $`\delta_n\sqrt n`$ makes it
vanish. The finite error divided by angle RSS therefore tends to
$`\sqrt2`$.

### An exact finite relative spectrum

When the differences at depth d are $`\pm e_d`$, the finite relative
spectrum also depends only on these magnitudes, not on the base angles
or the choices of signs. This is stronger than the differential statement.
Let $`U_\theta,U_\phi`$ be the two old prefix unitaries on
$`\mathcal S_d`$, and put
$`E_d=U_\theta^\dagger U_\phi`$ and
$`F=U_\theta^\dagger D`$, where D is the diagonal sign matrix.
With $`c=\cos e_d`$ and $`s=\sin e_d`$, the next relative unitary
on $`\mathcal S_d\oplus\mathcal M_d`$ is

```math
E_{d+1}=
\begin{pmatrix}cE_d&-sF\\sF^\dagger E_d&cI\end{pmatrix}.
```

The matrix F is unitary. Conjugating by $`\mathrm{diag}(I,F^\dagger)`$
gives

```math
\begin{pmatrix}cE_d&-sI\\sE_d&cI\end{pmatrix}.
```

Diagonalize the old unitary E_d. For each old eigenvalue z, the two new
eigenvalues are exactly the roots of

```math
\lambda^2-\cos(e_d)(z+1)\lambda+z=0.
```

Starting with the single eigenvalue $`E_0=1`$ recursively determines
the full relative spectrum, with multiplicities. The final active space
is the whole logical space. Its distance from identity, and therefore
$`\|W(\phi)-W(\theta)\|`$, depends only on the depth magnitudes.
This polynomial statement needs no ordering or unwrapping of eigenphases.

Compare with zero base angles and all-positive differences e_d. The
prepared first column then has overlap $`\prod_d\cos e_d`$ with
$`|0^n\rangle`$. Isospectrality consequently gives, at every base
tuple and every sign pattern,

```math
\|W(\phi)-W(\theta)\|
\ge\sqrt{2-2\prod_d\cos e_d}.
```

### A necessary precision bound for nearest angular grids

Consider the restricted architecture that rounds every angle to a nearest
point of an independently specified uniform angular grid of spacing
$`2\delta`$. Choose an interior midpoint target at each coordinate,
with both adjacent grid points permitted by its angle domain. Every
rounding error is then $`\pm\delta`$. The signs may be chosen jointly;
the finite spectral result makes the lower bound independent of them.
For $`0\le\delta\le\pi/2`$, an operator-error guarantee
$`0\lt\eta\lt1`$ therefore requires

```math
\cos^n\delta\ge1-\eta^2/2,\qquad
\delta^2\le\frac{-2\ln(1-\eta^2/2)}n,
```

where the second implication uses
$`\cos\delta\le\exp(-\delta^2/2)`$. Thus
$`\delta=O(\eta/\sqrt n)`$ at small fixed or varying accuracy.
If the grid spacing is $`\gamma2^{-m}`$ for a fixed angular range
constant gamma, the necessary precision is

```math
m\ge\log_2(1/\eta)+\tfrac12\log_2 n-O(1).
```

This is an actual worst-case lower bound for this nearest angular-grid
architecture, not merely a sufficient error certificate. It is not a
lower bound for the dyadic cosine/sine source, another quantizer, the
capped compiler, or T-depth. No gate count is inferred from the number
of grid bits.

The later exploratory material is preserved in
[Source leakage and the limits of angular error estimates](../research/depth/SOURCE_LEAKAGE_DIAGNOSTICS.md).
