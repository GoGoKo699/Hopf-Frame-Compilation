# Angle errors and coherent source leakage in complete Hopf frames

[Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Capped precision](AMORTIZED_DIRTY_LOOKUP.md#capping-the-source-precision) · [Current depth frontier](OPEN_PROBLEM.md)

Exact angle perturbations admit a dimension-independent root-sum-square
bound for the complete real Hopf frame. The current amplified source
layers do not admit that replacement for their full initialized-isometry
errors: a family of their actual shared-flag circuits accumulates leakage
linearly. These are different error models. Neither result changes the
current count/depth frontier or proves that fixed-accuracy frame depth
must contain a logarithm.

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

## 2. Two exact checks distinguish the error notions

A bound by only the maximum local error is impossible even for the first
column. Start from zero angles and perturb only prefix p equal to zero
at every depth by delta. The resulting state has overlap
$`(\cos\delta)^n`$ with $`|0^n\rangle`$, so its error is exactly

```math
\sqrt{2(1-\cos^n\delta)}.
```

Every layer error is $`2\sin(|\delta|/2)`$. Taking
$`\delta=1/n`$ makes their ratio grow as $`\sqrt n`$.

For a small complete-frame example, let n equal two and change all three
angles from zero to $`\pi/3`$. In basis $`00,01,10,11`$,

```math
W=\begin{pmatrix}
1/4&-\sqrt3/2&-\sqrt3/4&0\\
\sqrt3/4&1/2&-3/4&0\\
\sqrt3/4&0&1/4&-\sqrt3/2\\
3/4&0&\sqrt3/4&1/2
\end{pmatrix}.
```

Both layer operator errors equal one. The squared-error matrix has
eigenvalues $`(5\pm\sqrt{13})/4`$, each twice, giving
$`\|W-I\|=\sqrt{5+\sqrt{13}}/2\gt\sqrt2`$.
Its first-column error is only $`\sqrt{3/2}`$. Thus even unit-coefficient
RSS of layer operator errors is not valid for the full frame. Section 1
uses angle differences, not these chordal errors, and its proved constant.

## 3. Actual source layers can accumulate leakage linearly

We now use the literal two-flag Q and amplified layer
$`\mathcal A=-QRQ^\dagger RQ`$ from
[the operator-source proof, Sections 4–5](OPERATOR_SOURCE_COMPILER.md#4-two-flags-encode-a-suffix-controlled-rotation).
Here J initializes both flags, $`P=JJ^\dagger`$, and $`R=I-2P`$.
All mask queries and inverses are the actual circuits; selector and
suffix work returns exactly. The original source core and both flags
are reused between layers. In initialized-isometry formulas, ideal
logical operators are tensored with identity on all dirty work.

**Counterfamily.** Give every layer a common source width m and set

```math
\Delta=2^{2-m}\le(256n)^{-2},\qquad
\sigma_d=(-1)^{n-1-d}.
```

At every active row of depth d choose the same angle

```math
\theta_d=
\begin{cases}
\arcsin\sqrt{\Delta/2},&\sigma_d=-1,\\
\arcsin\sqrt{3\Delta/2},&\sigma_d=+1.
\end{cases}
```

These are canonical positive real Hopf angles. At each depth choose any
permitted rational estimates for its cosine and sine, and reuse those
estimates across its equal active rows in the existing certified encoding.
If $`\mathcal V=\mathcal A_{n-1}\cdots\mathcal A_0`$, then

```math
\epsilon_d:=\|\mathcal A_dJ-JL_d\|\le4\Delta,\qquad
\|\mathcal VJ-JW\|\ge\frac{n\Delta}{8}.
```

The norm includes flag leakage and arbitrary dirty inputs. Consequently
neither a dimension-independent maximum-error bound nor a
dimension-independent RSS bound in these actual local errors is possible:
the global error divided by $`\sqrt{\sum_d\epsilon_d^2}`$ is at
least $`\sqrt n/32`$.

### Rounding gives alternating radial defects

The coefficient grid spacing is Delta and the certified estimate error
is at most $`\Delta/8`$, so final coefficient error is at most
$`5\Delta/8`$. Write the encoded coefficients as c_d,s_d and set
$`\rho_d=1-c_d^2-s_d^2`$. For $`\Delta\le1/16`$, the cosine
drop in the first case lies between $`\Delta/4`$ and
$`4\Delta/15`$; in the second it lies between $`3\Delta/4`$ and
$`4\Delta/5`$. Every permitted estimate therefore rounds to

```math
c_d=1\quad(\sigma_d=-1),\qquad
c_d=1-\Delta\quad(\sigma_d=+1).
```

Together with $`|s_d-\sin\theta_d|\le5\Delta/8`$, this gives

```math
|\rho_d-\sigma_d\Delta/2|\le2\Delta^{3/2}.
```

Equation (24) of the source chapter gives the full local bound
$`\epsilon_d\le(5\sqrt2/2)\Delta\le4\Delta`$.

### The entire rejected space has the same zero-reference sign

The scalar source block can be written
$`\mathcal S_f=cI+X_aD_f`$, where

```math
D_f=(MN_f-N_fM)/2,\qquad D_f^\dagger=-D_f,\qquad
D_f^2=-(1-c^2)I.
```

For comparison only, replace the cosine scalar by I and the sine scalar
by $`X_aD_s/\sqrt{1-s_d^2}`$. The cosine change has norm
$`\sqrt{2(1-c_d)}\le\sqrt{2\Delta}`$; the sine change is at
most $`\sqrt2|s_d|\le2\sqrt\Delta`$. Hence the resulting
unitary $`Q_d^0`$ obeys $`\|Q_d-Q_d^0\|\le2\sqrt\Delta`$.
This comparison does not require the masks to converge as m changes.

Explicitly,

```math
Q_d^0=H_b\,\mathrm{diag}_b(I,X_aU)\,H_b,\qquad
U=(XZ_t)D_s/\sqrt{1-s_d^2},\qquad U^\dagger=U,\quad U^2=I.
```

Conjugation by $`C_a(U)`$ commutes with R and reduces this to the
two-flag matrix $`Q_{\rm base}=H_b\mathrm{CX}_{b\to a}H_b`$.
Direct multiplication gives
$`-Q_{\rm base}RQ_{\rm base}RQ_{\rm base}=-R`$.
Thus its amplified comparison is $`\mathcal A_d^0=-R`$ on the
whole logical/core space, not just the initialized subspace. Inactive
suffixes already have this exact zero-coefficient form. The three Q
occurrences in amplification now give the uniform bounds

```math
\|\mathcal A_d+R\|\le6\sqrt\Delta,\qquad
\|L_d-I\|\le2\sqrt\Delta.
```

The normalized comparison is an analytic unitary, not an extra circuit
or an uncharged operation in the compiler.

### A single rejected-flag matrix element witnesses accumulation

For the accepted block $`B_d=J^\dagger Q_dJ`$, exact amplification
algebra gives, on an active row,

```math
\mathcal A_dJ=\rho_d Q_dJ+2JB_d.
```

The component with flags $`(b=1,a=0)`$ is therefore
$`\rho_d(c_dI-s_dXZ_t)/2`$. Taking the logical-zero matrix element
kills $`XZ_t`$ and leaves $`\rho_dc_d/2`$ times identity on every
dirty input. The ideal frame has zero amplitude in these rejected flags.

Telescope the actual product against the ideal one. In each summand
replace the future physical stages by $`(-R)^{n-1-d}`$ and the ideal
past by I. Their respective norm errors are at most
$`6(n-1-d)\sqrt\Delta`$ and $`2d\sqrt\Delta`$. Multiplication by
the local error bound $`4\Delta`$ and summation bound the remainder by
$`16n(n-1)\Delta^{3/2}`$. Thus, for any normalized dirty state psi,
the matrix element from $`|00\rangle|0^n\rangle|\psi\rangle`$
to $`|b=1,a=0\rangle|0^n\rangle|\psi\rangle`$ is

```math
\sum_{d=0}^{n-1}\sigma_d\rho_dc_d/2+E,\qquad
|E|\le16n(n-1)\Delta^{3/2}.
```

Each displayed scalar differs from $`\Delta/4`$ by at most
$`(17/16)\Delta^{3/2}`$. Consequently

```math
\left|\text{amplitude}-\frac{n\Delta}{4}\right|
\le17n^2\Delta^{3/2}.
```

Since $`n\sqrt\Delta\le1/256`$, the amplitude is at least
$`n\Delta/8`$ in absolute value. This proves the counterfamily using
the complete physical rejected space, with no reset or postselection.

## 4. What this closes, and what it leaves open

The angle theorem applies to exact perturbed Hopf rotations. It does not
bound amplified-source leakage by the accepted block's angle error.
For the literal common-precision source family above, the usual
initialized-isometry certificate cannot be uniformly replaced by RSS.
Changing only its error estimate is insufficient; a circuit modification
must first establish a stronger full-input error structure.

The counterfamily is a composition of the existing common-m layer
gadgets. It is not a lower bound for the exact capped allocation, another
compiler, joint source synthesis, or fixed-accuracy T-depth. Delta shrinks
with n in the witness. It does not prove that every fixed-accuracy
implementation requires $`n\log n`$ depth or total assigned precision.

Even an angle-only RSS certificate with
$`e_d\le C2^{-m_d}`$ retains a logarithmic precision budget. Convexity
under $`2\sum_d C^2 2^{-2m_d}\le2^{-2L}`$ gives

```math
\frac1n\sum_d m_d\ge L+\frac12\log_2 n+\log_2 C+\frac12.
```

This concerns that sufficient certificate, not the true optimized error.
The finite midpoint argument in Section 1 separately makes the same
half-logarithmic order necessary for nearest uniform angular grids only.
Neither statement transfers that necessity to the current dyadic
cosine/sine implementation or to unrestricted T-depth.

The [flag-echo audit](HOPF_FLAG_ECHO.md) now supplies exact errors for the
two-half-angle square and the three other diagonal Pauli flag echoes.
They retain first-order radial error at generic angles; the two-flag-Z
echo has an exact same-mask exception. A separate
[phase-calibrated fixed-point filter](HOPF_RADIAL_FILTER.md) does suppress
the full radial error quadratically. Its charged phase words use the
same two flags, allowing the modified physical frame to combine angular
square-sum stability with an additive quadratic remainder. Its smaller
source-precision cap leaves the asymptotic count/depth frontier unchanged.

The [bounded checks](../tests/test_hopf_error_accumulation.py) cover small
ideal frames and a compressed representation of the exact source algebra.
They support matrix identities and the finite counterfamily. They are
not an emitted general source compiler, an asymptotic proof, or evidence
for a new depth frontier. The uniform bounds are the arguments above.
