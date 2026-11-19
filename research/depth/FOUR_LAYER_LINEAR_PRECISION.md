# A linear rank–precision obstruction for four native bands

For every fixed clean budget $`a`$, approximate the complete-frame witness
$`I\otimes R_y(\tfrac12\arccos(1/3))`$ by a four-native-T-layer circuit with full
initialized-isometry error at most $`2^{-L}`$. Write $`r`$ for the binary symplectic pairing
rank of the transported Pauli axes in its two middle T bands. The comparison circuit may use any
dirty width. For sufficiently large $`L`$,
```math
 r=\Omega_a(L).
```
The conclusion holds on each physical decomposition, hence on the minimum attainable rank. The
original R74 estimate gives $`\Omega_a(\sqrt L)`$ for this witness. Only its estimate of the
support of a nonzero cluster is changed below. The statement gives a growing resource
requirement within four layers. It supplies no growing lower bound on general frame depth at
fixed precision. R75 already excludes four zero-clean layers for sufficiently accurate
approximation. The rank improvement is particularly relevant with positive, fixed clean budgets.
[R74 §5.1](FOUR_LAYER_CROSSING_RANK.md#51-quantitative-rigidity-on-the-initialized-sector)

## Divisibility of phase functions on a binary cube

Put $`\omega=e^{i\pi/4}`$, $`\mathcal O=\mathbb Z[\omega]`$, and $`\pi=1-\omega`$. Then
$`\mathcal O/(\pi)=\mathbb F_2`$ and the field norm of $`\pi`$ is $`2`$ by evaluating the
cyclotomic polynomial $`X^4+1`$ at one. The $`\pi`$-valuation of $`2`$ is four, since the four
conjugates $`1-\omega^t`$ for odd $`t`$ are associates (their ratios are integral and have norm
one). A weighted cubic phase has exponent consisting modulo eight of a constant, arbitrary
linear coefficients, even quadratic coefficients, and cubic coefficients divisible by four. This
property survives an affine binary substitution: lift a parity to its integer multilinear
polynomial, whose degree-$`j`$ coefficients are divisible by $`2^{j-1}`$.

Expand $`u\mapsto\omega^{q(u)}`$ on a Boolean cube as the unique multilinear polynomial over
$`\mathcal O`$. Each coefficient on a set $`S`$ is divisible by $`\pi^{|S|}`$. Indeed an
exponent monomial $`c\prod_{i\in J}u_i`$ contributes the factor
$`1+(\omega^c-1)\prod_{i\in J}u_i`$. The valuation of its nonconstant coefficient is at least 1,
2 and 4 for $`|J|=1,2,3`$ respectively. When terms are multiplied, the size of the union of
their variable sets is at most the sum of the sizes. The conclusion persists after summation
with any algebraic-integer coefficients.

Consider any nonzero sum $`\alpha`$ of this type on a cube. Let $`v`$ be the minimum valuation
of its nonzero values on that cube. This equals the minimum valuation of its nonzero multilinear
coefficients, by evaluation and integral Möbius inversion. After division by $`\pi^v`$ and
reduction modulo $`\pi`$, the resulting polynomial is nonzero over $`\mathbb F_2`$ and has
degree at most $`v`$. A nonzero binary polynomial of degree $`j`$ takes value one on at least a
$`2^{-j}`$ fraction of the cube, by summing over a maximal monomial's free coordinates. Hence
$`\alpha`$ is nonzero on at least a $`2^{-v}`$ fraction. For any nonzero sample $`\alpha(u_0)`$
we also have
```math
 v\le \log_2\left|N_{\mathbb Q(\omega)/\mathbb Q}\big(\alpha(u_0)\big)\right|.
```
This last bound follows because $`\pi^v`$ divides that sample. It is independent of the number
of variables used to encode the central, commuting part of either layer.

## Replacement cluster estimate

Use the notation of R74 §5.1. There are $`K\ge2`$ signed boundary summands with coefficients
$`\pm2^{-b/2}`$, and $`m=2r`$ active doubled Choi coordinates. For the initialized identity
discrepancy $`K=2^{a+1}`$ and $`b=a`$. Within one cluster at scale $`t`$, the
untouched-coordinate marginals are supported on affine spaces of codimension
$`e_i\le\ell=2(K-1)t`$ in their common affine hull $`H`$. For one active output coordinate $`z`$
the total normalized amplitude in this cluster is $`2^{-\dim H/2}\lambda_z(u)`$. Each term of
$`\lambda_z`$ is its support indicator times $`\pm2^{(e_i-b)/2}`$ times a normalized conditional
amplitude $`v_i(u,z)`$.

Choose an output $`z_0`$ whose amplitude is nonzero at some point of $`H`$ and put
$`\alpha=(\sqrt2)^{b+2m}\lambda_{z_0}`$. A left term is an indicator times a weighted-cubic
phase with prefactor $`\pm(\sqrt2)^{e_i+2m-b_i}`$, where $`b_i\le m`$ is the constant affine
fiber dimension. A right term sums over $`b_i`$ fiber bits with prefactor
$`\pm(\sqrt2)^{e_i+m-b_i}`$. Each fixed path contributes a weighted cubic on affine parameters.
Every embedding $`\sigma`$ preserves normalization of the conditional fiber: it sends phases to
unit phases and the active Hadamards to themselves up to sign. Thus
$`|\sigma(\alpha)|\le K 2^{m+\ell/2}`$ throughout the hull.

Call a polynomial of filtration shift $`d`$ one whose multilinear coefficient on $`S`$ is
divisible by $`\pi^{\max(0,|S|-d)}`$. Products add their shifts since Boolean multilinearization
can only decrease the number of variables. A single affine parity bit or its complement has
shift at most one: for a nonconstant term of degree $`j`$ its integer coefficient is
$`\pm2^{j-1}`$. A conjunction of $`D`$ such tests has shift at most $`D`$.

At the selected output $`z_0`$, work on all of $`H`$. Extend each phase and affine fiber
parametrization from its support to the hull by an affine retraction. In a left term the
projected and active membership bits cost at most $`e_i+m`$ in total shift. Its prefactor
$`(\sqrt2)^{e_i+2m-b_i}`$ has $`\pi`$-valuation $`2(e_i+2m-b_i)\ge e_i+m`$, and hence pays all
these shifts. For a right term only the $`e_i`$ projected membership bits are needed. Its
prefactor has valuation $`2(e_i+m-b_i)\ge e_i`$. Thus $`\alpha`$ has filtration shift zero on
the entire hull.

Let $`v`$ again be the minimum valuation. The mod-$`\pi`$ polynomial shows that the fraction of
$`H`$ attaining precisely $`v`$ is at least $`2^{-v}`$. At each of those points the algebraic
norm satisfies
```math
 2^v \le |N(\alpha)|=|\alpha|^2 |\sigma(\alpha)|^2,\qquad
 |\sigma(\alpha)|^2\le K^2 2^{2m+\ell},
```
where $`\sigma`$ represents the conjugate pair of embeddings other than the identity and complex
conjugation. Since $`|\alpha|^2=2^{b+2m}|\lambda|^2`$, the squared coordinate is at least
$`2^{v-b-4m-\ell}/K^2`$. Its factor $`2^v`$ cancels the lower bound on the size of the
minimum-valuation stratum. We obtain
```math
 \|f_{\rm cluster}\|^2\ge 2^{-b-4m-\ell}/K^2 .
```

## Separating different clusters and transferring the witness

Set $`A=2(K-1)`$, $`C=4m+b+4\log_2K+1`$, $`t_0=0`$, and $`t_{j+1}=At_j+C`$. Among the $`K`$
successive pairs of partitions by projected affine support distance, there is a pair with the
same components. Across those components the sum of all interfering inner products is bounded by
$`K^2 2^{-t_{j+1}}`$, at most half the new nonzero-cluster bound. This follows from the original
affine-marginal overlap bound and the definitions of $`A,C`$. If the entire discrepancy is
nonzero it follows that
```math
 \|f\|^2\ge 2^{-T},\qquad T=C(1+A+\cdots+A^{K-1}).
```
For fixed $`K,b`$, this is linear in $`m`$. With $`m=2r`$, $`K=2^{a+1}`$ and $`b=a`$, the
original initialized preservation and witness argument applies: below $`\min(c_0,2^{-(T+2)/2})`$
every tested occupied Pauli would be preserved exactly. This contradicts the same nonnative
$`R_y`$ witness, where $`c_0=\sqrt{2-(2+\sqrt2)/\sqrt3}`$. Therefore $`T\ge2L-2`$ once
$`2^{-L}\lt c_0`$. As an explicit, conservative instance, $`a=1`$ gives $`K=4`$, $`A=6`$, and
$`T=259(8r+10)`$. The constants do not purport to be sharp at small precisions.
