# A robust linear T-count bound for the complete dirty indicator

[Native sharing requirement](NATIVE_INDICATOR_SHARING.md) ·
[Retained linear-count indicator](NONUNIFORM_DIRTY_INDICATOR.md) ·
[Complete approximate transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md)

Sharing one native T gate across many bank observables can be cheap.
The one-layer example below attains the earlier incidence requirement
with only p T gates, but has the wrong joint Pauli action. Checking the
whole Pauli family supplies a stronger count conclusion: every zero-clean
complete indicator at full operator error less than one over square root
two needs at least m-1 T gates, where m is two to the power p and p is at
least two. Arbitrary dirty width and arbitrary intermediate address/bank
motion are permitted.

This matches the order of the retained O(m) exact T-count. It does not
give a growing T-depth bound, a faster construction, or a change to R96's
complete-frame frontier. The remaining constructive question is how to
supply the complete joint action at that linear count with sublogarithmic
address T-depth.

## 1. The standard Pauli-subspace count argument

For a Q-qubit unitary V, let S(V) be the binary label space of Hermitian
Paulis P for which $`VPV^\dagger`$ is again a signed Hermitian Pauli.
Only the labels are taken modulo scalar; the conjugation itself retains
its literal sign. Products show that S(V) is a binary linear subspace.

If V has a native Clifford+T word containing q T/T-dagger gates, then

```math
\dim S(V)\ge2Q-q.                                      \qquad(1)
```

Here is a direct proof in the present full-space model. Initially all
2Q independent Pauli labels have Pauli images. A Clifford preserves that
property. Immediately before a physical T gate, restrict the currently
tracked subspace to those images commuting with the Z on its target.
This is at most one binary linear condition. Every retained image then
remains exactly the same Pauli after T or T-dagger. Induction leaves a
subspace of dimension at least 2Q-q. It is contained in S(V), even if
additional Pauli images reappear later.

Equation (1) is the familiar unitary stabilizer-nullity T-count method,
not a new general monotone. See Jiang and Wang, *Lower bound for the T
count via unitary stabilizer nullity*, Physical Review Applied 19,
034052 (2023), [Theorem 11](https://arxiv.org/pdf/2103.09999).
Their unitary nullity is $`2Q-\dim S(V)`$. The calculation below
specializes this method to the complete indicator and proves an explicit
constant-error extension using its Pauli-image gap. No continuity of
nullity for arbitrary targets is assumed.

## 2. The complete indicator's exact Pauli subspace

Let x have p bits, Y have $`m=2^p`$ bits, and W have w arbitrary
dirty bits. Put $`Q=p+m+w`$ and use the phase form

```math
D_p|x,Y,W\rangle=(-1)^{Y_x}|x,Y,W\rangle,
\qquad I_p=H_YD_pH_Y.
```

There is no initialized register. All systems may be correlated with
each other and with an external reference. Write an input Pauli's X mask
as (a,v,d) on x,Y,W, respectively. Its Z mask is unrestricted. Conjugation
by D_p leaves its basis permutation unchanged and multiplies its columns
by the Boolean derivative

```math
h_{a,v}(x,Y)=Y_x+Y_{x+a}+v_{x+a}.                      \qquad(2)
```

All additions in (2) are binary. A diagonal sign is a scalar Pauli
exactly when its Boolean exponent is affine. Therefore the input Pauli
belongs to S(D_p) exactly when h is affine in all x,Y variables.

If a is zero, (2) is simply v_x. It is affine precisely when the m-bit
truth table v is the table of an affine Boolean function on p bits.
These v form a space of dimension p+1.

If a is nonzero, choose any bank coordinate u. Its coefficient in (2)
is $`[x=u]+[x=u+a]`$. For p at least two this coefficient is one
at two addresses and zero at another address, so it is nonconstant.
An affine function of all x,Y variables must have constant coefficients
on its Y variables. The term $`v_{x+a}`$ contains no Y and cannot remove
this obstruction. Thus no Pauli with a nonzero address X mask belongs
to S(D_p).

All Q Z-mask directions are free, all w dirty X-mask directions are
free, and the allowed bank X masks have dimension p+1. Consequently

```math
\begin{aligned}
\dim S(D_p\otimes I_W)&=Q+w+p+1\\
                     &=m+2p+1+2w,\\
2Q-\dim S(D_p\otimes I_W)&=m-1.                    \qquad(3)
\end{aligned}
```

This explicitly retains both Pauli directions on every arbitrary dirty
wire. Increasing dirty width does not change the deficit. The condition
p at least two is necessary: the one-bit indicator is Clifford.

## 3. A uniform Pauli-image gap makes the count robust

For every input Hermitian Pauli P outside the space in (3), and every
signed Hermitian Pauli R on the entire physical register,

```math
\|D_pPD_p^\dagger-R\|_{\mathrm{op}}\ge\sqrt2.        \qquad(4)
```

This is a full operator-norm gap, independent of Q. To prove it, write
the target image as $`P D_h`$, where h is the nonaffine derivative
from (2). If R has a different X mask from P, their outputs on each
computational input have distinct destinations, giving column-distance
square two. If their X masks agree, multiply on the left by P-dagger.
The resulting Pauli $`P^\dagger R`$ is diagonal, and has the form

```math
\zeta Z^z,\qquad \zeta\in\{1,-1,i,-i\}.
```

For zeta equal to plus or minus i, every diagonal coefficient has distance
square two from the real sign D_h. For zeta equal to plus or minus one,
nonaffinity of h guarantees a disagreeing sign at some input, with
distance two. This proves (4). It includes the possible imaginary scalar
in the product of two Hermitian Paulis; discarding that scalar would
give an unjustified stronger gap.

Now suppose an actual native V satisfies

```math
\|V-D_p\otimes I_W\|_{\mathrm{op}}
 \le\tau\lt1/\sqrt2.                               \qquad(5)
```

For each input Pauli P in S(V), its image R is a signed Hermitian Pauli,
and unitary telescoping gives

```math
\|R-D_pPD_p^\dagger\|\le2\tau\lt\sqrt2.
```

Equation (4) therefore forces $`S(V)\subseteq S(D_p\otimes I_W)`$.
Combine (1) and (3) to obtain

```math
T(V)\ge m-1.                                          \qquad(6)
```

Hadamards on Y give the same statement for I_p without changing T-count
or T-depth. The argument also survives an optimized common scalar on V,
because that scalar cancels in every conjugation. No input-dependent
phase is allowed. The actual inverse has the same T-count and the same
norm error from the inverse target. Approximate return of x,Y,W is
already included in (5); no exact intermediate return is required.

Equation (6) is unrestricted within this zero-clean unitary model. It
does not assume a row decomposition, read-only address gates, Gaussian
interlayers, bounded interaction arity beyond the native alphabet, or
an exact basis permutation. It does not apply merely from accuracy on
initialized helper columns, and it is not a one-clean frame lower bound.

Together with the retained exact indicator, it establishes the order
$`\Theta_\tau(m)`$ for the standalone T-count at each fixed
$`\tau\lt1/\sqrt2`$. The exact upper construction simultaneously has
O(m) dirty width and Clifford count. Its logarithmic address T-depth
remains to be improved; a linear count alone supplies no growing depth
lower bound when linear width is available.

## 4. Full native incidence in one layer can still have the wrong action

The prior sharing theorem requires at least mp total active incidences
between T gates and the propagated individual bank observables. A literal
native word can meet that condition using only p T gates. Define

```math
P_j=Z_{x_j}\prod_{u=0}^{m-1}Z_{Y_u},\qquad
V_{\rm share}=\prod_{j=1}^p T_{P_j}.                  \qquad(7)
```

These p axes commute and are independent because each has its own address
Z coordinate. Let C be the product of the pm CNOTs from every Y_u into
every x_j. All these CNOTs commute. They preserve Y but temporarily
change x. The chronological word C, a parallel T on all x_j, then the
actual C adjoint implements (7), with the common ledger

```math
T=p,\qquad D_T=1,\qquad G=2pm,\qquad w=0.           \qquad(8)
```

For each bank probe $`X_{Y_u}`$, every one of the p T gates is active.
Its final Pauli expansion contains the entire affine coset generated
by the p independent P_j, with nonzero coefficients, so its affine
support dimension is exactly p. Thus the word has exactly mp individual
bank incidences and meets the earlier per-probe support condition.

Nevertheless, if $`h=\bigoplus_uY_u`$, its basis action is only

```math
V_{\rm share}|x,Y\rangle
 =\omega^{\sum_{j=1}^p(x_j\mathbin\oplus h)}|x,Y\rangle,
\qquad \omega=e^{i\pi/4}.                            \qquad(9)
```

It does not perform selected-row discrimination. For example, at x=0
take two different bank positions equal to one, including position zero.
Then h is zero, the word's phase is one, and the required phase is minus
one. Its literal distance from D_p is therefore two. More generally,
the joint Pauli-subspace test (3)--(6) rules out every sufficiently
accurate p-T replacement when $`p\lt m-1`$.

This example is a component identity and a boundary check on the sharing
criterion. It does not propose a syndrome-only compiler or repeat that
architecture's rank obstruction. It shows concretely why incidence alone
cannot be treated as a depth cost: even a source that satisfies all the
individual bank incidence inequalities may still have the wrong action
on products of those probes and on address Paulis.

One explicit missed joint observable is
$`P_*=X_{x_1}\cdots X_{x_p}X_{Y_0}`$. It commutes with every axis
in (7), and hence with the completed word. Its nonzero address X mask
puts it outside the target Pauli subspace in Section 2.

The [exact finite checks](../../tests/test_indicator_pauli_nullity.py)
verify the complete target Pauli classification on small native examples,
including arbitrary dirty directions and the imaginary-scalar gap case.
They also check (7)--(9), all physical columns and the actual inverse.
The all-width and constant-error conclusions follow from the analytic
arguments above, not from finite enumeration.
