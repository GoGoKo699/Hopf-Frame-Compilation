# Limits of two source-reuse routes

[Open endpoint](OPEN_PROBLEM.md) · [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Shared-source compiler](FAULT_TOLERANT_COMPILER.md)

The open endpoint remains
$`a=2`$, $`b=N+n+7`$, $`L=N`$, $`n\ge3`$, with
$`\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n))`$
by [conditional-suffix grouping](CONDITIONAL_SUFFIX_COMPILER.md).
Here $`\ell_*(n)=1+\log_2^*(n+2)`$, and $`\log_2^*`$ counts
base-two logarithms until the value is at most one.
This note gives two restrictions on proposed ways to reuse precision work.
Neither restriction is a lower bound for an unrestricted Hopf-frame
compiler, and neither changes the retained resource theorems.

## 1. A nilpotent carried source needs initialized dimension

The sufficient-clean compiler uses a truncated shift on an initialized
geometric source. Its useful property is a full-output relation
$`S_d g\approx 2^{-d/2}g`$, not only a small projected overlap.
Could an encoding of arbitrary dirty inputs supply the same property with
only two initialized qubits?

**Lemma 1.** Let $`E:\mathbb C^D\to\mathbb C^{rD}`$ be an isometry,
where $`D,r`$ are positive integers. Let $`K`$ be a nilpotent contraction
on the ambient space, and let $`0\lt\alpha\lt1`$. Then

```math
\varepsilon:=\|KE-\alpha E\|
\ge \frac{\alpha^r}{(1+\alpha)^{r-1}}.
```

There is no restriction on the encoding cost, the nilpotence index, or
the dimension $`D`$ of the arbitrary input.

*Proof.* Put $`A=\alpha I-K`$. The restriction of $`A`$ to the
$`D`$-dimensional range of $`E`$ has norm at most $`\varepsilon`$.
The singular-value min–max principle therefore gives at least $`D`$
singular values of $`A`$ no larger than $`\varepsilon`$. All remaining
singular values are at most $`\|A\|\le1+\alpha`$. Nilpotence gives

```math
\alpha^{rD}=|\det A|
\le \varepsilon^D(1+\alpha)^{(r-1)D}.
```

Taking $`D`$-th roots proves the bound. ∎

If $`a`$ initialized qubits are available, the ambient/input dimension
ratio is $`r=2^a`$. Every arbitrary logical input and dirty helper must
be included in $`D`$. Adding arbitrary dirty wires increases both
dimensions equally and does not improve the ratio. The lemma also
allows an arbitrary initial encoding $`E=UJ_a`$.

For the coefficient shift $`\alpha=1/2`$ used by the geometric source,

```math
\varepsilon\ge\frac32\,3^{-2^a}.
```

In particular, $`a=2`$ gives $`\varepsilon\ge1/54`$, regardless of
dirty width. Requiring $`\varepsilon\le2^{-L}`$ implies

```math
2^a\ge
\frac{L+\log_2(3/2)}{\log_2 3},
\qquad
a\ge\log_2 L-O(1).
```

Thus this nilpotent full-output source interface cannot replace a
precision-dependent initialized source by arbitrarily many dirty wires
while keeping the clean budget constant.

### The clean-width order is attainable for this interface

This is a dimension bound, not a claim of an optimal preparation circuit.
For $`r=2^a`$, take

```math
S_r=\sum_{k=0}^{r-2}|k\rangle\langle k+1|,
\qquad
K=S_r\otimes I_D,
\qquad E=|g\rangle\otimes I_D,
```

with

```math
g_k=
\sqrt{\frac{1-\alpha^2}{1-\alpha^{2r}}}\,\alpha^k.
```

Only the last coordinate contributes to the residual, so

```math
\|KE-\alpha E\|=
\alpha^r\sqrt{\frac{1-\alpha^2}{1-\alpha^{2r}}}.
```

At $`\alpha=1/2`$, choosing $`a=\lceil\log_2 L\rceil`$ suffices for
error at most $`2^{-L}`$, for integer $`L\ge1`$. Hence the
$`\Theta(\log L)`$ initialized-width order is matched as a
linear-algebra statement. This normalized vector is not asserted to
have an exact Clifford+T preparation. The native capped source in the
[shared-source proof](FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel)
also has this width order.

### Applicability boundaries

The lemma applies to the carried contraction actually satisfying the
displayed full-output relation. Its nilpotence must be proved.
A compression of a nilpotent operator need not be nilpotent.
Likewise, projecting flags in a unitary dilation does not automatically
give a nilpotent contraction.

The full residual kernel of the shared-source compiler contains an
identity component; it is not the nilpotent shift itself. Lemma 1
cannot be applied to that kernel without a new argument.
A known-zero logical suffix may provide initialized dimension within
a restricted sector, but its coherent inactive-sector action must then
be proved separately. Under the complete-frame contract, the logical
register is arbitrary and is not free initialized workspace.

The lemma leaves nonnilpotent sources, other operator representations,
and a directly constructed whole-frame block open.

## 2. Factoring out the source can move its cost into the masks

Use the operator-source notation

```math
R_k=\exp(i\pi Y_kX_{k+1}/8),\qquad
U_m=R_{m-2}\cdots R_0,\qquad
M_m=U_mX_0U_m^\dagger.
```

When a common source width is used, conjugating a complete scalar block
and its amplification by $`U_m^\dagger`$ replaces each source by $`X_0`$.
Since $`U_m`$ acts only on the dirty core, it commutes with the
flag reflections and logical operations. Adjacent outer basis changes
can then cancel between stages.

The programmed Pauli masks, however, become

```math
\widetilde P_f=U_m^\dagger P_fU_m.
```

They are no longer necessarily Clifford. The following valid single-bit
mask quantifies the issue.

**Lemma 2.** For $`0\le j\le m-2`$, let
$`D_j=U_m^\dagger Z_jU_m`$. Its minimum exact T-count is

```math
T_{\rm exact}(D_j)=2j.
```

Arbitrary returned clean stabilizer helpers and arbitrary returned dirty
helpers are allowed.

*Upper bound.* The source rotations with index above $`j`$ commute with
$`Z_j`$, while

```math
R_j^\dagger Z_jR_j
=B_j:=\frac{Z_j+X_jX_{j+1}}{\sqrt2}
=\mathrm{CNOT}_{j\to j+1}H_j\mathrm{CNOT}_{j\to j+1}.
```

Thus

```math
D_j=(R_0^\dagger\cdots R_{j-1}^\dagger)
B_j(R_{j-1}\cdots R_0).
```

The middle factor is Clifford, and each of the $`2j`$ remaining
Pauli rotations costs one T or T-dagger. Common scalars cancel.

*Lower bound.* With the Pauli probe $`P=Z_0`$, expansion of this
word gives

```math
c_j=2^{-m}\operatorname{Tr}(PD_jPD_j)=1-2^{-j}.
```

For $`j\ge1`$, exactly one term in the anticommuting Pauli expansion
of $`D_j`$ anticommutes with $`P`$: it is
$`-Y_0Z_1\cdots Z_{j-1}Y_j`$ with coefficient
$`2^{-(j+1)/2}`$. All other terms commute with $`P`$.
This proves the displayed correlation. The case $`j=0`$ follows
directly from $`B_0`$.

For $`j\ge1`$, the rational $`1-2^{-j}`$ has least
square-root-of-two denominator exponent $`2j`$. The
[returned-helper Pauli-transfer argument](OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers)
therefore requires at least $`2j`$ T gates. For $`j=0`$ the
Clifford upper bound is already optimal. ∎

The choice $`j=m-2`$ is allowed by the dyadic encoding, which fixes
only the final bit to zero. This single transformed mask costs
$`2m-4`$ T gates exactly.

### The obstruction persists at sufficiently fine approximation

Suppose a circuit $`V`$ with $`t`$ T or T-dagger gates obeys

```math
\|VJ-J(D_j\otimes I_s)\|\le\epsilon,
\qquad \epsilon\le2^{-j-2},
```

where $`s`$ arbitrary dirty helpers and any initialized stabilizer
helpers are included in the full-isometry contract. Then

```math
t\ge\max\{0,j-1\}.
```

Here is a proof that also accounts for approximate helper return.
After a Clifford change of initialization, take the clean helpers
to be zero. Set $`d=m+s`$, $`C=Z_0\otimes I_s`$,
$`A=C\otimes I_{\rm clean}`$, $`F=VJ`$, and

```math
c=2^{-d}\operatorname{Tr}(CF^\dagger AF).
```

The complete-isometry estimate gives
$`|c-(1-2^{-j})|\le2\epsilon`$.
The Pauli-transfer denominator argument gives
$`c\in(\sqrt2)^{-t}\mathbb Z[\sqrt2]`$.

Let $`\sigma`$ fix $`i`$ and send $`\sqrt2`$ to
$`-\sqrt2`$. Applied to a native circuit, it sends
$`H`$ to $`-H`$ and $`T^{\pm1}`$ to $`ZT^{\pm1}`$;
$`S`$ and CNOT are unchanged. Therefore $`\sigma(V)`$ is
also unitary. The initialization and Pauli probes above are fixed by
$`\sigma`$, so $`c'=\sigma(c)`$ is another normalized
Hermitian-channel correlation. In particular,
$`c,c'\in[-1,1]`$. No return promise for $`\sigma(V)`$ is needed.

If $`c\ne1`$, the nonzero quadratic algebraic norm satisfies

```math
|(1-c)(1-c')|\ge2^{-t},
\qquad
1-c\ge2^{-t-1}.
```

The first inequality follows by writing
$`1-c=z/(\sqrt2)^t`$ with nonzero
$`z\in\mathbb Z[\sqrt2]`$ and taking its integer norm.
On the other hand, the assumed accuracy gives

```math
2^{-j-1}\le1-c\le3\cdot2^{-j-1}.
```

Consequently $`t\ge j-\log_2 3`$, and integrality proves the claim.
At $`j=m-2`$ and $`\epsilon\le2^{-m}`$, this gives
$`t\ge m-3`$.

This is a lower bound for one transformed mask at the stated accuracy.
It is not permissible to add these primitive costs and infer a
whole-frame lower bound. Joint synthesis could avoid treating the masks
as independent operations.

## 3. Consequence for the research direction

These results rule out two particular shortcuts: replacing the
initialized nilpotent geometric source by a constant-clean dirty
encoding with the same full-output property, and treating all
source-conjugated programming masks as cheap Clifford operations.

The sufficient [whole-frame half-unitary block](OPEN_PROBLEM.md#a-sufficient-construction-to-seek)
remains a valid target. A successful construction may avoid intermediate
source return, may use a different carried operator, or may synthesize
the interleaved source and programming operations jointly.
The bounds for unrestricted complete-frame compilation remain unchanged.

## 4. Finite checks and evidence limits

Run:

```bash
python -m unittest discover -s tests -p 'test_source_reuse_limits.py'
```

The tests exercise the dimension inequality on small nilpotent
contractions and encoded subspaces, the explicit geometric construction,
the need for nilpotence, and the transformed-mask correlations.
These finite checks support indexing and assumption boundaries.
The dimension and T-count statements rest on the analytic proofs above;
the tests do not establish an unrestricted impossibility theorem or
literature priority.
