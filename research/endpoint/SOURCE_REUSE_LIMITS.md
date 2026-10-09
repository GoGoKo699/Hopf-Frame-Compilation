# Source reuse and structured residuals

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Open endpoint](../../docs/OPEN_PROBLEM.md) · [Operator-source compiler](../../docs/OPERATOR_SOURCE_COMPILER.md) · [Shared-source compiler](../../docs/FAULT_TOLERANT_COMPILER.md)

The open endpoint remains
$`a=2`$, $`b=N+n+7`$, $`L=N`$, $`n\ge3`$, with
$`2N-1\le T^\star_{F,\mathbb R}\le O(N\ell_*(n))`$
by the [literal-phase lower bound](../../docs/FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness)
and [conditional-suffix grouping](../../docs/CONDITIONAL_SUFFIX_COMPILER.md).
Here $`\ell_*(n)=1+\log_2^*(n+2)`$, and $`\log_2^*`$ counts
base-two logarithms until the value is at most one.
This note gives scoped restrictions on reusing precision work and encoding
the frame, an exact tree factorization of residual coefficients, and a
Boolean frame family for which exactly returned dirty work reduces T-count.
The factorization separates classical compression from the still-charged
task of coherent tree transport. These interface results do not establish
an unrestricted superlinear lower bound or a linear-T frame compiler.
The [tree-transport continuation](ENDPOINT_TREE_TRANSPORT.md) gives an
exact sparse representation, normalized unitary columns, height-independent
weighted norms, and a scoped obstruction to finite-order coarse corrections.

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
[shared-source proof](../../docs/FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel)
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
c_j=2^{-m}\,\mathrm{Tr}(PD_jPD_j)=1-2^{-j}.
```

For $`j\ge1`$, exactly one term in the anticommuting Pauli expansion
of $`D_j`$ anticommutes with $`P`$: it is
$`-Y_0Z_1\cdots Z_{j-1}Y_j`$ with coefficient
$`2^{-(j+1)/2}`$. All other terms commute with $`P`$.
This proves the displayed correlation. The case $`j=0`$ follows
directly from $`B_0`$.

For $`j\ge1`$, the rational $`1-2^{-j}`$ has least
square-root-of-two denominator exponent $`2j`$. The
[returned-helper Pauli-transfer argument](../../docs/OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers)
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
c=2^{-d}\,\mathrm{Tr}(CF^\dagger AF).
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

### The current paired source also has expensive transformed masks

Moving a dirty-only basis change entirely outside a word does not
improve its accepted-block accuracy. Write
$`Q=(I\otimes U)K(I\otimes U^\dagger)`$, let J initialize only
the signal flags, and put $`B=J^\dagger KJ`$. For a desired
logical operator E,

```math
\bigl\|J^\dagger QJ-E\otimes I_{\rm dirty}\bigr\|
=\bigl\|B-E\otimes I_{\rm dirty}\bigr\|.
```

The full-isometry identity is equally exact:

```math
\|QJ-J(E\otimes I_{\rm dirty})\|
=\|KJ-J(E\otimes I_{\rm dirty})\|.
```

The dirty-only U passes through J as an input-dirty unitary and commutes
with the returned-work target. Unitary invariance proves both equalities.
If the interior K were independently cheap, the outer basis changes could
be omitted; precision needed for the target must survive inside K.

The [one-clean compiler](../../docs/ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks)
uses a different, paired source. Its programmed masks have the same
limitation; no transfer of the preceding one-tail formula is assumed.
Let q be its precision parameter, $`m=q+1`$, and write its actual
loader as $`U_q`$, so $`M=U_qX_0U_q^\dagger`$. Choose any mask
$`P_g`$ supplied by its certified programming rule for
$`\theta=\pi/3`$, and set

```math
D_g=U_q^\dagger P_gU_q,\qquad
r=\frac{\beta}{2}=\frac{\sqrt5-1}{8},
\qquad \beta=\sin(\pi/10).
```

**Corollary.** For $`q\ge5`$, a native circuit approximating
$`D_g`$ with full-isometry error at most $`2^{-q}`$, including
return of arbitrary dirty and initialized stabilizer helpers, requires

```math
t\ge\max\{0,q/2-6\}
```

T or T-dagger gates.

*Proof.* The exact normalized Pauli correlation is

```math
2^{-m}\,\mathrm{Tr}(X_0D_g^\dagger X_0D_g)
=2^{-m}\,\mathrm{Tr}(MP_g^\dagger MP_g)=s_g.
```

The paired-source coefficient bounds give
$`|s_g-r|\le(15/8)2^{-q}`$. Use the same returned-helper
Pauli-transfer argument as above, now with probe $`X_0`$. The
approximating circuit's correlation c and its Galois conjugate obey

```math
c\in(\sqrt2)^{-t}\mathbb Z[\sqrt2],\qquad
c,c'\in[-1,1],\qquad
|c-r|\le\frac{31}{8}2^{-q}.
```

Put $`f(x)=16x^2+4x-1`$. Its two roots are outside
$`\mathbb Q(\sqrt2)`$, so $`f(c)\ne0`$. Since
$`f(c)\in2^{-t}\mathbb Z[\sqrt2]`$, its nonzero algebraic norm
has absolute value at least $`2^{-2t}`$. Meanwhile $`f(r)=0`$,
$`|f'|\le36`$ on $`[-1,1]`$, and $`|f(c')|\le21`$.
Consequently

```math
2^{-2t}\le|f(c)f(c')|
\lt140\cdot21\,2^{-q}\lt2^{12-q},
```

which proves the claim. The conjugate circuit is unitary regardless of
whether it returns its helpers. ∎

Factoring a common $`U_q`$ outside successive words therefore does
not make the intervening transformed masks free or constant-cost
Cliffords. This corollary prices one separately implemented mask.
Its costs cannot be added to obtain a depth-dependent compiler lower
bound.

### Exact hoisting of the paired source

The actual cancellation can be displayed without assigning independent
costs to its transformed masks. At one fixed precision let
$`M=UX_0U^\dagger`$ and let $`P_g`$ be the complete programmed mask.
The scalar word is exactly $`\mathcal S_g=U\mathcal K_gU^\dagger`$,
where

```math
\mathcal K_g
=H_r C_{r=0}(X_0)
 (U^\dagger P_gU)X_0(U^\dagger P_g^\dagger U)
 C_{r=1}(X_0)H_r.
```

All logical/signal routing and amplification reflections commute with
the unconditional dirty-only U. Thus a stream of real primitives using
this core and precision cancels its common exterior
$`U^\dagger U`$ pairs, including when the logical target changes.
The self-borrowed masks and core-borrowed predicate helpers have the
same complete operators as their external-work versions, so they
preserve this identity.

For the fixed mask f of the one-clean compiler, all tail rotations
commute with $`P_f`$: each joins Majoranas with equal f signs. Writing
$`U=U_{\rm tail}U_{\rm seed}`$ leaves

```math
U^\dagger P_fU=U_{\rm seed}^\dagger P_fU_{\rm seed},
```

which uses only the two seed rotations and their inverses. The arbitrary
programmed mask g retains $`U^\dagger P_gU`$. Each real primitive
contains five programmable scalar occurrences after amplification and
hence ten such transforms. Since U uses $`2q`$ T gates, the displayed
hoisted word for s primitives has source-rotation charge at most
$`(40s+4)q+O(s)`$, in addition to its queries and predicate gates.
This is a direct upper certificate for that word, not an additive lower
bound against a joint rewrite of the coupled transformed masks.

## 3. Tree generators compress the residual classically

The column-forest block stores coefficients indexed by a marker and an
ancestor depth. A single s-level group has $`\Theta(s2^s)`$ permitted
entries, although its two underlying frames have only $`O(2^s)`$ local
words. The extra entries are not independent classical parameters. The
following factorization records their exact dependence without dividing by
a subtree amplitude.

### Node overlaps and path products

Index a binary tree by heap labels: its root is 1, the children of v are
$`2v,2v+1`$, and its N leaves are $`N,\ldots,2N-1`$. For a frame
$`A\in\{C,W\}`$, let $`U_v^A`$ be the local two-by-two unitary. The
target words $`U_v^W`$ are real rotations; the actual coarse words
$`U_v^C`$ may be complex Clifford+T unitaries. Define the normalized
subtree vector and its prescribed complement by

```math
\begin{pmatrix}|q_v^A\rangle&|e_v^A\rangle\end{pmatrix}
=\begin{pmatrix}|q_{2v}^A\rangle&|q_{2v+1}^A\rangle\end{pmatrix}U_v^A.
```

The leaf vectors are the corresponding computational basis vectors.
The complete frame has first column $`q_1^A`$ and the complement
$`e_v^A`$ in its prescribed marker column. Below, the symbol $`\ast`$
denotes the first column and a node label v denotes that marker column;
this is only a common reindexing of the two frames.

**Lemma 3 — tree generators for the residual.** Set $`g_v=1`$ at every
leaf. At every internal node compute the two-by-two matrix

```math
K_v=(U_v^C)^\dagger
\begin{pmatrix}g_{2v}&0\\0&g_{2v+1}\end{pmatrix}U_v^W
=\begin{pmatrix}g_v&k_v\\h_v&d_v\end{pmatrix}.
```

These $`N-1`$ fixed-size matrices, together with the local words, determine
every entry of $`E=C^\dagger W-I`$.

For internal nodes v and u, let the path from v to its strict descendant u have edge
bits $`b_0,\ldots,b_{t-1}`$ and nodes
$`v_0=v,\ldots,v_t=u`$. Define

```math
\beta_A(v,u)=(U_{v_0}^A)_{b_0,1}
\prod_{j=1}^{t-1}(U_{v_j}^A)_{b_j,0}.
```

Let $`a_A(u)`$ be the product of first-column edge coefficients along
the root-to-u path, with $`a_A(1)=1`$. Then

```math
\begin{aligned}
E_{\ast,\ast}&=g_1-1,& E_{u,u}&=d_u-1,\\
E_{u,\ast}&=a_W(u)h_u,&
E_{\ast,u}&=\overline{a_C(u)}\,k_u,\\
E_{u,v}&=\beta_W(v,u)h_u,&
E_{v,u}&=\overline{\beta_C(v,u)}\,k_u
\quad(u\text{ strictly below }v).
\end{aligned}
```

Entries between incomparable internal nodes are zero.

*Proof.* The two child subtrees have disjoint supports. Their overlap
matrix in the pairs $`(q_{2v},q_{2v+1})`$ is therefore
$`\mathrm{diag}(g_{2v},g_{2v+1})`$. Applying the two local words gives
the displayed recurrence, with
$`g_v=\langle q_v^C|q_v^W\rangle`$,
$`h_v=\langle e_v^C|q_v^W\rangle`$,
$`k_v=\langle q_v^C|e_v^W\rangle`$, and
$`d_v=\langle e_v^C|e_v^W\rangle`$.
The restriction of an ancestor complement $`e_v^A`$ to u's subtree is
exactly $`\beta_A(v,u)q_u^A`$; the restriction of the root state is
$`a_A(u)q_u^A`$. Taking the corresponding inner products proves every
listed entry. Incomparable supports are disjoint. ∎

This is an $`O(N)`$-scalar description and an $`O(N)`$ count of fixed-size
classical recurrence operations. It is not an $`O(N)`$ bit-complexity
claim at precision L. An individual path product still has up to n
factors. All formulas remain valid when a branch amplitude is zero, so
the factorization retains the prescribed marker columns at singular
charts rather than reconstructing a frame from its first column alone.

### Two and three levels expose the target transport

For two levels, use the actual marker order $`0,1,2,3`$. Let the target
and coarse root angles be $`\alpha`$ and $`\gamma`$, respectively, and
let the left child have target-minus-coarse angle $`\delta`$. Other
child angles are arbitrary. For real coarse rotations, the four entries
are exactly

```math
\begin{aligned}
E_{1,0}&=\sin\delta\cos\alpha,&
E_{1,2}&=-\sin\delta\sin\alpha,\\
E_{0,1}&=-\cos\gamma\sin\delta,&
E_{2,1}&=\sin\gamma\sin\delta.
\end{aligned}
```

Replacing target transport by coarse transport in the first row pair
changes that two-entry vector by Euclidean norm

```math
2\left|\sin\delta\,
\sin\!\left(\frac{\alpha-\gamma}{2}\right)\right|.
```

The error can be second order in small coarse discrepancies, but it is
not identically zero. The real coarse specialization isolates the
algebra; the lemma itself allows actual complex native words.

Continuing the real-coarse specialization at three levels, let target root and left-child angles be
$`\alpha_0,\alpha_1`$, their coarse counterparts
$`\gamma_0,\gamma_1`$, and let the left-left internal node have angle
discrepancy $`\delta`$. Its marker is 1; the left-child and root markers
are 2 and 4. Regardless of the other local words,

```math
\begin{aligned}
E_{1,0}&=\sin\delta\cos\alpha_0\cos\alpha_1,\\
E_{1,2}&=-\sin\delta\sin\alpha_1,\\
E_{1,4}&=-\sin\delta\sin\alpha_0\cos\alpha_1,\\
E_{0,1}&=-\sin\delta\cos\gamma_0\cos\gamma_1,\\
E_{2,1}&=\sin\delta\sin\gamma_1,\\
E_{4,1}&=\sin\delta\sin\gamma_0\cos\gamma_1.
\end{aligned}
```

Thus the same local residual coefficient multiplies different target
path products in the downward entries and coarse path products in the
reverse entries. Deeper trees extend these products.

### What classical compression does and does not supply

The generator description identifies a concrete representation to exploit;
it does not supply a charged coherent evaluator for it. In particular:

- Loading all already-expanded ancestor coefficients through the existing
  table primitive still has $`\Theta(nN)`$ permitted rows for one full
  n-level group. The $`O(N)`$ generator count alone does not change that
  query circuit.
- Evaluating a product from loaded local generators requires coherent
  arithmetic or another explicit block construction, including its work,
  precision, inverse, and dirty-return costs. Classical preprocessing
  cannot stand in for this quantum operation.
- Preparing the downward path amplitudes by invoking the target subtree
  frame would call part of the frame that is being compiled. The target
  factors cannot simply be replaced by the coarse ones; the two-level
  witness displays the omitted term.

A giant group also retains a separate charged task: the current streamed
coarse interpreter uses $`w=O(n)`$ symbols at its addressed rows. Its
retained estimate is $`O(nN)`$, even if the residual table were compressed.
This estimate uses the uniform-star representation's fine coarse accuracy.
For the weighted representation, the
[retained borrowed compiler already supplies a cheap coarse frame](ENDPOINT_TREE_TRANSPORT.md#the-retained-borrowed-compiler-already-supplies-a-cheap-coarse-frame):
at constant coarse accuracy and the endpoint dirty allocation it gives
$`T(C)=O(\sqrt N)`$, $`G(C)=O(N)`$, with exactly returned helpers.
This prices C, not the weighted transport or its filters. An improved joint
construction must still charge every component. These are limitations of
the displayed routes, not a lower bound against another compiler.

The [explicit transport construction](ENDPOINT_TREE_TRANSPORT.md) now
realizes these path columns by local three-mode unitaries without a depth
register. Their Gram matrix is known exactly. After the h/k weights are
included, a finite-tree embedding estimate removes the height factor from
their operator norms. The standalone implementation still retains separate
precision-bearing stages, and a fixed-order perturbation around coarse
transport misses necessary mixed terms. Thus improving the endpoint requires
reducing the precision charge of an actual weighted-operator circuit,
including occupied flags and its literal inverse; classical compression
alone does not do this.
The [forward weighted-block continuation](WEIGHTED_TRANSPORT_BLOCK.md)
now provides a complete one-signal-flag dilation with native synthesis using
only borrowed work: $`T=O(N+nL)`$, $`G=O(NL)`$, and
$`b\ge L+n+7`$, with approximate core and amplification-signal return included in the error.
Its $`O(N\log N)`$ endpoint T-count still misses the linear target;
it is a component construction, not a complete frame compiler. The earlier
$`O(L\sqrt N)`$ route retains exact dirty return at its stated width.

## 4. Consequence for the research direction

Sections 1–2 rule out two particular shortcuts: replacing the
initialized nilpotent geometric source by a constant-clean dirty
encoding with the same full-output property, and treating all
source-conjugated programming masks as cheap Clifford operations.

The sufficient [whole-frame half-unitary block](../ROUTE_HISTORY.md#a-sufficient-construction-to-seek)
remains a valid target. A successful construction may avoid intermediate
source return, may use a different carried operator, or may synthesize
the interleaved source and programming operations jointly.
The best retained complete-frame bound remains unchanged. The
[affine residual assembly](RESIDUAL_ASSEMBLY.md) supplies an actual two-flag
word at $`O(N+nL)`$ T cost by merging the diagonal with the forward term.
This resolves its composition and flag-allocation step; reducing its repeated
precision cost remains a separate construction task.

The tree factorization above provides a compact classical starting point
for another construction, while keeping its coherent implementation as an
explicit unresolved task. A linear generator count is not yet a linear
T-count.

### Directly merging the scalar and rotation flags fails

The [two-flag rotation block](../../docs/OPERATOR_SOURCE_COMPILER.md#4-two-flags-encode-a-suffix-controlled-rotation)
uses separate flags for scalar symmetrization and cosine/sine selection.
Consider identifying these flags in the direct two-branch word. Fix one
address row, put $`J=XZ=-iY`$ on the logical target, and let the dirty
programmed sources $`N_c,N_s`$ encode the real coefficients $`c,s`$.
Define

```math
D_c=\frac{MN_c-N_cM}{2},\qquad
D_s=\frac{MN_s-N_sM}{2}.
```

The source identities give $`D_r^\dagger=-D_r`$ and
$`D_r^\dagger D_r=(1-r^2)I`$ for $`r=c,s`$. The merged word is a
Hadamard on the single flag, followed by branch actions $`MN_c`$ and
$`JN_sM`$, followed by a Hadamard. Its accepted block is therefore

```math
B_{\rm merge}
=\frac{MN_c+JN_sM}{2}
=\frac{cI+sJ+D_c-JD_s}{2}.
```

The unwanted dirty operators do not cancel. With
$`E=B_{\rm merge}-(cI+sJ)/2`$, tracing over the target eliminates the
cross terms because $`\mathrm{Tr}(J)=0`$:

```math
\frac12\mathrm{Tr}_{\rm target}(E^\dagger E)
=\frac{2-c^2-s^2}{4}I,
\qquad
\|E\|\ge\frac12\sqrt{2-c^2-s^2}.
```

The inequality follows because normalized partial trace is a unital
positive map. For true cosine and sine coefficients it gives
$`\|E\|\ge1/2`$; increasingly accurate digit tables retain this
constant defect. Even the exact identity-angle row fails under the usual
cubic amplification: $`c=1,s=0`$ gives the projector
$`B_{\rm merge}=(I-JD_s)/2`$. Its amplified accepted block is

```math
3B_{\rm merge}-4B_{\rm merge}B_{\rm merge}^\dagger B_{\rm merge}
=-B_{\rm merge},
```

which is not identity.

This calculation excludes only this direct merged-flag word. The
[one-clean compiler](../../docs/ONE_CLEAN_COMPILER.md) uses a different conjugated
scalar-source word and cancels the dirty terms through a second overlap.
Thus the diagnostic is not a one-clean lower bound; the constructive
replacement is now explicit.

### Moving a query past its address change is not cancellation

Let x be an address bit and z a dirty core bit. The valid addressed
Pauli mask $`P=\mathrm{CZ}_{x,z}`$ satisfies

```math
P R_y(\theta)_xP^\dagger
=\exp(-i\theta Y_xZ_z).
```

On dirty input $`|1\rangle_z`$ this applies
$`R_y(-\theta)`$, with error $`2|\sin\theta|`$ from the intended
$`R_y(\theta)`$. The dirty bit itself is returned. Thus delaying a
mask's inverse past a gate that changes its address can leave a logical
error even with exact dirty return. A fused word must unload before that
change or explicitly account for the conjugated operation.

Likewise, sharing a signal does not make accepted scalar blocks multiply.
For the current paired source's fixed mask,

```math
\mathcal S_f=\frac12 I+X_{\rm flag}D_f,\qquad
D_f^2=-\frac34I,\qquad
\langle0|\mathcal S_f^2|0\rangle=-\frac12I.
```

The product of its two accepted coefficients would instead be
$`I/4`$. This exact $`3/4`$ discrepancy is independent of q.
It excludes concatenation of these scalar words as a substitute for a
new fusion identity; it does not exclude the complete affine assembly
or another choice of rejected-space action.

### A scoped diagnostic for Pauli routing of rejected components

A separate proposed shortcut routes successive scalar-filter rejections
through the four states of two flags. For the common operator source M,
the [scalar block](../../docs/OPERATOR_SOURCE_COMPILER.md#3-dirty-programming-and-a-one-flag-scalar-block)
has the form

```math
\mathcal S_f=c_f I+X_{\rm flag}\otimes A_f,
\qquad A_f=\frac{MN_f-N_fM}{2}.
```

The source identities imply

```math
A_f^\dagger=-A_f,\qquad
A_f^\dagger A_f=(1-c_f^2)I,\qquad MA_fM=-A_f.
```

Suppose each rejection X is replaced by a Hermitian two-flag Pauli
$`P_j`$. Its flip vector records which flag bits it changes. To prevent
each one-step return to $`00`$, every flip vector must be nonzero; to
prevent every two-step return, they must be distinct. Three such vectors
exhaust the nonzero vectors in $`\mathbb F_2^2`$. Their sum is zero,
so

```math
\xi=\langle00|P_3P_2P_1|00\rangle,
\qquad |\xi|=1.
```

Pauli signs or diagonal Pauli phases do not remove this three-step return.
For an explicit coherent witness on three logical qubits, let
$`\Pi_j`$ project onto the respective ordered basis pairs
$`(0,4),(4,6),(6,7)`$, and let $`F_j`$ apply the literal
$`R_y(\pi/2)`$ on that pair and identity elsewhere. The actual unitaries

```math
Q_j=I_{\rm flags}\otimes B_j\otimes I_{\rm dirty}
+P_j\otimes F_j\Pi_j\otimes A_j,
\qquad B_j=F_j(I-\Pi_j+c_j\Pi_j)
```

are a rotation after a conditional scalar filter. The two-stage accepted
block is $`B_2B_1\otimes I`$, but the three-stage block is

```math
B_3B_2B_1\otimes I
+\xi\,(F_3\Pi_3F_2\Pi_2F_1\Pi_1)\otimes A_3A_2A_1.
```

The logical mixed factor equals $`|111\rangle\langle000|`$. Its dirty
factor is odd under conjugation by M, whereas every logical operator
tensored with the dirty identity is even. The odd projection
$`X\mapsto(X-(I\otimes M)X(I\otimes M))/2`$ is norm-contractive.
Consequently the three-stage block has distance at least

```math
\prod_{j=1}^3\sqrt{1-c_j^2}
```

from **every** operator of the form $`B\otimes I_{\rm dirty}`$.
For any source width $`m\geq3`$, the valid single-bit mask $`f=e_1`$
gives $`c_j=1/2`$ and the defect is at least $`3\sqrt3/8`$.
Increasing precision does not shrink this particular defect.

This diagnostic concerns Pauli routing of these scalar coefficient
filters. It is not a lower bound for arbitrary use of two clean flags,
and these filters are not the complete amplified half-unitary frame
blocks. More general history couplings can avoid this three-stage return;
doing so alone would still not combine the separately charged source
calls. The [native routing fixtures](../../tests/test_source_merge.py)
check the explicit word and the stated class boundary.

## 5. A shared conjugator does not close a branching fork

There is a native source-cancellation proposal for which the outer
conjugators really do cancel. Its failure is therefore more specific than
an invalid cancellation rule: the remaining full word has unwanted
rejected returns, and even arbitrary retuning of its masks cannot give a
generic normalized fork.

Use the paired source and scalar words of the
[one-clean construction](../../docs/ONE_CLEAN_COMPILER.md#3-conjugating-scalar-blocks-produces-a-rotation).
Its fixed mask has moment $`c=1/2`$. Write the corresponding
anti-Hermitian dirty operators as $`D_f,D_g`$, and put

```math
s=\langle M,N_g\rangle,\qquad
\tau=\langle N_f,N_g\rangle,\qquad r=s/2-\tau.
```

The brackets are real coefficient-vector inner products in the Majorana
basis. Thus $`\{D_f,D_g\}=2rI`$ and $`D_f^2=-3I/4`$.
Let Z be parity on the logical register. Each physical Hopf pair e differs
in one bit. Its symmetric pair swap $`X_e`$, extended by zero off that
pair, anticommutes with Z and obeys $`X_e^2=\Pi_e`$.

Define one common full-unitary conjugator and edge words by

```math
\begin{aligned}
A&=\frac12I+X_f ZD_f,\\
B_e&=I-\Pi_e+s_e\Pi_e+X_fX_eD_{g_e},\\
Q_e&=A^\dagger B_eA.
\end{aligned}
```

Here f denotes the signal bit when it appears in $`X_f`$, and the fixed
mask label when it appears in $`D_f`$. Dirty and logical identities are
implicit. A is the fixed native scalar word conjugated by flag-controlled
logical parity. Each B is a predicate-conditioned native scalar word with
the original flag-to-target CNOT routing; it is identity off its pair.
Queries are unloaded inside their scalar word before a later operation
changes the address. Thus the whole-space cancellation

```math
Q_k\cdots Q_1=A^\dagger B_k\cdots B_1A
```

is valid, including arbitrary dirty states and occupied signal ports.

Expanding one complete word gives

```math
Q_e=E_e+X_fX_eK_e,\qquad
E_e=I-\Pi_e+s_e\Pi_e+r_eX_eZ,
```

```math
K_e=D_{g_e}+2r_eD_f,\qquad
K_e^\dagger=-K_e,\qquad
K_e^\dagger K_e=(1-s_e^2-r_e^2)I.
```

Consider the physical two-level fork with basis labels
$`a=00,b=10,c=01,d=11`$: root pair $`(a,b)`$, left child
$`(a,c)`$, and right child $`(b,d)`$. Keeping every rejected return,

```math
\begin{aligned}
\langle0|Q_RQ_LQ_0|0\rangle
={}&E_RE_LE_0\otimes I\\
&+|c\rangle\langle b|\otimes K_LK_0
+|d\rangle\langle a|\otimes K_RK_0.
\end{aligned}
```

Indeed $`X_RX_L=0`$, $`X_LX_0=|c\rangle\langle b|`$,
and $`X_RE_LX_0=|d\rangle\langle a|`$. If
$`\kappa_e=\sqrt{1-s_e^2-r_e^2}`$, the difference from
$`E_RE_LE_0\otimes I`$ has norm exactly
$`\max\{\kappa_L\kappa_0,\kappa_R\kappa_0\}`$, because
its two input/output supports are orthogonal and each dirty product has
constant singular value. Under the original angle programming this tends
to $`1-\beta^2`$, where $`\beta=\sin(\pi/10)`$.

This norm is not automatically a distance from every logical operator:
the dirty products can contain scalar parts. Independently, the product
of ideal local accepted blocks has singular values
$`\beta,\beta,\beta^2,\beta^2`$; inactive local modes were identity.
It is not one uniformly normalized frame.

### Retuning the masks does not repair this word

The same word has a precision-independent obstruction even if its local
masks are chosen jointly rather than by the original angle rule.
The Clifford vectors
$`F=N_f-M/2`$ and $`G=N_g-sM`$ have squared lengths
$`3/4`$ and $`1-s^2`$, with inner product $`-r`$.
Cauchy--Schwarz therefore gives the feasible ellipse

```math
s^2+\frac43r^2\le1.
```

Let $`V=R_R(\theta_R)R_L(\theta_L)R_0(\theta_0)`$ be the real
target fork, with the conventional orientation on each ordered pair.
The right pair's parity reverses its encoded sine sign, which does not
change the ellipse. For any $`\gamma>0`$, define

```math
\epsilon=
\left\|\langle0|Q_RQ_LQ_0|0\rangle
-\gamma V\otimes I_{\rm dirty}\right\|.
```

Then every choice of the three masks satisfies

```math
\epsilon\ge
\frac{\gamma|\cos\theta_L|}{2}
\left(1-\sqrt{1-\frac{\sin^2\theta_0}{4}}\right).
```

To prove this, set $`t=\gamma\cos\theta_L`$,
$`v=(s_0,r_0)`$, and $`u=(\cos\theta_0,\sin\theta_0)`$.
The spectator entry $`(c,c)`$ and the row a restricted to columns a,b
contain no rejected-return term. They give, respectively,

```math
|s_L-t|\le\epsilon,\qquad
\|s_Lv-tu\|_2\le\epsilon.
```

Since $`\|v\|_2\le1`$, these imply
$`|t|\|v-u\|_2\le2\epsilon`$. The support function of the
feasible ellipse in the unit direction u is
$`\sqrt{1-\sin^2\theta_0/4}`$. Hence
$`\|v-u\|_2\ge1-\sqrt{1-\sin^2\theta_0/4}`$, as claimed.
The bound is trivial when t is zero and requires no division by it.

For $`\gamma=1/2`$ and $`\theta_0=\theta_L=\pi/4`$, the
lower bound exceeds 0.011, regardless of source precision; the right
angle can independently be nonzero. Multiplying the actual and intended
accepted blocks by the same native $`C^\dagger`$ preserves this error,
including when C is complex. This supplies no hardness result: even an
easy target can expose a bad word. It excludes this unamplified
shared-conjugator ansatz, not other jointly programmed words, completions,
or compilers.

Finally, outer cancellation alone does not give a global precision ledger.
Each remaining B still contains three source calls in the stated native
emission. Disjoint nodes may be batched by depth; that emission retains
a precision charge per depth. This is a cost of the displayed word, not
an additive lower bound under every possible circuit rewrite.

## 6. Changing source width without renewing its preparation

The grouped compiler uses unequal source widths. Factoring a source out of
each complete group therefore introduces a change of basis between groups,
even when its basis changes cancel at a fixed width. This section prices
that change and gives a different loader whose width transitions are cheap.
It resolves the source-boundary ledger for that loader, not the native cost
of the transformed group programs.

### The original loader has a costly bridge in the required direction

Retain $`U_m=R_{m-2}\cdots R_0`$ from Section 2, extending every
operator by identity to a common dirty pool. For $`M\geq m`$,

```math
U_MU_m^\dagger=R_{M-2}\cdots R_{m-1}
```

uses only the added rotations. But if a complete group has been written
$`G_m=U_mK_mU_m^\dagger`$, the bridge when width decreases from
M to m is $`U_m^\dagger U_M`$. These are different operators.
For $`m\geq2`$, already a one-bit decrease requires

```math
U_m^\dagger U_{m+1}=U_m^\dagger R_{m-1}U_m.
```

Use actual native loaders in the emitted word. In the one-T realization,
$`\widehat R_j=e^{-i\pi/8}R_j`$ and
$`\widehat U_m=e^{-i(m-1)\pi/8}U_m`$; hence the literal bridge
is $`e^{-i\pi/8}U_m^\dagger R_{m-1}U_m`$. Its scalar is retained,
and its displayed implementation has $`2m-1`$ T gates. This count is
optimal for this complete bridge, even allowing returned stabilizer and
arbitrary dirty helpers.

For the lower bound, write

```math
\begin{aligned}
U_m^\dagger(Y_{m-1}X_m)U_m
&=\sum_{k=0}^{m-1}a_kY_kZ_{k+1}\cdots Z_{m-1}X_m,\\
a_0^2=a_1^2&=2^{-(m-1)},\qquad
a_k^2=2^{-(m-k)}\ (k\geq1).
\end{aligned}
```

The normalized Pauli correlation of the bridge with probe $`Z_k`$ is
$`1-(1-1/\sqrt2)a_k^2`$. At k zero it is

```math
1-(1-1/\sqrt2)2^{-(m-1)}
=\frac{1+(2^{m-1}-1)\sqrt2}{(\sqrt2)^{2m-1}}.
```

Its numerator is not divisible by $`\sqrt2`$ in
$`\mathbb Z[\sqrt2]`$, so the returned-helper Pauli-transfer
argument of Section 2 gives the exact lower bound $`2m-1`$.

The obstruction also holds at the group's requested precision, without
requiring error smaller than the source's last coefficient. Suppose a
native circuit approximates this bridge with full-isometry error at most
$`2^{-L}`$, where $`L\geq6`$ and $`m\geq L-3`$. Choose
$`k=m-L+4`$, which lies between 1 and $`m-1`$. The target correlation
has gap $`16(1-1/\sqrt2)2^{-L}`$ below one. The actual correlation c
differs by at most $`2^{1-L}`$, so

```math
0\lt1-c\lt8\,2^{-L}.
```

For a circuit with t T gates, its Galois conjugate correlation lies in
$`[-1,1]`$ and the same algebraic-norm argument gives
$`1-c\geq2^{-t-1}`$. Therefore $`t\geq L-3`$. This prices one
specified bridge; it is not a compiler lower bound or a license to add
bridge costs after arbitrary circuit rewrites.

### A reverse-order star loader has the same coefficient grid

Keep the Majoranas $`\Gamma_j=Z_0\cdots Z_{j-1}X_j`$ and define
actual native rotations, for $`j\geq1`$, by

```math
q_j=e^{-i\pi/8}\exp(i\pi Y_0X_j/8),\qquad
F_j=\prod_{k=1}^{j}\mathrm{CNOT}_{k\to0},\qquad F_0=I,
```

```math
\begin{aligned}
p_j&=F_{j-1}q_jF_{j-1}^\dagger
=e^{-i\pi/8}\exp(i\pi Y_0Z_1\cdots Z_{j-1}X_j/8),\\
V_m&=p_1p_2\cdots p_{m-1}.
\end{aligned}
```

Products here are matrix products, with the rightmost factor acting first.
The phase in each q is literal: if
$`C_j=(S_0H_0)H_j`$, then

```math
q_j=C_j\mathrm{CNOT}_{0\to j}T_j^\dagger
\mathrm{CNOT}_{0\to j}C_j^\dagger.
```

This uses one T-dagger and eight Clifford gates. No controlled T gate or
initialized helper is hidden in the definition.

Conjugation successively splits the coefficient of $`\Gamma_0`$,
in the reverse index order. Thus

```math
\begin{aligned}
M_m^\star&=V_mX_0V_m^\dagger
=\sum_{j=0}^{m-1}\sqrt{w_j}\,\Gamma_j,\\
w_0&=2^{1-m},\qquad w_j=2^{j-m}\ (1\leq j\lt m).
\end{aligned}
```

These are the original geometric weights in reverse order. Reversing the
entire certified sign word, including the literal endpoint cases, preserves
the signed coefficient grid and its error bound. The programmed masks
remain $`P_f=\prod_jZ_j^{f_j}`$, implemented by the same addressed
dirty queries before the source basis is extracted. Coefficient tables
retain their row counts. Scalars in V cancel in every source and
controlled-source conjugation.

The long Pauli strings do not force a quadratic Clifford loader cost.
The CNOTs in the F words commute and square to identity, giving

```math
V_m=q_1\mathrm{CNOT}_{1\to0}q_2\mathrm{CNOT}_{2\to0}
\cdots\mathrm{CNOT}_{m-2\to0}q_{m-1}F_{m-2}^\dagger.
```

This literal word has $`m-1`$ T gates and $`10m-12`$ Clifford gates.
For a decrease $`M\gt m`$, the required bridge is now

```math
\begin{aligned}
V_m^\dagger V_M&=p_mp_{m+1}\cdots p_{M-1}\\
&=F_{m-1}q_m\mathrm{CNOT}_{m\to0}q_{m+1}\cdots
\mathrm{CNOT}_{M-2\to0}q_{M-1}F_{M-2}^\dagger.
\end{aligned}
```

It has $`M-m`$ T gates and $`8(M-m)+2M-4`$ Clifford gates in this
emission. An equal-width bridge is omitted. All these identities are on
the entire dirty space, including correlations with any other registers.

### Three groups and the general boundary ledger

Let $`G_g`$ denote the complete actual word for group g, with source
$`M_{m_g}^\star`$, and define $`K_g=V_{m_g}^\dagger G_gV_{m_g}`$.
For three groups executed in order, with $`m_1\geq m_2\geq m_3`$,

```math
G_3G_2G_1
=V_{m_3}K_3(V_{m_3}^\dagger V_{m_2})K_2
 (V_{m_2}^\dagger V_{m_1})K_1V_{m_1}^\dagger.
```

The same identity holds for any number of groups and for all occupied
flag and dirty inputs. The bridges use actual inverse words. There is no
projection, reset, or multiplication of accepted blocks in this equation.
For the retained monotone group widths, the total T count of the initial
loader, final loader, and bridges in the displayed emission is exactly

```math
(m_{\max}-1)+(m_{\min}-1)
 +\sum_g(m_g-m_{g+1})=2(m_{\max}-1)=O(L+n).
```

Their Clifford count is $`O(Rm_{\max})`$, where $`R\leq n`$ is
the number of groups. Since $`m_{\max}=L+O(n)`$ and $`L\geq6`$,
this is $`O((L+n)n)=O(NL)`$. The fixed deepest-layer tail remains
separate. These counts exclude the interior programs $`K_g`$.

Every basis change uses the existing dirty pool. At width m the source
acts on its first m wires; the other wires are available as arbitrary
dirty selectors under the existing per-group allocation. A bridge may
couple wires that change roles, and they may remain correlated. Exact
dirty-query contracts already permit such inputs. No wire is declared
clean by a change of role, and the whole-word error must still include
all final dirty return and reference correlations.

Within a group, logical maps and flag operations commute with V when their
complete action is identity on the dirty pool. If they borrow a source
wire, this commutation is justified only for their completed exact-return
primitive words, not for their individual physical gates. The actual
coarse word has that exact-return contract. This observation does not
make the coefficient queries commute with V: their masks become

```math
\widetilde P_f=V_m^\dagger P_fV_m.
```

Those transformed words must be implemented and charged. Assigning them
the old Pauli-query cost merely because the classical sign tables are
unchanged would omit the unresolved native work.

### A legal grouped coefficient still gives an expensive transformed mask

The issue is present for a table arising from the actual grouped residual,
not only for an arbitrary sign mask. Take a group of height $`s\geq1`$,
local dimension $`M=2^s`$, and padded term count
$`K=2^{\lceil\log_2(8s+12)\rceil}`$. Let its native coarse frame be
C equal to identity. Set every target local word to identity except the
root rotation $`R_y(\delta)`$, where

```math
\sin\delta=\frac{1}{9sK\sqrt M}.
```

The small positive choice obeys
$`\delta\leq2/(9sK\sqrt M)\lt1/(4sK\sqrt M)`$, so it meets
the grouped coarse-word allowance. The forward root-column entry at row
$`M/2`$ is $`\sin\delta`$. Its actual normalized scalar table entry is

```math
c_\ast=K\sqrt M\sin\delta=\frac1{9s}.
```

Let f be its certified source sign word after reversal. Its exact moment
$`c_f`$ satisfies $`|c_f-c_\ast|\leq(5/2)2^{-m}`$. For
$`D=V_m^\dagger P_fV_m`$, the normalized Pauli correlation with
$`X_0`$ is exactly $`c_f`$.

Suppose a native circuit approximates this D within full-isometry error
$`2^{-m}`$, including returned dirty and initialized stabilizer helpers.
Its correlation c then satisfies

```math
|c-c_\ast|\leq(9/2)2^{-m},\qquad
c\in(\sqrt2)^{-t}\mathbb Z[\sqrt2],\qquad c,c'\in[-1,1],
```

where $`c'`$ is its Galois conjugate. The number $`9sc-1`$ is nonzero:
a rational in this dyadic quadratic ring cannot have the odd denominator
of $`1/(9s)`$. Its nonzero algebraic norm therefore gives

```math
\begin{aligned}
2^{-t}&\leq|(9sc-1)(9sc'-1)|\\
&\leq\frac{81s}{2}(9s+1)2^{-m}
\leq405s^2 2^{-m}\lt512s^2 2^{-m}.
\end{aligned}
```

Consequently

```math
t\geq\max\{0,\ m-2\log_2s-9\}.
```

This bounds one separately implemented transformed mask at the stated
accuracy. The target witness itself is a simple one-rotation frame; it is
not a hard instance for all compilers. No sum of these mask bounds is an
unrestricted frame lower bound. A jointly synthesized interior program
could avoid exposing these masks as separate subroutines.

The reverse-order loader therefore supplies the proposed cheap source
boundaries and preserves the coefficient encoding, live width, and
Clifford budget. It does not yet supply cheap transformed group bodies.
The next native construction must address those bodies jointly or change
their interface; the generic endpoint frontier is unchanged.

## 7. A flag-correlated source boundary and its query cost

A correlated flag can carry a nontrivial source code even when the dirty
input is arbitrary. The following candidate has cheap changes of source
width, but its simplest query-and-renewal interface does not remove the
precision charge. This is a test of that interface, not a restriction on
all flag- or logical-correlated boundaries.

Use the **original chain source**, rather than the star-source eigenbasis
of Section 6. Let $`L_m=\widehat R_{m-2}\cdots\widehat R_0`$ be
its literal native loader, retaining the scalar phases of its actual
Pauli T words, and put $`M_m=L_mX_0L_m^\dagger`$. Embed every width
in the same maximum-width dirty register. Define

```math
\begin{aligned}
J_X|\psi\rangle
&=\frac{|0\rangle_f|\psi\rangle+|1\rangle_fX_0|\psi\rangle}{\sqrt2},\\
E_m&=(I_f\otimes L_m)J_X.
\end{aligned}
```

Preparing this isometry from a zero flag costs $`m-1`$ T gates for
the literal loader, plus Clifford operations; final decoding uses its
charged actual inverse. Only f is initialized. All dirty amplitudes and
reference correlations are retained. Its projector and encoded source action are

```math
\Pi_m=E_mE_m^\dagger=\frac{I+X_f\otimes M_m}{2},
\qquad
(I_f\otimes M_m)E_m=(X_f\otimes I)E_m.
```

Thus M can be replaced by a flag X **while the input is in this code**.
That qualification is essential.

### Width changes are cheap on this correlated boundary

For $`k\lt m`$, the actual loader words give

```math
L_kL_m^\dagger
=\widehat R_{k-1}^\dagger\cdots\widehat R_{m-2}^\dagger,
\qquad
(I_f\otimes L_kL_m^\dagger)E_m=E_k.
```

The transition costs $`m-k`$ T gates and $`O(m-k)`$ Clifford gates;
its actual inverse handles growth. The equalities retain every literal
loader phase. This is the physical code transition $`L_kL_m^\dagger`$,
not the oppositely oriented eigenbasis bridge considered in Section 6.

After shrinking to k, the projector is identity on the released tail.
In fact any operator Z supported on that tail obeys
$`(I_f\otimes Z)E_k=E_kZ`$. Those wires are again arbitrary dirty
inputs, not initialized bits: they may remain entangled with the retained
core or an external reference. The identity permits a tail to be reassigned
as returned dirty workspace, but does not prove compatibility with a
particular group program.

### A valid scalar mask changes the code

For a Pauli mask P, let $`N=PM_mP^\dagger`$. Applying P to the dirty
register changes the code projector to
$`(I+X_f\otimes N)/2`$. If $`\{M_m,N\}=2cI`$, direct multiplication
gives the full-input leakage identity

```math
\left[(I-\Pi_m)(I_f\otimes P)E_m\right]^\dagger
\left[(I-\Pi_m)(I_f\otimes P)E_m\right]
=\frac{1-c}{2}I.
```

For every $`m\ge5`$, the allowed mask $`P=Z_1Z_2Z_3`$ programs
$`c=1/8`$: it flips source weights $`1/4,1/8,1/16`$, with the final
mask bit zero. This row is realized by a one-level grouped dictionary:
take native C equal to identity and a real root angle satisfying
$`\sin\theta=1/(8K\sqrt2)`$. Its coarse error is at most
$`2\sin\theta=\delta_c`$, while the special forward coefficient is
$`K\sqrt2\sin\theta=1/8`$. The mask's leakage norm is exactly
$`\sqrt7/4`$, independently of precision. In particular,

```math
\left\|\bigl(I_f\otimes M_m-X_f\otimes I\bigr)
       (I_f\otimes P)E_m\right\|=\sqrt7/2.
```

The previously cheap flag substitution therefore fails after the query.
One may instead carry the changed code, but subsequent operations must
then respect that different boundary. Cheap width transport alone does
not supply this program-dependent closure.

### Full syndrome renewal recharges precision

A natural proposed repair uses the second clean flag b to coherently
extract the syndrome of $`S_m=X_f\otimes M_m`$. Write
$`\Pi_\pm=(I\pm S_m)/2`$. The complete-input extraction target is

```math
\mathcal E|\psi\rangle
=|0\rangle_b\Pi_+|\psi\rangle
 +|1\rangle_b\Pi_-|\psi\rangle,
```

where psi is arbitrary on the carried flag and dirty register. An actual
extractor V with $`VJ_b=\mathcal E`$ obeys

```math
V^\dagger Z_bVJ_b=J_bS_m.
```

Stripping the free $`X_f`$ thus implements $`M_m`$ and returns the
extraction flag. The source's exact returned-helper minimum
$`T(M_m)=2m-4`$ implies

```math
T(V)\ge m-2.
```

A concrete extractor is $`H_bC_b(S_m)H_b`$; the controlled M inside
it has its explicitly charged native source implementation. The inequality
concerns the full-syndrome interface, not every possible way to carry or
change the code.

The obstruction also applies at the compilation accuracy. If
$`\|VJ_b-\mathcal E\|\le\epsilon`$, then, using the **actual** inverse,

```math
\|V^\dagger Z_bVJ_b-J_bS_m\|\le2\epsilon.
```

This follows by $`Z_b\mathcal E=\mathcal E S_m`$ and two applications
of the extraction error. For $`L\ge6`$, $`m\ge L+2`$, and
$`\epsilon\le2^{-L}`$, it gives

```math
T(V)\ge\frac{L-4}{2}.
```

For completeness, put $`j=L-3`$ and probe the implemented M with
$`Z_j`$. The exact source correlation is $`1-2^{-j}`$. The doubled
extractor error changes that correlation by at most $`4\epsilon`$,
which is at most $`2^{-j-1}`$. Hence the native correlation q satisfies

```math
2^{-j-1}\le1-q\le3\,2^{-j-1}.
```

If the resulting source circuit has t T gates, its Pauli-transfer
coefficient lies in $`(\sqrt2)^{-t}\mathbb Z[\sqrt2]`$. Its Galois
conjugate is another unitary-channel correlation in $`[-1,1]`$, including
initialized and arbitrary dirty helpers. The nonzero algebraic norm gives
$`1-q\ge2^{-t-1}`$, so $`t\ge j-1=L-4`$. The circuit uses V and
its actual inverse, proving the bound above.

This candidate consequently resolves source-width transport but not
source reuse through a scalar program. A transition that explicitly
renews this full syndrome pays another length-L cost. A more general
boundary may avoid renewal, correlate the source with logical data, or
defer individual group action; none is excluded by these identities.

## 8. Small products: compress before synthesis

Small examples can identify an interface before a general construction is
available. Begin with a positive case, then change its address or logical
support. These checks distinguish an actual product identity from a
cancellation that exists only for separately accepted blocks.

### A positive small example: multiply first within one two-mode space

Write a determinant-one one-qubit unitary as

```math
U(q)=wI-i(xX+yY+zZ),\qquad q=(w,\mathbf v),\qquad
w^2+\|\mathbf v\|^2=1.
```

The complete product stays in these same four real coordinates:

```math
q_2\star q_1=
(w_2w_1-\mathbf v_2\cdot\mathbf v_1,
 w_2\mathbf v_1+w_1\mathbf v_2+\mathbf v_2\times\mathbf v_1).
```

For example, let $`R_j(\theta)=e^{-i\theta\sigma_j}`$ and abbreviate
$`c_\alpha=\cos\alpha`$, $`s_\alpha=\sin\alpha`$. Two noncommuting
rotations give

```math
R_x(\alpha)R_y(\beta)
=c_\alpha c_\beta I-i(s_\alpha c_\beta X
 +c_\alpha s_\beta Y+s_\alpha s_\beta Z).
```

Left-multiplying by $`R_z(\gamma)`$ gives the quaternion

```math
\begin{pmatrix}
c_\gamma c_\alpha c_\beta-s_\gamma s_\alpha s_\beta\\
c_\gamma s_\alpha c_\beta-s_\gamma c_\alpha s_\beta\\
c_\gamma c_\alpha s_\beta+s_\gamma s_\alpha c_\beta\\
c_\gamma s_\alpha s_\beta+s_\gamma c_\alpha c_\beta
\end{pmatrix}.
```

Any further factor updates these four numbers; it does not enlarge the
logical support. These are exact matrix identities with the literal phase
retained. This is why noncommutation alone need not cause repeated precision
cost.

If every factor has the same unchanged k-bit address, classically compute
the complete product at each of its $`S=2^k`$ rows. The retained certified
Euler procedure then approximates each entire product by three rotations.
The [one-target compiler](../../docs/ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
already supplies $`T=O(S+L)`$ and
$`G=O(SL)`$ with its stated sufficient allocation
$`a=1,b\ge L+k+9`$. Its complete-isometry estimate includes all logical
inputs, dirty return, and references. This is a corollary of that compiler,
not a new workspace theorem. The source core is reused through a constant
number of calls, independent of the original word length; this does not
prove a single literal preparation and unpreparation.

Certified preprocessing must approximate the whole product to the required
matrix tolerance before allocating the fixed number of synthesis errors.
Its work and working precision may grow with the original word length.
No quantum-controlled unknown angles or free coherent evaluator are assumed.

### A changed address invalidates rowwise multiplication

Let $`A_j=\sum_x|x\rangle\langle x|\otimes U_{j,x}`$ and let C
change the address. The complete action is

```math
A_2(C\otimes I)A_1
=\sum_{x,y}C_{yx}|y\rangle\langle x|\otimes U_{2,y}U_{1,x}.
```

Multiplying only $`U_{2,x}U_{1,x}`$ at each old row loses the terms
with different x and y. A single address Hadamard already exposes this
failure. Classical product compression remains valid on the full joint
logical space, but its program is then no longer the original unchanged address table. This is the next condition to test when extending the
positive example.

### The first change of logical support

On three modes, two overlapping pair rotations instead give

```math
R_{bc}(\beta)R_{ab}(\alpha)=
\begin{pmatrix}
c_\alpha&-s_\alpha&0\\
c_\beta s_\alpha&c_\beta c_\alpha&-s_\beta\\
s_\beta s_\alpha&s_\beta c_\alpha&c_\beta
\end{pmatrix}.
```

The mixed path amplitude $`s_\beta s_\alpha`$ is already present. With
$`A_{ij}=|i\rangle\langle j|-|j\rangle\langle i|`$,
$`[A_{ab},A_{bc}]=A_{ac}`$: the two edges generate the three-dimensional
Lie algebra of SO(3). A fixed three-mode space still has constant-size
Euler coordinates, but adding a connected fourth mode gives SO(4), of
dimension six. Repeated path commutators on the connected N-mode Hopf tree
generate every $`A_{ij}`$, giving the Lie algebra of SO(N).

This last statement concerns arbitrary repeated products. The prescribed
once-per-edge Hopf frame still has only $`N-1`$ parameters and its retained
linear generator description. It is not a gate lower bound. The lesson for
the next small example is precise: a useful compression must retain a
bounded-cost program as the active modes grow, not only as more rotations
are multiplied on one fixed target.

The [small tree Cayley construction](ENDPOINT_TREE_TRANSPORT.md#6-small-products-suggest-a-cayley-representation)
continues this product-first experiment on four and eight logical modes.
It preserves the complete residual in a coupled recursive representation;
its generic native source-and-program implementation remains unproved.
The [four-mode benchmark](ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark)
now prices a real four-mode row with two retained one-target programs and
a borrowed logical spectator. The eight-mode extension retains a
controlled root coupling; independent factorization does not remove it.
The [joint changing-target word](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
works in the original tree basis and shares one fixed scalar conjugator.
Its eight-/sixteen-mode source-appearance reduction survives complete
signal and dirty inputs, with two final parity boundaries. Hoisting the
loader and simplifying the fixed mask exposes the same leading precision
cost in the comparison word: the remaining programmable masks still
carry a precision charge per layer.

## 9. Dirty programs and phase-gradient sources

The following identities distinguish complete dirty-input operations from
prepared-source promises. They do not exclude joint native synthesis.
The broader literature comparison is in
[Related Work, Section 19](../RELATED_WORK.md#19-unary-phase-source-reuse-and-linear-t-depth-3-october-2026).

### Whole-word XOR cancellation and a cyclic alternative

Suppose a unitary family indexed by bit strings satisfies
$`U_{z\oplus s}U_z^\dagger=V_s`$ independently of arbitrary dirty z.
Taking $`z=0`$ and substituting back gives

```math
U_s=V_sU_0,\qquad V_{z\oplus s}=V_sV_z.
```

Thus the $`V_s`$ must be commuting involutions; conversely that condition
suffices. This is a literal-phase identity. In particular, if the two echo
outputs $`U_sU_0^\dagger`$ and $`U_0U_s^\dagger`$ both approximate V
within $`\epsilon`$,
then $`\|V-V^\dagger\|\le2\epsilon`$. For $`V=R_y(\theta)`$ this
requires $`\epsilon\ge|\sin\theta|`$. An involution alphabet alone is
insufficient: $`U_{z_1,z_2}=H^{z_2}X^{z_1}`$ and $`s=(1,0)`$ produce
X or Z depending on the dirty second bit. The retained
[reflection interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
instead completes each symbol's cancellation before starting the next.

A cyclic character does allow a complete dirty phase echo. For
$`M=2^m`$, define $`P_M|z\rangle=e^{2\pi iz/M}|z\rangle`$ and
$`A_f|x,z\rangle=|x,z+f(x)\bmod M\rangle`$. Then

```math
A_f^\dagger(I\otimes P_M)A_f(I\otimes P_M^\dagger)
=\sum_x e^{2\pi if(x)/M}|x\rangle\langle x|\otimes I.
```

Evaluating on $`|x,z\rangle`$ proves the identity, hence also dirty return
under arbitrary reference entanglement. Native costs are
$`2T(A_f)+2T(P_M)`$; actual adjoints and complete-isometry primitive errors
$`\epsilon_A,\epsilon_P`$ give total error at most
$`2\epsilon_A+2\epsilon_P`$. Both operations still need charged circuits.
XOR lookup cannot replace modular addition: for $`M=4,f=1`$ its phase is
i on even z and $`-i`$ on odd z. A cheap dirty phase-gradient unitary,
a cheap dirty addition table, and complete-frame amortization remain
separate missing constructions.

### A prepared tuple is not an arbitrary dirty register

Zhang, Tan, Kothari, Gosset, and Gidney,
[*Quantum circuit compilation with constant overhead*](https://arxiv.org/html/2609.39092v1),
Section 5, give unitary $`O(m)`$-gate phase-gradient preparation at error
$`2^{-m}`$. Their construction prepares a coefficient register from
$`|0^B\rangle`$ and a tuple register $`|+^M\rangle`$, with
$`B=O(m)`$ and $`M=128m`$ (Eqs. 46–54). These are initialized registers.

There is a useful operator extension: remove the initial Hadamards on the
data while retaining those prepared registers. Equations 60–64 then give
a block approximating $`\gamma e^{-iH}`$ for
$`H=\pi\sum_{j=1}^m2^{-j}Z_j`$, uniformly on every data input. Indeed,
the selected Paulis have norm one, so the preparation, tuple-tail, and
Taylor-tail bounds apply to an arbitrary data vector and reference.
Constant-round oblivious amplification and correction of the known phase
of $`\gamma`$ preserve linear T-count with that initialized work. This
does not fit two clean qubits; at $`m=N`$ the tuple alone exceeds the
total physical width.

The obstruction to a direct dirty substitution is explicit. Let F project
onto tuple strings with at least one 1. Although
$`\|(I-F)|+^M\rangle\|=2^{-M/2}`$, one has $`\|I-F\|=1`$.
Even valid tuples fail to give the required average. At Taylor order one,
the orthogonal valid states

```math
|a\rangle=|1\rangle|+^{M-1}\rangle,\qquad
|b\rangle=|01\rangle|+^{M-2}\rangle
```

select $`Z_1`$ and $`Z_2`$, respectively. On data $`|10\cdots0\rangle`$,
SELECT sends $`(|a\rangle+|b\rangle)/\sqrt2`$ to the orthogonal dirty
state $`(-|a\rangle+|b\rangle)/\sqrt2`$. Its distance from any common
data-only scalar phase, on these two dirty inputs, is at least $`\sqrt2`$.
The matrix element in the prepared tuple state is an average; the
tuple-controlled unitary is not that average tensored with identity.

### A Majorana first moment does not supply independent powers

Let $`p_j\ge0`$ sum to one, let $`\Gamma_j`$ be anticommuting Hermitian
involutions, and let $`P_j`$ be commuting Hermitian involutions on disjoint
data wires. Set

```math
M=\sum_j\sqrt{p_j}\Gamma_j\otimes I,\qquad
N=\sum_j\sqrt{p_j}\Gamma_j\otimes P_j,\qquad
H_0=\sum_jp_jP_j.
```

Pairing the cross terms proves
$`M^2=N^2=I`$ and $`(MN+NM)/2=I\otimes H_0`$ on the entire source
space. For $`p_j=2^{-j}`$ ($`1\le j\le m`$), a capped identity tail
$`p_{m+1}=2^{-m}`$ with $`P_{m+1}=I`$, and $`P_j=Z_j`$ for
$`1\le j\le m`$, the geometric Majorana loader
and a data-controlled Pauli mask give an $`O(m)`$-T first-moment block.
The mask is Clifford because its source-Pauli exponents are linear
functions of the data bits.

Suppressing the source identity in $`H_0`$, the actual one-flag word
$`H_f\,\mathrm{diag}(MN,NM)\,H_f`$ is

```math
S=I_f\otimes H_0+X_f\otimes D,\qquad
D=(MN-NM)/2,\qquad [D,H_0]=0,\quad D^2=H_0^2-I.
```

Reusing its occupied flag therefore gives

```math
\langle0|S^2|0\rangle=2H_0^2-I,\qquad
\langle0|S^k|0\rangle=T_k(H_0),
```

where $`T_k`$ is the Chebyshev polynomial, rather than the independent
tuple moment $`H_0^k`$. For $`m=2`$,
$`H_0=Z_1/2+Z_2/4+I/4`$ vanishes on $`|10\rangle`$; the second
compressed moment is $`-I`$ there instead of zero. Signal processing can
use these Chebyshev moments, but a degree-K word then charges K source
calls in the direct implementation. At the Taylor degree
$`K=\Theta(m/\log(m+2))`$, this is $`O(mK)`$ T gates before phase
synthesis, not a single $`O(m)`$ charge. No source-fusion identity is
supplied by the first-moment equation.

There are two further interface conditions. The paired-Majorana realization
uses at least $`\lceil(m+1)/2\rceil`$ source qubits in addition to the
m-bit data bank; at $`m=N`$ this misses the endpoint width for large n.
Identifying source and data would invalidate their assumed commutation.
More generally, for nonzero real weights and logical Hermitian coefficients,

```math
S=\sum_j a_j A_j\otimes\Gamma_j
\quad\Longrightarrow\quad
S^2=\sum_j a_j^2A_j^2\otimes I
+\sum_{j\lt k}a_ja_k[A_j,A_k]\otimes\Gamma_j\Gamma_k.
```

For the independent paired Majoranas, the identity and distinct grade-two
monomials are linearly independent. Hence this linear Hermitian source is
unitary exactly when $`\sum_j a_j^2A_j^2=I`$ and all logical coefficients
commute. Two overlapping tree-edge swaps A,B obey
$`\|[A,B]\|=\sqrt3`$: their products are opposite three-cycles on their
three-dimensional support. Thus
$`S=(A\otimes\Gamma_0+B\otimes\Gamma_1)/\sqrt2`$ has
$`\|S^2-I\|=\sqrt3/2`$. Higher Majorana grades are outside this linear
ansatz; their full-word structure is treated below. None of these
restrictions is an unrestricted compiler lower bound.

### Phase halving requires its promised catalyst

Appendix B, Proposition 16 of the same paper assumes the catalyst
$`|\chi_+\rangle=R(\theta)|+\rangle`$, where
$`R(\theta)=\mathrm{diag}(1,e^{i\theta})`$. The measurement-free word
with both Toffolis retained has, on an arbitrary catalyst c and a zero
helper, the complete effective action

```math
\begin{aligned}
U_\theta={}&I_c\otimes|00\rangle\langle00|\\
&+R_c(2\theta)X_c\otimes
 (|01\rangle\langle01|+|10\rangle\langle10|)\\
&+e^{2i\theta}I_c\otimes|11\rangle\langle11|.
\end{aligned}
```

The helper returns to zero. Define

```math
Q_\theta=e^{-i\theta}R(2\theta)X
=\cos\theta X+\sin\theta Y.
```

The word is the desired
$`R(\theta)\otimes R(\theta)`$ followed by $`Q_\theta`$ on c,
controlled on odd data parity. The promised catalyst is its +1
eigenvector. Its orthogonal −1 eigenvector produces an extra
$`Z\otimes Z`$ on the data. In particular,

```math
|0\rangle_c|01\rangle|0\rangle_a
\longmapsto e^{2i\theta}|1\rangle_c|01\rangle|0\rangle_a,
```

at distance $`\sqrt2`$ from the requested dirty-return output
$`e^{i\theta}|0\rangle_c|01\rangle|0\rangle_a`$. Correcting the
controlled $`Q_\theta`$ requires a further angle-dependent operation.
Replacing the paper's measurement-based four-T uncomputation by exact
Toffolis changes a constant cost, not this catalyst requirement.

A prepared approximate catalyst remains useful: if a full ideal sequence
returns $`|\chi_+\rangle`$ and implements W, replacing that catalyst by
a normalized state within $`\delta`$ changes the complete action,
relative to returning the actual supplied state, by at most $`2\delta`$.
Compare the initial and final ideal catalyst through the full unitary.
This bound need not accumulate per catalytic use, but it presupposes a
legal preparation and does not apply to an arbitrary dirty qubit.

## 10. Higher-grade Spin identities on arbitrary dirty work

Put $`q=N/2`$ and use the N paired Jordan--Wigner Majoranas
$`\Gamma_{2r}=Z_{\lt r}X_r`$ and $`\Gamma_{2r+1}=Z_{\lt r}Y_r`$ on q
dirty qubits. A coordinate-plane rotation has a half-angle Spin lift
$`U_e=\exp(\pm\theta_e\Gamma_u\Gamma_v/2)`$, with orientation fixed
so its ordered tree product U satisfies

```math
U\Gamma_jU^\dagger=\sum_i W_{ij}\Gamma_i.
```

The following identities concern the full dirty Hilbert space; they use
no selected spinor state. The displayed lift still contains $`N-1`$
independently specified rotations, so its native synthesis cost must be
accounted for separately.

### Covariance and an exact native reflection

Define the isometry from spinor space to logical-index times spinor space,

```math
A|\psi\rangle=\frac1{\sqrt N}\sum_j|j\rangle\Gamma_j|\psi\rangle.
\qquad A^\dagger A=I,\qquad (W\otimes U)A=AU.
```

The first identity uses orthogonality of the index labels; covariance
follows by substituting the conjugation equation and using
$`\sum_jW_{kj}W_{ij}=\delta_{ki}`$. Let
$`S=\sum_j|j\rangle\langle j|\otimes\Gamma_j`$ and
$`T=S(H^{\otimes n}\otimes I)`$. Then

```math
R=2AA^\dagger-I
=T(2|0^n\rangle\langle0^n|\otimes I-I)T^\dagger
```

is a native reflection on **every** logical and dirty input. This formula
does not require the logical register to start in zero.

Two exact [dirty-word queries](../../docs/ONE_CLEAN_COMPILER.md) implement
the X and Z masks of S, using at most n returned dirty selectors.
The S gate on the index parity bit supplies the literal i in $`Y=iXZ`$.
The logical zero reflection uses the exact borrowed-MCX construction
with $`O(n^2)`$ T gates, borrowing already returned work. Its scalar
minus sign, when needed, is the exact Clifford word XZXZ. Thus R has
$`O(N)`$ T-count, $`O(Nq)`$ elementary Clifford count, and uses
$`q+n`$ dirty qubits with no added clean qubit. All selectors return
exactly. The image of A occupies only a $`1/N`$ fraction of the
logical-spinor space, so covariance on this image is not a full-input
implementation of W.

### Low-rank reflection products require N factors

Let d denote the full dirty input dimension. For any projectors
$`P_1,\ldots,P_r`$ of rank at most d on the ambient physical space,

```math
\mathrm{rank}\left[\prod_{j=r}^{1}(I-2P_j)-I\right]\le rd.
```

Indeed $`AB-I=A(B-I)+(A-I)`$ makes the rank subadditive, and each
factor differs from identity by rank at most d. The statement includes
arbitrary conjugated rank-d reflections; any conjugating unitary may
depend on the target. For the preceding Majorana isometry, $`AA^\dagger`$
has rank d, including identity factors on additional dirty work. The same
bound applies if its opposite reflection convention
$`2AA^\dagger-I`$ is used an even number of times.

Choose the admitted Hopf frame with every upper angle zero and all
final-depth angles $`\pi/4`$. It is a direct sum of $`N/2`$ real
rotations and satisfies $`\sigma_{\min}(W-I)=2\sin(\pi/8)`$.
Let J embed its Nd-dimensional logical/dirty input into any clean-work
extension. If the reflection product V satisfies

```math
\|VJ-J(W\otimes I_d)\|\lt2\sin(\pi/8),
```

then $`(V-I)J`$ is injective. A nonzero kernel vector would give
approximation error at least that smallest singular value. Therefore
$`Nd\le\mathrm{rank}(V-I)\le rd`$, so $`r\ge N`$.
This is a restriction on products of these low-rank reflections; it
does not apply to arbitrary interleaved operations or projectors whose
rank grows beyond d.

### The reflection-only interface loses sensitivity

Direct expansion gives

```math
A^\dagger(W\otimes I)A
=\frac{\mathrm{tr}(W)I+
 \sum_{j\lt k}(W_{jk}-W_{kj})\Gamma_j\Gamma_k}{N}.
```

For a single plane rotation of angle theta all overlap singular values
equal

```math
s_\theta=\sqrt{1-\frac{4(N-2)(1-\cos\theta)}{N^2}}.
```

Indeed the overlap is
$`[(N-2+2\cos\theta)I\pm2\sin\theta\Gamma_a\Gamma_b]/N`$;
squaring its singular values gives the formula. Covariance implies

```math
R_U=(I\otimes U)R(I\otimes U^\dagger)
=(W^{\mathsf T}\otimes I)R(W\otimes I),\qquad
\|R_U-R\|=\frac{4\sqrt{(N-2)(1-\cos\theta)}}{N}.
```

The norm follows from the principal angles of the two equal-rank
projector ranges: reflection distance is $`2\sqrt{1-s_\theta^2}`$.
At $`\theta=\pi`$ it is $`\Theta(N^{-1/2})`$, whereas
$`\|W-I\|=2`$. A unitary hybrid therefore requires
$`\Omega(\sqrt N)`$ queries to distinguish the two target isometries
with fixed error below one if the only angle-dependent operation is
$`R_U`$. Controlled queries, inverses and fixed auxiliary encodings obey
the same hybrid bound. Direct U queries and angle-dependent interlayers
are outside this reflection-only restriction.

### Matched dirty frames force commuting logical labels

Suppose Clifford programming sends a family of anticommuting dirty
Majoranas to the matched form
$`Q_j=P_j\otimes\Gamma_j`$, with Hermitian logical Paulis
$`P_j`$ and the same dirty $`\Gamma_j`$. For distinct j and k,

```math
0=\{Q_j,Q_k\}
=[P_j,P_k]\otimes\Gamma_j\Gamma_k.
```

The dirty product is invertible, so all logical labels commute. For
$`A_0=\sum_j a_j I\otimes\Gamma_j`$ and
$`B_0=\sum_j b_jP_j\otimes\Gamma_j`$, the scalar-overlap identity is

```math
\frac{A_0B_0+B_0A_0}{2}
=\sum_j a_jb_jP_j\otimes I.
```

It therefore lies in a commuting logical Pauli algebra. A permutation of
the matched dirty frame gives the same result after relabelling. The
restriction is the matched form itself; different dirty frames or more
general extraction words are outside this argument.

### A fixed Clifford organizes the complete Choi expansion

For each tree edge let $`P_e=i\Gamma_u\Gamma_v`$ and absorb fixed
orientation signs into its angle. The binary Pauli labels of these
$`N-1`$ bivectors are independent. In fact, the N Majorana labels have
symplectic Gram matrix $`I+\mathbf1\mathbf1^{\mathsf T}`$ over
$`\mathbb F_2`$, whose square is I because N is even. They form a
basis. The tree incidence columns are independent by leaf elimination;
multiplying them by this basis gives the independent edge labels.

Expand the actual ordered word

```math
U=\prod_e\bigl(\cos(\theta_e/2)I+i\sin(\theta_e/2)P_e\bigr).
```

Every subset of edges has a distinct Pauli label. In canonical
$`X^xZ^z`$ order its coefficient phase is a linear power of i times
a quadratic sign in the subset bits: multiplication contributes
$`(-1)^{z\cdot x'}`$. Thus a fixed, angle-independent Clifford E on
the $`2q=N`$ Choi qubits satisfies

```math
|U\rangle\!\rangle
=E\left[\bigotimes_e
  \bigl(\cos(\theta_e/2)|0\rangle+\sin(\theta_e/2)|1\rangle\bigr)
  \otimes|0\rangle\right],
```

where $`|U\rangle\!\rangle`$ is normalized vectorization. Construct E
by the S/CZ phase function, an invertible CNOT map extending the edge
labels to a basis, and the Clifford Bell-basis encoding. This accounts
for all Majorana grades and their literal phases.

This is a Choi-state identity, not a physical pre/post Clifford
factorization of U: E may mix the two Choi halves. Using it as a compiler
requires both a charged synthesis of the independent fine-angle factors
and a full-input extraction with the allowed clean/dirty return.
Neither prepared Choi qubits nor measurement-based gate teleportation is
part of the endpoint contract. The
[bounded-diagonal route](BOUNDED_DIAGONAL_FACTORIZATION.md) and
[tree Cayley reduction](TREE_CAYLEY_REDUCTION.md) provide alternative
global representations without assuming those resources.

### Actual source words share a Spin algebra before addressing

The paired-source primitive admits a different, exact joint representation
that retains its rejected action. Here use source precision q and
$`m=q+1`$ core qubits, and write $`D_f=MF`$, $`D_g=MG`$ as in the
[one-clean primitive](../../docs/ONE_CLEAN_COMPILER.md#3-conjugating-scalar-blocks-produces-a-rotation).
The source vectors F and G are perpendicular to M. For a Hermitian Pauli
K on separate arbitrary dirty work, define

```math
C_K=|0\rangle\langle0|_r\otimes I
    +|1\rangle\langle1|_r\otimes K.
```

Dressing both routing axes $`X_rZ_t`$ and $`X_rX_t`$ by K conjugates
the complete primitive Q by $`C_K`$. Since $`C_K`$ commutes with its
four signal reflections, the actual amplified word obeys

```math
Q_K=C_KQC_K,\qquad V_K=C_KVC_K.
```

It also commutes with the ideal logical operation tensored with work
identities. Thus the full-operator error, approximate dirty/reference
return, and exact inactive-predicate identity are unchanged. This does
not require intermediate signal or core return.

Choose n anticommuting Hermitian Pauli factors $`K_j`$ on
$`\lceil n/2\rceil`$ separate dirty qubits. If
$`F_1,\ldots,F_{2m-1}`$ is an orthonormal Majorana-vector basis
perpendicular to M, then these vectors together with

```math
\Lambda_{j,\alpha}=X_r\sigma_{j,\alpha}K_jM,
\qquad \alpha\in\{x,y,z\},
```

are mutually anticommuting Hermitian involutions. Pauli anticommutation
supplies the sign within each target, the K factors supply it between
targets, and M supplies it against every F. An unconditioned, address-independent
scalar route is a scalar plus a simple bivector in this one Clifford
algebra, hence a Spin rotation. Signal reflection conjugation preserves
the bivectors, and its four occurrences cancel in the five-call word.
Also $`\Lambda_{j,x}\Lambda_{j,y}=iZ_j`$, with cyclic analogues.
Consequently changing logical targets or rotation axes preserves this
common Spin representation for an unconditioned, address-independent
stream of the actual regular coins. Logical predicates are outside this
closure statement, although the exact dirty-gauge identity still applies
to them. This is a product representation; synthesizing
the resulting orthogonal matrix still needs a native cost certificate.

### Unequal addresses generate a term outside that Spin algebra

Let x be an address bit and t a distinct target. For the regular coin
$`Q(u)=R_y(\pi/2)R_x(\arctan u)R_z(\arctan u)`$, put

```math
B=|0\rangle\langle0|_x\otimes Q(u)
 +|1\rangle\langle1|_x\otimes Q(v),
\qquad u,v\in(-\sqrt3,\sqrt3).
```

Direct multiplication gives

```math
Q(u)ZQ(u)^\dagger
=-\frac{2u}{1+u^2}Y+\frac{u^2-1}{1+u^2}Z
=:\boldsymbol\nu(u)\cdot\boldsymbol\sigma.
```

The conjugate of $`\Lambda_{t,z}`$ has a component outside the
Clifford-vector span equal to

```math
\frac12 Z_xX_rK_tM
 [\boldsymbol\nu(u)-\boldsymbol\nu(v)]\cdot\boldsymbol\sigma_t.
```

It is orthogonal to that span in normalized Hilbert--Schmidt inner
product: it has nonidentity Pauli support on two logical targets, whereas
each original vector has support on at most one. Its norm in that inner
product is

```math
d(u,v)=\frac{|u-v|}{\sqrt{(1+u^2)(1+v^2)}}.
```

Any G in this particular Spin group maps vectors back to their span.
The estimate $`\|G-B\|\le\epsilon`$ would imply conjugation error
at most $`2\epsilon`$ in operator norm, hence also in normalized
Hilbert--Schmidt norm. Therefore

```math
\epsilon\ge d(u,v)/2,
\qquad (u,v)=(1/4,5/4)\ \Longrightarrow\
\epsilon\ge8/\sqrt{697}.
```

This restricts full-operator approximation inside the specified fusion
group. It does not restrict all native circuits, target-dependent
algebras, complete tree products, or other completions required only on
an initialized signal sector.

The K dressing uses Clifford gates only. At common q, the same stream of
$`2n`$ real banks still has the source charge
$`(80n+4)q+O(n)`$ from Section 2, plus queries and predicates, and
uses $`q+1+\lceil n/2\rceil`$ dirty qubits besides the supplied signal.
Thus its displayed cost remains $`O(N+nq+n^4)`$ T gates. The exact
common algebra neither discards rejected outputs nor removes this charge.

## 11. A single Clifford cannot carry the full frame through an encoding

Consider a physical word $`B_\theta^\dagger C_\theta B_\theta`$, where
$`C_\theta`$ is Clifford and the encoder $`B_\theta`$ may depend on every
angle. Allow arbitrary encoder cost and structure. On
$`M=n+b+a`$ physical qubits, with $`a\le2`$ initialized qubits, the required
complete-input action would give an isometry
$`E_\theta=B_\theta J_a`$ satisfying

```math
\|C_\theta E_\theta-E_\theta(W_\theta\otimes I_b)\|
\le\eta,
\qquad \eta=2^{-N},\qquad N=2^n.
```

**Proposition.** If such an isometry and literal physical Clifford exist
for every complete real Hopf frame, then

```math
\boxed{3M\ge\frac N2(N-n-11)-7.}
```

Consequently this architecture requires
$`M\ge N^2/6-O(N\log N)`$ and cannot fit
$`M\le N+2n+9`$ for $`n\ge5`$. The encoding may vary discontinuously with
the target. Neither exact invariance nor exact dirty return is assumed.

### Literal Clifford spectra have few distinct types

Let $`G_M`$ be the finite group generated by the literal H, S, and CNOT
matrices on $`M`$ qubits, and let $`k(G)`$ count conjugacy classes of a
finite group. Its symplectic quotient is $`\mathrm{Sp}(2M,2)`$; the
kernel is the Pauli group with eighth-root scalar phases and has size
$`8\cdot4^M`$. In particular the scalar phase is not discarded.

For any finite-group quotient $`G\to Q`$ with kernel $`K`$,
$`k(G)\le|K|k(Q)`$. Choose one representative from each conjugacy class
of $`Q`$. Every element of $`G`$ can be conjugated so that its image is
one of those representatives, and there are only $`|K|`$ elements over
each representative. Since conjugate unitaries have the same spectrum,
the number of distinct literal spectra of $`G_M`$ is at most

```math
8\cdot4^M k(\mathrm{Sp}(2M,2))
\le121.6\,2^{3M}\lt2^{3M+7}.
```

The imported inequality is
[Fulman–Guralnick, Theorem 3.13(2)](https://arxiv.org/pdf/0902.2238):
for even $`q`$, $`k(\mathrm{Sp}(2m,q))\le15.2q^m`$. Their theorem
applies to every rank; no limit with fixed rank or large field is used.

### Approximate eigenvectors force actual spectral multiplicity

Let $`E:\mathbb C^D\to\mathbb C^K`$ be an isometry and let $`C`$ be
unitary. If $`|\lambda|=1`$ and
$`\|(C-\lambda I)E\|\le\eta`$, the spectral projection of $`C`$ onto
$`|z-\lambda|\lt r`$ has rank at least $`D`$ for every $`r\gt\eta`$.
Otherwise its composition with $`E`$ has a nonzero kernel. For a unit
vector $`v`$ in that kernel, spectral calculus gives
$`\|(C-\lambda I)Ev\|\ge r`$, a contradiction.
This argument concerns an approximate eigenspace and does not replace it
by an exactly invariant one.

### Final-depth frames supply too many spectra

Set every earlier angle to zero. The complete operator is then

```math
W=\bigoplus_{j=0}^{N/2-1}R_y(\theta_j),
\qquad \mathrm{spec}(R_y(\theta_j))
=\{e^{i\theta_j},e^{-i\theta_j}\}.
```

Put $`p=N/2`$ and $`Q=2^{N-4}`$. Choose the positive-angle grid

```math
\gamma_j=\frac\pi4+\frac{(j-1)\pi}{2Q},
\qquad 1\le j\le Q.
```

Its positive eigenphases have pairwise chord distance at least
$`1/Q=16\eta`$. Each $`p`$-element subset of this grid, assigned in
increasing order to the final-depth angles, gives a distinct target
spectrum. There are $`\binom Qp`$ such spectra.

Every selected positive eigenphase of $`W\otimes I_b`$ has an eigenspace
of dimension $`2^b`$. The preceding multiplicity argument applied to its
image under $`E_\theta`$ forces at least $`2^b`$ eigenvalues of
$`C_\theta`$, counted with multiplicity, in its radius-$`2\eta`$ disk.
The grid disks are disjoint. A single physical Clifford spectrum therefore
supports at most

```math
\frac{2^M}{2^b}=2^{n+a}\le4N
```

grid phases, and can serve at most $`\binom{4N}p`$ of the selected target
spectra. Counting literal Clifford spectra gives

```math
\binom{Q}{N/2}
\le2^{3M+7}\binom{4N}{N/2}
\le2^{3M+7+4N}.
```

On the other hand, $`\binom Qp\ge(Q/p)^p`$, whose base-two logarithm is
$`(N/2)(N-n-3)`$. Combining the bounds proves the proposition. ∎

The zero upper angles are admitted singular state charts: the full frame
still retains every final-depth rotation. The same count holds with
strictly positive upper angles. Set them all to
$`0\lt\delta\le\eta/n`$; telescoping the upper layers changes the full
operator by at most $`\eta`$. An $`\eta`$ implementation is then a
$`2\eta`$ intertwiner for the final-depth-only target. Radius-$`4\eta`$
spectral disks still remain disjoint on this grid, giving the same bound.

For general accuracy and clean count, a grid of
$`Q=\Theta(1/\eta)`$ gives the necessary inequality

```math
\binom{Q}{N/2}
\le2^{3M+7}\binom{2^aN}{N/2},
```

whenever $`Q\ge N/2`$. Thus

```math
3M\ge\frac N2
\bigl(\log_2(1/\eta)-a-\log_2N-O(1)\bigr)-7.
```

This restriction comes from the spectral multiplicities required by
arbitrary dirty inputs. A fixed one-qubit Clifford conjugation makes this
final-depth-only family diagonal, so the
[packed diagonal compiler](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames)
already handles it at linear T-count. It therefore does not identify an
unrestricted hard family. Several interleaved
non-Clifford operations need not have Clifford spectra, and remain outside
the proposition. The Spin covariance, reflection, and Choi identities in
Section 10 have different interfaces and retain their separate statements.

### Exact return does not erase the value of dirty work

The [complete-return rigidity result](../../docs/FAULT_TOLERANT_COMPILER.md#104-rigidity-of-approximate-workspace-return)
can force an approximate native implementation to return all its work
exactly at a sufficiently small T-count. Exact return still does not
permit removing the dirty wires without increasing that count.

Let $`S=2^{n-1}=N/2`$ and let $`f`$ be any Boolean table on
$`n-1`$ bits. Set all earlier angles to zero and use
$`\theta_{n-1,p}=\pi f(p)`$. Since $`R_y(\pi)=-I_2`$ literally,
the resulting full Hopf frame is

```math
W_f=D_f\otimes I_2,
\qquad D_f=\mathrm{diag}_p((-1)^{f(p)}).
```

These are $`2^S`$ distinct literal targets. The exact
[dirty-bank SelectSwap query](../../docs/OPERATOR_SOURCE_COMPILER.md#7-trading-additional-dirty-banks-for-lookup-cost)
implements the phase-free XOR oracle for this one-bit table at cost

```math
T=O(S/\lambda+\lambda),
\qquad b\le\lambda+n-1
```

in additional arbitrary banks and selectors, for a power of two
$`1\le\lambda\le S`$. The separate query-output qubit can be one supplied
clean qubit prepared in $`|-\rangle`$. Phase kickback, followed by the
inverse of that Clifford preparation, implements $`W_f`$ exactly.
The query returns every dirty bank and selector exactly on its full input
space, so this construction also returns the clean qubit and preserves all
external-reference correlations. Choosing $`\lambda`$ within a factor of
two of $`\sqrt S`$ gives

```math
T=O(\sqrt N),\qquad a=1,\qquad b=O(\sqrt N+n),
```

within the permitted allocation for all sufficiently large $`n`$.

By contrast, the
[fixed-width coherent word count](../../docs/FAULT_TOLERANT_COMPILER.md#102-fixed-width-coherent-counting)
bounds the number of literal native circuits on $`q=n+2`$ wires with at
most $`t`$ T or T-dagger gates by

```math
2^{2q^2+3q+5+(2q+1)t}.
```

Allowing two clean wires already includes every circuit with fewer clean
wires by padding identities. One exact circuit cannot implement two
different $`W_f`$ on that initialized-input subspace. Covering all
$`2^{N/2}`$ Boolean targets consequently requires, in the worst case,

```math
t\ge\frac{N/2-(2q^2+3q+5)}{2q+1}
=\Omega(N/n).
```

This exceeds any fixed constant times $`\sqrt N`$ for sufficiently large
$`n`$. Hence an exact-return circuit cannot in general be compressed to
$`n+2`$ wires at unchanged T-count, even within the prescribed Hopf family.
The rigidity theorem concerns return, while this comparison concerns the
cost of deleting returned workspace; neither implies a superlinear lower
bound in the allowed dirty-work model.

<a id="9-finite-checks-and-evidence-limits"></a>
<a id="10-finite-checks-and-evidence-limits"></a>
<a id="11-finite-checks-and-evidence-limits"></a>

## 12. Finite checks and evidence limits

Run:

```bash
python -m unittest discover -s tests -p 'test_source_reuse_limits.py'
python -m unittest discover -s tests -p 'test_tree_residual_structure.py'
python -m unittest discover -s tests -p 'test_source_merge.py'
python -m unittest tests.test_shared_conjugator_merge
python -m unittest tests.test_precision_carry
python -m unittest tests.test_correlated_precision_carry
python -m unittest tests.test_small_product_compilation
python -m unittest tests.test_endpoint_structural_limits
python -m unittest tests.test_endpoint_refinement_routes
```

The tests exercise the dimension inequality on small nilpotent
contractions and encoded subspaces, the explicit geometric construction,
the need for nilpotence, and the transformed-mask correlations.
They also check the literal paired-source hoist, its fixed-tail
cancellation, the low-rank reflection kernel witness, and the matched-frame
commutation constraint. These finite checks support indexing and assumption
boundaries.
The common-conjugator fixture additionally checks every dirty and occupied
signal column of the literal native fork, its target-half-block failure,
and its source counts. Small mask sweeps check the ellipse and row
constraints; the precision-independent retuning bound is proved above.
The width-transition fixtures check actual native loader phases, bridge
orientation and gate counts, complete multi-group identities, and the
legal grouped coefficient witness. They do not establish a cheaper
implementation of the transformed group bodies. The flag-correlated
fixtures check literal chain phases, three unequal widths, arbitrary
released tails, a realizable scalar query, and full syndrome extraction
with the actual inverse and nonzero approximation leakage. They do not
implement a complete group in the proposed code. The small-product checks
compare quaternion multiplication with complete addressed matrices, retain
noncommuting factors and mixed path amplitudes, and reject rowwise
compression after an address change. They are examples of the existing
fixed-address compiler, not a new general precision-sharing theorem.
The dimension and T-count statements rest on the analytic proofs above;
the tests do not establish an unrestricted impossibility theorem or
literature priority.

The [endpoint structural fixtures](../../tests/test_endpoint_structural_limits.py)
also check the noncommuting dirty-word counterexample, the valid-tuple
SELECT disturbance and its clean compression, the actual Majorana word's
first and second moments, and the complete phase-halving action on prepared
and arbitrary catalysts. These checks verify the stated
interfaces; they do not convert prepared-source promises into dirty ones.

The [refinement-route fixtures](../../tests/test_endpoint_refinement_routes.py)
also check the common Clifford-vector algebra, the complete dirty-gauge
identity, and the unequal-address conjugation defect. The identities above
are full-operator statements; the scoped Spin obstruction is not an
unrestricted T-count lower bound.

The [tree-residual fixtures](../../tests/test_tree_residual_structure.py)
independently compare addressed circuits, recursive subtree vectors, and
the scalar-generator reconstruction for two- and three-level frames.
They include actual complex native words, general complex SU(2) words,
singular target angles, and the explicit path-product witnesses. These
small checks validate the identities and marker conventions; they do not
construct a lower-cost coherent coefficient evaluator.
