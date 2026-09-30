# Source reuse and structured residuals

[Open endpoint](OPEN_PROBLEM.md) · [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Shared-source compiler](FAULT_TOLERANT_COMPILER.md)

The open endpoint remains
$`a=2`$, $`b=N+n+7`$, $`L=N`$, $`n\ge3`$, with
$`\Omega(N)\le T^\star_{F,\mathbb R}\le O(N\ell_*(n))`$
by [conditional-suffix grouping](CONDITIONAL_SUFFIX_COMPILER.md).
Here $`\ell_*(n)=1+\log_2^*(n+2)`$, and $`\log_2^*`$ counts
base-two logarithms until the value is at most one.
This note gives scoped restrictions on proposed ways to reuse precision work
and an exact tree factorization of the residual coefficient data. The
factorization separates classical compression from the still-charged task
of coherent tree transport. None of these statements is a lower bound for
an unrestricted Hopf-frame compiler or changes the retained resource bounds.
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
transport misses necessary mixed terms. Thus the remaining target is a
jointly charged block for the weighted operator, including its occupied
flags and literal inverse, rather than classical compression alone.

## 4. Consequence for the research direction

These results rule out two particular shortcuts: replacing the
initialized nilpotent geometric source by a constant-clean dirty
encoding with the same full-output property, and treating all
source-conjugated programming masks as cheap Clifford operations.

The sufficient [whole-frame half-unitary block](OPEN_PROBLEM.md#a-sufficient-construction-to-seek)
remains a valid target. A successful construction may avoid intermediate
source return, may use a different carried operator, or may synthesize
the interleaved source and programming operations jointly.
The bounds for unrestricted complete-frame compilation remain unchanged.

The tree factorization above provides a compact classical starting point
for another construction, while keeping its coherent implementation as an
explicit unresolved task. A linear generator count is not yet a linear
T-count.

### A scoped diagnostic for Pauli routing of rejected components

A separate proposed shortcut routes successive scalar-filter rejections
through the four states of two flags. For the common operator source M,
the [scalar block](OPERATOR_SOURCE_COMPILER.md#3-dirty-programming-and-a-one-flag-scalar-block)
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
calls. The [native routing fixtures](../tests/test_source_merge.py)
check the explicit word and the stated class boundary.

## 5. Finite checks and evidence limits

Run:

```bash
python -m unittest discover -s tests -p 'test_source_reuse_limits.py'
python -m unittest discover -s tests -p 'test_tree_residual_structure.py'
python -m unittest discover -s tests -p 'test_source_merge.py'
```

The tests exercise the dimension inequality on small nilpotent
contractions and encoded subspaces, the explicit geometric construction,
the need for nilpotence, and the transformed-mask correlations.
These finite checks support indexing and assumption boundaries.
The dimension and T-count statements rest on the analytic proofs above;
the tests do not establish an unrestricted impossibility theorem or
literature priority.

The [tree-residual fixtures](../tests/test_tree_residual_structure.py)
independently compare addressed circuits, recursive subtree vectors, and
the scalar-generator reconstruction for two- and three-level frames.
They include actual complex native words, general complex SU(2) words,
singular target angles, and the explicit path-product witnesses. These
small checks validate the identities and marker conventions; they do not
construct a lower-cost coherent coefficient evaluator.
