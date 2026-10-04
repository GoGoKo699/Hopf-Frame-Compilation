# One precision charge for sparse nested tree updates

**Proved restricted case.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Research status](../../docs/OPEN_PROBLEM.md) · [Prefix-free updates](ANTICHAIN_COMPILER.md) · [Conditional workspace](../../docs/CONDITIONAL_SUFFIX_COMPILER.md)

Arbitrarily deep nesting can still admit a single precision charge. When
the ancestor closure of the changed nodes is small, the nonnative part of
the tree acts on a small set of computational basis states. An exact
permutation packs those states into a smaller logical register. The
remaining logical bits supply conditional initialized workspace for one
dense residual block.

The resulting promised-family bound is

```math
a=1,\qquad b\ge L+n+7,\qquad
T=O(N+L),\qquad G=O(NL).
```

In particular, arbitrary changes along one root-to-leaf path satisfy this
bound, regardless of the path length. Unlike the prefix-free theorem,
this construction uses one initialized compiler qubit. It does not give
the same bound for a general frame: the dimension of the changed part,
rather than nesting alone, is the restriction.

## 1. Promise and theorem

Use the parameters and literal native baseline of the
[prefix-free update theorem](ANTICHAIN_COMPILER.md#1-target-and-native-baseline-promise).
Thus $`N=2^n`$, $`n\ge1`$, $`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
Each depth-d baseline word $`U_v^C`$ is an actual determinant-one
Clifford+T word of length at most $`c_0(n-d+1)`$, with fixed
$`c_0`$. Local gates act on the ordered pair

```math
|v0^{n-d}\rangle,\qquad |v1\,0^{n-d-1}\rangle.
```

Let S be a nonempty, ancestor-closed set of internal nodes. Every strict
ancestor of a node in S also belongs to S. The target has arbitrary
effectively specified SU(2) local words on S and agrees literally with
the supplied baseline outside S. One may take S to be the ancestor
closure of the nodes that change; the parameter below counts that
closure, including its unchanged ancestors. Set

```math
t=|S|,\qquad m=t+1,\qquad
s=\lceil\log_2m\rceil,\qquad M=2^s.
```

**Theorem.** If

```math
2s+32\le n,
```

then there is a coherent Clifford+T circuit $`\widetilde W`$ using one
initialized compiler qubit and at most $`b=L+n+7`$ arbitrary dirty
qubits such that

```math
\|\widetilde WJ_1-J_1(W\otimes I_b)\|\le\eta,
\qquad T=O(N+L),\qquad G=O(NL).
```

Here $`J_1`$ initializes only the external compiler qubit. The norm
covers every logical and dirty input, arbitrary references, initialized
leakage, and dirty-work return. It is an isometry guarantee on that
initialized input, not a full-operator guarantee for an arbitrary value
of the external qubit. Larger work allocations may leave unused wires
untouched. The constant 32 is a sufficient allocation, not an optimized
threshold.

**Path corollary.** If all changes lie on one root-to-leaf path, the same
bound holds for every n, without the displayed sufficient-size condition.
This includes arbitrary angles at every internal node of that path.

There are no measurements, resets, supplied resource states, or uncharged
quantum oracles. Certified coefficient evaluation and classical circuit
and table preparation are separate preprocessing costs. No bound on that
classical running time is asserted.

## 2. Factor the native forest and identify the active modes

Let $`Q_S`$ be the product of the target gates in S, in the prescribed
shallow-to-deep order. Let F contain the baseline gates outside S in
their original relative order. A node outside S cannot be an ancestor
of a node in S. Consequently, moving an outside gate left past later
S-gates exchanges only incomparable gates, whose logical supports are
disjoint. Thus the exact complete-operator identity is

```math
W=FQ_S.
```

The native mask construction in
[the prefix-free proof](ANTICHAIN_COMPILER.md#exact-compilation-of-the-baseline-and-its-masks)
implements F with $`O(N)`$ T and Clifford gates and at most
$`n+2`$ arbitrary dirty wires. Its helpers return exactly. F is
compiled from the supplied local native words; it is not a free oracle.

Define the computational basis labels

```math
\mathcal B_S=\{0^n\}\ \cup
\{v1\,0^{n-|v|-1}:v\in S\}.
```

These are m distinct labels: the rightmost one in a nonzero marker
identifies its node. The right endpoint of every S-gate is its own marker.
Its left endpoint is either $`0^n`$ or the marker belonging to the
ancestor at the last one in its prefix. Ancestor closure puts that
ancestor in S. Hence

```math
Q_S=I\quad\text{on}\quad
\mathrm{span}\{|x\rangle:x\notin\mathcal B_S\}.
```

This support statement concerns the complete operator, including all
frame columns. It does not retain only a prepared state or its first
column. It also holds for complex, noncommuting local SU(2) words.

## 3. Pack the support with a charged exact permutation

Choose a computational basis permutation P taking the ordered labels in
$`\mathcal B_S`$ to the integers $`0,\ldots,m-1`$. It can be
constructed with at most m transpositions: for each source label in turn,
swap its current image with its desired image. Earlier desired images
are distinct and remain fixed. P may move other basis states, whose
temporary motion cancels upon applying its actual inverse.

Each transposition of distinct n-bit labels a and b has an exact
$`O(n^2)`$ implementation with one returned dirty helper. First use
X gates to send a to zero. Choose a pivot where $`a\oplus b`$ is
one, and use CNOTs from that pivot to the other nonzero positions. This
invertible affine circuit sends b to the pivot's one-hot label. A
negative-controlled X on the pivot, conditioned on every other logical
bit being zero, swaps just those two labels. Undo the affine circuit.
The retained borrowed-MCX construction charges the middle operation and
returns its helper on all inputs. Thus

```math
T(P),G(P)=O(mn^2).
```

Write the packed logical register as r upper bits and s lower data bits,
where $`r=n-s`$. There is an M-dimensional unitary U such that

```math
P Q_S P^\dagger
=|0^r\rangle\langle0^r|\otimes U
+(I-|0^r\rangle\langle0^r|)\otimes I_M.
```

Inside U, the unused labels $`m,\ldots,M-1`$ are also unchanged.
The added identity block introduces no determinant or scalar-phase
adjustment.

Compute $`h=[\text{upper }r\text{ bits}=0]`$ into the external
clean qubit. This costs $`O(n^2)`$ native gates and one returned dirty
helper. On $`h=1`$, the upper logical bits are zero and can be used
as private initialized work. On $`h=0`$, every complete subroutine
below is exact identity on its full input space, including arbitrary
private bits. This is the
[conditional-suffix contract](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally).

## 4. An exactly native coarse unitary and a dense column dictionary

Fix

```math
K=4M,\qquad \delta=\frac{1}{8K\sqrt M}.
```

Build an actual native coarse circuit $`C_0`$ for U with
$`\|C_0-U\|\le\delta`$. U is the product of t known two-level
SU(2) gates on the packed labels. Approximate each local word within
$`\delta/t`$, with literal phase-calibrated native length
$`O(\log M)`$, and use exact permutations and sector selection to
address its pair. The three Euler rotations may be approximated with the
retained phase-cancelling words; their common synthesis scalars are not
dropped. This gives an exactly unitary native $`C_0`$, with exactly
returned helpers and total cost $`t\,\mathrm{poly}(n)`$.

Its h-conditioned implementation uses the retained determinant-one
involution rewrite and dirty predicate echo, as in the exact native mask
construction. Each addressed reflection is completed before the next is
started. This charges literal controls and preserves the inactive sector
exactly; arbitrary controlled T gates are not assumed to be free.

Define

```math
V=C_0^\dagger U,\qquad E=V-I_M,\qquad \|E\|\le\delta.
```

For each column j, split the entries of E into four nonnegative phase
components, using $`\omega\in\{1,-1,i,-i\}`$:

```math
E_{ij}=\sum_\omega\omega a_{ij\omega},
\qquad a_{ij\omega}\ge0,\qquad
c_{ij\omega}=K\sqrt M\,a_{ij\omega}\le\frac18.
```

The split takes positive and negative parts of the real and imaginary
components. Certified interval evaluation and clipping give these
coefficients to any prescribed finite tolerance; exact sign or zero
tests are unnecessary. Evaluate each scaled coefficient within
$`2^{-q}/16`$, retaining the interval $`[0,1/4]`$ by clipping.
The actual coarse circuit is fixed before evaluating E.

The K term labels are the pair $`(j,\omega)`$, using s column bits
and two phase bits. For each j define the rank-one column map

```math
D_j=H^{\otimes s}X^j|j\rangle\langle j|
=\frac1{\sqrt M}\sum_{i=0}^{M-1}|i\rangle\langle j|.
```

Its unitary dilation uses one distinct atom rejection flag. Starting
that flag at zero, toggle it by the predicate $`[\mathrm{data}\ne j]`$,
then apply $`X^j`$ to the data and finally $`H^{\otimes s}`$.
The accepted block is exactly $`D_j`$. A comparison uses XORs,
a borrowed multi-controlled operation, and their actual inverses; its
temporary comparison work is erased before changing the data word.

With $`\mathrm{diag}(c_{\cdot j\omega})`$ acting after
this map, the exact dictionary identity is

```math
\frac1K\sum_{j,\omega}
\omega\mathrm{diag}(c_{\cdot j\omega})D_j=E.
```

The scalar table has $`KM=4M^2`$ entries. It is dense only within
the packed active register.

## 5. One scalar source, separate rejection flags, and full error

Use a separate initialized scalar flag for the
[dirty scalar-source block](../../docs/OPERATOR_SOURCE_COMPILER.md#3-dirty-programming-and-a-one-flag-scalar-block).
Its q arbitrary dirty core bits and k dirty query selectors give a scalar
accepted coefficient on the entire dirty space. Its native rounding error
is at most $`(5/2)2^{-q}`$ per coefficient. Including the certified
evaluation tolerance above gives

```math
|\widehat c_{ij\omega}-c_{ij\omega}|\lt3\,2^{-q}.
```

For each term, perform the atom dilation first, then the scalar block
addressed by the **current output data** i and the unchanged label
$`(j,\omega)`$. The scalar block acts on a different rejection flag
and leaves both its data address and the atom flag unchanged. Consequently
the accepted product is

```math
\mathrm{diag}(\widehat c_{\cdot j\omega})D_j
\otimes I_{\rm dirty}.
```

Every actual scalar query and its inverse completes while i and the term
label remain unchanged. The original input j is supplied by the label;
the query does not try to recover an overwritten data address. The two
separate rejection flags justify the displayed projected product. No
intermediate rejected amplitude is discarded.

Apply the literal phase outside the scalar/atom product. For binary phase
bits u and v, use $`\omega=(-1)^u i^v`$. Its h- and mode-conditioned
implementation is charged exactly. In particular, a controlled S phase
$`i^F`$ is obtained by toggling a dirty bit z by F, applying S,
untoggling z, and applying $`S^\dagger`$. This gives
$`i^F(-1)^{zF}`$; an exact multiple-controlled Z corrects the last
factor. One additional returned MCX helper suffices. Controlled Z is
already a multiple-controlled phase. These operations do not assume
initialized phase scratch.

Prepare a mode bit f and the K-label register uniformly. On $`f=0`$
apply identity, and on $`f=1`$ apply the phase, atom, and scalar
SELECT; undo the preparations. Let Q be this actual half-block word.
All preparations, atoms, and complete scalar words are h-conditioned.
With J initializing the private upper bits, its active accepted block B
satisfies

```math
B=J^\dagger QJ
=\frac12\left(I_M+
\frac1K\sum_{j,\omega}
\omega\mathrm{diag}(\widehat c_{\cdot j\omega})D_j\right)
\otimes I_{\rm dirty}.
```

Since $`\|D_j\|=1`$, each rounded atom has error less than
$`3\,2^{-q}`$. Averaging gives

```math
\|2B-V\otimes I_{\rm dirty}\|\lt3\,2^{-q}.
```

There is only one shared scalar block in SELECT. It is addressed by
$`(i,j,\omega,h,f)`$, so its queries have
$`k=2s+4`$ bits. Inactive query rows contain the zero **XOR word**.
Every controlled source is implemented by unconditional source basis
changes around its controlled central X; those changes cancel when the
control is inactive. This makes each complete h-inactive scalar word
identity even on rejected scalar-flag inputs and arbitrary dirty states.
The source is not copied once per column or term.

Let $`R_h`$ be the identity on $`h=0`$ and, on $`h=1`$, the
reflection about the complement of the all-zero private upper register.
It is a charged $`O(n^2)`$ borrowed multi-controlled phase. Define

```math
\mathcal A_h=Z_h Q R_h Q^\dagger R_h Q.
```

On $`h=1`$, $`Z_h`$ supplies the literal minus sign in the
[normalization-two amplification lemma](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error).
On $`h=0`$, the entire word is identity. The lemma bounds the
complete active isometry, including rejected components, by

```math
\|\mathcal A_hJ-J(V\otimes I_{\rm dirty})\|
\lt12\,2^{-q}.
```

Apply the actual h-conditioned coarse $`C_0`$, then reverse the
original computation of h. Ideal outputs return all private logical bits
to zero before that uncomputation. Actual leakage is propagated by
unitaries, so neither step increases its norm. The inactive sector remains
exact identity. Choose

```math
q=L+6,\qquad
12\,2^{-q}=\frac3{16}\,2^{-L}\le\eta.
```

Finally apply $`P^\dagger`$ and F. They are actual exact circuits
on the full logical and dirty space. Using
$`W=FP^\dagger(PQ_SP^\dagger)P`$ proves the theorem's complete
initialized-input error bound. All inverses above are actual circuit
inverses.

## 6. Workspace and charged costs

Only h is initialized externally. On its active sector the private upper
register uses

| Private role | Logical zero bits |
|---|---:|
| Column label j | s |
| Literal phase label | 2 |
| Identity/correction mode f | 1 |
| Atom rejection flag | 1 |
| Scalar-source flag | 1 |
| Total | $`s+5`$ |

No initialized q-bit output word or endpoint register is present. The
remaining upper bits stay zero on ideal active inputs.

Reserve eight additional arbitrary dirty helpers, separate from the core
and selectors. The multi-controlled predicate, comparison, and reflection
circuits need one returned MCX helper at a time. A literal conditional
S phase uses one dirty echo bit and that MCX helper. If a conditional T
phase occurs in a fixed control decomposition, the retained nested T/S
echo uses two dirty echo bits and one MCX helper. The coarse reflection
interpreter additionally uses a dirty predicate bit, returned before the
next reflection. Fixed controlled H gates use a constant native
conjugation of controlled Z, with opposite synthesis scalars cancelling;
source central controls use the same returned MCX helper. These steps are
sequential. Even reserving their four distinct echo/predicate roles,
one MCX helper, and three spare helpers uses at most eight wires. None is
assumed zero, and none is also counted as a private flag or a table selector.

Thus the dense step fits

```math
\begin{aligned}
b_{\rm dense}
&\le q+k+8\\
&=L+2s+18\\
&\le L+n+7,
\end{aligned}
```

and the logical private requirement fits because
$`r=n-s\ge s+32\ge s+5`$. The native forest, support permutation,
coarse circuit, and predicate computations use the same pool sequentially;
their helpers return exactly on arbitrary inputs. The stronger sufficient
slack 32 avoids optimizing the fixed control allowance.

There are a constant number of complete source and query calls, including
amplification. With $`2^k=16M^2`$, their costs are

```math
T=O(M^2+L+\mathrm{poly}(n)),
\qquad
G=O(M^2L+\mathrm{poly}(n)).
```

Add the native forest $`O(N)`$, exact packing $`O(mn^2)`$, and
the coarse cost $`t\,\mathrm{poly}(n)`$. The hypothesis gives
$`M^2\le2^{n-32}`$ and $`t\lt M\le\sqrt N`$.
For every fixed-degree polynomial p, $`\sqrt N\,p(n)=O(N)`$,
with an absolute constant. Consequently the complete circuit satisfies

```math
T=O(N+L),\qquad G=O(NL).
```

This is an asymptotic statement with sufficient fixed constants. It does
not claim a practical advantage at the allocation threshold, an optimal
Clifford count, or a T-depth improvement.

## 7. Arbitrarily long paths and the remaining restriction

For changes on one root-to-leaf path, its ancestor closure has
$`t\le n`$. Thus $`s\le\lceil\log_2(n+1)\rceil`$, and
the sufficient condition holds for every sufficiently large n. The n
that fail it form an absolute finite set. For those cases, implement the
at most n target local gates in the ancestor closure in their original
tree order. Synthesize its changed words to error $`\eta/n`$, using
literal native sector echoes and exactly returned borrowed helpers; copy
its unchanged native words exactly. Then apply the exact off-S forest F
from Section 2. Their total cost is
$`O(N+n(L+\log n)\mathrm{poly}(n))`$. Since n is bounded in
this fallback, it is $`O(N+L)`$ with an absolute constant and fits
the stated workspace. The empty update is just the supplied baseline.
This proves the path corollary for all n.

For a real Hopf specialization, take the branching native baseline
$`U_v^C=R_y(\pi/4)=HZ`$ and arbitrary real rotations at the
promised changed nodes. The construction approximates the complete real
frame, so its existing fixed-parameter
[QBP substitution guarantee](../../docs/QBP_APPROXIMATION.md) applies. A generic
complex native baseline instead gives a structured SU(2) tree operator;
it is not automatically a prescribed real Hopf frame or the separate
phase-dressed magnitude family.

At $`L=N`$ the promised family therefore has linear T count with one
clean compiler qubit and $`b=N+n+7`$ dirty qubits. It permits
arbitrary-depth nested changes and sparse branching. The sparse dimension
condition is essential to this proof: a general changed tree has
$`t=N-1`$, hence $`s=n`$, and supplies no packed zero suffix.
Neither this theorem nor the prefix-free update theorem closes that
unrestricted endpoint. No matching lower bound for the promised family
is claimed.
