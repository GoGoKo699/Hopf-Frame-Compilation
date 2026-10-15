# Exact depth for deterministic Hopf frames

[Homogeneous magnitudes](HOMOGENEOUS_MAGNITUDE_DEPTH.md) · [Preserved-prefix recursion](RECURSIVE_TERMINAL_DEPTH.md) · [One-clean query and fallback contracts](ONE_CLEAN_DEPTH.md#2-complete-query-replacements-and-returned-predicates)

**Local result R84.** A complete real Hopf frame whose local rotations
are signed quarter turns has an exact compiler with one clean qubit and
additive T-depth $`O(\log(n+2)(1+\log^*n))`$. The orientation and sign
may vary independently at every node. A tagged descendant source is an
exact signed computational string; one classical last-deviation query
then erases its tag. The protected logical suffix supplies the short
clean tag, and recursion completes every marker column.

This is a corner-family result. It does not compile an arbitrary
in-place tree automorphism for free, extend the classical tag eraser to
positive overlapping descendants, or improve the unrestricted R79
bound. Zero incoming amplitudes do not remove the prescribed columns
inside the corresponding subtrees.

## 1. Exact statement

Let $`n\ge1`$ and $`N=2^n`$. At every node $`v`$ the prescribed rotation
has the form

```math
R_v=\begin{pmatrix}c_v&-s_v\\s_v&c_v\end{pmatrix},\qquad
(c_v,s_v)\in\{(1,0),(-1,0),(0,1),(0,-1)\}.
```

For every integer dirty allocation $`b\ge17(n+13)`$, there is one native
Clifford+T circuit $`V`$ with one external initialized qubit and at most
$`b`$ arbitrary dirty qubits such that

```math
VJ_1=J_1(W\otimes I_b),
```

```math
\begin{aligned}
T&=O\!\left(\sqrt N+\frac Nb\right),\\
G&=O(N),\\
D_T&=O\!\left(\frac N{b^2}
  +\log_2(n+2)(1+\log^*n)\right).
\end{aligned}
```

Here $`W`$ is the entire prescribed Hopf frame and $`J_1`$ initializes
only the sole external signal. The equality includes literal signs,
every logical column, exact work return, arbitrary dirty inputs, and
arbitrary reference correlations. All constants are absolute. Unused
supplied wires remain untouched. As elsewhere, $`G`$ counts elementary
Clifford gates, including gates between T layers; Clifford depth is not
part of the T-depth objective.

At $`b=\Theta(\sqrt N)`$ the same circuit therefore has
$`T=O(\sqrt N)`$, $`G=O(N)`$, and
$`D_T=O(\log(n+2)(1+\log^*n))`$. There is no accuracy parameter because
the compiler is exact. The promise on the four possible local pairs
makes their identification effective: certified evaluations to absolute
error less than $`1/4`$ distinguish all four cases. No undecidable sign
or zero test is requested. Classical preprocessing is terminating,
without an efficiency claim. No QRAM, measurement, reset, supplied
precision state, or free coherent evaluator is used.

## 2. Signed leaves and the last-deviation inverse

For each node define $`r_v\in\{0,1\}`$ as its unique nonzero scaling
child and $`\sigma_v\in\{1,-1\}`$ as that coefficient. Then

```math
R_v=\sigma_v(XZ)^{r_v}.
```

The operator order is literal: chronological $`Z`$ followed by $`X`$
gives $`XZ`$. For a local subtree of height $`h`$ beneath a preserved
outside prefix $`P`$, write its root input as $`0^h`$ and its marker
inputs as

```math
j=p1\,0^{h-t},\qquad |p|=t-1,\qquad 1\le t\le h.
```

The root output follows the favored child $`r_v`$ at every step, with
sign equal to the product of $`\sigma_v`$ along that path. For a marker,
all earlier layers are inactive because the marker bit has not yet been
visited. At its node $`Pp`$, the column $`R_{Pp}|1\rangle`$ chooses the
other child, with coefficient $`(-1)^{r_{Pp}}\sigma_{Pp}`$. Every later
layer follows the favored child. Thus its output leaf and sign are

```math
z(P,p,t)=p\,(1-r_{Pp})\,u(P,p,t),
```

```math
\chi(P,p,t)=(-1)^{r_{Pp}}\sigma_{Pp}
   \prod_{v\text{ on the later favored path}}\sigma_v.
```

The string $`u`$ has length $`h-t`$. Empty products are one. For $`t=0`$
use the root leaf and root sign just described, with empty preserved
local prefix. Crucially, marker signs have no factor from ancestors
above $`Pp`$. These formulas therefore include markers below every zero
incoming edge.

For an arbitrary output leaf $`x\in\{0,1\}^h`$, define

```math
\delta_j(P,x)=x_j\mathbin\oplus r_{P x_1\cdots x_{j-1}},
```

```math
d(P,x)=\max\bigl(\{j:\delta_j(P,x)=1\}\cup\{0\}\bigr).
```

The root output has no deviation. The output of marker $`p1\,0^{h-t}`$
has a deviation at position $`t`$ and follows favored children after
that position, so its last deviation is exactly $`t`$. Its preceding
prefix is still $`p`$. Hence the unsigned leaf map $`\pi_P`$ obeys

```math
d(P,\pi_P(j))=\mathrm{lastone}(j).
```

Conversely, the last deviation of any leaf reconstructs its unique
marker, or the root if there is no deviation. This proves the signed
permutation description for every column. Computing the finite table
$`d`$ is classical preprocessing. A coherent XOR query to that table
will be priced explicitly; it is not an in-place inverse permutation
oracle.

## 3. Exact selected involutions and signed string sources

Let $`H`$ be a preserved activity bit, $`T`$ a short tag, and
$`X_1,\ldots,X_h`$ the logical data. They are disjoint registers. Put

```math
s=\lceil\log_2(h+1)\rceil.
```

The selected involution banks of
[R80, section 3](HOMOGENEOUS_MAGNITUDE_DEPTH.md#3-actual-one-qubit-words-and-selected-involution-banks)
apply a one-qubit involution on positions selected by a Boolean mask of
$`H,T`$. A completed mask query XORs its values into arbitrary dirty
outputs; a controlled-involution bank, the actual inverse query, and
the same controlled bank cancel the dirty contribution. The address
$`H,T`$ is read only and is disjoint from every target. Private
non-Clifford supports permit parallel selected operations; shared
address interactions are Clifford gates only.

Use these banks both for the single-bit mask $`H[t=j]`$, which selects
$`X_j`$, and for the suffix mask $`H[t\le j]`$, which selects a Hadamard
on $`X_j`$. Tag zero selects the whole suffix bank and no marker-clear
bit. Padded tags above $`h`$ select nothing. Controlled Hadamard is an
exact constant-size Clifford+T word, as in R80; the selected-Hadamard
bank requires no angle approximation. The completed banks have
polynomial local count and dirty width and $`O(\log(h+2))`$ T-depth.
They are exactly identity on $`H=0`$ even for arbitrary occupied tags
and dirty inputs.

For a valid tag and preserved label $`(P,p,t)`$, write $`z_{\mathrm{suf}}`$
for the prescribed output suffix from section 2 and write
$`\chi=(-1)^\kappa`$. On variable suffix bits $`y`$, define

```math
f(P,p,t,y)=\kappa(P,p,t)
    \mathbin\oplus\bigl(z_{\mathrm{suf}}(P,p,t)\cdot y\bigr).
```

The dot product is over $`\mathbb F_2`$. If $`\mathcal H_t`$ denotes
Hadamards on precisely that suffix, the exact source identity is

```math
\mathcal H_t\,\mathrm{diag}((-1)^f)\,\mathcal H_t
  =\chi(P,p,t)X^{z_{\mathrm{suf}}(P,p,t)}.
```

It is an operator identity on the entire variable suffix. In particular,
it maps its zero input to the required signed computational string.
The prefix and tag are unchanged throughout. The sign term is part of
the same Boolean phase; it is not discarded as a global phase of a
single column.

To implement the middle phase, use one arbitrary dirty bit $`a`$ and
a completed XOR query $`Q_{Hf}`$ whose address is $`P,T,X,H`$. The
chronological word

```math
Q_{Hf},\quad Z_a,\quad Q_{Hf}^{\dagger},\quad Z_a
```

is exactly multiplication by $`(-1)^{Hf}`$ and returns $`a`$ and all
query helpers on every input. Assign zero phase to padded tags. On
$`H=0`$ both selected Hadamard banks and the completed phase word are
identity. The two copies of $`\mathcal H_t`$ are actual inverses.

There is only one long data register. At fixed $`t`$, its prefix bits
encode $`p`$ and its remaining bits encode $`y`$. If the outside prefix
has $`a_{\mathrm{ext}}`$ bits, the source phase table consequently has
at most

```math
C(h+1)2^{a_{\mathrm{ext}}+h}
```

rows, including the activity bit and tag padding. No second copy of
$`X`$ or initialized long index is hidden in this estimate.

## 4. A complete active stage

Split the remaining $`m`$ logical bits into $`h=m-r`$ early bits $`X`$
and a protected suffix $`Z`$ of $`r`$ bits. The previously completed
prefix $`P`$ has length $`n-m`$ and remains read only during this stage.
The sole initialized external signal first saves

```math
H=[Z=0^r].
```

Use the exact returned-two-helper conjunction and its actual inverse
from [R72, section 2](ONE_CLEAN_DEPTH.md#2-complete-query-replacements-and-returned-predicates).
Both arbitrary helpers are borrowed between completed queries. Choose
$`r\ge s`$, and let $`T`$ occupy $`s`$ designated bits inside $`Z`$.
On the active branch $`T=0`$; on the inactive branch it may be arbitrary.
The remaining bits of $`Z`$ stay untouched.

The early stage has the following chronological word:

1. Apply $`F_H`$, the exact XOR of $`H\,\mathrm{lastone}(X)`$ into $`T`$.
2. Clear the marker bit using the selected $`X_j`$ mask $`H[t=j]`$.
3. Apply the signed string source of section 3, guarded by $`H`$.
4. Apply the completed output query $`T\mathrel{\oplus}=H\,d(P,X)`$.
5. Reverse the actual saved-predicate circuit for $`H`$.

The [last-one XOR word](HOMOGENEOUS_MAGNITUDE_DEPTH.md#4-exact-last-one-xor-with-arbitrary-dirty-work)
uses read-only conjunctions and dirty echoes. It preserves $`X,H`$,
works on arbitrary initial $`T`$, costs $`O((h+2)^4)`$ gates and dirty
width, has $`O(\log(h+2))`$ T-depth, and is exactly identity when
$`H=0`$.

On the active branch, the first two operations leave the prefix of the
input marker in place and set its entire source suffix to zero. The
source then produces exactly the signed output of that marker, or the
root. Section 2 shows that the output query XORs precisely the existing
tag value, returning $`T=0`$. Its read-only address is $`P,X,H`$,
disjoint from its $`s`$ output bits. Thus $`Z`$ is exactly restored and
the final saved-predicate inverse returns the sole clean signal to zero.

On $`H=0`$, each completed operation is exactly identity on the full
occupied space, including $`T`$ and all dirty helpers. Hence the same
saved-predicate inverse is valid there as well. Linearity proves the
complete initialized-isometry equality, including superpositions of
active and inactive inputs and entangled dirty systems. This argument
uses no approximate work return or amplification.

The only global-address query appearances are two completed one-bit
XOR queries in the phase echo and one completed short-output decoder
query. Their multiplicity is constant. Packing and selected banks
involve only the local logical bits or the short tag; no sequence of
spectral comparisons or long controlled source network occurs.

## 5. Local height and global address ledger

Put

```math
a_{\mathrm{ext}}=n-m,\qquad h=m-r,\qquad
\overline Q=2^{a_{\mathrm{ext}}+h}=N2^{-r},\qquad K=m+2.
```

For a table with $`R`$ rows, $`q`$ output bits, and a disjoint arbitrary
dirty pool $`B`$, the retained exact sharp rectangular query has

```math
\begin{aligned}
T_Q&=O(q+\sqrt{Rq}+Rq/B),\\
G_Q&=O(Rq),\\
D_Q&=O\!\left(q+Rq/B^2
 +\log_2(\lceil\log_2R\rceil+2)\right),
\end{aligned}
```

provided $`B\ge16c_I(\lceil\log_2R\rceil+2)`$. The constant $`c_I`$
is the absolute helper constant of the
[nonuniform rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
There is no multiplicative polynomial in global address length.

| Completed operation | Rows or local work | Output width |
| --- | --- | --- |
| Signed source phase | $`O((h+1)\overline Q)`$ | one bit; two query appearances |
| Last-deviation decoder | $`O(\overline Q)`$ | $`s=O(\log(h+2))`$ |
| Last-one packing | $`O((h+2)^4)`$ local gates and dirty width | short tag |
| Selected marker clear and suffix Hadamards | polynomial local masks | private parallel outputs |
| Saved suffix predicate and its inverse | $`O(r+1)`$ exact Toffolis | one signal; two returned helpers |

The short-tag selected banks fit the same $`O(K^4)`$ local gate and
width envelope, by the explicit R80 bounds. Their depth is
$`O(\log(K))`$. Only the tag needs conditional initialization. Taking
$`B=\Theta(b)`$ gives the sufficient stage ledger

```math
\begin{aligned}
T_m&\le C\!\left(\sqrt{K\overline Q}
       +\frac{K\overline Q}{b}+K^4\right),\\
G_m&\le C(K\overline Q+K^4),\\
D_m&\le C\!\left(\log_2(n+2)+\log_2(K)
       +\frac{K\overline Q}{b^2}\right).
\end{aligned}
```

All constants are absolute. The decoder's $`s`$ output term contributes
only $`O(\log K)`$ depth. The phase-table and decoder row-output volumes
are each $`O(K\overline Q)`$. All query outputs and helpers return
before the next completed query; the common large query pool can be
reused. Reserve the separate $`O(K^4)`$ local dirty support to justify
parallel selected banks without overlapping non-Clifford supports.

## 6. Protected recursion and the simultaneous sums

Choose an absolute constant $`A\ge16`$, enlarged if necessary to fit
the short tag, and put

```math
r(m)=\lceil A\log_2(m+2)\rceil.
```

Fix an absolute cutoff $`m_*`$ large enough that $`r(m)\le m/2`$ for
$`m\gt m_*`$. Starting with $`m_0=n`$, apply section 4 and then recurse
on the suffix of height $`m_{i+1}=r(m_i)`$. The completed early bits
become part of the preserved external prefix. The residual rotations
are precisely the prescribed subtree rotations indexed by that prefix;
they still satisfy the same corner promise. No conditioning on root
amplitude and no removal of zero-weight subtrees is performed.

For $`m\le m_*`$, use the exact direct layers of section 7. There are
only constantly many terminal layers. Their prefix tables may still
have $`O(N)`$ rows: the global query is priced, rather than treated as a
constant-size subproblem. They contribute $`O(\sqrt N+N/b)`$ T gates,
$`O(N)`$ Clifford gates, and $`O(N/b^2+\log(n+2))`$ T-depth when the
shallow allocation applies.

Enlarge the fixed cutoff so that

```math
r(r(m))\le\log_2m,\qquad r(m)+2\le\sqrt{m+2}
```

above it. The first inequality gives $`O(1+\log^*n)`$ stages; the
second makes $`\log_2(m_i+2)`$ decrease geometrically. The heights also
decrease at least by a factor of two before the cutoff. With
$`K_i=m_i+2`$ these facts give

```math
\begin{aligned}
\sum_i\sqrt{K_i}\,2^{-r(m_i)/2}&=O(1),\\
\sum_i K_i\,2^{-r(m_i)}&=O(1),\\
\sum_i K_i^4&=O((n+2)^4),\\
\sum_i\log_2(K_i)&=O(\log_2(n+2)).
\end{aligned}
```

For example, the first two summands are bounded by
$`K_i^{(1-A)/2}`$ and $`K_i^{1-A}`$; reversing the height sequence bounds
their sums by convergent geometric series. The polynomial standalone
cost is absorbed by $`O(\sqrt N)`$ for T gates and $`O(N)`$ for Clifford
gates. These are uniform absolute asymptotic bounds, with the finitely
many small heights handled by enlarging the same constants.

Summing the stage ledger and the terminal layers gives

```math
T=O(\sqrt N+N/b),\qquad G=O(N),
```

```math
D_T=O\!\left(N/b^2+\log_2(n+2)(1+\log^*n)\right).
```

The global-prefix query term $`\log(n+2)`$ occurs at each recursive
stage; it has not been replaced by a local-height term. This is the
remaining iterated-logarithm factor.

For sufficiently large $`n`$, set $`b_c=\sqrt N/n`$. If $`b\ge b_c`$,
reserve a fixed fraction of $`b`$ for the large query pool and another
fraction for all local polynomial support. The exponential lower bound
on $`b`$ exceeds $`C(n+2)^4`$ and the shallow query's absolute helper
threshold for every stage and terminal layer. The largest source
address has at most $`n+O(\log(n+2))`$ bits. The conditional clean tag
occupies the protected logical suffix and does not consume an external
initialized wire. This proves the recursive branch with at most the
supplied dirty allocation. Section 7 covers smaller widths and the
finite remaining heights without changing the literal threshold.

## 7. An exact fallback at every literal width

The fallback also uses only the finite corner data. Write
$`\sigma_v=(-1)^{\kappa_v}`$, so that

```math
R_v=(-1)^{\kappa_v}(XZ)^{r_v}.
```

At logical position $`j`$, save in the sole clean signal the predicate
that all unvisited lower bits are zero. The node prefix, this suffix,
and the signal are preserved throughout the selected operation; the
target bit $`X_j`$ is excluded from every table address.

Use a dirty output bit to query $`Hr_v`$. A selected-$`Z`$ involution
echo and a separate selected-$`X`$ involution echo apply the chronological
word $`Z^{Hr_v},X^{Hr_v}`$ and return the dirty bit. Do not echo $`XZ`$
as if it were an involution. Add the exact Boolean phase echo for
$`H\kappa_v`$ and then reverse the actual suffix-predicate circuit.
This is exactly the desired guarded corner rotation, with the correct
literal sign. Each layer is identity on its inactive space, and all
work returns exactly.

There are constantly many completed one-output queries per layer. At
prefix length $`d=j-1`$, they have $`O(2^d)`$ rows including the activity
bit. Reserve one arbitrary output bit and reuse the rest of the dirty
allocation for the queries; borrow the two returned predicate helpers
between queries. At $`b\ge17(n+13)`$ this fits the routed query condition
$`B\ge4(d+3)`$ at every layer. It needs no precision core. Geometric
sums of table sizes and their square roots, together with the
$`O(n^2)`$ total predicate count, give

```math
T=O(\sqrt N+N/b),\qquad G=O(N).
```

If the absolute shallow threshold $`B\ge16c_I(n+2)`$ fits, the summed
depth is

```math
D_T=O(N/b^2+n\log_2(n+2)).
```

If it does not fit, use the routed schedule from R72. It replaces each
query's logarithmic address-depth term by its address length and gives
$`O(N/b^2+n^2)`$ depth. Failure of the shallow threshold implies
$`b=O(n+2)`$, with an absolute constant depending only on $`c_I`$.
The bounded sequence $`n^2(n+2)^2/2^n`$ then absorbs the $`n^2`$ term
into $`O(N/b^2)`$. Thus the preceding displayed fallback bound holds
at every stated literal allocation.

When $`b\lt b_c=\sqrt N/n`$, the term $`N/b^2`$ exceeds $`n^2`$ and
absorbs the fallback's serial depth. For the finitely many heights below
the cutoff needed by section 6, use the same exact fallback at every
eligible width and enlarge the absolute constant multiplying
$`\log(n+2)(1+\log^*n)`$. This completes one simultaneous exact theorem
at $`b\ge17(n+13)`$. No approximate generic fallback is substituted
into the exact equality.

## 8. Scope of the decoder improvement

The useful new interface is a constant-query exact tag eraser for
signed deterministic descendants, including markers below zero incoming
branches. It saves the repeated spectral comparisons needed when those
descendants overlap in the computational basis. The source table has a
single global address and a short tag, and every query is a native dirty
word with its actual inverse and all helpers returned.

For a positive interior bias, different marker columns can have
nonzero coefficients on the same computational leaf. A leaf-controlled
XOR cannot erase their distinct tags there: one output leaf cannot name
two different input depths. A general extension must create a basis in
which the depth supports separate, or supply another interference
mechanism. The exact corner theorem does not grant that basis change or
an arbitrary in-place tree permutation without a native resource proof.

The bounded fixture in
[`tests/test_orientation_depth.py`](../../tests/test_orientation_depth.py)
checks the packed-tag chronology against an independent prescribed frame,
including negative corner pairs, nonlinear orientations, root and marker
signs, physical dirty phase/mask echoes, exact last-deviation labels, and
occupied inactive tag inputs. The file's five tests also exercise the
affine construction R83. These finite checks support the literal algebra
and work-return claims; they do not establish the asymptotic resource
bounds, implement a scalable compiler, or test an unrestricted family.
