# Conditional suffix workspace and star residuals

[Open endpoint](OPEN_PROBLEM.md) · [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Grouped residuals](FAULT_TOLERANT_COMPILER.md#6-grouped-dictionaries-and-an-exactly-clean-coarse-frame)

[Section 10](#10-the-grouped-bounds-need-only-one-external-clean-qubit) extends
both count bounds to one external clean qubit. Sections 1–9 retain the
original two-clean construction and supply the group proof it reuses.

The logical suffix of a Hopf group is zero on the sector where that group
acts. It can supply temporary initialized workspace when every actual
subroutine is identity on the inactive sector. A decomposition of the
residual into columns indexed by ancestor depth makes the required private
work logarithmic in the group height. Several layers then share one
precision charge, and successive group heights can grow exponentially.

**Theorem.** Let $`n\geq1`$, $`N=2^n`$, $`0\lt\eta\leq1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. Define

```math
\ell_*(n)=1+\log_2^*(n+2),
```

where $`\log_2^*x`$ is the number of successive base-two logarithms
needed to obtain a value at most one. At every allocation

```math
a\geq2,\qquad b\geq L+n+7,
```

every prescribed complete real Hopf frame has a coherent Clifford+T
compiler satisfying

```math
\|VJ_a-J_a(W\otimes I_b)\|\leq\eta,
\qquad
T=O\bigl(N+L\ell_*(n)\bigr),\qquad G=O(NL).
```

Only two of the available clean qubits are used. At the selected endpoint
$`L=N`$, $`a=2`$, $`b=N+n+7`$, this gives
$`O(N\ell_*(n))=O(N\log_2^*N)`$ T gates, improving the preceding
$`O(N\log(n+2))`$ grouping bound. It does not prove $`O(N)`$ T count
or an optimal T-depth tradeoff.

The error includes initialized-work leakage and the return of arbitrary
dirty inputs, together with any references. There are no measurements,
resets, supplied resource states, or uncharged quantum oracles. The fixed
thresholds below need not give a practical saving for small systems.
Certified coefficient evaluation and table construction are separate
classical preprocessing costs. Angles must be effectively specified so
that the finite residual entries admit certified rational evaluation.
The theorem here concerns the prescribed real frame; no additional
complex-frame theorem is inferred from the use of complex coarse words.

## 1. An active suffix supplies workspace conditionally

Consider a group of depths $`[e-s,e)`$, and write

```math
r=n-e,\qquad p=e-s=n-r-s.
```

The first $`p`$ logical bits are a preserved prefix, the next $`s`$ bits are
the local system, and the final $`r`$ bits are the suffix. The ideal group is
identity unless the entire suffix is zero. On that active sector it is a
direct sum, indexed by the prefix, of complete $`s`$-level frames.

Compute $`h=[\mathrm{suffix}=0]`$ into the first external clean qubit using
the exact borrowed-MCX construction. One local system bit can be borrowed
and returned for this computation; its value need not be known. The cost is
$`O(r^2)`$ Toffolis. On the sector $`h=1`$, the $`r`$ suffix bits are now
available as zero-initialized private work. On $`h=0`$ they are arbitrary.

Every actual group subroutine below preserves $`h`$ and is identity on
$`h=0`$, on the whole input space of its other registers. Its active action
returns the suffix to zero up to the stated group error. Finally, reverse
the actual predicate computation. The ideal output has the same suffix
predicate as the input, so this last unitary returns $`h`$ to zero without
increasing the error. This reasoning does not assume that the actual suffix
is exactly zero after an approximate subroutine.

### Exact conditioning without a Clifford-count penalty

Conditioning a whole-word XOR query is done by adding $`h`$ to its address
and giving every inactive row the zero word. Its T-count increases by a
constant factor; its many leaf CNOTs are still Clifford gates. The same
rule applies when there is an additional identity/correction mode bit.

A controlled operator source has the form

```math
U_m\,C_F(X_0)\,U_m^\dagger,
```

where $`F`$ is the conjunction of a bounded number of flag literals. The
central multiple-controlled X costs a constant number of T gates, and
the unconditional source basis changes cancel exactly when $`F=0`$.

Other native gates needing a bounded number of controls also have constant
exact cost with a constant number of returned dirty helpers. For completeness,
the controlled-phase case does not require assuming that inactive suffix
scratch is clean. If $`z`$ is a dirty bit and $`F`$ is a Boolean conjunction
not depending on it, the chronological phase echo

```math
X_F^{(z)},\ T_z,\ X_F^{(z)},\ T_z^\dagger
```

returns $`z`$ and contributes the phase
$`\exp(i\pi(F-2zF)/4)`$. Correct its residual factor by the phase
$`i^{zF}`$. A second dirty echo using S implements that phase up to
$`(-1)^{wzF}`$, and this last correction is an exact multiple-controlled Z.
Only a bounded number of controls and returned dirty helpers occur. An
S-phase uses the last two steps of the same construction. A controlled H
uses $`H=R_y(\pi/8)ZR_y(-\pi/8)`$: the two fixed rotations have constant
Clifford+T cost up to opposite common scalars, which cancel in this
conjugation. They also cancel when the central controlled Z is inactive.
X, Z, and CNOT controls reduce to exact multiple-controlled X or Z.
All such helpers are disjoint from the controls and are returned before
the next operator-source call. A fixed number of core or selector wires
can serve this purpose; the ledger below also allows a separate fixed
helper reservation.

The same conditioning method applies to fixed coarse program words.
It applies to complete subroutines; their individual native gates need not
be identity on the inactive sector.

## 2. The residual support is a union of column forests

Consider an active group of s consecutive layers. Its local dimension is
$`M=2^s`$, and the unchanged external prefix is x. Its target $`W_g`$
is a direct sum of complete s-level frames. Let C be a coarse, exactly
unitary Clifford+T group with the same addressed layer structure and
exactly returned work. Its one-qubit words may be complex. The
[grouped-support lemma](FAULT_TOLERANT_COMPILER.md#6-grouped-dictionaries-and-an-exactly-clean-coarse-frame)
still applies, because the column supports depend on the addressed tree
structure rather than on real-valued entries.

For a nonzero local marker

```math
v=z\,1\,0^{s-j-1},\qquad 0\leq j\lt s,
```

its column is supported inside the aligned interval with prefix z and length
$`D_j=2^{s-j}`$. Two permitted marker intervals overlap exactly when their tree
nodes are comparable. Column zero has the whole local register as its
permitted support. Special angle choices may introduce additional zeros.
Therefore every
permitted entry of

```math
E_x=C_x^\dagger W_{g,x}-I
```

is diagonal, connects zero with another marker, or connects a nonzero
marker with one of its strict descendants.

The local indices in the interval of v are not all descendants of v.
The first index $`z0^{s-j}`$ is zero or an ancestor marker, and v itself
is the midpoint. Exclude both from the downward entries assigned to v.
All other indices in that interval are its strict descendants. Assign
zero-to-nonzero entries to a separate column type, and assign all
diagonal entries separately.

This is a disjoint partition. The number of strict downward pairs is

```math
(M-1)+\sum_{j=0}^{s-1}2^j(2^{s-j}-2)
=(s-1)M+1.
```

Including both orientations and the M diagonal entries gives
$`(2s-1)2^s+2`$, exactly the permitted-pair count of the grouped-support
lemma.

### Uniform column maps without endpoint registers

There are $`s+1`$ column types, indexed by
$`\nu\in\{\ast,0,\ldots,s-1\}`$. Define

```math
\Pi_\ast=|0^s\rangle\langle0^s|,
\qquad U_\ast=H^{\otimes s},\qquad D_\ast=2^s,
```

and, for $`0\leq j\lt s`$,

```math
\Pi_j=\sum_{z\in\{0,1\}^j}
|z1\,0^{s-j-1}\rangle\langle z1\,0^{s-j-1}|,
\qquad
U_j=(I_{2^j}\otimes H^{\otimes(s-j)})X_j,
\qquad D_j=2^{s-j}.
```

The marker-bit X acts before the Hadamards. Thus $`U_j\Pi_j`$ maps each
selected marker to the **positive** uniform superposition over its
interval, of amplitude $`1/\sqrt{D_j}`$. It preserves the interval's
prefix. The special type does the same for column zero on the whole
local register.

A diagonal filter following $`U_j\Pi_j`$ can independently assign every
entry in these columns. Its table is zero at the interval's first index
and at v, so it retains only strict descendants. The special filter is
zero at output zero. No quantum register stores an endpoint word: the
input predicate and the output's existing local bits specify the two
endpoints together with the depth label.

## 3. A small coefficient table and a streamed coarse circuit

Use four nonnegative phase components with
$`\omega\in\{1,-1,i,-i\}`$. There are four diagonal types, four
forward types for each $`\nu`$, and four reverse types for each
$`\nu`$. Pad their total to

```math
K=2^{\lceil\log_2(8s+12)\rceil}=O(s).
```

For a forward type, $`a_{x,u,\ell}`$ is the nonnegative component of
$`(E_x)_{u,v}`$ with phase $`\omega_\ell`$. For a reverse type, it is
the corresponding component of $`(E_x)_{v,u}`$. The endpoint v is
zero for the special type and the midpoint of u's aligned interval for
a depth-j type. The structural exclusions above receive literal zero
coefficients. Diagonal types use the components of $`(E_x)_{u,u}`$.

Choose the actual coarse words with one-qubit error at most

```math
\delta_c=\frac{1}{4sK\,2^{s/2}}.
```

Unitary telescoping over the s local layers gives

```math
\|E_x\|\leq s\delta_c
\leq\frac{1}{4K\,2^{s/2}}.
```

Each nonnegative split coefficient is bounded by the same quantity.
Define the real scalar-source tables by

```math
c_{x,u,\ell}=
\begin{cases}
K a_{x,u,\ell},&\text{diagonal type},\\
K\sqrt{D_\nu}\,a_{x,u,\ell},&\text{column or reverse type}.
\end{cases}
```

Every coefficient lies in $`[0,1/4]`$. Certified enclosures and clipping
at zero produce the phase components without an exact sign test.
Multiplication by $`\sqrt{D_\nu}`$ is classical coefficient evaluation,
not a supplied quantum amplitude. The rounded table uses the same
operator-source rule as the other compiler, including exact encodings
for structural zeros.

The phase-calibrated coarse words have length
$`w=O(\log(1/\delta_c))=O(s)`$. Storing such a whole word would defeat
the new clean-space bound. Instead, fix a common constant-alphabet
schedule and stream one symbol at a time. At each schedule position,
query its constant number of bits into a reusable constant-size clean
buffer, interpret the symbol, and erase those bits by the actual query
inverse. The query address consists of the unchanged logical prefix and
preceding local bits. The target bit is not part of that address, so
interpretation does not obstruct unloading. A reusable predicate flag
restricts the word to the appropriate local zero suffix; this flag is
computed and erased exactly. Identity padding uses an empty symbol.

Whole-word dirty lookup on a constant-size output costs $`O(2^d)`$
T gates at absolute depth d. Summing over word positions and depths in
the group gives

```math
T(C)=O(s2^e+\mathrm{poly}(n)),
\qquad G(C)=O(s2^e+\mathrm{poly}(n)).
```

The polynomial has a fixed degree; it accounts for literal controlled
word gates, suffix predicates, and returned scratch. The constant-size
buffer and predicate flag can reuse the residual block's private work.
Every coarse word is an actual native circuit with exactly returned
work, not an approximate-return operator-source implementation.
Its precise circuit is fixed before its residual is evaluated.

## 4. A common source implements all column types

Use the second external clean qubit as the scalar flag. Inside the active
suffix reserve a mode f, a $`\log_2K`$-bit term label, a distinct atom
rejection flag, and a constant number of temporary bits. The streamed
coarse interpreter reuses these temporary bits. There are fixed
constants A and B such that all initialized private work inside the
suffix occupies at most

```math
A\log_2(s+2)+B
```

wires. Endpoint registers and an s-bit program register are absent.

For every column type, construct an actual unitary $`\mathcal D_\nu`$ whose
zero-atom-flag block is $`U_\nu\Pi_\nu`$: toggle that flag by the
negation of the marker predicate, then apply $`U_\nu`$. Erase predicate
scratch before changing the local word. A depth selected by the term
label needs only reversible comparisons, marker predicates, and
conditionally applied Hadamards and X gates. They can be unrolled with
a polynomial number of gates in s and returned dirty helpers. This
unrolling does not call the precision source once per depth. Even the
all-zero local predicate can borrow a core or selector wire; it does not
require an additional logical bit outside its controls. The
predicates with many controls are charged as such; they are not treated
as bounded-control gates of constant cost.

Let $`\mathcal S_c`$ be the scalar-source unitary with its distinct
zero flag and m dirty core bits. Its accepted coefficient satisfies

```math
\widehat c_{x,u,\ell}I_{\rm dirty},
\qquad
|\widehat c_{x,u,\ell}-c_{x,u,\ell}|
\leq\epsilon_m:=\frac54\,2^{1-m}.
```

The query is addressed by the external prefix x, the **current local
logical word** u, and the term label. Its output is the arbitrary dirty
core itself. No initialized m-bit output word is present.

For a forward column type, define the real-base unitary

```math
T_\nu(c)=\mathcal S_c\mathcal D_\nu.
```

The scalar operation leaves the local word unchanged and acts on a
different rejection flag. With J initializing these two flags,

```math
J^\dagger T_\nu(c)J
=\mathrm{diag}(\widehat c_{x,u,\ell})U_\nu\Pi_\nu.
```

A reverse type uses the **actual inverse** $`T_\nu(c)^\dagger`$,
with the independently specified reverse coefficient table. Its block
is the adjoint of the displayed real-base block. The desired literal
phase $`\omega_\ell`$ is applied outside the forward or inverse base
word; it is not inadvertently conjugated for a reverse type. Diagonal
types use only $`\mathcal S_c`$ and their phase. Padding has a zero
coefficient and an identity base map.

The term label is preserved. Its direction, diagonal/column class, and
phase can be decoded into a fixed number of temporary bits by reversible
comparisons and erased afterward; the depth remains part of the label.
This allows the padded K-label encoding above without requiring separate
full-size label fields. In a common SELECT circuit, apply the selected $`\mathcal D_\nu`$ before the scalar operation on
forward branches, and the selected $`\mathcal D_\nu^\dagger`$ after the inverse
scalar operation on reverse branches. Diagonal branches need no atom
operation. **Two shared scalar blocks suffice**, one forward and one
inverse, independently of the number of depth labels. Their additional
source controls involve only h, f, the direction, and the scalar flag.
The large coefficient tables include the term label in their addresses;
no depth-controlled copy of each native source gate is introduced.

Every coefficient query and its inverse finishes before a local logical
word is changed by a subsequent atom operation. In the reverse word the
scalar query and its inverse occur first, followed by $`\mathcal D_\nu^\dagger`$.
Thus the actual inverse does not assume that a changed address still
contains its old value. Dirty selectors return after each complete query.

Let $`A_{x\ell}`$ denote the unrounded accepted atom, including its phase.
The scaling in Section 3 and the disjoint support partition give

```math
\frac1K\sum_{\ell=0}^{K-1}A_{x\ell}=E_x.
```

For a forward type, rounding changes its atom by at most

```math
\left\|\mathrm{diag}(\widehat c-c)
U_\nu\Pi_\nu\right\|\leq\epsilon_m.
```

The same bound holds for reverse types by taking adjoints, and for
diagonal types directly. In particular, neither the column size nor
the factor $`\sqrt{D_\nu}`$ amplifies the rounding error.

Prepare f and the term label uniformly, use identity for $`f=0`$ and
this SELECT for $`f=1`$, then undo the uniform preparations. All actual
subroutines are h-conditioned as in Section 1. Projecting the private
work and scalar flag to zero gives the exact accepted block

```math
B=\frac12\left(I+\frac1K\sum_\ell\widehat A_{x\ell}\right)
\otimes I_{\rm dirty}.
```

Therefore, uniformly over the coherent prefix direct sum,

```math
\|2B-(C^\dagger W_g\otimes I_{\rm dirty})\|
\leq\frac1K\sum_\ell
\|\widehat A_{x\ell}-A_{x\ell}\|\leq\epsilon_m.
```

The distinct scalar and atom flags justify their projected product.
No intermediate rejected component has been discarded.

## 5. Amplification and the complete group contract

Let $`R_h`$ be identity on $`h=0`$ and, on $`h=1`$, the reflection
$`I-2JJ^\dagger`$, with J now initializing the entire private suffix
and the external scalar flag. Testing unused active suffix bits is
harmless because they remain zero. The reflection costs $`O(r^2)`$
Toffolis with a returned borrowed local bit. Set

```math
\mathcal A_h=Z_h Q R_h Q^\dagger R_h Q,
```

where Q is the actual half-block circuit of Section 4. On $`h=0`$ all
factors are identity. On $`h=1`$, $`Z_h=-1`$ supplies the literal
leading minus sign required by normalization-two amplification.
The [complete-isometry amplification bound](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
gives an error at most

```math
4\epsilon_m=10\,2^{-m}
```

relative to $`C^\dagger W_g`$, including every rejected component and
all dirty references. Apply the actual coarse C afterward, then reverse
the original computation of h. On ideal private-zero columns C returns
its work exactly. Earlier leakage is propagated by actual unitaries,
so neither this application nor the predicate uncomputation increases
the group error. The result approximates the intended group on every
logical input and is exactly identity on the original inactive sector.

This argument permits the same two external clean qubits and dirty pool
to be reused between groups. It does not assume that actual intermediate
work has been reset or exactly returned.

## 6. Exponentially growing groups and the explicit workspace ledger

Choose fixed constants $`C_1\geq4`$ and $`r_0`$ sufficiently large for
the clean and dirty ledgers below. If $`n\leq r_0`$, use the existing
operator-source compiler for all layers. Otherwise reserve the deepest
$`r_0`$ layers for that compiler, and form the remaining groups backward
from the deep end. With r already reserved suffix wires, take

```math
s=\min\{n-r,\;2^{\lfloor r/C_1\rfloor}\},
\qquad e=n-r,
\qquad m=L+\lfloor r/4\rfloor+8,
```

then replace r by $`r+s`$. This specifies a partition, not the physical
execution order. Execute the groups in the prescribed shallow-to-deep
order, followed by the deepest reserved layers.

Choose $`C_1`$ large enough relative to the constants A and B in the
private-work bound. Since $`\log_2(s+2)\leq r/C_1+2`$, increasing the
fixed $`r_0`$ if necessary ensures that the active suffix contains all
private clean work.

A coefficient query has at most

```math
k=e+\log_2K+O(1)=n-r+O(\log(s+2))
```

address bits, including the identity/correction mode and h. The exact
dirty traversal needs k selectors. Other bounded dirty helpers contribute
a fixed number of wires. Increasing $`C_1`$ and $`r_0`$ once gives

```math
m+k+O(1)\leq L+n+7.
```

Indeed the left side is at most
$`L+n-3r/4+r/C_1+O(1)`$. The fixed additive eight in m and the
helper reservation are absorbed by the growing unused part of r.
No active private suffix bit is also counted as an arbitrary dirty
core or selector bit. Other queries and the streamed coarse interpreter
can reuse the same returned selector/helper pool. Larger dirty budgets
may be left unused.

Until the final capped group, r grows at least exponentially in r after
a fixed initial threshold. More explicitly, for sufficiently large r,
two uncapped updates dominate $`r\mapsto2^r`$: after the first update
$`r'\geq2^{r/C_1-1}`$, and then
$`2^{r'/C_1-1}\geq2^r`$. Consequently the number R of groups is

```math
R=O(1+\log_2^*(n+2))=O(\ell_*(n)).
```

### Error budget, including the reserved layers

For a deepest layer with r lower logical wires, use the old source width
$`m=L+r+5`$. Its explicit error bound, before conservative rounding,
is $`10\sqrt2\,2^{-m}`$. All reserved layers together therefore
contribute less than

```math
\frac{5\sqrt2}{8}\,2^{-L}.
```

A new group's error is at most
$`10\,2^{-L-\lfloor r/4\rfloor-8}`$. Distinct groups have distinct
integer r values. Choosing $`r_0`$ to be a sufficiently large multiple
of four gives

```math
\sum_{\rm groups}10\,2^{-L-\lfloor r/4\rfloor-8}
\leq\frac{80}{256}\,2^{-L-r_0/4}
\lt\left(1-\frac{5\sqrt2}{8}\right)2^{-L}.
```

The sum follows by grouping all integers $`r\geq r_0`$ into blocks of
four. Full-isometry telescoping, in the physical prescribed order, now
gives total error less than $`2^{-L}\leq\eta`$. Ideal preceding
circuits return their work; actual earlier errors are propagated
unitarily. This controls all initialized leakage and dirty/reference
return jointly.

## 7. T and Clifford counts

The common coefficient table has $`O(K2^e)=O(s2^e)`$ rows and m
output bits. Each SELECT uses only a constant number of scalar-source
calls. The inverse calls and three amplification appearances change
only constants. Including coarse streaming, predicates, and reflections,
a group satisfies

```math
T_g=O(s2^e+m+\mathrm{poly}(n)),
\qquad
G_g=O(s2^e m+m+\mathrm{poly}(n)).
```

The polynomials have fixed degrees. They include unrolling all selected
marker predicates and variable Hadamards, so no depth-selection work is
omitted.

Since $`s\leq2^{r/C_1}`$, $`e=n-r`$, and the r values are distinct,

```math
\sum_gs2^e
\leq N\sum_{r\geq r_0}2^{-(1-1/C_1)r}=O(N).
```

Also $`\sum_gm=O(LR+nR)=O(L\ell_*(n)+N)`$. All polynomial
routing and predicate overheads, even summed over at most n groups, are
$`O(N)`$ because any fixed polynomial in n is $`O(2^n)`$.
The reserved tail contributes $`O(N+L)`$, with constants depending on
the fixed $`r_0`$. These estimates establish
$`T=O(N+L\ell_*(n))`$ uniformly for every $`L\geq6`$.

For Clifford gates, the weighted table sum is essential:

```math
\begin{aligned}
\sum_g s2^e m
&\leq N\sum_{r\geq r_0}
2^{-(1-1/C_1)r}(L+r/4+8)\\
&=O(NL).
\end{aligned}
```

The remaining source costs $`O(L\ell_*(n)+N)`$ are also $`O(NL)`$;
the streamed coarse tables, polynomial overheads, and reserved tail obey
the same bound. Thus $`G=O(NL)`$ without assuming $`L\geq n`$.

## 8. Additional dirty banks give a width tradeoff

Put $`B_0=L+n+7`$. If $`b\geq2B_0`$, the same construction gives

```math
T=O\!\left(\sqrt{NL}+L\ell_*(n)+\frac{NL}{b}\right),
\qquad G=O(NL),\qquad a\geq2.
```

The proof still uses only two clean qubits. Reserve $`B_0`$ dirty wires
for the core, selectors, and returned helpers proved sufficient above.
The disjoint additional pool has size
$`K_{\rm bank}=b-B_0\geq b/2`$. In particular, every core word of
length $`m\leq B_0`$ fits in that bank pool.

For a coefficient table with $`Q_g=\Theta(s2^e)`$ padded rows and an
m-bit word, choose a power-of-two bank count $`\mu`$ within a constant
factor below

```math
\max\left\{1,\min\left(Q_g,\sqrt{Q_g/m},K_{\rm bank}/m\right)\right\}.
```

The [exact whole-word dirty SelectSwap query](OPERATOR_SOURCE_COMPILER.md#7-trading-additional-dirty-banks-for-lookup-cost)
then costs

```math
T_{\rm query}=O\!\left(\sqrt{Q_gm}+m+\frac{Q_gm}{b}\right),
\qquad G_{\rm query}=O(Q_gm).
```

Its $`\mu m`$ bank bits are disjoint from the core and selectors. The
query restores all banks and selectors exactly on arbitrary inputs.
The complete inactive table has zero rows, so the complete query is
identity there even though its individual loader and routing operations
need not be. This retains the inactive-sector contract without controlling
every leaf CNOT. The selected forward and inverse scalar blocks still
call only a constant number of coefficient queries per group.

The constant-size coarse symbol buffers also admit this banked query.
For a symbol table with $`S_j=\Theta(2^{p+j})`$ rows, its cost is
$`O(\sqrt{S_j}+S_j/b+1)`$ T gates. Streaming $`O(s)`$ positions at
each of the s local depths therefore costs

```math
T_{\rm coarse,g}
=O\!\left(s2^{e/2}+\frac{s2^e}{b}+\mathrm{poly}(n)\right),
\qquad
G_{\rm coarse,g}=O(s2^e+\mathrm{poly}(n)).
```

The same base selectors and returned helpers suffice. Banks can be reused
between a coarse symbol query and a coefficient query.

For $`C_1>2`$, the distinct-r weighted sums obey

```math
\begin{aligned}
\sum_g\sqrt{Q_gm_g}
&=O\!\left(\sqrt N\sum_{r\geq r_0}
2^{-(1-1/C_1)r/2}\sqrt{L+r+8}\right)
=O(\sqrt{NL}),\\
\sum_g Q_gm_g&=O(NL),\\
\sum_g s2^{e/2}
&\leq\sqrt N\sum_{r\geq r_0}2^{-(1/2-1/C_1)r}
=O(\sqrt N),\\
\sum_g s2^e&=O(N).
\end{aligned}
```

After the fixed initial threshold, uncapped suffix lengths at least double;
hence $`\sum_g r=O(n)`$ and
$`\sum_gm_g=O(L\ell_*(n)+n)`$. The fixed-degree polynomial overheads
are $`O(\sqrt N)`$ after summation and are absorbed by
$`O(\sqrt{NL})`$. The fixed number of reserved deepest layers use the
existing banked operator-source compiler, contributing
$`O(\sqrt{NL}+L+NL/b)`$ T gates and $`O(NL)`$ Clifford gates.
Combining these costs proves the banked bound.

In the **two-clean budget** $`a=2`$, total width
$`q=n+2+b=\Theta(b)`$ in this regime. The [retained worst-case lower
bound](FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage) is $`\Omega(\sqrt{NL}+L+NL/b)`$. Consequently the new upper
bound is matching whenever either

```math
L\ell_*(n)^2\leq N
\qquad\text{or}\qquad
b\leq N/\ell_*(n),
```

because $`L\ell_*(n)`$ is then absorbed by $`\sqrt{NL}`$ or
$`NL/b`$, respectively. These are sufficient matching regimes, not a
claim that they exhaust the possibilities. In particular, they do not
close the selected high-precision constant-clean endpoint.

## 9. Verification and scope

The [conditional-suffix tests](../tests/test_conditional_suffix_compiler.py)
exercise the support partition, native scalar-source blocks, separate
atom flags, actual inverse words, literal complex phases, amplification,
and the complete dirty-input contract in small instances. Finite checks
support the indexing and circuit identities. The asymptotic result rests
on the dimension-independent construction and resource ledger above;
small numerical matrices do not prove it by themselves.

This construction uses the existing operator source, dirty queries, and
normalization-two amplification. Its changes are the ancestor-column
residual representation, the streamed constant-size coarse-program
buffer, and the conditional suffix allocation. They leave only an
iterated-logarithmic multiplicity on the precision charge. The unrestricted
matching two-clean frontier, including whether the selected endpoint
admits $`O(N)`$ T gates, remains open. No optimal T-depth, optimal
Clifford-count, or literature-priority claim is established by this note.

## 10. The grouped bounds need only one external clean qubit

The [one-clean rotation primitive](ONE_CLEAN_COMPILER.md) removes the
second external clean qubit from the grouped construction. Under the same
complete-frame and certified-input assumptions, put $`B_0=L+n+7`$.
For every $`a\geq1`$ and $`b\geq B_0`$, there is a compiler with

```math
\|VJ_a-J_a(W\otimes I_b)\|\leq\eta,
\qquad
T=O\bigl(N+L\ell_*(n)\bigr),\qquad G=O(NL).
```

For $`b\geq2B_0`$, the same one-clean construction gives

```math
T=O\!\left(\sqrt{NL}+L\ell_*(n)+\frac{NL}{b}\right),
\qquad G=O(NL).
```

Only one external clean qubit is used; further clean qubits may be left
untouched. The changes are a private-work relocation and a different
compiler for the fixed deepest tail. The column forests, complex coarse
words, four phase components, and shared scalar blocks of Sections 2–5
remain unchanged.

### Group workspace and conditioning

Keep the sole external clean qubit as $`h=[\mathrm{suffix}=0]`$.
Move the scalar flag that Section 4 placed on the second external qubit
into the active suffix, alongside the distinct atom flag, mode, label,
and temporary work. The private initialized reservation is now

```math
w_g\leq A\log_2(s+2)+B+1.
```

Choose the same fixed $`C_1`$ sufficiently large and increase the fixed
multiple-of-four threshold $`r_0`$ until $`w_g\leq r`$ for every
group. This consumes one additional suffix bit, not an additional external
qubit. On $`h=1`$ every private bit is initialized; on $`h=0`$ the complete
conditioned subroutines remain exactly identity on arbitrary suffix and
dirty inputs. The scalar and atom flags remain distinct, so the accepted
product and actual reverse words in Section 4 are still valid.

The active reflection in Section 5 now tests the entire private suffix,
including its scalar flag. It still costs $`O(r^2)`$ Toffolis using a
returned borrowed local bit. The external $`Z_h`$ supplies the same literal
amplification sign. Applying the actual coarse circuit and reversing the
predicate computation gives the same group error $`10\,2^{-m_g}`$.
The argument includes suffix leakage and arbitrary dirty references.

The dirty core and query address are unchanged:

```math
m_g=L+\lfloor r/4\rfloor+8,
\qquad
k_g=n-r+\log_2K+O(1),
\qquad
m_g+k_g+O(1)\leq B_0.
```

The last inequality follows from the same
$`L+n-3r/4+r/C_1+O(1)`$ estimate after increasing $`r_0`$.
No suffix bit counted as private clean work is reused as a dirty core,
selector, or helper. The source controls have bounded size as before;
relocating their scalar control does not introduce a depth-dependent
source multiplicity. The coefficient tables still have $`O(s2^e)`$
rows and $`m_g`$ output bits. The coarse interpreter reuses returned
private work and the same dirty selector pool.

### The one-clean tail fits the exact dirty reservation

For a reserved layer at depth d, let $`r=n-d-1`$ and use the new
primitive with precision parameter and source-core width

```math
q_d=L+n-d+5=L+r+6,
\qquad m_d=q_d+1.
```

Its one clean signal flag is the sole external clean qubit. It uses d
dirty table selectors and one dedicated arbitrary helper, so its complete
dirty reservation is exactly

```math
m_d+d+1=L+n+7=B_0.
```

The suffix remains logical data in these tail stages. Conditional source
centers and the actual inverse words give exact inactive action; no suffix
bit is assumed initialized. Each stage has complete-isometry error
strictly below $`30\,2^{-q_d}`$. Consequently all reserved tail layers
contribute less than

```math
30\sum_{r\geq0}2^{-L-r-6}
=\frac{15}{16}\,2^{-L}.
```

If $`n\leq r_0`$, compile all layers this way, so no group or extra
private-work assumption is needed. Otherwise the group errors still sum
to at most

```math
\frac{80}{256}\,2^{-L-r_0/4}
\lt\frac1{16}\,2^{-L}
```

when the chosen threshold also satisfies $`r_0\geq12`$. Thus the total
error is strictly below $`2^{-L}\leq\eta`$. Full-isometry telescoping
uses the prescribed shallow-to-deep execution order. It propagates earlier
leakage by actual unitaries and does not assume an intermediate reset of
the external flag, private suffix, or dirty core.

### Resource sums and additional dirty banks

The group counts of Sections 7–8 are unchanged. The one-clean primitive
uses a constant number of source and coefficient-query calls per tail
stage. For the fixed number of reserved depths, its unbanked counts sum
to $`O(N+L)`$ T gates and $`O(NL)`$ Clifford gates, including the
charged suffix predicates. The finitely bounded case $`n\leq r_0`$
obeys the same estimates with constants depending only on $`r_0`$.

When $`b\geq2B_0`$, reserve the same disjoint base pool of $`B_0`$
dirty wires and use the remaining $`b-B_0\geq b/2`$ for word banks.
Both $`m_g`$ and $`m_d`$ fit in the base reservation and hence in the
additional bank pool. A tail table has $`S_d=2^d`$ rows; applying the
existing exact banked query to each of its constantly many scalar tables
costs

```math
O\!\left(\sqrt{S_dm_d}+m_d+\frac{S_dm_d}{b}\right)
```

T gates and $`O(S_dm_d)`$ Clifford gates. Banks and selectors return
exactly on arbitrary inputs. The two Pauli-mask tables are queried
sequentially and reuse those returned banks; their output is the dirty
core itself, not a fresh initialized program word. Actual inverse calls
use the same reservation. The fixed tail therefore contributes
$`O(\sqrt{NL}+L+NL/b)`$ T gates and $`O(NL)`$ Clifford gates.
Combining it with the unchanged group sums proves both displayed bounds.

For the fixed one-clean budget $`a=1`$, total width is
$`n+1+b=\Theta(b)`$ in the banked regime. The existing unrestricted
lower bound therefore makes the banked upper bound matching under the
same sufficient conditions $`L\ell_*(n)^2\leq N`$ or
$`b\leq N/\ell_*(n)`$. This extends the clean-workspace range; it
does not remove the iterated logarithm at $`L=N`$ or establish an optimal
T-depth bound.
