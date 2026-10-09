# Assembling the complete residual with two clean flags

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Research status](../../docs/OPEN_PROBLEM.md) · [Tree transport](ENDPOINT_TREE_TRANSPORT.md) · [Weighted blocks](WEIGHTED_TRANSPORT_BLOCK.md)

The diagonal and forward transport can be compiled together as one
affine tree map. Its dilation uses one signal flag. A second flag selects
between this dilation and the reverse weighted dilation, giving the
complete residual at normalization two. All native synthesis work is
borrowed, including when both clean flags are occupied.

This closes the composition step for the weighted-tree construction.
With $`N=2^n`$, $`n\ge1`$, $`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$, the resulting
complete-frame compiler has

```math
a=2,\qquad b\ge L+n+7,\qquad
T=O(N+nL),\qquad G=O(NL).
```

Its error includes initialized-flag leakage and approximate return of all
borrowed work, on arbitrary logical and dirty inputs with references.
At $`L=N`$ this is $`O(N\log N)`$ T gates. It does not improve the
retained $`O(N\ell_*(n))`$ complete-frame endpoint bound. The construction gives an explicit assembly within two clean flags; a separate
precision charge remains at every tree depth.

[Section 10](#10-collective-correction-and-its-precision-interfaces) separates
three requirements on a collective alternative: independent-stage error
budgets, orthogonal branch averaging, and correction of coarse work
leakage. The branch-rank restriction leaves the actual low-rank Hopf atoms
open, and the polar alignment formula is not a native compiler.

## 1. Use one algebraic target for both branches

Fix a dyadic constant $`0\lt\varepsilon_0\le1/64`$, independently
of n and eta, and put

```math
\alpha=4\varepsilon_0,\qquad s=2-\alpha.
```

Use the retained borrowed compiler to produce an actual native coarse
frame C at accuracy $`\varepsilon_0/2`$, with its per-depth error
allocation. Its external helpers return exactly, so its physical action
is $`C\otimes I_b`$ on all borrowed inputs. The
[coarse-frame argument](ENDPOINT_TREE_TRANSPORT.md#the-retained-borrowed-compiler-already-supplies-a-cheap-coarse-frame)
also bounds every subtree discrepancy by $`\varepsilon_0/2`$.
At the width used here its cost is $`T(C)=O(N)`$ and $`G(C)=O(N)`$;
at the selected endpoint the sharper retained estimate is
$`T(C)=O(\sqrt N)`$.

Replace each requested real local rotation by an algebraic real rotation
with operator error at most

```math
\delta=\min\left\{\frac{\varepsilon_0}{4n},
                         \frac{\eta}{64n}\right\}.
```

Certified rational stereographic coordinates on the unit circle, choosing
between its two charts, supply finite choices. Call the resulting
complete frame $`W'`$. The complete depth layers telescope to give
$`\|W'-W\|\le n\delta\le\eta/64`$, and every subtree of
$`W'`$ is within $`3\varepsilon_0/4`$ of its corresponding coarse
subtree. The native C and the algebraic $`W'`$ therefore have algebraic
overlap generators $`g_v,h_v,k_v,d_v`$ with

```math
|g_v|,|d_v|\le1,\qquad
\epsilon_v=1-|g_v|^2\le\varepsilon_0^2.
```

All zero tests, square roots, and completion choices below can be made
with algebraic zero tests and certified enclosures. Classical preprocessing
and its bit complexity remain separate from the quantum gate counts.
Completions need not vary continuously across zero defects: the circuit
approximates the chosen complete algebraic unitary.

Use the common heap indexing from the transport chapter and define

```math
\begin{aligned}
A&=\mathrm{diag}(g_1,\{d_v\}_{v=1}^{N-1}),\\
F&=\iota D_h\mathcal P_{W'},&
R&=\mathcal P_C^\dagger D_k\iota^\dagger,& S&=A+F.
\end{aligned}
```

The exact identity and retained weighted norm bounds are

```math
C^\dagger W'=S+R,\qquad
\|F\|,\|R\|\le2\varepsilon_0,\qquad
\|S\|\le1+2\varepsilon_0\lt s.
```

Using the same $`W'`$ for both branches makes their sum an exact unitary.
There is no independent rounding of residual coefficients whose sum would
then need a separate unitarity argument.

## 2. A scalar recursion dilates the affine forward map

Write $`u_0,u_1`$ for the columns of the local word $`U_v^{W'}`$.
For logical input amplitudes $`x_\ast,x_v`$, set

```math
a_1=x_\ast,\qquad
a_{2v+b}=u_0[b]a_v+u_1[b]x_v.
```

The desired outputs of S are $`g_1x_\ast`$ at the root and
$`h_va_v+d_vx_v`$ at each marker v. The added diagonal term must be
included in the dilation itself; this construction does not obtain S by
treating the earlier forward block as an opaque subroutine.

For the subtree rooted at v, write its stop-output map as
$`[t_v\ \ M_v]`$: the first column multiplies its incoming flow,
and the remaining columns multiply its own and descendant marker inputs.
Here $`\|t_v\|^2=\epsilon_v`$ by defect telescoping, while
$`M_v`$ is a marker diagonal plus restricted forward transport. Hence
$`\|M_v\|\le\mu:=1+2\varepsilon_0\lt s`$. Define

```math
\begin{aligned}
\rho_v
&=\sup_x\left(\|t_v+M_vx\|^2-s^2\|x\|^2\right)\\
&=t_v^\dagger(I-M_vM_v^\dagger/s^2)^{-1}t_v
\le\frac{\epsilon_v}{1-\mu^2/s^2}
\le\frac32\epsilon_v.
\end{aligned}
```

This is a strictly concave quadratic maximization. Completing the square
gives the inverse identity and proves finiteness before any recursive
calculation. As in the
[weighted-block recursion](WEIGHTED_TRANSPORT_BLOCK.md#2-a-scalar-recursion-supplies-the-missing-orthogonality),
these are classical Schur-complement operations, not calls to a quantum
matrix-inverse oracle.

Set leaf values to zero. The bottom-up scalar recurrence is

```math
\begin{gathered}
R_v=\mathrm{diag}(\rho_{2v},\rho_{2v+1}),\qquad
A_v=u_0^\dagger R_vu_0,\quad
B_v=u_1^\dagger R_vu_1,\quad
C_v=u_0^\dagger R_vu_1,\\
\beta_v=s^2-|d_v|^2-B_v,\qquad
E_v=\overline{h_v}d_v+C_v,\qquad
\rho_v=|h_v|^2+A_v+\frac{|E_v|^2}{\beta_v}.
\end{gathered}
```

Indeed, eliminate all descendant marker inputs in the variational
definition. The remaining marker amplitude z contributes

```math
|h_v+d_vz|^2+\sum_b\rho_{2v+b}|u_0[b]+u_1[b]z|^2-s^2|z|^2.
```

Maximizing this scalar quadratic gives the displayed recurrence. Uniformly,

```math
\beta_v\ge s^2-1-\tfrac32\varepsilon_0^2>\frac{11}{4},\qquad
|g_1|^2+\rho_1\le1+\tfrac32\varepsilon_0^2\lt s^2.
```

### Complete local unitaries

At a nonterminal node, order semantic output modes as stop, left
continuation, right continuation, rejection. For $`\rho_v>0`$,
prescribe the first two columns

```math
p_v=\frac1{\sqrt{\rho_v}}
\begin{pmatrix}
h_v\\ \sqrt{\rho_{2v}}u_0[0]\\
\sqrt{\rho_{2v+1}}u_0[1]\\ -\overline{E_v}/\sqrt{\beta_v}
\end{pmatrix},\qquad
q_v=\frac1s
\begin{pmatrix}
d_v\\ \sqrt{\rho_{2v}}u_1[0]\\
\sqrt{\rho_{2v+1}}u_1[1]\\ \sqrt{\beta_v}
\end{pmatrix}.
```

The recurrence gives unit norms and $`p_v^\dagger q_v=0`$.
Complete these columns by Gram–Schmidt on the standard basis, in order,
skipping zero projections and taking positive normalizing norms. Phase
the last added column to make the determinant one. If $`\rho_v=0`$,
first choose an algebraic unit vector orthogonal to the still-defined
$`q_v`$, then complete in the same way. The stop basis vector is not
generally orthogonal to $`q_v`$ because its first entry is now
$`d_v/s`$. The accepted incoming continuation amplitude at this node is
zero, so the first-column choice does not affect the required block.

At terminal nodes, $`\beta_v=s^2-|d_v|^2`$. Put
$`\omega_v=h_v/|h_v|`$ when h is nonzero and $`\omega_v=1`$
otherwise. Use the complete two-mode unitary

```math
U_v=\frac1s
\begin{pmatrix}
\omega_v\sqrt{\beta_v}&d_v\\
-\overline{d_v}&\overline{\omega_v}\sqrt{\beta_v}
\end{pmatrix}.
```

It has determinant one. For nonzero h the recurrence gives
$`\rho_v=|h_v|^2s^2/\beta_v`$, so its first stop coefficient is
$`h_v/\sqrt{\rho_v}`$, as required. Equivalently, this phases only
the rejection row of the two-column formula. It introduces no missing
relative scalar phase. If h is zero, the incoming continuation is zero.

## 3. Exactly one signal flag supplies the modes

Use heap logical labels $`0,1,\ldots,N-1`$, where 0 is the root,
and write a physical mode as $`(f,v)`$ with signal bit f. Initially
the signal-zero half contains all N logical input amplitudes.

On slots $`(0,0),(1,0),(1,1)`$, apply an SU(3) root factor with
first column

```math
\begin{pmatrix}
g_1/s\\
\sqrt{1-(|g_1|^2+\rho_1)/s^2}\\
\sqrt{\rho_1}/s
\end{pmatrix}.
```

The outputs are root stop, root rejection, and root continuation.
Complete this column algebraically and fix the determinant using a
completion column. Set $`c_v=(1,v)`$ for every internal node.
Process nodes from shallow to deep. At nonterminal v the physical slots
are $`c_v,(0,v),(1,2v),(1,2v+1)`$. Use the four-mode factor above
with output rows reordered as $`(0,3,1,2)`$: stop remains in
$`c_v`$, rejection goes to $`(0,v)`$, and the two continuations
go to their child slots. This row permutation is even. Terminal nodes
use $`c_v,(0,v)`$ with stop and rejection in that order.

Root mixing consumes two previously unused modes; the
$`N/2-1`$ nonterminal nodes consume two each. The total is exactly N,
including $`n=1`$. At completion, the root stop is at $`(0,0)`$
and each marker stop is at $`(1,v)`$. Flip the signal iff
$`v\ne0`$ to gather all stops into signal zero and all rejections
into signal one.

Inductively, the continuation at v has amplitude
$`\sqrt{\rho_v}a_v/s`$. Its stop is
$`(h_va_v+d_vx_v)/s`$, and its child continuations preserve this
invariant. At zero rho, nonnegativity of the recurrence gives
$`h_v=0`$ and $`\sqrt{R_v}u_0=0`$, so the same conclusion holds.
The complete unitary $`U_S`$ therefore obeys

```math
(\langle0|_f\otimes I)U_S(|0\rangle_f\otimes I)=S/s.
```

Every rejected-input column is specified by the actual factors. No signal
sector is reset, and the circuit is valid on arbitrary signal inputs.

For the reverse term, apply the earlier weighted-forward construction to
the tree C with weights $`\overline{k_v}`$ and normalization alpha.
Its accepted block is $`\iota D_{\overline k}\mathcal P_C/\alpha`$.
The actual inverse circuit defines $`U_R`$, whose accepted block is
$`R/\alpha`$. The row-defect identity for C and k supplies the same
norm and denominator bounds as the earlier h construction.

## 4. Native synthesis, including external controls

The fixed SU(4) elimination and Euler templates in
[weighted-block Section 5](WEIGHTED_TRANSPORT_BLOCK.md#5-a-native-implementation-by-depth-batching)
apply to the new factors. Their native compilation uses the
[borrowed-signal bound](WEIGHTED_TRANSPORT_BLOCK.md#6-a-faster-implementation-using-a-borrowed-signal),
which approximates the full unitary and has exact identity action whenever
its unchanged predicate is inactive. The occupied dilation signal is data;
a distinct arbitrary borrowed wire supplies the synthesis signal.

For clarity, the affine S routing changes are explicit. At nonterminal depth d,
write $`v=2^d+x`$. Partition the register into signal f,
$`n-d-2`$ leading logical bits, two mode bits A and B, and the d-bit
suffix x. Here A names a bit in this packing description, separately from
the diagonal operator above. An A-controlled right cyclic shift of the last
$`d+1`$ logical bits, using d Fredkins, packs the four slots into
three-bit labels $`101,001,110,111`$ with common suffix x.
The existing fixed permutation $`\pi=(0\ 4\ 7\ 3\ 6\ 2\ 5)`$,
fixing 1, sends these labels to $`000,001,010,011`$.
This works also at depth zero, where the cyclic shift is identity.
The packed leading bits form an unchanged zero predicate; two mode bits
carry the SU(4) factor, with x the address. Root SU(3) embeds in SU(4)
on the signal and low logical bit, fixing $`(0,1)`$, with the upper
$`n-1`$ logical bits zero. Terminal factors now use at most three
Euler templates rather than a single diagonal rotation.

Each four-mode factor needs at most 18 fixed Euler templates. Take an
absolute constant $`J\ge40`$ covering all templates in any depth,
including root mixing. To achieve full-operator error at most
$`\eta/64`$ for each complete dilation, choose a sufficiently large
absolute constant C and put

```math
q_d=L+n-d+C,\qquad K=C-4.
```

For $`d\ge K`$, fix K address bits as sector literals and use the
remaining $`k=d-K`$ free address bits. The dirty reservation is exactly

```math
(q_d+1)+k+1+1=L+n+7:
```

the source core, free-address selectors, one predicate helper, and one
borrowed amplification signal. The fixed number of sector circuits act
on disjoint invariant address sectors, with exact full-space identity on
the others, so their errors take a maximum. For $`d\lt K`$, direct
borrowed-sector echo compiles the bounded number of rows with
$`O(L+n)`$ native word length and $`O(n^2)`$ predicate cost. Increasing
C changes only construction constants and never the stated dirty threshold.

To select either branch using another logical bit c, add the literal
$`c=0`$ or $`c=1`$ to every template predicate. It introduces no
initialized helper, free address bit, or controlled-T synthesis assumption.
Packing and its inverse may stay unconditional around a selected template:
the inactive action is exactly $`P^\dagger IP=I`$ on the full space.
For the affine S branch, controlled gathering uses X on f under its
selector literal, followed by X on f under that literal and $`v=0`$.
This costs $`O(n^2)`$ exact Toffolis with an arbitrary helper.
The matching selector literal is $`c=0`$ for S and $`c=1`$ for R;
negative controls are implemented by paired X gates on c.

The reverse branch retains the earlier weighted-forward allocation before
taking its actual inverse. In that allocation the root continuation is
$`(0,0)`$ and $`(1,1)`$ is padding, so its depth-zero packing and
gathering remain those of weighted-block Section 5. Specifically, its
gathering flips f when $`v\ge2`$ and then cycles
$`(0,0)\mapsto(0,1)\mapsto(1,1)\mapsto(0,0)`$.
Condition each toggle in this retained permutation on the reverse-branch
selector literal, including both toggles implementing the three-cycle.
This also costs $`O(n^2)`$ with an arbitrary helper. The new affine
gathering must not be substituted for it. The controlled reverse branch
is the actual inverse of this entire controlled weighted-forward circuit,
including its own gathering and root mixing.

For either branch, physical marker reindexings before and after the entire
selected heap circuit may remain unconditional and cancel in the inactive
branch. They retain their charged $`O(n^4)`$ cost.

There are $`O(n^2)`$ packing gates and $`O(n^4)`$ charged physical
reindexing gates. The deep template sums give, for each selected branch,

```math
\begin{aligned}
T&=O\!\left(\sum_{d=0}^{n-1}(2^d+L+n-d+n^2)+n^4\right)
  =O(N+nL),\\
G&=O\!\left(\sum_{d=0}^{n-1}2^d(L+n-d)+nL+n^4\right)
  =O(NL).
\end{aligned}
```

Here $`n^4=O(2^n)`$ with an absolute constant. These are full-operator
approximations on the logical register, both occupied flags, and all dirty
work. Core and borrowed-signal return errors are included in the norm.
The borrowed pool is reused sequentially by the two branches and their
actual inverses.

## 5. Two-term selection and normalization-two amplification

The two initialized flags are the selector c and signal f. Choose the
determinant-one selector preparation P with first column
$`(\sqrt{s/2},\sqrt{\alpha/2})^{\mathsf T}`$. Define

```math
\begin{aligned}
\mathrm{SELECT}
&=|0\rangle\!\langle0|_c\otimes U_S
 +|1\rangle\!\langle1|_c\otimes U_R,\\
Q&=(P^\dagger\otimes I)\mathrm{SELECT}(P\otimes I),\\
J_2&=|00\rangle_{cf}\otimes I.
\end{aligned}
```

Because $`s+\alpha=2`$, its accepted block is exactly

```math
J_2^\dagger QJ_2
=\frac{s}{2}\frac{S}{s}
 +\frac{\alpha}{2}\frac{R}{\alpha}
=\frac{C^\dagger W'}2.
```

P is a fixed one-qubit rotation. Synthesize it to error at most
$`\eta/128`$ with $`O(L)`$ native gates and use its actual inverse.
A scalar phase in its synthesized word cancels between preparation and
unpreparation. Compile each selected dilation to full-operator error at
most $`\eta/64`$ as above. Both preserve c exactly and are exact
identity on their inactive c sectors. Their SELECT error is therefore at
most $`\eta/64`$, even on arbitrary shared dirty inputs, rather than
requiring initialized work between calls.

The two initialized flags have the following lifetimes. Initialization
occurs once, before the complete amplified circuit; no row below performs
a reset or assumes that work from the preceding row is clean.

| Stage | Selector c | Dilation signal f | Borrowed pool |
|---|---|---|---|
| Initial input | Zero | Zero | Arbitrary, including reference correlations |
| P in one Q call | Prepared superposition | Unchanged | Unchanged |
| Selected S dilation | Preserved | May be occupied | Source, selectors, helper, and synthesis signal |
| Selected R dilation | Preserved | May be occupied | Same pool, reused sequentially |
| Actual inverse of P | May remain occupied | Unchanged | Unchanged |
| Amplification reflections and actual inverse Q | Both flag sectors allowed | Both flag sectors allowed | Arbitrary input to every call; no intermediate reset |
| Final unconditional C | Unchanged | Unchanged | Coarse helpers return exactly |

The final complete-isometry bound below controls the return of both flags
to zero and of the entire borrowed pool to its input state.

Let $`\widetilde Q`$ be the actual native circuit, extending ideal Q
by identity on the dirty pool. Then

```math
\|\widetilde Q-Q\otimes I_b\|\le\frac{\eta}{32},\qquad
\left\|2J_2^\dagger\widetilde QJ_2
 -(C^\dagger W')\otimes I_b\right\|\le\frac{\eta}{16}.
```

Use the retained
[robust normalization-two lemma](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
with $`\mathcal R=I-2J_2J_2^\dagger`$ and the actual native word

```math
\widetilde{\mathcal A}
=-\widetilde Q\mathcal R\widetilde Q^\dagger
                    \mathcal R\widetilde Q.
```

The reflection tests only the two flags and is Clifford; the leading minus
sign is a literal Clifford scalar. Two forward calls and one actual inverse
are charged. The lemma controls the complete output, including rejected
flag sectors and borrowed-work disturbance:

```math
\left\|\widetilde{\mathcal A}J_2
 -J_2\bigl((C^\dagger W')\otimes I_b\bigr)\right\|
\le\frac{\eta}{4}.
```

Finally apply the actual unconditional C, which returns its borrowed
helpers exactly on every input. With
$`V=(I_{cf}\otimes C\otimes I_b)\widetilde{\mathcal A}`$,

```math
\|VJ_2-J_2(W\otimes I_b)\|
\le\frac{\eta}{4}+\frac{\eta}{64}
\lt\eta.
```

This is the complete-frame contract, so the retained QBP conclusions
apply. There are no measurements, resets, supplied catalysts, uncharged
quantum table oracles, or ignored reference systems.

## 6. What this assembly resolves

The two clean flags suffice for the whole residual, not only for one
weighted term. The diagonal is incorporated into the forward tree's
recursion, the reverse branch is an actual inverse, and amplification
reuses the same selector and signal. At the present per-depth precision
cost, this provides $`O(N+nL)`$ T gates and $`O(NL)`$ Clifford gates.

The abstract assembly needs only accepted-block estimates for actual
unitary selected branches. It does not require approximation to the
particular Gram–Schmidt completions chosen above. Suppose the branches'
signal-zero compressions, including all arbitrary dirty inputs, differ
from $`(S/s)\otimes I_b`$ and $`(R/\alpha)\otimes I_b`$ by at
most $`\delta_S,\delta_R`$. With selector-preparation error
$`\delta_P`$, the doubled accepted-block error against
$`(C^\dagger W')\otimes I_b`$ is at most

```math
\zeta\le2\max\{\delta_S,\delta_R\}+4\delta_P.
```

The preparation estimate uses the unitarity of the actual selected words;
the branch estimate is their weighted accepted-block sum. The same robust
amplification lemma includes all rejected-output and work-return errors.
Comparison directly with W adds at most $`n\delta`$ to zeta instead.
The native construction in Section 4 supplies stronger full-operator
estimates, but future branch implementations may choose different
rejected-space completions.

Thus, if improved native selected implementations have costs
$`T_S,T_R`$ at the required constant-fraction accepted-block accuracy
and the same clean/dirty allocation, the total cost is

```math
T=O(T_S+T_R+N+L).
```

The $`O(N)`$ allowance covers coarse compilation and explicit routing;
preparation costs $`O(L)`$, and selection/amplification use a constant
number of branch calls. This is a composition statement with specified
controls and work lifetimes. A faster implementation of the old isolated
F block alone does not automatically price the new affine S block, nor
does an uncontrolled black-box gate count automatically price SELECT.
An endpoint improvement now requires an improved implementation of these
structured selected branches with accepted-block error controlled on all
logical and dirty inputs. The linear-T endpoint and optimal T-depth remain
open.

## 7. A bounded audit of fusion across tree depths

The affine and reverse-adjoint branches share an injection-and-stopping
tree interface. This section tests whether merging adjacent depths removes
their repeated precision charge. The complete local unitaries do close under
merging, including arbitrary inputs in the modes used for rejection and
continuation. The resulting hierarchy retains linear classical storage.
Those facts do not yet supply a cheaper Clifford+T implementation.

### A parent and its two children give an exact ten-mode word

Use the physical row order stop, rejection, left continuation, right
continuation for a nonterminal local unitary. Write its first two rows as
$`H_v`$ and its continuation rows as $`\ell_0,\ell_1`$, so

```math
U_v=\begin{pmatrix}H_v\\ \ell_0\\ \ell_1\end{pmatrix},
\qquad H_v\in\mathbb C^{2\times4}.
```

Its four input modes are incoming continuation, marker, and the two modes
allocated to its children. For each child b, write its complete four-mode
unitary as $`U_b=[\tau_b\ \ Z_b]`$, separating the incoming
continuation column from its three other input columns. Applying the parent
and then both children gives the exact matrix

```math
\mathcal M_v=
\begin{pmatrix}
H_v&0&0\\
\tau_0\ell_0&Z_0&0\\
\tau_1\ell_1&0&Z_1
\end{pmatrix}.
```

The column groups have sizes 4, 3, 3; the row groups have sizes 2, 4, 4.
Thus this is a ten-mode unitary, with all ten input amplitudes arbitrary.
It is the actual product of the three local factors, with fixed input and
output orderings. No rejected sector has been projected out. In particular,
the $`Z_b`$ columns cannot be omitted when this word occurs inside
selection or amplification. The modes here are basis states of the existing
signal/logical register, not ten initialized qubits.
The displayed input and output orderings differ, so its determinant need
not be one even when all three physical local factors have determinant one.
Restoring the common physical basis includes the corresponding permutation;
no scalar phase is discarded.

### The hierarchy closes with one coupling column per edge

For a closed subtree containing $`m_v`$ internal nodes, the same
invariant is

```math
V_v=[\tau_v\ \ Z_v]\in U(2m_v).
```

Its first input is the incoming continuation; the other
$`2m_v-1`$ inputs are arbitrary. Its outputs are all stops and
rejections in that subtree. A terminal subtree uses its two-mode factor.
For a nonterminal node, substitute the complete child matrices
$`V_b=[\tau_b\ \ Z_b]`$ into the displayed merge formula.
The result has dimension $`2(1+m_0+m_1)=2m_v`$; its first column
is the new $`\tau_v`$.

To check unitarity directly, each child satisfies
$`\tau_b^\dagger\tau_b=1`$,
$`\tau_b^\dagger Z_b=0`$, and $`Z_b^\dagger Z_b=I`$.
The parent block therefore contributes $`U_v^\dagger U_v=I_4`$
to the four parent-input columns, while cross terms with child-private
columns vanish. This proves the invariant on the whole input space.
Root mixing and final gathering remain the explicit outer operations
already charged above.

Each child coupling $`\tau_b\ell_b`$ has rank at most one.
Keeping the constant-size local matrices and child pointers therefore
uses $`O(N)`$ classical generators. There is no need to materialize
every expanding $`\tau_v,Z_v`$ array. Rank one on an individual
tree edge does not mean rank one across a whole depth cut. On initialized
accepted inputs, the continuations at depth d have matrix

```math
\mathcal C_d
=\mathrm{diag}\!\left(\frac{\sqrt{\rho_v}}{\sigma}\right)
 W_d,\qquad v\text{ at depth }d,
```

where $`\sigma=s`$ or alpha and $`W_d`$ is the complete
$`2^d`$-dimensional truncated propagation frame on the root and
strictly higher marker inputs. Its rank equals the number of positive
$`\rho_v`$, generically $`2^d`$. On arbitrary full input modes,
the frontier rows are rows of a unitary and have rank $`2^d`$.
All of these modes already fit in the existing register. This rank count
does not require extra clean qubits and is not a circuit lower bound.

For comparison, eagerly expanding a generic h-level affine accepted map
into entries gives

```math
1+\sum_{j=0}^{h-1}(j+2)2^j=h2^h+1
```

nonzero entries. A depth-j marker row contains its root contribution, its
j strict ancestors, and its own diagonal. With zero root and diagonal,
as in the reverse-adjoint map, the count is
$`(h-1)2^h+1`$. This diagnoses the extra factor of h in an eager
entrywise representation. It does not apply to the linear-size hierarchy
just proved or establish how many coherently queried rows another compiler
must use.
The counts are attained within both branch families: take all coarse words
to be $`R_y(\pi/4)`$ and all target words to be
$`R_y(\pi/4+\vartheta)`$ for a sufficiently small positive perturbation
whose rotation matrix has algebraic entries. Identical child subtrees give nonzero h, k, and d, and both
propagation columns have nonzero entries. Choosing
$`\vartheta\le\varepsilon_0 2^{-n-4}`$ also meets the retained
per-depth coarse allowances.

### Fixed permutations do not make the completed band an unchanged-address table

An open band of g nonterminal levels has

```math
D_g=2(2^g-1)+2^g=3\cdot2^g-2
```

output modes: two stops/rejections per processed node and a continuation
for every frontier node. This count agrees with the ten-mode case at
$`g=2`$. There are legal affine inputs for which the band's incoming
continuation column is nonzero on all $`D_g`$ outputs.

For example, choose native $`C=I`$ and the same positive real algebraic
rotation $`R_y(\theta)`$ at every target node, with
$`\theta\le\varepsilon_0 2^{-n-4}`$. Positive rational
stereographic coordinates supply such rotations. This choice meets even
the retained per-depth coarse error allowances. Equal child subtrees give
equal child rho values and $`C_v=u_0^\dagger R_vu_1=0`$ in the
scalar recurrence. Here $`h_v,d_v,E_v`$ are positive, so each
nonterminal incoming continuation column has four nonzero entries.
Every path through the open band consequently contributes a nonzero
amplitude to its own output, without competing paths to cancel it.

The bipartite support graph of the full band matrix is therefore connected:
this one column meets every output row, and unitarity makes every other
input column meet some row. Independent input and output permutations
preserve connected components. They cannot turn this fixed completed word
into a direct sum of smaller blocks. In particular, for fixed r mode bits,
permutation-only repacking cannot represent all band heights as one
unchanged-address table with blocks of size at most $`2^r`$.

This statement concerns exact permutation-only repacking of the specified
complete word. It supplies no approximation lower bound, excludes no
nonpermutation change of basis, and does not constrain a different unitary
completion with the same accepted block.

### A strict external norm margin does not damp the internal continuation

There is also a concrete obstruction to replacing this tree transport by
a bounded-depth local polynomial merely because its accepted map is a
strict contraction. Fix a rational $`0\lt z\le\varepsilon_0/8`$,
independently of n, and set

```math
\theta=2\arctan z,\qquad
c=\frac{1-z^2}{1+z^2},\qquad
t=\frac{2z}{1+z^2}.
```

Take every coarse local word to be $`R_y(\pi/4)`$. Let target
$`W'`$ use the same words except at terminal internal nodes, where
the word is $`R_y(\pi/4+\theta)`$. The coarse frame is an actual
native circuit: $`R_y(\pi/4)=HZ`$, so each depth uses selected H and
Z gates. The retained exact identity $`H=VZV^\dagger`$ with
$`V=SHTHS^\dagger`$ reduces each selected H to a selected Z; the
borrowed multi-controlled-X construction supplies these selected Z gates
with $`O(n^2)`$ gates per depth. The total is $`O(n^3)`$ native
gates. No claim that this entire coarse frame is Clifford is needed.
The target is algebraic, and it can itself be the requested W. Moreover,

```math
W'=(I\otimes R_y(\theta))C,\qquad
\|W'-C\|=2\sin(\theta/2)\le\theta\le\varepsilon_0/4.
```

The same bound holds for every subtree. Terminal overlaps have
$`g_v=d_v=c`$, $`h_v=t`$, and $`k_v=-t`$.
At all higher nodes the overlap matrix is $`cI`$, so
$`g_v=d_v=c`$ and $`h_v=k_v=0`$.

Let H be the bottom-depth row block of
$`\mathcal P_C=\mathcal P_{W'}`$, and let
$`T=\iota_{\rm bottom}H`$. This is a partial isometry from the
root and strictly higher marker inputs to the bottom marker outputs.
Its initial and final projections partition the logical space, and
$`T^2=0`$. The residual pieces are

```math
F=tT,\qquad R=-tT^\dagger,\qquad S=cI+tT,
\qquad C^\dagger W'=cI+t(T-T^\dagger).
```

In this family $`J=T-T^\dagger`$ satisfies $`J^\dagger=-J`$ and
$`J^2=-I`$, so the complete residual is $`e^{\theta J}`$. This
special relation uses both forward and reverse pieces.

Thus $`\|S/s\|\lt0.54`$ and
$`\|R^\dagger/\alpha\|=t/\alpha\le1/16`$: both external
norm margins are independent of height. Nevertheless, the affine recursion
has the same positive message at every internal node,

```math
\rho_v=\frac{t^2s^2}{s^2-c^2}.
```

At the bottom this is the terminal formula. Above it, equal child messages,
orthogonal local columns, and zero h give $`A_v=B_v=\rho_v`$
and $`C_v=E_v=0`$, proving the claim inductively.
The reverse-adjoint weighted recursion likewise has $`\rho_v=t^2`$
at every node. In either case, the nonterminal incoming continuation column,
in semantic stop/left/right/rejection order, is exactly

```math
p_v=(0,u_0[0],u_0[1],0)^{\mathsf T}.
```

There is no stop or rejection attenuation before the bottom. In the
transport notation of the earlier chapter, the rescaled internal shift is

```math
\mathcal B_\rho
=D_{\sqrt\rho}\mathcal B D_{\sqrt\rho}^{-1}
=\mathcal B,\qquad
\|\mathcal B^{n-1}|1\rangle\|=1.
```

For $`n\ge2`$ it has norm one and nilpotency index n. Nilpotence
does not give a height-independent geometric decay of its successive
powers.

Consequently a replacement of the precise form

```math
\iota D_h\,p_D(\mathcal B_\rho)\mathcal G,
\qquad \deg p_D=D\lt n-1,
```

has zero root-to-bottom output, whereas the true forward term sends that
input to a bottom vector of norm t. Its normalized error is at least
$`t/s`$ for the affine branch and $`t/\alpha`$ for the
reverse-adjoint branch. Adding the affine diagonal cannot restore this
missing bottom component. The same support argument applies to sums of
words containing at most D factors from
$`\mathcal B_\rho,\mathcal B_\rho^\dagger`$, interleaved only
with operators block diagonal in tree depth: each shift factor changes
depth by at most one.

This is a restriction on the stated local-depth word class. Global logical
gates, tree shortcuts, another completion, or a jointly synthesized native
word are outside it. Indeed, the witness itself has the simple global
one-angle formula through C displayed above. It is not a hard endpoint
instance and proves no additive $`nL`$ T-count lower bound.

The bounded audit therefore retains the closed linear-size hierarchy and
rejects two specific shortcuts: exact permutation-only compression to a
fixed-size unchanged-address block, and constant-degree local truncation
justified only by the external norm margin. Keeping the original scatterers
still gives $`T=O(N+nL)`$, $`G=O(NL)`$. No T-count bound has
improved in this audit; the linear-T endpoint remains open.

## 8. A coupled completion and its native cost

The coupled unitary residual admits a simpler complete boundary contract
than separate forward and reverse contractions. This elementary completion
is useful for testing a recursive word: it has one signal flag, constant
normalization, and a continuous zero-defect limit. It does not by itself
give a cheaper native compiler.

For any unitary V, define

```math
\mathcal D(V)=
\begin{pmatrix}
V/2&-\sqrt3 I/2\\
\sqrt3 I/2&V^\dagger/2
\end{pmatrix}.
```

Direct multiplication proves unitarity on both signal sectors. Its
accepted block is V/2, and for any two unitaries V and U,

```math
\|\mathcal D(V)-\mathcal D(U)\|=\frac12\|V-U\|.
```

In particular the complete word, including its rejected columns, approaches
the fixed $`\mathcal D(I)`$ as the residual approaches identity. Tensor
every block with dirty identity to include arbitrary borrowed inputs and
references. No clean history or additional logical padding is required.

### Exact branching on the whole input space

Let $`\mathcal R_v=C_v^\dagger W_v`$ be the complete residual of
a subtree. Embed its local root words $`U_v^C,U_v^W`$ on the two
child-root modes, acting as identity on all other subtree modes. In this
common basis the complete-frame recursion is

```math
K_v=\mathcal R_{2v}\oplus\mathcal R_{2v+1},\qquad
\mathcal R_v=(U_v^C)^\dagger K_vU_v^W.
```

This uses all child columns, not only their root overlaps. Put
$`E_v=(U_v^C)^\dagger U_v^W`$. Then the exact signal/logical word is

```math
\begin{aligned}
\mathcal D(\mathcal R_v)
={}&\mathrm{diag}(I,E_v^\dagger)
 (I_f\otimes(U_v^C)^\dagger)\,
 \mathcal D(K_v)\\
&\quad\cdot (I_f\otimes U_v^C)
 \mathrm{diag}(E_v,I).
\end{aligned}
```

Here $`\mathcal D(K_v)`$ is the coherent direct sum of the two
complete child words on the same signal f. To verify the equation, write
$`A_v=(U_v^C)^\dagger K_vU_v^C`$. The four resulting blocks are
$`A_vE_v/2`$, $`-\sqrt3 I/2`$, $`\sqrt3 E_v^\dagger E_v/2`$,
and $`E_v^\dagger A_v^\dagger/2`$. Thus every occupied rejected
input is included and the normalization remains two at every height.
Leaf residuals are identity.

There is also a stable approximate boundary. Suppose the actual whole
child words have errors at most $`\delta_0,\delta_1`$, and replace
$`E_v`$ by an actual unitary $`\widehat E_v`$, using its actual
inverse in the other wrapper. With the coarse wrappers exact, the merged
whole-word error is at most

```math
\max\{\delta_0,\delta_1\}
+\frac12\|\widehat E_v-E_v\|.
```

First replace the children by their exact direct sum, which costs the
maximum child error under unitary wrappers. The remaining word is exactly
$`\mathcal D(A_v\widehat E_v)`$; apply the norm identity above.
Independent unrelated approximations to the two wrappers do not enjoy
this identity. The same argument applies on a common arbitrary dirty
space if each error includes return of that space.

### One signal preparation does not remove target precision

Unfolding this recursion gives the exact word

```math
\mathcal D(C^\dagger W)
=\mathrm{diag}(C^\dagger,W^\dagger)
 \mathcal D(I_N)
 \mathrm{diag}(W,C).
```

The middle factor is one signal rotation. Approximating it once to error
$`\epsilon`$ changes the whole word by at most $`\epsilon`$.
However, both W and its actual inverse still appear in the wrappers.
They cannot be supplied as target-frame oracles. At the local level the
same cost is the fine-precision $`E_v`$ in every parent wrapper.
Native substitution preserves the algebraic contract, but each emitted
target program, its controls, and its inverse must be charged. Using the
existing layerwise primitives retains their repeated precision cost.
Invoking the grouped compiler does not establish an improvement over its
known bound; a conditional invocation still needs a valid work allocation
and explicit literal-control accounting. Counting only the one middle
rotation proves no new T-count estimate.

This answers the boundary part of the branching question. The missing
part is a native implementation that avoids reintroducing the complete
target program, not an additional proof of normalization or mode closure.

### Anchoring at the zero-defect completion leaves a mixed term

A more economical proposed word would reuse $`\mathcal D(I)`$ as
a fixed connector between independently encoded child and parent
residuals. For arbitrary unitaries A and B, its accepted block obeys

```math
\left[\mathcal D(A)\mathcal D(I)^\dagger\mathcal D(B)\right]_{00}
-\frac{AB}{2}
=-\frac38(A-I)(B-I).
```

This follows by multiplying the three complete two-by-two block matrices;
no rejected block is discarded. In the tree merge take
$`A=(U_v^C)^\dagger K_vU_v^C`$ and $`B=E_v`$. Their product is
exactly the coupled residual $`\mathcal R_v`$, so the coupled
unitarity relations do not in general cancel the displayed mixed term.
For a fixed nonzero parent and independently changed children its norm
can stay positive as synthesis precision increases.

If A and B change disjoint invariant subspaces, the mixed term vanishes
and the full anchored word equals $`\mathcal D(AB)`$. This special
case does not cover generic comparable tree updates. Small coarse-frame
error can make the mixed term small, but a fixed coarse tolerance does
not make it arbitrarily accurate. The identity concerns this anchored
word; it is not a lower bound on all choices of completion or fusion.

### The full missing repair acts on at most four modes

The mixed term also gives a bounded next target. In a parent merge B
acts only on the two parent-root modes. Put $`X=A-I`$, $`Y=B-I`$,
$`c=1/2`$, $`t=\sqrt3/2`$, and
$`F=\mathcal D(A)\mathcal D(I)^\dagger\mathcal D(B)`$.
The complete error, including rejected ports, factors as

```math
F-\mathcal D(AB)
=c\begin{pmatrix}tX\\cX^\dagger\end{pmatrix}
  \begin{pmatrix}-tY&cY^\dagger\end{pmatrix}
-c\begin{pmatrix}0&0\\0&Y^\dagger X^\dagger\end{pmatrix}.
```

Because B is unitary, Y is normal and
$`\mathrm{range}(Y)=\mathrm{range}(Y^\dagger)`$.
If $`r=\mathrm{rank}(B-I)`$, the first summand has rank at
most r and the second at most r. Consequently the accepted error has
rank at most two and the complete error has rank at most four at a
tree parent. No signal sector has been discarded in the latter bound.
These ranks are on signal/logical space. After tensoring with dirty
identity, the Hilbert-space rank multiplies by the dirty dimension;
the repair acts trivially on that factor.

The exact left repair

```math
\mathcal K=\mathcal D(AB)F^\dagger
```

is unitary and satisfies $`\mathrm{rank}(\mathcal K-I)\le4`$.
A unitary minus identity is normal, so the repair is identity on the
orthogonal complement of its at-most-four-dimensional moving subspace.
Moreover, $`\det\mathcal D(V)=1`$: factor it as
$`\mathrm{diag}(V,I)\mathcal D(I)\mathrm{diag}(I,V^\dagger)`$.
Thus F and the repair also have determinant one. The restriction to its
moving subspace is special unitary, with no discarded scalar phase.
If $`|r_0\rangle,|r_1\rangle`$ are the two parent-root basis modes,
that moving subspace is contained in the span of

```math
\begin{pmatrix}t(A-I)|r_b\rangle\\c(A^\dagger-I)|r_b\rangle\end{pmatrix},
\qquad
\begin{pmatrix}0\\|r_b\rangle\end{pmatrix},
\qquad b\in\{0,1\}.
```

These are four generally transported vectors, not four computational
basis labels or four free clean qubits. The definition of the repair is
an existence identity; implementing $`\mathcal D(AB)`$ by its target
wrappers would repeat the unresolved problem. A useful native repair
must make both A and its adjoint's boundary vectors accessible with a
charged circuit and maintain their return promises. It must also prove
that recursive child calls and precision work do not multiply at every
parent. The rank bound alone supplies neither operation nor resource
bound. The next section tests this full four-mode repair, including its
rejected-space action and its recursive call ledger.

The [coupled-merge checks](../../tests/test_coupled_residual_merge.py)
audit the same recursion at a fork and height three, using real targets
and actual complex native coarse words, all four signal blocks, actual
inverses, and reference-correlated spectators. They also test the mixed
term and the stability estimate. These are bounded operator checks, not
an elementary endpoint emitter. The
[shared-source fork audit](SOURCE_REUSE_LIMITS.md#5-a-shared-conjugator-does-not-close-a-branching-fork)
separately tests a fully native source-cancellation candidate.

## 9. A repair word without an ill-conditioned transported basis

The four-mode repair has an explicit commutator implementation using
the whole child word and two conditional parent words. It does not need
to orthonormalize small difference vectors or invert an overlap gap.
Its actual inverse calls cancel against the anchored word. After that
cancellation, however, the circuit is exactly the fine-precision parent
wrapper construction above. This resolves the access and call-count
questions for this particular repair, without improving the endpoint bound.

### Whole-port commutator and support

Retain the parent unitaries A and B from Section 8, and define

```math
\begin{aligned}
J&=\mathcal D(I),& U&=\mathcal D(A)J^\dagger,\\
R_B&=\mathrm{diag}(I,B^\dagger),& S_B&=\mathrm{diag}(B,I).
\end{aligned}
```

The local completion satisfies $`\mathcal D(B)=R_BJS_B`$.
Consequently the exact left repair is the complete unitary word

```math
\begin{aligned}
\mathcal K&=R_BU R_B^\dagger U^\dagger\\
&=R_B\mathcal D(A)J^\dagger R_B^\dagger J\mathcal D(A)^\dagger.
\end{aligned}
```

Indeed, for the anchored word $`F=U\mathcal D(B)`$,

```math
\begin{aligned}
\mathcal KF
&=R_BU R_B^\dagger U^\dagger U R_BJS_B\\
&=R_B\mathcal D(A)S_B
=\mathcal D(AB).
\end{aligned}
```

All cancellations are between complete unitaries and their actual
inverses; they hold on occupied rejected ports as well as accepted ones.
The repair uses the existing signal f and the same logical register.
It introduces no initialized modes, measurements, or resets.

Let P project onto the two lower-signal parent-root modes. The unitary
$`R_B`$ acts as identity outside P. The commutator therefore acts as
identity outside

```math
\mathcal M=\mathrm{range}(P)+U\,\mathrm{range}(P),
\qquad \dim\mathcal M\le4.
```

Both P and its transported image have orthonormal root frames. With
$`c=1/2`$ and $`t=\sqrt3/2`$, direct multiplication gives

```math
U|1,r_b\rangle
=|1,r_b\rangle+
c\begin{pmatrix}
t(A-I)|r_b\rangle\\
c(A^\dagger-I)|r_b\rangle
\end{pmatrix}.
```

Thus this support is exactly the span used in Section 8, after redundant
vectors are removed. The two frames can coincide or intersect; the word
requires no lower bound on their angle and no division by a small
singular value. If A or B is identity, the repair word is identity.
Tensoring every operation with dirty identity preserves this argument
for arbitrary borrowed inputs and reference correlations.

### Exact cancellation also applies to actual approximate words

Suppose $`\widehat V`$ is an actual child word and
$`\widehat J`$ is an actual connector word, on the same physical
registers. The child need not have an exactly canonical block matrix.
Let $`\widehat R,\widehat S`$ be actual local wrapper words, and
**construct** the local completion as

```math
\widehat D_B=\widehat R\widehat J\widehat S,
\qquad
\widehat U=\widehat V\widehat J^\dagger.
```

Using literal inverses in
$`\widehat K=\widehat R\widehat U\widehat R^\dagger\widehat U^\dagger`$
and $`\widehat F=\widehat U\widehat D_B`$ gives the exact identity

```math
\widehat K\widehat F
=\widehat R\widehat V\widehat S.
```

It holds on the full physical space, including arbitrary dirty inputs
and any work leakage, without assuming that a helper has been freshly
returned between the displayed factors. The identity does not apply to
an independently compiled local completion merely because it has the
same accepted block. Its complete word must have the matched construction
above. Approximation and return errors of the surviving three factors
still require their usual complete-input estimates.

Emitting the repair and anchored word independently would use the
child word three times: twice forward and once inverse. Expanding and
cancelling the actual words leaves exactly one child call, with no
child inverse. It also removes all occurrences of the connector at that
parent. Applying the same cancellation at every node leaves one copy of
each child program. The identical connector at all leaves is a direct
sum on one shared signal, hence one global signal rotation. This is a
valid call-count reduction; it is not a reduction in the precision cost
of the surviving local wrappers.

### Repair sensitivity does not relax the parent precision

The commutator provides quantitative bounds without choosing a basis for
its moving subspace. For exact unitaries A and B,

```math
\|\mathcal K(A,B)-I\|
\le\|A-I\|\,\|B-I\|.
```

Here $`\|U-I\|=\|A-I\|/2`$ and
$`\|R_B-I\|=\|B-I\|`$; expand the commutator using
$`[R_B,U]=[R_B-I,U-I]`$. For unitary approximations
$`\widehat A,\widehat B`$, a useful asymmetric stability bound is

```math
\begin{aligned}
\|\mathcal K(A,B)-\mathcal K(\widehat A,\widehat B)\|
\le{}&\|B-I\|\,\|A-\widehat A\|\\
&+\|\widehat A-I\|\,\|B-\widehat B\|.
\end{aligned}
```

For the first change, subtract the two conjugations of $`R_B^\dagger-I`$
by U, costing $`2\|R_B-I\|\|U-\widehat U\|`$.
For the second, subtract the two conjugations of
$`\widehat U-I`$ by the parent wrappers, costing
$`2\|\widehat U-I\|\|R_B-R_{\widehat B}\|`$.
The factor-one-half completion identity gives the displayed constants.

This attenuation concerns the repair alone. The complete repaired word
with the parent replaced consistently by $`\widehat B`$ is
$`\mathcal D(A\widehat B)`$. In particular,

```math
\|\mathcal D(A\widehat B)-\mathcal D(AB)\|
=\frac12\|\widehat B-B\|.
```

There is no small factor $`\|A-I\|`$ in this final error. With
arbitrary whole child-word error delta, the matched-inverse argument
still gives $`\delta+\|\widehat B-B\|/2`$. Independently
approximating unmatched wrapper circuits gives only the separately
charged error estimates. The small norm of the repair therefore does
not justify using coarse precision for the parent in the complete word.

### Native implementation and the remaining ledger

At a tree parent, $`B=(U_v^C)^\dagger U_v^W`$ acts on the two
existing child-root modes. Its positive word applies the requested real
$`U_v^W`$ first and the actual native $`(U_v^C)^\dagger`$ second;
its inverse reverses those actual words. The surviving wrappers select
this pair on $`f=0`$ and its inverse on $`f=1`$. The occupied signal
is an unchanged predicate, not an initialized synthesis helper.

The [native controlled templates](#4-native-synthesis-including-external-controls)
implement the real local rotations using a separate borrowed synthesis
signal and the same dirty pool. Literal native coarse words, their
controls, pair routing, predicates, and actual inverses are all charged.
The coarse controls use the retained
[determinant-one reflection interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter),
not an assumption that arbitrary controlled T gates are native.
At a fixed depth the local pairs are disjoint addressed SU(2) blocks;
the two signal values require a constant number of selected templates.
This uses the existing full-operator inactive-sector contract, rather
than assuming a controlled isometry is available at its uncontrolled
cost. No extra clean qubit is supplied by the two-mode root support.

The retained layerwise allocation consequently implements the reduced
word with

```math
\begin{aligned}
T&=O\!\left(L+\sum_{d=0}^{n-1}(2^d+L+n-d+n^2)+n^4\right)
  =O(N+nL),\\
G&=O(NL),\qquad b\ge L+n+7,
\end{aligned}
```

within the stated two-clean budget. The connector is synthesized once;
the fine-precision local tables still occur at every depth. Their complete
word errors include dirty return and reference correlations, and the
actual inverse of a compiled local word is used wherever required.
At $`L=N`$ this construction gives $`O(N\log N)`$, weaker than
the retained grouped bound. Its improvement over a literal uncancelled
repair recursion only removes redundant calls; it reproduces the direct
wrapper construction's resource bottleneck.

The tested route therefore stops here: the commutator supplies a
degeneracy-safe repair and exact inverse cancellation, but no new way
to share the surviving parent precision. A further endpoint advance
requires a different native word or a separate global precision ledger
for those wrappers. This conclusion concerns the displayed construction,
not all low-rank repairs or all complete-frame compilers.

## 10. Collective correction and its precision interfaces

The following results distinguish a collective correction from independently
returned stages. They hold on the complete logical and dirty input space,
including references. None is an unrestricted compiler lower bound.

### Independent-stage precision budgets

Suppose R stages use certified errors $`a_g2^{-\ell_g}`$, with
$`a_g\gt0`$, and their only composition certificate is

```math
\sum_{g=1}^R a_g2^{-\ell_g}\le2^{-L}.
```

AM–GM gives

```math
\sum_g\ell_g\ge R(L+\log_2R)+\sum_g\log_2a_g.
```

Equality in the continuous relaxation assigns error $`2^{-L}/R`$ to
every stage. More generally, with positive precision weights $`w_g`$ and
$`w=\sum_gw_g`$, the minimizing allocation is
$`\epsilon_g=2^{-L}w_g/w`$, giving

```math
\min\sum_gw_g\ell_g
=wL+\sum_gw_g\log_2(a_gw/w_g).
```

This follows either from a Lagrange multiplier or weighted AM–GM. Integer
precision and workspace constraints can only increase this relaxed minimum.
For prefactors and weights bounded above and below by positive absolute
constants, it retains a declared
$`\Omega(RL)`$ precision charge. Correlated errors and joint circuit
rewrites are outside this certificate.

The [conditional-suffix compiler](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#6-exponentially-growing-groups-and-the-explicit-workspace-ledger)
already makes its total table cost $`O(N)`$; its remaining charge is
one length-L source allowance per group. Reallocating these certified
errors therefore does not remove its growing group multiplicity. Nor
does a small active sector reduce operator error: for every nonzero
projector P,

```math
\|(I-P)+e^{i\alpha}P-I\|=|e^{i\alpha}-1|,
```

independently of its rank. A logical input supported in that sector sees
the full discrepancy.

### Orthogonal scalar branches and the low-rank Hopf escape

**Lemma.** Let $`E:\mathbb C^D\to\mathbb C^{\rho D}`$ be any
isometry and let $`P_1,\ldots,P_K`$ be pairwise orthogonal projectors.
They may have unequal ranks and need not sum to identity. If
$`p_j\gt0`$ and

```math
\|E^\dagger P_jE-p_jI_D\|\lt p_j
\qquad\text{for every }j,
```

then $`K\le\rho`$.

*Proof.* Each compression is positive definite, hence has rank D. Thus
$`\mathrm{rank}(P_j)\ge D`$, and orthogonality gives
$`\rho D\ge\sum_j\mathrm{rank}(P_j)\ge KD`$. ∎

For a complete family $`\sum_jP_j=I`$, this gives an explicit uniform
averaging counterexample. Suppose a preparation, sign selection and actual
unpreparation promise

```math
E^\dagger\left(\sum_jf_jP_j\right)E
\approx\frac1K\left(\sum_jf_j\right)I_D
\qquad\text{for every }f\in\{1,-1\}^K.
```

If $`K\gt\rho`$, some compression $`E^\dagger P_jE`$ has a
kernel vector v. Choose only $`f_j=-1`$. On v the actual compression
is identity, whereas the requested average is $`1-2/K`$. Its worst-case
accepted-block error, and therefore its full-output error, is at least
$`2/K`$.

With a supplied a-qubit clean register, $`\rho=2^a`$: every arbitrary
logical and dirty input is included in D. Additional dirty wires multiply
both dimensions equally. Two clean qubits therefore cannot supply this
uniform scalar averaging interface at error below $`2/K`$ for more than
four orthogonal branches, even with an arbitrary entangling encoding. On
an active r-zero-suffix sector the corresponding cap is $`2^{r+2}`$.
This explains the clean-label
constraint of that interface; it is not a restriction on all coefficient
encodings. Anticommuting source terms are not orthogonal label projectors.

The actual [column-forest atoms](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#2-the-residual-support-is-a-union-of-column-forests)
leave a more substantial escape. For a local group dimension $`M=2^s`$,
the forward atoms have the form $`D_jU_j\Pi_j`$. The marker-depth
projectors, together with the zero-marker projector, partition its input:

```math
\sum_j\mathrm{rank}(D_jU_j\Pi_j)
\le\sum_j\mathrm{rank}(\Pi_j)=M.
```

The reverse atoms have the same rank sum; the unchanged prefix multiplies
both sums by its dimension. This statement concerns one set of forward or
reverse coefficient atoms; splitting coefficients into four phase classes
only multiplies the bound by four. Requiring a separate uniform scalar
label replaces each partial atom by a full-rank branch. The lemma does
**not** imply that a joint encoding of the actual Hopf atoms needs
$`\log s`$ clean bits. Their summed partial ranks are compatible with a
constant initialized-dimension ratio, without proving such an encoding.

Inferring the depth from the input marker does not immediately supply it.
The positive column map $`\sum_jU_j\Pi_j`$ is not an isometry:
nested uniform interval columns overlap, and the zero and root-marker
columns coincide. Balanced Haar signs make these columns orthogonal with
the already charged fixed mixer, but the following coefficient filter
still depends on both the original ancestor marker and the current output.
A native joint filter that preserves and erases this correlation remains
missing. The rank lemma does not exclude it.

### A logical correction cannot erase coarse work leakage

Section 1 uses an actual coarse C whose work returns exactly. Replacing it
by a cheaper approximate-return circuit changes the correction problem.
Let J insert the supplied clean work, put $`P=JJ^\dagger`$, and let
$`\mathcal C`$ denote an actual physical coarse circuit. Write
$`U=W\otimes I_b`$ for the desired logical and dirty action.

If a unitary correction F preserves P, then

```math
\|(I-P)F\mathcal CJ\|=\|(I-P)\mathcal CJ\|,
\qquad
\|F\mathcal CJ-JU\|\ge\|(I-P)\mathcal CJ\|.
```

Indeed F is block diagonal with respect to P, so it preserves the norm of
the rejected component. This includes every logical-only correction. For
example, $`\mathcal C=R_y(\alpha)_{\rm flag}\otimes U`$ with a
zero-initialized flag has leakage $`|\sin\alpha|`$ on every input,
which no logical-only left correction changes.

There is also a robust version. If $`\|(I-P)FP\|\le\epsilon`$,
the equal-rank projection identity gives
$`\|F^\dagger PF-P\|=\|(I-P)FP\|`$. Hence

```math
\begin{aligned}
\|F\mathcal CJ-JU\|
&\ge\|(I-F^\dagger PF)\mathcal CJ\|\\
&\ge\|(I-P)\mathcal CJ\|-\epsilon.
\end{aligned}
```

A fine correction whose own initialized-input leakage is at most epsilon
cannot silently erase a much larger coarse leakage. A correction designed
to act on that occupied rejected space has a different interface.

### Exact polar alignment and the unpriced native step

Such a full-space correction exists algebraically. Define

```math
\begin{aligned}
Q&=\mathcal CP\mathcal C^\dagger,\qquad D=P-Q,\\
Z&=PQ+(I-P)(I-Q)=\frac{I+R_PR_Q}{2},\\
R_P&=2P-I,\qquad R_Q=\mathcal CR_P\mathcal C^\dagger.
\end{aligned}
```

Direct multiplication gives $`Z^\dagger Z=I-D^2`$. If
$`\|D\|\lt1`$, the polar unitary

```math
F_{\rm align}=Z(I-D^2)^{-1/2}
```

maps the range of Q onto the range of P. To see the intertwining, use
$`ZQ=PZ`$ and note that $`D^2`$ commutes with both projectors.
Thus $`F_{\rm align}\mathcal CJ`$ returns initialized work exactly;
its accepted action is the polar part of $`J^\dagger\mathcal CJ`$.
This compression can still act nontrivially on dirty work. Its remaining
logical discrepancy and dirty-return error relative to U have not been
corrected.

If $`\|\mathcal CJ-JU\|\le\delta`$, the projection difference
is at most $`2\delta`$. The binomial coefficients
$`c_j={2j\choose j}/4^j\le1`$ give the polynomial

```math
F_k=Z\sum_{j=0}^kc_jD^{2j},\qquad
\|F_k-F_{\rm align}\|
\le\frac{(2\delta)^{2k+2}}{1-4\delta^2},
\qquad2\delta\lt1.
```

The estimate follows by summing the omitted geometric tail and using
$`\|Z\|\le1`$. It is a whole-space matrix estimate, including every
dirty input. **The polynomial $`F_k`$ is not a supplied unitary circuit.**
A block encoding, normalization, amplification, phase synthesis and work
return all require native constructions and charges. The reflections have
explicit ingredients: $`R_P`$ tests clean work and $`R_Q`$ uses
$`\mathcal C`$ and its actual inverse. A one-clean coarse compiler
leaves the second supplied clean wire available as a signal, but this
observation alone does not implement the polynomial.

The direct repeated-query route retains the precision problem. With R
source-bearing groups, a certified coarse accuracy $`\delta=2^{-L_0}`$
and $`L_0\asymp L/R`$ give the retained coarse upper bound
$`O(N+RL_0)=O(N)`$ at $`L=N`$. Making the displayed tail certificate
at most $`2^{-L}`$ requires $`k=\Omega(L/L_0)=\Omega(R)`$
when $`L_0`$ grows. Pricing a degree-$`O(k)`$ reflection polynomial
by $`O(k)`$ complete coarse calls gives the formal query ledger
$`O(k[N+RL_0])=O(NR)`$ at $`k=\Theta(R)`$, before its remaining
operations are priced. This is not a native compiler or a query lower
bound. A useful collective alternative must jointly realize the alignment
and the remaining logical and dirty correction without separately paying
for each complete coarse word.
