# Assembling the complete residual with two clean flags

[Research status](OPEN_PROBLEM.md) · [Tree transport](ENDPOINT_TREE_TRANSPORT.md) · [Weighted blocks](WEIGHTED_TRANSPORT_BLOCK.md)

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
retained $`O(N\ell_*(n))`$ complete-frame endpoint bound. The useful
new ingredient is an explicit assembly within two clean flags; a separate
precision charge remains at every tree depth.

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
$`|h_v+d_vz|^2+\sum_b\rho_{2v+b}|u_0[b]+u_1[b]z|^2-s^2|z|^2`$.
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
[robust normalization-two lemma](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
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
