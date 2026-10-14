# Disjoint unitary blocks and reflection sources

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Borrowed native words](../../docs/BORROWED_WORKSPACE_COMPILER.md) · [One-clean phases](../../docs/ONE_CLEAN_COMPILER.md)

A direct sum of special-unitary blocks with a linear sum of squared dimensions has a linear-T compiler, even when its entries require exponentially fine accuracy. The construction first emits an exactly logical coarse word, then corrects its residual through two reflection sources and literal amplitude amplification. The reflection source also gives a reusable one-clean operation on an arbitrary mode register.

## 1. Disjoint unitary blocks

Fix constants $`\kappa,c\gt0`$, let $`N=2^n`$ with $`n\ge1`$, and set $`\epsilon=2^{-N}/c`$. Let $`V`$ be a computably specified direct sum of $`SU(m_j)`$ matrices on known disjoint subsets of the $`N`$ logical labels, with identity elsewhere and

```math
\sum_j m_j^2\le\kappa N.
```

There is an effective coherent Clifford+T word $`\widehat V`$ satisfying

```math
T=O_{\kappa,c}(N),\qquad
\|\widehat VJ-J(V\otimes I_{\rm dirty})\|\le\epsilon,
```

with two clean qubits and at most $`N+n+7`$ arbitrary dirty qubits. Here $`J`$ initializes only those two clean qubits. The norm includes every logical column, both returned clean qubits, arbitrary dirty inputs, and reference systems. Every inverse below is the actual inverse of its emitted native word. All literal phases are retained.

The constant is independent of the block entries and their conditioning. Classical certified evaluation and finite search are allowed preprocessing; every coherent table lookup is charged. The determinant-one condition concerns the exact coarse logical words: arbitrary block phases must not be discarded or assumed to have exactly block-supported native approximants.

The proof uses the existing [dirty reflection interpreter](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter), [logical-sector echo](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector), and [literal diagonal primitive](../../docs/ONE_CLEAN_COMPILER.md#literal-diagonal-phases-use-the-same-flag). The new routing, source, and composition arguments are supplied below.

### Precision-parametrized high-precision form

For the same fixed $`\kappa,c\gt0`$, there is a threshold
$`n_0(\kappa,c)`$ such that, whenever

```math
n\ge n_0(\kappa,c),\qquad L\in\mathbb Z,\quad L\ge6,
\qquad L\gt N/n,
```

the same target family has a two-clean compiler with

```math
T=O_{\kappa,c}(N+L),\qquad b\le L+n+7,\qquad
\|\widehat VJ-J(V\otimes I_b)\|\le2^{-L}/c.
```

The constant and threshold are independent of L and the target entries.
The full-input, literal-phase, actual-inverse, and certified-input
conventions are unchanged. The [precision replacement](#precision-replacement-in-the-high-regime)
below proves this form using the same routing and coarse words. The
finite-dimension fallback for the original $`L=N`$ theorem is not used
to assert a bound uniform in arbitrary L.

With additional dirty word banks, the
[banked high-precision interface](#banked-high-precision-interfaces)
improves the count to
$`O_{\kappa,c}(\sqrt{NL}+L+NL/b)`$ at
$`b\ge2(L+n+7)`$, with $`O_{\kappa,c}(NL)`$ elementary Clifford
gates on that same circuit. The corresponding reflection and state
interfaces retain their single clean compiler flag.

## 2. Exact routing and coarse words

### Dirty-only routing and Boolean phases

On $`M=2^k`$ basis labels, every permutation has a Beneš network of $`2k-1`$ addressed-X matching stages. For completeness, form the degree-two bipartite graph connecting input pairs to output pairs under the desired permutation. Alternately color each cycle; the colors route the two channels through separate half-size networks. Input/output pair swaps and recursion give the stated depth. Each matching is one Boolean reflection table. With $`\lambda\asymp\sqrt M`$ dirty banks, the interpreter costs $`O(\sqrt M)`$ T gates per stage. Thus

```math
T_{\rm perm}=O(k\sqrt M),\qquad
b_{\rm perm}=O(\sqrt M+k).
```

The permutation is exact and phase free, uses no initialized workspace, and returns every borrowed wire on its entire input space.

Known Boolean center phases also admit exact dirty-only implementation. On two borrowed bits let

```math
X_4|z\rangle=|z+1\pmod4\rangle,\qquad
Z_4|z\rangle=i^z|z\rangle.
```

The gate $`Z_4`$ is S on the low bit and Z on the high bit. A conditional increment $`X_4^{f(x)}`$ consists chronologically of a high-bit toggle controlled by the old low bit and $`f(x)`$, followed by a low-bit toggle controlled by $`f(x)`$. These are two exact Boolean reflection tables, each priced by the same dirty interpreter. The literal commutator

```math
Z_4 X_4^{f(x)}Z_4^\dagger X_4^{-f(x)}=i^{f(x)}I
```

returns both borrowed bits exactly, including on entangled inputs. Squaring the clock gives $`(-1)^{f(x)}`$. Consequently all center S/Z operations below cost $`O(\sqrt M+k^2)`$ T gates, with no additional clean flag or uncharged phase oracle.

### Exactly logical determinant-one approximants

For either $`R_y`$ or $`R_z`$, X conjugation reverses the rotation. If $`Q_x`$ approximates its quarter-angle up to phase, define the actual word

```math
C_x=XQ_x^\dagger XQ_x,\qquad D_x=C_x^2.
```

Then $`\det C_x=1`$ and $`XC_xX=C_x^\dagger`$ exactly. The unwanted scalar of $`Q_x`$ cancels, and $`D_x`$ approximates the desired literal rotation. Three such factors approximate any SU(2) gate. At inverse-polynomial accuracy in $`N`$, polynomial-length words in $`n`$ suffice; no sharper synthesis-length estimate is needed here. The interpreter emits the resulting table as an exact logical direct sum tensored with identity on its borrowed work. Approximation concerns only the logical block entries.

A small addressed table can be restricted to a larger register's fixed prefix without expanding its row count. Let $`h`$ be that prefix predicate, let $`z`$ be one borrowed bit, let $`T_h`$ toggle $`z`$ by $`h`$, and put $`Y_z=C_zX`$. If $`G_z`$ is the controlled-$`z`$ multiplexor for $`D_x`$, then

```math
Y_z(T_hG_zT_h)G_z^\dagger Y_z=D_x^{h}.
```

Indeed the middle exponent on the original sector $`z`$ is $`(z\mathbin\oplus h)-z`$; conjugation by $`X^z`$ removes its sign. Predicate toggles use the existing exact borrowed-MCX construction. All inactive sectors are literally identity, including occupied inputs to every borrowed bit. The table address consists only of the low bits $`x`$.

## 3. Coarse synthesis with a linear block-entry budget

### Constant-size reduction of oversized blocks

Ignore identity $`SU(1)`$ blocks. There are at most $`64\kappa`$ blocks with $`m\gt\sqrt N/8`$. For sufficiently large $`n`$, partition each such block into chunks of size at most $`d=\lfloor\sqrt N/32\rfloor`$. The number of chunks is bounded in terms of $`\kappa`$ alone.

Block QR eliminates a lower tile by an SU operation on its two row chunks. The stacked tile has at most the upper chunk's number of columns, leaving a nonzero output complement in the two-chunk space; a phase adjustment there makes the row operation determinant one without changing the elimination. Earlier zero tiles stay zero. A block-upper-triangular unitary is block diagonal. Its chunk determinants have product one, so paired diagonal SU(2) phases between chunk pivots transfer their phases to a reference pivot. The remaining chunk factors are individually SU. This gives a bounded number, depending only on $`\kappa`$, of factors supported on at most $`\sqrt N/8`$ labels.

This existence argument does not require an algorithm to decide exact QR pivots. Parameterize each factor by a fixed adjacent-rotation mesh and enumerate finite angle grids. Density and a product perturbation bound guarantee a product strictly within any requested tolerance; certified matrix comparison then terminates. Allocate at most $`\epsilon/8`$ to this preliminary approximation. Every grid factor still has exact determinant one.

Pack the original small blocks into a bounded number $`H_\kappa`$ of bins satisfying

```math
\sum_{j\ {
m in\ bin}}(m_j^2+m_j)\le N/2.
```

Each large-block factor may be its own bin. Since each remaining item is at most $`N/4`$ for sufficiently large $`n`$, greedy packing uses only $`O_\kappa(1)`$ bins. All subsequent constants may depend on this fixed bound.

### Compressed adjacent-row meshes

For one bin, group blocks by $`B=2^{\lceil\log_2m\rceil}`$ and pack them into aligned size-$`B`$ regions. Since $`m\ge2`$, the bin budget implies $`\sum m\le N/4`$ and $`\sum B\lt N/2`$. Padding therefore fits among inactive labels. Disable all gates touching the added padding exactly.

An adjacent-row SU(2) QR rotation eliminating row $`i`$ in column $`j`$ is assigned time

```math
t=B-i+2j,\qquad 0\le j\lt m-1,\quad j+1\le i\lt m.
```

Pairs at one time are disjoint and have the same parity. The lower entry in that column and the preceding-column entry on row $`i-1`$ were eliminated earlier. Thus this is a valid schedule of depth at most $`2B-3`$. The remaining determinant-one diagonal is a product of adjacent-pair Rz phases; cumulative diagonal phases give their parameters, and two matchings implement them.

Even pairs use one target bit. Odd pairs use a cyclic increment of the block suffix, the same matching, and the inverse increment; the wraparound pair is disabled. These increments have polynomial-$`n`$ exact borrowed cost. A separate arbitrary permutation is not paid at every mesh step.

Pack the entire size class into one dyadic region of length $`R_B`$ less than twice its padded size. These regions fit because their total length is less than $`N`$. The compressed-prefix predicate above makes the table length $`S_B=R_B/2`$. Since every nontrivial $`m`$ in this class exceeds $`B/2`$,

```math
S_B=O(\kappa N/B).
```

Choose $`\lambda\asymp\sqrt{S_B}`$. Inverse-polynomial accuracy per mesh gate gives polynomial-$`n`$ word length, including all guards. Its $`O(B)`$ stages cost

```math
T_B=O_\kappa\!\left(\mathrm{poly}(n)\sqrt{NB}\right).
```

Here and below “polynomial” denotes a fixed polynomial, independent of target entries. Summing over dyadic $`B=O(\sqrt N)`$ gives $`O_\kappa(\mathrm{poly}(n)N^{3/4})=O_\kappa(N)`$. One-time pack/unpack permutations also cost $`O(N)`$. The emitted coarse word $`V_0`$ has exact prescribed block support, exact identity on all padding, and exact arbitrary-dirty return, with

```math
\|V_0-V\|\le\frac1{4\sqrt{m_{\max}}}.
```

The peak dirty allocation is $`O(\sqrt N+\mathrm{poly}(n))`$, and no initialized helper is used. With $`R=V_0^\dagger V`$, Cauchy–Schwarz now bounds every row and column absolute sum by

```math
\sum_i|R_{ij}|,\ \sum_j|R_{ij}|\le
1+\sqrt{m_{\max}}\|R-I\|\le5/4.
```

## 4. Disjoint reflection and state sources

This section is also an independent reusable interface. Let $`M\le2N`$ be a dyadic mode dimension on at most $`n+1`$ arbitrary data qubits. Let real unit vectors $`v_a`$ have known mutually disjoint supports. Suppose their supports, each rounded to a dyadic size, have total padded size at most $`M`$. Then the complete mode unitary

```math
R_v=I-2\sum_a|v_a\rangle\langle v_a|
```

has a one-clean compiler using four fine diagonal calls, exact charged interlayers, and $`O(N)`$ T gates at accuracy $`2^{-N}`$ times any fixed constant. It acts on **all mode inputs** and ideally fixes every added padding coordinate. The one initialized wire is the diagonal compiler flag; all mode qubits may be occupied or entangled. The source requires no second initialized mode flag.

### A quantitative two-diagonal local state chart

Pack the supports into aligned dyadic blocks, with a designated center $`e`$ in each block. Singleton reflections are just the exact center Z phases already priced above; no SU(1) state preparation is assumed. On each remaining block of dimension $`B\ge2`$, let $`U`$ be its Walsh transform, $`r=U^Te`$, $`s=\min_j|r_j|=B^{-1/2}`$, and

```math
A(\phi)=U\,\mathrm{diag}(e^{i\phi})U^T,\qquad
D=I+(i-1)|e\rangle\langle e|,
\qquad F=A(\phi)D A(\psi)D^\dagger.
```

Write $`L=U\,\mathrm{diag}(r)`$ and parameterize $`\phi=L^{-1}b`$, $`\psi=L^{-1}c`$, with real $`b,c`$ and $`c_e=0`$. At zero, the derivative of $`Fe`$ is $`ib+c`$. Use the imaginary parts of all amplitudes and the real noncenter amplitudes as local sphere coordinates; the coordinate derivative is the identity.

For parameter vectors $`x=(b,c_{\ne e})`$, the operator norm of the second derivative is at most $`2s^{-2}`$: each diagonal derivative is bounded by $`s^{-1}`$, and the two phase directions contribute a factor at most $`\sqrt2`$ each. On $`\|x\|\le\rho=s^2/8`$ the Jacobian is within $`1/4`$ of identity. For target coordinates $`\|y\|\le\rho/2`$, the map

```math
x\longmapsto x+y-f(x)
```

is a contraction of that ball. Certified iteration computes its solution. The center's real part remains positive there, so matching the coordinates gives literal equality of the normalized states. In particular $`\|v'-e\|\le1/(16B)`$ suffices. No state phase is discarded.

### Exactly logical coarse preparation and full reflection

Prepare each real $`v_a`$ coarsely by a block Hopf word, with crossing-block edges disabled. The dirty interpreter emits an exactly logical native word $`V_0`$, with exact work return and

```math
\|V_0^\dagger v_a-e_a\|\le1/(32B_{\max}).
```

A finite angle-grid search supplies the preparation even at zero subtree masses. Its depth is at most $`\log_2M`$ and its required word precision is polynomial in $`n`$; the total coarse T cost is $`O(\mathrm{poly}(n)\sqrt M)=O(N)`$. Exact coarse words may be complex. Apply the local chart to the actual complex vector $`V_0^\dagger v_a`$, using its exact computable entries. At the ideal-diagonal level,

```math
R_v=V_0F R_{\rm centers}F^\dagger V_0^\dagger.
```

This is equality on the **whole mode space**, not only on the center columns. It therefore fixes the padded coordinates because the target vectors vanish there. There are four variable diagonal slots. Variable-size block Walsh transforms are exact addressed-H reflection tables, one per bit; their total cost is $`O(\log M\sqrt M)`$. Center S/Z phases use the dirty clock echo in Section 2. Packing permutations and every actual inverse are included.

The [one-clean literal diagonal](../../docs/ONE_CLEAN_COMPILER.md#literal-diagonal-phases-use-the-same-flag) has full initialized-isometry error below $`30\,2^{-q}`$. Its [self-borrowed whole-word queries](../../docs/OPERATOR_SOURCE_COMPILER.md#self-borrowed-whole-word-queries) eliminate external selectors when $`q+1\ge2\log_2M`$. Sources are unconditional here; any center predicate work occurs in separate exact subroutines. Thus each fine call uses the same $`q+1`$ arbitrary dirty core, and the reflection error is below $`120\,2^{-q}`$. At $`q=N+O(1)`$, the coarse and fine reservations fit $`N+n+7`$ for sufficiently large $`n`$. Finite smaller dimensions use direct synthesis.

For a real unit state $`z_a`$ in each support, a full-space state extension can be chosen with one such reflection and a fixed center sign. When $`1+(z_a)_e`$ is bounded below, use

```math
V_z=R_{(e+z_a)/\|e+z_a\|}\,R_e,
\qquad V_ze=z_a.
```

The alternative $`R_{(e-z_a)/\|e-z_a\|}`$ works when $`1-(z_a)_e`$ is bounded below. Overlapping certified tests choose one of the two charts without deciding equality. Both choices fix padding on their whole mode space, and require only four fine diagonals. A complex reflection follows by conjugating its real magnitude reflection with two literal phase diagonals, for six fine calls. At zero entries, compute rational coefficient approximants and assign phase zero to rational zero; the reflection perturbation is controlled in vector norm, not by continuity of an argument function.

**Precision-parametrized source interface.** For every fixed $`c\gt0`$,
the same real reflection and real state-extension constructions have
error at most $`2^{-L}/c`$, T-count $`O_c(N+L)`$, and at most
$`L+n+7`$ dirty qubits when $`L\ge6`$, $`L\gt N/n`$, and n is
above a fixed threshold depending only on c. They still use one clean
compiler flag and four fine diagonal calls; the complex reflection uses
six. This is a complete initialized-isometry guarantee on every mode
input, with all clean leakage included. Choose each fine precision
$`q=L+O_c(1)`$, reserve a fixed share of the error for certified chart
and phase data, and use the [shared-pool estimate](#precision-replacement-in-the-high-regime)
below. The coarse word, packing, and interlayers are unchanged. No
state-only contract or new initialized mode register is introduced.

## 5. Two sources for the coarse residual

For each bin compute a rational complex matrix $`\widehat R`$ on every active block with

```math
\|\widehat R-R\|\le\eta,\qquad
\eta=\epsilon/(48H_\kappa),
```

and certify that all row and column absolute sums remain below $`3/2`$. The strict $`5/4`$ margin makes this terminating. Norm accuracy can be obtained from entry accuracy $`\eta/m`$, using $`N+O(n+\log H_\kappa)`$ classical bits. Exact zero tests concern only rational entries of $`\widehat R`$.

Use the first clean flag as a mode bit, giving $`2N`$ modes. Original active inputs are column anchors $`c_j`$. Inside flag one allocate $`m^2`$ pair modes $`(i,j)`$ and $`m`$ row anchors $`r_i`$ per block. The bin budget leaves sufficient space. Define real magnitude states

```math
\begin{aligned}
p_j&=\sqrt{1-\tfrac12\sum_i|\widehat R_{ij}|}\,c_j
       +\sum_i\sqrt{|\widehat R_{ij}|/2}\,|i,j\rangle,\\
q_i&=\sqrt{1-\tfrac12\sum_j|\widehat R_{ij}|}\,r_i
       +\sum_j\sqrt{|\widehat R_{ij}|/2}\,|i,j\rangle.
\end{aligned}
```

Column-source supports are mutually disjoint; row-source supports are separately disjoint. Each family's total unpadded size is $`\sum(m^2+m)\le N/2`$, so dyadic padding is below $`N`$ and fits in the $`2N`$ modes. Householder normals $`(c_j-p_j)/\|c_j-p_j\|`$ and $`(r_i-q_i)/\|r_i-q_i\|`$ are uniformly well conditioned: the residual is near unitary, so its row/column absolute sums are bounded away from zero as well as below two. The reflection interface gives full unitary sources $`P,Q`$ mapping the anchors to these states.

Let $`D_{\rm pair}`$ have phase $`e^{-i\arg\widehat R_{ij}}`$ on a nonzero pair entry and phase one at a zero entry or any other mode. Let $`B`$ exchange column anchors with the corresponding row anchors. Then

```math
U=B^\dagger Q^\dagger D_{\rm pair}^\dagger P
```

has the initialized mode-zero block

```math
A=\bigoplus_j(\widehat R_j/2)\ \oplus\ I_{\rm inactive}.
```

The pair phase has the indicated conjugation so that its contribution is $`\widehat R_{ij}/2`$, including complex signs. All inactive original labels are fixed by the ideal word. Its two reflection families use eight fine diagonals, and its pair phase uses one.

## 6. Complete recovery and workspace ledger

Let $`Z_f`$ be the mode-flag sign. The literal word

```math
\mathcal A=-UZ_fU^\dagger Z_fU
```

has accepted block $`3A-4AA^\dagger A`$. If $`UJ=[A;B_f]`$, its rejected block is exactly

```math
B_f(I-4A^\dagger A).
```

On an active block these are

```math
\frac32\widehat R-\frac12\widehat R\widehat R^\dagger\widehat R,
\qquad B_f(I-\widehat R^\dagger\widehat R).
```

For $`\eta\le1`$, a cubic telescope bounds the accepted error by $`5\eta`$ and the rejected norm by $`3\eta`$. The complete initialized-isometry error is therefore below $`6\eta`$. There is no square-root leakage estimate or postselection. On inactive inputs the polynomial is $`-I`$ with zero rejection; one known diagonal sign corrects it to $`I`$. The leading minus is also literal, for example $`XZXZ`$ on any qubit. Reattach the exact coarse $`V_0`$ on the left.

Each bin has at most $`28`$ fine diagonal calls: nine in each of three forward/inverse overlap calls, plus the inactive sign. The latter could instead use the exact Boolean phase lemma. Put $`C_\kappa=28H_\kappa`$ and choose

```math
q=N+\left\lceil\log_2(120cC_\kappa)\right\rceil.
```

Then the full fine-word hybrid contributes less than $`\epsilon/4`$. Rational residual recovery contributes less than $`\epsilon/8`$ in total, and oversized-block preprocessing at most $`\epsilon/8`$. The resulting error is below $`\epsilon`$ with room for certified source-table rounding.

The only initialized wires are the mode flag and the one-clean diagonal compiler flag. The mode flag is ordinary data during fine calls; their largest address has $`n+1`$ bits. Their shared reservation is $`q+1`$ dirty qubits, with self-borrowing valid once $`q+1\ge2(n+1)`$. Exact routing, coarse words, and center phases reuse this pool and return it exactly. For sufficiently large $`n`$,

```math
q+1\le N+n+7.
```

The complete source reflections are controlled by their full-space unitary contract. No fine word is assumed ideal on an occupied compiler flag, and no flag is reset between operations. A unitary hybrid compares each actual circuit with its ideal action only on the initialized subspace; earlier leakage is propagated by the later actual unitaries. Actual inverse words satisfy the corresponding inverse-target error bound.

All finitely many smaller dimensions have a uniform finite fallback. Direct two-level SU(2) decomposition reduces the target to Ry/Rz factors; determinant phases of intermediate chunks are transferred between pivots as above. Each controlled two-level factor uses the literal determinant-one echo and borrowed logical-sector construction from Section 2. Certified finite native nets at the fixed required tolerance cover the compact target class, without optimizing away a global phase or initializing another register. The available dirty pool suffices for direct sector-by-sector synthesis in each such fixed dimension. Taking the maximum T-to-$`N`$ ratio over these finitely many dimensions enlarges only the fixed constant. The trivial $`n=0`$ special-unitary target is identity. This completes the stated full-contract theorem.

### Precision replacement in the high regime

For the precision-parametrized forms, set $`\epsilon_L=2^{-L}/c`$.
All exact routing estimates, bin capacities, coarse precision
requirements, and local chart radii depend on N and the fixed
block-entry budget, not L. Replace
only the requested preliminary and residual accuracies by their
$`\epsilon_L`$ versions, including
$`\eta_L=\epsilon_L/(48H_\kappa)`$. Certified rational residual
entries then need $`L+O_{\kappa,c}(n)`$ classical bits. These are
classical computations, with no new coherent table oracle.

The same at most $`C_\kappa=28H_\kappa`$ fine slots use

```math
q=L+\left\lceil\log_2(120cC_\kappa)\right\rceil.
```

Their total native initialized-isometry error is below
$`\epsilon_L/4`$. Request the certified diagonal data so that each
ideal diagonal differs from its computed unitary diagonal by at most
$`\epsilon_L/(4C_\kappa)`$ in operator norm, allocating at most
$`\epsilon_L/4`$ to all such data errors. Certified evaluation of the
local charts and phases terminates at this tolerance; their required
classical accuracy does not change the native fine precision q. The
residual recovery and oversized-block preliminary approximation each
contribute at most $`\epsilon_L/8`$ as before. These allocations sum
to less than $`\epsilon_L`$. For a standalone reflection or state
extension use its four (or six) fine slots and the same argument, with
fixed fractions for input/chart evaluation and native realization.

Every fine call costs $`O(N+q)=O_{\kappa,c}(N+L)`$. There are only
constantly many calls. Coarse synthesis and routing still cost
$`O_\kappa(\mathrm{poly}(n)N^{3/4})=O_\kappa(N)`$ for blocks,
and $`O(\mathrm{poly}(n)\sqrt N)=O(N)`$ for the source interface.
Their peak dirty requirement remains bounded by
$`K_\kappa(\sqrt N+p(n))`$ for a fixed polynomial p and constant
$`K_\kappa`$. Consequently, uniformly in
$`L\gt N/n`$, increasing the fixed threshold on n ensures

```math
b_{\rm coarse}\le K_\kappa(\sqrt N+p(n))\le N/n\lt L,
\qquad 2(n+1)\le q+1\le L+n+7.
```

The first inequality fits every coarse reservation into the dirty pool;
the second permits self-borrowing for every fine address and fits its
entire core into that same pool. The reservations are reused between
completed subroutines, not added. The one-clean source statement has
the same estimate, with its constants depending only on c.

All error comparisons remain full-isometry comparisons. Later actual
unitaries propagate earlier leakage without amplification, and actual
inverse words inherit the inverse-target bound. In particular no step
assumes that the shared clean flag is exactly zero after a previous
approximate call. This proves the high-regime forms for every L above
the stated threshold, without a finite-net argument at variable L.

### Banked high-precision interfaces

Fix $`\kappa,c\gt0`$. After increasing the threshold on n by an
amount depending only on these constants, the high-precision block
theorem above strengthens to

```math
\begin{gathered}
L>N/n,\qquad L\ge6,\qquad B_0=L+n+7,\qquad b\ge2B_0,\\
a=2,\qquad
T=O_{\kappa,c}\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O_{\kappa,c}(NL).
\end{gathered}
```

The error is at most $`2^{-L}/c`$ in the same complete
initialized-isometry norm. The real reflection and state-extension
interfaces of Section 4 have these counts with $`a=1`$ and constants
depending only on c; the complex reflection also retains one clean
flag. The ambient mode dimension is still at most $`2N`$, with at
most $`n+1`$ occupied address qubits. These are separate resource
statements for the same respective target families, not a claim that
the block compiler itself needs only one clean qubit.

**Banked fine calls without an extra selector reservation.** Retain the
same constant fine precision $`q=L+c_0`$, with $`c_0`$ fixed by
the number of fine calls and their error allocation. Put $`m=q+1`$
and reserve the same $`B_0`$-wire base pool. By increasing the fixed
threshold on n, this pool fits the complete m-bit dirty core and the
fixed number of separate predicate helpers. The slack
$`B_0-m=n+6-c_0`$ can also reserve one further arbitrary signal
wire outside every core, selector set, and bank, if needed by a
subsequent real-source construction. Also
$`m\ge2(n+1)`$ since $`L>N/n`$. For a mask on a k-bit address,
$`k\le n+1`$, split the core into two balanced parts A and B,
each of size at least k. Implement the query into A using k bits of B
as selectors, then implement the query into B using k bits of A.
Each completed query returns those selectors exactly on every input;
the second therefore accepts the A output of the first. This is the
[self-borrowed query](../../docs/OPERATOR_SOURCE_COMPILER.md#self-borrowed-whole-word-queries)
with the
[exact dirty word-bank query](../../docs/OPERATOR_SOURCE_COMPILER.md#7-trading-additional-dirty-banks-for-lookup-cost)
used for each part.

The additional bank pool has size $`K=b-B_0\ge b/2\ge B_0\ge m`$.
For a table with $`S=2^k`$ rows choose a power of two $`\mu`$
within a factor two below

```math
\max\!\left\{1,\min\!\left(S,\sqrt{S/m},K/m\right)\right\}.
```

Only $`\mu m`$ bank wires are reserved. The A and B queries use
their corresponding portions of those banks sequentially; their
selectors belong to the other core part, never to a bank being loaded
or routed. Their combined complete operation is the original m-bit
XOR query tensored with identity on the bank pool. Its counts are

```math
T_{\rm mask}=O\!\left(\sqrt{Sm}+m+\frac{Sm}{b}\right),
\qquad G_{\rm mask}=O(Sm).
```

For $`k=0`$ the mask is a fixed Pauli word and needs no bank.
Conjugation by core Hadamards gives the same Z mask with the same
orders of cost. Every mask is unchanged as a full unitary, not merely
on the initialized flag columns. The original source words, literal
amplification signs, actual inverses, and complete inactive-sector
guards therefore remain valid. Long source loaders remain
unconditional; only their central gates acquire predicate controls.
All bank and selector work returns before the next source gate.

Since $`S\le2N`$ and $`m=L+O_{\kappa,c}(1)`$, a fine call,
including its $`O(m)`$ source gates and polynomial-in-n predicates,
costs

```math
T_{\rm fine}=O_{\kappa,c}\!\left(\sqrt{NL}+L+\frac{NL}{b}
                                      +\mathrm{poly}(n)\right),
\qquad
G_{\rm fine}=O_{\kappa,c}(NL+\mathrm{poly}(n)).
```

There are only constantly many fine calls. Their precision and errors
are unchanged, so their full-isometry hybrid is exactly the one above.
The extra banks are returned spectators in each complete query, and
are valid even when correlated with the core, mode flags, and an
external reference.

**Coarse Clifford and T counts.** The coarse word, routing, Walsh
interlayers and Boolean center phases do not use endpoint-precision
tables. They retain the dirty interpreter and the polynomial-in-n
word lengths of Sections 2--4. An S-row native-symbol table has
$`O(S)`$ elementary Clifford cost per symbol, including the dirty
traversal and its actual reverse; its one-bit bank routes have no
larger Clifford cost. For a block-size class B, there are $`O(B)`$
mesh stages and $`S_B=O_\kappa(N/B)`$ rows per stage. Thus its
Clifford count is $`O_\kappa(N\mathrm{poly}(n))`$. There
are at most n such size classes and a fixed number of bins, which
only enlarges the fixed polynomial. The one-time permutations have
$`O(n)`$ matching tables each of size $`O(N)`$. The source-only
coarse Hopf word and its variable-size Walsh transforms have at most
$`O(n)`$ layers of at most $`2N`$ rows. Exact arithmetic, center
phases, short native gates and predicate echoes contribute only the
same bound. Consequently all coarse/interlayer Clifford work is
$`O_{\kappa,c}(N\mathrm{poly}(n))`$.

Their already proved total T count is
$`O_{\kappa,c}(N^{3/4}\mathrm{poly}(n))`$; source-only
interfaces have the stronger square-root bound. In the present regime,

```math
\sqrt{NL}>\frac{N}{\sqrt n},\qquad
N^{3/4}\mathrm{poly}(n)=O(\sqrt{NL}),\qquad
N\mathrm{poly}(n)=O(NL).
```

Each comparison is uniform over every $`L>N/n`$ after increasing one
fixed n threshold. The coarse dirty reservation remains
$`O_{\kappa,c}(\sqrt N+\mathrm{poly}(n))\lt L`$, as proved
above, and reuses the base pool between completed fine calls. No
bank or mode coordinate is simultaneously counted as this work.
Adding the constant fine calls proves both displayed gate bounds at
the literal bank threshold. All approximation and actual-inverse
guarantees follow unchanged from the high-precision proof. The same
Clifford ledger applies to the original high-precision interfaces at
$`b\ge B_0`$: their unbanked fine masks also have $`O(Sm)`$
Clifford cost, while all coarse words are unchanged. Thus those
interfaces have $`T=O_{\kappa,c}(N+L)`$ and
$`G=O_{\kappa,c}(NL)`$ on the same circuits and at their original
clean allocations. ∎

## 7. Verification

[Focused tests](../../tests/test_disjoint_block_compiler.py) check the literal dirty clock echo on every input column, the adjacent-row mesh schedule, the local chart's literal complex-state phase, complex pair-source orientation, exact inactive labels, and the complete accepted/rejected recovery identities. These finite matrix checks do not establish native T counts or asymptotic source coverage; those follow from the charged constructions above.
