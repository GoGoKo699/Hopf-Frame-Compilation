# Borrowed conjugate signals and recovery sectors

[One-clean source](../../docs/ONE_CLEAN_COMPILER.md) · [Rank-weighted blocks](RANK_WEIGHTED_BOUNDARY_COMPILER.md) · [Endpoint contract](../../docs/OPEN_PROBLEM.md)

The literal phase gadget can borrow its signal when the complete high-level
operation is real. It approximates conjugate phases on the two signal
sectors, not the same phase on both. This gives zero-clean real reflection
banks and one-clean real block compilers. The remaining mode/recovery flag
has a different role: its initialization cannot be removed from an
accepted-block guarantee alone.

The fixed-precision statements use $`N=2^n`$, $`n\ge3`$, and
$`\eta=2^{-N}`$; the high-precision forms below use
$`\epsilon=2^{-L}`$. Only T-count is optimized here. Arbitrary
dirty/reference inputs, literal phases, and leakage are included in every
stated norm. No measurements, resets, supplied sources, or uncharged
coherent operations are allowed. The source-lifting interface identifies
which signal can be borrowed in the
[growing-source compiler](GROWING_SOURCE_COMPILER.md), while retaining its
initialized anchor marker and original logical recovery space.

## 1. The actual finite phase word

Use the paired-Majorana source, sign programming, and certified digit rule
of [the one-clean construction](../../docs/ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks).
For one unchanged address, its scalar words are

```math
\mathcal S_f=cI+X_rF,\quad \mathcal S_g=sI+X_rG,\quad c=1/2,
\quad F^\dagger=-F,\quad G^\dagger=-G,\quad
FG+GF=2(cs-\tau)I.
```

The literal phase word, using actual source inverses, is

```math
Q=\mathcal S_f^\dagger S_r\mathcal S_gS_r^\dagger\mathcal S_f
 =\begin{pmatrix}zI&B\\-B&\overline z I\end{pmatrix},
\quad z=s+i(\tau-s/2),\quad B=-i(G/4+FGF).
```

Expansion of the two Pauli sectors gives the formula. The source identities
imply $`B^\dagger=B`$ and $`B^2=(1-|z|^2)I`$. Put $`x=|z|^2`$ and

```math
f(x)=5-20x+16x^2,\qquad g(x)=1-12x+16x^2.
```

With $`Z=Z_r`$, the actual five-call word is

```math
V=QZQ^\dagger ZQZQ^\dagger ZQ
 =\begin{pmatrix}f(x)zI&g(x)B\\-g(x)B&f(x)\overline zI\end{pmatrix}.
```

For $`D_\phi=\mathrm{diag}(e^{i\phi},e^{-i\phi})\otimes I`$ and
$`\alpha=f(x)z-e^{i\phi}`$,

```math
(V-D_\phi)^\dagger(V-D_\phi)
 =\bigl(|\alpha|^2+(1-x)g(x)^2\bigr)I.
```

Thus full operator error equals the flag-zero isometry error. The existing
certified programming estimate yields $`\|V-D_\phi\|\lt30\,2^{-q}`$.
Taking the direct sum over addresses gives $`D\oplus\overline D`$ on an
arbitrary signal input. Replacing this target by $`D\otimes I_r`$ is invalid.
The dirty core and both signal sectors have been retained throughout.

## 2. Lift the whole high-level word

Define

```math
\mathcal P_r(U)=|0\rangle\langle0|_r\otimes U
 +|1\rangle\langle1|_r\otimes\overline U.
```

It preserves products and actual adjoints. For a real final target,
$`\mathcal P_r(W)=I_r\otimes W`$. Exact real gates lift unchanged, while

```math
\mathcal P_r(S_x)=S_x\mathrm{CZ}_{r,x},\qquad
\mathcal P_r(T_x)=T_r^\dagger\mathrm{CNOT}_{x\to r}
 T_r\mathrm{CNOT}_{x\to r}.
```

The latter costs two T/T-dagger gates; its basis phase is
$`\exp(i\pi[(r\mathbin\oplus x)-r]/4)`$. Clifford interlayers remain
Clifford. Lift every interlayer consistently, and use Section 1 directly
for variable diagonal slots. Do not lift an already compiled clean-signal
gadget while pretending its original initialization disappeared.
A full-unitary hybrid bounds the total native error by the sum of slot
errors on arbitrary signal, data, dirty, and reference inputs.

### A real initialized target is sufficient

Let J be a real computational-basis initialization isometry, including
identity on every arbitrary dirty input. Let W be a real unitary on its
input space, and let V be any full unitary satisfying
$`\|VJ-JW\|\le\delta`$. Extend $`\mathcal P_r`$ to rectangular
matrices by the same two-sector formula. Then

```math
\begin{aligned}
J_r&=I_r\otimes J,\\
\mathcal P_r(V)J_r-J_r(I_r\otimes W)
&=|0\rangle\langle0|_r\otimes(VJ-JW)
 +|1\rangle\langle1|_r\otimes\overline{VJ-JW},\\
\|\mathcal P_r(V)J_r-J_r(I_r\otimes W)\|
&=\|VJ-JW\|.
\end{aligned}
```

Complex conjugation preserves the singular values. This identity applies
to arbitrary superpositions and entanglement of r, not only to its basis
inputs. The full V, its coarse approximants, and its action on the
complement of J need not be real. Only J and the specified initialized
target W have to be real. Tensoring an external reference preserves the
same operator norm.

Suppose V is an ideal high-level word with literal diagonal slots and
otherwise native Clifford+T interlayers. Remove the clean phase compiler
signal before forming that word: it is not a data wire or a recovery
coordinate of V. Lift every native interlayer, including all exact coarse
words and actual inverses, and replace each lifted diagonal slot directly
by the finite word of Section 1. If slot j has full-operator error
$`\xi_j`$, the emitted word has initialized error at most
$`\delta+\sum_j\xi_j`$. The interlayer lift at most doubles its
T-count, introduces no initialized wire, and needs only the shared r.
Intermediate T-lift gates may change r; the completed lift implements the
stated two-sector operator exactly. The full-operator hybrid therefore
remains valid even when earlier finite slots have mixed its sectors.

For an unchanged real predicate h, use the
[completed conditional phase word](../../docs/ONE_CLEAN_COMPILER.md#5-exact-conditioning-and-the-resource-ledger).
On $`h=0`$ its actual source loaders cancel with their inverses and
its four signal reflections cancel, including on occupied r. On
$`h=1`$ Section 1 gives the conjugate-phase estimate. Its full target
is thus the lifted conditional diagonal, with zero inactive error.
Lifting a completed guarded interlayer likewise preserves its exact
inactive identity. An identity entry in a phase table alone does not
supply this guard.

### Lift an overlap before its recovery

This form applies when the original fine signal was also tested by an
outer reflection. Let $`J_A`$ initialize the genuine logical recovery
coordinates and act as identity on arbitrary dirty inputs. Start with
ideal high-level sources P and Q, with their fine diagonal slots still
uncompiled and no phase-signal wire, such that

```math
\|J_A^\dagger Q^\dagger PJ_A-W/2\|\le\delta_0,
\qquad W=\overline W.
```

Reserve r outside every source core for the entire construction. In
particular lift the exact cheap scalar words and their replay on the
same physical core, as well as the fine, attenuation, coarse, routing,
and tag operations. Identical original replay words remain identical
after lifting; their rejected dirty operators are conjugated on the
same sector, not replaced by another word with the same scalar block.
The exact lifted sources have accepted overlap
$`\mathcal P_r(J_A^\dagger Q^\dagger PJ_A)`$ on
$`J_r=I_r\otimes J_A`$. Its distance from
$`(I_r\otimes W)/2`$ is at most $`\delta_0`$.

Let $`\widehat P_r,\widehat Q_r`$ be the emitted sources. If all
their fine replacements have total full-operator errors
$`\sigma_P,\sigma_Q`$, respectively, the native unitary
$`\widehat U_r=\widehat Q_r^\dagger\widehat P_r`$ obeys

```math
\left\|J_r^\dagger\widehat U_rJ_r
 -\frac{I_r\otimes W}{2}\right\|
\le\delta_0+\sigma_P+\sigma_Q=:\delta.
```

Certified-data errors may be included in $`\delta_0`$ or in these
slot budgets, but are charged only once. Define

```math
R_r=I-2J_rJ_r^\dagger
   =I_r\otimes(I-2J_AJ_A^\dagger).
```

The [normalization-two bound](../../docs/FAULT_TOLERANT_COMPILER.md#4-from-a-half-unitary-block-to-the-complete-output)
then gives full initialized-output error at most $`12\delta`$ for
$`-\widehat U_rR_r\widehat U_r^\dagger R_r\widehat U_r`$,
including arbitrary r and all actual fine leakage. The reflection does
not test r. It still tests every coordinate genuinely initialized by
$`J_A`$; removing those tests is not licensed by the lift.

For the growing source, $`J_A`$ initializes the original logical
anchor suffix in the active sector. The external marker h remains
initialized and records membership in that sector. Its literal
conditional leading sign and recovery use

```math
Z_h\widehat U_{r,h}R_{r,h}\widehat U_{r,h}^{\dagger}
 R_{r,h}\widehat U_{r,h},\qquad
R_{r,h}=I-2\Pi_{h=1,\,x\in\mathcal A}\otimes I_r.
```

Every completed source and reflection is identity on $`h=0`$.
After recovery, reversing the actual anchor predicate returns h on the
ideal output; its leakage is covered by the same full-isometry bound.
The removed c was the fine compiler signal, whereas the original
logical suffix remains the recovery space. Section 3's unconstrained
complementary-sector action therefore creates no obstruction: both
sectors of the borrowed r are already inside the accepted input space.

### Zero-clean real reflection and state-extension banks

Let real unit vectors have mutually disjoint known computational supports
in a dyadic mode space $`M\le4N`$. Their full reflection
$`I-2\sum_a v_av_a^T`$ has a compiler with

```math
\|\widehat R-R\otimes I_b\|\le\eta/c,\qquad
T=O_c(N),\qquad a=0,\quad b\le N+n+7
```

for fixed $`c\gt0`$. Apply the lift to the complete real reflection word
in [the disjoint-source proof](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources),
including its possibly complex coarse word and phase chart. This does not
assume that its intermediate matrices are real.

Dyadically padded support sizes have sum below $`2M`$. Dyadic items of
size at most M and total at most $`2M`$ pack into two M-capacity bins by
largest-first packing. Pack, apply, and unpack each bin. All padding is
arbitrary data: the ideal full reflection fixes it. There are at most eight
fine diagonal appearances, with the same $`O(\mathrm{poly}(n)\sqrt M)=O(N)`$
coarse/routing cost. The dirty reservation is $`q+2`$: one borrowed
conjugate signal and $`q+1`$ source-core wires. Choosing

```math
q=N+\lceil\log_2(480c)\rceil
```

leaves strict margin for coefficient selection; self-borrowing needs
$`q+1\ge2\log_2 M`$. These conditions and $`q+2\le N+n+7`$ hold for
sufficiently large n. The finite remaining dimensions use the existing
direct literal-synthesis fallback. A real state extension is a reflection
and a known center sign, so it has the same resource order and full-input
contract. Dirty return here is approximate in operator norm, not exact.

### One-clean real block corollary

The [rank-weighted compiler](RANK_WEIGHTED_BOUNDARY_COMPILER.md) has a fixed
number $`F_\kappa`$ of fine slots. Before compiling those slots, lift its
entire ideal high-level word. Retain the initialized mode flag and borrow
the phase signal. Since the embedding and final real target are real, the
ideal initialized-input error is the direct sum of an error and its
conjugate, with equal norm. At
$`q=N+\lceil\log_2(120cF_\kappa)\rceil`$, the fine hybrid is at most
$`\eta/(4c)`$. The original strict representation margins handle the
other errors. The reservation is again $`q+2`$ dirty wires. Therefore

```math
U_j\in SO(m_j),\quad \mathrm{rank}(U_j-I)\le r_j,\quad
\sum_jm_jr_j\le\kappa N
\quad\Longrightarrow\quad
T=O_{\kappa,c}(N),\quad a=1,\quad b\le N+n+7.
```

This includes the squared-size real-block theorem and its promised fusion
and matching corollaries. Complex final targets retain different conjugate
sectors. The genuine mode/recovery flag has not been removed.

### Precision-parametrized real blocks and sources

Fix $`\kappa,c\gt0`$. There is a threshold $`n_0(\kappa,c)`$
such that, whenever $`n\gt n_0`$, $`L\ge6`$, and $`L\gt N/n^2`$,
the following full-contract bound holds for known disjoint real blocks:

```math
U_j\in\mathrm{SO}(m_j),\qquad \sum_jm_j^2\le\kappa N
\quad\Longrightarrow\quad
\begin{gathered}
T=O_{\kappa,c}(N+L),\qquad a=1,\qquad b=L+n+7,\\
\|\widehat UJ_1-J_1(U\otimes I_b)\|\le 2^{-L}/c.
\end{gathered}
```

Use the
[precision-parametrized disjoint-block word](DISJOINT_BLOCK_COMPILER.md#precision-parametrized-high-precision-form)
before compiling its phase slots. Its single mode/recovery flag stays
initialized. Its complete target U is real, though its coarse word and
residual need not be. Apply the initialized-target lifting identity
above to this complete high-level word. The original certified-data
and representation errors retain their norms. A fixed number
$`F_\kappa`$ of fine slots now uses the Section 1 full-operator
primitive at

```math
q=L+\left\lceil\log_2\max\{1,120cF_\kappa\}\right\rceil.
```

Their native hybrid is below $`2^{-L}/(4c)`$, the same allocation
as in the disjoint-block high-precision proof. Together with its
certified-data and representation budgets, the total is at most
$`2^{-L}/c`$. All native interlayers, actual inverses, and exact
predicate signs are lifted; there is no assumption that the residual
or an intermediate source is real. The arbitrary signal r is restored
approximately as part of the final full-isometry error.

The number of slots is independent of L. Each costs
$`O(N+q)=O_{\kappa,c}(N+L)`$; lifting at most doubles the existing
$`O_\kappa(N)`$ coarse/routing T-count. The fine reservation is
$`q+2`$ dirty wires: $`q+1`$ for its core and one separate r.
Its address has at most $`n+1`$ bits, so self-borrowing needs
$`q+1\ge2(n+1)`$. All other coarse work is
$`O_\kappa(\sqrt N+\mathrm{poly}(n))`$ and is reused in separate
stages. Since $`L\gt N/n^2`$, a single larger threshold makes that
coarse work, r, and any fixed returned helper fit within $`L+n+7`$,
and ensures the fine self-borrowing inequality. This threshold is
independent of L. No finite-net argument at unbounded precision is
used.

In the same regime, a real disjoint reflection or state-extension bank
from the earlier source interface, on at most $`4N`$ modes, has
$`T=O_c(N+L)`$, zero clean qubits, $`b=L+n+7`$, and full-operator
error at most $`2^{-L}/c`$. Its final full-mode target is real, so
the lift acts as identity on r in that target. The fixed number of
fine slots, packing into two bins, exact treatment of arbitrary
padding, and source charts are unchanged; only
$`q=N+O_c(1)`$ becomes $`q=L+O_c(1)`$. At most $`n+2`$ address
bits and the same $`q+2`$ reservation fit after the same type of
threshold increase. These interfaces assert T-count and literal
workspace; no additional Clifford-count or depth bound is inferred.

### Banked real interfaces

The independently charged
[banked high-precision interfaces](DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces)
also admit this signal replacement. In their high regime and for
$`b\ge2(L+n+7)`$, the real squared-size block family above has

```math
a=1,\qquad
T=O_{\kappa,c}\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),
\qquad G=O_{\kappa,c}(NL),
```

with the same error $`2^{-L}/c`$. Real reflection and state-extension
banks have these counts with zero clean qubits and constants depending
only on c.

For this statement retain r inside the base pool and outside the fine
core, selectors, and additional banks. The banked proof reserves this
extra constant wire explicitly after increasing its fixed threshold.
Each complete banked query is exactly the original XOR or Z mask on
every input, tensored with identity on its banks. It therefore leaves
the Section 1 full-operator phase calculation unchanged. Lift all
remaining interlayers as above. An S lift adds one Clifford gate; a
T lift uses two T gates and two CNOTs. Thus these lifts send
$`T\mapsto O(T)`$ and $`G\mapsto O(G+T)`$, preserving the displayed
Clifford bound as well as the banked T-count. Their error allocations
are unchanged. For the source interface on up to $`4N`$ modes,
there are at most $`n+2`$ address bits and $`4N`$ table rows;
the same base slack and self-borrowing argument apply after increasing
the threshold, and the extra factors are constant. This banked ledger
does not imply any T-depth bound.

For complete real frames, the separate
[high-precision depth proof](../depth/HIGH_PRECISION_DEPTH.md#2-exact-replacement-at-a-completed-mask-boundary)
uses the real full-operator identity to resynthesize completed XOR/Z
masks without touching r. Together with a separately bounded remainder,
this supplies $`D_T=O(L+NL/b^2)`$ when $`L\gt N/n^2`$ and
$`b\ge17(L+n+7)`$, retaining optimal counts and at most two clean flags.

## 3. Recovery must be checked on both sectors

For an actual unitary with equal sector dimensions, write

```math
U=\begin{pmatrix}A&B\\C&D\end{pmatrix},\qquad
Z_f=\mathrm{diag}(I,-I).
```

Full multiplication gives

```math
\mathcal R(U)=-UZ_fU^\dagger Z_fU
=\begin{pmatrix}
3A-4AA^\dagger A&B(I-4D^\dagger D)\\
C(I-4A^\dagger A)&3D-4DD^\dagger D
\end{pmatrix}.
```

If $`A=E/2`$ with E unitary, unitarity forces
$`BB^\dagger=C^\dagger C=3I/4`$ and
$`D^\dagger D=DD^\dagger=I/4`$. Hence
$`\mathcal R(U)=E\oplus F`$, where $`F=2D`$ is an otherwise unconstrained
unitary. Every pair occurs, for example in

```math
U_{E,F}=\begin{pmatrix}E/2&-\sqrt3 I/2\\\sqrt3 FE/2&F/2\end{pmatrix}.
```

If E and F are real, so is this source. Reality alone does not make F=E.
For approximate balanced blocks,
$`\|A-E/2\|\le\epsilon_0`$, $`\|D-E/2\|\le\epsilon_1`$, a product
telescope gives accepted errors at most $`10\epsilon_i`$ and rejected
errors at most $`6\epsilon_i`$. Thus

```math
\|\mathcal R(U)-(E\oplus E)\|
\le\sqrt{136}\sqrt{\epsilon_0^2+\epsilon_1^2}
\le17\max(\epsilon_0,\epsilon_1).
```

This is a sufficient two-sector source test, not a construction or an
additional requirement on the original two-clean endpoint.

## 4. Exact native regression against deleting the recovery flag

On two dirty core bits use real anticommuting involutions
$`\Gamma_0=X\otimes I`$, $`\Gamma_1=Z\otimes X`$,
$`\Gamma_2=Z\otimes Z`$, and

```math
M=\Gamma_0,\qquad N=\Gamma_0/2+\Gamma_1/\sqrt2+\Gamma_2/2.
```

Two native Majorana-plane conjugations construct N from M with four
T/T-dagger appearances; their scalar phases cancel. The actual scalar word

```math
S_f=H_f C_{f=0}(M)N C_{f=1}(M)H_f
=I/2+X_f(MN-NM)/2
```

is real, satisfies $`Z_fS_fZ_f=S_f^\dagger`$ and $`S_f^3=-I`$, and recovers
to identity. Let x be a logical bit and set $`U=\mathrm{CZ}_{f,x}S_f`$.
The controlled-Z fixes flag zero pointwise and commutes with $`Z_f`$, so

```math
J_0^\dagger UJ_0=I/2,\qquad
\mathcal R(U)=\mathrm{CZ}_{f,x}\otimes I_{\rm core},\qquad
\|\mathcal R(U)-I\|=2.
```

Both source and recovery are real. Another conjugate lift cannot repair
this difference. [The generator](../../compiler_robust_hopf/conjugate_recovery.py)
emits the complete chronological H/S/T/CNOT word, expanding X, Z and CZ.
Its source has 51 gates and four T/T-dagger gates; recovery has 169 and
twelve. These unoptimized counts are a diagnostic, not a lower bound.

[Regression tests](../../tests/test_conjugate_recovery.py) multiply this
16-dimensional word exactly using the repository's retained
$`\mathbb Q(i,\sqrt2)`$ arithmetic. They also test the full finite phase
word, both conjugate sectors, noncommuting block recovery, and perturbations
mixing flags and dirty data. Finite numerical tests do not certify an
endpoint-scale arbitrary-input compiler. Sharing physical helper wires
does not make repeated fine-source calls free.
