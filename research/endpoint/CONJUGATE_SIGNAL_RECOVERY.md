# Borrowed conjugate signals and recovery sectors

[One-clean source](../../docs/ONE_CLEAN_COMPILER.md) · [Rank-weighted blocks](RANK_WEIGHTED_BOUNDARY_COMPILER.md) · [Endpoint contract](../../docs/OPEN_PROBLEM.md)

The literal phase gadget can borrow its signal when the complete high-level
operation is real. It approximates conjugate phases on the two signal
sectors, not the same phase on both. This gives zero-clean real reflection
banks and one-clean real block compilers. The remaining mode/recovery flag
has a different role: its initialization cannot be removed from an
accepted-block guarantee alone.

Throughout, $`N=2^n`$, $`n\ge3`$, and $`\eta=2^{-N}`$. Only T-count is
optimized. Arbitrary dirty/reference inputs, literal phases, and leakage
are included in every stated norm. No measurements, resets, supplied
sources, or uncharged coherent operations are allowed. These are primitive
and promised-family results, not a resolution of the unrestricted endpoint.

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
