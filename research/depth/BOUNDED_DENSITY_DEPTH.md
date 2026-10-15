# Bounded-density compilation against a positive structured reference

[Retained gauge identity](../endpoint/PRECISION_GAUGE_COMPILER.md) ·
[Native sources](COHERENT_DESCENDANT_LOADER.md) ·
[Common-flag core](HOMOGENEOUS_MAGNITUDE_DEPTH.md) ·
[Affine-oriented core](AFFINE_ORIENTATION_DEPTH.md) ·
[Recursive terminal](RECURSIVE_TERMINAL_DEPTH.md)

**Local result R82.** A positive density bounded above and below relative to
an effectively supplied positive product distribution, or a positive
affine-oriented reference from Section 12, admits a
complete-frame compiler of fixed-accuracy T-depth $`O(\log n\,(1+\log^*n))`$
at square-root dirty width. The reference is finite classical input to
preprocessing, not a supplied quantum state. Its depth-dependent branch
probabilities may be arbitrarily biased. The constants depend on the
relative density bounds and accuracy, not on those biases or the supplied
affine masks.

The proof first treats the uniform reference, then derives the
product-reference version and its affine-oriented extension in Section 12.
The latter uses the explicitly priced [R83 core](AFFINE_ORIENTATION_DEPTH.md),
not an arbitrary reference-frame compiler. The main mechanism is a bounded
positive correction to a charged common-flag frame, implemented by sparse
row/column sources and constant-normalization recovery. The generic
complete-frame frontier remains [R79](UNIFORM_SUBTREE_DEPTH.md).

All logarithms below are binary except in elementary trigonometric
expressions.

## 1. Fixed-accuracy theorem and precise scope

Fix an integer $`L\ge6`$ and known constants $`0\lt c\le1\le C`$. Let $`N=2^n`$ and let
$`W`$ be the complete real Hopf frame prescribed by certified real node
pairs $`(\cos\theta_v,\sin\theta_v)`$. Assume the root state $`\psi`$ has
the global density promise

```math
c\le w_x:=N|\psi_x|^2\le C\qquad\text{for every leaf }x.
```

The compiler need not decide the promise. Classical preprocessing
terminates; no efficiency bound is imposed on it. For every integer
$`b\ge17(L+n+7)`$, there is one native Clifford+$`T`$ circuit $`V`$ using
one external initialized qubit and at most $`b`$ arbitrary dirty qubits such
that

```math
\|VJ_1-J_1(W\otimes I_b)\|\le2^{-L},
```

```math
T=O_{L,c,C}(\sqrt N+N/b),\qquad G=O_{L,c,C}(N),
```

```math
D_T=O_{L,c,C}\!\left(N/b^2+\log_2(n+2)[1+\log^*(n+2)]\right).
```

All bounds refer to the same circuit. The norm covers all logical columns,
literal phases, work leakage and return, arbitrary dirty states, and
external references. Supplied unused wires are untouched. No QRAM,
measurement, reset, free coherent evaluator, supplied precision state, or
uncharged complete-frame oracle is used. At $`b=\Theta(\sqrt N)`$, the
additive depth is $`o(\log^2n)`$.

This is a supporting-family theorem, not a generic improvement. It allows
node-dependent magnitudes under a complementary bounded-density promise. The
uniform-reference and level-homogeneous families overlap; neither generally
contains the other. The product-reference corollary below unifies their
strictly positive parts; its affine-oriented extension admits the concrete
reference family of R83 as well. The same displayed resource and norm
bounds apply to these supplied-reference promises. Constants are permitted
to depend on $`L,c,C`$; there is no all-precision claim.

## 2. Signs and uniform conditional density bounds

The
[retained sign gauge](../endpoint/PRECISION_GAUGE_COMPILER.md#1-contract-and-effective-normalization)
gives exactly

```math
W=D_{\mathrm{sign}}W_{\mathrm{positive}}S_{\mathrm{sign}}.
```

At node $`v`$, write
$`e_0=\mathrm{sgn}(\cos\theta_v),\ e_1=\mathrm{sgn}(\sin\theta_v)`$. The
output leaf sign is the product of edge signs along that path. The
state-input sign is one, and the marker-column sign at $`v`$ is the product
of ancestor edge signs times $`e_0(v)e_1(v)`$. This is the sign-only
specialization of the retained gauge, not a new identity.

Here every local sine and cosine is nonzero. Indeed a child-to-parent
probability ratio is at least $`c/(2C)`$, because child and parent masses
are the corresponding original leaf-density sums. Thus each magnitude is at
least $`\sqrt{c/(2C)}`$. Certified evaluation determines all signs in finite
time, without equality tests. The sign tables are exact finite Boolean data.

Each sign diagonal is implemented by an arbitrary-dirty XOR query $`Q_f`$
into one arbitrary dirty bit $`d`$, followed chronologically by

```math
Q_f,\quad Z_d,\quad Q_f^\dagger,\quad Z_d.
```

The exact phase is $`(-1)^f`$ and all dirty work returns. The retained
sharp/routed query splice gives, at every eligible original width,
$`T=O(\sqrt N+N/b)`$, $`G=O(N)`$, $`D=O(N/b^2+\log_2(n+2))`$; it needs no
clean phase bit. The sign words run sequentially with the positive circuit
and reuse its dirty pool.

For the positive compiler, use the fixed common bounds

```math
w_-=c/C,\qquad w_+=C/c.
```

Every conditional or coarsened subtree density is an average of the ORIGINAL
$`w_x`$ over its descendant block, divided by the average of ORIGINAL
$`w_x`$ over its current prefix block. Hence it always lies in
$`[w_-,w_+]`$. These bounds do not square at successive recursion stages.

In particular, if a preserved prefix $`P`$ leaves $`m`$ active bottom bits,
and $`r`$ bottom bits $`Z`$ are protected, the upper $`h=m-r`$ stage has,
for each $`P`$, a normalized $`h`$-bit distribution $`q_P`$ obtained by
averaging over $`Z`$. Its density $`2^hq_P`$ lies in $`[w_-,w_+]`$. The
exact target stage is

```math
\sum_P|P\rangle\langle P|\otimes
\bigl[W_{P,h}\otimes\Pi_Z+I\otimes(I-\Pi_Z)\bigr],
\qquad \Pi_Z=|0^r\rangle\langle0^r|.
```

This is the ordinary retained complete-frame early/tail factorization.

## 3. Retained gauge identity in the balanced reference frame

Fix one $`h`$-bit positive stage, $`Q=2^h`$, with density $`w_x=Qq_x`$ in
$`[w_-,w_+]`$. Let $`H_{\mathrm{bal}}`$ denote the complete balanced
Haar/Hopf frame (not $`H^{\otimes h}`$). For each node $`v`$ let $`I_v`$ be
its leaf interval, $`l_v`$ and $`r_v`$ the means of $`w`$ on its two
children, and set

```math
g_v=\sqrt{(l_v+r_v)/2},\qquad
\delta_v=\frac{r_v-l_v}{r_v+l_v},\qquad
\tau_v=\frac{2\sqrt{l_vr_v}}{r_v+l_v},
```

```math
C_v=\frac1{g_v\tau_v},\qquad A_v=\frac{\delta_v}{g_v\tau_v}.
```

The state-input diagonal entry of $`C`$ is one. Let $`U_{\mathrm{scal}}`$
have the uniform subtree vector $`u_v=\mathbf1_{I_v}/\sqrt{|I_v|}`$ in
marker column $`\lambda(v)`$, and zero state column. Put
$`D=\mathrm{diag}(\sqrt{w_x})`$. Then exactly

```math
W_q=D(H_{\rm bal}C-U_{\rm scal}A).\qquad (1)
```

This is
[the retained gauge identity](../endpoint/PRECISION_GAUGE_COMPILER.md),
specialized to balanced coarse pairs. The useful new argument below is its
dimension-independent positive-correction norm and sparse native overlap,
not identity (1).

For a direct verification, let $`v_v`$ be the normalized positive weighted
scaling vector on $`I_v`$ and $`e_v`$ the desired positive-frame marker
vector. Then

```math
De_v^{\mathrm{bal}}=g_v(\delta_vv_v+\tau_ve_v),\qquad Du_v=g_vv_v.
```

Their displayed linear combination proves every marker column of (1);
$`Du_{\mathrm{root}}`$ is the root state. Also

```math
\|D\|\le d_0:=\sqrt{w_+},\qquad
C_v^2=\frac{l_v+r_v}{2l_vr_v}\le\frac1{w_-},\qquad
\|C\|\le c_0:=\frac1{\sqrt{w_-}}.
```

## 4. A dimension-independent entrywise-positive correction norm

Let $`\Delta_v=(r_v-l_v)/2`$. Direct substitution gives
$`A_v^2\le\Delta_v^2/w_-^3`$. For every subtree $`I_u`$, finite Haar
orthogonality gives the exact variance identity

```math
\sum_{I_v\subseteq I_u}|I_v|\Delta_v^2
=\sum_{x\in I_u}\bigl(w_x-\mathrm{mean}_{I_u}w\bigr)^2
\le(w_+-w_-)^2|I_u|.
```

Thus $`\lambda_v=|A_v|^2|I_v|`$ satisfies

```math
\sum_{I_v\subseteq I_u}\lambda_v\le B|I_u|,\qquad B=\frac{(w_+-w_-)^2}{w_-^3}.
```

For completeness, this finite Carleson estimate implies

```math
\sum_v\lambda_v|\mathrm{mean}_{I_v}f|^2\le B\|M_d|f|\|_2^2,
```

by integrating level sets: the maximal dyadic intervals on which the mean of
$`|f|`$ exceeds $`t`$ are disjoint, and the sum of $`\lambda`$ below them is
at most $`B`$ times their total size. The finite dyadic maximal function
obeys $`\|M_df\|_2\le2\|f\|_2`$. Indeed the maximal disjoint intervals give

```math
t\,|\{M_d|f|>t\}|\le\sum_{x:M_d|f|(x)>t}|f(x)|.
```

Integrating twice in $`t`$ gives
$`\|M_d|f|\|_2^2\le2\langle|f|,M_d|f|\rangle`$, and Cauchy-Schwarz proves
the bound. Therefore

```math
\|U_{\mathrm{scal}}|A|\|\le K_0:=2(w_+-w_-)w_-^{-3/2}.\qquad (2)
```

The entrywise absolute matrix of $`U_{\mathrm{scal}}A`$ is
$`U_{\mathrm{scal}}|A|`$. This is why a fixed normalization is possible
despite the retained generic full-Haar entrywise barrier. Split continuously

```math
A=A_+-A_-,\qquad A_+=\max(A,0),\qquad A_-=\max(-A,0),
```

```math
K_+=U_{\mathrm{scal}}A_+,\qquad K_-=U_{\mathrm{scal}}A_-.
```

No exact-zero or exact-sign decision is needed for this split. $`K_+`$ and
$`K_-`$ are nonnegative matrices of norm at most $`K_0`$, including zero
state columns.

## 5. Certified Schur weights and the sparse source interface

Choose a fixed computable $`\alpha>K_0`$, for example
$`\alpha=\max(1,2K_0)`$, with $`\alpha=1`$ when $`K_0=0`$. A rational strict
upper bound may instead be used. For either $`K=K_+`$ or $`K_-`$, define

```math
B_K=\begin{pmatrix}0&K\\K^{\mathsf T}&0\end{pmatrix},\qquad \binom uv=(I-B_K/\alpha)^{-1}\mathbf1.
```

Since $`\|B_K\|\lt\alpha`$, the positive Neumann series converges. All
$`u_x,v_j\ge1`$ and

```math
Kv=\alpha(u-1),\qquad K^{\mathsf T}u=\alpha(v-1).
```

This is certified computable from the input with a known strict norm gap.
Neither a Perron eigenvector nor an equality decision is required. Large
intermediate classical numbers do not imply large quantum registers: the
circuit prepares normalized vectors to a stated absolute isometry error.

For a supported pair $`x\in I_v`$, $`j=\lambda(v)`$, put

```math
r_{xj}=\sqrt{\frac{K_{xj}u_x}{\alpha v_j}},\qquad l_{xj}=\sqrt{\frac{K_{xj}v_j}{\alpha u_x}}.
```

These same continuous formulas give zero when $`K_{xj}=0`$; no zero test is
needed. Then

```math
\sum_xr_{xj}^2=1-1/v_j,\qquad \sum_jl_{xj}^2=1-1/u_x,\qquad l_{xj}r_{xj}=K_{xj}/\alpha.
```

Use a common overlap space consisting of the logical $`h`$-bit $`X`$, a
depth tag $`T`$ of $`s=\lceil\log_2(h+1)\rceil`$ bits, two failure bits
$`F`$, and $`O(p)`$ fresh source work. No second $`h`$-bit initialized
leaf/index register is copied. Tag $`t=d+1`$ together with the depth-$`d`$
prefix of $`x`$ identifies the unique marker ancestor $`j`$.

The ideal column source $`R`$ has

```math
R|j\rangle=\sum_xr_{xj}|x,t(j),F=00\rangle
+\frac1{\sqrt{v_j}}|\mathrm{anchor}(j),t(j),F=01\rangle.
```

For state input $`j=0`$, the whole column goes to $`X=0,T=0,F=01`$. For a
marker at depth $`d`$, $`\mathrm{anchor}(j)`$ is its $`d`$-bit prefix
followed by zeros. Distinct marker columns have distinct
$`(\text{prefix},\text{tag})`$ sectors; the root failure has tag zero. Thus
$`R`$ is an isometry.

The ideal row source $`L`$ has

```math
L|x\rangle=\sum_{v:x\in I_v}l_{x,\lambda(v)}|x,t(v),F=00\rangle
+\frac1{\sqrt{u_x}}|x,T=0,F=10\rangle.
```

It preserves $`X`$, and hence is an isometry. Row/column failure supports
are orthogonal. Consequently

```math
L^\dagger R=K/\alpha.\qquad (3)
```

The external prefix $`P`$, any branch enable, and all dirty/reference
systems are preserved.

### Native realization, not supplied sources

$`R`$ first uses exact R78 marker packing, conditionally XORing last-one tag
$`t(X)`$ into $`T`$ and then clearing the marker bit selected by $`T`$. On
initialized $`T`$ this leaves the marker prefix, zero suffix, and its depth
tag. The retained generic real-vector descendant loader prepares the full
normalized vector on suffix $`X`$ plus $`F`$. Include the failure coordinate
in that vector; it remains normalized even when the accepted column
vanishes. The prefix and $`T`$ remain unchanged. Padded unused tags are
assigned fixed normalized targets and arbitrary unitary completions; no
active whole-macro identity on unused tags is assumed. Exact inactive
identity comes from the effective enable, not from those padding choices.

$`L`$ preserves $`X`$ and $`P`$ and applies the same generic real-vector
source construction to the short target $`T`$ plus $`F`$, including its
row-failure coordinate. The geometric Boolean-phase expansion in the
[native source proof](COHERENT_DESCENDANT_LOADER.md) applies to arbitrary
normalized real vectors; its contract is not limited to a probability
product with no failure coordinate.

Each source uses a constant number of native Boolean phase/calibration
queries, their actual inverses, small gradient/coefficient preparation
words, and the inner fixed amplification. Its conditionally initialized
source bank is $`O(p)`$,
$`p=1+\lceil\log_2(1/\delta_{\mathrm{src}})\rceil`$. All dirty pools are
arbitrary and returned exactly by completed query words. Every completed
source is literal identity when its effective enable is zero, including on
occupied $`T,F`$ and source work. On the active initialized input it
approximates its ideal isometry to $`\delta_{\mathrm{src}}`$, with all
source scratch included in the norm.

If $`U_R,U_L`$ are those actual source unitaries on the same padded bank and
$`J`$ initializes that bank, then

```math
\|J^\dagger U_L^\dagger U_RJ-K/\alpha\|\le2\delta_{\mathrm{src}}.
```

No full-space action of $`U_R/U_L`$ on occupied source scratch, or exact
return of approximate source work, is being assumed.

## 6. The balanced factor is an explicitly guarded native word

The balanced factor $`H_{\mathrm{bal}}`$ uses one separate clean $`s`$-bit
tag $`T_H`$ inside the already protected $`Z`$. It does not call a recursive
whole-frame oracle and requires no additional external clean qubit.

For a preserved effective enable $`E`$ (outer saved $`H`$ and, when needed,
a constant-size LCU branch predicate), let $`F_E`$ XOR $`E`$ times the
last-one position of $`X`$ into arbitrary $`T_H`$. Its read-only
conjunction/dirty-echo construction has $`\mathrm{poly}(h)`$ cost and dirty
width, $`O(\log_2(h+2))`$ T-depth, and is exactly identity on $`E=0`$ for
every input.

At each bit let $`U=HZ`$, a balanced real rotation. Write $`U=JZ`$ with
$`J=VZV^\dagger=H`$ and $`V=(SH)T(SH)^\dagger`$, an exact native
single-qubit word. Let $`B=U^{\otimes h}`$. Let $`S_E(T_H)`$ apply $`U`$ on
precisely positions $`j\ge t`$ when $`E=1`$. It is implemented by a selected
$`Z`$ echo followed by a selected $`J`$ echo, each with a read-only
short-tag mask query into arbitrary dirty bits. Controlled $`J`$ is
$`V^\dagger,\mathrm{CZ},V`$ chronologically. These are full-input involution
echoes, not an invalid echo of noninvolutory $`U`$. Effective $`E`$ enters
every query predicate. Their costs are $`\mathrm{poly}(h)`$, depth
$`O(\log_2(h+2))`$; no large Clifford network receives generic controls.

The exact chronological five-stage word is

```math
F_E,\quad S_E(T_H),\quad B^\dagger,\quad F_E,\quad B.\qquad (4)
```

On $`E=1,T_H=0`$, it implements every column of $`H_{\mathrm{bal}}`$ and
returns $`T_H=0`$. The last-one marker remains identifiable after
$`B^\dagger`$, so the second $`F_E`$ erases its tag. On $`E=0`$, both
$`F_E`$ and $`S_E`$ are identity on the full occupied input space, and
$`B^\dagger B`$ cancels. Thus (4) has literal inactive identity even on
arbitrary tag and dirty inputs. Branch controls are read-only inputs of the
predicates; their constant number changes only constants. This is the exact
balanced specialization of the
[common-flag construction](HOMOGENEOUS_MAGNITUDE_DEPTH.md#3-actual-one-qubit-words-and-selected-involution-banks).

## 7. A complete guarded constant-normalization early circuit

Save $`H=[Z=0]`$ in the sole external initialized qubit, using the retained
exact arbitrary-dirty predicate. Do not retest $`Z`$ while it holds
protected work. All following macros preserve $`H`$ and $`P`$. Their
effective enables always include $`H`$, including internal source predicates
and the scalar recovery phase.

Reserve disjoint protected factor banks:

* an attenuation flag and source work for $`C/c_0`$;
* a clean balanced tag $`T_H`$;
* one padded source/tag/failure bank for each $`K_+`$ and $`K_-`$;
* a separate attenuation flag and source work for $`D/d_0`$;
* two LCU branch bits and a fresh damping flag with any preparation work.

The number of these banks is constant, and their total size is
$`O_{c,C}(p+s)`$. Sequential query dirty pools are reused. Distinct
product-factor banks are essential: an occupied rejected $`C`$ or $`K`$ bank
must not be returned by the later $`D`$ factor.

The diagonal attenuators are native real-vector sources on a success flag,
preserving the leaf or marker and $`P`$. Their ideal initialized columns are

```math
\begin{aligned}
|j\rangle|0\rangle&\longmapsto|j\rangle\left[\frac{C_j}{c_0}|0\rangle+\sqrt{1-(C_j/c_0)^2}|1\rangle\right],\\
|x\rangle|0\rangle&\longmapsto|x\rangle\left[\frac{D_x}{d_0}|0\rangle+\sqrt{1-(D_x/d_0)^2}|1\rangle\right].
\end{aligned}
```

Use the generic source contract at tolerance $`\delta_{\mathrm{src}}`$, not
an unpriced table-controlled rotation. Their private source banks are
included in the corresponding factor's acceptance projector.

Short branch enables do not require generic control of a large circuit.
Large Boolean tables simply add the constant-size enable to their read-only
address, a constant row factor. In a coefficient-rotation echo, the product
of selected CZ gates is the phase $`E`$ times the parity of its targets:
compute that parity by Clifford CNOTs, apply one constant-arity
controlled-$`Z`$, and uncompute. Similarly the simultaneous $`X`$ mask used
in the gradient echo is Clifford-conjugate to one constant-arity
controlled-$`X`$ on a pivot. These have constant extra T-depth and
polynomial local Clifford count. All inactive cancellations are identities
for arbitrary occupied work, because the effective $`E`$ includes the outer
$`H`$ in every query and every echoed middle mask.

Let $`\beta_0=c_0+2\alpha`$. Prepare the constant-size branch state

```math
\sqrt{c_0/\beta_0}|0\rangle+\sqrt{\alpha/\beta_0}|1\rangle+\sqrt{\alpha/\beta_0}|2\rangle.
```

Use actual native PREP and its actual inverse. It may be unconditionally
applied: when $`H=0`$ the entire SELECT is identity, so PREP followed by its
actual inverse cancels on arbitrary occupied branch inputs. Approximate PREP
is compared in full operator norm to the desired fixed finite-dimensional
preparation. Branch 3 is identity with zero ideal preparation amplitude.

SELECT has branches

| Branch | Guarded operation |
| --- | --- |
| 0 | The $`C/c_0`$ source, then balanced factor (4) |
| 1 | $`U_{L,+}^\dagger U_{R,+}`$ with conditional scalar minus |
| 2 | $`U_{L,-}^\dagger U_{R,-}`$ |

The scalar minus on branch 1 is an exact controlled phase including $`H`$.
Each branch acts as identity on $`H=0`$ regardless of its occupied private
banks. All branch actual inverses are the reverse native words. Compression
onto all-zero private work after
$`\mathrm{PREP}/\mathrm{SELECT}/\mathrm{PREP}^\dagger`$ is, in the exact
comparison circuit,

```math
\frac{H_{\mathrm{bal}}C-K_++K_-}{\beta_0}.
```

Apply $`D/d_0`$ using its separate bank. Compression onto the union of banks
is then exactly

```math
W_q/\beta,\qquad\beta=d_0\beta_0.\qquad (5)
```

The factorization of compression follows because each new product factor
leaves the earlier factor banks untouched. Arbitrary rejected histories
cannot return to the all-zero product anchor through a later factor.
Constant many source/PREP errors give a dimension-independent block error
$`\epsilon\le C_{c,C}\delta_{\mathrm{src}}`$. This argument is a block
estimate only; full-isometry recovery is supplied next.

Choose a known rational upper bound $`\overline\beta\ge\beta`$ and an
integer $`k>\overline\beta`$. Put

```math
\theta=\frac\pi{4k+2},\qquad a=\sin\theta,\qquad\rho=\beta a\lt1.
```

Add the fresh damping contraction $`\rho`$, including its actual native
error. The complete actual block-encoding word $`U`$ has

```math
\|J^\dagger UJ-aW_q\|\le\epsilon,
```

where $`J`$ initializes ALL its protected work. It is exact identity on
$`H=0`$. The number of source/PREP/SELECT occurrences is constant depending
only on $`c,C`$.

### Full-isometry recovery with a linear error constant

Here is a direct proof that a block estimate suffices, without assuming
approximate work returns exactly. Work on the active $`H=1`$ space and write
$`A=J^\dagger UJ`$, $`C_{\mathrm{rej}}=(I-JJ^\dagger)UJ`$. For $`\epsilon`$
small enough depending only on $`a`$, $`I-A^\dagger A`$ is bounded below by
a positive constant. Thus

```math
C_{\mathrm{rej}}=C_0\sqrt{I-A^\dagger A}
```

with $`C_0`$ an isometry into the rejection space. Define exact ideal
columns

```math
X_*=JaW_q+C_0\sqrt{1-a^2}.
```

They are an isometry. Since $`\|A-aW_q\|\le\epsilon`$,

```math
\|A^\dagger A-a^2I\|\le(2a+\epsilon)\epsilon,
```

and the square-root function has a constant Lipschitz bound on the relevant
positive spectral interval. Consequently $`\|UJ-X_*\|\le C_a\epsilon`$. The
retained
[dimension-free isometry-extension lemma](SUBLINEAR_T_DEPTH.md#4-from-a-source-error-to-an-honest-unitary-comparison)
provides a comparison unitary $`U_*`$ with

```math
U_*J=X_*,\qquad\|U-U_*\|\le6C_a\epsilon.
```

The correction is mathematical only and is never implemented. Take it
sectorwise, preserving $`P,H`$; on $`H=0`$ choose $`U_*=U=I`$. Dirty systems
tensor throughout.

Let $`R_H=I-2HJJ^\dagger`$, i.e. the phase flips precisely $`H=1`$ and ALL
protected-work bits zero. It never tests $`P`$ or the logical $`h`$-bit
$`X`$. Let $`Z_H`$ be the literal phase -1 on $`H=1`$ and +1 on $`H=0`$. The
physical recovery word is

```math
F(U)=[Z_HUR_HU^\dagger R_H]^kU.\qquad (6)
```

On $`H=0`$ it is exactly identity on the complete occupied-work space. On
$`H=1`$ its accepted amplitude is $`\sin((2k+1)\theta)=1`$ for the ideal
$`U_*`$, and its rejected amplitude vanishes. The factor $`Z_H`$ corrects
the otherwise present $`(-1)^k`$ recovery phase. Thus $`F(U_*)J=JW_q`$
literally. The word has exactly $`2k+1`$ forward/inverse $`U`$ occurrences.
Unit norm telescoping gives

```math
\|F(U)J-JW_q\|\le(2k+1)6C_a\epsilon.\qquad (7)
```

Every inverse in the physical circuit is the actual inverse. The all-work
reflection is an exact returned conjunction phase with $`O(\log_2(r+2))`$
depth and $`\mathrm{poly}(r)`$ cost, not a free clean-ancilla test. All
source banks, branch bits, both failures, balanced tag, and damping work
belong to its anchor.

Finally reverse the actual initial saved-predicate word. The ideal stage has
restored $`Z`$ and hence $`H`$; its actual inverse propagates every leakage
component unitarily, including any final $`H`$ error. Set
$`\delta_{\mathrm{src}}=2^{-l}/C_\kappa`$ for a fixed sufficiently large
$`C_\kappa`$ depending only on $`c,C`$. Then the complete stage has full
initialized-isometry error at most $`2^{-l}/16`$. This contract permits
reuse of the same external clean qubit at the next stage, with errors
handled by telescoping and no reset.

## 8. Native local-height/global-prefix ledger

Let $`a_{\mathrm{ext}}=n-m`$, $`h=m-r`$,
$`\overline Q=2^{a_{\mathrm{ext}}+h}=N2^{-r}`$, and $`K=m+l+2`$. Source
precision satisfies $`p=O_{c,C}(l+1)`$. All initialized work is
$`O_{c,C}(l+\log_2(m+2))`$; include the constant number of disjoint factor
banks from section 7.

For $`R`$ table rows, $`q`$ outputs and a disjoint dirty query pool $`B`$,
the exact sharp rectangular query has

```math
T_Q=O(q+\sqrt{Rq}+Rq/B),\qquad G_Q=O(Rq),
```

```math
D_Q=O\!\left(q+Rq/B^2+\log_2(\lceil\log_2R\rceil+2)\right),
```

provided $`B\ge16c_I(\lceil\log_2R\rceil+2)`$. This is the
[sharp rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
In particular it has no multiplicative polynomial in the global address
length.

Every large source query includes the preserved external prefix where
needed. Sparse overlap avoids a second $`h`$-bit index. The following
sufficient table bounds include padding and all actual inverses:

| Completed operation | Sufficient rows | Outputs / multiplicity |
| --- | --- | --- |
| Column $`R`$ phase/calibration | $`C_\kappa(h+1)\overline Qp`$ | at most $`Cp`$; constant many |
| Row $`L`$ phase/calibration | $`C_\kappa(h+1)\overline Qp`$ | at most $`Cp`$; constant many |
| $`D,C`$ attenuation phase/calibration | $`C_\kappa\overline Qp`$ | at most $`Cp`$; constant many |
| Source all-work anchor | $`C_\kappa(h+1)\overline Qp\delta_{\mathrm{src}}^{-2}`$ | one output; constant many |
| Exact marker packing | $`2^h`$ rows and $`O(\log h)`$ outputs | constant many |
| Local suffix/tag masks | $`O(h+1)`$ rows and $`O(h)`$ outputs | parallel-output short query |
| Balanced last-one predicates and selected involutions | no large global table | $`\mathrm{poly}(h)`$ count/width, $`O(\log h)`$ depth |
| Outer all-work phase and saved predicates | no data table | $`\mathrm{poly}(r)`$ count/width, $`O(\log r)`$ depth |
| Coefficient/gradient/PREP/damping | no large data table | $`\mathrm{poly}(l)`$ gates, $`O(l)`$ depth |

For row $`L`$, $`X`$ and $`P`$ are the read-only address and only the short
depth/failure target is prepared. For column $`R`$, tag/prefix identify its
subtree and the descendant $`X`$ supplies the only long target. Therefore
the first two row bounds are $`\overline Q`$ times a polynomial in local
height and precision, not $`\overline Q^2`$. An anchor need not include
$`P`$, but doing so is a valid harmless upper bound. Its factor
$`\delta_{\mathrm{src}}^{-2}`$ comes from
$`2^u\le C\delta_{\mathrm{src}}^{-2}p`$ for one source bank, not from
exponentiating the union of all LCU banks. The outer union-of-banks
reflection is the separately charged conjunction phase.

There are only $`C_{c,C}`$ completed source appearances including outer (6),
whereas R78 has $`O(\log m)`$ sequential comparisons. This constant
multiplicity is the depth improvement. Wide local masks use the retained
parallel-output query, so their $`h`$ outputs do not create a serial
$`h`$-depth term. Large calibration output widths are $`O(l)`$.

To make the polynomial explicit, $`p\le C_\kappa K`$ and
$`\delta_{\mathrm{src}}^{-2}\le C_\kappa4^l`$ give maximum calibration
volume $`C_\kappa K^3\overline Q`$ and maximum individual anchor volume
$`C_\kappa K^2 4^l\overline Q`$. Constant multiplicity leaves their summed
volumes at most $`C_\kappa K^3 4^l\overline Q`$ and summed square-root terms
at most $`C_\kappa K^2 2^l\sqrt{\overline Q}`$. The retained source
preparations, short-mask/predicate helpers and tag packing fit the same
local degree-four envelope used in
[R81 Section 3](RECURSIVE_TERMINAL_DEPTH.md#3-native-costs-depend-polynomially-on-the-local-height);
the new balanced last-one word is independently $`O((h+2)^4)`$, and the
outer anchor/saved predicates are $`O((r+2)^4)`$. Constant branch
parity/fanout guards and constant-dimensional PREP add lower-degree costs.
Since $`r\le m`$ on every compiled stage, these are all bounded by
$`C_\kappa K^4 4^l`$, including the conservative standalone source
allowance. Thus the following degree-ten envelope has ample explicit slack
and covers every listed operation:

```math
\begin{aligned}
T_i&\le C_\kappa\left[K^{10}2^l\sqrt{N2^{-r}}+\frac{K^{10}4^lN2^{-r}}b+K^{10}4^l\right],\\
G_i&\le C_\kappa\left[K^{10}4^lN2^{-r}+K^{10}4^l\right],\\
D_{T,i}&\le C_\kappa\left[\log_2(n+l+2)+l+\log_2(m+2)+\frac{K^{10}4^lN2^{-r}}{b^2}\right].
\end{aligned}\qquad (8)
```

A disjoint $`B=\Theta(b)`$ pool is reused sequentially for large queries;
all nonquery dirty reservations have a fixed polynomial bound in $`m+l`$.
The displayed standalone term deliberately covers small preparations,
guards, last-one circuitry, and predicates. It cannot be discarded. In (8),
the global logarithmic query term occurs only a constant number of times per
stage. There is no multiplicative $`\log m`$ decoder factor.

## 9. Recursive schedule, error and count sums

Choose a fixed sufficiently large $`A\ge40`$, enlarged for the total
conditional-clean footprint and the degree-ten envelope in (8). Let
$`\kappa`$ be a fixed integer large enough for the stage error constants.
Start $`m_0=n`$. At each nonfinal stage choose

```math
l_i=L+2\lceil\log_2(m_i+2)\rceil+\kappa,\qquad K_i=m_i+l_i+2,
```

```math
r_i=\lceil A[l_i+\log_2K_i]\rceil,\qquad m_{i+1}=r_i.\qquad (9)
```

There is a fixed stopping height $`m_*=m_*(L,c,C,A)`$ so that above it
$`r(m)\le m/2`$ and $`r(m)=O_{L,c,C}(\log_2(m+2))`$. Stop at $`m_t\le m_*`$.
Every earlier stage has $`h_i\ge m_i/2`$ and a protected bank large enough
for all fields in section 7. The early/tail factorization leaves precisely
the original $`m_{i+1}`$ bottom layers, selected by the now preserved
external prefix. The density guarantee in section 2 therefore applies afresh
with the SAME $`w_-,w_+`$.

Stage errors are at most

```math
\epsilon_i\le\frac{2^{-l_i}}{16}\le\frac{2^{-L}}{16\cdot2^\kappa(m_i+2)^2}.
```

Since the decreasing heights have $`r_i\le m_i/2`$, their reversed sequence
grows geometrically. Thus $`\sum_i(m_i+2)^{-2}`$ is bounded by an absolute
constant. Choose $`\kappa`$ to make the stage sum below $`2^{-L}/16`$. For
the final at most $`m_*`$ original layers, use the retained R72
smaller-error bottom with $`q_j=L+j+7`$ and $`m_j=L+j+8`$, as in R81. Its
full error is below $`15\cdot2^{-L}/64`$. All exact sign words add no error.
Full-isometry telescoping gives total below $`2^{-L}`$ (with substantial
reserve).

The protected-bank suppression gives, for $`A\ge40`$ and the degree-ten
envelope,

```math
K_i^{10}2^{l_i}2^{-r_i/2}\le K_i^{-10},\qquad K_i^{10}4^{l_i}2^{-r_i}\le K_i^{-30}.
```

Their sums are finite because the heights grow geometrically in reverse.
Hence all large tables in (8) sum to

```math
T=O_{c,C}(\sqrt N+N/b),\qquad G=O_{c,C}(N),\qquad D_{\mathrm{tables}}=O_{c,C}(N/b^2).
```

The standalone terms obey $`4^{l_i}\le C_{L,c,C}(m_i+2)^4`$ and
$`K_i=O_{L,c,C}(m_i+2)`$. Therefore

```math
\sum_iK_i^{10}4^{l_i}\le C_{L,c,C}(n+2)^{14}=O_{L,c,C}(\sqrt N).
```

The corresponding Clifford bound is weaker and fits $`O(N)`$. There is no
unpaid $`\log^*n`$ factor in the count budgets. The finite bottom has

```math
T_{\mathrm{bottom}}=O_{L,c,C}(\sqrt N+N/b+n^2),\qquad G_{\mathrm{bottom}}=O_{L,c,C}(N),
```

```math
D_{\mathrm{bottom}}=O_{L,c,C}(N/b^2+\log_2(n+2)),
```

and $`n^2=O(\sqrt N)`$.

The number $`t`$ of nonfinal stages is $`O_{L,c,C}(1+\log^*(n+2))`$. One
elementary justification is that for sufficiently large $`m`$ the
composition of two updates $`f(f(m))`$, with $`f(m)\le A'(1+\log m)`$, is at
most $`\log m`$. Thus two stages perform at least one ordinary logarithmic
reduction until a fixed threshold, after which only a fixed number remain.
Also

```math
l_i=O_{L,c,C}(1+\log_2(m_i+2))=O_{L,c,C}(\log_2(n+2)).
```

Consequently summing the additive depth in (8) gives

```math
D_{\mathrm{additive}}=O_{L,c,C}\!\left(\log_2(n+2)[1+\log^*(n+2)]\right).
```

The two exact sign phases and constant-height bottom obey the same final
bound.

## 10. Every dirty width at the original threshold

For sufficiently large $`n`$ (cutoff depending only on fixed $`L,c,C`$ and
native constants), split at $`b_c=\sqrt N/n`$.

If $`b\ge b_c`$, reserve $`B=b/4`$ for each sequential large query and at
most $`b/4`$ for all local polynomial-in-$`(n+l_0)`$ dirty helpers.
Exponential dominance gives room for these local pools and ensures

```math
B\ge16c_I(\lceil\log_2R\rceil+2)
```

for every query, because all addresses have at most
$`n+O_{L,c,C}(\log_2(n+2))`$ bits. Output registers, preserved $`P,X`$, and
protected initialized work are not charged as part of $`B`$. The final
bottom fits its retained shallow schedule as well. All sources/predicates
return their query pools before reuse. The protected work consists of
original logical bits $`Z`$, not additional initialized ancillas; only $`H`$
is external.

If $`b\lt b_c`$, use the retained R79 compiler on the ORIGINAL target and
ORIGINAL $`L`$. It already fits $`b\ge17(L+n+7)`$, has the required count
and error, and at fixed $`L`$ has additive depth $`O_L(\log_2^2(n+2))`$.
Since $`N/b^2>n^2`$, that additive term is absorbed by the width term. For
the finitely many $`n`$ below the large-dimension cutoff, use R79 at every
width; its finite additive depth is absorbed in the claimed bound with
constants depending on $`L,c,C`$.

Thus no larger precision reservation or unknown indicator constant is
silently required at the literal threshold $`17(L+n+7)`$. The proof selects
one complete branch circuit, not incompatible count-optimal and
depth-optimal circuits.

## 11. The full initialized-isometry contract

The construction and its actual inverse have the same resource bounds. If
the forward word satisfies $`\|VJ_1-J_1W\|\le\eta`$, then

```math
V^\dagger J_1-J_1W^\dagger
=V^\dagger(J_1W-VJ_1)W^\dagger,
```

so its actual inverse has the same error, including all work and arbitrary
reference correlations. Composition between recursive stages and with the
exact signs always uses this initialized-isometry estimate. No intermediate
return is promoted from approximate to exact.

## 12. Positive product-reference corollary

A stronger version has the same asymptotic resources and proof. Suppose the
classical input also supplies a certified strictly positive product
probability

```math
q_x=\prod_{j=1}^nq_j(x_j),\qquad
q_j(0)=a_j^2>0,\quad q_j(1)=b_j^2>0,\quad a_j^2+b_j^2=1,
```

and promises, for $`\mu_x=|\psi_x|^2`$,

```math
c\le\mu_x/q_x\le C.
```

There is no uniform lower bound on $`q_x`$ or its one-bit factors. The
constants in the theorem depend only on $`L,c,C`$, not on the biases or the
smallest $`q`$ atom. The reference $`q`$ is supplied effectively; searching
for such a reference is not asserted. Uniform $`q`$ is the theorem above.
Taking $`\mu=q`$ includes every strictly positive level-homogeneous
magnitude input. Taking bounded nonproduct perturbations of $`q`$ allows
independent node-dependent magnitudes. Reference zeros are excluded from
this corollary.

Use the complete depth-homogeneous Hopf frame $`W_q`$ as reference. For a
node $`v`$ let $`a,b`$ be its reference child amplitudes, let $`l`$, $`r`$
be the child conditional $`q`$-means of $`w=\mu/q`$, and put

```math
g^2=a^2l+b^2r,\qquad C_v=\frac g{\sqrt{lr}},\qquad A_v=\frac{ab(r-l)}{g\sqrt{lr}}.
```

Let $`U_q`$ have the positive reference scaling vector

```math
(u_v^q)_x=\sqrt{q_x/q(I_v)}\,\mathbf1_{x\in I_v}
```

in marker column $`v`$, and zero state column. Then the retained general
gauge identity is

```math
W_\mu=D(W_qC-U_qA),\qquad D=\mathrm{diag}\sqrt w.
```

For verification, $`De_v^q=g(\alpha_vv_v+\beta_ve_v)`$, with

```math
\alpha_v=\frac{ab(r-l)}{g^2},\qquad\beta_v=\frac{\sqrt{lr}}{g^2}.
```

Hence $`C=1/(g\beta)`$ and $`A=\alpha/(g\beta)`$, giving the stated formula.
The bounds remain

```math
C_v^2=a^2/r+b^2/l\le1/w_-,\qquad A_v^2\le\frac{a^2b^2(r-l)^2}{w_-^3}.
```

Finite $`q`$-weighted variance decomposes into the increments

```math
\sum_{I_v\subseteq I_u}q(I_v)a_v^2b_v^2(r_v-l_v)^2
=\sum_{x\in I_u}q_x\bigl(w_x-\mathbb E_q[w\mid I_u]\bigr)^2.
```

It follows that $`\lambda_v=A_v^2q(I_v)`$ has weighted Carleson constant
$`B=(w_+-w_-)^2/w_-^3`$. For a leaf vector $`f`$ in ordinary Euclidean norm,
set $`F_x=f_x/\sqrt{q_x}`$. Then

```math
\|(U_q|A|)^\dagger f\|^2
=\sum_v\lambda_v|\mathbb E_q[F\mid I_v]|^2
\le4B\sum_xq_x|F_x|^2
=4B\|f\|^2.
```

The finite weighted maximal inequality has the same elementary level-set
proof: replace interval cardinality by $`q`$ mass throughout. This argument
has no dependence on a minimum $`q`$ atom. Thus exactly the same $`\alpha`$
and Schur-resolvent row/column formulas apply to $`K_+=U_qA_+`$ and
$`K_-=U_qA_-`$. Their support is still the same subtree interval, so the
single-leaf/depth overlap, failure padding, table sizes, and native source
multiplicities are unchanged.

For each early reference of height $`h`$, compile the common single-qubit
factors into actual native $`U_j=J_jZ`$, $`J_j=V_jZV_j^\dagger`$, with
per-factor operator error at most $`\delta_{\mathrm{src}}/(Ch)`$. The
common-flag five-stage word (4), with $`B=\bigotimes_jU_j`$ and the selected
$`J_j/Z`$ echoes, implements its perturbed reference frame exactly and
returns the clean tag exactly. The target reference frame differs by at most
$`O(\delta_{\mathrm{src}})`$: each logical depth layer is a direct sum of
the same one-qubit factor, and $`h`$-layer telescoping adds the per-factor
errors. This is charged as another block error in section 7; the exact $`q`$
remains in all classical $`D,C,K`$ tables.

Ancilla-free synthesis of each $`V_j`$ to this tolerance has T-depth
$`O(\log(h/\delta_{\mathrm{src}}))`$ and count
$`O(\log(h/\delta_{\mathrm{src}}))`$; all $`j`$ factors are parallel inside
each bank. Actual adjoints and selected involution echoes preserve literal
inactive identity even though the synthesized $`U_j`$ need not be exactly
real. The reference factor therefore adds $`\mathrm{poly}(h+l)`$ count/width
and $`O(\log_2(h+2)+l)`$ depth, already covered in (8). It is not an
unpriced $`W_q`$ oracle, and no exact tag-return claim is made for an
arbitrary unrelated approximate frame compiler.

Finally, every recursive conditional ratio is

```math
\frac{\mathbb E_q[w\mid\text{current prefix plus active block}]}{\mathbb E_q[w\mid\text{current prefix}]},
```

using ORIGINAL $`w=\mu/q`$ and ORIGINAL product $`q`$. It lies in
$`[c/C,C/c]`$ at every stage, and the conditional reference remains a
product of the remaining depth factors. The same common bounds therefore
apply without repeated squaring. Exact target signs are certified computable
because every target amplitude and local branch coefficient is nonzero;
unlike the uniform-reference special case, no fixed lower magnitude is
needed for finite termination. The sign queries and all-width fallback are
unchanged.

This proves the product-reference version with the same stage, error, count,
and literal-width bounds.

### Affine-oriented reference extension

The same R82 bounds hold for the following larger, concretely supplied
reference family. At global position j, supply certified strictly positive
$`a_j,b_j`$ with $`a_j^2+b_j^2=1`$, and finite affine orientation data

```math
r_j(v)=\rho_j+\sum_{k\lt j}A_{jk}v_k
\quad\text{over }\mathbb F_2.
```

The positive reference tree has node pair $`(a_j,b_j)`$ when
$`r_j(v)=0`$ and $`(b_j,a_j)`$ otherwise. Let q be its root probability
distribution and $`W_q`$ its complete positive frame. Promise, as before,

```math
c\le\mu_x/q_x\le C,\qquad \mu_x=|\psi_x|^2.
```

This reference is supplied by its certified magnitudes and binary affine
masks; no search for a suitable reference is asserted. The constants depend
only on fixed $`L,c,C`$, not on the masks, branch biases, or smallest q atom.
Reference zeros remain excluded. Product references are included by taking
every orientation mask zero. The complete native reference factor is the
[direct shared-tag core of R83](AFFINE_ORIENTATION_DEPTH.md#6-a-complete-early-frame-using-one-shared-tag),
whose exact leaf/marker maps and guarded implementation are proved there.
This extension does not import an arbitrary reference-frame oracle.

**Conditioning and coarsening.** Let P be a preserved prefix of length a,
let m bottom bits remain, protect r of them, and put $`h=m-r`$. Define
the normalized h-bit distributions by descendant masses:

```math
\mu_{P,h}(x)=\frac{\mu(Px)}{\mu(P)},\qquad
q_{P,h}(x)=\frac{q(Px)}{q(P)}.
```

With the original global $`w=\mu/q`$,

```math
\frac{\mu_{P,h}(x)}{q_{P,h}(x)}
=\frac{\mathbb E_q[w\mid Px]}{\mathbb E_q[w\mid P]}
\in[c/C,C/c].
```

These are the same fixed $`w_-,w_+`$ used throughout the proof. The local
reference orientation at position j is

```math
r_{a+j}(Pu)=\rho_{a+j}
+\sum_{k\le a}A_{a+j,k}P_k
+\sum_{k\lt j}A_{a+j,a+k}u_k.
```

It remains affine in the local prefix u, with constants affine in unchanged
P. Summing descendants leaves the earlier normalized branch probabilities
unchanged. Thus every conditional, coarsened reference has exactly the
outside-prefix interface supported by the direct R83 core. Strict
positivity and the lower density-ratio bound also make every target branch
nonzero, so the exact target sign tables are certified computable without
a uniform minimum-amplitude assumption.

**The weighted correction is unchanged.** In the preceding gauge formula,
let $`a_v,b_v`$ be the actual reference child amplitudes at the node,
including its orientation, and let $`l_v,r_v`$ be the two child conditional
means of the local density ratio. The two-child calculation still gives

```math
W_\mu=D[W_qC-U_qA],\qquad
C_v=\frac{g_v}{\sqrt{l_vr_v}},\qquad
A_v=\frac{a_vb_v(r_v-l_v)}{g_v\sqrt{l_vr_v}},
\qquad g_v^2=a_v^2l_v+b_v^2r_v.
```

The bounds on D,C,A and the finite q-weighted variance identity use only
these conditional probabilities, not independence of the reference bits.
Therefore $`\lambda_v=A_v^2q(I_v)`$ has the same Carleson constant
$`(w_+-w_-)^2/w_-^3`$, and
$`\|U_q|A|\|\le2(w_+-w_-)w_-^{-3/2}`$. The positive Schur resolvent,
orthogonal failure flags, and sparse leaf/depth overlap are unchanged.
All q-dependent row and column vectors are still prepared by the charged
generic source. The reference changes their classical table contents,
not their coherent label count or support.

**One dedicated tag, with no inner clean signal.** Replace the reference
factor inside Section 7's SELECT by the direct local R83 core. It uses the
same dedicated short tag inside the outer protected logical bank Z, first
for its signed marker permutation and then for its common-flag word. Both
parts return the tag exactly. It uses no recursive whole-frame call or
additional initialized activity flag. Its effective enable is the saved H
and a constant-size branch condition; the bounded-arity native guards in
R83 include that condition directly and remain identity on occupied work
when disabled. All attenuation, Schur-source, branch, and damping banks
remain disjoint as in Section 7.

For local tolerance $`\delta_{\rm src}`$, use the same actual native
single-qubit factors in the reference source and its common conjugator.
The exact native leaf/marker maps surround that word. It returns its tag
and dirty work exactly and differs from $`W_{q_{P,h}}`$ by
$`O(\delta_{\rm src})`$ in operator norm, uniformly over P. Keep the exact
q in D,C,K and all source tables. Since $`\|C/c_0\|\le1`$, this reference
approximation adds only $`O(\delta_{\rm src})`$ to the accepted-block
error. The same fixed-round recovery and linear full-isometry comparison
then apply. There is no assumption of exact return by the other
approximate source factors.

**The additional standalone cost.** R83's short-tag-selected quadratic
phases may contain individual outside-prefix bits. They have addresses
$`(T,E,z_j,P_i)`$, so their physical T-depth is local even though their
number depends polynomially on the original n. A sufficient per-reference
call allowance is

```math
T_{\rm ref},G_{\rm ref}
=O\!\left((n+2)^4+(h+1)[l+\log_2(h+2)]\right),\qquad
w_{\rm ref}\le C_1(n+2)^4,
```

```math
D_{T,\rm ref}=O(\log_2(h+2)+l).
```

Here $`C_1`$ is an absolute native width constant, and
$`\delta_{\rm src}=2^{-l}/C_\kappa`$ as in Section 7; fixed density
constants are absorbed in the displayed count and depth constants.
Single-qubit synthesis is ancilla-free. The shared tag still requires only
$`O(\log(h+2))`$ conditional clean bits.

Append these standalone terms to (8). They are not purely local-height
costs. The large source/calibration/anchor volumes and their constant
multiplicity are unchanged. Under Section 9's recursion the additional
count sum is a fixed polynomial in n, hence
$`O_{L,c,C}(\sqrt N)`$ T gates and $`O_{L,c,C}(N)`$ Clifford gates.
Its depth is already covered by the additive $`l+\log(m+2)`$ in (8).
At $`b\ge\sqrt N/n`$, the polynomial dirty pool fits beside the large-query
pools above a fixed dimension cutoff. Reuse it only after actual return.
Below that width, and in the finite small-dimension range, Section 10's
original-precision R79 fallback is unchanged. The real literal threshold
therefore remains $`17(L+n+7)`$, with the same one-clean full-norm theorem.

The direct R83 theorem is sharper when the target is itself the supplied
affine-oriented reference: it gives logarithmic depth. This corollary
allows arbitrary node-dependent perturbations within the fixed density
ratios and retains R82's iterated-logarithmic factor. For arbitrary
orientations the leaf map need not be affine and the marker correction
need not be a polynomial collection of short-address quadratic phases;
the priced reference interface above cannot be omitted.


## 13. Prescribed complex phases on the same one clean qubit

For the prescribed phase-dressed completion $`D_\phi W_{\mathbb R}`$, assume
any of the real reference promises above. Run the real theorem at precision
$`L+1`$ and the retained
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at that same precision. Both fit

```math
b\ge17(L+n+8).
```

At fixed $`L`$ the diagonal has

```math
T_D=O_L(\sqrt N+N/b),\qquad G_D=O_L(N),\qquad
D_{T,D}=O_L(N/b^2+\log_2(n+2)).
```

Reuse the sole external signal and the same dirty pool. With
$`R=W_{\mathbb R}\otimes I_b`$ and $`D=D_\phi\otimes I_b`$, the actual
composition obeys

```math
V_DV_RJ_1-J_1DR
=V_D(V_RJ_1-J_1R)+(V_DJ_1-J_1D)R.
```

Its norm is at most $`2\cdot2^{-(L+1)}=2^{-L}`$. Earlier actual leakage is
propagated by the next actual unitary; only the ideal comparison term uses
initialized work. Counts and depths add, while the width requirement is
their maximum. Thus the same R82 resource bounds hold, with the displayed
complex threshold. This is the prescribed phase-dressed magnitude
completion, not arbitrary complex-unitary synthesis. No supplied quantum
reference or phase state is used.

## 14. Verification boundary and scope

The proof combines retained signed-gauge algebra, the native real-vector
source, the common-flag selected-involution core, sharp dirty queries, the
retained isometry extension, and the retained early/tail factorization. Its
additions are the finite Carleson normalization, sparse Schur row/column
source assembly, constant-normalization complete early LCU, and recursive
relative-density ledger.

The bounded checks in
[`tests/test_weighted_haar_depth.py`](../../tests/test_weighted_haar_depth.py)
cover finite variance/Carleson identities, uniform and vanishing correction
columns, independent weighted-frame reconstruction, padded row/column
overlaps, complete three-branch LCU with separate product flags, literal
amplification, and the actual inverse of a perturbed source. These support
the algebraic and full initialized-space identities. The asymptotic
resources follow from the native table-volume, physical-support, recursion,
and all-width arguments above; finite fixtures do not establish them or
constitute an emitted large native compiler.

The affine-oriented extension additionally uses the complete native
leaf/marker and shared-tag identities proved in
[R83](AFFINE_ORIENTATION_DEPTH.md). Its new outside-prefix standalone
polynomial costs and direct reference-core interface are charged in
Section 12; the large weighted-source table ledger is unchanged.

The density promise excludes arbitrarily concentrated magnitude tables
unless they remain within fixed factors of a supplied positive product or
affine-oriented reference. No claim of an efficient classical test for that
promise or an efficient search for the reference is made. The unrestricted growing-depth
lower-bound and generic fixed-accuracy improvement remain separate open
tasks.
