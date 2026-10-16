# General frame depth by product-coin dilation

[Resource frontier](../../docs/OPEN_PROBLEM.md#large-width-t-depth) ·
[Geometric recovery](GEOMETRIC_BAND_ORIENTATION_DEPTH.md) ·
[Native family sources](COHERENT_DESCENDANT_LOADER.md) ·
[Selected routing](SELECTED_DIRTY_ROUTING.md) ·
[Recursive terminal](RECURSIVE_TERMINAL_DEPTH.md) ·
[Session outcome](INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)

**Local result R87.** At fixed accuracy, arbitrary node-dependent
magnitudes admit the complete-frame depth
$`O_L(N/b^2+\ell\log_2(\ell+2))`$, where $`N=2^n`$ and
$`\ell=\log_2(n+2)`$, with one external clean qubit, the retained
T-count and Clifford-count budgets, and the original literal widths.
The new local decoder represents each branch decision by its own
uniform coin word. After a tensor Hadamard, the last nonzero coin
block gives the marker depth exactly. A paid source and its actual
inverse turn this identity into fine-tag erasure.

This is the unique proof home for the product-coin mechanism.
The [interval-dilation theorem R86](INTERVAL_DILATION_DEPTH.md)
supplies the finite global rounding and the outer geometric allocation.
Its interval-overlap amplification is replaced here. No efficient
classical preprocessing, scalable native emitter, precision-uniform
strengthening, or optimal-depth assertion is made.

## 1. Theorem and inherited interface

Fix integers $`n\ge1,L\ge6`$, and put $`\eta=2^{-L}`$. For every
complete prescribed real frame W and integer dirty width
$`b\ge17(L+n+7)`$, there is a Clifford+T circuit V with one
external initialized qubit and at most b arbitrary dirty qubits such that

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\!\left(N/b^2+\ell\log_2(\ell+2)\right).
\qquad(1)
```

The same statement holds for the prescribed complex magnitude
completion, with supplied leaf phases, at $`b\ge17(L+n+8)`$.
The bounds hold on the same circuit. The norm includes every logical
column, literal phase, leakage and work return, arbitrary dirty input,
and external references. There is no measurement, reset, supplied
precision state, QRAM, or free coherent evaluator.

R86's [finite rounding](INTERVAL_DILATION_DEPTH.md#2-finite-rounding-and-the-atom-embedding)
selects one globally consistent positive magnitude tree with

```math
\tau=\left\lceil2\log_2(n/\eta)\right\rceil+12,\qquad K=2^\tau,
\qquad p_v=k_v/K,\quad 1\le k_v\le K-1.
\qquad(2)
```

Its certified rational approximation, nearest-integer tie rule, clipping,
and finite sign choice require no exact-zero or exact-sign test on the
original coefficients. Its full operator error is at most
$`2n\,2^{-\tau/2}\le\eta/32`$. Every prescribed node is rounded,
including nodes below an originally zero incoming state amplitude.
The rounded signs have exact input and output diagonal gauges.
Every source, projector, local tree, and terminal factor below uses this
same rounded tuple. Section 7 restores the gauges and supplied phases.

The inherited early circuit preserves an outer prefix u of length a,
a coarse bin R with endpoints $`(a,d]`$, and a saved activity bit H.
Its selected descendant inverse exposes a complete q-bit local frame,
$`q=d-a`$, and a zero logical suffix. Write its scaling projectors
as $`A_j^{u,q}`$ and its depth bands as

```math
\Pi_0^{u,q}=A_0^{u,q},\qquad
\Pi_j^{u,q}=A_j^{u,q}-A_{j-1}^{u,q}\quad(1\le j\le q).
\qquad(3)
```

The old fine tag is $`t=a+j`$. Only the first bin includes the
local root $`j=0`$, and then $`a=0`$. It suffices to erase this
tag uniformly on the whole orthogonal sum of bands and preserved
prefixes, while restoring the data and exposed work. A separate
compiler for the entire local frame is not required.

## 2. A product of coin words

Fix a positive rounded q-level local tree. Use q coin blocks
$`C_1,\ldots,C_q`$, each of tau bits. Put $`r=q\tau`$ and
identify a block with an integer in $`\{0,\ldots,K-1\}`$.
For a leaf $`x=x_1\cdots x_q`$, define

```math
B_j(x)=
\begin{cases}
\{0,\ldots,k_{x_{\lt j}}-1\},&x_j=0,\\
\{k_{x_{\lt j}},\ldots,K-1\},&x_j=1,
\end{cases}
\qquad b_j(x)=|B_j(x)|.
```

The preserved outer prefix u is implicit in the node index.
Let

```math
\mathcal R_x=\prod_{j=1}^q B_j(x),\qquad
m_x=\prod_{j=1}^q b_j(x),\qquad
E|x\rangle=m_x^{-1/2}\sum_{z\in\mathcal R_x}|z\rangle.
\qquad(4)
```

These rectangles partition all $`K^q`$ coin strings. Indeed, a
coin string z determines a unique leaf $`\chi_u(z)`$ recursively:
at step j, choose zero exactly when
$`z_j\lt k_{u,x_{\lt j}}`$. Conversely, those comparisons are precisely
the conditions defining $`\mathcal R_x`$. Thus the rectangles are
disjoint, nonempty, and exhaustive; $`E^\dagger E=I`$.
The original leaf probability is $`m_x/K^q`$.

For a prefix v of length e, write $`\mathcal R_v`$ for the product
of its first e constrained coin sets, and $`|\mathcal R_v\rangle`$
for its normalized indicator. Let $`\rho_v`$ be the normalized
positive logical scaling vector below v. Then

```math
E\rho_v=|\mathcal R_v\rangle\otimes
|+\rangle^{\otimes\tau(q-e)}.
\qquad(5)
```

To check every amplitude, put $`m_v=\prod_{j\le e}b_j(v)`$.
The logical coefficient of a descendant leaf x in $`\rho_v`$ is
$`\sqrt{m_x/(m_vK^{q-e})}`$. Multiplication by the coefficient
$`m_x^{-1/2}`$ in (4) leaves
$`(m_vK^{q-e})^{-1/2}`$ on every coin string with decoded prefix v.
Those descendant rectangles partition the unrestricted later coins,
proving (5) without a tensor-product promise on the original probabilities.

At an internal prefix v of length $`j-1`$, the prescribed positive
marker has literal signs

```math
w_v=-\sqrt{1-p_v}\,\rho_{v0}+\sqrt{p_v}\,\rho_{v1}.
```

Consequently

```math
Ew_v=|\mathcal R_v\rangle\otimes\omega_v
       \otimes|+\rangle^{\otimes\tau(q-j)},
\qquad(6)
```

where

```math
\omega_v=
-\sqrt{\frac{1-p_v}{k_v}}\sum_{z\lt k_v}|z\rangle
+\sqrt{\frac{p_v}{K-k_v}}\sum_{z\ge k_v}|z\rangle.
```

The disjoint supports give $`\|\omega_v\|=1`$, and

```math
\langle+|^{\otimes\tau}\omega_v
=-\sqrt{p_v(1-p_v)}+\sqrt{p_v(1-p_v)}=0.
\qquad(7)
```

The root maps to $`|+\rangle^{\otimes r}`$. Equations (5)–(7)
retain the full frame, including every marker and its phase.

## 3. Exact spectral depth from the last nonzero block

Let $`\mathcal H=H_2^{\otimes r}`$, where $`H_2`$ is the
ordinary one-qubit Hadamard, and define the classical function

```math
\lambda(z)=\max\bigl(\{j:z_j\ne0\}\cup\{0\}\bigr).
```

Let $`Q_j`$ project onto coin strings with $`\lambda(z)=j`$.
The later uniform coin blocks in (6) become zero. Equation (7)
means that the distinguished j-th block has no zero component after
Hadamard. Earlier blocks may be in arbitrary superposition. Therefore

```math
\mathcal H E\rho_\varnothing=|0^r\rangle,\qquad
\mathcal H E\,\mathrm{Ran}\Pi_j^{u,q}
\subseteq\mathrm{Ran}Q_j,\quad 0\le j\le q.
\qquad(8)
```

This is an inclusion, not an assertion that the local image fills the
much larger coin sector. Linearity proves it for all superpositions in
each band, and orthogonal sectors handle superpositions across bands.
The whole tau-bit block is tested for nonzero; the last nonzero
individual bit is not itself the logical depth.

Let $`X_c`$ denote bitwise XOR of the binary word c into the fine
tag, and set

```math
\mathcal L_a=\sum_z|z\rangle\langle z|\otimes X_{a+\lambda(z)},
\qquad
\mathcal F_a=\sum_{j=0}^q\Pi_j^{u,q}\otimes X_{a+j}.
```

Equation (8) gives the exact intertwining identity

```math
\mathcal L_a(\mathcal H E\otimes I_T)
=(\mathcal H E\otimes I_T)\mathcal F_a.
\qquad(9)
```

Thus the chronological word

```math
E,\quad \mathcal H,\quad
T\mathrel{\oplus}=a+\lambda(C),\quad
\mathcal H,\quad E^\dagger
\qquad(10)
```

implements the full local spectral-label XOR on its initialized input.
On the consistent old tag $`a+j`$ it erases T and restores the
local vector with phase +1. Its inverse embedding is applied to the
same image after the tag operation; no action of $`E^\dagger`$
outside that image is assumed. The first bin's root has a zero tag.
On a later bin the reached band has $`j\ge1`$; the freely defined
local root case $`a+0`$ is not substituted for a global root.

## 4. Paid physical embedding and actual inverse

Use a q-bit local leaf bank X, an r-bit coin bank C, and a separate
short source bank A. A real family source, preserving u, R, H, and X,
prepares $`|\mathcal R_x\rangle`$ into initialized C and returns A
within its source-isometry tolerance. Follow it with the exact query

```math
X\mathrel{\oplus}=H\chi_u(C).
\qquad(11)
```

On the ideal source image, (11) clears X because
$`\chi_u(z)=x`$ throughout $`\mathcal R_x`$. The physical
embedding is therefore

```math
\overline E|\psi\rangle
=|0^q\rangle_X(E|\psi\rangle)_C|0\rangle_A,
```

with all preserved controls and arbitrary dirty systems unchanged.
The finite rectangle amplitudes are algebraic. The total membership
function is obtained by terminating classical preprocessing on the
finite coin domain. Its evaluation on the quantum register is a charged
truth-table query, not a sequence of free adaptive comparisons.

For the next two estimates restrict to the active H-one sector.
Let J initialize C and A, and let U be this actual source/query word.
The [family-source theorem](COHERENT_DESCENDANT_LOADER.md), with
returned dirty queries, gives

```math
\|UJ-\overline E\|\le\delta_E.
\qquad(12)
```

This estimate is uniform over the other preserved labels and arbitrary
dirty/reference inputs. The physical inverse is exactly the reverse
word: membership query followed by the actual source inverse.

Let $`\mathcal M=\mathcal H\mathcal L_a\mathcal H`$ on the
physical coin and tag banks, identity elsewhere. From (9),
$`\mathcal M\overline E=\overline E\mathcal F_a`$. Hence

```math
\begin{aligned}
\|U^\dagger\mathcal M UJ-J\mathcal F_a\|
&\le\|U^\dagger\mathcal M(UJ-\overline E)\|\\
&\quad+\|U^\dagger\overline E\mathcal F_a-J\mathcal F_a\|\\
&\le2\delta_E.
\end{aligned}
\qquad(13)
```

This full initialized-isometry bound includes all source and coin
leakage and final work return. It holds for arbitrary tag input as
well as the intended consistent tag. It uses no full-space agreement
between U and a guessed source completion. Approximate output is
never reset or treated as an exactly initialized input.

Guard the middle XOR by H. Each completed family source and membership
query is H-zero identity on arbitrary occupied work. The two
Hadamard layers may be unguarded: on H zero their middle word is
identity, so they cancel exactly. For variable r they may act on the
entire fixed physical coin bank. The selected query and source ignore
inactive padding; the Hadamards on that padding commute through and
cancel. Parked dirty values outside the active banks are untouched.
This proves literal inactive identity without testing arbitrary
padding for zero.

## 5. Selected registers, queries, and local cost

Retain the [R86 thin bins](INTERVAL_DILATION_DEPTH.md#thin-bins-and-actual-donor-space).
First reserve disjoint short metadata and
$`m_0=\lceil c_L\ell\rceil\ge4(\tau+2)`$ free logical zeros.
The metadata contains the old fine tag, coarse label, search and
threshold records, source work, and saved-predicate interfaces; it is
not counted as free donor space. Only afterwards choose the independent
protected height $`M=\lceil A_L\ell\rceil`$ and put $`h=n-M`$.

Starting backwards at d=h, set

```math
v=h-d+m_0,\qquad
q=\min\left\{d,\left\lfloor\frac{v}{2(\tau+2)}\right\rfloor\right\},
\qquad a=d-q,\qquad r=q\tau.
\qquad(14)
```

The selected descendant unloader exposes h-d zero early bits.
Together with the untouched reserve these give v donor bits, with

```math
v\ge2(r+2q),\qquad a+q+r=d+r\le h+m_0.
\qquad(15)
```

One r-bit coin bank suffices for (10). The retained extra donor
margin safely accommodates the selected layout. No second atom bank,
square-block output, damping flag, or repeated overlap-amplification
anchor is needed. Constant-factor thinner bins from the design would
also work; (14) is the already proved allocation used here.

The [selected-routing lemma](SELECTED_DIRTY_ROUTING.md) moves local
data into a fixed external leaf bank and donor zeros into the active
coin prefix. Former arbitrary bank values are parked in the original
local/donor slots. One completed layout surrounds the entire word
$`U^\dagger\mathcal M U`$; its actual inverse restores the parked
dirty systems, including correlations with references. Every source
anchor tests only the selected initialized fields.

The relevant large query addresses are:

| Operation | Relevant data bits | Target or action |
|---|---:|---|
| Rectangle family source | $`a+q+r`$ | Common source grammar on active C |
| Coin-to-leaf membership | $`a+r`$ | q-bit XOR into X |
| Fine-tag erasure | r | XOR $`a+\lambda(C)`$ into T |

Short bin, enable, expansion, and calibration fields contribute fixed
polynomial factors. The last query also reads the bin, which determines
a, q, and r. Its target T is disjoint from every relevant address.
The membership target is the fixed external X bank; all preserved
prefix addresses remain in the original early register. Thus X is
disjoint from the union of every bin's prefix/coin address sets.
T and source-program targets are separately reserved and are also
disjoint from the full relevant address union for their query positions.
All truth tables, invalid-sector extensions, and masks are finite
preprocessing data.

Coherently variable target sizes use one common family-source grammar,
with selected active prefixes and returned short masks, exactly as in
[R86's coherent-variable-size interface](INTERVAL_DILATION_DEPTH.md#coherent-variable-sizes).
Whole sources sharing quantum targets are not run in parallel.
Instead each Boolean-query position uses the
[ragged selected-query lemma](SELECTED_DIRTY_ROUTING.md#2-parallel-selected-tables-with-different-relevant-addresses):
each bin reads its own actual address, private simultaneous dirty
pools are summed, and shared outputs use the complete
query/CNOT/unquery/CNOT echo. All non-Clifford supports are private.
No table is padded with the union of every physical public bank and
then charged a smaller exponent.

The retained bin geometry gives

```math
p=O_L(\ell^2),\qquad k_{\rm bin}=\lceil\log_2p\rceil
=O_L(\log(\ell+2)),\qquad
\mathcal Q:=\sum_i2^{a_i+q_i+r_i}\le p\,2^{h+m_0}.
\qquad(16)
```

Here $`q_{\max}=O_L(n/\ell)`$ and $`r_{\max}=O(n)`$.
The source in (12) has a fixed-L tolerance, so its expansion and
calibration precision are fixed-L constants. The existing charged
source and query bounds give, for some fixed-degree polynomial,

```math
T_{\rm loc},w_{\rm loc}\le P_L(n)\sqrt{\mathcal Q},\qquad
G_{\rm loc}\le P_L(n)\mathcal Q,\qquad
D_{\rm loc}=O_L(\ell).
\qquad(17)
```

This includes both source appearances, membership, the fine query,
returned predicates, and all actual inverses. The two tensor Hadamard
layers have zero T-depth and O(n) charged Clifford gates. The boundaries
of the episode add the selected-router depth
$`O_L(\ell\log\log(\ell+4))`$. Source precision need not pay for
one approximation per coin or tree level: (12) is one entire family
source estimate.

## 6. Complete early word and simultaneous global budget

Use the same rounded early tagged source as R85/R86. For increasing
endpoints $`0=d_0\lt\cdots\lt d_p=h`$, compute the coarse
label $`R\mathrel{\oplus}=Hg(T)`$ from the old fine tag, where

```math
g(t)=\min\{i:0\le i\lt p,\ t\le d_{i+1}\}.
```

The first bin includes the root. Compute the selected endpoint from R,
apply its actual descendant inverse, establish the selected layout,
run (10) with the physical U and its actual inverse, reverse the layout,
and reload the descendants. Uncompute the endpoint using unchanged R.
All these controls are independent of T except for the middle XOR.
Equations (8)–(13) clear T on the entire ideal source image and
restore the early data. Only the coarse label remains.

Erase R by the retained coherent binary search, XOR, and actual reverse
search from [R85 Section 7](GEOMETRIC_BAND_ORIENTATION_DEPTH.md#7-complete-recovery-word-phases-and-operator-error).
Its comparison is $`g(t)\le i`$ exactly when
$`t\le d_{i+1}`$. Padded thresholds use h and are always true.
The transported reflection includes its initialized source-work
projector. In particular the endpoint-h reflection is
$`I-2(A_h\otimes|0\rangle\langle0|_A)`$, which is minus identity
on initialized A, not an identity reflection. The old fine tag and old
coarse label are not read by the search. Its
$`k_{\rm bin}`$ rounds cost $`O_L(\ell\log(\ell+2))`$ depth.

Choose the completed source tolerances

```math
\delta_0=\delta_U=\delta_E=\eta/128,\qquad
\delta_R=\frac{\eta}{512(k_{\rm bin}+1)}.
```

The same ideal-source-image hybrid now uses (13) for the local word:

```math
\epsilon_{\rm early}
\le\delta_0+2\delta_U+2\delta_E+8k_{\rm bin}\delta_R
\lt\eta/8.
\qquad(18)
```

Orthogonal bin and prefix sectors take a maximum, not a sum over
columns. Ideal reached states supply the exposed zeros; actual earlier
leakage propagates through unitaries. On H zero the local sandwich
is exact identity, the threshold calculations meet their inverses,
and coarse search meets its inverse with its middle XOR disabled.
Save H from the protected zero predicate before using protected work,
then reverse that actual predicate computation last. H returns on the
ideal output and the existing error is preserved on actual leakage.

By (15)–(17), the complete early word has a fixed polynomial degree c
such that

```math
T_E,w_E\le C_L(n+2)^c\left(2^{(h+m_0)/2}+1\right),\qquad
G_E\le C_L(n+2)^c\left(2^{h+m_0}+1\right),
```

```math
D_E=O_L(\ell\log(\ell+2)).
\qquad(19)
```

The outer sources and comparisons have the same envelopes as R86.
Their finest source precision is
$`O_L(\log(k_{\rm bin}+2))`$; all exponential precision costs
remain fixed polynomials in n. The changed local word uses the same
largest address as R86, fewer source appearances, and one extra
finite last-block query whose volume is smaller. Hence it introduces
no new exponential factor. Simultaneous private query pools and
external banks belong to $`w_E`$; serially returned pools are reused.

Fix c, $`c_L`$, and all metadata before choosing
$`A_L>c_L+2c+20`$, increased if needed for that already fixed
metadata. Since $`2^{h+m_0}=N2^{-(M-m_0)}`$,

```math
T_E=O_L(\sqrt N),\qquad G_E=O_L(N),\qquad
w_E=o(\sqrt N/n).
\qquad(20)
```

The independent protected margin suppresses the table volume.
Its wires have not also enlarged the source address.

## 7. Terminal, literal widths, signs, and phases

Use this wide construction beyond a fixed L-dependent dimension cutoff
and when $`b\ge\sqrt N/n`$. Equation (20) fits the entire early
dirty peak into a fixed fraction of b. The same rounded frame has
the retained exact early/terminal factorization, with terminal height
$`M=O_L(\ell)`$. The [R81 terminal theorem](RECURSIVE_TERMINAL_DEPTH.md),
at error $`\eta/4`$, gives

```math
T_T=O_L(\sqrt N+N/b),\qquad G_T=O_L(N),\qquad
D_T^{\rm terminal}
=O_L\!\left(N/b^2+\ell\log(M+2)\right).
\qquad(21)
```

Its reservation fits beyond the chosen cutoff, and it reuses the
returned dirty pool. Initialized-isometry composition with actual
unitaries does not assume perfect actual intermediate work return.

The exact finite sign gauges are retained Boolean phase tables, with
$`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, and
$`D_T=O_L(N/b^2+\ell)`$. The prescribed complex leaf diagonal is
implemented at error $`\eta/4`$ by the retained one-clean diagonal
compiler with the same budgets. Its phases are included at every
leaf, including leaves whose original state-column amplitude was zero.
The total error is at most

```math
\eta/32+\eta/8+\eta/4+\eta/4\lt\eta.
```

Equations (19)–(21) prove (1) in the wide branch. If
$`b\lt\sqrt N/n`$, use the original
[R72 one-clean compiler](ONE_CLEAN_DEPTH.md) at the original total
accuracy eta. Its depth $`O_L(N/b^2+nL)`$ is
$`O_L(N/b^2)`$, since $`N/b^2>n^2`$. It retains the required
counts and literal widths $`17(L+n+7)`$ and $`17(L+n+8)`$.
The same original-accuracy compiler covers every dimension below the
fixed cutoff, absorbed by L-dependent constants. No stricter internal
tolerance is substituted into this fallback reservation.

This completes the general fixed-accuracy theorem. R79 remains the
absolute all-precision result. At square-root dirty width, the current
general upper depth is $`O_L(\ell\log(\ell+2))`$; a
dimension-growing unrestricted lower bound and optimal general depth
remain open.

## 8. Evidence and relation to the factorization design

The [bounded coin tests](../../tests/test_coin_dilation_depth.py)
check the heterogeneous rectangle partition, every scaling vector and
literal marker, the last-nonzero-block support, and the physical
source/membership/Hadamard/XOR/actual-inverse word with leakage and
references. Finite completions test the interfaces; they are not a
native emitter and do not prove the asymptotic bounds.

The [session record](INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
also retains the independent
[weighted-chain factorization](WEIGHTED_CHAIN_FACTORIZATION.md).
That route changes a complete atom basis through a short sequence of
coefficient stages with bounded absolute norm; its separate native
implementation obligations remain explicit.
The present theorem needs only spectral-depth readout on the embedded
logical image. Its product-coin representation supplies that readout
directly, so the interval-overlap normalization and repeated outer
amplification do not enter its circuit.
