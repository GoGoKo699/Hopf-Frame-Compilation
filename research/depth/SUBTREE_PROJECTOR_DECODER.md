# Polylogarithmic fixed-accuracy depth from subtree projectors

[One-clean all-precision theorem](ONE_CLEAN_DEPTH.md) · [Coherent descendant loader](COHERENT_DESCENDANT_LOADER.md) · [Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md) · [Resource statement](../../docs/OPEN_PROBLEM.md)

**Local result R78.** The complete early frame has nested subtree scaling
spaces. Native reflections about these spaces recover a column's depth
label coherently in logarithmically many queries. Computing that label,
XORing it into the source tag, and reversing the computation erases the
tag without amplitude amplification. The original descendant loader
supplies each reflection; the target frame is never supplied as an oracle.
The resulting one-clean compiler has polylogarithmic T-depth at fixed
accuracy and retains the optimal T-count and linear Clifford budget.

## 1. Fixed-accuracy theorem and scope

Fix $`0\lt\eta\le1/64`$, and set

```math
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad
N=2^n,\qquad \ell_n=\log_2(n+2).
```

For every integer $`n\ge1`$ and the indicated dirty allocation, each
prescribed operator W has one coherent Clifford+T circuit V:

| Prescribed family | Required arbitrary dirty width |
|---|---|
| Complete real Hopf frame | $`b\ge17(L+n+7)`$ |
| Complex magnitude frame $`D_\phi W_{\mathbb R}`$ | $`b\ge17(L+n+8)`$ |

There is only one external initialized qubit. With its initialization
isometry $`J_1`$, the same circuit satisfies

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\mathrm{op}}\le\eta,
```

```math
T=O_\eta\!\left(\sqrt N+\frac Nb\right),\qquad
G=O_\eta(N),\qquad
D_T=O_\eta\!\left(\frac N{b^2}+\ell_n^2\right).
```

Constants may depend on the fixed accuracy, but not on n, b, or the angle
tuple. The bound includes every logical column, literal phase, all work
leakage and return, arbitrary dirty inputs, and external references.
Unused supplied wires remain untouched. T-depth permits arbitrary
Clifford interlayers; G counts their elementary gates.

In particular, at $`b=\lceil C_b\sqrt N\rceil`$ for any fixed positive
$`C_b`$ and sufficiently large n,

```math
T=O_\eta(\sqrt N),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(\log^2(n+2)).
```

Choose a sufficiently large fixed $`C_b`$ if this specialization must
meet the literal width threshold also in every finite dimension.
The inherited worst-case lower bounds match T at every eligible width.
They match depth throughout each nonempty interval

```math
17B_0\le b\le\frac{\sqrt N}{\ell_n}:
\qquad T^\star=\Theta_\eta(N/b),\qquad
D_T^\star=\Theta_\eta(N/b^2),
```

where $`B_0`$ is the family-specific base in the table. Integer-width
nonemptiness means $`17B_0\le\lfloor\sqrt N/\ell_n\rfloor`$.
Optimal depth beyond this interval remains open.

The inputs have the same terminating certified evaluation interface as
R72, including signed, unequal, and singular angles and prescribed leaf
phases at zero amplitudes. Classical preprocessing may be expensive;
no running-time bound is added. Every query, source preparation, actual
inverse, guard, and reflection is charged. There is no QRAM, measurement,
reset, supplied precision state, or free coherent evaluator. This is a
fixed-accuracy refinement, not a precision-uniform replacement for R72.

## 2. The literal nested subtree spaces

For an h-bit frame, let a depth-d node have binary prefix p. Define full
h-bit vectors, including the prefix, recursively:

```math
v_{h,p}=|p\rangle,
```

```math
\begin{aligned}
v_{d,p}&=c_{d,p}v_{d+1,p0}+s_{d,p}v_{d+1,p1},\\
e_{d,p}&=-s_{d,p}v_{d+1,p0}+c_{d,p}v_{d+1,p1},
\end{aligned}\qquad
c_{d,p}=\cos\theta_{d,p},\quad s_{d,p}=\sin\theta_{d,p}.
```

These are the [literal scaling and marker vectors](../../docs/HOPF_INTERFACE.md#marker-columns).
Their definition never divides by an incoming amplitude. Disjoint child
supports make both vectors normalized and orthogonal even at singular
angles. At each node they are an orthonormal change of basis of its two
child scaling vectors. Consequently

```math
A_d=\sum_{|p|=d}|v_{d,p}\rangle\langle v_{d,p}|,\qquad
A_{d+1}=A_d+\sum_{|p|=d}|e_{d,p}\rangle\langle e_{d,p}|.
```

Thus $`A_0\le A_1\le\cdots\le A_h=I`$. The mutually orthogonal
projectors

```math
\Pi_0=A_0,\qquad \Pi_{d+1}=A_{d+1}-A_d
```

sum to I. They project onto the state column and the depth-d marker
columns, respectively. In particular,

```math
A_d=\sum_{t=0}^d\Pi_t,\qquad
(I-2A_d)x_t=(-1)^{[t\le d]}x_t
\quad (x_t\in\mathrm{Ran}\,\Pi_t).
```

The initial source below preserves the literal signs of the individual
marker vectors. The subsequent projector operations act as scalars on
each entire band, so they do not alter those signs or the orientation
inside a band. For padded thresholds $`d\ge h`$, define $`A_d=I`$.

## 3. Native reflections from a source-isometry guarantee

Use the [original descendant loader](COHERENT_DESCENDANT_LOADER.md#1-source-contract-and-disjoint-descendant-sectors),
with local gradient preparation on each call. Let D be a preserved
$`s=\lceil\log_2(h+1)\rceil`$-bit threshold register. Its effective
depth is $`d'=\min(d,h)`$. The target at label (d,p) is the normalized
subtree scaling vector $`v_{d',p}`$. This is a legal effectively
specified real-vector family for the source contract, including the
scalar +1 target when $`d'=h`$.

Let A contain **all** fresh initialized source work: its program,
gradient, expansion index, and attenuation bit. No retained gradient
is used. Define the anchor and desired-output isometries within D=d by

```math
J_d|p\rangle=|p0^{h-d'}\rangle_X|0\rangle_A,\qquad
B_d|p\rangle=v_{d',p}|0\rangle_A.
```

Use the source's existing saved activity flag H. On its active block,
the native loader $`\widehat U`$ preserves D and the prefix exactly,
returns its arbitrary dirty pool exactly, and satisfies

```math
\|\widehat U_dJ_d-B_d\|\le\delta.
```

Different prefixes and selectors occupy orthogonal invariant sectors,
so the coherent error is a maximum, with no number-of-columns factor.
The source tag and search record are spectators for this loader; H is
read but unchanged. The following source and projector formulas refer
to H=1. On H=0 the inherited completed loader is exactly identity on
its full input space. No additional activity flag is supplied.

Put $`P_d=J_dJ_d^\dagger`$ and $`V_d=\widehat U_dJ_d`$. Since
both $`V_d`$ and $`B_d`$ are isometries,

```math
\begin{aligned}
\|V_dV_d^\dagger-B_dB_d^\dagger\|
&\le\|(V_d-B_d)V_d^\dagger\|
    +\|B_d(V_d^\dagger-B_d^\dagger)\|\\
&\le2\delta.
\end{aligned}
```

Therefore the actual conjugated reflection has the **full operator-norm**
bound

```math
\left\|\widehat U_d(I-2P_d)\widehat U_d^\dagger
-\left[I-2\left(A_{d'}\otimes|0\rangle\langle0|_A\right)\right]
\right\|\le4\delta.
\qquad\text{(1)}
```

The ideal comparison is identity off the indicated output subspace; it
is not asserted to be $`(I-2A_{d'})\otimes I_A`$ on occupied A.
Equation (1) includes source leakage and arbitrary unitary completion
on the other source inputs. No polar completion is required. The
direct sum over D preserves the same bound, as does tensoring dirty
systems and references. The dagger is the actual reversed native word.

### Exact controls and the paid central predicate

Let H be the saved activity flag and c the current search control. With
$`P_D=\sum_d|d\rangle\langle d|\otimes P_d`$, use

```math
\widehat C=
\widehat U\left(I-2[H=1,c=1]P_D\right)\widehat U^\dagger.
```

Chronological execution is $`\widehat U^\dagger`$, central phase,
$`\widehat U`$. Neither exterior word changes H or c. When either
control is zero the actual exterior words cancel exactly, on arbitrary
occupied A and other registers, with literal phase +1. Controlled
synthesis of the loader is unnecessary. Replacing the transported
projector by $`A_{d'}\otimes|0\rangle\langle0|_A`$ defines
$`C_*`$ with $`\|\widehat C-C_*\|\le4\delta`$ globally.

The central phase is an exact arbitrary-dirty Boolean phase for

```math
f(D,X,A,H,c)=Hc\,[A=0]\,[X_{d'+1:h}=0].
```

The suffix test is vacuous when $`d'=h`$. The zero test excludes the
prefix, D, the occupied source tag, search record, H, c, and all dirty
helpers. The [dirty phase echo](SHALLOW_DIRTY_QUERIES.md#2-a-boolean-phase-with-arbitrary-dirty-work)
compiles the whole predicate. If A has u bits its table has at most

```math
2^{h+s+u+2}=O((h+1)2^h2^u)
```

rows, one output, and returned arbitrary dirty work. This table, rather
than an unpriced long zero test, is included in Section 6. Padded scalar
loaders need not be exact identities; their source guarantee already
proves (1) for $`d\ge h`$.

## 4. Reversible spectral search and tag erasure

Let Y and D be fresh s-bit registers. For $`j=s-1,\ldots,0`$,
perform the following chronological stage:

1. XOR into D the threshold
   $`d_j(Y)=\sum_{k>j}2^kY_k+(2^j-1)`$.
   Copy the higher bits by CNOT, leave bit j zero, and flip the lower bits.
2. Apply a Hadamard to $`Y_j`$, the controlled reflection
   $`\widehat C`$ with $`c=Y_j`$, then a Hadamard to $`Y_j`$.
3. Apply $`\mathrm{CNOT}(H;Y_j)`$ and reverse the threshold preparation.

Threshold preparation reads only higher Y bits, which this stage does
not modify. It needs no carries, comparator, or stored midpoint history.
The final CNOT corrects a literal sign: on active initialized A the
Hadamard sandwich computes $`[t\le d_j]`$, while the desired bit is
its complement. Equivalently, a controlled active minus sign converts
$`I-2A_d`$ to $`2A_d-I`$ before the sandwich.

On $`x_t\in\mathrm{Ran}\,\Pi_t`$ and the already correct higher
bits of Y, the threshold is the last integer in the lower half of the
current binary interval. Hence $`[t\le d_j]=1-t_j`$. Starting from
Y=D=A=0, the exact comparison search $`Q_*`$ satisfies

```math
Q_*\,x_t|0\rangle_Y|0\rangle_D|0\rangle_A
=x_t|t\rangle_Y|0\rangle_D|0\rangle_A
\quad(H=1).
\qquad\text{(2)}
```

There is no residual scalar phase. When a padded threshold exceeds h,
membership is true and the returned bit is zero, as required. Every
comparison stage preserves the data vector within its band and returns
D and A. At H=0 each actual stage is identity even on occupied Y,D,A:
the reflection cancels, the two Hadamards cancel, the CNOT is disabled,
and threshold preparation meets its actual inverse.

Let T denote the old s-bit source tag, which the search never reads.
Use one exact Toffoli per bit for the activity-guarded XOR

```math
C_H:\quad T\longmapsto T\mathbin\oplus HY,
\qquad \widehat E=\widehat Q^\dagger C_H\widehat Q.
```

On H=1 and initialized Y,D,A the ideal comparison has the spectral-XOR action

```math
E_*=\sum_{t=0}^h\Pi_t\otimes X_t,
\qquad X_t|z\rangle_T=|z\mathbin\oplus t\rangle_T,

\qquad\text{(3)}
```

with Y,D,A returned. This displayed formula describes the initialized
subspace; $`Q_*^\dagger C_HQ_*`$ specifies its full unitary extension.
Uncomputation uses X and Y, never the erased T. The guard on the middle
XOR is essential for inactive occupied tags and search records. At H=0,
the whole actual word is exactly identity also by
$`\widehat Q^\dagger I\widehat Q=I`$.

There are s reflections in the forward search and s in its actual
inverse. Equation (1), exact intervening arithmetic, and unitary
telescoping imply the global bounds

```math
\|\widehat Q-Q_*\|\le4s\delta,\qquad
\|\widehat E-E_*\|\le8s\delta.
\qquad\text{(4)}
```

In (4), $`E_*`$ denotes the full comparison just specified. Full-norm
comparison permits arbitrary earlier leakage into A; no step assumes
that actual work has returned exactly to zero.
The separately approximated native projectors need not be nested or
commute. Only the exact comparison projectors have those properties;
the global error estimate compares the entire actual search to them.

## 5. Complete early frame with one external clean qubit

Choose a protected terminal block of M original logical bits Z, with
$`h=n-M\ge1`$, and let X be the first h bits. The exact product of
the first h addressed layers is

```math
W_E=W_h\otimes\Pi_Z+I_X\otimes(I_Z-\Pi_Z),\qquad
\Pi_Z=|0^M\rangle\langle0^M|,\qquad W=W_TW_E.
```

This is the complete-operator factorization proved in
[R66, Section 2](SUBLINEAR_T_DEPTH.md#2-protect-one-terminal-logical-block).
Before using Z, compute $`H=[Z=0^M]`$ into the sole external clean
qubit using an exact returned-helper predicate. Do not retest Z while
it is work. On H=1 the registers T,Y,D,A can occupy distinct parts of Z;
on H=0 they are arbitrary.

Use the [packed tagged source](SUBLINEAR_T_DEPTH.md#3-compress-the-complete-columns)
at tolerance $`\delta_0`$. Its ideal action on every input column is

```math
R|j\rangle=W_h|j\rangle_X|t(j)\rangle_T|0\rangle_A,
\quad t(0)=0,\quad t(p1\,0^{h-d-1})=d+1.
```

It first performs the exact permutations
$`T\mathrel{\oplus}=Ht(X)`$ and
$`X\mathrel{\oplus}=Hm(T)`$, where
$`m(d+1)=2^{h-d-1}`$, $`m(0)=0`$, and unused tags have zero mask.
The descendant loader then prepares the literal signed column from
its prefix anchor. Its two tables have Q rows with s outputs and
$`2^s`$ rows with h outputs, respectively, where $`Q=2^h`$.
All subsequent source tables have $`Q\,\mathrm{poly}(n)`$ rows.
No second initialized h-bit index is copied. The whole source is
exactly identity at H=0 on its complete occupied input space and has
active operator-isometry error at most $`\delta_0`$.

Every ideal source column lies in its indicated $`\Pi_{t(j)}`$.
Equations (2)–(3) therefore erase T, while retaining the entire signed
frame column and returning Y,D,A. Applying the actual predicate
inverse after $`\widehat E\widehat R`$ gives an early compiler with

```math
\|V_EJ_1-J_1(W_E\otimes I_b)\|
\le\delta_0+8s\delta.
\qquad\text{(5)}
```

Indeed the ideal word restores Z, so the final predicate inverse
returns H. The same actual inverse propagates all errors unitarily,
including any nonzero final H. The inactive branch is exact identity;
orthogonality of the saved activity sectors gives a maximum, not a sum
over logical columns. Dirty systems and references tensor throughout.

Write $`\eta_L=2^{-L}\le\eta`$. Choose the common tolerance

```math
\delta_0=\delta=\frac{\eta_L}{64(s+1)},\qquad
p=1+\lceil\log_2(1/\delta)\rceil.
```

Then (5) is at most $`\eta_L/8`$, and at fixed accuracy
$`p=O_\eta(\log\log(n+4))`$. The simultaneous conditional clean
allocation is exactly three s-bit records T,Y,D, one reused
$`O(p)`$-bit source bank A, and at most a fixed number of additional
bits already charged by the source interface. Only H is external.

## 6. Native resource ledger and protected-block choice

This section uses the conservative query bound, so it does not rely
on an improvement hidden in a query constant. For an R-row, m-output
table with $`a=\lceil\log_2R\rceil`$, the complete arbitrary-dirty
query has

```math
T,w=O(m(a+2)^3\sqrt R),\qquad
G=O(mR+m(a+2)^3\sqrt R),\qquad D_T=O(\log_2(a+2)).
```

Here w includes outputs when necessary by enlarging the same bound.
The [source ledger](COHERENT_DESCENDANT_LOADER.md#6-native-resource-interface)
has expansion length $`J=2^q=O(p)`$, two r-bit program/gradient
registers with $`r=\log_2(1/\delta)+O(1)`$, and one attenuation bit.
Thus its complete fresh bank has

```math
u=2r+q+1=O(p),\qquad 2^u\le C\delta^{-2}p.
```

Fixed extra bits change only C. Calibration has a constant number of
appearances; its individual rotation tolerance is a constant multiple
of $`\delta`$. Dividing preparation error among r parallel rotations
changes their synthesis counts, not r or this register bound.

| Completed operation | Sufficient rows | Output width or other cost |
|---|---|---|
| Descendant phase/calibration queries | $`C(h+1)Qp`$ | At most $`Cp`$ outputs |
| Central all-work anchor phase | $`C(h+1)Qp\delta^{-2}`$ | One output |
| Variable-suffix Hadamard masks | $`C(h+1)`$ | h outputs, constant-depth private controlled-H bank |
| Initial tag packing | Q and $`2^s`$ | s and h outputs, respectively |
| Local gradient/coefficient preparations | No large table | Polynomial in p gates; depth $`O(p)`$ |
| Search threshold arithmetic | No table | $`O(s^2)`$ Clifford gates |
| Guarded tag XOR | No table | s exact Toffolis, conservatively serial depth $`O(s)`$ |

For example the first row contributes at most
$`Cp(h+p+2)^3\sqrt{(h+1)Qp}`$ to T and w. The second contributes
at most $`C(h+p+2)^3\delta^{-1}\sqrt{(h+1)Qp}`$.
The short masks cost at most
$`Ch(\log_2(h+2)+2)^3\sqrt{h+1}`$; native controlled-H banks,
returned predicates, and small preparations have smaller polynomial
cost. The same substitutions bound G by the table volume plus these
costs. Consequently a safe envelope for one source or reflected-source
macro is

```math
\begin{aligned}
T,w&\le C(h+p+2)^6\delta^{-1}(\sqrt Q+1),\\
G&\le C(h+p+2)^6\delta^{-2}(Q+1),\\
D_T&\le C[\log_2(h+2)+p].
\end{aligned}
\qquad\text{(6)}
```

Every query pool returns exactly after its completed macro. Its peak
width is reused, not added over the sequence. Formula (6) includes
actual inverses and local gradient preparation on each call. The
central reflection tests that gradient as part of A; no shared-gradient
invariant is assumed.

There is one initial tagged source, $`4s`$ descendant-source appearances
and $`2s`$ central phases in the decoder, and the exact arithmetic above.
The initial and final saved-predicate words cost polynomially in M
and have $`O(\log_2(M+2))`$ T-depth. With fixed accuracy a deliberately
loose complete early envelope is

```math
\begin{aligned}
T_E,w_E&\le C_\eta(n+2)^{10}(\sqrt Q+1),\\
G_E&\le C_\eta(n+2)^{10}(Q+1),\\
D_{T,E}&=O_\eta(\ell_n^2).
\end{aligned}
\qquad\text{(7)}
```

Indeed $`s,\delta^{-1},p\le C_\eta(n+2)`$ in the coarse count
bound, while the sharper $`s=O(\ell_n)`$ and
$`p=O_\eta(\log\log(n+4))`$ give the depth. The degree ten
majorant covers all rows of the table and both predicate boundaries.

Choose a fixed constant $`A_0\ge24`$, increased if necessary for
the simultaneous clean-work constant, and put

```math
M=\lceil A_0\ell_n\rceil,\qquad h=n-M,\qquad
Q=2^h\le N(n+2)^{-A_0}.
```

For sufficiently large n, $`h\ge n/2`$ and T,Y,D,A fit inside M.
Equations (6)–(7) imply

```math
\begin{aligned}
T_E,w_E&\le C_\eta\left[\frac{\sqrt N}{(n+2)^2}+(n+2)^{10}\right],\\
G_E&\le C_\eta\left[\frac{N}{(n+2)^{14}}+(n+2)^{10}\right].
\end{aligned}
\qquad\text{(8)}
```

This is why the query overhead is affordable: all large tables address
the h-bit early frame, not all n original bits. The protected block
has only logarithmic length, allowing the direct terminal construction.

## 7. Terminal layers and the all-width splice

The remaining M original layers use the
[one-clean smaller-error bottom construction](ONE_CLEAN_DEPTH.md#the-final-logical-layers-at-smaller-error),
with $`R_0=M`$. No shared terminal stream is needed. For
$`j=1,\ldots,M`$ counted upward from the bottom, its parameters are

```math
q_j=L+j+7,\qquad m_j=L+j+8,\qquad Q_j=N2^{-j}.
```

The core, original selectors and helper reserve exactly
$`B'_0=L+n+9`$ dirty bits. At the supplied threshold,
$`B=b-B'_0\ge b/2`$. If $`b/2\ge16c_I(n+2)`$, the shallow
rectangular query schedule gives

```math
\begin{aligned}
T_T&=O\!\left(\sqrt{NL}+NL/b+ML+n^2\right),\\
G_T&=O(NL),\\
D_{T,T}&=O\!\left(NL/b^2+M[L+M+\ell_n]\right).
\end{aligned}
\qquad\text{(9)}
```

Its error is less than
$`\sum_{j=1}^M30\,2^{-(L+j+7)}\lt15\eta_L/64\lt\eta_L/4`$.
The geometric sums of $`Q_jm_j`$ and $`\sqrt{Q_jm_j}`$ give
the first two lines of (9); the source and predicate depths sum to
the displayed last line. At fixed L and $`M=O(\ell_n)`$, this is
the claimed count and $`O_\eta(N/b^2+\ell_n^2)`$ depth. In
particular $`n^2=O(\sqrt N)`$ absorbs its polynomial count overhead.

Reuse H as the terminal signal after its early predicate inverse.
If $`R_E=W_E\otimes I_b`$ and $`R_T=W_T\otimes I_b`$, then

```math
V_TV_EJ_1-J_1R_TR_E
=V_T(V_EJ_1-J_1R_E)+(V_TJ_1-J_1R_T)R_E.
```

This bounds total error by less than $`3\eta_L/8`$. It propagates
actual earlier leakage, with no reset or assumed exact clean return.
Early and terminal dirty pools are sequential and may use the same
supplied wires.

To cover every eligible b, set

```math
b_c=\frac{\sqrt N}{n+2}.
```

Choose a fixed cutoff $`n_0(\eta)`$ large enough that, for
$`n\ge n_0`$, the above logical reservations fit, (8) gives
$`w_E\le b_c/2`$, and
$`b_c/2\ge\max\{L+n+9,16c_I(n+2)\}`$.
The exponential growth compared with every displayed polynomial
establishes all these conditions simultaneously.

- **If $`b\ge b_c`$**, use the new early circuit and (9). Their
  peak dirty use fits b, and all stated error, count and depth bounds
  hold on this one complete circuit.
- **If $`b\lt b_c`$**, use the
  [original one-clean all-layer fallback](ONE_CLEAN_DEPTH.md#3-a-one-clean-fallback-at-every-width).
  It has depth $`O(NL/b^2+n(L+n))`$, count
  $`O(\sqrt{NL}+nL+NL/b)`$, and $`G=O(NL)`$ at the literal
  threshold. Here $`N/b^2\gt(n+2)^2`$ absorbs the additive
  depth at fixed L; $`nL=O_\eta(\sqrt N)`$ absorbs its count
  overhead. Its error is less than $`15\eta_L/16`$.
- **If $`n\lt n_0(\eta)`$**, use that same fallback at every
  eligible b. Its finite additive depth is bounded by
  $`C_\eta\ell_n^2`$, uniformly in b and the angle tuple.
  The same count estimate holds. This is an actual native compiler
  fallback, not a finite net or an assertion of exact synthesis.

These cases prove the real theorem at every n and literal eligible
allocation. No efficient numerical crossover is claimed.

## 8. Complex phases, lower bounds, and limitations

For $`D_\phi W_{\mathbb R}`$, use the real theorem with error
$`2^{-(L+1)}`$ and compose the
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at the same half-error. Both fit $`b\ge17(L+n+8)`$. Its additional
resources at fixed accuracy are

```math
T_D=O_\eta(\sqrt N+N/b),\qquad G_D=O_\eta(N),\qquad
D_{T,D}=O_\eta(N/b^2+\ell_n).
```

The full initialized-isometry hybrid from Section 7 reuses the same
single external signal and dirty pool. Errors add to at most
$`2^{-L}\le\eta`$; counts and depths add within the theorem's
bounds. The specified leaf phases are retained literally, including
at zero-amplitude leaves. This proves the complex statement.

The [inherited complete-frame and diagonal lower bounds](ONE_CLEAN_DEPTH.md#8-matching-intervals-and-evidence)
give at fixed accuracy

```math
T^\star=\Omega_\eta(\sqrt N+N/b),\qquad
D_T^\star=\Omega_\eta(N/b^2),
```

because physical width $`n+b+1=\Theta(b)`$. At
$`b\le\sqrt N/\ell_n`$, the width term absorbs the additive
$`\ell_n^2`$ depth and $`N/b`$ absorbs $`\sqrt N`$.
This proves the stated matching interval, without a growing lower
bound at square-root width.

The mechanism does not improve the entry-indexed normalization in
[R67](ENTRY_INDEXED_ENCODING_LIMIT.md). It replaces that overlap
recovery architecture with target-dependent, natively compiled
projector queries. It also does not infer cheap recovery from a
general compressor definition: Sections 3–7 construct and price the
entire decoder. The count theorem, all-precision R72, and the separate
R74–R77 constant-layer obstructions retain their scopes.

## 9. Evidence and proof dependencies

The proof uses the literal addressed-frame identities, the original
coherent descendant source, complete arbitrary-dirty shallow queries,
returned predicates, and the one-clean bottom/fallback/diagonal words.
The new dependencies are the nested-projector identity, full-mode
source-to-reflection comparison, reversible spectral search, and the
fixed-accuracy width splice proved here. It does not use the shared
gradient, shared terminal stream, or old outer amplitude amplification.

The bounded checks in
[`test_subtree_projector_decoder.py`](../../tests/test_subtree_projector_decoder.py)
audit the new algebra and its failure-sensitive phase/work conventions.
Finite matrices and verification-only unitary completions are not
native loader implementations or proofs of asymptotic bounds. The
analytic native interfaces above supply those bounds. No general
scalable emitter, efficient preprocessing, practical crossover, or
external novelty claim follows from this result.
