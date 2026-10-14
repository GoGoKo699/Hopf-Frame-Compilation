# Coherent signed descendant preparation with conditional workspace

This chapter supplies the state-source part of the
[sublinear-depth theorem](SUBLINEAR_T_DEPTH.md). Its native
resource statement uses the separately charged
[shallow dirty-query construction](SHALLOW_DIRTY_QUERIES.md). It is not
an unrestricted T-depth lower bound or, by itself, a complete-frame
depth theorem.

[Complete frame and marker convention](../../docs/HOPF_INTERFACE.md) ·
[Conditional-work contract](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally) ·
[Existing literal routing and predicates](../../docs/T_DEPTH_COMPILER.md)

The source is based on the Boolean-phase state-preparation method of
[Gosset, Kothari, and Wu, *Quantum State Preparation with Optimal
T-Count*, arXiv:2411.04790v3, Lemma 3.2 and Section 3.3](https://arxiv.org/html/2411.04790v3#S3).
The changes below make the construction coherent in a preserved source
label, replace its Boolean-oracle implementation by returned dirty
queries, and give a strictly certifiable residual recursion. No
classical preprocessing running-time bound is asserted. Every quantum
query, inverse, small rotation, guard, and reflection is charged.

## 1. Source contract and disjoint descendant sectors

Let an h-bit data register have computational dimension $`Q=2^h`$.
A preserved tag k specifies a prefix length $`d(k)`$. Within tag k,
the first $`d(k)`$ data bits form the preserved prefix p; the remaining
$`s(k)=h-d(k)`$ bits are the local target register. For every allowed
pair $`a=(k,p)`$, let $`|\psi_a\rangle`$ be an effectively specified
real unit vector on those local bits. Its signs and literal phase are
part of the target.

The source receives local input $`|0^{s(k)}\rangle`$ and a separate
zero scratch register A. It preserves k and p exactly. At fixed k,
different p values therefore occupy disjoint invariant sectors, and
different tags are orthogonal. The desired column action is

```math
|k\rangle|p\rangle|0^{s(k)}\rangle|0^A\rangle
\longmapsto
|k\rangle|p\rangle|\psi_{k,p}\rangle|0^A\rangle.
```

In the frame application, tag zero denotes the state column. A marker
at depth d receives tag $`d+1`$ and retains its d-bit node prefix.
The signed descendant vector is its complete marker column restricted
to that subtree. The marker-to-tag permutation is separate from this
source; the source does not copy an h-bit input index into clean work.

An external activity bit z is computed before the reserved suffix is
used as scratch. On $`z=1`$, the designated source scratch is zero.
On $`z=0`$, it is arbitrary, possibly entangled with the data and
references. Every completed native macro below is exactly identity
on $`z=0`$ on its full input space. The activity bit is never modified.
All supplied dirty query work is arbitrary on both sectors and returns
exactly after each completed query. Approximate source return, including
leakage into A, is included in the source-isometry error.

The scalar case $`s(k)=0`$ is included. A real unit vector is then
$`+1`$ or $`-1`$. Its phase is a Boolean phase on the preserved label;
there are no local Hadamards. No nonexistent target bit is borrowed.

## 2. A strictly certifiable geometric expansion

The inherited coarse-state lemma says that every real unit vector u on
s qubits admits Boolean diagonal matrices $`B_1,B_2`$ such that

```math
|\phi\rangle=B_2H^{\otimes s}B_1H^{\otimes s}|0^s\rangle,
\qquad \langle u|\phi\rangle\ge\frac1{\sqrt2}.
```

The following fixed-coefficient recursion avoids normalizing a zero
residual or deciding whether a computable residual is exactly zero.
Put

```math
\gamma=\frac35,\qquad \beta=\frac45,\qquad r_0=\psi.
```

At step j, enumerate Boolean pairs and accept any resulting unit vector
$`\phi_j`$ with a certified strict inequality

```math
\|r_j-\gamma\beta^j\phi_j\|\lt\beta^{j+1}.
```

Then set $`r_{j+1}=r_j-\gamma\beta^j\phi_j`$. The search dovetails
all finitely many candidates with increasing arithmetic precision; it
does not wait to decide one possibly false or equal comparison.

**Existence and termination.** Suppose $`\|r_j\|\le\beta^j`$.
If the residual is nonzero, the coarse lemma supplies overlap at least
$`1/\sqrt2\gt2/3`$ with its normalized direction. Writing
$`x=\|r_j\|/\beta^j\in[0,1]`$, that candidate satisfies

```math
\begin{aligned}
\frac{\|r_j-\gamma\beta^j\phi_j\|^2}{\beta^{2j}}
&\le x^2+\frac9{25}-\frac45x\\
&\le\max\left\{\frac9{25},\frac{14}{25}\right\}
 =\frac{14}{25}\lt\frac{16}{25}=\beta^2.
\end{aligned}
```

The quadratic is convex, so its maximum on the interval is at an
endpoint. At zero residual every candidate instead has norm
$`\gamma\beta^j\lt\beta^{j+1}`$. A strictly passing candidate
therefore always exists and eventually receives a finite certificate.
Induction gives

```math
\left\|\psi-\gamma\sum_{j=0}^{J-1}\beta^j\phi_j\right\|
\lt\beta^J.
```

Choose a common power of two J for every source label, with
$`2\beta^J\le\epsilon_{\rm data}`$. Put

```math
S_a=\sum_{j=0}^{J-1}\beta^j\phi_{a,j},\qquad
\Phi_a=\frac{S_a}{\|S_a\|},\qquad
A_J=\sum_{j=0}^{J-1}\beta^j=\frac{1-\beta^J}{1-\beta}.
```

The sum is nonzero, since $`\|\gamma S_a-\psi_a\|\lt\beta^J\lt1`$.
Normalization changes the vector by at most its original error, giving

```math
\|\Phi_a-\psi_a\|\le2\beta^J\le\epsilon_{\rm data}.
```

Here $`J=O(\log(1/\epsilon_{\rm data}))`$ and the index uses
$`q=\log_2J=O(\log\log(1/\epsilon_{\rm data}))`$ bits. The
Boolean phase tables and finite geometric sums are computable. In
particular, the subsequent calibration does not require a sign or
zero test on an unknown analytic amplitude.

## 3. Constant-query flagging and literal inner amplification

The normalized coefficient state factors into one-qubit real rotations:

```math
|\chi\rangle
=\frac1{\sqrt{A_J}}\sum_{j=0}^{J-1}\beta^{j/2}|j\rangle
=\bigotimes_{\ell=0}^{q-1}
\frac{|0\rangle+\beta^{2^\ell/2}|1\rangle}
     {\sqrt{1+\beta^{2^\ell}}}.
```

Prepare this index state; apply the local Hadamards, the Boolean phase
$`B_{1,a}^{(j)}`$, the local Hadamards, and $`B_{2,a}^{(j)}`$;
then reverse the coefficient preparation. All controlled Boolean
phases combine into two tables on the index, tag, and data. The
accepted zero-index component is

```math
\frac{S_a}{A_J}=\xi_a\Phi_a,
\qquad \xi_a=\frac{\|S_a\|}{A_J}.
```

The triangle inequality gives $`\xi_a\le1`$, and the residual bound
gives a uniform lower bound:

```math
\xi_a\ge
\frac{(1-\beta^J)/\gamma}{(1-\beta^J)/(1-\beta)}
=\frac{1-\beta}{\gamma}=\frac13.
```

Add one zero attenuation bit. A positive real Y rotation, controlled by
the preserved label a, maps its zero input to

```math
\frac{\sin(\pi/10)}{\xi_a}|0\rangle
+\sqrt{1-\frac{\sin^2(\pi/10)}{\xi_a^2}}|1\rangle.
```

This is well conditioned because $`\sin(\pi/10)\lt1/3`$. The
rotation is label dependent; it cannot be replaced by one universal
angle. Its angle is computed from the finite Boolean states and the
fixed rational geometric weights.

Let $`V_a`$ denote this ideal flagging word, including attenuation,
and let A now denote its index and attenuation bits. Its column is

```math
V_a|0^A0^{s(k)}\rangle
=\sin(\pi/10)|0^A\rangle|\Phi_a\rangle+|\tau_a\rangle,
\qquad (|0^A\rangle\langle0^A|\otimes I)\tau_a=0.
```

Define within each preserved sector

```math
R_{0,a}=2|0^A0^{s(k)}\rangle\langle0^A0^{s(k)}|-I,
\qquad R_A=2|0^A\rangle\langle0^A|\otimes I-I,
\qquad R_{1,a}=V_aR_{0,a}V_a^\dagger.
```

Two amplitude-amplification rounds give, with rightmost factors acting
first,

```math
(R_{1,a}R_A)^2V_a|0^A0^{s(k)}\rangle
=|0^A\rangle|\Phi_a\rangle.
```

Indeed the usual amplification step is $`-R_{1,a}R_A`$, and its
minus sign squares to one. The final success amplitude is
$`\sin(5\pi/10)=1`$, with the displayed positive phase. The
expanded word contains five occurrences of V or its actual inverse
and a constant number of reflections. Thus the construction has a
constant number of Boolean phase calls, not J sequential calls.

The reflected zero input consists only of the local target and source
index/attenuation registers. The preserved prefix is not tested for
zero. For a zero-length target it tests A alone.

## 4. Native guarded operations and their return

### Boolean phases and variable suffix Hadamards

An exact arbitrary-output XOR query $`Q_f`$ and one arbitrary dirty
bit c implement the Boolean phase by the chronological echo

```math
Q_f,\ Z_c,\ Q_f^\dagger,\ Z_c.
```

Its phase is $`(-1)^{f}`$ independently of c. Include z in the query
address and set every $`z=0`$ row to zero. The completed word is then
exactly identity on the inactive sector, including arbitrary source
scratch, dirty inputs, and references.

For a local Hadamard bank, query the mask
$`m_i=z[i\ge d(k)]`$ into arbitrary dirty bits $`c_i`$. Apply the
controlled-H gates $`\mathrm{c}H(c_i;t_i)`$ on the distinct target
wires, reverse the mask query, and repeat the same controlled-H bank.
The dirty offsets cancel because H is an involution and the different
target operations commute. The mask depends only on z and k, which
are preserved. The completed bank applies exactly the selected suffix
Hadamards and is identity at $`z=0`$ on all inputs.

The exact controlled-H implementation and cancellation of its scalar
phases are given in the
[conditional geometric source](../../docs/CONDITIONAL_GEOMETRIC_SOURCE.md#native-phases-and-an-explicit-reservation).
All gates in a bank have disjoint control-target pairs, so its non-query
T-depth is constant. The dirty masks are distinct from targets and all
live source/program registers.

### Guarded real rotations and the calibration table

Every continuous operation above is a real rotation
$`R_y(\theta)=e^{-i\theta Y}`$. Let p be a preserved program bit,
t its target, and $`Z_c=\mathrm{CCZ}(z,p;t)`$. If a native word V
approximates $`R_y(\theta/2)`$ with literal operator error e, use the
chronological sequence

```math
V_t,\ Z_c,\ V_t^\dagger,\ Z_c.
```

On $`z=p=1`$, its ideal matrix is

```math
Z R_y(-\theta/2) Z R_y(\theta/2)=R_y(\theta).
```

Its error is at most $`2e`$. If either control is zero, the actual
word is exactly $`V^\dagger V=I`$, on every target input. Its
inverse is the actual reversed native word. For a fixed coefficient
rotation, replace CCZ by $`\mathrm{CZ}(z,t)`$. No controlled-T
primitive, initialized phase state, or row-dependent discarded scalar
is assumed.

For each a, choose an r-bit approximation to the calibration angle,
with circular error at most $`2\pi2^{-r}`$. Load its binary digits
into a separate conditionally zero r-bit program register, apply the
preceding echo to every angle bit, and reverse the actual query.
The address depends only on the preserved tag and prefix, not on the
rotated attenuation bit. The program is used only as controls.
Consequently it returns exactly to zero on active initialized inputs;
the supplied dirty query pool returns exactly on all inputs. At
$`z=0`$, the query and every echo are exactly identity, even if the
program is occupied.

If each half-angle word has error $`e_{\rm syn}`$, the complete
active initialized-program error is at most

```math
2\pi2^{-r}+2r e_{\rm syn}.
```

For target error $`\epsilon_{\rm rot}`$, take
$`r=O(\log(1/\epsilon_{\rm rot}))`$ and
$`e_{\rm syn}=O(\epsilon_{\rm rot}/r)`$. The small-gate
contribution to T, G, and T-depth is

```math
O\!\left(r\log(r/\epsilon_{\rm rot})\right).
```

The fixed coefficient rotations have smaller polynomial-logarithmic
cost and need no calibration program. Their native synthesis uses the
retained [one-qubit synthesis premise](../../docs/FAULT_TOLERANT_COMPILER.md#primary-references).
The same physical program register is reused only after its exact return.

On an occupied active program input, load/rotate/unload implements the
XOR-offset program. It is not asserted to equal the desired calibrated
rotation tensored with identity on that program. The initialized-source
error estimate and the linked full-norm completion respect this
distinction.

### Reflections and exact inactive phases

Every zero-state reflection is a Boolean phase depending on the
preserved tag, relevant local target bits, and source index. It can
therefore use the exact guarded Boolean phase query. The native word
for $`I-2z\Pi`$ uses the Boolean function $`z[\Pi]`$. When the
ideal formula calls for $`2\Pi-I`$, include $`Z_z`$ to change only
the active scalar sign. The resulting macro is still exactly identity
at $`z=0`$. Both the scalar signs in the inner amplification and all
actual inverses are retained.

Fixed zero tests may alternatively use the returned-helper
[logarithmic-depth predicate](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
Variable suffix zero tests can be folded into the Boolean table. Neither
construction assumes a clean dirty helper or a zero preserved prefix.

## 5. Coherent error and full-unitary comparison

Let $`A_{\rm all}`$ include the index, attenuation, and calibration
program registers. The active calibration load/rotate/unload preserves
every computational program value exactly: the address is unchanged,
program bits are controls only, and the actual query inverse removes
the loaded word. Thus its zero-program subspace is exactly invariant
for every target/index input and under actual forward and inverse
words. The other source macros do not change that register. Native
rotation error comparisons and unitary telescoping may therefore be
performed on this invariant subspace, without asserting the desired
row rotation on occupied programs.

Choose the data and native rotation errors so that the completed native
source has, for every allowed a,

```math
\|v_a-|\psi_a\rangle|0^{A_{\rm all}}\rangle\|\le\delta,
```

where $`v_a`$ is its actual normalized output, including any leakage
into index, attenuation, and program scratch. A constant number of
native flagging words and reflections occur, so ordinary unitary
telescoping bounds the native error by a constant multiple of the
allocated rotation errors. The data error contributes separately.

The native circuit preserves tag and prefix exactly. All complete dirty
query words return their dirty registers exactly, so the mode-space
unitary is tensored with identity on the supplied dirty pool. Hence
orthogonality of the preserved sectors gives the full source-isometry
bound

```math
\left\|\widetilde UJ_{\rm src}
      -J_{\rm out}\Psi\right\|\le\delta.
```

This bound holds for every superposition of source anchors, arbitrary
dirty input, and arbitrary references. There is no factor equal to the
number of columns. Source scratch return is part of the same norm.

The dimension-free isometry-extension lemma in the
[complete depth construction](SUBLINEAR_T_DEPTH.md#4-from-a-source-error-to-an-honest-unitary-comparison)
converts this actual source-isometry bound into a full-mode comparison
for outer recovery. Its comparison unitary may act on leaked source
scratch; no implementation of that mathematical correction is assumed.
In particular, the source-isometry estimate is not misidentified with a
specified full-space action on occupied program registers.

## 6. Native resource interface

Put $`p=1+\lceil\log_2(1/\delta)\rceil`$. The source uses a
common $`J=O(p)`$ expansion length, an $`O(\log p)`$-bit index,
one attenuation bit, and a reused $`O(p)`$-bit calibration program.
Its additional conditionally initialized footprint is $`O(p)`$;
the tag itself is separate. All large query banks and mask outputs are
arbitrary dirty qubits.

For the h-bit descendant family, each tag has at most Q data strings.
With $`O(h+1)`$ tags and J expansion indices, a sufficient padded
Boolean address space has

```math
S=O((h+1)Qp)
```

rows. Calibration output width is $`O(p)`$. The variable-suffix
Hadamard masks have only $`O(h+1)`$ short-tag rows and h outputs;
their costs are polynomial in h and are separately included. Tables
for invalid tags may be set to identity without changing any source
column.

There are a constant number of large phase queries, calibration queries,
their actual inverses, and reflection queries. Denote by
$`T_{\rm q},G_{\rm q},D_{\rm q},b_{\rm q}`$ the charged resources
of the exact [shallow dirty-query primitive](SHALLOW_DIRTY_QUERIES.md)
for these padded row and output sizes. The source has

```math
\begin{aligned}
T&=O(T_{\rm q}+\mathrm{poly}(h+p)),\\
G&=O(G_{\rm q}+\mathrm{poly}(h+p)),\\
D_T&=O(D_{\rm q}+p^2+\log(h+2)),\\
b&=O(b_{\rm q}+\mathrm{poly}(h+p)).
\end{aligned}
```

The polynomial terms include short-tag mask queries, exact returned
predicate helpers, and the native small rotations. For an R-row,
m-output query, the conservative standalone query theorem gives
$`T_{\rm q},b_{\rm q}=O(m(r+2)^3\sqrt R)`$,
$`G_{\rm q}=O(mR+m(r+2)^3\sqrt R)`$, and
$`D_{\rm q}=O(\log(r+2))`$, where $`r=\lceil\log_2R\rceil`$.
Substituting the actual row and output sizes above yields a fixed
polynomial P such that

```math
T,b\le P(h,p)2^{h/2},\qquad
G\le P(h,p)2^h,\qquad
D_T=O(\log(h+2)+p^2).
```

The additional conditional-zero scratch is $`O(p)`$ beyond the
separately stored $`O(\log(h+2))`$-bit tag. At
$`\delta^{-1}=\mathrm{poly}(n)`$ and $`h\le n`$, the complete
source depth is therefore $`O(\log^2 n)`$. The final early/tail
allocation must still absorb the displayed polynomial count and width
factors; they are not omitted from the native ledger.

The complete source proof uses no measurement, postselection, reset,
QRAM, supplied magic state, or catalyst. Approximate native rotations
are not required to preserve exact real amplitudes; their effect and
all resulting work leakage are covered by the isometry estimate.
