# Literal leaf phases preserve the all-precision depth bound

[Real-frame theorem, R70](ALL_PRECISION_DEPTH.md) · [Packed literal diagonal](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames) · [Complex magnitude convention](../../docs/HOPF_INTERFACE.md#6-phase-dressed-complex-magnitude-frame)

**Local result R71.** The prescribed phase-dressed complex magnitude
frame has the same all-precision count/depth bounds as the real frame,
with one additional unit in the precision reservation. The new local
step is a native depth schedule for the existing packed literal diagonal.
It replaces complete masks by exact rectangular queries and runs a
constant number of linear-precision source words serially. A routed
indicator fallback retains the literal dirty-width threshold when the
shallow indicator's sufficient width constant does not fit.

The [one-clean refinement](ONE_CLEAN_DEPTH.md), R72, retains these
bounds and literal threshold with a different native diagonal word.
This chapter preserves the two-clean packed construction and its proof.

## 1. Complete complex-frame theorem

Let $`n\ge1`$ and $`L\ge6`$ be integers, put $`N=2^n`$, and define

```math
F(n,L)=\min\left\{nL,(\sqrt n+L)\log_2(n+L+2)\right\}.
```

For the prescribed real magnitude frame and independently specified
leaf phases, write

```math
W_{\mathbb C,\mathrm{mag}}=D_\phi W_{\mathbb R},\qquad
D_\phi=\sum_{x=0}^{N-1}e^{i\phi_x}|x\rangle\langle x|.
```

Every integer allocation

```math
a=2,\qquad b\ge17(L+n+8)
```

admits one coherent Clifford+T circuit V such that

```math
\|VJ_2-J_2(W_{\mathbb C,\mathrm{mag}}\otimes I_b)\|\le2^{-L},
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+F(n,L)\right).
```

All constants are absolute and all three resources belong to the same
circuit. The embedding $`J_2`$ initializes only the two external flags;
logical and dirty inputs are unrestricted. The norm includes every
logical column, literal common and relative phases, clean or logical-work
leakage, dirty return, and arbitrary external reference correlations.
Unused supplied dirty wires remain untouched. G counts elementary
Clifford gates, including those between T layers.

Throughout the high-precision range there is also the stronger bound

```math
L\gt N/n^2
\quad\Longrightarrow\quad
D_T=O\!\left(\frac{NL}{b^2}+L\right),
```

with the same count, Clifford-count, width, and error guarantees.
The T-count is worst-case optimal at every eligible allocation. Section 6
states the matching depth interval; neither additive depth term is
claimed necessary at unrestricted width.

Real angles and leaf phases (or the corresponding real and imaginary
diagonal entries) admit terminating certified evaluation. Phases at
zero-amplitude leaves are independently specified as part of the operator.
No efficient classical
preprocessing bound is asserted. The target is the specified
phase-dressed magnitude completion, including singular and unequal
tuples, rather than an arbitrary complex unitary.
No measurement, reset, QRAM, supplied precision state, or uncharged
coherent evaluator is used.

## 2. The existing packed diagonal and its exact native word

We first prove the following standalone statement. For every integer
$`\ell\ge6`$ and $`b\ge17(\ell+n+7)`$, a prescribed literal diagonal
has a two-clean circuit $`V_D`$ with

```math
\|V_DJ_2-J_2(D_\phi\otimes I_b)\|\le2^{-\ell},
```

```math
T_D=O\!\left(\sqrt{N\ell}+\ell+\frac{N\ell}{b}\right),\qquad
G_D=O(N\ell),
```

```math
D_{T,D}=O\!\left(\frac{N\ell}{b^2}+\ell+\log_2(n+2)\right).
```

Its source, coefficient rule, and recovery are exactly the
[packed diagonal construction](../../docs/OPERATOR_SOURCE_COMPILER.md#8-literal-diagonal-unitaries-and-phase-dressed-frames).
Only the complete query implementations change. Set

```math
q=\left\lceil\frac{\ell+4}{2}\right\rceil,\qquad
r=n+1,\qquad Q=2^r=2N.
```

Thus $`q=\Theta(\ell)`$ and $`q\le\ell`$; the source precision
does not grow with n. On q arbitrary dirty core qubits the packed
source uses $`2q`$ paired Majoranas and an exact involution

```math
M=UX_0U^\dagger=\sum_{j=0}^{2q-1}a_j\gamma_j,\qquad
a_j^2=\begin{cases}
2^{-(j+1)},&j\lt2q-1,\\
2^{-(2q-1)},&j=2q-1.
\end{cases}
```

Every adjacent splitting rotation in U is one T gate conjugated by a
fixed Clifford, up to a scalar that cancels in the actual conjugation.
The source and its centrally controlled version therefore have at most
$`4q-2`$ T gates and $`O(q)`$ Cliffords, with no initialized core or
source-private helper. Their serial T-depth is $`O(q)`$.

For a sign word f, the inherited literal mask has the order
$`P_f=X^{u(f)}Z^{v(f)}`$ and satisfies

```math
N_f=P_fMP_f^\dagger,\qquad
\frac{MN_f+N_fM}{2}=c_fI,\qquad
c_f=\sum_ja_j^2(-1)^{f_j}.
```

These identities hold on the entire dirty input space. For each
address x, certified rounding chooses cosine and sine tables with
coefficient error at most $`(5/4)2^{1-2q}`$. Write the two clean
flags as $`a_0,a_1`$, with $`a_1`$ selecting the cosine or sine row.
The actual scalar word and literal diagonal block are

```math
\mathcal S_f
 =H_{a_0}C_{a_0=0}(M)N_fC_{a_0=1}(M)H_{a_0},\qquad
Q_D=H_{a_1}S_{a_1}\mathcal S_{f(a_1,x)}H_{a_1}.
```

Products here act rightmost first. The query address is the preserved
tuple $`(a_1,x)`$. The S gate supplies the literal factor i on the
sine branch, giving

```math
B_D:=J_2^\dagger Q_DJ_2
 =\frac12\,\mathrm{diag}(c_x+i s_x)\otimes I_b,
```

```math
\zeta:=\|2B_D-D_\phi\otimes I_b\|
 \le\frac{5\sqrt2}{4}\,2^{1-2q}\le\frac14.
```

Use the actual normalization-two word

```math
R=I-2J_2J_2^\dagger,\qquad
V_D=-Q_DRQ_D^\dagger RQ_D.
```

R tests only the two flags and is Clifford. The leading minus is the
literal Clifford scalar $`XZXZ=-I`$. The
[full-isometry amplification estimate](../../docs/OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
then gives

```math
\|V_DJ_2-J_2(D_\phi\otimes I_b)\|
 \le4\zeta\le5\sqrt2\,2^{1-2q}
 \le2^{4-2q}\le2^{-\ell}.
```

In particular, core disturbance and rejected flag amplitudes are
included; no accepted scalar is silently normalized.

Each $`Q_D`$ contains three source appearances and four complete
q-bit XOR/Z queries: two for $`P_f`$ and two for its actual inverse.
Consequently $`V_D`$ contains exactly twelve such queries and nine
source appearances, counting its actual inverse call. The other flag
gates, both reflections, and the leading sign are Clifford. No suffix
predicate or long controlled source loader is present.

## 3. Complete masks on a disjoint dirty pool

Reserve

```math
B_0=\ell+n+7,\qquad B=b-B_0.
```

The original packed base $`q+n+1`$ fits in $`B_0`$. It contains the
source core and original selectors, which may remain unused after query
replacement. The supplied threshold implies

```math
B\ge\frac{16}{17}b\ge16B_0.
```

All new indicator outputs, indicator helpers, high-address traversal
selectors, and the controlled-bilinear helper lie in this additional
pool B. They are disjoint from the core, both flags, the logical
address, and every base reservation. Different complete queries reuse
the pool only after its exact return.

For either mask table F, the replacement has the full-space action

```math
|y,z,w\rangle\longmapsto|y,z\mathbin\oplus F(y),w\rangle,
```

with arbitrary output z and helper w. A Z mask is the complete XOR
query conjugated by Hadamards on all q core outputs. Implement
$`P_f=X^{u(f)}Z^{v(f)}`$ by executing the Z mask before the X mask;
its inverse is the actual reversed word. The mask need not be
Hermitian, and its address-dependent Pauli scalar is retained.

All original table rows are unchanged, including rows reached on
occupied flags. No individual query is assumed to receive an
initialized output or helper. The complete replacement equals the old
mask tensored with identity on B, so it preserves the scalar identity,
every rejected block, and the preceding error proof. Query helpers
return exactly on arbitrary inputs and references. The source core's
final return is governed by the full-isometry estimate above.

### Sufficient width for the nonuniform indicator

Let $`c_I\ge1`$ be the fixed constant in the
[nonuniform rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
If

```math
B\ge16c_I(r+2),
```

that exact query, with output width q, gives

```math
T_Q=O\!\left(q+\sqrt{Qq}+\frac{Qq}{B}\right),\qquad
G_Q=O(Qq),
```

```math
D_{T,Q}=O\!\left(q+\frac{Qq}{B^2}+\log_2(r+2)\right).
```

The indicator constant is not assumed to fit at every allocation
$`b\ge17B_0`$. The following explicit fallback covers the rest.

### Routed fallback at the literal threshold

The retained
[routed indicator](../../docs/PARALLEL_DIRTY_LOOKUP.md#1-an-exact-dirty-indicator-by-conjugated-routing)
on s address bits uses exactly its $`2^s`$ arbitrary output wires,
no additional helpers, $`O(2^s)`$ T and Clifford gates, and at most
$`8s`$ T layers. Use it in the existing
[selected bilinear query and two-indicator echo](../../docs/BLOCKED_BILINEAR_LOOKUP.md#2-select-a-bilinear-block-with-dirty-traversal).
Choose

```math
J=\mathrm{pow2floor}\,
 \min\left\{Q,\sqrt{Qq},B/4\right\},\qquad
H=\min\{J,Q/J\},\qquad K=Q/(HJ),\qquad p=\log_2K.
```

Here pow2floor rounds down to a power of two. All H, J, and K are
positive integer powers of two, $`H\le J`$, and $`HJ\le Q`$.
The two indicator outputs use at most $`H+J\le2J\le B/2`$ wires.
The high-address selectors and one bilinear helper use
$`p+1\le r+1=n+2\le B/2`$ more. Thus the literal pool B suffices
without a hidden indicator width constant or initialized copy bank.

The four-call bilinear echo, including both indicator inverses, has

```math
T_Q=O(KqH+H+J),\qquad G_Q=O(Qq+H+J),
```

```math
D_{T,Q}=O(Kq+\log_2H+\log_2J).
```

Power-of-two rounding gives

```math
KqH=Qq/J=O\!\left(q+\sqrt{Qq}+Qq/B\right).
```

Also $`H+J\le2J\le2\min\{Q,\sqrt{Qq}\}`$. If
$`J\ge\sqrt Q`$, then $`H=Q/J`$ and $`K=1`$. Otherwise
$`H=J`$; a count cap gives $`J\gt\sqrt Q/2`$, while a width cap
gives $`J\gt B/8`$. In either case

```math
K=O(1+Q/B^2).
```

The routed fallback therefore has the same count bounds and

```math
D_{T,Q}=O\!\left(q+\frac{Qq}{B^2}+r\right).
```

All coordinate changes, routes, selected bilinear words, traversal
selectors, and their actual inverses are charged by the linked native
identities. No output-bank routing or extra initialized register is
added. As in the shallow branch, the entire query is an exact
phase-free XOR on every input.

## 4. Uniform standalone diagonal resources

On the shallow branch, twelve queries, nine serial sources, and the
fixed Clifford flag words immediately give

```math
T_D=O\!\left(\sqrt{N\ell}+\ell+\frac{N\ell}{b}\right),\qquad
G_D=O(N\ell),
```

```math
D_{T,D}=O\!\left(\frac{N\ell}{b^2}+\ell+\log_2(n+2)\right).
```

The H conjugations for Z masks add only $`O(q)`$ Cliffords per query.
The serial source depth is bounded using its own $`O(q)`$ native word,
not by inferring depth from the entire diagonal's T-count.

On the routed branch, failure of the shallow width condition implies

```math
b\le\frac{17}{16}B\lt17c_I(n+3).
```

Consequently its additional r term is absorbed by the width penalty:

```math
\frac{r b^2}{N\ell}
\lt\frac{(17c_I)^2(n+1)(n+3)^2}{2^n\ell}=O(1).
```

This is uniform for every $`n\ge1,\ell\ge6`$: the numerator is
a fixed polynomial in n and the denominator contains $`2^n`$.
For example $`(n+3)^3/2^n\le32`$ for every positive integer n.
Thus $`r=O(N\ell/b^2)`$ on precisely this fallback branch.
It gives the asserted diagonal depth with the same optimal counts,
the same literal allocation, and absolute constants. No large-dimension
cutoff or precision-dependent finite exception is needed.

Both branches implement the same completed mask operators. Hence the
error bound in Section 2 is unchanged, including all rejected sectors,
dirty core disturbance, and arbitrary references. This proves the
standalone diagonal theorem. ∎

## 5. Half-error composition on the same two flags

For the complex frame set $`\ell=L+1`$. Apply the
[real-frame theorem R70](ALL_PRECISION_DEPTH.md#1-the-uniform-theorem)
at error $`2^{-\ell}`$ and the standalone diagonal theorem at the
same error. Both fit exactly the supplied condition

```math
b\ge17(\ell+n+7)=17(L+n+8).
```

Extend either circuit by identity on unused supplied wires, and use
the same two flags and the same b dirty inputs. Let
$`R_0=W_{\mathbb R}\otimes I_b`$ and $`D_0=D_\phi\otimes I_b`$.
Execute $`V_R`$ followed by $`V_D`$. The actual unitary circuits obey

```math
V_DV_RJ_2-J_2D_0R_0
=V_D(V_RJ_2-J_2R_0)+(V_DJ_2-J_2D_0)R_0.
```

Taking norms gives at most $`2\cdot2^{-\ell}=2^{-L}`$.
Earlier actual work leakage is propagated by $`V_D`$; only the ideal
comparison term uses initialized flags and unchanged dirty inputs.
There is no intermediate reset, reinitialization, or fresh work promise.
The identity remains valid after tensoring any reference.

The actual inverse retains the same full-isometry guarantee. With
$`V=V_DV_R`$ and $`W_0=D_0R_0`$,

```math
V^\dagger J_2-J_2W_0^\dagger
=V^\dagger(J_2W_0-VJ_2)W_0^\dagger.
```

Its native count, depth, and workspace equal those of the forward word.

Counts and depths add, while the workspace requirement is their
maximum. Since $`\ell=\Theta(L)`$, the T and Clifford counts have
the form claimed in Section 1. For the depth, uniformly for
$`n\ge1,L\ge6`$,

```math
F(n,L+1)=O(F(n,L)),\qquad
F(n,L)\ge L,\qquad F(n,L)\ge\log_2(n+2).
```

The first assertion follows by bounding each of the two entries of
the minimum after the shift $`L\mapsto L+1`$ by an absolute
multiple of its original value. The other two hold for both entries
of the minimum. Thus the extra diagonal depth is absorbed into
$`O(NL/b^2+F(n,L))`$ on the very circuit supplying the count bounds.

If $`L\gt N/n^2`$, then $`\ell\gt N/n^2`$, so R70 supplies
real depth $`O(N\ell/b^2+\ell)`$. Also

```math
\log_2(n+2)\le n+1=O(N/n^2+1)=O(L).
```

Here boundedness of $`n^3/2^n`$ gives the middle relation uniformly
in n. The diagonal adds only $`O(NL/b^2+L)`$. This proves the
stronger high-precision clause without shifting its stated condition
from L to $`L+1`$. ∎

## 6. Optimal count, matching depth, and scope

Set all magnitude angles to zero. The canonical complex magnitude
frame then contains every literal diagonal $`D_\phi`$. The
[existing diagonal-subfamily lower bound](../../docs/FAULT_TOLERANT_COMPILER.md#92-literal-diagonals-and-the-complex-magnitude-frame)
and its fixed-width counting consequence apply to the total physical
width $`n+b+2=\Theta(b)`$, giving

```math
T^\star=\Omega\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
D_T^\star=\Omega(NL/b^2).
```

Thus the construction attains the full worst-case T-count frontier
at every eligible allocation. This lower-bound argument uses the
diagonal family within the canonical angle domain; it does not
require identifying the entire unrestricted real-angle family with
a subfamily of that domain.

On every nonempty interval

```math
17(L+n+8)\le b\le\sqrt{\frac{NL}{F(n,L)}},
```

the same circuit has simultaneous matching orders

```math
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

The upper endpoint absorbs F into $`NL/b^2`$. Since $`F(n,L)\ge L`$
and $`F(n,L)\ge1`$, it also implies $`b\le\sqrt N\le N`$ and
$`b\le\sqrt{NL}`$, absorbing both L and $`\sqrt{NL}`$ into
$`NL/b`$. For integer n and L, the displayed interval is nonempty
exactly when

```math
289(L+n+8)^2F(n,L)\le NL.
```

As in R70, the high-precision upper bound supplies no matching
unrestricted lower bound for an additive L term. The result retains
the prescribed magnitude-frame completion and its leaf phases. It
does not compile arbitrary complex unitaries, reduce the clean budget
to one, assert practical crossover estimates, or supply a general
native emitter.

The native source, literal diagonal, and query identities already have
their own bounded evidence in their linked proof homes. This result
changes their schedule and composes their established full-isometry
contracts; it introduces no new native gate identity, source state,
or uncharged evaluator. No additional finite fixture is required to
justify the symbolic resource and width inequalities above.
