# Optimal count and all-precision depth with one clean qubit

[One-clean native compiler](../../docs/ONE_CLEAN_COMPILER.md) · [Real two-clean theorem](ALL_PRECISION_DEPTH.md) · [Phase-dressed two-clean theorem](COMPLEX_FRAME_DEPTH.md) · [Resource statement](../../docs/OPEN_PROBLEM.md)

**Local result R72.** One external initialized qubit suffices for the
simultaneous count and depth bounds of R70 and R71, at their original
dirty-width thresholds. The early compressed-column circuit already uses
one external flag. The terminal padding-rejection bit can occupy protected
logical work instead of a second external flag. A new width split supplies
the uniform fallback using the existing one-clean layerwise primitive.
No clean-to-dirty substitution or inference of depth from an optimal
total count is used.

The later [uniform subtree-depth theorem R79](UNIFORM_SUBTREE_DEPTH.md)
improves the all-precision additive depth at the same family-specific
widths. This chapter retains the baseline theorem and
the one-clean terminal, fallback, diagonal, and high-precision proofs
used by R79. The [R78 decoder](SUBTREE_PROJECTOR_DECODER.md) remains
the proof home for its nested-projector construction.

## 1. The one-clean theorem

Let $`n\ge1`$ and $`L\ge6`$ be integers, and define

```math
N=2^n,\qquad \eta_L=2^{-L},\qquad
F(n,L)=\min\{nL,(\sqrt n+L)\log_2(n+L+2)\}.
```

For either prescribed family in the following table, every integer
allocation $`b\ge17B_0`$ admits one coherent Clifford+T circuit V
using one external clean qubit and at most b arbitrary dirty qubits:

| Prescribed operator W | Base reservation $`B_0`$ |
|---|---|
| Complete real Hopf frame | $`L+n+7`$ |
| Complex magnitude frame $`D_\phi W_{\mathbb R}`$ | $`L+n+8`$ |

Writing $`J_1`$ for initialization of that sole external qubit,

```math
\|VJ_1-J_1(W\otimes I_b)\|\le2^{-L},
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+F(n,L)\right).
```

All constants are absolute and all bounds hold on the same circuit.
The norm includes every logical input, literal phases, work leakage and
return, arbitrary dirty inputs, and arbitrary reference correlations.
Unused supplied wires remain untouched. G counts elementary Clifford
gates, including those between T layers. T-depth permits arbitrary
intervening Clifford circuits.

For both families there is the stronger clause

```math
L\gt N/n^2\quad\Longrightarrow\quad
D_T=O\!\left(\frac{NL}{b^2}+L\right),
```

with the same simultaneous count, Clifford, width, and error guarantees.
The T-count is worst-case optimal at every eligible width. Section 8
states the matching depth intervals; an unrestricted depth optimum is
not inferred.

Angles and leaf phases admit terminating certified evaluation, including
singular and unequal tuples. Every leaf phase is prescribed, including
phases at zero-amplitude leaves. Classical preprocessing is terminating
without an efficiency claim. Every coherent query, source preparation,
actual inverse, guard, and recovery is charged. There is no measurement,
reset, QRAM, supplied precision state, or uncharged coherent evaluator.
The complex family is the specified phase-dressed magnitude completion,
not arbitrary complex-unitary synthesis.

## 2. Complete query replacements and returned predicates

We use the exact one-clean real primitive from
[Sections 1–5 of its proof](../../docs/ONE_CLEAN_COMPILER.md#1-contract-and-addressed-primitive).
At precision q it has a $`q+1`$-qubit arbitrary dirty core, one
initialized signal, a constant number of complete mask queries and
source words, and full-isometry error less than $`30\,2^{-q}`$.
Its source costs $`O(q)`$ T and Clifford gates and has a serial
$`O(q)`$ T-depth schedule. Long logical guards occur only on the
central X of each source conjugation; the loaders and masks remain
unconditional. The completed primitive is exactly identity whenever
its unchanged suffix predicate is false, even on occupied signal and
core inputs. Its five-call recovery uses four $`-Z`$ reflections;
physical Z gates preserve the literal scalar because there are four.

For a table with $`Q=2^r`$ rows and m output bits, retain the entire
mask operator, on every input:

```math
Q_f|x,y,z\rangle=|x,y\mathbin\oplus f(x),z\rangle.
```

The output y and all helpers z may be arbitrary and entangled. The
additional dirty query pool B is disjoint from the core, address,
signal, and reserved base. Z masks are complete XOR queries conjugated
by Hadamards on all m output wires. Keep the literal ordered Pauli
mask $`P_f=X^{u(f)}Z^{v(f)}`$ and its actual inverse. Replacing a
completed mask by the same operator tensored with identity on B
preserves every rejected block, inactive cancellation, and error bound.
No individual half-query or occupied table row is changed.

### Routed schedule at every literal allocation

The [routed rectangular allocation](COMPLEX_FRAME_DEPTH.md#routed-fallback-at-the-literal-threshold)
works for any m and r when $`B\ge4(r+2)`$. Its parameters are

```math
J=\mathrm{pow2floor}\min\{Q,\sqrt{Qm},B/4\},\qquad
H=\min\{J,Q/J\},\qquad K=Q/(HJ).
```

Here pow2floor rounds down to a power of two. The two arbitrary
indicator outputs occupy $`H+J\le B/2`$ wires; the high-address
selectors and one bilinear helper occupy
$`\log_2K+1\le r+1\le B/2`$ more. No initialized bank is used.
The exact four-call echo, with routed indicators and actual inverses,
has

```math
T_Q=O\!\left(m+\sqrt{Qm}+\frac{Qm}{B}\right),\qquad
G_Q=O(Qm),\qquad
D_{T,Q}=O\!\left(m+\frac{Qm}{B^2}+r\right).
```

Indeed its raw costs are $`O(KmH+H+J)`$ T gates and
$`O(Km+\log_2H+\log_2J)`$ T-depth. Power-of-two rounding gives

```math
KmH=Qm/J=O(m+\sqrt{Qm}+Qm/B),\qquad K=O(1+Q/B^2).
```

If $`Q=1`$, emit the fixed Clifford mask
without indicators. All Clifford coordinate changes and traversal
words remain included in $`G_Q`$.

### Shallow schedule when its absolute helper constant fits

Let $`c_I\ge1`$ be the absolute constant of the
[nonuniform rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
When $`B\ge16c_I(r+2)`$, the same complete query has the same
T and Clifford counts and

```math
D_{T,Q}=O\!\left(m+\frac{Qm}{B^2}+\log_2(r+2)\right).
```

The sufficient condition is checked where this schedule is selected;
its helper constant is not assumed to fit at the literal threshold.

For source guards, use the
[returned two-helper conjunction](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
A conjunction on p controls costs $`O(p+1)`$ exact Toffolis and
$`O(\log_2(p+2))`$ Toffoli depth, hence the same orders in native
T-count, Clifford count, and T-depth. It restores two arbitrary helpers.
Borrow both from B between completed queries, when every query helper
has returned. They are distinct from controls, the core target, and
the signal. The original reserved predicate helper may stay idle.
No predicate borrows a live query value, and no source loader becomes
a sequence of controlled T gates.

## 3. A one-clean fallback at every width

At layer $`d=n-j`$, $`1\le j\le n`$, use exactly the original
one-clean precision allocation

```math
q_j=L+j+5,\qquad m_j=q_j+1=L+j+6,\qquad Q_j=N2^{-j}.
```

The core, original $`n-j`$ selectors, and original helper reserve

```math
m_j+(n-j)+1=L+n+7=B_0.
```

Set $`B=b-B_0\ge16b/17\ge16B_0`$. This fits the routed query
for every layer and supplies both returned predicate helpers. The
query replacements preserve the original layer error and exact
inactive identity. Full-isometry telescoping gives

```math
\sum_{j=1}^n30\,2^{-q_j}
\lt\frac{15}{16}\,2^{-L}.
```

Earlier actual signal/core leakage is propagated by later unitaries;
only comparison terms use the ideal initialized signal and returned
core. The argument tensors with arbitrary references.

The weighted sums satisfy

```math
\sum_jQ_jm_j=O(NL),\qquad
\sum_j\sqrt{Q_jm_j}=O(\sqrt{NL}),\qquad
\sum_jm_j=O(nL+n^2).
```

The first two follow from convergence of the geometric sums in j,
uniformly for $`L\ge6`$. The serial sources, constant query
multiplicities, and central guards therefore give one circuit with

```math
T=O\!\left(\sqrt{NL}+nL+\frac{NL}{b}\right),\qquad G=O(NL),
```

```math
D_T=O\!\left(\frac{NL}{b^2}+n(L+n)\right).
```

Here the $`O(n^2)`$ count of source overhead and guards is absorbed
using $`n^2=O(\sqrt N)`$. The routed indicators contribute at most
$`O(n^2)`$ depth, already displayed. The sum of predicate depths
is $`O(n\log_2(n+2))=O(n^2)`$. No cap on the layer precisions is
needed. When $`L\le N/n^2`$, the remaining $`nL`$ count is at
most $`\sqrt{NL}`$, so this fallback already has optimal count.

### The final logical layers at smaller error

For only the deepest $`R_0\le n`$ layers use precision parameter
$`L+2`$, so

```math
q_j=L+j+7,\qquad m_j=L+j+8,\qquad 1\le j\le R_0.
```

Their base is $`B'_0=L+n+9`$ and their error is less than
$`15\eta_L/64\lt\eta_L/4`$. The original supplied threshold
still leaves

```math
B=b-B'_0\ge16b/17-2\ge b/2.
```

Whenever $`b/2\ge16c_I(n+2)`$, use the shallow query schedule.
The bottom circuit then has

```math
T_{\rm bottom}
=O\!\left(\sqrt{NL}+\frac{NL}{b}+R_0L+n^2\right),\qquad
G_{\rm bottom}=O(NL),
```

```math
D_{T,\rm bottom}
=O\!\left(\frac{NL}{b^2}
       +R_0[L+R_0+\log_2(n+2)]\right).
```

These are ungrouped one-clean native layers. In every refined branch
below, $`L\le N/n^2`$, so $`R_0L\le nL\le\sqrt{NL}`$;
their count is optimal there. The single external qubit is reused as
their signal after the preceding protected-work circuit, with a
full-isometry hybrid and no reset.

## 4. Relocate the terminal padding flag

The early circuit in
[R66, Section 10](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
already uses only one external initialized flag H. Its gradient,
program, addition register, tag and both attenuation fields are in
the protected logical suffix Z. Their simultaneous width is
$`O(L+\log_2(n+2))`$. Its
[shared-gradient reflection](SHARED_GRADIENT_RECOVERY.md)
excludes the retained gradient and includes the returned addition
register exactly as proved there. No change to this circuit or its
reflection is made.

The [binary terminal stream](SHARED_BINARY_TERMINAL.md#5-global-padding-and-full-occupied-space-inactive-identity)
and [whole-frame stream](WHOLE_FRAME_BINARY_DEPTH.md) save
$`H=[Z=0^{R_0}]`$ in the sole external initialized flag. Their
second external flag is used only for padding rejection. Choose one
distinct coordinate p inside Z for that rejection bit. It is outside
the main source, padding sampler, preparation scratch, private failures
F, retained history, and coarse/predicate program. It is touched only
by its padding flip and the complete zero-work reflection.

All other live initialized fields have width at most $`C_Z\lambda`$,
where $`\lambda=\log_2(n+L+2)`$ for the whole-frame stream and
$`\lambda=\log_2(n+2)`$ for the linear-precision terminal stream.
Enlarge the fixed bank constant to $`A_0\ge\max\{24,C_Z+3\}`$.
Then $`R_0=\lceil A_0\lambda\rceil`$ supplies those fields and
p simultaneously whenever a nonempty stream is used. This is an
explicit extra protected logical bit, not an initialized dirty input.
The source registers and their scratch remain logarithmic in their
caps; the early gradient has its separate, larger reservation above.

On H one, p is initially zero. Keep the exact chronological word

1. Prepare the padding sampler $`B_{\rm pad}`$ unconditionally.
2. Prepare the main source $`P_s`$ unconditionally.
3. Run the H-guarded stream and retained history.
4. Apply the actual $`P_s^\dagger`$.
5. Flip p iff H is one and the padding mode is zero.
6. Apply the actual $`B_{\rm pad}^\dagger`$.

The saved H is never recomputed from Z while Z holds work. On H
zero, the stream, history and p flip are identity on every occupied
input. The preparations therefore cancel in reverse nested order,
giving the complete identity

```math
U_H\big|_{H=0}
=B_{\rm pad}^\dagger P_s^\dagger P_sB_{\rm pad}=I.
```

This uses no initialized-work premise on the inactive sector. On H
one let $`J_0`$ initialize all Z, including p and unused protected
coordinates, but exclude H and arbitrary dirty inputs. The unchanged
accepted block is

```math
B=J_0^\dagger U_HJ_0=\widetilde w V_0,\qquad
\|2B-W_{\rm stream}\|\le2\xi+2\epsilon_0.
```

Use the physical reflection $`R_H`$ that is identity on H zero and
$`I-2J_0J_0^\dagger`$ on H one. The accepted projector includes p;
it does not demand that occupied p be reset on rejected paths. All
source inversions have occurred at this boundary. The actual recovery

```math
\mathcal A_H=Z_HU_HR_HU_H^\dagger R_HU_H
```

uses $`Z_H=I-2|1\rangle\langle1|_H`$ and retains its literal
active-sector sign and full-isometry error

```math
8(\xi+\epsilon_0)\le3\eta_L/256\lt\eta_L/16.
```

On H zero every factor is identity. Finally reverse the original
H predicate; its unitary action propagates the error, while ideal
outputs return Z and H. The arbitrary dirty/reference contract follows
on both sectors. The p flip is a constant-size controlled gate and
the reflection gains at most one control; their native bounds and
returned-helper allocations stay within the existing ledgers.

## 5. Refined circuits above a new width boundary

Put

```math
\ell=\log_2(n+2),\qquad \lambda=\log_2(n+L+2),\qquad
b_c=\sqrt{\frac{NL}{n(L+n)}}.
```

The uniform lower estimate

```math
b_c^2\Big/\frac{N}{n^2}=\frac{nL}{L+n}\ge\frac67
```

holds for every $`n\ge1,L\ge6`$. Thus
$`b_c\ge\sqrt{6/7}\,\sqrt N/n`$. The following cutoffs are
absolute, chosen once after the native constants, independently of L,
b, and the prescribed frame.

### Linear precision: early columns and a terminal stream

Use R66's constants $`A\ge\max\{16,C_s+1\}`$,
$`c=1/(8A)`$, and
$`M=\lceil A[L+\lambda]\rceil`$. For $`6\le L\le cn`$
and sufficiently large n, its protected early bank fits and
$`M\le n/2`$. The already proved early dirty-width envelope is

```math
w_E\le C_1n^{-2}\sqrt N+C_2n^6(2^{n/4}+1).
```

The terminal stream, with $`R_0=\min\{M,\lceil A_0\ell\rceil\}`$,
has the [binary terminal ledger](SHARED_BINARY_TERMINAL.md#7-simultaneous-resources)

```math
w_{\rm stream}\le C_3\sqrt N(n+2)^{-2}+C_4(n+2)^{10}.
```

Here its source cap is $`J=O(L+\ell)=O(n)`$, and
$`2^{-R_0/2}`$ suppresses its stated degree-ten table overhead
whenever the stream is nonempty, since $`A_0\ge24`$.
Increase the same cutoff so that each
displayed width is at most
$`\tfrac14\sqrt{6/7}\,\sqrt N/n`$, and so that

```math
\tfrac12\sqrt{6/7}\,\sqrt N/n\ge16c_I(n+2).
```

These inequalities eventually hold by exponential dominance. The
additional polynomial returned-helper pools may be included in the
fixed constants. Hence both protected circuits and the shallow
bottom queries fit every $`b\ge b_c`$. They run sequentially;
their work reservations are maximized, not added.

Keep R66's source and query counts. Section 4 moves only one padding
bit, so the early and stream budgets remain $`O(\sqrt{NL})`$ T
and $`O(NL)`$ Clifford gates. The early depth is
$`O(L+\sqrt n\,\ell)`$, and the terminal stream depth is
$`O((L+\ell)\ell)`$. Section 3's one-clean bottom at
$`R_0=O(\ell)`$ has the same terminal order, since its extra
$`R_0^2+R_0\ell`$ is $`O(\ell^2)`$. If $`R_0=M`$ there
is no middle stream, and this bottom estimate applies directly.
Thus one complete circuit has

```math
T=O(\sqrt{NL}+NL/b),\qquad G=O(NL),\qquad
D_T=O\!\left(NL/b^2+(\sqrt n+L)\ell\right).
```

Its full error is at most
$`9\eta_L/64+3\eta_L/256+15\eta_L/64\lt\eta_L`$.
For a missing middle stream omit its error. Each ideal intermediate
embedding initializes only the same external flag; actual earlier
leakage is propagated by later unitaries.

### Larger precision covered by the whole-frame stream

Use the relocated stream of Section 4 with the absolute enlarged
$`A_0`$ and $`\kappa=1/(4A_0)`$. For
$`6\le L\le2^{\kappa n}`$ and sufficiently large n, the
protected bank $`R_0=\lceil A_0\lambda\rceil\le n/2`$ exists.
The proof of R68 gives, with $`x=n+L+2`$ and
$`u=n+\lceil\log_2J\rceil+20=O(n)`$,

```math
u^6\sqrt J\,2^{-R_0/2}\le Cx^{-11/2},\qquad
u^6J\,2^{-R_0}\le Cx^{-17}.
```

Consequently its dirty-width envelope is
$`C_5\sqrt N(n+2)^{-11/2}+C_6n^6`$, which fits at
$`b\ge b_c`$ after the same type of absolute cutoff. Its native
count and Clifford bounds remain $`O(\sqrt{NL})`$ and $`O(NL)`$,
with upper depth $`O(L+n\lambda)`$.

Apply Section 3's one-clean bottom. In this precision range
$`\lambda=O(n)`$, $`R_0=O(\lambda)`$, and $`\ell\le\lambda`$,
so

```math
R_0[L+R_0+\ell]=O(L\lambda+\lambda^2)
=O((n+L)\lambda).
```

In the later selection this branch also has $`L\le N/n^2`$,
which absorbs the bottom's $`R_0L`$ count. It therefore supplies

```math
T=O(\sqrt{NL}+NL/b),\qquad G=O(NL),\qquad
D_T=O\!\left(NL/b^2+(n+L)\lambda\right).
```

The upper and bottom full errors total less than
$`3\eta_L/256+15\eta_L/64\lt\eta_L`$.
This argument does not presume an overlap with the high-precision
range. A larger $`A_0`$ can decrease kappa; Section 6 covers the
intervening precision range separately.

## 6. Select one real-frame circuit at every precision

Fix a single absolute dimension cutoff $`n_0`$ meeting all the
conditions above, R69's high-precision conditions, and
$`2^{\kappa n}\ge n`$ for $`n\ge n_0`$.

For $`n\lt n_0`$, use the original unbanked one-clean layerwise
compiler at the original base $`L+n+7`$. Its
$`T=O(N+nL)=O(L)`$ and $`G=O(NL)`$ bounds are uniform in
unbounded L after one enlargement depending only on the fixed
$`n_0`$. A serial native schedule has $`D_T\le T=O(L)`$.
Since $`F(n,L)\ge L`$, this proves both desired clauses in these
dimensions at every eligible width. This is an actual all-precision
circuit, not a finite-accuracy net.

Henceforth $`n\ge n_0`$. If $`L\gt N/n^2`$, use the actual
large-dimension circuit of
[R69](HIGH_PRECISION_DEPTH.md#1-high-precision-theorem).
It already uses the one-clean growing-source compiler and leaves
its stated second flag entirely untouched. Delete that spectator.
Every rescheduled complete mask is real on its full input space, so
R69's replacement after the conjugate-signal lift remains valid.
The borrowed signal stays in its dirty base, and no new T layer is
routed through it. This supplies the stronger
$`O(NL/b^2+L)`$ depth and optimal count. The preceding bounded-n
branch explicitly replaces R69's two-clean fallback.

All remaining precisions satisfy $`L\le N/n^2`$, so Section 3's
fallback has optimal count. If $`b\lt b_c`$, select that fallback:

```math
n(L+n)\lt NL/b^2.
```

Its entire additive depth is absorbed by the width term, at the
literal lower threshold $`b\ge17(L+n+7)`$.

For $`b\ge b_c`$, select as follows.

| Precision | Selected one-clean circuit | Additive depth |
|---|---|---|
| $`L\le cn`$ | Early columns and terminal stream, Section 5 | $`(\sqrt n+L)\ell`$ |
| $`cn\lt L\le2^{\kappa n}`$ | Whole-frame stream and bottom, Section 5 | $`(n+L)\lambda`$ |
| $`L\gt cn`$ and $`L\gt2^{\kappa n}`$ | Layerwise fallback, Section 3 | $`n(L+n)`$ |

Each displayed additive term is $`O(F(n,L))`$ on its selected range:

- In the first branch $`\ell\le\lambda`$. The ratio
  $`(\sqrt n+L)\ell/(nL)`$ is the sum of
  $`\ell/(L\sqrt n)`$ and $`\ell/n`$, both uniformly bounded.
  Thus the additive depth is bounded by an absolute multiple of both
  entries of F.
- In the second branch $`n+L\le(1+c^{-1})L`$ and
  $`\lambda=O(n)`$. Its depth is at most an absolute multiple
  of $`(\sqrt n+L)\lambda`$ and of $`nL`$.
- In the last branch $`L\gt2^{\kappa n}\ge n`$, so
  $`n(L+n)\le2nL`$, while
  $`nL\lt\kappa^{-1}L\log_2L\le\kappa^{-1}(\sqrt n+L)\lambda`$.

These cases exhaust every n, L, and eligible b. Every case selects
one complete circuit supplying all its counts, depth, allocation and
full-isometry error. This proves the real theorem. ∎

## 7. A one-clean diagonal and complex composition

For integer $`\ell\ge6`$, the
[one-clean literal diagonal](../../docs/ONE_CLEAN_COMPILER.md#literal-diagonal-phases-use-the-same-flag)
uses

```math
q=\ell+5,\qquad m=q+1=\ell+6,\qquad Q=N,\qquad r=n.
```

Its base is exactly $`m+n+1=\ell+n+7`$. Unlike the packed
two-clean diagonal of R71, its table address contains only x. The
fixed scalar word A has no addressed mask. With a denoting the sole
initialized signal, its actual phase word is

```math
Q_{\rm phase}=A^\dagger(S_a\mathcal S_gS_a^\dagger)A,
\qquad
\langle0|Q_{\rm phase}|0\rangle
=s+i(\tau-s/2).
```

Program $`s=\beta\cos\phi_x`$ and
$`\tau=s/2+\beta\sin\phi_x`$, where
$`\beta=\sin(\pi/10)`$. Use the original certified sign rule,
paired two-tail source, and five-call amplification. The four
$`R=-Z_a`$ reflection signs cancel literally. The full diagonal
error is less than $`30\,2^{-q}=(15/16)2^{-\ell}`$, including
signal leakage and dirty core return.

Each phase word has nine source appearances and four complete m-bit
XOR/Z queries. Five forward/inverse appearances give forty-five
sources and twenty queries. All source depth is $`O(q)`$; with
no suffix predicate their central flag-controlled X gates are
Clifford. Fixed masks and routing S gates are also Clifford.
These are constant multiplicities, different from R71's packed word.

Reserve $`B_0=\ell+n+7`$ and
$`B=b-B_0\ge16b/17`$. Use Section 2's shallow queries if
$`B\ge16c_I(n+2)`$, and its routed schedule otherwise. In
either case the diagonal has

```math
T_D=O\!\left(\sqrt{N\ell}+\ell+N\ell/b\right),\qquad
G_D=O(N\ell).
```

The shallow depth is $`O(N\ell/b^2+\ell+\log_2(n+2))`$.
On the routed branch failure of that width condition gives
$`b\lt17c_I(n+2)`$, and hence

```math
\frac{nb^2}{N\ell}
\lt\frac{(17c_I)^2n(n+2)^2}{2^n\ell}=O(1).
```

Its extra n depth is therefore absorbed by $`N\ell/b^2`$.
This proves the standalone one-clean diagonal statement

```math
b\ge17(\ell+n+7),\qquad
D_{T,D}=O\!\left(N\ell/b^2+\ell+\log_2(n+2)\right).
```

For the complex magnitude frame set $`\ell=L+1`$. The real
theorem at this precision and the standalone diagonal both fit
$`b\ge17(L+n+8)`$. Extend each by identity on unused dirty
wires and reuse the same sole external signal. With
$`R_0=W_{\mathbb R}\otimes I_b`$ and $`D_0=D_\phi\otimes I_b`$,

```math
V_DV_RJ_1-J_1D_0R_0
=V_D(V_RJ_1-J_1R_0)+(V_DJ_1-J_1D_0)R_0.
```

Its norm is at most $`2\cdot2^{-\ell}=2^{-L}`$. The first
term propagates all actual earlier leakage; no intermediate reset
or exact return is assumed. Counts and depths add; work is reused.
Uniformly for $`n\ge1,L\ge6`$,

```math
F(n,L+1)=O(F(n,L)),\qquad
F(n,L)\ge L,\qquad F(n,L)\ge\log_2(n+2).
```

These follow by bounding each entry of the minimum separately. They
absorb the extra diagonal depth and the half-error precision shift.
The real theorem is applied at its shifted precision in every case,
including any crossing of its internal branch boundaries.

If $`L\gt N/n^2`$, then $`\ell\gt N/n^2`$ and the real factor has
depth $`O(N\ell/b^2+\ell)`$. Boundedness of $`n^3/2^n`$ gives

```math
\log_2(n+2)\le n+1=O(N/n^2+1)=O(L).
```

The diagonal therefore preserves the high-precision
clause with its stated condition on L. This proves the complex
theorem. Actual inverses preserve the same error because, for
$`W_0=W\otimes I_b`$,

```math
V^\dagger J_1-J_1W_0^\dagger
=V^\dagger(J_1W_0-VJ_1)W_0^\dagger.
```

Their native counts, depth, and workspace equal those of the forward
word. ∎

## 8. Matching intervals and evidence

The physical width is $`n+b+1=\Theta(b)`$. The
[inherited complete-frame lower bounds](../../docs/FAULT_TOLERANT_COMPILER.md#10-matching-lower-bounds-and-their-lineage)
and [literal-diagonal subfamily bound](../../docs/FAULT_TOLERANT_COMPILER.md#92-literal-diagonals-and-the-complex-magnitude-frame)
therefore imply

```math
T^\star=\Omega(\sqrt{NL}+L+NL/b),\qquad
D_T^\star=\Omega(NL/b^2).
```

The diagonal family is contained in the complex magnitude family by
setting its magnitude angles to zero. Thus count is optimal for
each family at every eligible width. For either base from Section 1,
every nonempty interval

```math
17B_0\le b\le\sqrt{NL/F(n,L)}
```

has simultaneous worst-case matching

```math
T^\star=\Theta(NL/b),\qquad D_T^\star=\Theta(NL/b^2).
```

Its integer-width nonemptiness is exactly $`289B_0^2F(n,L)\le NL`$.
The upper endpoint absorbs F into the width penalty; since
$`F\ge L\ge6`$, it also absorbs L and $`\sqrt{NL}`$ into
$`NL/b`$. The high-precision depth is $`O(L)`$ at every eligible
width, as in R70; no unrestricted matching lower bound for L follows.
Optimal depth remains an open problem.

The proof changes register allocation, completed query schedules and
the width splice. Its native one-clean scalar/amplification identities,
source words, exact dirty queries, predicates and stream history retain
their existing proof homes and bounded checks. No additional finite
fixture proves or is needed for the asymptotic allocation inequalities.
In particular, existing two-clean complete-frame fixtures are not
claimed to emit this new one-clean schedule. R70 and R71 remain valid
at their stated two-clean allocations. Result C's narrower base and
banked width thresholds, the sharp two-clean complex endpoint, and
the borrowed-signal real count theorem retain their own scopes. No
general native emitter, efficient preprocessing guarantee, or practical
crossover is supplied.
