# A shared binary geometric source for terminal layers

[Precision-uniform frame theorem](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision) · [Grouped terminal construction](GROUPED_TERMINAL_DEPTH.md) · [Shared-source compression](../../docs/FAULT_TOLERANT_COMPILER.md#8-shared-private-work-history-and-the-final-normalization)

One logarithmic-width binary geometric source can serve every middle
terminal layer. A retained failure counter supplies the intermediate
accepted-path projections, so native source preparation is paid only
once per complete stream. The residual of a single addressed qubit has
sixteen fixed matrix-unit/phase types; its coefficient queries need no
full endpoint register. A protected logarithmic terminal bank then
suffices, leaving only logarithmically many deepest layers to the
original arbitrary-dirty source.

The proof uses the retained
[exact capped binary source](../../docs/FAULT_TOLERANT_COMPILER.md#3-an-exact-geometric-source-of-logarithmic-width),
[contractive-shift kernel](../../docs/FAULT_TOLERANT_COMPILER.md#7-a-reusable-source-and-the-local-correction-kernel),
and [nonwrapping-history compression](../../docs/FAULT_TOLERANT_COMPILER.md#8-shared-private-work-history-and-the-final-normalization).
Literal activity controls use the
[bounded-control native lemma](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#exact-conditioning-without-a-clifford-count-penalty).
The coarse interpreter uses
[whole-word prefetch](GROUPED_TERMINAL_DEPTH.md#2-prefetch-the-exact-coarse-native-word),
and data access uses the established
[shallow dirty query](SHALLOW_DIRTY_QUERIES.md#1-parallel-output-bits-with-a-read-only-address).
All source words and their actual inverses are charged. No measurements,
resets, supplied precision states, or uncharged coherent evaluators occur.

## 1. Target and protected bank

Put $`N=2^n`$, $`\ell=\log_2(n+2)`$, $`p=L+\ell`$, and $`\eta_L=2^{-L}`$.
Keep the current early/terminal split

```math
6\le L\le cn,\qquad
M=\left\lceil A[L+\log_2(n+L+2)]\right\rceil=O(p).
```

Choose an absolute dimension cutoff $`n_0`$ large enough that the
existing protected-suffix split has $`M\le n/2`$ for every
$`n\ge n_0`$ and $`6\le L\le cn`$. It will be increased
below by further absolute thresholds. All constants are independent of
n, L, b, and the prescribed frame. Choose a sufficiently large $`A_0`$,
and set

```math
R_0=\min\{M,\lceil A_0\ell\rceil\},\qquad R=M-R_0.
```

Leave the deepest $`R_0`$ layers to the retained capped compiler. If
$`R=0`$, there is no stream; then $`M\le\lceil A_0\ell\rceil`$ and the
old bottom depth is already $`O(NL/b^2+p\ell)`$.

For $`R\ge1`$, save $`H=[Z=0^{R_0}]`$ into the first external clean
flag, where Z is the bottom logical bank. On H one, Z supplies
initialized work for the entire stream of R middle single layers.
No later local suffix predicate tests Z while it holds that work.
The original H computation is reversed only after full middle recovery.

**Terminal lemma.** For every $`n\ge n_0`$ in the displayed
precision range and every allocation

```math
a=2,\quad b\ge17(L+n+7),\quad b\ge\sqrt{N/n},
```

the final M layers have one circuit with complete-isometry error below
$`\eta_L/2`$ and simultaneous bounds

```math
T_T=O(\sqrt{NL}+NL/b),\quad G_T=O(NL),\quad
D_{T,T}=O(NL/b^2+p\ell).
```
 The unrestricted
parallel queries are used only on this large-width branch. The
retained low-width compiler remains a separate main-theorem fallback.

## 2. Sixteen residual types, with no endpoint register

At stage i let x be the unchanged logical prefix, t its target qubit,
and $`G_i`$ the all-zero predicate on logical bits below t and above Z.
This predicate excludes t and Z. Compute it into separate zero
scratch with direct H control before the stage and apply its actual
inverse afterward. On H one its value is the stated predicate; the
empty lower range gives $`G_i=1`$. On H zero the predicate word
is identity even when its target and helpers are occupied. Exact
returned-helper conjunctions charge this computation explicitly.

Choose actual, phase-calibrated native words $`C_{i,x}`$ for the desired
real rotations. Their addressed coarse action has direct H and $`G_i`$
controls. On H one, write $`W_i`$ for the desired full layer and $`C_i`$
for this exact coarse logical circuit. Its residual has the dictionary

```math
C_i^\dagger W_i-I
=\sum_{\lambda=(u,v,\omega)}
a_{i,\lambda}(x,G_i)\,\omega|u\rangle\langle v|_t,\qquad K=16.
```

Here $`u,v\in\{0,1\}`$ and $`\omega\in\{1,-1,i,-i\}`$. Use the
nonnegative positive/negative real and imaginary parts of the four
matrix entries. All coefficients are literally zero when $`G_i=0`$.
They are diagonal on preserved prefix/lower bits and commute with
the target matrix units. The only endpoint label is the pair of
one-bit values u,v; there is no n-bit endpoint or uniform prefix label.

For R at least one define

```math
q_\pi=\lceil\log_2(8R)\rceil,\quad
\pi=2^{-q_\pi},\quad \gamma=(1-\pi)^R,
```

```math
J=2^{\lceil\log_2(L+\lceil\ell\rceil+20)\rceil},\quad
M_s=2J,\quad
\varepsilon=\frac{\pi}{(1-\pi)KJ},\quad
\delta_c=\varepsilon/4.
```

Thus $`R\pi\le1/8`$ and $`7/8\le\gamma\le1`$. Full-operator native
synthesis at accuracy $`\delta_c`$ gives words of length

```math
w=O(\log(R+1)+\log J)=O(\log(p+2)).
```

Define the residual from
these actual words, retaining their literal phases. Each coefficient
satisfies $`0\le a_{i,\lambda}\le\delta_c=\varepsilon/4`$.

For each finite classical table row, take a nonnegative rational
lower enclosure within $`\varepsilon2^{1-J}`$ of a, then round down
to that grid. This gives Boolean digits with

```math
\widehat a_{i,\lambda}
=\varepsilon\sum_{k=0}^{J-1}f_{i,\lambda,k}(x,G_i)2^{-k},
\quad
0\le a_{i,\lambda}-\widehat a_{i,\lambda}
\lt2\varepsilon2^{1-J}.
```

Clipping certified rational enclosures at zero handles positive and
negative parts without an exact sign or exact-zero decision. Inactive
and structural-zero rows have literal zero digits. Coarse accuracy
only bounds these coefficients; the final target of the correction
is exactly $`C_i^\dagger W_i`$, so no separate coarse error is added.

## 3. Paid geometric source and the complete local kernel

The retained exact source preparation $`P_s`$ has initialized column

```math
g=\sum_{j=0}^{M_s-2}2^{-(j+1)/2}|j\rangle
 +2^{-(M_s-1)/2}|M_s-1\rangle.
```

It has $`O(\log M_s)`$ initialized width and $`O(M_s)=O(p)`$ gates
and conservative T-depth. Its actual inverse has the same cost.
Its preparation scratch is separate from retained failures/history
and returns exactly on every source input when initialized to zero.

For $`0\le d\lt M_s`$, let

```math
S_d=\sum_{j=d}^{M_s-1}|j-d\rangle\langle j|,\qquad
\|S_dg-2^{-d/2}g\|\le e_s
=\sqrt{4-2\sqrt2}\,2^{-M_s/2}\lt2\,2^{-J}.
```

The d=0 error is zero. All displacements are $`d=2k\le2J-2`$.
The exact dilation XORs $`[j\lt2k]`$ into an underflow flag, then
subtracts $`2k`$ modulo $`M_s`$. The flag is retained. Controlled
comparison and ripple subtraction on $`O(\log J)`$ bits cost
$`O(\log J)`$ count and depth, with separately returned exact scratch.
Their actions and actual inverses are defined on every source/flag input.

One shared failure register F contains rare-mode bits, their decoded
rare flag, the four-bit term label, the log J digit label, the loaded
digit, and distinct underflow, atom-rejection, and digit-zero flags.
It excludes the geometric source, history, H, global padding fields,
and every exact coarse/arithmetic/preparation scratch register.

On the initialized active F column, $`B_i`$ creates uniform rare,
term, and digit labels, decodes the all-zero rare word, and loads
$`f_{i,\lambda,k}(x,G_i)`$ into the digit. Uniform Hadamards and the
rare decoding are H-guarded. The query address is
$`(H,x,G_i,\lambda,k)`$, with zero H-zero rows; it excludes t and
the source. Consequently $`B_i`$ is identity on the full H-zero space.

The rare event is not treated as a bounded literal for free. Decode the
all-zero predicate on its $`q_\pi=O(\log(p+2))`$ bits into the
reserved rare flag using the retained exact returned-helper conjunction,
then erase it from the unchanged rare bits in the actual $`B_i`$
inverse. Its count and dirty width have the safe cubic bound and its
T-depth is $`O(\log(q_\pi+2))`$. Uniform H-controlled Hadamards
have their literal native implementations; their private target T gates
may run in parallel while shared-H CZ gates stay Clifford.

SELECT is identity in common mode. In rare mode with digit zero it
flips the separate digit-zero failure flag. In rare mode with digit
one it applies the downward source dilation and an independent
matrix-unit dilation. The latter XORs $`[t\ne v]`$ into its atom flag,
then flips t iff $`u\ne v`$. Its zero-flag block is $`|u\rangle\langle v|`$.
Apply the literal phase $`\omega`$ on that selected branch.
Both off-diagonal orientations use this same downward shift;
neither is constructed by inverting the source operation.

All consumed operations retain direct H and bounded rare/digit
controls. Literal controlled phases use the retained dirty-helper
echo lemma, with residual corrections, not relative-phase substitutions.
Adding these bounded controls to ripple arithmetic changes costs only
by an absolute factor.

SELECT preserves x, $`G_i`$, all labels, the rare flag, and the digit.
The actual $`B_i^\dagger`$ therefore unloads the digit even after t
changes and on failed branches. Dirty query work returns exactly.
The independent source and atom flags justify the product of their
accepted contractions without discarding intermediate rejection.

Let $`K_i=B_i^\dagger\mathrm{SELECT}_iB_i`$. On H one its F-zero
block is the exact contraction

```math
T_i=(1-\pi)I+\frac{\pi}{KJ}
\sum_{\lambda,k}S_{2k}\otimes
 f_{i,\lambda,k}(x,G_i)\,\omega|u\rangle\langle v|_t.
```

In particular $`\|T_i\|\le1`$. On $`G_i=0`$ all rare terms reject,
so its block is $`(1-\pi)I`$, not I. The digit-zero rejection is
not $`G_i`$-guarded: bypassing a locally inactive layer would make
global normalization depend on the input.

On H zero the complete $`K_i`$ is literal identity for arbitrary
occupied F/source/scratch/dirty inputs. On H-one failed inputs only
unitarity and exact scratch/query return are required.

## 4. Local and stream accepted-column estimates

Write $`E_g=|g\rangle\otimes I`$ and
$`D_i=(1-\pi)C_i^\dagger W_i`$. Source replacement and digit rounding give

```math
\|T_iE_g-E_gD_i\|
\le\pi e_s+\frac{2\pi}{J}2^{1-J}.
```

The source term is at most $`\pi e_s`$ because there are KJ
norm-at-most-one summands weighted by $`\pi/(KJ)`$. The second term
is $`(1-\pi)K`$ times the coefficient rounding bound. Thus no number
of prefixes or target entries multiplies this estimate. Both $`T_i`$
and $`D_i`$ are contractions, and $`\|D_i\|=1-\pi`$.

After $`K_i`$ execute the actual coarse $`C_i`$. Its whole native word
is prefetched into a separate zero program buffer, using address
$`(H,x)`$ and zero H-zero rows. Every symbol has direct H and $`G_i`$
controls. The target is excluded from both addresses. Actual unload
and the $`G_i`$ inverse return this buffer, predicate, and all helpers
on every logical input, regardless of retained F/source/history.
Its exact scratch is not F. On that exact-zero scratch subspace
the coarse word is $`C_i\otimes I_{F,\mathrm{source},\mathrm{history}}`$.

Between stages increment a separate history iff $`F\ne0`$, always
with direct original H control. For $`R\ge2`$ use
$`q_h=\lceil\log_2R\rceil`$ bits; for $`R=1`$ omit the counter.
There are R-1 possible increments, so no path from zero wraps to zero.
Every increment uses returned exact scratch, separately from F/history.

In detail, compute the nonzero predicate of F into an exact temporary
outside F, apply the H-controlled modular counter increment using that
temporary as enable, and reverse the predicate from the unchanged F.
The test has $`O(\log(p+2))`$ controls and uses the same
returned-helper conjunction bound. Its temporary and arithmetic work
return before the next stage. No many-controlled history operation is
assigned constant cost.


The final F-zero/history-zero block inserts F-zero projectors between
every stage. Indeed a path with any recorded nonzero F ends with
counter in $`\{1,\ldots,R-1\}`$; all other paths have intermediate
F zero. Coarse words commute with these F projectors, and later
kernels do not alter the history except for the specified increments.
The source is never projected between stages. The exact accepted block is

```math
V_{\rm acc}=(I\otimes C_R)T_R\cdots(I\otimes C_1)T_1.
```

Contractive telescoping on the source column yields

```math
\|V_{\rm acc}E_g-E_g\gamma W_{\rm mid}\|\le\xi,\qquad
W_{\rm mid}=W_R\cdots W_1,
```

```math
\xi\le R\pi\left(e_s+\frac2J2^{1-J}\right)
\le\frac38\,2^{-J}\le\eta_L/1024.
```

The bound uses J at least L+20 and J at least 6. Source disturbance
on earlier accepted steps is propagated contractively, not reset.
Apply the actual $`P_s^\dagger`$ after the whole stream. It uses
separate preparation scratch and cannot touch or clear retained
failures/history. Projecting source zero at this final boundary gives
$`\|V_0-\gamma W_{\rm mid}\|\le\xi`$.

## 5. Global padding and full occupied-space inactive identity

Choose

```math
w_0=\frac1{2\gamma}\in[1/2,4/7],\qquad
M_p=2^{\lceil\log_2(L+12)\rceil}.
```

A second independent capped binary source, followed by the exact
binary-digit sampler, prepares a padding mode with one probability
$`\widetilde w`$ satisfying

```math
0\le w_0-\widetilde w\le2^{1-M_p}
\le\eta_L/2048=:\epsilon_0.
```

Here $`w_0`$ is an exactly specified rational number, since $`\pi`$ is dyadic.
Its finite digits are computable exactly. Denote the actual source
and digit preparation by $`B_{\rm pad}`$. It has O(L) count/depth
and O(log L) initialized width, with its own exact scratch.
It does not touch F or history.

Define the complete half-block $`U_H`$ in chronological order:

1. Apply $`B_{\rm pad}`$ unconditionally.
2. Apply $`P_s`$ unconditionally.
3. Execute the H-guarded stream and history.
4. Apply the actual $`P_s^\dagger`$.
5. Flip a separate padding-rejection flag iff H is one and padding
   mode is zero.
6. Apply the actual $`B_{\rm pad}^\dagger`$.

The stream is not controlled by the padding mode. On H one, the
initialized-work block is exactly

```math
B=\widetilde w V_0,\qquad
\|2B-W_{\rm mid}\|\le2\xi+2\epsilon_0.
```

On H zero the entire interior, including history and padding, is
identity on arbitrary occupied inputs. The actual preparations cancel
in reverse nested order:
$`B_{\rm pad}^\dagger P_s^\dagger P_sB_{\rm pad}=I`$.
This does not assume clean source preparation scratch on H zero.
No inactive-enable property of the Gray preparation is used.
Both native preparation/inverse pairs are charged.

For simplicity allow separate exact scratch for main preparation,
padding preparation, arithmetic/history, and coarse program. Each
is O(log p). No inverse preparation borrows F or history, and
each preparation scratch remains exactly zero on every active
failed path while it is idle.

## 6. Full-isometry recovery and return

On H one let $`J_0`$ initialize all Z and the second external flag.
The saved H itself is excluded. Define $`R_H`$ to be identity on
H zero and $`I-2J_0J_0^\dagger`$ on H one. The zero test includes
source, padding, F, labels, history, and approximation work. It
never tests arbitrary dirty inputs. Testing all Z is valid: actual
source preparations have been inverted at the $`U_H`$ boundary.

The exact word

```math
\mathcal A_H=Z_H U_H R_H U_H^\dagger R_H U_H
```

has complete-isometry error

```math
8(\xi+\epsilon_0)\le3\eta_L/256\lt\eta_L/16
```

on the active initialized columns, by the retained normalization-two
lemma. Its inverse is the actual inverse of the whole finite word,
including failed branches, sources, coarse words, and history.
On H zero all factors are identity; the sign is $`Z_H`$, never an
unconditional minus sign.

Reverse the original H predicate computation last. The ideal middle
target preserves the original bottom-zero sector and returns its work.
The predicate inverse propagates the established error by unitarity;
it does not presume an actual reset or exact return after approximation.
This gives the full middle contract with arbitrary dirty references.

## 7. Simultaneous resources

Every initialized field fits $`C_Z\log_2(n+2)`$ protected bits:
two logarithmic-width sources; their separate exact scratch; rare,
term, and digit labels; constant atom/failure flags; log R history;
padding mode; arithmetic/history exact scratch; and the separate
O(log p) coarse buffer/predicate. The first external flag saves H,
and the second can hold padding rejection. Choose $`A_0\ge C_Z+2`$.

A stage has O(ell) depth for its constant number of shallow queries,
O(log J+log(R+1)) arithmetic/history depth, O(ell) local-predicate
depth, and O(log p) interpreted coarse symbols. Thus stage depth
is O(ell). Boundary source and padding preparation/inverses cost O(p)
per full stream. Global recovery has only three stream appearances,
one reversed; global reflections and H predicates add O(ell). Hence

```math
D_{T,\mathrm{middle}}=O(p+R\ell)=O(p\ell).
```

This uses no free reflection about g: the only recovery reflections
are zero-work reflections after actual source inversion, and history
implements the intermediate product.

At absolute depth d, a digit query has at most C J2^d rows and
constant output width; the coarse query has at most C 2^d rows
and O(log p) outputs. Their addresses exclude the target and have
at most O(n) bits. Here $`d\le n-R_0-1`$. There are a constant
number of completed queries per stage, including inverse calls and
the constant recovery multiplicity.

The conservative shallow-query lemma gives T/dirty width
$`O(m(q+2)^3\sqrt Q)`$ and Clifford count
$`O(mQ+m(q+2)^3\sqrt Q)`$. Uniformly J,R,p are O(n), so safe
common majorants for the full middle word are

```math
T_{\rm middle},w_{\rm middle}
\le C(n+2)^{10}\bigl(\sqrt N\,2^{-R_0/2}+1\bigr),
```

```math
G_{\rm middle}
\le C(n+2)^{10}\bigl(N2^{-R_0}+1\bigr).
```

They include all local source, arithmetic, native symbols, predicates,
history, reflections, and returned helpers.

There are at most $`O(n)`$ stages and long queries, source caps
$`O(n)`$, and at most $`O(n^2)`$ native symbol or predicate
occurrences. Charging the safe cubic count/dirty-work bound to every
local conjunction and the exact linear source and ripple words fits a
fixed polynomial of degree at most six; degree ten safely bounds the
combined terms. The short fixed padding tables fit that same polynomial
allowance. This is a native count bound, not an unspecified oracle cost.

Peak dirty work is reused
only after complete return. No initialized Z field is double-counted
as arbitrary dirty work.

Choose $`A_0\ge24`$ too. Increasing an absolute dimension threshold gives

```math
w_{\rm middle}=o(\sqrt{N/n}),\qquad
T_{\rm middle}=O(\sqrt N),\qquad G_{\rm middle}=O(N),
```

uniformly over $`6\le L\le cn`$. Thus it fits the supplied large-width
allocation and retains the optimal exponential count budgets.

## 8. Bottom and terminal composition

For the deepest $`R_0`$ layers retain exactly

```math
m_j=L+6+\min\{j,\lceil\log_2(8n)\rceil\},\qquad Q_j=4N2^{-j}.
```

The [retained capped layerwise allocation](SUBLINEAR_T_DEPTH.md#7-the-final-m-layers-retain-the-count-budgets)
has $`m_j\le L+6+j`$ and $`m_j=O(p)`$. Convergent
geometric sums, with absolute constants for $`L\ge6`$, give

```math
\sum_{j=1}^{R_0}Q_jm_j=O(NL),\qquad
\sum_{j=1}^{R_0}\sqrt{Q_jm_j}=O(\sqrt{NL}),\qquad
\sum_{j=1}^{R_0}m_j=O(R_0p).
```

Reserve its actual base $`B_0'=L+n+9`$ and additional returned
pool $`B=b-B_0'`$. Increasing the same absolute dimension cutoff,
$`L\le cn`$ and $`b\ge\sqrt{N/n}`$ imply
$`B\ge b/2`$ and $`B\ge16c_I(n+3)`$, where
$`c_I`$ is the absolute
[nonuniform-indicator constant](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
Thus its count-efficient rectangular queries fit this pool with no
precision-dependent hidden threshold. The weighted sums prove

```math
T_{\rm bottom}=O(\sqrt{NL}+NL/b),\quad
G_{\rm bottom}=O(NL),\quad
D_{T,\rm bottom}=O(NL/b^2+pR_0)=O(NL/b^2+p\ell).
```

Native source and returned predicate counts are polynomial in n in
this precision range and fit the same exponential count budgets above
the absolute cutoff. Their depth is already included in $`pR_0`$.

The bottom full-input error is below $`\eta_L/4`$. Middle and bottom
dirty pools are used sequentially, so their maximum reservation fits.
Full-isometry telescoping, with the actual bottom acting on middle
leakage, gives terminal error below
$`\eta_L/16+\eta_L/4\lt\eta_L/2`$ and the stated resource bounds.

The unchanged early segment can now be composed with this terminal
lemma on the large-width branch. Since $`p\ell=L\ell+\ell^2`$
and $`\ell^2=O(\sqrt n\,\ell)`$, their total depth is

```math
D_T=O(NL/b^2+\sqrt n\,\ell+L\ell).
```

The [main theorem](SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
handles the all-width and finite-dimension fallbacks separately.
This terminal lemma does not require its shallow query pools to fit
the literal minimum dirty allocation outside the stated large-width
branch.
