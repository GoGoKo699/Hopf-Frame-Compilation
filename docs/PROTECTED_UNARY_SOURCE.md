# One protected unary source across all early groups

[Unary source and literal shifts](UNARY_PHASE_GRADIENT.md) · [Uniform low-precision schedule](UNIFORM_PRECISION_DEPTH.md) · [Improved tail queries](NONUNIFORM_DIRTY_INDICATOR.md) · [Current frontier](OPEN_PROBLEM.md)

A fixed terminal logical block can hold the unary phase source throughout
all early groups. One initialized flag records that block's original zero
sector. The second flag records each group's activity without inspecting
the modified source block again. The source is prepared once and its
actual inverse is applied once before the late tail.

The preparation itself is unconditional. On the initial nonzero-block
sector, every intervening group is identity and the actual boundary
words cancel. Thus there is no unpriced controlled preparation or third
initialized flag. Source preparation error is charged twice over the
entire early segment, including final return of its initial-zero flag.

This removes the repeated early source-boundary depth allowance. Logical
stages, group activity predicates, and program loading still have their
separate linear upper allowances. The complete-frame frontier is unchanged.
The later [windowed-predicate refinement](WINDOWED_GROUP_PREDICATES.md)
also makes the activity-predicate allowance sublinear by reserving a small
cache in this bank. Logical stages and program queries remain separate.

## 1. Registers, flags, and the guarded interface

Use the established parameters $`n\ge1`$, $`N=2^n`$,
$`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
Split the logical register into A and a fixed terminal block B of
$`s_*`$ bits. B contains the q-bit unary source core, its binary index
and preparation work, the shift-output word, and the convolution
helper pool. All group programs, selected row words, selector trees,
enables, and their private work are outside B. The two external flags
H and h start at zero. Arbitrary dirty helpers and references are allowed.

First apply an exact predicate word $`C_*`$ computing
$`H=[B=0^{s_*}]`$, with two separately reserved returned dirty helpers.
Apply the existing native preparation U unconditionally on B. On B
initially zero it prepares a state in the unary source subspace exactly;
all non-core B work is zero after preparation. The approximate phase
words affect its amplitudes, not this exact encoding/cleanup property.
On a nonzero initial B, no property of U's output is assumed.

At a group with local targets $`t_0,\ldots,t_{g-1}`$, let Z be its
outer logical suffix with B removed. Compute into the initially zero h

```math
h=H[Z=0].
```

Use an exact multi-control predicate with the same two returned dirty
helpers. Its original controls exclude B. Load the h-enabled program,
run the existing unary stages and incremental selector schedule, return
all local work, apply the actual program inverse, and reverse this exact
h-predicate word. No group retests B as zero.

On H equal to one and h equal to one, every reserved group-work bit
inside Z was zero. The selected one-hot word and its inverse are exact
for every source-core input, as is the cyclic shift's temporary-word
return. Convolution and preparation helpers in B remain zero. The group
changes only its logical targets and the source core; it restores every
bit of Z. This remains true when the source is entangled with logical
inputs or contains the full preparation error. Hence h uncomputes to
zero exactly at every group boundary.

On a group with h equal to zero, the existing completed stage is identity
on arbitrary source and work inputs. The program's inactive row is zero,
and the actual selector and enable words cancel. Thus h returns exactly
also on inactive groups, without treating dirty copies as predicates.
In particular, on $`Hh=00`$ every complete group is identity on
arbitrary B and nonbank work. This last statement requires the initialized
second flag h; it is not a claim about arbitrary inputs with $`H=0,h=1`$.

After the last early group, apply the actual $`U^\dagger`$ on B,
then the actual $`C_*^\dagger`$. Only afterward begin the ordinary
late tail, whose targets may eventually enter B.

## 2. One global initialized-isometry error estimate

All early ideal layers are identity when the original terminal B is
nonzero, since every such layer requires its full outer suffix to be
zero. Therefore their prescribed ideal word has the form

```math
W_E=G_A\otimes\Pi_0+I_A\otimes(I_B-\Pi_0),
\qquad \Pi_0=|0^{s_*}\rangle\langle0^{s_*}|_B.
```

Here $`G_A`$ is the early grid-angle word restricted to B equal to
zero, including every A-input column; it need not act on every A bit.
Define embeddings with the external flag order H,h,

```math
E_0v=|0^{s_*}\rangle_B|10\rangle_{Hh}\otimes v,
\qquad
E_\phi v=|\phi\rangle_B|10\rangle_{Hh}\otimes v,
```

where v includes arbitrary A, dirty, and reference inputs, and
$`\phi=\lambda\Phi_q\otimes0_{\rm noncore}`$. The fixed scalar
lambda is the preparation's phase convention, independent of v. The
existing native approximation supplies

```math
\|UE_0-E_\phi\|\le\delta.
```

Let V contain all complete early group words, including their predicates,
program loads/unloads, and selector cleanup. Each ideal-source shift
has its literal eigenvalue. The exact interface of Section 1 gives

```math
VE_\phi=E_\phi G_A.
```

This identity does not assume a valid program was supplied externally:
each complete group includes its exact query and actual inverse. It
holds on all logical columns with all borrowed dirty inputs.

Since V and U are unitary,

```math
\begin{aligned}
\|U^\dagger VUE_0-E_0G_A\|
&\le\|U^\dagger V(UE_0-E_\phi)\|
   +\|(U^\dagger E_\phi-E_0)G_A\|\\
&\le2\delta.
\end{aligned}
```

On the orthogonal initial nonzero-B sector, computing $`C_*`$ leaves
$`Hh=00`$. V is identity there even after U acts on arbitrary B, so
$`U^\dagger VU=I`$ exactly on that sector. H is preserved until the
last predicate inverse. Combining both sectors and applying the same
actual $`C_*^\dagger`$ to both compared outputs proves

```math
F=C_*^\dagger U^\dagger VU C_*,\qquad
\|FJ_2-J_2(W_E\otimes I_{\rm dirty})\|\le2\delta.
```

Only the ideal output is known to restore B to zero before the final
predicate inverse. Actual source leakage and any residual H are included
in this norm. There is no assumed exact return of the approximate source,
no ideal-inverse substitution, and no error factor from the number of
groups. Relative phases between the two H sectors are preserved; lambda
cancels against the actual source inverse.

## 3. A disjoint static bank and a revised sufficient cutoff

Put $`q=2^\ell`$, $`\rho=\log_2 3`$, and $`R=q^\rho=3^\ell`$.
Choose one absolute constant $`A\ge1`$ large enough for both of the
following separate reservations inherited from the unary construction:

| Register family | Sufficient conditional-zero logical width |
|---|---:|
| Static protected bank: source, preparation/index work, shift buffer, and $`5R`$ convolution pool | $`A(q^\rho+q\ell+\ell+1)`$ |
| Dynamic group work: full program, selected row, private row-selection terms, prefix nodes/copies, and suffix enables | $`A(q2^g+g)`$ |

The second family is entirely outside B. Define

```math
s_*=\left\lceil A(q^\rho+q\ell+\ell+1)\right\rceil,
\qquad K=16s_*.
```

At remaining height $`k>K`$, take

```math
g=\left\lfloor\log_2\frac{k}{16Aq}\right\rfloor.
```

Use the same global phase modulus as the unary theorem, so $`q^2\ge n`$.
Then $`g\le\ell`$, and the cutoff ensures a positive group height
and $`g\ge\lfloor(\rho-1)\ell\rfloor`$. Moreover,

```math
s_*\le k/16,\qquad Aq2^g\le k/16,\qquad
Ag\le s_*\le k/16,\qquad g\le k/16.
```

Thus the static bank and all dynamic work together use at most
$`3k/16`$ logical bits, whereas the outer suffix has
$`k-g\ge15k/16`$ bits. Dynamic work fits outside the fixed B.
No conditional source bit is also counted as an arbitrary dirty query
helper. The two external flags and two dirty predicate helpers are
separately reserved.

The last early group can cross the cutoff by at most $`\ell`$ bits.
Its remaining height $`K'`$ therefore satisfies
$`K'\gt K-\ell\ge15s_*`$ and $`K'\le K`$ when the early
segment is nonempty. In particular no early target enters B, and the
source inverse and H erasure occur before any tail target can enter it.
If the early segment is empty, use the retained compiler without this
protected-bank construction.

## 4. Source error, depth, and the remaining complete-frame ledger

Choose q as before to be the least power of two with
$`q\ge4\pi\sqrt{2n}/\eta`$, but now set

```math
\delta=\eta/8,\qquad
\tau=1+\log_2((\ell+1)/\delta)
=O(L+\log\log(n+2)).
```

The phase words each have error at most $`\delta/\ell`$. Since the
source is prepared only once, Section 2 charges $`2\delta=\eta/4`$
for its entire use and final return. Early angle rounding still costs
at most eta divided by four under the complete-frame angular bound.
The old capped tail at precision parameter $`L'=L+2`$ still costs
at most eta divided by four. Exact table operations add no error.
Unitarity composes the tail with actual source/flag leakage without a
reset, giving total error at most eta.

The two global source words cost

```math
T=O(q+\ell\tau),\qquad G=O(q\ell+\ell\tau),\qquad
D_T=O(\ell+\tau).
```

The initial B-zero predicate and its actual inverse cost
$`O(s_*)`$ T and Clifford gates and $`O(\log(s_*+1))`$ T-depth,
with two returned dirty helpers. The entire source-boundary contribution
is consequently

```math
D_{T,\rm boundary}=O(\ell+\tau+\log(s_*+1))=O(L+\log(n+2)).
```

Each group, excluding its query and predicate pairs, now costs only

```math
T=O(q2^g+gR),\qquad G=O(q2^g+gqR),\qquad D_T=O(g).
```

Its h-predicate has $`k-g-s_*+1`$ controls, including H. Thus it
still costs $`O(k)`$ gates and $`O(\log(k+2))`$ T-depth.
The program has $`w=q(2^g-1)\le k/(16A)`$ bits. Its actual
prefetch/unload and private dirty-pool allocation retain the existing
separate query ledger. They do not borrow B or a live local program.

In the uniform range $`6\le L\le\log_2(n+2)/16`$,
$`q=O((n+2)^{9/16})`$ and

```math
K=O\!\left((n+2)^{9\rho/16}
              +(n+2)^{9/16}\log(n+2)\right),\qquad
K[L+\log(n+2)]=o(n).
```

These estimates have absolute constants. The early group count remains
$`O(n/\log(n+2))`$. All local gate counts remain polynomial in n
uniformly in this range, hence fit the same $`O(\sqrt{NL})`$ T
and $`O(NL)`$ Clifford budgets. The early parallel-query pools obey
the same exponentially suppressed bound, since every early group starts
above K; their widths fit the existing sufficiently large-width branch.
The improved nonuniform-indicator tail remains
$`O(NL/B_{\rm extra}^2+K[L+\log(n+2)])`$ in depth whenever
its sufficient extra-pool reservation fits.

At fixed accuracy and sufficient square-root dirty width, the updated
ledger is therefore:

| Contribution | Depth after this refinement |
|---|---:|
| Global source preparation/conversion/inverse and B-zero predicate pair | $`O(\log(n+2))`$ total |
| Early target-dependent shifts and incremental selectors | $`O(n)`$ |
| Early h-predicates and their inverses | $`O(n)`$ |
| Early program prefetch and unload | $`O(n)`$ |
| Improved late tail | $`o(n)`$ |

The source-boundary row is the improvement. The other three early
linear allowances are not lower bounds, but they still determine the
current complete-frame upper ledger. No sublinear complete-frame depth,
new matching interval, or high-precision endpoint bound is asserted.
The established full-frame theorems retain their literal width thresholds:
use this refinement only where its stated bank, cutoff, and query
reservations fit, and retain their proved fallback otherwise.

## 5. Proof and evidence boundary

The guarded group action, exact local work return on arbitrary source
inputs, and global two-boundary norm comparison are analytic statements.
The four checks in
[test_protected_unary_source.py](../tests/test_protected_unary_source.py)
audit the following bounded interfaces:

- The literal native $`q=4`$ preparation on its full six-wire bank,
  including its initial-zero column and the initial bank predicate.
- Unequal group heights $`(1,2)`$ and $`(2,1)`$ on every initial
  logical/bank input with both flags zero, and exact interior identity
  on $`Hh=00`$ with arbitrary prepared bank inputs.
- Analytic unitary source perturbations outside the unary code, the
  global $`2\delta`$ bound, exact group-h return, actual final H
  leakage, reference correlations, and inactive boundary cancellation.
- Negative cases that retest the modified physical bank or replace the
  actual preparation inverse by an ideal inverse.

The small source preparation is built from a native emitted word; its
actual inverse is evaluated as that matrix's adjoint. Predicates, signed shifts, and group
actions are reduced exact operators; its perturbations are analytic
unitary fixtures. It does not emit loaded-program queries or their
unloading, a native complete group, or a scalable complete-frame compiler.
The checks do not replace the symbolic reservation, uniform cutoff, or
complete-input error proof.
