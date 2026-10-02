# A parallel geometric source on a conditional zero suffix

[Conditional suffix workspace](CONDITIONAL_SUFFIX_COMPILER.md) · [Operator source](OPERATOR_SOURCE_COMPILER.md) · [Exact source-depth scope](SOURCE_T_DEPTH.md) · [Current frontier](OPEN_PROBLEM.md)

A logical suffix that is zero on the active sector can supply a shallow
geometric source. The construction uses the same two external clean
flags as the original compiler. Its additional initialized work consists
of explicitly reserved logical suffix bits, and is available only after
computing the suffix predicate. On the inactive sector the completed
source operations are exactly identity on every input, including arbitrary
states of this work.

The geometric preparation and its zero-state reflections have
$`O(\log(m+1))`$ T-depth, $`O(m)`$ T and Clifford count, and at most
$`7m`$ conditionally clean suffix bits. This is a different source
interface from the arbitrary-dirty Majorana operator. It does not
resynthesize that operator, contradict its restricted depth lower bound,
or make its source core clean for free.

This improves the source-and-reflection part of the depth ledger.
Addressed queries and outer suffix predicates retain their separately
charged costs, so no improved complete-frame depth frontier follows here.

## 1. An exact parallel first-one preparation

Let $`m\ge2`$ and $`r=m-1`$. Use r coin wires
$`x_0,\ldots,x_{r-1}`$ and an endpoint wire y as the m-wire source
core. All are initially zero on the active sector. For
$`j=1,\ldots,r`$, write

```math
g_j(x)=\bigvee_{i=0}^{j-1}x_i.
```

An exact reversible prefix routine C writes these bits into r distinct
zero output wires and returns all its internal work to zero. Its action
on that initialized-work subspace is

```math
C:\ |x\rangle|0^r\rangle|0^{w}\rangle
\longmapsto |x\rangle|g_1(x),\ldots,g_r(x)\rangle|0^{w}\rangle.
```

Here is a self-contained linear-size implementation. In a balanced
binary tree with r leaves, compute the OR of each interval into one
fresh zero wire. There are $`r-1`$ internal interval wires. A down-sweep
computes the OR strictly preceding each subtree: copy the parent's
prefix to its left child, and OR that prefix with the left interval
value to obtain the right child's prefix. Reserve at most $`2r-1`$
prefix-node wires, including the zero root prefix. Copy the leaf-prefix
values for leaf indices $`1,\ldots,r-1`$ and the root interval value
to the r designated g outputs, and
reverse both sweeps. For r equal to one, copy $`x_0`$ directly to
$`g_1`$.

Computing $`u\lor v`$ into a zero target uses two CNOTs and one
Toffoli. At each sweep level, copy the controls of its OR gates into
private zero wires, perform the disjoint Toffolis, and uncopy those
controls. At most $`2r`$ reusable copy wires suffice. Consequently C
uses at most $`4(r-1)`$ Toffolis, $`O(r)`$ Clifford gates, and
at most $`16\lceil\log_2r\rceil`$ T-depth. Each Toffoli uses
the literal seven-T, four-layer phase schedule from
[the native batch proof](T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers),
specialized to one target pair and conjugated by target Hadamards.

The compute-copy-uncompute structure is essential. Internal interval
ORs need not remain invariant under the operations below. They are
already erased before any coin is changed; only the designated full
prefix outputs remain live.

Define V by this chronological sequence:

1. Apply H to every coin wire.
2. Apply C.
3. For $`j=1,\ldots,r-1`$, apply H to $`x_j`$ controlled by
   its distinct stored bit $`g_j`$. These controlled-H gates act on
   disjoint pairs and run together.
4. Toggle y by $`1\oplus g_r`$, using X and CNOT.
5. Apply the actual inverse of C.

### Why the work can be erased

On a computational basis string x, let its first one occur at i. Every
earlier coin has false prefix control and remains zero. Coin i itself
also has false prefix control and remains one. Only later coins are
Hadamard transformed. Therefore every prefix-OR value remains unchanged
on the complete resulting superposition. The all-zero string is also
unchanged by these controlled Hadamards. Equivalently, each controlled
Hadamard commutes with every prefix-OR Boolean projector.

The endpoint toggle does not alter a coin. Thus the stored g bits still
equal the prefix functions of the current coin register, and the last
actual inverse erases them. This argument holds for every coin input,
not just the uniform state produced in step 1. V preserves zero prefix
and internal work on the entire source-core input space.

Starting with the zero core, the initially uniform coins have probability
$`2^{-(j+1)}`$ of first one at j. The controlled Hadamards transform
the subsequent uniform tail to zero, with no phase. The endpoint wire
is toggled only for the all-zero coin string. Hence

```math
V|0^m\rangle|0^{\rm aux}\rangle
=\left(\sum_{j=0}^{m-1}a_j|e_j\rangle\right)|0^{\rm aux}\rangle,
\qquad
a_j=\begin{cases}
2^{-(j+1)/2},&j\lt m-1,\\
2^{-(m-1)/2},&j=m-1.
\end{cases}
```

The states $`|e_j\rangle`$ are the m one-hot core strings. These are
exactly the original squared geometric weights, including the duplicated
last weight.

### Native phases and an explicit reservation

Controlled-H is literal, with two T gates and no helper. Put
$`D=SH`$ and $`V_y=DTD^\dagger`$. Then
$`V_y ZV_y^\dagger=H`$, so the chronological word
$`V_y^\dagger,\mathrm{CZ},V_y`$ implements controlled-H.
The common scalar in $`V_y`$ cancels against its actual inverse.
All controlled-H pairs in step 3 use distinct wires.

The core, prefix outputs, interval nodes, prefix nodes, and reusable
copy pool occupy at most

```math
m+r+(r-1)+(2r-1)+2r\le7m
```

conditionally zero wires. The preparation has the conservative ledger

```math
T(V)\le58(m-2),\qquad G(V)=O(m),\qquad
D_T(V)\le32\lceil\log_2(m-1)\rceil+2.
```

For m equal to two, V is Clifford. The displayed depth bound is only an
upper bound. Neither this preparation nor its work is an arbitrary-dirty
implementation of the Majorana involution.

## 2. Programming the same coefficient table

For an m-bit row f, let
$`P_f=\prod_j Z_j^{f_j}`$ act on the one-hot core. Its initialized
compression is

```math
\langle0^m|V^\dagger P_fV|0^m\rangle
=\sum_j a_j^2(-1)^{f_j}=c_f.
```

The zero auxiliary work is implicit and returned exactly. Thus the
existing signed-dyadic encoding, including zero and both endpoints,
applies without increasing its m-bit program width.

For a preserved address x, an external active flag h, and a branch flag
b, the whole-word mask is the exact operation

```math
P_{h f(b,x)}=H^{\otimes m}Q_{h f(b,x)}H^{\otimes m},
```

where the existing arbitrary-input XOR query targets the m core wires.
There is no extra clean program register. The query address is
$`(h,b,x)`$, its output width is m, and its inactive h-zero rows are
the zero word. For a depth-d prefix this is at most $`4\cdot2^d`$
rows, within the original whole-word query ledger. All lookup selectors,
additional banks, and their borrowed work obey the existing exact-return
contract.

## 3. Two external flags suffice for the conditional source

For a Hopf layer at depth d, let its logical suffix have length
$`\ell=n-d-1`$. Compute $`h=[\mathrm{suffix}=0]`$ into the first
external zero flag, and reserve at most $`7m`$ suffix bits for the
preparation above when $`\ell\ge7m`$. The second external zero flag is b.
On h equal to one all reserved suffix bits are zero. On h equal to zero
they remain arbitrary, including entanglement with other inputs.

Set $`K=XZ=-iY`$ on the rotation target and define

```math
Q_h=H_b\,C_{h=b=1}(K)\,V^\dagger P_{h f(b,x)}V\,H_b.
```

The controls h and b are preserved between the outer Hadamards. On the
active initialized core/b subspace, its accepted block is
$`(cI+sK)/2`$ as before. The controlled K uses literal controlled-Z
followed chronologically by controlled-X, with both controls h,b; these
are constant-size exact native words.

On h equal to zero, the mask and controlled K are identity. The actual
$`V^\dagger V`$ cancels, even on arbitrary source and auxiliary inputs,
and the outer Hadamards cancel. Thus $`Q_h=I`$ on that entire sector.
No prefix helper is assumed clean there.

For amplification use

```math
R_h=I-2[\,h=1,\ b=0,\ {\rm core}=0^m\,],\qquad
A_h=Z_h Q_h R_h Q_h^\dagger R_h Q_h.
```

The factor $`Z_h`$ supplies the literal minus on the active sector.
On the inactive sector every completed factor is identity. On the active
sector this is the original normalization-two amplification, with a
zero-core block rather than an arbitrary-dirty scalar block. Prefix
helpers stay zero for every core input under V and its inverse, so the
reflection need not test them.

The reflection can use the existing two-dirty-helper exact multi-control
construction, with $`O(m)`$ T/Clifford count and
$`O(\log(m+1))`$ T-depth. Its two helpers are distinct from the
logical suffix, core, and external flags. The original dirty source-core
reservation is now unused by V and has at least two wires, so it can
supply these helpers between completed queries. The outer suffix
predicate uses the same two returned helpers at a different time. Its
$`O(\ell)`$ T count and $`O(\log(\ell+1))`$ T-depth are separate costs.

The two-helper multi-control construction is the existing imported
[Khattar–Gidney primitive](SOURCE_MAP.md), not a new reflection synthesis.
The whole A_h uses six V/V-adjoint occurrences and three whole-word
queries, plus constant target and flag gates and two reflections.

## 4. Complete-input error and exact inactive return

On the active sector, certified coefficient programming supplies the
same $`\zeta=\|(c,s)-(\cos\theta,\sin\theta)\|_2`$ as the
original source. The existing robust amplification proof gives the
full initialized-isometry error at most $`4\zeta`$ when
$`\zeta\le1/4`$, as ensured by the capped allocation below, including core
leakage. Prefix computation work is returned exactly on that sector;
the source core returns to zero within the same error bound.
The scalar block and reflection argument are restricted to the invariant
active subspace with zero prefix work. They do not assert a scalar
compression for arbitrary helper states on h equal to one.

The inactive sector is exactly identity on arbitrary suffix and dirty
inputs. Combining the two sectors gives an operator-norm estimate on
their coherent direct sum. Finally reverse the actual outer suffix
predicate. The ideal output has the original suffix and predicate value,
so this last unitary restores h in the ideal comparison and does not
increase the error. It is not assumed that the approximate output core
has returned exactly before this inverse. All statements extend to
arbitrary references.
When composing layers, unitary telescoping handles preceding leakage;
it does not assume that later calls receive freshly initialized work.

The original additive layer-error allocation therefore remains valid.
No radial filter is needed here: the native pi-over-three phase words
of that optional variant would reintroduce $`O(m)`$ precision depth.
There are no measurements, resets, extra external clean flags, or
supplied phase states.

## 5. What this improves in the depth ledger

Use the unfiltered additive capped assignment

```math
m_d=L+4+\min\{n-d,\lceil\log_2(8n)\rceil\}.
```

At layers with $`n-d-1\ge7m_d`$, replace the source and its local
reflections by this construction. Retain the original dirty source on
the remaining layers. Every source still has an m-bit program table;
the eligible replacement has $`O(m)`$ native count and
$`O(\log(m+1))`$ source/reflection T-depth. The unchanged original
dirty reservation suffices, and no conditionally clean suffix bit is
also counted as an arbitrary dirty query bank.

Write $`M=L+\log_2(n+2)`$. There are at most
$`O(\min(n,M))`$ ineligible late layers, since every assigned m is
$`O(M)`$. The source-and-local-reflection contribution is consequently

```math
T_{\rm source}=O(nM),\qquad G_{\rm source}=O(nM),\qquad
D_{T,{\rm source}}=
O\!\left(n\log(M+2)+\min(n,M)M\right).
```

When M is too large for a useful conditional suffix, this is only a
valid upper bound; the unchanged all-dirty source remains available.
At fixed L the displayed depth becomes

```math
O\!\left(n\log\log(n+2)+\log^2(n+2)\right)
=O\!\left(n\log\log(n+2)\right)
```

asymptotically. Its T count remains $`O(n\log(n+2))`$.

This ledger excludes addressed queries and the outer suffix predicate.
The current query schedule alone can still contribute
$`O(n\log(n+2))`$, and recomputing the outer predicate for each
eligible layer has that same depth upper bound. The complete-frame
count/depth frontier therefore remains unchanged. Sharing those costs
across a group would require a separate full-input circuit proof, with
its live program and workspace reservations charged.

## 6. Evidence boundary

The [bounded checks](../tests/test_conditional_geometric_source.py) audit
the first-one preparation, its native controlled-H phase, actual inverse,
programmed scalar block, and conditional amplification on small inputs.
They distinguish zero-helper source isometries from inactive arbitrary
helper return. Large-size prefix scheduling, imported reflection costs,
and the complete-frame ledger are analytic statements, not inferred
from those finite matrices.
