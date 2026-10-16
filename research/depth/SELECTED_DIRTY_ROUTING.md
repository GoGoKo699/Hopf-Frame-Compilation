# Selected routing and variable-size clean windows

[Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md) · [Native Fredkin schedule](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers) · [Interval-dilation compiler](INTERVAL_DILATION_DEPTH.md)

This supporting lemma supplies coherent physical routing, selected query
families, and the clean-window geometry used by the interval-dilation
construction. It does not assume a free selected permutation or charge a
variable-size table at the size of its largest padded physical address.
All dirty-return identities below hold on the complete input space,
including external reference correlations. T-depth permits intervening
Clifford circuits, whose elementary gate counts remain charged.

## 1. Exact routing selected by a short preserved label

Let R be a preserved k-bit register, H a preserved activity bit, and let
$`\pi_0,\ldots,\pi_{p-1}`$ be prescribed permutations of w target
wires, with $`k=\lceil\log_2\max\{p,2\}\rceil`$. Every target is
disjoint from R and H. Padded selector values specify the identity. There
is an exact Clifford+T word implementing

```math
|R,H\rangle|x\rangle|d\rangle
\longmapsto |R,H\rangle\pi_R^H|x\rangle|d\rangle
```

with no initialized routing work, literal phase one, and simultaneous
conservative bounds

```math
T=O\!\left(w\log(w+2)(k+3)^3\sqrt p\right),
```

```math
G=O\!\left(w\log(w+2)
  \left[p+(k+3)^3\sqrt p\right]\right),
\qquad
D_T=O\!\left(\log(w+2)\log(k+3)\right),
```

```math
w_{\mathrm{dirty,peak}}
=O\!\left(w(k+3)^3\sqrt p\right).
```

The dirty-width bound includes all padding, controls, indicators and
helper pools, but excludes the original target and public-address
registers. Finite zero-width cases are identities.

### A fixed network and its switch configurations

Pad to $`W=2^{\lceil\log_2\max\{w,2\}\rceil}`$ sites with arbitrary
dirty wires. Extend every prescribed permutation to fix the padding
sites. These sites may move during routing; only their final return is
required. In particular padding introduces no clean-work assumption.

A fixed Beneš network has $`2\log_2 W-1`$ stages, each a matching of
at most W/2 possible swaps. Recursively it has a first matching pairing
adjacent inputs, a last matching pairing adjacent outputs, and two
disjoint W/2-site subnetworks between them. The two members of each
outer pair feed the two subnetworks. The W=2 base case is one switch.

Here is a configuration proof. Make a bipartite multigraph with input
pairs on the left and output pairs on the right. Each input site gives
an edge to the pair containing its prescribed output. Every vertex has
degree two, so each component is an even cycle, including a two-edge
cycle of parallel edges. Alternately color each cycle's edges zero and
one. The first switches send the two colors into separate subnetworks;
the last switches send them to their prescribed output sites. Each
color now specifies a permutation of W/2 sites, which is routed
recursively. The stage recurrence is $`L(W)=L(W/2)+2`$.

Apply this finite classical construction separately to every selector
value. The switch settings $`\sigma_{s,j}(R)`$ are therefore ordinary
finite truth tables of the short unchanged selector. The active setting
is $`f_{s,j}(R,H)=H\sigma_{s,j}(R)`$. No coherent classical routing
algorithm is assumed. Fixed network wiring can be represented by fixed
physical pair lists; if emitted as wire permutations it uses only
Clifford SWAP/CNOT gates, with $`O(W\log W)`$ count included above.

### One native stage with arbitrary dirty controls

Allocate one arbitrary dirty control $`c_j`$ per switch. All controls
are disjoint from the routed sites, padding, R and H. The
[parallel-output query](SHALLOW_DIRTY_QUERIES.md#1-parallel-output-bits-with-a-read-only-address)
toggles $`c_j`$ by $`f_{s,j}(R,H)`$, using a separate arbitrary dirty
indicator and private helper pool for each output. Call this word Q.
The public address occurs only as Clifford CNOT controls. Every
non-Clifford support is private to its output query.

Let $`S_j`$ swap switch j's two sites and put
$`C=\prod_j S_j^{c_j}`$. Execute chronologically

```math
Q,\ C,\ Q^\dagger,\ C.
```

The controls and address are preserved by C. The pairs are disjoint and
each swap is a literal involution, so the net operation is

```math
\prod_j S_j^{c_j\oplus f_{s,j}}S_j^{c_j}
=\prod_j S_j^{f_{s,j}}.
```

Every dirty control and helper returns. This operator identity also
holds when the controls, routed data and external references are
entangled.

Each Fredkin consists of two CNOTs surrounding a Toffoli. The exact
[phase-polynomial schedule](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
implements a bank of Fredkins in at most four T layers, including
literal phases. Here the switch controls as well as target pairs are
disjoint, so those schedules run in parallel directly. A relative-phase
Toffoli is not substituted. A query's private helpers exclude every
routed wire, padding wire, control other than its declared output, and
every preserved public-address wire. Shared-address Clifford gates may
serialize between T layers; their counts are summed.

The stage query has $`2^{k+1}=O(p)`$ rows and O(W) output bits. The
conservative query bounds give

```math
T_s,w_s=O(W(k+3)^3\sqrt p),\qquad
G_s=O(Wp+W(k+3)^3\sqrt p),\qquad
D_s=O(\log(k+3)).
```

The two Fredkin banks add O(W) count and constant depth. Execute all
network stages sequentially, reusing each fully returned pool. This
proves the simultaneous bounds. The peak width is not multiplied by
the stage count. Use the actual reversed native word for the inverse.
When H is zero, every completed stage is literally the identity even
on occupied inputs. At the end, all padding is fixed by the extended
prescribed permutation and therefore returned exactly.

## 2. Parallel selected tables with different relevant addresses

Physical padding is not part of a truth-table address unless the
function actually depends on it. The following interface permits a
selected family to use different relevant physical subsets without
serializing the large queries.

For each bin i let $`A_i`$ be an ordered set of $`s_i`$ public input
wires, preserved throughout the query. Let
$`f_{i,j}(A_i,R,H)`$ be the jth output bit, extended by zero unless
$`H=1,R=i`$. The sets $`A_i`$ may overlap. Every final output target
and every query's private work must be disjoint from the **union** of
all public address wires, including R and H.

First query each $`f_{i,j}`$ into its own arbitrary dirty bit
$`c_{i,j}`$. All per-bin, per-output queries run in parallel T-depth:
their public wires occur only as read-only Clifford controls and their
non-Clifford supports are disjoint. Pad shorter native schedules with
identity layers. The total elementary Clifford count is the sum of
the separate counts. Denote the parallel word by Q.

For a shared final output word z, let C be the Clifford bank
$`z_j\mathrel{\oplus}=c_{i,j}`$ for every i,j. The chronological word

```math
Q,\ C,\ Q^\dagger,\ C
```

returns every private dirty bit and maps

```math
z_j\longmapsto z_j\oplus\bigoplus_i f_{i,j}(A_i,R,H).
```

The only active summand is the selected bin. C preserves every public
address because the final targets were excluded from their union.
This exclusion is essential: merely being disjoint from the selected
bin's address is insufficient for this parallel schedule.

For a Boolean phase, use a separate dirty phase target for each bin,
perform Q, apply Z to those targets, reverse Q, and repeat those Zs.
The resulting phase is
$`(-1)^{\bigoplus_i f_i(A_i,R,H)}`$. These identities hold on every
input, return dirty work exactly, and give literal identity for H=0.
No initialized flag holding $`[R=i]`$ is required.

If $`Q_i=2^{s_i+k+1}`$ bounds the ith table's rows and m is the
maximum output length, the conservative simultaneous bounds are

```math
T,w=O\!\left(m\sum_i(s_i+k+3)^3\sqrt{Q_i}\right),
```

```math
G=O\!\left(m\sum_i
   [Q_i+(s_i+k+3)^3\sqrt{Q_i}]\right),
\qquad
D_T=O\!\left(\log(\max_i s_i+k+3)\right).
```

Extra preserved indices simply contribute their actual address bits to
$`s_i`$. The helper widths are summed over simultaneous bins; they
cannot be reused within this parallel layer. Completed layers may reuse
the returned pools. This interface applies to the elementary Boolean
queries inside a source compiler, whose phase/program targets are
separate from its row and leaf inputs. It does not assert that whole
source unitaries sharing quantum targets can run in parallel.

## 3. Thin-bin geometry for an atom register of size q tau

Put $`\ell=\log_2(n+2)`$, choose an integer
$`1\leq\tau\leq C_\tau\ell`$, and set

```math
\kappa=\tau+2,\qquad m_0\geq4\kappa,
\qquad m_0=O(\ell).
```

The constant in $`m_0=c_0\ell+O(1)`$ is fixed by this requirement
before choosing a larger total protected suffix M. Write $`h=n-M`$.
This construction is used when $`h\geq1`$; a caller handles the
remaining finite-size range separately. Of the protected M suffix
wires, exactly $`m_0`$ form the free clean reserve. All persistent
metadata has a separate reservation in the other $`M-m_0`$ wires.

Construct endpoints backward, starting with b=h. At each step put

```math
v=h-b+m_0,\qquad
q=\min\!\left\{b,\left\lfloor\frac{v}{2\kappa}\right\rfloor\right\},
\qquad a=b-q,\qquad r=\tau q,
```

record the bin $`(a,b]`$, and continue with b=a until zero. Since
$`v\geq m_0\geq4\kappa`$, every nonempty step has $`q\geq1`$.
The capacity and address inequalities are

```math
v\geq2\kappa q=2(r+2q),
\qquad
a+q+r=b+r\leq h+m_0.
```

After unloading the selected endpoint b, the ideal early data tail
$`X_{b+1},\ldots,X_h`$ is zero. Together with the free protected
reserve it supplies exactly v exposed clean wires. Thus r atom wires,
two q-wire leaf/work fields and generous spare capacity fit without
borrowing any persistent metadata. If the local wrapper only needs
$`r+q`$ clean wires, the second q field can remain unused.

Every nonterminal step obeys
$`q\geq v/(4\kappa)`$, so v grows by at least a factor
$`1+1/(4\kappa)`$. Therefore the number p of bins satisfies

```math
p\leq1+
 \frac{\log(1+h/m_0)}{\log(1+1/(4\kappa))}
=O\!\left(\kappa\log(1+h/m_0)\right)=O(\ell^2).
```

The coarse selector needs only $`O(\log\ell)`$ bits. Moreover,
$`q_{\max}=O(n/\ell)`$ when $`\tau=\Theta(\ell)`$ and
$`r_{\max}=O(n)`$; more generally
$`q_{\max}+r_{\max}=O(n)`$. Fixed physical banks for all leaf and
atom fields therefore have O(n) wires.

The full arbitrary overlap table has preserved prefix a, input leaf
q, and atom index r. Its relevant address is bounded by
$`a+q+r\leq h+m_0`$. The interval-membership table used to erase a
leaf record needs prefix a and atom index r, excluding its output
leaf register. A source for the inverse interval embedding can use
prefix a, its fresh leaf register q, and atom register r. Separate
calls use their corresponding relevant subsets. No such table is
charged at $`h+q_{\max}+r_{\max}`$, and no arbitrary table depending
jointly on two leaf registers and the atom register is supplied by
this inequality.

## 4. A fixed selected layout for an entire local word

Allocate fixed **external arbitrary dirty** banks
$`X_{\mathrm{loc}}`$, Y and Z of lengths
$`q_{\max},q_{\max},r_{\max}`$, plus an optional q-maximal bank
if the wrapper needs the spare clean leaf/work field. They are disjoint
from the original early data X, the protected suffix, all metadata and
all query helper pools.

For bin i, prescribe a wire permutation with these completed effects:

1. Swap $`X_{a_i+1:b_i}`$ into
   $`X_{\mathrm{loc}}[1:q_i]`$, parking that bank's old arbitrary
   contents in the original local positions.
2. Choose fixed, distinct slots from
   $`X_{b_i+1:h}`$ and the free protected reserve, and swap their
   zeros into $`Y[1:q_i]`$ and $`Z[1:r_i]`$, as well as the
   optional extra active field if used. The old arbitrary bank contents
   are parked in those chosen slots.
3. Fix the preserved original prefix $`X_{1:a_i}`$, all metadata,
   and the inactive padding of each external bank.

The capacity proved above supplies all required slots. Their choice
and this permutation are finite classical data for bin i. Section 1
implements the selected permutation, guarded by H, including its
intermediate movement of any network padding. Only the **completed**
router is used as the layout boundary.

Keep this placement unchanged throughout the entire local word,
including every block-encoding call, actual inverse call,
amplification reflection, local fine-label XOR and its reverse
transport. The logical input is always
$`X_{\mathrm{loc}}[1:q_i]`$. Row/atom queries read their relevant
active physical bank prefixes and the preserved original prefix.
Section 2 compiles different bins' elementary tables in parallel.
All final query targets must still be disjoint from the union of that
call's public address wires; a separate program or phase target
satisfies this requirement.

The body must act as identity on parked dirty contents and inactive
padding for each selected bin. Its per-bin table functions ignore
these values, and a native query never modifies its public addresses.
Other constituents, such as variable-length balanced transforms or
selected rotations, must also supply this complete-space padding
interface. Merely approximating a map on fully zero-padded maximum-size
banks would not suffice.

After the complete local episode, apply the actual inverse router.
It restores every parked arbitrary bank value and places the local
logical output and active-work residuals back in their original slots.
The transported dirty bank values are untouched by the body, so their
return is exact; any initialized-work approximation error belongs to
the caller's operator-norm bound on the exposed slots. The statement
includes coherent superpositions of bins and dirty reference
correlations.

It is sufficient to route once around an amplified local compiler,
or once around the larger fine-erasure episode, rather than once per
amplification call. This requires one fixed selected placement for
both a word and its actual inverse, and the same placement for every
success reflection. For $`p=O(\ell^2)`$ and O(n) routed sites, the
two boundary routers have polynomial count and dirty width, no clean
work, and

```math
D_{T,\mathrm{route}}
=O\!\left(\ell\log\log(\ell+4)\right).
```

### Activity and occupied work

The router itself is exactly identity for H=0. A whole local
fine-erasure word may have the form
$`\widehat C^\dagger L_H\widehat C`$, where
$`L_0=I`$ and $`\widehat C`$ preserves H. It is then exactly identity
on the entire H-zero sector, even if every nominal work wire is
occupied: the two occurrences of $`\widehat C`$ are actual inverses.
No initialized-source promise is applied on that sector. This outer
cancellation can be used even if an intermediate local compiler has
no separate H-zero identity. Any different wrapper must establish its
own complete-space inactive action.

## 5. Active-only success reflections and metadata

An amplification success anchor tests precisely the initialized fields
of the selected local wrapper. It excludes the logical output
$`X_{\mathrm{loc}}[1:q_i]`$, the preserved prefix, the old fine tag,
parked dirty values and all inactive bank padding. It includes the
active Y and Z fields, any active optional work field, and every
source, damping, overlap and balanced-transform ancillary field that
the wrapper promises to return. In a wrapper ending with an X/Y swap,
the tested Y is the field containing the old input after that swap.

For each bin i form the structural conjunction

```math
g_i=[H=1]\,[R=i]\,
 [\text{every active ancillary wire for bin }i\text{ equals its anchor value}].
```

Use the
[read-only conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
with private arbitrary dirty targets and helper pools for every bin,
and the parallel Boolean-phase echo of Section 2. At most one g_i is
active. Each conjunction has O(n) literals, so these reflections have
$`O(\log(n+2))`$ T-depth and polynomial count and dirty width. This
is a structural conjunction, not an arbitrary truth table on the
union of all ancillary bits. No tested wire can be reused as a dirty
helper in its own or another simultaneous predicate. The entire union
of tested and preserved public wires is excluded from every private
helper pool.

Persistent metadata is reserved **outside** the free $`m_0`$ bank.
For the interval application it includes the old fine tag, endpoint
and coarse records, the balanced-transform tag, and the constant
number of fresh source/overlap/damping fields. In particular the
balanced tag has length $`O(\log(r_{\max}+2))=O(\ell)`$, not
merely $`O(\log\ell)`$. Source precision inverse-polynomial in ell
uses $`O(\log\ell)`$ bits per constant-size source bank. A fixed
$`O(\ell)`$ reservation covers these fields. The local source and
amplification proof remains responsible for specifying its complete
list; this lemma supplies its physical placement and exact predicate.

## 6. Table suppression and simultaneous count margin

Suppose there are polynomially many native elementary source/query
calls and output bits, and every extra source index has polynomially
many values. By Section 3, each per-bin large table then has at most

```math
\mathrm{poly}(n)\,2^{h+m_0}
=\mathrm{poly}(n)\,N2^{-(M-m_0)}
```

rows. Section 2 sums, rather than serializes, the p bins' private
query pools. Since $`p=O(\ell^2)`$, the total native envelope for
all large queries has some fixed exponent d satisfying

```math
T_{\mathrm{early}},w_{\mathrm{early}}
\leq O\!\left(n^d\sqrt N\,2^{-(M-m_0)/2}\right),
\qquad
G_{\mathrm{early}}
\leq O\!\left(n^dN2^{-(M-m_0)}
             +n^d\sqrt N\,2^{-(M-m_0)/2}\right).
```

All independent simultaneous pools are included in d. Returned pools
are reused between calls. Routing, structural anchors, fixed external
banks and small selected operations add only polynomial counts and
width. Their actual serial T-depth is charged separately.

First fix $`m_0=c_0\ell+O(1)`$ from the clean-work requirement and
fix the native polynomial envelope. Then choose

```math
M=m_0+\lceil A_{\mathrm{slack}}\ell\rceil,
\qquad A_{\mathrm{slack}}>2d+4,
```

increasing the constant further if needed for the separate metadata
reservation. This gives early T-count and peak dirty width
$`o(\sqrt N/n)`$ and Clifford count O(N). It leaves M logarithmic
and therefore compatible with the retained logarithmic-tail compiler.
Increasing M does not also increase $`c_0`$: that distinction is
required for the suppression estimate. The main compiler must still
check its actual number of source calls, its source errors and its
tail resources on the same physical circuit.

These lemmas establish routing and resource interfaces. The interval
normalization bound, rectangular-to-square local block, amplification
error and whole-frame recovery are proved in the
[interval-dilation construction](INTERVAL_DILATION_DEPTH.md).
