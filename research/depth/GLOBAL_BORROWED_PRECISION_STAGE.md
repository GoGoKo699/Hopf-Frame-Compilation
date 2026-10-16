# A global early stage with borrowed precision work

[Coarse frame and sparse correction](COARSE_FLAG_SPARSE_CORRECTION.md) ·
[Small-workspace sources](SMALL_WORKSPACE_DESCENDANT_SOURCE.md) ·
[Causal compact flags](CAUSAL_COMPACT_COIN_DEPTH.md) ·
[Precision-bridge design](BORROWED_PRECISION_BRIDGE_PROMPT.md)

This chapter prices
the small-workspace early word with its actual unchanged global prefix.
Its clean reservation separates geometric metadata from requested
accuracy. It is not a complete-frame depth theorem by itself.

## 1. Uniform stage contract

Let the logical register consist of an unchanged prefix P of length
$`n-m`$, an a-bit early register X, and an r-bit bottom register Z,
where $`a=m-r`$. For each P prescribe a real a-bit Hopf frame
$`W_P`$, with all its signs and singular pairs. The target factor is

```math
E=\sum_P|P\rangle\langle P|\otimes
 \left(W_P\otimes|0^r\rangle\langle0^r|
       +I\otimes(I-|0^r\rangle\langle0^r|)\right).
```

For integer stage precision $`l\ge6`$, put

```math
K=m+l+2,\qquad t=\lceil\log_2(m+2)\rceil,\qquad
q_l=\lceil\log_2(l+2)\rceil,\qquad \chi=q_l+t^3.
\qquad(1)
```

There are absolute constants C, $`d_0`$, d and $`d_w`$ fixed by
the native constructions below. After fixing them and every clean-bank
coefficient, choose an absolute A sufficiently large, set
$`r=\lceil A\chi\rceil`$, and use this stage only when
$`r\lt m/2`$. It has one external initialized qubit H and at most b
arbitrary dirty qubits. Its complete initialized-isometry error is

```math
\|VJ_1-J_1(E\otimes I_b)\|_{\mathrm{op}}
\le2^{-l}/16.                                         \qquad(2)
```

Every table with a long local-data address, including P whenever its
data depend on P, has address length at most
$`n-r+C\chi`$. A constant number of such table positions are
sequential; at each position the number of simultaneous one-bit
instances is at most $`CK^{d_0}`$. Set

```math
\overline Q=2^{n-r+C\chi}.
```

A sufficient simultaneous allocation is

```math
b\ge 8W_{\mathrm{loc}}+CK^{d_w}(n+\chi+2),\qquad
W_{\mathrm{loc}}\le2^{C\chi}.                      \qquad(3)
```

Under (3), the same physical word satisfies

```math
\begin{aligned}
T(V)&\le CK^d\left(\sqrt{\overline Q}
                         +\overline Q/b+1\right)
                         +2^{C\chi},\\
G(V)&\le CK^d(\overline Q+1)+2^{C\chi},\\
D_T(V)&\le C\left[l+m+\log_2(n+\chi+2)
                         +K^d\overline Q/b^2\right].
\end{aligned}                                         \qquad(4)
```

Here $`W_{\mathrm{loc}}`$ includes every standalone source core,
external leaf/coin bank, local router, output, and returned helper
outside the private pools used for the global queries. An explicit
larger allocation within the displayed majorant may be chosen before
b. These are sufficient component widths, not a replacement for a
complete theorem's literal minimum. The initialization in (2) contains
only H; the bottom zeros used by the body are original logical inputs
on its active sector.

## 2. Geometry and metadata are independent of final precision

Use the coarse mathematical frame $`W_P^0`$ of the
[coarse/sparse proof, Section 2](COARSE_FLAG_SPARSE_CORRECTION.md#2-a-coarse-exact-tree-and-a-small-positive-residual-norm).
Each two-dimensional pair is rounded with error at most
$`\epsilon_0=[64(m+1)^2]^{-1}`$. Choose a single integer
$`\tau_0\le C t`$ for the dyadic probability precision, independent
of l and P. The clipping and finite sign procedure treats every node,
including nodes below a zero incoming amplitude. The exact matrix
$`D_P=W_P-W_P^0`$ has at most m+1 entries per row and satisfies

```math
\|\,|D_P|\|\le m\sqrt{m+1}\epsilon_0\lt1/16.
\qquad(5)
```

The approximation in (5) is mathematical. The physical implementation
of $`W_P^0`$ below is at final tolerance, not at $`\epsilon_0`$.

Choose a free donor reserve $`m_{\mathrm{free}}=\lceil c t\rceil`$
with $`m_{\mathrm{free}}\ge4(\tau_0+2)`$. It is disjoint from
all live source and metadata fields. The retained backward bin rule is

```math
q_i\le\frac{a-d_i+m_{\mathrm{free}}}{2(\tau_0+2)},
\qquad p\le C\tau_0\log_2(m+2)\le Ct^2.              \qquad(6)
```

All geometric quantities here are functions of m, r and the coarse
rounded tree; none requires accuracy l. The r dependence only limits
a and supplies sufficient free zeros. It does not enter a metadata
word as an arbitrary r-bit string.

For a bin of height $`q_i`$, enlarge its geometric exponent
$`s_i=q_i\tau_0`$ to

```math
s_i'=2^{\lceil\log_2(s_i+1)\rceil}-1,
\qquad s_i\le s_i'\le2s_i+1.                         \qquad(7)
```

Every original dyadic alias acceptance probability still has this
denominator. Its fixed coin now has $`2^{q_i}(s_i'+1)`$ atoms.
Each geometric weight class has at least $`2^{q_i}`$ atoms, and
the fixed coin uses $`q_i+O(t)`$ bits.

For each fixed P choose the shear catalogues separately, uniform only
over that P's original internal prefixes and fiber-cardinality groups.
There are at most $`2^m\mathrm{poly}(m+2)`$ such families.
At a fixed split and shear stage, each remaining family is fixed
conditional on the earlier choices. A random shear succeeds with
probability at least one half. Thus $`O(m+t)`$ choices at that
position suffice simultaneously by the retained union bound. There
are at most $`\mathrm{poly}(m+2)`$ stages and splits. The
selected group, split, stage and seed records therefore have O(t)
bits per block. They do not have $`O(\log n)`$ or O(l) bits.
The chosen matrices may depend on P, which remains part of every
physical matrix-table address.

The custom atom register $`Q_i`$ has $`q_i+O(t)`$ bits, as in
R89. Its weight classes have density at least $`(m+2)^{-C}`$.
The dense-XOR-matching proof, unioned over at most
$`2^m\mathrm{poly}(m+2)`$ original-prefix/weight families,
still needs only $`\mathrm{poly}(m+2)`$ shifts. Matching
indices consequently have O(t) bits per block. The shifts are chosen
separately for each P and their query addresses include P. Their
indices live outside Q, so the catalogue size does not feed back into
the density bound. Rebuild these matchings for the padded coins (7).

The full custom register has

```math
|Q|=a+O(pt)=a+O(t^3).                                \qquad(8)
```

The high-erasure and two matching-selector functions use causal
finite decoding, tabulated before the quantum word is emitted.
They are each one paid global query, not a physical loop through
prefixes. Catalogue selection is also finite preprocessing. No
common catalogue across all external P is requested.

## 3. Constant many ordinary source banks carry all precision indices

Use the [small-workspace descendant source](SMALL_WORKSPACE_DESCENDANT_SOURCE.md)
for the initial tagged source, selected endpoint sources, fine
rectangle source, joint rank source, and both sparse residual sources.
There are only constantly many chronological source positions. Set
every completed source tolerance to a sufficiently small fixed
multiple of $`2^{-l}`$, as fixed in Section 7. Its expansion length
is $`J=O(l+1)`$ and its initialized bank consists of

```math
\log_2J+1=O(q_l)
```

index and attenuation bits. A constant number of these banks is
reserved, including disjoint coarse and residual factor banks.

The calibration signal and its $`l+O(1)`$-bit precision core are
arbitrary dirty wires. They are not tested by any initialized source
anchor. The real borrowed-signal calibration has a full-unitary
estimate on all of these wires. Its scalar loaders are unconditional;
only the constant-arity central operations are enabled. Replacing
its mask queries by their exact R91 words preserves its complete
operator, actual inverse, inactive identity and literal phases.
There is no precision-sized initialized program, gradient or adder.

For the joint rank source, use exactly R89's
[metadata-selected target-set extension](CAUSAL_COMPACT_COIN_DEPTH.md#5-one-joint-source-on-metadata-selected-rank-windows).
At fixed metadata $`\mu`$, its target is the rank-window union
$`S(\mu)`$ in Q, and its preserved label is
$`(P,\mu,Q\setminus S(\mu))`$. Its source phases read one Q
register and one logarithmic expansion index. Variable Hadamard
masks and zero-anchor tests select only the rank windows. This does
not create a second a-bit source-address register. Source masks
preserve the metadata exactly; their targets exclude its positions.

Apply source completion on the entire physical initialization domain,
including arbitrary dirty work and external references. For actual
source U and exact initialized source isometry $`\Psi`$, it gives
an analysis-only $`U_*`$ with

```math
\|U-U_*\|\le6\delta,\qquad U_*J_{\mathrm{src}}=\Psi,
\qquad
U_*J_{\mathrm{src}}J_{\mathrm{src}}^\dagger U_*^\dagger
=\Psi\Psi^\dagger.                                  \qquad(9)
```

The last projector includes identity on every borrowed precision
wire. Source anchors exclude invalid padded tags and include all
their own initialized source fields. No exact actual dirty-core
return between source calls is assumed; earlier leakage propagates
through later full unitaries.

## 4. Exact reference words and their occupied-sector cancellation

The geometric part of a padded coin has power-of-two dimension
$`M_i=s_i'+1`$. The
[recursive Gray-path geometric source](../../docs/FAULT_TOLERANT_COMPILER.md#3-an-exact-geometric-source-of-logarithmic-width)
prepares it exactly, with T-count and Clifford count O($`M_i`$),
T-depth O($`M_i`$), and O($`\log_2M_i`$) initialized scratch.
Its scratch returns exactly for every coin input. The uniform-column
part uses Hadamards only. Equation (6), together with
$`m_{\mathrm{free}}\le r`$, gives

```math
s_i\le m/2,\qquad \max_i M_i\le m+2,
\qquad \sum_iO(\log_2M_i)=O(t^3).                    \qquad(10)
```

Run these reference words on disjoint coin and private scratch
registers, in parallel. Their depth is O(m), total count is at most
$`Cpt(m+2)`$ after a harmless enlargement, and no l-dependent
field occurs in their initialized banks. Include every private
scratch bit in the corresponding augmented nonzero-block test.
The root and all coarse marker bands then obey R89's full-band
identity, including padding and arbitrary band superpositions.

The reference transforms are unconditional. Their full word is the
actual reference inverse, the enabled coarse-label XOR, and the
actual forward reference transform. If the enable is zero, the
middle XOR is identity on the entire input space, so the surrounding
actual words cancel even with occupied coin or scratch registers.
This supplies exact inactive identity without parallel controlled
reference circuits sharing an unpriced initialized control bank.

The simultaneous conditional allocation can be fixed as follows.

| Register or purpose | Initialized size | Lifetime |
|---|---:|---|
| Code padding, fixed hash records and rank-window metadata | $`O(t^3)`$ | Compact embedding through its inverse |
| Causal matching-index register | $`O(t^3)`$ | Returned on the reached matched image |
| Exact reference scratch, one private bank per bin | $`O(t^3)`$ | Reference inverse through reference forward |
| Ordinary source indices and attenuation fields | $`O(q_l)`$ | Constant many separate source positions/banks |
| Fine/coarse tags, selected endpoints, router parameters | $`O(t)`$ per constant position | Computed and actually uncomputed |
| Residual depth tag, two failure coordinates and LCU branch | $`O(t)`$ | Combined block through recovery |
| Free fine-decoder donor zeros | $`O(t)`$ | Selected fine layout through its inverse |

All live rows of this table sum to $`C\chi`$ after including
disjoint coarse and residual factor banks. Source labels are separate
from source scratch. A protected logical bit appears in at most one
simultaneous allocation. The remaining Z bits are idle padding.

## 5. Global table addresses and physical scheduling

The retained coarse word is chronological: initial tagged source,
coarse-label computation, endpoint inverse, selected fine layout,
one rectangle-source/fine-label sandwich, inverse layout, endpoint
reload, endpoint uncomputation, compact compressor, one joint rank
source, causal matching, exact-reference coarse readout, and the
actual inverse compact embedding. Signs are separate exact phase
queries. Sparse correction and normalization-two recovery add only
constantly many repetitions of these positions.

At each Boolean-query position, use private read-only-address query
instances and the complete query/CNOT/unquery/CNOT output echo.
Whole sources sharing targets are not parallelized. Ragged selected
bins retain their own actual relevant addresses and share only
Clifford public controls. Their output union is physically disjoint
from the union of addresses, as in the retained selected-query lemma.

Write $`t_P=n-m`$. The following envelopes count the actual P
address even for short-catalogue functions.

| Table position | Relevant address bound | Separate target/action |
|---|---:|---|
| Initial tagged and endpoint source phases | $`t_P+a+O(t+q_l)`$ | Boolean phase / scalar mask |
| Fine rectangle source, bin i | $`t_P+\alpha_i+q_i+q_i\tau_0+O(t+q_l)`$ | Selected coin source |
| Fine membership, bin i | $`t_P+\alpha_i+q_i\tau_0+O(t)`$ | Fixed external leaf bank |
| Fine-label readout | $`q_i\tau_0+O(t)`$ | Fine tag |
| All compact metadata | $`t_P+a+O(t)`$ | Fixed metadata positions |
| Local shear catalogue | $`t_P+\lvert\mathrm{High}_i\rvert+O(t)`$ | Separate Low bank |
| Global high erasure | $`t_P+a+O(t^3)`$ | Separate High banks |
| Joint rank-source phases | $`t_P+\lvert Q\rvert+O(q_l)`$ | Separate source fields |
| Matching indices, first and last queries | $`t_P+\lvert Q\rvert+O(1)`$ | Separate matching register |
| Matching shifts | $`t_P+O(t)`$ | Disjoint $`Q_i`$ targets |
| Augmented coarse-block readout | $`a+O(t^3)`$ | Old coarse tag |
| Residual row and column source phases | $`t_P+a+O(t+q_l)`$ | Source fields / short tag |
| Ordinary source-anchor phase | Its source label/target union plus $`O(t+q_l)`$ | Boolean phase |
| Exact sign gauges | $`t_P+a+O(1)`$ | Boolean phase |

For the fine rows, $`\alpha_i=d_i-q_i`$ is the internal preceding
prefix length. The bin geometry proves

```math
\alpha_i+q_i+q_i\tau_0\le a+m_{\mathrm{free}}.
\qquad(11)
```

The actual endpoint-unloaded suffix supplies the corresponding
selected zeros. The contiguous interval routers park arbitrary old
external-bank values in separate original local/donor positions and
restore them through their actual inverses. No inactive padded bank
is included merely because a physical wire is present. Conversely,
no actual P or source-index field is omitted from its table price.

Equations (8) and (11) show that every row above has at most
$`n-r+C\chi`$ bits. Scalar calibration mask tables have
$`O(l+1)`$ output bits but their address is the unchanged source
label, not the precision core. The calibration core is a target/work
system, never a public address bank. Inner source reflections test
only the selected source-zero input and O($`q_l`$) source fields.
Transported anchors use (9) and never append an O(l) dirty-core test.

The all-work reflection in normalization-two recovery tests the
O($`\chi`$) work actually used by the combined LCU block. It does
not test idle Z padding or any arbitrary dirty pool. The original
$`[Z=0^r]`$ predicate is separately computed and inverted by a
returned-work conjunction on the local r bits. Its local
O($`\mathrm{poly}(m)`$) resources are charged; its r input
bits are not appended to any P-dependent source/table address.

All scalar source phases/calibrations have the same constant native
chronological grammar. One fixed polynomial bounds their local
output lengths, one-bit instances, router arrays and metadata fields.
For example p is O($`t^2`$), each output has O($`m+l+\chi`$)
bits, and active-stage $`\chi\le m`$; increasing one fixed
exponent covers their products and the retained mask schedules.
No factor polynomial in global n occurs in the instance count.

For shared H or P, non-Clifford gates are confined to disjoint
private query supports by the R88/R91 read-only interface. Their
Clifford control gates may serialize between T layers and are counted.
Local shear and shift tables consequently run in parallel across
bins. Contiguous routers have their proved O($`\log(m+2)`$)
depth with disjoint banks. The ordinary coefficient rotations use
their shared-H half-angle echoes on disjoint targets. The only
parallel whole source words are the unconditional, disjoint exact
reference words of Section 4. This is a physical schedule, not an
inference from commuting ideal projectors.

## 6. Native count, width and depth ledger

Fix all address constants and the finite polynomial instance exponent
$`d_0`$ in Section 5 before choosing A. At each global-query
position there are at most $`I=CK^{d_0}`$ one-bit instances.
Rows of each instance are at most $`\overline Q`$, after
increasing its fixed exponent coefficient. Give each instance a
private pool

```math
B_j=\left\lfloor\frac{b}{C_1K^{d_0}}\right\rfloor,
\qquad B_j\ge C_Q(s_j+2),
\qquad s_j\le n-r+C\chi.                             \qquad(12)
```

The second term in (3), with a sufficiently large fixed $`d_w`$,
ensures every individual minimum. The pools, private dirty outputs,
logical/data banks, calibration cores and routing helpers are
disjoint. After its completed exact word returns them, a pool may
be reused at the next sequential position. Allocating a fixed
fraction of b to all simultaneous $`B_j`$ leaves a fixed fraction
for the local reservations in (3).

Apply the read-only sharp query to every instance:

```math
T_j\le C(1+\sqrt{\overline Q}+\overline Q/B_j),\quad
G_j\le C(\overline Q+1),\quad
D_j\le C\left[1+\log_2(s_j+2)+\overline Q/B_j^2\right].
```

The constant number of chronological appearances, including every
actual inverse and output echo, yields the query terms of (4) for
one fixed $`d\ge3d_0`$. Multiplying Clifford gates by parallel
output counts is included in $`K^d`$; their depth is not charged
as extra T-depth.

Every table with a long local-data address is assigned to this same
query family, EVEN WHEN P IS EMPTY or the table happens not to depend
on P. This includes source phases and calibration masks, the fine
rectangle and membership functions, compact metadata and high
erasure, joint rank source, matching endpoints, the augmented block
readout, signed packing and both residual sources. Their actual
address bounds are precisely those in Section 5, so each is covered
by $`n-r+C\chi`$. Calling one of these tables "local" does not
place it in the standalone bill. Adding them to the family changes
only its fixed instance and sequential-position constants.

The only remaining standalone finite tables read short tags,
metadata, source indices, or routing parameters: at most
$`C\chi`$ bits, with a fixed polynomial number of outputs.
Give them separate sufficient square-root pools. R91 gives depth
O($`\log(\chi+2)`$), with count and simultaneous width at most
$`K^{C_2}2^{C_3\chi}`$. Since $`\log_2K\le C\chi`$,
one absolute enlargement of C absorbs the polynomial factor into
$`2^{C\chi}`$. Alternatively these short tables may also be
included in the global family; the displayed separate majorant
makes their simultaneous reservation explicit.

The other standalone contributions are explicit: each borrowed
calibration core and its native source loaders use O(l) wires,
O(l) depth and fixed-polynomial count; an ordinary coefficient bank
has O($`q_l[l+\log(q_l+2)]`$) count and O(l) depth; exact
reference words contribute O(m) depth and polynomial local count;
selected routers and saved-predicate conjunctions have polynomial
local count/width and logarithmic depth. Only constantly many
ordinary sources and layouts are sequential. These contributions
give $`W_{\mathrm{loc}}`$ and the standalone term in (4).
In particular, there is no residual $`2^{Cm}`$ standalone term.
A less precise $`2^{C(m+\chi)}`$ bound is valid but is not needed.

Crucially, this derivation uses neither an unspecified
$`2^{P(m+l)}`$ local compiler nor a table on l arbitrary program
bits. The l-sized quantum systems are dirty precision targets; all
precision-bearing source address bits number O($`q_l`$).

## 7. Complete error and return

The exact compact embedding and its causal matchings carry every
coarse band into the fixed-reference flag. Source completion (9)
justifies each transported anchor on the entire physical space.
The augmented exact reference test of Section 4 then erases the
coarse label with literal phase +1. The usual forward/middle/actual-
inverse inequality gives a constant multiple of the completed
ordinary-source errors. There is no sum over p approximate reference
sources because every reference here is exact. There is no sum over
P because all its preserved sectors are orthogonal. Hence the
complete coarse word has full initialized-isometry error
$`C_0\delta`$ from $`W_P^0\otimes I_D`$.

For the sparse correction retain the Schur weights and source pair
from the [coarse/sparse construction, Section 6](COARSE_FLAG_SPARSE_CORRECTION.md#6-native-sparse-residual-sources-and-normalization-two-recovery).
Equation (5) supplies the strict Schur norm gap. Its column source
keeps root tag zero live across all leaves, and its row source uses
the at most a+1 ancestor/root tags. Distinct failure coordinates are
orthogonal, and their exact source overlap is $`D_P`$, with its
literal signs. The two sources now use Section 3's borrowed-precision
implementation. They add only $`C_1\delta`$ compression error,
including their dirty-core leakage.

Disjoint conditional factor banks, one Hadamard LCU branch and its
actual inverse give the combined all-work-zero block

```math
\frac12(W_P^0+D_P)=\frac12W_P.                       \qquad(13)
```

Every SELECT enable includes H and its branch condition. Exact
enabled query/echo macros and constant-arity returned conjunctions
supply these conditions; a large physical word is not generically
controlled. On H zero SELECT is identity, and the PREP pair cancels
for arbitrary occupied branch and work inputs.

Let J initialize precisely the complete used conditional-work bank,
including both factor banks, branch, tags and reference scratch, and
let $`R_H=I-2HJJ^\dagger`$. It excludes P, X, idle padding and
all dirty systems. With rightmost factors acting first, use

```math
Z_H U R_H U^\dagger R_H U,\qquad Z_H=(-1)^H.           \qquad(14)
```

The retained normalization-two recovery has an absolute linear
block-error-to-full-isometry bound at success amplitude 1/2.
The sign $`Z_H`$ supplies the literal desired phase. Choose
$`\delta=2^{-l}/C_*`$ after the fixed number of source
appearances and this recovery constant are fixed. This changes the
source precision by an absolute additive constant only and yields
(2). All scratch leakage and approximate dirty return are included.
Every inverse in (14) is the actual native inverse, not an inverse
of a guessed ideal extension.

Finally invert the original saved predicate. On the ideal active
output all used work is zero and idle Z remained zero; H therefore
returns to its initialized value. On an original inactive input the
whole body is literal identity on its occupied work, including the
unconditional reference sandwiches, so the predicate inverse also
returns H exactly there. Earlier active leakage passes through this
last full unitary without enlargement. This proves (2) on all logical
columns, arbitrary dirty inputs and external references.

## 8. Intended use and scope

The separated reservation (1) contracts an input height O(l) first
to O($`\log^3(l+2)`$), then to O($`\log(l+2)`$), after an
absolute cutoff. If $`r\ge m/2`$, omitting the stage gives the
same corresponding next-height inequality with a factor two. Both
running stages have O(l) reference depth. These observations do not
by themselves absorb (3)–(4) into a global complete-frame budget;
that requires the separate first-stage, terminal-catalogue, error,
literal-width and complex-phase assembly.

The construction permits terminating nonuniform preprocessing and
uses no measurement, reset, QRAM, supplied gradient or free coherent
evaluator. It asserts no efficient preprocessing, practical crossover,
scalable native emitter or general lower bound.

## 9. Fast-reference variant

The [coherent capped-reference construction](FAST_CAPPED_GEOMETRIC_REFERENCE.md)
can replace Section 4's full exact reference sources. Keep every
mathematical coin, atom weight, alias assignment, hash and matching
unchanged. Only the native reference column is replaced by a shorter
exact geometric cap padded into the SAME full coin register. At
reference error $`\epsilon_i=2^{-l}/[C_*(p+1)]`$, it has
O($`l+\log_2(p+1)`$) native depth and the same O($`t^3`$)
total initialized scratch bound.

For completeness, if original and short cap dimensions are powers
of two $`R\lt M`$, their overlap and column distance are

```math
\langle g_M,g_R\rangle=1-(2-\sqrt2)2^{-R},\qquad
\|g_M-g_R\|^2=(4-2\sqrt2)2^{-R}.
```

Taking R to be the next power of two above
$`2[l+\lceil\log_2(p+1)\rceil+c]+2`$, capped at M,
provides the stated tolerance. The reference component constructs
an analysis-only full-space unitary comparison at this same error.
Its augmented zero-anchor projector is the original reference
projector, on the entire coin and ALL private scratch. Tensoring
these comparisons and using the actual inverse sandwich adds at
most $`2\sum_i\epsilon_i`$ to the coarse-readout error.
This is a full operator estimate, so it also propagates arbitrary
preceding source leakage and dirty/reference entanglement.

Fix c and $`C_*`$ after the finite recovery factors in Section 7.
The references still run unconditionally on disjoint physical
registers, and the middle enabled XOR still makes their sandwich
literal identity on every occupied inactive input. The original
high coin bits are retained in every augmented nonzero test; they
are not discarded because the actual prepared cap has zero high
bits. All matching and density arguments in Section 2 remain exact.

The replacement has polynomial count and no reference truth table.
Sections 5–6 therefore retain their SAME global query envelope,
standalone $`2^{C\chi}`$ bill, minimum pools and dirty-width
allocation. Its depth version of (4) is

```math
D_T(V)\le C\left[l+\log_2(n+\chi+2)
                  +K^d\overline Q/b^2\right].
\qquad(15)
```

Equation (2) and all count bounds are unchanged. Constants are fixed
before A as before. In particular (15) is available at m=n: all
large local-data tables remain within the globally suppressed
$`2^{n-r+C\chi}`$ envelope even though the external P is then
empty. Complete-frame range and terminal absorption still belong
to a separate assembly.
