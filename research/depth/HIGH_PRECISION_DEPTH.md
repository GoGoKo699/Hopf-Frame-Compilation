# A same-circuit depth schedule at high precision

[Optimal banked count frontier](../endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier) · [Conjugate-signal interfaces](../endpoint/CONJUGATE_SIGNAL_RECOVERY.md#banked-real-interfaces) · [Rectangular dirty query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence)

At sufficiently high precision, the existing growing-source compiler can
retain its optimal count while its complete precision masks use a separate
rectangular-query schedule. The masks are real full-space operators, so
they can be resynthesized after the conjugate-signal lift without sending
their native T layers through one shared signal. Every remaining word has
a separately bounded native count that is only linear in the precision
in this regime. This proves a depth bound on the same actual circuit.
Its stable result identifier is **R69**; the
[all-precision theorem](ALL_PRECISION_DEPTH.md) combines it with the
lower-precision constructions.

## 1. High-precision theorem

There is an absolute constant C such that, for every integer $`n\ge1`$,
$`N=2^n`$, integer $`L\ge6`$ satisfying $`L\gt N/n^2`$, and allocation

```math
a=2,\qquad b\ge17(L+n+7),
```

every prescribed complete real Hopf frame W has one native Clifford+T
circuit V with at most b arbitrary dirty qubits and

```math
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-L},
```

```math
T\le C\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G\le CNL,\qquad
D_T\le C\left(L+\frac{NL}{b^2}\right).
```

Here $`J_2`$ initializes only the two external clean flags. The norm
includes every logical input, work leakage, arbitrary dirty inputs and
external references, with literal phase. Above one absolute dimension
cutoff the construction uses the existing one-clean compiler and leaves
the second flag untouched. Bounded dimensions use the already proved
two-clean simultaneous count/depth compiler. No measurement, reset,
initialized query bank, or new approximation is introduced.

The inherited construction is the complete real-frame
[growing-source compiler](../endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
with its
[optimal banked count ledger](../endpoint/GROWING_SOURCE_COMPILER.md#12-optimal-banked-count-frontier),
[disjoint block/source interfaces](../endpoint/DISJOINT_BLOCK_COMPILER.md#banked-high-precision-interfaces),
and conjugate-signal lift. The replacement is the
[nonuniform rectangular query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).

## 2. Exact replacement at a completed mask boundary

Let F be the original classical m-bit table on $`Q=2^k`$ addresses.
Its completed XOR mask is the full-space real permutation

```math
Q_F|x,y,z\rangle=|x,y\mathbin\oplus F(x),z\rangle,
```

where z contains every query helper and bank. Both y and z may be
arbitrary and correlated with all other registers. The new rectangular
query implements exactly this operator, including its literal phase.

Retain every row of the original table. In particular, do not zero rows
on an inactive or occupied logical sector unless the original table
already has those zero rows. In the growing-source construction masks
are generally unconditional: only central source gates receive the long
arena predicate, and the masks cancel with their actual inverses on the
inactive sector. Changing those rows would not be licensed by the
original source proof.

The borrowed signal r is outside the query address, output, and helper
pool. The conjugate-sector map is

```math
\Pi_r(U)=|0\rangle\langle0|_r\otimes U
 +|1\rangle\langle1|_r\otimes\overline U.
```

Since the completed mask is real,

```math
\Pi_r(Q_F\otimes I_{\rm helpers})
 =I_r\otimes Q_F\otimes I_{\rm helpers}.
```

Thus, after lifting the high-level word, replace a completed lifted mask
by the new real query that leaves r untouched. Do not lift the individual
T gates of the new query: they would share r, and their native depth
schedule would require a different proof. Equality of the completed
operators on the same physical space already licenses the direct
synthesis. It remains exact when preceding finite phase slots have mixed
the r sectors.

An addressed Z mask is $`H_{\rm out}Q_FH_{\rm out}`$, with one H on
each output wire. It is real as well, and its replacement adds only
those Clifford layers. The paired source's literal mask
$`P_f=X^{x(f)}Z^{z(f)}`$ is real even when non-Hermitian. Keep that
order and its actual inverse $`Z^{z(f)}X^{x(f)}`$; no canonical-Pauli
phase is discarded. Replace the complete full-output mask, not an
internal half-query of the old self-borrowed implementation.

On a common physical base and additional B-bit bank, the old and new
completed masks equal the same base operator tensored with $`I_B`$.
Unused helpers in either implementation are spectators. This full-unitary
equality preserves the entire rejected operator of each cheap word, its
mixed gauge, and its actual column replay on the same physical core.
Fix the new actual row word before computing the gauge, then replay that
exact word and core in the column, as in the inherited construction.

The actual inverse of a replaced word is its literal reversed native
word. All inherited address-mutation boundaries remain unchanged.
Short metadata is uncomputed before a carrier endpoint changes.
Conditional metadata echoes use the original extra short control bit
as an address bit, adding only a factor two to Q. Their long arena guard
and returned helper are unchanged. No long precision loader is promoted
to a gate-by-gate controlled word.

## 3. A disjoint pool with the literal dirty allocation

Reserve the inherited base

```math
B_0=L+n+7,
```

including the arbitrary signal r, the currently live precision core,
and the other simultaneous guard/metadata reservations. The inherited
high-precision proof ensures, after one absolute cutoff:

- Fine cores have $`L+O(1)`$ wires and leave room for r and fixed helpers.
- Cheap, coefficient, and refinement cores are smaller than $`L/2`$.
- Their remaining live helper reservations are polynomial in n and fit
  the rest of the base.
- Coarse/routing words use $`O(\sqrt N+\mathrm{poly}(n))`$ dirty wires,
  also inside the base, in their separate stages.
- Original logical arena coordinates are never counted as dirty storage.

Use only the additional pool

```math
B=b-B_0\ge\frac{16}{17}b\ge16B_0
```

for the new query's indicators, high-address selectors, and controlled
bilinear helper. This pool is disjoint from the output core, all public
address and guard bits, r, and every other live field. It is reused only
after a completed query returns it exactly. No part of a source core
serves as selector for these new masks.

Every replaced query has $`k\le C_a(n+1)`$ address bits for an absolute
$`C_a`$. Coefficient tables have $`O((N/2^d)P_d)`$ rows with
$`P_d=O((d+1)^2)`$, fine/source banks have at most $`4N`$ rows, and
bounded metadata has no larger row order. Long precision occurs in the
output width, not in a quantum address field. The remainder of a large
precision exponent remains offline mask data. Increase the absolute
cutoff until

```math
\frac{N}{n^2}\ge c_I\bigl(C_a(n+1)+2\bigr),
```

where $`c_I`$ is the nonuniform indicator constant. Then
$`L\gt N/n^2`$ and $`B\ge16L`$ imply the exact sufficient query
condition $`B\ge16c_I(k+2)`$. All helper constants are covered by this
same absolute cutoff, independently of L, b, and W.

The rectangular lemma gives simultaneously on each completed mask

```math
T_j=O\!\left(m_j+\sqrt{Q_jm_j}+\frac{Q_jm_j}{B}\right),\qquad
G_j=O(Q_jm_j),
```

```math
D_j=O\!\left(m_j+\frac{Q_jm_j}{B^2}+\log(k_j+2)\right).
```

For $`Q_j=1`$, emit the fixed Pauli word; these larger bounds still
hold. The additional H layers for Z masks cost $`O(m_j)`$ Cliffords,
already included because $`Q_j\ge1`$. All repeated masks, actual
inverses, and constant amplification multiplicities below are included
in the query list.

## 4. Which queries are rescheduled

The rescheduled list contains exactly the following interfaces.

1. **Growing coefficient masks at each height d.** They satisfy

   ```math
   P_d\asymp(d+1)^2,\qquad S_d=N/2^d,\qquad Q_d=O(S_dP_d),
   ```

   ```math
   m_d=O\!\left(L/P_d+P_d^3\log(P_d+2)\right).
   ```

   There are only a constant number of whole masks per kernel, including
   the scalar constituents and literal X/Z pair. Power/type labels are
   classical table addresses, not separate scalar-source calls.

2. **Cheap row words and their actual column replays.** They have
   $`Q_d=O(S_d)`$ and $`m_d=O(L/P_d+\log(P_d+2))`$. Their source
   and central-guard schedules are unchanged; the cheap rejected
   operator is preserved exactly.

3. **Fine source, attenuation, tag, replay, and disjoint bottom-block
   phase banks.** There are a constant number of banks, each with a
   constant number of phase-gadget masks, $`Q=O(N)`$ and $`m=O(L)`$.
   Their constants depend on the fixed original anchor size and fixed
   error margins, hence are absolute and independent of n and L.
   This includes the real block/source conjugate-signal interfaces,
   not merely the growing upper kernels.

4. **Bounded multiplier/type metadata at each height.** This includes
   the preserved node-type bit loaded once before its entire splitter
   word and unloaded afterward. Row charge is no larger than
   $`O(S_dP_d)`$, output width is polynomial in d, and aggregate
   unweighted output width and call count are polynomial in n. Short
   prefix echoes add only a constant number of these query calls.

The complete coarse mesh, packing, signed-normalization permutation,
center-phase routing, and coarse reflection words are left untouched.
Some contain many native-symbol table lookups; no claim is made that
all queries in the old compiler have polynomial multiplicity. Their
entire separate native count bound is used in Section 6 instead.

The strictly halving precision-refinement words contain no long table.
Their native gates remain untouched. Their precision widths sum to
$`O(L)`$, with additional controls/arithmetic of fixed polynomial size
in n.

## 5. Uniform sums of the rescheduled queries

The inherited native ledger gives, for the preceding list including
constant replays and outer recovery appearances,

```math
A_1:=\sum_jQ_jm_j=O(NL),\qquad
A_2:=\sum_j\sqrt{Q_jm_j}=O(\sqrt{NL}),
```

```math
A_3:=\sum_jm_j=O(L+P(n)),\qquad
A_4:=\sum_j\log(k_j+2)=O(P(n)),
```

for one fixed polynomial P with absolute coefficients. Here and below,
a larger fixed polynomial may absorb finitely many inherited polynomial
native ledgers; it never depends on L or the target data.

For coefficient masks, the row/width product is bounded by

```math
Q_dm_d\le CN2^{-d}\bigl[L+P_d^4\log(P_d+2)\bigr].
```

Both its geometric sum and its square-root geometric sum converge,
giving $`O(NL)`$ and $`O(\sqrt{NL})`$ uniformly for $`L\ge6`$.
Its unweighted widths sum to

```math
O\!\left(L\sum_d\frac1{P_d}
 +\sum_dP_d^3\log(P_d+2)\right)=O(L+P(n)).
```

Cheap words satisfy the same estimates with smaller row/width products.
The finitely many global fine banks contribute $`O(NL)`$,
$`O(\sqrt{NL})`$, and $`O(L)`$. Metadata contributes geometric
row products times fixed polynomials in d, hence $`O(N)`$,
$`O(\sqrt N)`$, and polynomial unweighted cost. The number of these
masks is polynomial in n and $`k_j=O(n)`$, giving $`A_4`$. No factor
proportional to L is introduced into the number of masks.

Sequential composition of the new query schedules suffices:

```math
T_{\rm queries}=O\!\left(\sqrt{NL}+L+P(n)+\frac{NL}{B}\right),
\qquad G_{\rm queries}=O(NL),
```

```math
D_{\rm queries}=O\!\left(L+P(n)+\frac{NL}{B^2}\right).
```

No overlap between separate queries is assumed. Private T-layer
parallelism is used only inside each proved rectangular query.

## 6. Every remaining native word is affordable in depth

The unchanged source loaders and central predicates contribute
$`O(\sum_jm_j+P(n))=O(L+P(n))`$ T gates. Each loader's long word
can run serially. Conjugate lifting costs at most two T gates per
original T gate, preserving this bound even if the shared signal forces
serialization. Strictly halving refinements add $`O(L+P(n))`$ T
gates. Bounded-control metadata consumption and splitter interpretation
add only $`P(n)`$.

The independently proved disjoint-block/coarse ledger gives all remaining
routing, coarse meshes, signed normalization, packing, coarse state-source
words, center-phase echoes, and actual reversals total T-count at most

```math
C_*N^{3/4}P_*(n)
```

and Clifford count at most $`C_*NP_*(n)`$, for a fixed polynomial
$`P_*`$ and absolute $`C_*`$. Source-only subinterfaces have a smaller
square-root bound, which is not needed. These words retain their old
native synthesis and exact controls. Lifting doubles T at most and
adds $`O(T)`$ Cliffords. Original-anchor and bottom recoveries multiply
the already listed calls only by fixed constants.

This uses $`D_T\le T`$ only for this separately bounded remainder;
it does not infer whole-circuit depth from optimal count. Increase the
absolute cutoff until

```math
P(n)+C_*N^{3/4}P_*(n)\le C_0N/n^2\lt C_0L,
```

```math
P(n)+C_*P_*(n)\le C_1N/n^2\lt C_1L.
```

Such a cutoff exists because $`N=2^n`$ and all polynomial degrees
and constants were fixed in the inherited native constructions. It is
independent of L. The remainder therefore has T-count and serial T-depth
$`O(L)`$, and Clifford count $`O(NL)`$. This charges every lifted
gate, predicate, phase correction, and short controlled word without
assuming their internal parallelism.

Using $`B\ge(16/17)b`$ in Section 5 and adding this remainder gives
on the same actual circuit

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad D_T=O\!\left(L+\frac{NL}{b^2}\right).
```

## 7. Error, actual inverses, and finite dimensions

Every changed macro is an exact full-space replacement on the same
base and bank registers. All source precisions, table values, slot
errors, certified-data budgets, full arena guards, mixed gauges, and
source rejection operators retain the inherited construction's values.
Original-anchor recovery excludes the arbitrary signal r and uses its
literal guarded sign; it is unchanged. Neither arbitrary dirty work nor
r is projected or reset. Actual source/column inverses use the reversed
new native word. The complete-isometry error is the existing bound,
at most $`2^{-L}`$, and tensoring any external reference preserves it.
The spare second external clean flag remains untouched above the cutoff.

For n below the single absolute cutoff, use the
[two-clean uniform simultaneous-depth compiler](../../docs/UNIFORM_PRECISION_DEPTH.md#1-the-uniform-theorem-and-its-stronger-low-precision-corollary)
at the same literal $`b\ge17(L+n+7)`$. It has

```math
T=O\!\left(\sqrt{NL}+\frac{NL}{b}+nL\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+nL\right).
```

Since n is bounded by that fixed cutoff, $`nL=O(L)`$ with one absolute
enlarged constant, for every unbounded L. This supplies all three bounds
on the same fallback circuit and retains the full error contract. It is
not a finite-net argument and does not infer fallback depth from count.
This proves the theorem. ∎

The theorem is an upper bound in the stated high-precision range. It
does not assert a matching depth lower bound or a low-precision
complete-frame depth guarantee.
