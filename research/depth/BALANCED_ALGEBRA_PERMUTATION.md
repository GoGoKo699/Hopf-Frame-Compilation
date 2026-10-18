# Balanced native compilation of the HRT algebra permutations

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Historical HRT source and ledger](NATIVE_CLIFFORD_TRANSPORT.md#7-third-level-membership-alone-does-not-bound-tableau-degree) ·
[Returned bilinear query R91](PARALLEL_BILINEAR_WORD_QUERY.md)

For every integer $`r\ge1`$, the HRT permutation has a quadratic inverse
and forward Boolean degree r on $`w=2^r-1`$ arbitrary input qubits. At
$`r=1`$ it is the identity and needs no gates or helpers. Its hierarchy membership does
not supply a shallow native source. This proof supplies a complete
zero-clean native source and its actual inverse with

```math
T\lt48\,3^r,\qquad G\lt60\,3^r,\qquad
D_T\le40\lceil\log_2r\rceil,\qquad
B_{\mathrm{extra}}\lt2\,3^r.                              \qquad(1)
```

Here B counts additional arbitrary dirty qubits, excluding the w input
qubits. Every helper returns exactly, and the error and leakage are zero.
The same simultaneous bounds apply to the actual inverse. These are
conservative bounds, not optimal T-count or Clifford-count claims.
They improve the paid source depth relative to the cardinality-batched
baseline below, at the same asymptotic count and large dirty allocation.
They do not improve the general R96 frame frontier or produce a
bounded-depth, unbounded-degree query family.

The permutation family is due to He, Robitaille and Tan,
[*Characterization of permutation gates in the third level of the
Clifford hierarchy*](https://arxiv.org/html/2510.04993v1), Definition 5.1,
Proposition 5.2 and Lemma 5.7; see also
[the published source](https://doi.org/10.4230/LIPIcs.TQC.2026.2).
The native ingredients here are the retained R91 phase derivative and
literal dirty-register echoes. No external synthesis-priority claim is
made. The source algebra and balanced native construction are proved
below without an uncharged inverse or initialized coefficient bank.

## 1. Algebra and the in-place multiplier

Let $`s\ge1`$ and let
$`R_s=\mathbb F_2[e_1,\ldots,e_s]/(e_1^2,\ldots,e_s^2)`$, including
its unit. Its coefficient space has dimension $`2^s`$. Every element b
with zero constant coefficient satisfies $`b^2=0`$, since cross terms
cancel in characteristic two and every nonconstant basis monomial has
square zero. Input b is arbitrary in that augmentation ideal.
There are $`N\ge1`$ arbitrary input banks, each $`a_\ell\in R_s`$.
We require the complete permutation $`a_\ell\mapsto(1+b)a_\ell`$,
preserving b, with all helpers
returned. These banks can represent the coefficients of a square-zero
ideal supported on fixed left labels. No a-bank is initialized.

Choose one right generator e and write

```math
R_s=R_t\oplus eR_t,\quad t=s-1,\qquad
a_\ell=x_\ell+ey_\ell,\quad b=b_0+ev,\qquad
m=2^t,\quad K=3^t.
```

The constant coefficient of $`b_0`$ is structurally zero; it is omitted from
all native linear forms. It is not represented by a qubit assumed zero.
v is arbitrary. Put $`u=1+b_0`$, so $`u^2=1`$. The desired map is

```math
x_\ell\longmapsto ux_\ell,\qquad
y_\ell\longmapsto uy_\ell+vx_\ell.
```

## 2. Two shears and two copies

For every bank index $`\ell`$, perform these four operations chronologically:

```math
x_\ell\mathrel{+}=b_0y_\ell,\quad
y_\ell\mathrel{+}=x_\ell,\quad
x_\ell\mathrel{+}=b_0y_\ell,\quad
y_\ell\mathrel{+}=x_\ell.
```

All additions and products are in $`R_t`$. A source bank is never a target
of the same shear. Write $`M=M_{b_0}`$ for multiplication by $`b_0`$;
then $`M^2=0`$. Starting from the original x,y, the successively changed
register values are

```math
x+My,\quad y+x+My,\quad (I+M)x,\quad
(I+M)y+Mx.
```

Thus the incomplete diagonal transform is exactly

```math
x_\ell=ux_\ell^{\mathrm{original}},\qquad
y_\ell=uy_\ell^{\mathrm{original}}+b_0x_\ell^{\mathrm{original}}.
```

Keep the displayed extra term until the shared correction. The two
multiplication shears can run over all N banks simultaneously; their
source $`b_0`$ is shared but unchanged.

## 3. A returned dirty correction and five native stages

Borrow one arbitrary d-bank in $`R_t`$, shared by all N pairs. Define

```math
Q:\ d\mathrel{+}=b_0+v+vb_0,\qquad
E:\ y_\ell\mathrel{+}=x_\ell d\quad\text{for every }\ell.
```

Execute chronologically $`Q,E,Q^\dagger,E^\dagger`$.
Q is one completed bilinear XOR plus copies of the b0 and v banks;
the structurally absent constant coefficient of b0 is not copied.
E consists of N simultaneous completed bilinear XORs. Their operators
are involutions, but the final occurrences are emitted as the actual
reversed native inverse words. This echo restores d and adds exactly
$`x_\ell(b_0+v+vb_0)`$ to each $`y_\ell`$. Since $`b_0^2=0`$,

```math
b_0+v+vb_0=(b_0+v)u,\qquad
ux_\ell^{\mathrm{original}}(b_0+v)u
 =(b_0+v)x_\ell^{\mathrm{original}}.
```

The b0 term cancels the retained term from section 2. The full result is

```math
(ux_\ell^{\mathrm{original}},
uy_\ell^{\mathrm{original}}+vx_\ell^{\mathrm{original}}),
```

with no initialized helper and no input-dependent phase.

There are six completed bilinear calls, counting the two diagonal
shears separately and the two E calls each as a simultaneous bank family.
Q commutes with all four operations in section 2: it uses only b0,v,d,
whereas those operations change only the x,y banks. Compile Q's bilinear
part together with the first diagonal shear in one shared native phase
word. Their output registers x,d are disjoint, and every common source
wire belongs to the unchanged b0 bank. Use distinct private phase wires
for every term in this combined word. Its output Hadamards and Clifford
copies are also charged. The five successive non-Clifford stages are

```math
(A_1\text{ together with }Q),\quad A_2,\quad E,\quad
Q^\dagger,\quad E^\dagger,
```

with the first y-copy after A1 and the second after A2. The copies inside
Q follow its bilinear word and can be placed before the first y-copy.
Each stage has at most eight T layers. This schedule does not share a
physical T target between concurrent calls. It implements the same exact
word because the earlier placement of Q is a literal commutation.

Compared with separately removing the b0 term before the correction,
this construction removes one N-bank multiplication family and overlaps
one further call. It improves the paid source constants, while the number
of balanced tree levels and the general frame frontier are unchanged.

## 4. Literal native bilinear word and exact counts

For one nonzero term $`Y_iX_j`$ in output $`Z_k`$, use one private arbitrary
wire h. Let $`\omega=e^{i\pi/4}`$ and let F have phase
$`\omega^{p(h;Y_i,Z_k)}`$, where

```math
p(h;y,z)=h-(h\oplus y)-(h\oplus z)+(h\oplus y\oplus z)
        =(4h-2)yz.
```

A Gray walk emits F using four T/T† gates and four CNOTs:

```math
T_h,\ \mathrm{CX}_{y,h},\ T_h^\dagger,\ \mathrm{CX}_{z,h},\quad
T_h,\ \mathrm{CX}_{y,h},\ T_h^\dagger,\ \mathrm{CX}_{z,h}.
```

For all nonzero multiplication tensor terms in parallel, emit
chronologically $`F^\dagger,L,F,L^\dagger`$, where L toggles each
private h by its corresponding $`X_j`$. Its exact phase is
$`(-1)^{\sum Y_iZ_kX_j}`$; all private wires return.
Hadamards on each output Z before and after convert this to the literal
bilinear XOR. If there are k nonzero tensor terms and q output bits,

```math
T=8k,\qquad G\le10k+2q,\qquad D_T\le8,\qquad
B_{\mathrm{private}}=k.
```

G counts H and CNOT separately. $`T^\dagger`$ is charged as one T gate. The native
word and its literal reverse use only H, CNOT, T and $`T^\dagger`$, so their gate counts
and depth bounds agree. A skipped structural zero term requires no wire.
All F copies use distinct private T targets; shared input/output controls
only add elementary CNOT count and Clifford depth, which remains free.
The source labels Y and X may overlap across different terms: all source
qubits are preserved, every private h is distinct, and the output and
private registers are disjoint from those sources. Thus the lemma also
applies directly to a quadratic shear on one shared source register.

A full product in $`R_t`$ has $`K=3^t`$ nonzero terms, one for each disjoint
ordered pair of subsets. A product with $`b_0`$ has $`K-m`$ terms, since its
constant coefficient is omitted. The two diagonal shears contribute
$`2N(K-m)`$, Q and $`Q^\dagger`$ contribute $`2(K-m)`$, and
E and $`E^\dagger`$ contribute $`2NK`$. Hence, with S counting terms
across all six completed calls,

```math
\begin{aligned}
S&=(4N+2)K-(2N+2)m,\\
T&\le8S,\qquad D_T\le40,\\
G&\le10S+10Nm+8m-2\\
 &= (40N+20)K-(10N+12)m-2.
\end{aligned}
```

The Clifford remainder comprises $`4Nm`$ H gates for the two diagonal
bilinear calls, $`4m`$ for Q and its inverse, $`4Nm`$ for E and its
inverse, and $`2Nm`$ CNOTs for the two diagonal copies. The b0 and v
copies in Q and its inverse cost $`2((m-1)+m)=4m-2`$ CNOTs.
These are upper counts when some structural-zero products or outputs
are omitted. Combining two calls in one stage preserves this gate-count
sum and their literal phases.

The parallel first stage needs $`(N+1)(K-m)`$ private phase wires; an E
stage needs NK. All later stages reuse that pool only after its previous
stage has returned it. The arbitrary d-bank lies outside this pool, so

```math
B_{\mathrm{extra}}=m+\max\{NK,(N+1)(K-m)\}
 \le(N+1)K\le2NK,\qquad a_{\mathrm{clean}}=0.
```

The first inequality uses $`m\le K`$ and N at least one. The actual
reversed native circuit reverses the five stages, with the same counts,
depth and simultaneous allocation. It returns every dirty bit on every
input; no independently synthesized substitute inverse is used.

No approximation, discarded branch, relative phase, or leakage occurs.
The proof is equality on the entire computational basis, including all
arbitrary helper columns, and hence equality on arbitrary superpositions
and external references.

## 5. The limited determinant obstruction and its escape

If whole a-banks were treated as rank-one free $`R_s`$ modules, bilinear
bank shears and copies have determinant one over $`R_s`$. The desired
$`\mathrm{diag}(1+b,1,\ldots,1)`$ has determinant $`1+b`$, so for
nonzero b that restricted representation
cannot return every helper. Splitting one generator changes the scalar
ring to $`R_t`$ and doubles the module rank. The desired matrix is

```math
\begin{pmatrix}u&0\\v&u\end{pmatrix},
```

with determinant $`u^2=1`$. The explicit construction above
therefore escapes exactly that representation-specific obstruction.
It is not a lower bound for the unrestricted native model.

## 6. Balanced HRT construction

The input c has one coefficient $`c_S`$ for every nonempty subset
$`S\subseteq[r]`$.
All these coefficients are arbitrary logical data, not initialized work.
The native target is the complete permutation

```math
F_r|c\rangle=|f(c)\rangle,\qquad
f(c)=\prod_{\varnothing\ne S\subseteq[r]}(1+c_Se_S)-1.
                                                               \qquad(2)
```

The source's quadratic inverse is

```math
k(a)_S=a_S+
\sum_{\substack{\{I,J\}:I,J\ne\varnothing\\I\sqcup J=S}}
a_Ia_J.                                                       \qquad(3)
```

The sum is over unordered pairs. The
[retained algebra proof](NATIVE_CLIFFORD_TRANSPORT.md#7-third-level-membership-alone-does-not-bound-tableau-degree)
establishes $`k=f^{-1}`$ on every coefficient vector. The construction
below directly supplies f; its actual reversed native word supplies k.

Partition the original c coefficients by their minimum subset label:

```math
v_i=\sum_{\min S=i}c_Se_S.
```

All products between two basis monomials in one $`v_i`$ vanish, so the HRT
map is $`f(c)=\prod_{i=1}^r(1+v_i)-1`$. Maintain each contiguous interval
as its product-minus-one, encoded in the physical coefficient banks whose
minimum labels belong to that interval. If left child is a and right
child b, their merged value is $`a+b+ab`$. The update $`a\mapsto a+ab`$,
with b unchanged,
has precisely that encoding. It is an invertible complete native update,
not an initialized accumulator.

For a node whose left interval starts at l and has length L, put
$`j=l+L`$ and $`s=r-j+1`$. Each nonempty subset of the L left labels
indexes one full a-bank over variables $`j,\ldots,r`$. Thus
$`N=2^L-1`$. The right child b
belongs to the augmentation ideal of that suffix algebra. Its coefficients
whose minimum labels lie after the right child's endpoint are structural
zeros and must be omitted, not replaced by unrelated physical qubits.
Using the full suffix tensor counts is a safe upper bound.

More explicitly, the physical coefficient qubit for an a-bank indexed
by a nonempty $`A\subseteq\{l,\ldots,j-1\}`$ and suffix mask
$`B\subseteq\{j,\ldots,r\}`$ is the original qubit $`c_{A\cup B}`$,
now holding its current transformed coefficient. These are all existing
data qubits. For b, only masks whose minimum lies in the right child's
interval have a physical source qubit. A requested b coefficient outside
that interval is the formal zero polynomial: omit its associated copy
or product term from the native word. Splitting by generator $`e_j`$
likewise supplies an existing source for each nonzero $`b_0`$ or v
coefficient. No omitted coefficient is represented by an arbitrary
qubit treated as zero. The d-bank and phase-helper pool are separately
allocated arbitrary qubits and return for every value.

Different nodes at a tree level operate on disjoint physical coefficient
banks classified by their minimum labels. Shared formal generator names
do not imply shared qubits. Give these nodes disjoint private helper pools.
A balanced tree has height $`\lceil\log_2r\rceil`$; each level has depth
at most 40 under the five-stage schedule.

For $`t=s-1=r-l-L`$,

```math
NK\lt2^L3^{r-l-L}=3^{r-l}(2/3)^L.
```

For fixed l, its internal ancestor nodes have distinct positive L, so
the sum of this expression is at most $`2\cdot3^{r-l}`$. Summing over
$`l=1,\ldots,r-1`$ gives $`\sum_{\mathrm{nodes}}NK\lt3^r`$.
Since $`m\le K`$ and $`N\ge1`$,

```math
\begin{aligned}
T&\le48\sum_{\mathrm{nodes}}NK\lt48\,3^r,\\
G&\le60\sum_{\mathrm{nodes}}NK\lt60\,3^r,\\
B_{\mathrm{extra}}&\le2\sum_{\mathrm{nodes}}NK\lt2\,3^r,\\
D_T&\le40\lceil\log_2r\rceil,\qquad a_{\mathrm{clean}}=0.
\end{aligned}
```

These are deliberately conservative full-circuit bounds, establishing
(1). Allocating private helpers for every node already fits the stated
width; reusing returned pools can only reduce it. The actual
reversed circuit implements $`k=f^{-1}`$ with the same bounds and returns
all dirty work. This improves the paid HRT source depth to logarithmic in
r, with $`O(3^r)`$ count and dirty width. It does not give constant source
depth or improve the R96 complete-frame frontier.

## 7. A returned query and the fair baseline

A literal returned query may conjugate a target CNOT by F and its actual
inverse: chronologically $`F,\mathrm{CX}_{[r],z},F^\dagger`$ toggles
z by $`f(c)_{[r]}`$, restores all c and helpers, and has T-depth at most
$`80\lceil\log_2r\rceil`$, twice the
source counts plus one CNOT. Its Boolean degree is exactly r: on
singleton-supported c the top coefficient is the product of the r
singleton coordinates, and no nonzero algebra product has more than r
nonconstant factors. This is a
matching mechanism stress test for logarithmic degree-growth statements,
not a bounded-depth unbounded-degree counterexample.

This query treats **all w coefficient qubits as its address**. Restricting
the nonsingleton c coefficients to zero would algebraically yield an
r-input AND, but those zeros are additional clean promises. No r-input
AND with the remaining coefficient bank arbitrary and returned follows
from this construction. The native source itself acts on every c column.

For the baseline comparison assume $`r\ge2`$; the $`r=1`$ identity was
handled above. The old serial source ledger is not a fair depth baseline. Equation (3)
can be compiled by decreasing target cardinality, as in the retained
Toffoli staircase. At a fixed target cardinality s, every control lies
on a strict subset, so no target in this batch is a control in it. The
entire batch is therefore a sparse quadratic XOR word and admits the
eight-layer returned construction of section 4. Let

```math
M_r=\frac{3^r-2^{r+1}+1}{2},\qquad
M_{r,s}=\binom rs(2^{s-1}-1),\qquad 2\le s\le r.
                                                               \qquad(4)
```

These count all unordered disjoint nonempty pairs, and those whose
union has size s, respectively. With one private dirty wire per term
within a cardinality batch, and a pool reused between batches, the
quadratic inverse and actual forward inverse have the complete bounds

```math
T=8M_r,\qquad G\le10M_r+2(2^r-r-1),\qquad
D_T\le8(r-1),\qquad
B=\max_{2\le s\le r}M_{r,s},\qquad a_{\mathrm{clean}}=0.
                                                               \qquad(5)
```

The output Hadamards occur once per nonsingleton coefficient at each
end of its own batch. Every private helper returns before the next
batch. For an arbitrary pool size B at least one, split the terms of
each cardinality batch into groups of at most B, keeping the output
Hadamards outside all its groups. Counts in (5) are unchanged, while

```math
D_T\le8\sum_{s=2}^r\left\lceil M_{r,s}/B\right\rceil
\le8(M_r/B+r-1).                                          \qquad(6)
```

This comparison charges the fairly simplified parity walk and includes
all helper return and inverses. At large dirty width the new balanced
source improves the dependence on r from linear to logarithmic, with
T and G still $`O(3^r)`$. It does not assert better constants for small
r: either complete circuit may be selected. Since $`r=\log_2(w+1)`$,
the new depth is $`O(\log\log(w+1))`$ in the source's input width.

## 8. What the construction resolves and what remains

The first unpaid operation in the balanced product idea was in-place
multiplication by $`1+b`$, not the out-of-place bilinear query. Sections
1–4 pay that operation. The whole-bank determinant obstruction in
section 5 fails after splitting one generator, exactly as the explicit
native construction demonstrates. It must not be promoted to an
unrestricted native lower bound.

The source depth still grows as $`\log r`$. Thus this family does not
refute a possible native depth bound of the form
$`\deg f\le2^{O(D_T)}`$, nor prove such a bound. The next constructive
question is whether the balanced merge depth can be reduced while
retaining the complete physical coefficient operator and its returned
dirty banks. A simultaneous logarithmic degree bound must also allow
this family. An arbitrary-table restored-address query with the R98
count and workspace budgets remains a separate missing construction.

## 9. Finite evidence

The seven [exact checks](../../tests/test_balanced_algebra_permutation.py)
include every column of the five-qubit native multiplier and its actual
adjoint, both disjoint and shared-source sparse phase words on arbitrary
dirty inputs, the shared correction bank across multiple a-banks, and every
physical column of the concurrent first-stage products with a shared b0 control.
Classical coefficient tests verify the balanced permutation and its
reversed tree against (2)–(3), exhaustively through $`r=3`$ and on seeded
inputs at $`r=4,5`$. Symbolic checks audit counts and simultaneous helper
allocation through $`r=64`$. These checks do not substitute for the native
source and work-return proof above.
