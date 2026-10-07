# Retained-source fusion and stabilizer completion

**Status:** exact full-frame factors and a shallow completion; ordered transport remains linear in group height. These refinements do not improve Result D. Sections 1–8 refer to the [unary-source proof](../../docs/UNARY_PHASE_GRADIENT.md); original numbering is retained below.

## 9. Retained-source fusion without a larger modulus

The fixed four-mode magic basis and eight-mode Bell projector in
[the existing local analysis](../endpoint/ENDPOINT_TREE_TRANSPORT.md#7-a-native-four-mode-benchmark)
also give exact identities for the retained unary source. Balanced
Laurent factors keep the original modulus q and preserve every source
input. The four- and eight-mode words use two and three selected-shift
rounds respectively. Section 10 batches their stabilizer completion for
growing groups; the remaining transport retains its linear-depth allowance.

### One cyclic source and balanced four-mode factors

Let S be the physical cyclic permutation of the q source wires. It
satisfies $`S^q=I`$. Use the commutative coefficient ring

```math
\mathcal R_q=\mathbb Z[i,1/2][z,z^{-1}]/(z^q-1),\qquad z\mapsto S,
```

where adjoints conjugate i and send $`z\mapsto z^{-1}`$. Define

```math
c_a=\frac{z^a+z^{-a}}2,\qquad
s_a=\frac{z^{-a}-z^a}{2i},\qquad
G_{ij}(a)\big|_{i,j}=\begin{pmatrix}c_a&-s_a\\s_a&c_a\end{pmatrix}.
```

On source character k, $`z=\omega_q^{-k}`$, these are the real
rotation entries for angle $`2\pi ka/q`$. The unary one-hot orbit
contains all q characters, not just the prepared character $`k=1`$.
Polynomial identities modulo $`z^q-1`$ therefore act identically in
every representation of S, including arbitrary physical source states
outside that orbit and their reference correlations.

Use root-most-significant-bit mode order and
$`W_4=G_{23}(c)G_{01}(b)G_{02}(a)`$. Let M have the fixed columns
$`[\Phi^+,i\Psi^+,i\Phi^-,\Psi^-]`$ from the linked magic-basis
analysis, and put $`\Pi_P^\pm=(I\pm P)/2`$. Its Clifford entries
use the scalar extension containing $`1/\sqrt2`$; the conjugated
Laurent coefficients below lie in $`\mathcal R_q`$ itself. Define

```math
A_r=\Pi_Z^++z^a\Pi_Z^-,\qquad
B_r=z^{-a}\Pi_Z^++\Pi_Z^-,
```

```math
A_c=\Pi_X^++z^{b+c}\Pi_X^-,\qquad
B_c=z^{-b}\Pi_X^++z^{-c}\Pi_X^-.
```

Direct multiplication gives the exact identity

```math
MW_4M^\dagger=(A_cA_r)\otimes_{\mathcal R_q}(B_cB_r).
```

The tensor product concerns the two logical factors: its entries multiply
Laurent operators on the same source. It does not allocate independent
sources. The two factor determinants are $`z^{a+b+c}`$ and
$`z^{-(a+b+c)}`$; their phases cancel literally. Dropping an unmatched
scalar from either factor changes the retained-source word. Consistent
paired normalizations can cancel, but their separate half-powers need not
belong to the original ring; the balanced factors avoid those choices.

Implement the paired root and child factors jointly. Their exponent
lists, in the Z and X product bases respectively, are

```math
(-a,0,0,a),\qquad(-b,-c,c,b).
```

Each list is one selected source shift with a fixed Clifford boundary.
The selected program uses the original one-hot labels, their reversed
indices, or the fixed zero-shift label. In particular the apparent
$`b+c`$ inside an individual factor requires no summed-label query or
unpriced convolution. Computing and erasing each constant-size selected
row uses $`O(q)`$ work and gates by the existing mask interface. Two
rounds give the full four-mode word, retaining all logical columns.

### Eight modes use one joint root shift at the original modulus

Conjugate the two four-mode children by $`K=I\otimes M`$. With
$`P_{\rm Bell}=|\Phi^+\rangle\langle\Phi^+|`$, the root becomes

```math
U=I+(c_a-1)I\otimes P_{\rm Bell}-i s_aY\otimes P_{\rm Bell}.
```

Its four commuting Pauli terms are
$`Q_0=YII`$, $`Q_1=YXX`$, $`Q_2=-YYY`$, and $`Q_3=YZZ`$.
Diagonalize the independent commuting generators $`YII,IXX,IZZ`$
with a fixed Clifford. If their signs are $`y,x,z_B\in\{1,-1\}`$,
the four terms have joint signs $`(y,yx,yxz_B,yz_B)`$. The total source
exponent is consequently

```math
\frac{ay(1+x)(1+z_B)}4\in\{0,a,-a\}.
```

Thus the root is one signed selected shift, using the original a and
the original q. There are no four independent quarter-angle sources.
The factor one fourth expresses a projector whose joint eigenvalues
are zero or one; it does not demand a fractional source index.
The fixed Clifford and its actual inverse preserve literal phase.

The two children have disjoint logical branches. Select the appropriate
child's label in each of their two paired rounds, so their combined
action occupies two slots on the shared source. Along with the root,
the eight-mode word therefore has three slots. This construction closes
the stated local word; it does not turn the child/root family into
globally commuting independent factors.

### An exact Bell-stabilizing factor with its native cost

For child $`s\in\{0,1\}`$, write its balanced factors as
$`C_s=A_s\otimes_{\mathcal R_q}B_s`$, and define

```math
V_s=A_sB_s^{\mathsf T},\qquad
K_s=B_s^{-\mathsf T}\otimes_{\mathcal R_q}B_s,\qquad
C_s=(V_s\otimes I)K_s.
```

Transpose acts only on logical two-by-two indices; it neither conjugates
nor inverts the source variable. The Bell identity gives
$`K_s(|\Phi^+\rangle\otimes|\psi\rangle)=|\Phi^+\rangle\otimes|\psi\rangle`$
for every source input. Put
$`\mathcal C=\mathrm{diag}(C_0,C_1)`$,
$`\mathcal V=\mathrm{diag}(V_0,V_1)\otimes I`$, and
$`\mathcal K=\mathrm{diag}(K_0,K_1)`$.
The branch-controlled $`\mathcal K`$ commutes with the root U even
for unequal children: it is identity on the Bell sector, and U is
identity on its complement. Source coefficients commute throughout.
Consequently, on every logical and physical source input,

```math
\mathcal C U=\mathcal VU\mathcal K,\qquad
\mathcal C U\mathcal C^\dagger=\mathcal VU\mathcal V^\dagger.
```

These factors have native schedules if the program also preloads the
derived labels $`b+c`$ and $`b-c`$ modulo q for each child. This is a
constant-factor expansion of the bounded four/eight-mode table. Charge
its exact query, actual unload, and logical storage; no quantum sum
lookup or convolution is free. The selected readouts still cost O(q).

Indeed $`V=A_c(A_rB_r)B_c`$ has three rounds: B_c in the X basis,
the Z exponents $`(-a,a)`$, and A_c in the X basis with derived label
$`b+c`$. K has Z exponents $`(0,a,-a,0)`$ and X exponents
$`(0,b-c,c-b,0)`$, hence two rounds. Multiplexing the child branch
keeps these slot counts. The full product $`\mathcal VU\mathcal K`$
uses six rounds against three for $`\mathcal C U`$; the conjugated
word uses seven against five. These are conservative upper schedules,
not minimum-round bounds. The exact extraction supplies no saving by
itself, and it does not make a generic recursively transported block
a single-source scalar rotation.

There is a limited transport identity on the incoming root-times-Bell
subspace. Let $`R_a(z)`$ be the two-by-two rotation defined above, and put

```math
F=\mathrm{diag}(V_0,V_1)(R_a(z)\otimes I_2)\in\mathrm{SU}(4,\mathcal R_q),
```

```math
\bigl[\mathcal VU-(F\otimes I_2)\bigr]\,(I_2\otimes P_{\rm Bell})=0.
```

Here F acts on the root and first Bell qubit. This is a boundary
isometry, not a replacement for the full U on the
Bell-orthogonal sector. Its logical block coefficients need not commute.
Legal child triples $`(0,1,0)`$ and $`(1,0,0)`$ give respectively
$`V_0=z^{-1}\Pi_X^++z\Pi_X^-`$ and
$`V_1=z^{-1}\Pi_Z^++z\Pi_Z^-`$. At $`q=8,k=1`$ their commutator
is $`iY`$, of norm one. A scalar two-by-two factorization therefore
cannot simply be iterated on these matrix-valued blocks; a new
composition rule and its full native ledger would be needed.

### Guard, resources, and the limit of the local simplification

At each slot the completed selected shift is diagonal in its current
logical selector basis. Compute its selected row, apply the guarded
shift, and erase the row before changing that basis. On the active
sector all selected-program work is conditionally zero, and the shift
returns its temporary word for every source input. On $`h=0`$ the
central guarded primitive is identity on arbitrary work. Actual selector
inverses and Clifford boundaries then cancel. These identities introduce
no new source error or supplied initialized state.

For growing groups, the retained compiler still gives the conservative
slot recurrence $`d_g=d_{g-1}+1`$, hence $`d_g=g`$. Its existing
ledger, excluding query/predicate pairs and source boundaries, is

```math
T=O(q2^g+gR),\qquad G=O(q2^g+gqR),\qquad D_T=O(g),
```

```math
w=O(q2^g+R+q\log q+g),\qquad R=3^\ell,\quad q=2^\ell.
```

This is the retained schedule and its existing selector accounting.
One may conjugate its full word by explicitly charged fixed Cliffords,
or restore the original logical basis before using those selectors.
Root diagonalization alone does not implement a new recursively
transformed selector-cache schedule. The displayed four/eight-mode
readouts have constant-size Clifford boundaries and use existing native
shift primitives; no larger-modulus compiler is assumed.

A single fixed-basis diagonal-shift family would commute for all program
choices, because its source coefficients commute. But
$`W_4(1,0,0)=G_{02}(1)`$ and $`W_4(0,1,0)=G_{01}(1)`$ do not
commute: their product difference has entry $`s_1^2`$ at position
$`(1,2)`$. At $`q=8,k=1`$ that entry is one half. This excludes
only a single-round template with the same fixed basis for every row.
It is not a minimum-round theorem for arbitrary circuits or a lower
bound on the complete-frame depth. The local forms above retain the
current linear group-depth allowance and prove no new frame theorem.

The seven [retained-source checks](../../tests/test_retained_source_fusion.py)
include exact Gaussian-dyadic Laurent identities for all 512 q=8 label
triples, their eight character evaluations, and small full physical q=4
source banks. They also check distinct eight-mode children, determinant
phases, integer root exponents, inactive Clifford cancellation, Bell
extraction, and Bell-only transport. Wrong-transpose, omitted-factor,
and full-space transport replacements fail. Native
Clifford matrices and numerical character/physical-bank evaluations are
separate from the exact ring checks. The fixture emits no scalable
selected-program, private-work, or full-group native circuit; the resource
and arbitrary-work return claims rely on the charged primitives in the linked unary-source proof and the accounting here.

## 10. Two source shifts for the stabilizer completion

For an even-height group $`g=2m`$, the stabilizer factors from Section 9
can be implemented together with two selected source shifts and
$`O(\log(g+2))`$ T-depth. The row predicates are computed before any
basis change and retained until the completion returns to the original
basis. This improves the growing completion component; the available
ordered transport schedule remains a separate $`O(g)`$ cost.

### Factor the full group without discarding columns

Index logical pairs by $`j=0,\ldots,m-1`$ from top to bottom. At pair j,
let u be the upper $`2j`$ bits. Its two-level word $`G_j`$ is multiplexed
over all such prefixes and controlled by all later pairs being zero.
Its row parameters
$`a_{j,u},b_{j,u},c_{j,u}`$ may also depend on the preserved external
prefix. Section 9 gives, on the active pair,

```math
G_j=T_j\overline K_j,\qquad
T_j=M_j^\dagger(V_{j,u}\otimes I)M_j,\qquad
\overline K_j=M_j^\dagger K_{j,u}M_j,
```

with both factors extended by identity outside the same prefix/suffix
controls. Since K fixes the Bell vector pointwise for every source,
$`\overline K_j`$ fixes its own 00 state pointwise. Thus
$`\overline K_j-I`$ is supported on

```math
Q_j=[\text{pair }j\ne00][\text{later pairs}=0].
```

For $`i\lt j`$, the support of $`G_i-I`$ requires pair j to be zero,
so it is orthogonal to $`Q_j`$. Hence $`[\overline K_j,G_i]=0`$.
The Q_j are also mutually orthogonal, so the stabilizer factors commute.
Moving them past the earlier completed words proves

```math
W_g=G_{m-1}\cdots G_0
=T_{m-1}\cdots T_0\,\mathcal K_g,\qquad
\mathcal K_g=\prod_{j=0}^{m-1}\overline K_j.
```

This is a full logical/source identity, with the completion applied
first. It does not discard a complementary column or restrict the
incoming logical state to a vacuum sector.

### Cache original row predicates before the global basis changes

Refine the support into
$`Q_{j,u}=[\text{upper}=u]Q_j`$. These projectors are mutually
orthogonal and partition every nonzero computational-basis string. There
are $`\sum_j4^j=(2^g-1)/3`$ rows. Cache each Boolean value using

```math
[\text{upper}=u,\text{later}=0]
\ \oplus\
[\text{upper}=u,\text{own pair}=00,\text{later}=0].
```

Compute both conjunctions with private copied literals and balanced
AND trees, copy their roots into the zero cache bit, and reverse the
actual trees and copies. An empty conjunction is a literal X. Parallel
trees over all rows use $`O(g2^g)`$ T and Clifford gates and private
work, with $`O(\log(g+2))`$ T-depth. On active zero-work inputs the
tree scratch returns to zero, leaving only the row cache. No all-input
tensor-return claim for this scratch is needed on the inactive sector.

Let $`M_{\rm all}=\bigotimes_jM_j`$ and $`H_{\rm all}=H^{\otimes g}`$.
With the original row cache retained, execute chronologically:

1. Apply $`M_{\rm all}`$.
2. Compute the q-bit selected program from the cached row and that row's
   own current pair, using Z-basis exponents $`(0,a,-a,0)`$. Apply the
   guarded source shift and erase the selected program.
3. Apply $`H_{\rm all}`$. Compute the selected program using the same
   cached row and its current pair, with exponents
   $`(0,b-c,c-b,0)`$. Apply the guarded shift and erase this program.
4. Apply $`H_{\rm all}`$, then $`M_{\rm all}^\dagger`$, and finally
   the actual inverse of the complete row-cache word.

The selected program defaults to $`e_0`$ if no cached row is active.
Each row readout uses its cached bit and only its own two current logical
bits. It never rereads the upper prefix after M or H has changed basis.
For example, the selected word can be formed as $`e_0`$ XOR the
row-and-pair masks times $`e_a\oplus e_0`$ for their respective labels.
The constant-degree private-copy masks cost $`O(q2^g)`$ gates/work
and constant T-depth by the existing selected-program primitive.

On a valid cached row $`(j,u)`$, each completed original-basis axis
acts only on pair j, preserves its 00-orthogonal complement, and leaves
the upper prefix and later zero pairs fixed; the other pairs' basis
changes cancel. Thus it preserves that Q row and, by orthogonality of
the row-controlled branches, all Q rows. Equivalently, the conjugated
rows remain valid through the two axes, including coherent branch and
source superpositions. The final M inverse restores the original basis,
so the actual cache inverse clears every active cache bit and scratch.
This is a valid-cache-embedding statement, not commutation of an
unconditional pair operation with every row projector on arbitrary caches.

On $`h=0`$ both native shifts are identity for arbitrary program, source,
and private-work inputs. Their surrounding program words cancel as
actual inverses, followed by the H and M boundaries. Thus the entire
middle is identity even if the initial cache computation entangled
arbitrary scratch with the logical input. Its actual inverse then
cancels that computation. This proves literal full inactive identity;
it never interprets an inactive dirty cache as a valid predicate.

### Charge the program, workspace, and remaining transport

The derived $`b-c\bmod q`$ one-hot rows are compiled into the classical
table and loaded with the original rows. Their exact quantum query,
actual unload, and $`O(q2^g)`$ logical program storage are charged.
All row caches, tree/copy scratch, selected programs, and source-shift
work are disjoint while live; returned private pools may then be reused.

Before absorption, add $`O(g2^g)`$ for the row cache to the count/work
bounds below. The condition $`g\le q`$ absorbs this cost into the
selected-program budget and includes the inherited range $`g\le\log_2q`$.
With $`R=3^\ell`$ and $`q=2^\ell`$, the completion alone has

```math
T,w_{\rm extra}=O(q2^g+R),\qquad
G=O(q2^g+qR),\qquad D_T=O(\log(g+2)).
```

These costs exclude the program query pair, outer h-predicates, and
source preparation/return. The q-bit source core and inherited
$`O(q\log q)`$ preparation/index bank are separate reservations.
The total uses the same asymptotic conditional group workspace with a
possibly larger constant; no previous literal allocation constant or
new full-frame width threshold is asserted. All-source exactness adds
no approximation error.

The available transport schedule uses three shifts per pair with the
charged derived rows of Section 9. Each block retains its original
lower-zero guard and changes only its own pair. Consuming suffix enables
before those targets change, then extending prefix selectors from their
final values, extends the [incremental selector schedule](../../docs/GROUPED_PROGRAM_PREFETCH.md#10-amortized-local-selectors-and-suffix-enables)
with constant overhead per pair. The factorized whole-group word has $`3m+2`$ slots
and $`O(g)`$ depth, distinct from the original compiler's g slots.
The full-group and complete-frame asymptotic depth are unchanged.
The new four-target, q=8 fixture checks five exact row sectors plus the
vacuum, both axes' Q invariance, the serial/batched completion and full
group factorization. Recomputing rows in the transient basis or omitting
the upper pair's suffix guard fails; the latter can preserve the first
column while corrupting complementary columns. These are exact ring
checks, not an arbitrary-cache simulation or a scalable native emitter.
