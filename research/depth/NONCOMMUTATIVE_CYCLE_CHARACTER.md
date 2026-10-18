# A noncommutative cycle character with complete dirty return

[Commutative character obstruction](NONLINEAR_CHARACTER_TRANSPORT.md) ·
[Native query baseline R91](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Conditional frame transfer R98](QUERY_DEPTH_FRAME_TRANSFER.md)

A noncommutative nilpotent coefficient algebra supports a degree-d Boolean
character on only $`d^2`$ independently arbitrary work bits. Its character
commutator is an exact full-space identity, and the construction below pays
for a controlled source, read-only-address transport, both actual inverses,
and arbitrary-output XOR. Thus the logarithmic character-degree bound for
commutative coefficient algebras does not extend to this representation.

The character commutator is a representation escape; its demonstrated
output is a single conjunction. The explicit parallel native word has
$`T,G=O(d^3)`$, extra dirty width $`d^2+d`$, and
$`D_T\le64(d-1)`$. The retained strictly read-only conjunction already
has logarithmic depth with smaller polynomial count and width; the ordinary
completed-circuit MCX baseline is smaller still. Section 6 gives their
distinct contracts. No arbitrary-table reader, sublogarithmic depth law,
or R91/R96/R98 improvement follows. Section 5d separately gives the
structured inverse-word query and cycle source logarithmic T-depth, with
cubic count and arbitrary dirty width. Its matrix inputs change internally,
and the full conjunction's separate transport cost remains unchanged.

## 1. The algebra and a nonlinear Boolean character

Fix $`d\ge1`$. Indices i are in $`\mathbb Z/d\mathbb Z`$.
Let $`A_d`$ have basis

```math
e_{i,k},\qquad 0\le i\lt d,\quad1\le k\le d,
```

with multiplication

```math
e_{i,k}e_{j,l}=
\begin{cases}
e_{i,k+l},&j=i+k\pmod d,\quad k+l\le d,\\
0,&\text{otherwise}.
\end{cases}                                               \qquad(1)
```

These are positive-length directed paths around an oriented d-cycle,
with paths longer than d set to zero. Concatenation proves associativity;
the truncation is an ideal. Hence $`\dim A_d=d^2`$ and
$`A_d^{d+1}=0`$. Every basis element squares to zero: a nonzero square
would require $`k=0\pmod d`$, hence k equal to d, but its resulting
length 2d is truncated.

The unitization gives the finite group $`G_d=1+A_d`$. Each arbitrary
$`d^2`$-bit word represents exactly one unit

```math
u_z=1+\sum_{i,k}z_{i,k}e_{i,k}.
```

There is no initialized coordinate, valid-code promise, padding sector,
or supplied eigenstate. Inverses are finite geometric series.

Distinguish the edge from vertex zero to vertex one. Let $`\nu(i,k)`$
be its number of occurrences in the path $`e_{i,k}`$, modulo two, and set

```math
\delta(e_{i,k})=\nu(i,k)e_{i,k},\qquad\delta(1)=0,
\qquad \tau(z)=\sum_i[e_{i,d}]z.                         \qquad(2)
```

Concatenation adds edge counts, so delta is a derivation. The functional
tau is a cyclic trace: $`\tau(vw)=\tau(wv)`$. For basis paths whose
lengths do not sum to d, both evaluations vanish. If their lengths sum
to d, a valid concatenation closes around the cycle exactly when its
cyclic rotation does, and both traces equal one. The identity extends
linearly and to the unitization with $`\tau(1)=0`$. In particular,
$`\tau(rsr^{-1})=\tau(s)`$ for every unit r.

Define

```math
\chi(u)=\tau(\delta(u)u^{-1})\quad\text{in }\mathbb F_2.
                                                               \qquad(3)
```

The product rule gives

```math
\delta(uv)(uv)^{-1}
 =\delta(u)u^{-1}+u\delta(v)v^{-1}u^{-1}.
```

Applying tau proves $`\chi(uv)=\chi(u)+\chi(v)`$. Thus chi is a
Boolean group character even though the group is noncommutative.

For the read-only address $`x\in\mathbb F_2^d`$, put

```math
N_x=\sum_{i=0}^{d-1}x_i e_{i,1},\qquad a_x=1+N_x.
```

Here $`\delta(a_x)=x_0e_{0,1}`$ and
$`a_x^{-1}=\sum_{r=0}^dN_x^r`$. Only the term with r equal to
$`d-1`$ reaches length d after multiplication by $`\delta(a_x)`$.
The unique continuation of the distinguished edge traverses every edge
once, so

```math
\chi(a_x)=\tau(x_0e_{0,1}N_x^{d-1})
         =\prod_{i=0}^{d-1}x_i.                           \qquad(4)
```

This includes even d and d equal to one. At d equal to seven, for example,
the degree is seven on 49 coefficients, exceeding
$`1+\lfloor\log_2 49\rfloor=6`$. The
[commutative-algebra theorem](NONLINEAR_CHARACTER_TRANSPORT.md#the-same-width-barrier-for-every-commutative-nilpotent-algebra)
therefore requires its commutativity hypothesis.

## 2. Exact return, controlled sources, and arbitrary outputs

Let $`P_x|u\rangle=|a_xu\rangle`$ and
$`D_\chi|u\rangle=(-1)^{\chi(u)}|u\rangle`$. Then

```math
D_\chi P_xD_\chi^\dagger P_x^\dagger
       =(-1)^{\chi(a_x)}I_W.                              \qquad(5)
```

The rightmost transport sends u to $`a_x^{-1}u`$. The two source phases
sum to $`\chi(a_x^{-1}u)+\chi(u)=\chi(a_x)`$, and the final work
word is u. This proof uses the character law, not commutativity of the
transport group.

For a separate arbitrary output bit h use the controlled source

```math
D_h|h,u\rangle=(-1)^{h\chi(u)}|h,u\rangle.
```

The analogous commutator has phase $`(-1)^{h\prod_i x_i}`$ and
identity action on all coefficient work. Conjugate by Hadamards on h.
The complete native word, in chronological order, is

```math
H_h,\quad P_x^\dagger,\quad D_h^\dagger,\quad
P_x,\quad D_h,\quad H_h.                                \qquad(6)
```

It implements $`h\mapsto h\oplus\prod_i x_i`$ for arbitrary h.
Every dagger denotes the actual reversed native inverse word. Equality
holds on every physical input column, hence on arbitrary superpositions
and with every external reference system. No phase is discarded.

## 3. An in-place native character source

Order basis paths by increasing length k, then increasing start i. Every
unit has a unique ordered factorization into $`1+c_{i,k}e_{i,k}`$.
A reversible triangular word F computes those factor bits on the same
$`d^2`$ arbitrary wires.

At pivot $`p=e_{i,k}`$, let its current coefficient be $`c=z_{i,k}`$.
Conceptually cancel it on the **left** of the residual unit by
$`1+cp`$, its own inverse because $`p^2=0`$. Retain c physically as
the factor bit. For each path $`q=e_{i+k,l}`$ **later than p in the
basis order**, with $`k+l\le d`$, execute

```math
z_{i,k+l}\longmapsto z_{i,k+l}+z_{i,k}z_{i+k,l}.           \qquad(7)
```

Previously processed physical coordinates store factor bits, not residual
coefficients. Their residual coefficients are zero, so they must be
excluded from (7). Every target has greater length than the pivot and
is unprocessed. Within one pivot, sources start at $`i+k`$ and targets
at i; these starts differ whenever an update exists. Thus no updated
target becomes a source for that pivot, and all three Toffoli wires are
distinct. The updates implement the conceptual residual multiplication
exactly, while storing the extracted c instead of its canceled zero.
Induction proves F on all inputs and reconstructs

```math
u=(1+c_1p_1)(1+c_2p_2)\cdots(1+c_wp_w)
```

in the chosen order. The actual reversed native word is its full inverse.

For a basis generator, (3) gives
$`\chi(1+ce_{i,k})=0`$ for $`k\lt d`$ and c for k equal to d.
Consequently $`\chi(u)=\sum_i c_{i,d}`$. Compute F, apply CZ from h
to each of the d top factor bits, then execute the actual inverse of F.
This realizes $`D_h`$ on every coefficient and output input. No T gate
is placed on h.

Put $`r=\lfloor d/2\rfloor`$. The lexicographic F has exactly

```math
N_d=d\left\lfloor\frac{(d-1)^2}{4}\right\rfloor
             +rd-\frac{r(r+1)}2                         \qquad(8)
```

Toffolis. Sources with $`l\gt k`$ contribute
$`d\sum_k\max\{d-2k,0\}`$, the first term. At l equal to k,
the source is later exactly when $`i+k\lt d`$, giving $`d-k`$
choices for each $`1\le k\le r`$ and the remaining sum.

Use the retained literal 7-T, 12-elementary-Clifford, at-most-five-T-layer
[Toffoli word](RETURNED_DIAGONAL_ECHO_DEGREE.md#sharp-cubic-example-and-a-literal-native-ledger).
Implement each CZ as H, CNOT from h, H. One controlled source has

```math
T_D=14N_d,\qquad G_D=24N_d+3d,\qquad d_D\le10N_d.         \qquad(9)
```

Both F and its actual inverse are included. The source's actual inverse
has equal counts in this H/CNOT/T alphabet.

## 4. Paid read-only-address transport and exact ledger

Left multiplication by $`1+N_x`$ has coefficient formula

```math
z_{i,1}\longmapsto z_{i,1}+x_i,\qquad
z_{i,k}\longmapsto z_{i,k}+x_i z_{i+1,k-1}\quad(2\le k\le d).
                                                               \qquad(10)
```

Visit lengths from d down to two and execute the product toggles, then
apply the d edge CNOTs. Every source coefficient is still its old value,
so the word is precisely $`P_x`$. A second forward call is generally
not its inverse; the inverse is the actual reversed word.

There are $`R_d=d(d-1)`$ product toggles. For each, the
[R91 private-dirty phase echo](PARALLEL_BILINEAR_WORD_QUERY.md#1-a-phase-derivative-on-an-arbitrary-dirty-wire),
conjugated by two target Hadamards, supplies an exact Toffoli with
eight T/T-dagger gates, 18 CNOTs, two Hadamards, and eight T layers.
Its separate arbitrary helper returns after the occurrence. Reuse one
such helper sequentially. Every original address use is a CNOT control;
original address wires are never T targets or CNOT targets. Thus

```math
T_P=8R_d,\qquad G_P=20R_d+d,\qquad d_P\le8R_d.             \qquad(11)
```

The actual inverse has equal counts, since the only affine loads here
are CNOTs. Combining (6), (9), and (11) gives this exact sequential ledger:

| Resource | Lexicographic sequential native word |
|---|---:|
| T/T-dagger count | $`28N_d+16R_d`$ |
| Elementary Clifford count | $`48N_d+40R_d+8d+2`$ |
| T-depth | at most $`20N_d+16R_d`$ |
| Extra arbitrary work | $`d^2+1`$ bits |
| Clean work | 0 |
| Error / leakage | 0 / 0 |
| Returned work | Exact identity on every coefficient and helper input |

The original address and arbitrary output are excluded from extra work.
Both output Hadamards, all sources and controls, and every actual inverse
are charged. These are upper bounds for a specified word, not an optimized
conjunction compiler.

If replacement native controlled sources and transports have full operator
errors $`\delta_D`$ and $`\delta_P`$, including identity return of
all extra arbitrary helpers, their actual inverses have the same errors.
Unitary telescoping bounds the completed word's error by
$`2\delta_D+2\delta_P`$. Accuracy of an uncontrolled source does not
imply this controlled-source promise. If the completed word also
approximates a prescribed Boolean XOR query at error epsilon and
$`\varepsilon+2\delta_D+2\delta_P\lt2`$, that Boolean query must
equal the conjunction, since distinct Boolean XOR queries have full
operator distance two. This is a forward source-error ledger, not a
constraint on arbitrary fused implementations.

## 5. A separate parallel schedule with linear T-depth

An alternate pivot ordering improves the construction to linear depth.
Keep increasing length k, but group the pivots of each length into at
most three matching classes of the multigraph with edges
$`\{i,i+k\}`$. For $`k\lt d`$ its components are cycles, including
two parallel edges for a directed two-cycle, so three edge colors suffice.
Within a matching class, different pivots have disjoint start pairs.
All coefficient wires used by their families are therefore disjoint.
The k equal to d pivots have no updates.

At one pivot c, the updates have form $`\mathrm{CCX}(c,s_j;t_j)`$,
where the pairs $`(s_j,t_j)`$ are disjoint and exclude c. Conjugate by H
on every target, and combine the phase identities

```math
4cst\equiv c+s+t-(c\oplus s)-(c\oplus t)-(s\oplus t)
                       +(c\oplus s\oplus t)\pmod8.       \qquad(12)
```

For m pairs, merge the m common-control terms into $`T_c^m`$.
Modulo eight, this is a bounded Clifford power and at most one T gate.
Four physical T layers suffice:

1. The individual s and t phases and the merged c phase.
2. The negative $`c\oplus s`$ and $`c\oplus t`$ phases, on distinct
   s and t wires.
3. The negative $`s\oplus t`$ phases on t.
4. The positive $`c\oplus s\oplus t`$ phases on t.

Each parity has its explicit CNOT load and unload. At each layer all T
targets are distinct; shared-control CNOTs are serialized at the layer
boundary and counted. This uses no initialized helper. Matching classes
run their pivot families together, giving at most $`12(d-1)`$ T layers
for F, and at most $`24(d-1)`$ for the controlled source.

The new within-length order changes the same-length contribution to
$`N_d`$. Its total is still $`O(d^3)`$: the $`l\gt k`$ term is
unchanged, and all same-length contributions are at most
$`d\lfloor d/2\rfloor`$. Use this new F and its actual inverse
consistently; do not combine its schedule with the exact lexicographic
count in (8). Merged Clifford powers and their literal inverses add
only $`O(d^2)`$ elementary Clifford gates.

For transport, reserve d separate arbitrary helpers. At each descending
length, process all d product toggles simultaneously in eight T layers
using the private-dirty echo. Their source and target coefficient wires
are distinct at that length, and every original address interaction is
still a Clifford control. The helpers return after each batch. Transport
therefore has depth at most $`8(d-1)`$.

The separate parallel construction has the simultaneous ledger

```math
T=O(d^3),\qquad G=O(d^3),\qquad D_T\le64(d-1),
\qquad B=d^2+d,\qquad\text{clean work}=0.                 \qquad(13)
```

Its error and leakage are zero, and all coefficient and private work
return exactly. The control, source, inverse, and reference guarantees
are those of the sequential construction. This alternate schedule gives
no sublogarithmic depth result.

## 5a. The cycle source is a unitriangular inverse word query

The full source has a stronger interpretation than its conjunction
restriction. Cut the directed cycle between vertices $`d-1`$ and zero,
and divide its coefficients into a strictly upper-triangular matrix U
and a lower-triangular matrix L:

```math
U_{ij}=z_{i,j-i}\quad(i\lt j),\qquad
L_{ij}=z_{i,d+j-i}\quad(i\ge j).
```

U contains the nonwrapping paths; L contains the wrapping paths,
including the d full cycles on its diagonal. Together they contain
exactly the original $`d^2`$ coefficient bits. Put

```math
R(U)=(I+U)^{-1}=I+U+\cdots+U^{d-1},\qquad
q=\frac{d(d-1)}2,\qquad r=\frac{d(d+1)}2.
```

The character is the inverse pairing

```math
\chi(1+z)=\mathrm{tr}(R(U)L)
 =\sum_iL_{ii}+\sum_{i\lt j}L_{ji}R(U)_{ij}.              \qquad(13a)
```

To prove it, one may distinguish the wrap edge in the derivation (2).
This leaves chi unchanged: every distinguished edge gives a character
that is zero on $`1+e_{i,k}`$ for $`k\lt d`$ and one for k equal to d;
the ordered factorization shows these basis units generate the group.
In the expansion of $`\delta(u)u^{-1}`$, the differentiated path now
wraps once. A contributing length-d closed path can have no further
wrap, so all subsequent factors follow nonwrapping paths. A wrapping
path from j to i therefore contributes $`L_{ji}R(U)_{ij}`$, including
the length-d path when i equals j. Summing those terms proves (13a).

Let $`\mathcal R(U)`$ be the r-bit lower-triangular word whose
$`(j,i)`$ entry is $`R(U)_{ij}`$ for $`i\le j`$. Its diagonal
entries are the known bit one. Equation (13a) makes the source

```math
D_\chi|U,L\rangle=(-1)^{\langle\mathcal R(U),L\rangle}|U,L\rangle.
```

With Hadamards on all r lower-coordinate wires, the exact identity is

```math
H_LD_\chi H_L:
  |U,L\rangle\longmapsto|U,L\oplus\mathcal R(U)\rangle.  \qquad(13b)
```

Thus the uncontrolled cycle source and the complete arbitrary-output
unitriangular inverse word query are Clifford equivalent. In either
direction, T-count, T-depth, extra-helper width and full operator error
are unchanged; elementary Clifford count increases by exactly $`2r`$
Hadamards. All output values and arbitrary helper inputs are retained,
so the statement includes arbitrary external references. This is not
an initialized matrix-inversion oracle or a matrix-element encoding.

The matrix input U is preserved by the **completed** source or query.
The equivalence preserves native read-only use of U only if the supplied
implementation already restricts it to Clifford controls throughout.
The factorization F in Section 3 does not establish that stronger
property for these coefficient wires. Calling U an address does not
make those native operations read-only. The original conjunction word
(6) separately retains its proved read-only property for x.

### A single inverse entry needs no controlled source

For $`i\le j`$, let $`C_{h,ji}`$ be the CNOT from a separate
arbitrary bit h to $`L_{ji}`$. Execute chronologically

```math
D_\chi,\quad C_{h,ji},\quad D_\chi^\dagger,\quad C_{h,ji}.
```

The phase difference is $`hR(U)_{ij}`$ and L returns exactly.
Conjugating this word by $`H_h`$ gives the complete XOR query
$`h\mapsto h\oplus R(U)_{ij}`$. Suppose a supplied zero-clean
native source has T-count $`T_D`$, T-depth $`d_D`$, elementary
Clifford counts $`G_D,G_{D^\dagger}`$ for its forward word and
actual inverse, and $`W_D`$ extra arbitrary helpers. This entry
query has

```math
T\le2T_D,\qquad G\le G_D+G_{D^\dagger}+4,
\qquad D_T\le2d_D,\qquad B=r+W_D.                       \qquad(13c)
```

B excludes the q matrix-input bits and the one output h. It includes
the whole arbitrary L register. The only source control is the explicit
CNOT translation; no controlled source or prepared one-hot L is assumed.
Full source error $`\delta_D`$, including helper identity action,
gives completed query error at most $`2\delta_D`$. Native read-only
use of U again depends on the supplied source.

### A fully charged controlled source from an uncontrolled source

A different wrapper realizes the **whole controlled source**
$`(-1)^{h\chi(U,L)}`$ while preserving both U and L. Reserve a fresh
arbitrary r-bit lower array E, separate from U, L and h. The same source
on U,E has phase $`\langle\mathcal R(U),E\rangle`$. Let

```math
S_h:(h,L,E)\longmapsto(h,L,E\oplus hL),
```

where hL means componentwise multiplication by h. Execute chronologically

```math
D_\chi(U,E),\quad S_h,\quad D_\chi(U,E)^\dagger,
\quad S_h^\dagger.                                     \qquad(13d)
```

The difference between the two source phases is
$`h\langle\mathcal R(U),L\rangle`$, and E returns. This is exact
identity on the h-zero sector for arbitrary U,L,E and helpers. It is
not an inference that an uncontrolled source can be controlled for free.

For each of the r products $`E_k\mapsto E_k+hL_k`$, use the
private-dirty echo from Section 4, with its own arbitrary helper.
All r products run in eight T layers; the helpers are distinct and
every use of h and L is a Clifford control. One shear has exactly
$`8r`$ T gates and $`20r`$ elementary Cliffords, and its actual
inverse has the same counts. For the supplied source ledger above,
the controlled wrapper has

```math
\begin{aligned}
T&\le2T_D+16r,\\
G&\le G_D+G_{D^\dagger}+40r,\\
D_T&\le2d_D+16,\\
B_{\rm extra}&=r+\max\{W_D,r\}.
\end{aligned}                                          \qquad(13e)
```

Here extra work excludes the caller's q-bit U, r-bit L and control h.
It includes E and a separate pool large enough for either the source's
$`W_D`$ helpers or the shear's r private helpers. That pool is reused
only after each completed source or shear has returned it. E cannot
be merged with the pool: it carries its live arbitrary value through
both source positions. There are zero clean helpers. The literal
inverse Clifford count remains separate because a new source may use
native S gates whose inverse has a different elementary count.

If the supplied uncontrolled source approximates
$`D_\chi\otimes I_{W_D}`$ on its complete physical space within
$`\delta_D`$, the actual-inverse wrapper error is at most
$`2\delta_D`$. The shears here are exact. This provides a controlled
source accuracy guarantee by an explicit circuit; it does not assume
one from the uncontrolled guarantee alone. All identities hold with
arbitrary reference correlations. For an approximate source, helper
return is correspondingly approximate; reusing the pool remains valid
because its source guarantee is on the full physical space, including
any residual helper correlations. Hadamards on L also convert (13d)
to the controlled inverse-word query, adding $`2r`$ Clifford gates.

### What this duality does not supply

Neither (13b) nor (13e) alone reduces the source depth.
It identifies the exact owed operation: a faster source is a faster
complete unitriangular inverse-word query in another basis. Ordinary
clean-scratch matrix algorithms do not automatically realize it.

For example, block inversion gives

```math
\begin{pmatrix}A&B\\0&C\end{pmatrix}^{-1}
 =\begin{pmatrix}A^{-1}&A^{-1}BC^{-1}\\0&C^{-1}\end{pmatrix}.
```

After inverting the diagonal blocks in place, one still owes the
reversible action $`B\mapsto A^{-1}BC^{-1}`$ on its arbitrary existing
block. A product-accumulation shear into a second dirty array is not
that in-place action. Alternatively, inverse blocks accumulated into
dirty arrays retain unknown offsets; removing them by the ordinary
four-corner product echo repeats the recursive queries. Those identities
alone give no logarithmic-depth recurrence. This is the same distinction
between a cheap doubled shear and a paid in-place map exposed by the
[Gaussian inverse analysis](GAUSSIAN_INVERSE_EXTRACTION.md#2-the-cheap-gaussian-shear-does-not-supply-the-inverse-map).
The source equivalence alone supplies neither a faster complete circuit
nor an unrestricted depth lower bound. Section 5d supplies the additional
native construction by retaining helper counter-actions until the final
readout and actual uncompute.

## 5b. Exact inversion by commuting quadratic pivots

There is a full-input in-place realization of the inverse needed in
Section 5a whose completed quadratic pivots commute. This gives an exact
path identity and a cubic-count, linear-depth native source. It does not
reduce the established asymptotic depth or establish native read-only U.

Use the strictly upper triangular U of Section 5a, with indices from zero
through d minus one. For each internal vertex $`1\le k\le d-2`$, define

```math
[T_k(U)]_{ij}=
\begin{cases}
U_{ij}+U_{ik}U_{kj},&i\lt k\lt j,\\
U_{ij},&\text{otherwise}.
\end{cases}                                               \qquad(13f)
```

The source coordinates lie in column k above its diagonal and row k to
its right. Neither set meets the target rectangle $`i\lt k\lt j`$.
Thus all updates in one pivot use unchanged sources, commute as Boolean
Toffolis, and make $`T_k`$ an involution on every upper matrix word.

### Path formula and commutation

For a set S of internal vertices, apply each pivot in S once in any order.
The resulting coordinate is

```math
[T_S(U)]_{ij}
 =\sum_{\substack{i=v_0\lt v_1\lt\cdots\lt v_m=j\\
                 \{v_1,\ldots,v_{m-1}\}\subseteq S}}
       \prod_{a=0}^{m-1}U_{v_a,v_{a+1}},                   \qquad(13g)
```

where every edge coefficient on the right is an original input bit.
The direct-edge path is included, and all arithmetic is over
$`\mathbb F_2`$. For empty S this is U. Suppose the formula holds for S
and add a pivot k outside S. Existing paths avoid k. The product of the
current entries at $`(i,k)`$ and $`(k,j)`$ enumerates exactly the paths
through k whose other internal vertices are in S: their prefix vertices
are strictly below k and their suffix vertices strictly above k. The two
pieces have no common internal vertex, and every path through k has one
such decomposition. The added product therefore supplies precisely the
missing paths, proving the formula by induction for every pivot order.

In particular, the completed maps $`T_k`$ commute pairwise. Together
with involutivity, they satisfy $`T_ST_{S'}=T_{S\mathbin\triangle S'}`$.
For all internal vertices the path sum is the finite geometric series:

```math
\mathrm{Inv}(U)
 :=\left(\prod_{k=1}^{d-2}T_k\right)(U)
 =U+U^2+\cdots+U^{d-1}=(I+U)^{-1}-I.                  \qquad(13h)
```

For d equal to one or two the product is empty and the same formula
holds. In particular $`\mathrm{Inv}^2=I`$ as a permutation of
upper matrix words; this is distinct from an inverse word query into a
separate arbitrary output array.

### A literal zero-clean native schedule

Pivot k has $`m_k=k(d-1-k)`$ product updates. Give every update its own
arbitrary helper and use the private-dirty echo from Section 4. Its
relative-phase subwords act on the corresponding row coefficient,
target coefficient and private helper; its translation is the CNOT
from the column coefficient to that helper. Across a pivot the source
sets and target set are disjoint, the targets are distinct, and all
phase targets are private helpers. Each of the four relative-phase
bands therefore groups into one physical T layer, including when row
or column sources are shared. The two relative-phase words and the two
translations implement all $`m_k`$ products in eight T layers, with
$`8m_k`$ T gates and $`20m_k`$ elementary Cliffords. Every helper returns
exactly for every input; no initialized copy of a row, column or product
is used.

Execute these pivots sequentially, in any fixed order, reusing the
helper pool only after each pivot returns it. Since
$`\sum_k m_k=\binom d3`$ and
$`\max_k m_k=\lfloor(d-1)^2/4\rfloor`$, the resulting literal native
word V has

```math
T_V=8\binom d3,\qquad G_V=20\binom d3,\qquad
D_T(V)\le8\max\{d-2,0\},\qquad
W_V=\left\lfloor\frac{(d-1)^2}{4}\right\rfloor.          \qquad(13i)
```

Empty sums and the helper maximum are zero when d is at most two.
Here $`G`$ uses the same elementary Clifford convention as Section 4,
and $`W_V`$ counts extra arbitrary helpers, excluding the
$`q=d(d-1)/2`$ matrix bits. The exact full-space action is

```math
V|U,e\rangle=|\mathrm{Inv}(U),e\rangle
\quad\text{for every }U\in\mathbb F_2^q,
\ e\in\mathbb F_2^{W_V}.                               \qquad(13j)
```

All phases are literal plus one. By linearity this includes arbitrary
superpositions and external references. Its actual gate-by-gate inverse
has the same T, G and depth ledger and implements the same ideal
involution. An implementation must still use that actual inverse when
a surrounding word calls for $`V^\dagger`$; commutation of the completed
pivots does not license dropping an inverse or a native relative phase.
The approximation error here is zero.

For completeness, this supplies a second uncontrolled cycle source.
Execute V on U, then CZ between $`U_{ij}`$ and $`L_{ji}`$ for all
$`i\lt j`$ and Z on each diagonal $`L_{ii}`$, then the actual
$`V^\dagger`$. Equation (13h) and the inverse-pairing identity (13a)
give exactly $`(-1)^{\chi(U,L)}`$, returning U, L and all helpers.
With CZ implemented by H, CNOT, H and forward Z by $`S^2`$, its ledger is

```math
T_D=16\binom d3,\quad
G_D=40\binom d3+3q+2d,\quad
D_T(D)\le16\max\{d-2,0\},\quad W_D=W_V.                \qquad(13k)
```

The $`q+r=d^2`$ U,L wires coexist with the helper pool. In an elementary
H,S,CNOT convention with $`S^\dagger=S^3`$, the actual inverse has
$`G_{D^\dagger}=40\binom d3+3q+6d`$; its T-count and depth agree with
the forward word. This keeps literal inverse accounting distinct even
though the ideal phase source is self-adjoint. Applying the controlled
wrapper (13d) still owes all of (13e), including its fresh r-bit array
and the simultaneous allocation $`r+\max\{W_D,r\}`$; (13k) is not a
free controlled source.

### Why this does not give parallel inversion

Commutation here concerns completed pivots, not every Toffoli or native
subword across pivots. For d equal to four, put
```math
\begin{aligned}
a&=U_{01},& b&=U_{12},& c&=U_{23},\\
p&=U_{02},& q'&=U_{13},& s&=U_{03}.
\end{aligned}
```
The two pivots update

```math
T_1:\ (p,s)\mapsto(p+ab,s+aq'),\qquad
T_2:\ (q',s)\mapsto(q'+bc,s+pc).
```

Either order gives $`s\mapsto s+aq'+pc+abc`$. Individual updates to
p and updates controlled by p do not commute; the other cross-pair
produces the matching cubic term in the reverse pivot order. A single
simultaneous evaluation of every quadratic update on original U would
instead give $`U+U^2`$, missing $`U^3`$ on a three-edge chain. Thus the
path identity alone does not put overlapping native pivots in one T
band or give a balanced-block depth recurrence. The demonstrated bound
in (13i), and hence (13k), remains linear in d.

Moreover V changes matrix-coordinate wires by nonlinear updates.
The resulting source returns U only at completion; it does not keep U
restricted to native Clifford controls. The word-query equivalence
(13b) therefore inherits no strict read-only-U guarantee from this
implementation. The original read-only conjunction circuit and its
fair comparisons in Section 6 are unchanged. This identity provides
neither a faster complete reader nor an unrestricted native depth
lower bound.

## 5c. The isolated dirty in-place block boundary

A balanced block continuation exposes a precise helper-return failure.
The exact block formula and the scoped array-shear invariant below
identify it without giving a new native-depth upper or lower bound.
This is the historical boundary for a returned in-place subroutine;
Section 5d bypasses that interface for the completed query.

### What a contiguous block of commuting pivots does

Partition the vertices into consecutive sets L,K,R, and write the
strictly upper input as

```math
U=\begin{pmatrix}P&B&D\\0&A&C\\0&0&Q\end{pmatrix},
\qquad R_A=(I+A)^{-1}.
```

The product of exactly the pivots whose vertices lie in K has the
complete coordinate action

```math
(A,B,C,D)\longmapsto
(R_A-I,\ B R_A,\ R_A C,\ D+B R_A C),                 \qquad(13l)
```

and leaves P,Q unchanged. No further input coefficient is promised zero. To prove
this, use the existing partial-path formula: a path internal to K gives
the inverse series $`R_A`$; a path from L to K has one initial B edge; a
path from K to R has one terminal C edge; and a path from L to R either
uses its direct D edge or one B edge, a K-internal path, and one C edge.
No internal vertex in L or R is permitted, so P,Q play no role.

This is also a direct all-input involution. After one application the
new internal inverse series is $`R_A^{-1}=I+A`$. Applying (13l) twice
restores B,C and adds $`B R_A C`$ twice to D. Thus (13l) is a full permutation identity,
including arbitrary external references, rather than an initialized
Schur-complement computation.

The important missing operation is not the cross-block XOR alone.
Even after inverting A, (13l) still changes its arbitrary incident blocks
by $`B\mapsto BR_A`$ and $`C\mapsto R_AC`$. A product accumulator into a
second dirty array does not implement these in-place changes.

### The two-sided XOR has a complete native ledger

For unit upper matrices A of order a and C of order c, define the
invertible linear operator M on arbitrary a-by-c arrays by $`M(B)=ABC`$.
The XOR shear $`E\mapsto E+M(B)`$ has an exact zero-clean native implementation.
Their unit diagonals are fixed constants, not initialized physical
registers: only the strict upper coefficients are input wires. Reserve
an arbitrary a-by-c F, and execute chronologically

```math
P:F\mathrel{+}=AB,\quad
Q:E\mathrel{+}=FC,\quad P^\dagger,\quad Q^\dagger.  \qquad(13m)
```

The two Q uses add $`(F+AB)C`$ and FC, so E gains ABC and F returns. A,C,B
are preserved throughout, and A,C,B enter the product shears only as
native Clifford controls. All phases are plus one and every private
phase helper returns on every input. The last two words are the actual
native inverses, not independently chosen replacements.

Put

```math
\tau_A=c\binom a2,\qquad\tau_C=a\binom c2.
```

P contains ac Clifford identity loads and $`\tau_A`$ nonlinear products;
Q similarly has ac identity loads and $`\tau_C`$ products. The retained
private-helper phase schedule permits shared targets: conjugate the
whole commuting product family by one Hadamard per target, and use a
distinct arbitrary phase helper per nonlinear product. Each family has
eight T layers. Using a harmless upper bound of twenty elementary
Cliffords per product, (13m) has

```math
T_S=16(\tau_A+\tau_C),\quad
G_S\le40(\tau_A+\tau_C)+4ac,\quad D_S\le32,\quad
W_S=ac+\max(\tau_A,\tau_C).                         \qquad(13n)
```

$`W_S`$ excludes A,C,B,E and includes F plus the largest private phase
pool. That pool is reused only after each complete shear returns it;
F cannot be merged with it. The bounds remain valid at small dimensions
where a nonlinear family is empty. This is a complete out-of-place XOR
interface, not an in-place multiplication result.

### Three shears leave a counter-action on arbitrary work

Suppose both $`S_M`$ and $`S_{M^{-1}}`$ are available with all their own helpers
returned. On two arbitrary arrays B,E, execute

```math
E\mathrel{+}=MB,\quad
B\mathrel{+}=M^{-1}E,\quad
E\mathrel{+}=MB,\quad \mathrm{SWAP}(B,E).       \qquad(13o)
```

The first three operations map $`(B,E)`$ to $`(M^{-1}E,MB)`$; hence the exact
completed result is

```math
(B,E)\longmapsto(MB,M^{-1}E).                         \qquad(13p)
```

No phase is missing, but the arbitrary helper has not returned. Setting
E to zero would hide this failure; arbitrary E has no such promise.
Returning E now requires the in-place map $`E\mapsto ME`$, the same type of
operation that the proposed conversion was supposed to supply.

For balanced inversion, the live diagonal blocks after their recursive
inversions are A,C, and the desired off-diagonal map is precisely ABC.
The middle shear in (13o) requires the original $`A^{-1},C^{-1}`$. They are
not free second copies. If native in-place diagonal inversion words
$`V_A,V_C`$ are supplied, one honest implementation of that middle shear
uses their actual inverses in parallel, then the polynomial XOR shear
with the restored diagonal coefficients, then $`V_A,V_C`$ in parallel.

Write the simultaneous diagonal costs as
$`T_I=T_A+T_C`$, $`D_I=\max(D_A,D_C)`$, and $`W_I=W_A+W_C`$.
Keeping the actual forward/inverse Clifford counts separate, (13o) therefore obeys

```math
\begin{aligned}
T&\le3T_S+2T_I,\\
G&\le3G_S+G_A+G_{A^\dagger}+G_C+G_{C^\dagger}+3ac,\\
D_T&\le3D_S+2D_I,\\
W_{\rm extra}&=ac+\max(W_S,W_I).
\end{aligned}                                         \qquad(13q)
```

The ac term outside the maximum is E, which stays live; B is the caller's
output block. The two diagonal work pools coexist when those inversions
run in parallel. The returned S work and diagonal work may otherwise
share a pool. The first and last shears use the live A,C; only the middle
one changes them internally. Thus (13q) is an explicitly paid circuit
for (13p), not a ledger for an imaginary returned-helper primitive.
The diagonal inversion calls also prevent this implementation from
establishing native read-only use of the matrix-input coefficients.

Including the initial diagonal inversions, this particular black-box
schedule has three sequential diagonal-inversion positions along its
critical path. It therefore does not even justify the proposed
$`D(d)=D(d/2)+O(1)`$ recurrence: its call accounting would read
$`3D(d/2)+O(1)`$, before the still-missing helper return. This is not a
complete inversion upper bound, nor a lower bound on other schedules.

If all building blocks are exact, all these statements are exact. More
generally, full-space errors add across the three shear calls and the
two diagonal-inversion positions. A valid bound, including the actual
inverse errors, is $`3\varepsilon_S+2(\varepsilon_A+\varepsilon_C)`$. It
compares to (13p), not to the unimplemented $`(MB,E)`$ map.

### A scoped invariant for the array-only repair

There is a precise reason that adding more of the same kind of
whole-array shear need not repair (13p). Fix a nonidentity unipotent
binary operator M with nilpotence index $`\nu`$ for $`M-I`$. Let

```math
\mathcal R=\mathbb F_2[t]/(t^\nu),\qquad t\mapsto M-I.
```

This is a faithful representation: the minimal polynomial of $`M-I`$ is
$`t^\nu`$, so no nonzero residue of degree below $`\nu`$ acts as zero. Grant any
number s of arbitrary equal-size arrays and any operations

```math
E_i\mathrel{+}=p(M)E_j\quad(i\ne j),
```

for arbitrary polynomials p, together with whole-array swaps and
whole-array CNOTs. The coefficients are frozen. Each operation is an
s-by-s matrix over the commutative ring $`\mathcal R`$ and has determinant one.
For swaps, the usual determinant sign is also one in characteristic
two. Thus every finite word in this move set has determinant one.

The sought full-return action $`(MB,E_2,\ldots,E_s)`$ would instead have array
matrix $`\mathrm{diag}(1+t,1,\ldots,1)`$, whose determinant is
$`1+t\ne1`$. Faithfulness means that equality on all physical input columns would imply equality
of these ring matrices, a contradiction. In comparison, the paired
action $`\mathrm{diag}(M,M^{-1})`$ in (13p) has determinant one, as it must.

This excludes a returned-helper conversion using only these uniform
polynomial array shears, even if every polynomial-coefficient shear is
granted at no cost. It is a limitation of that proposed interface, not
a gate-count or T-depth lower bound. Within-array coordinate changes,
coefficient mutations outside the fixed coefficient algebra, nonuniform
coordinate operations, and Hadamards are all outside this move set.
The existing native pivot word uses such coordinate operations and is
not obstructed. Keeping helper counter-actions across a larger algorithm
with different coefficient operators is also not ruled out here.

This attempted balanced recurrence therefore stops at (13p). Section 5d
instead keeps the helper action live until after a readout, with one
diagonal-inversion position per level and a final actual uncompute. That
completed-query construction leaves the isolated returned-helper move
set of this invariant. The linear-depth, smaller-work construction in
Section 5b and both conjunction comparisons in Section 6 remain valid.
Even a faster structured inverse-word query needs a separate arbitrary-table
encoding and the complete R98 ledger before supplying a general reader.

## 5d. A logarithmic-depth complete inverse-word query

The two-shear block step below computes the desired inverse coordinates
independently of its arbitrary helpers, while leaving an invertible
counter-action on those helpers. Keeping that counter-action until an
actual compute/readout/uncompute completes gives a valid full-return
inverse-word query. It avoids the failed three-child-call recurrence in
Section 5c. It does not give returned-work in-place inversion by itself,
or native read-only use of the matrix coefficients.
Here $`n=d\ge1`$ is the matrix order of the cycle representation.

### The two-shear block identity

For a unit upper input matrix, write

```math
X=\begin{pmatrix}A&B\\0&C\end{pmatrix},\qquad
N(E)=AEC.
```

The matrices A,C are unit upper of orders a,c; their unit diagonals are
constants, not initialized qubits. B and a separate helper E are arbitrary
a-by-c arrays. Execute, in chronological order:

1. $`B\mathrel{+}=AEC`$, using the original diagonal coefficients.
2. Recursively invert A and C in parallel, allowing private helper
   garbage but requiring their principal inverse coordinates to be
   independent of all those helper inputs.
3. $`E\mathrel{+}=A^{-1}BC^{-1}`$, using the now-inverted diagonal coefficients.
4. Swap B and E entrywise.

The first shear gives $`B_1=B_0+N(E_0)`$. The second gives

```math
E_1=E_0+N^{-1}(B_1)=N^{-1}(B_0).
```

Thus the completed principal block is the correct inverse and the
parent helper has the explicit action

```math
B_{\rm out}=A^{-1}B_0C^{-1},\qquad
E_{\rm out}=B_0+A E_0 C.                              \qquad(13s)
```

The child helpers undergo their own recursive counter-actions. Neither
child computation touches B or E. In particular, the second shear never
needs the original coefficient blocks again. There is only one recursive
inversion position, rather than an inversion, its reversal, and another
inversion.

### Exact all-input computation, with garbage kept live

Recursively split every vertex interval into balanced consecutive halves.
At every internal node v, reserve a distinct $`a_v\times c_v`$ arbitrary
array $`E_v`$. The total number of these permanent block-helper bits is

```math
\sum_v a_vc_v=\binom n2=:q.
```

Every original coefficient appears in exactly one off-diagonal block $`B_v`$.
For fixed original X, the helper action at that node is
$`E_v\mapsto B_v+A_vE_vC_v`$, with the original diagonal submatrices
$`A_v,C_v`$.
Each is an invertible affine map. Temporary native work is returned by
every completed product shear. Starting with $`V_1=I`$, the block identity
and induction therefore give the exact native word

```math
\begin{aligned}
V_n|U,e\rangle&=|\mathrm{Inv}(U),g_U(e)\rangle,\\
\mathrm{Inv}(U)&=(I+U)^{-1}-I.
\end{aligned}                                        \qquad(13t)
```

on every principal strict upper word U and every helper input e, with
literal phase plus one. Here $`g_U`$ is a permutation on the entire helper
fiber; permanent helper components have the node formula above, and
unused temporary components return. No helper is assumed initialized.
This is an all-column unitary statement, hence includes arbitrary
external references. The helpers in (13t) deliberately have not returned.

The actual inverse of $`V_n`$ reverses the native gate order and adjoints
every gate. It restores exactly the input U,e from each output in (13t).
An independently chosen inversion routine is not substituted for that inverse.

### Native resources and simultaneous work

Use the paid two-sided XOR shear from Section 5c for each of the two
parent shears. For a node with sizes a,c, put

```math
\tau=c\binom a2+a\binom c2
     =\binom{a+c}3-\binom a3-\binom c3.
```

One shear has $`T=16\tau`$, $`G\le40\tau+4ac`$, T-depth at most 32, and
extra arbitrary work

```math
W_S=ac+\max\{c\binom a2,a\binom c2\}.
```

This work includes its arbitrary product accumulator and private phase
helpers. Shared product targets use a single enclosing Hadamard per
target, with distinct phase helpers for the constituent products. All
source matrices enter that shear only as native Clifford controls.
The diagonal blocks are changed by their recursive inversions between
the two shears; the complete $`V_n`$ is not native read-only in U.

The two shears plus ac swaps therefore cost $`32\tau`$ T gates and at most
$`80\tau+11ac`$ elementary Cliffords. Summing over the recursion tree gives

```math
T_V=32\binom n3,\qquad
G_V\le80\binom n3+11\binom n2.                       \qquad(13u)
```

The actual inverse has the same counts: $`V_n`$ uses only H,CNOT,T and
$`T^\dagger`$, so every Clifford in the native word is self-adjoint.

The two children use disjoint principal blocks and private work pools,
so their T bands run in parallel. Parent shears precede and follow both
children. Nodes of size at most two have no T gates. Hence, with

```math
h_n=\max\{\lceil\log_2 n\rceil-1,0\},
\qquad D_T(V_n)\le64h_n.                             \qquad(13v)
```

For the simultaneous work allocation, set $`W(1)=0`$ and, for the balanced
split $`a=\lfloor n/2\rfloor`$, $`c=\lceil n/2\rceil`$, reserve the parent E
separately from a pool used in three successive stages: the first shear,
the parallel children, and the second shear. This gives the exact reservation

```math
W(n)=ac+\max\{W(a)+W(c),W_S\}.                       \qquad(13w)
```

The pool may contain child garbage at the second shear. This is valid:
the shear's arbitrary-work contract returns that current state exactly,
including its correlations, and $`V_n`$'s later actual inverse will undo the
child garbage. B and the parent E remain outside that pool throughout.

For the balanced split, (13w) has the closed form

```math
W(n)=\frac{ac(c+3)}2.                                \qquad(13x)
```

One way to check local dominance is to inductively use
$`W(m)\le m^2(m+7)/16`$. For $`c\in\{a,a+1\}`$ and $`a\ge1`$,

```math
\frac{a^2(a+7)+c^2(c+7)}{16}
 \le\frac{ac(c+1)}2=W_S.
```

The claimed closed form also satisfies that induction bound by the two
parity cases for m. Thus the shear pool dominates the simultaneous child
pools, proving (13x) for $`n\ge2`$. The base case $`n=1`$ is separate.
At $`n=2`$, the principal inverse is already the input, but this literal schedule still runs its
Clifford-only two-shear recursion and reserves $`W(2)=2`$. It deliberately
keeps that helper counter-action, rather than replacing $`V_2`$ by identity.
In particular, $`W(n)\ge q`$: all permanent node arrays coexist, but may
also serve as returned scratch at another recursion stage. There are zero
clean bits.

### A complete inverse-word query

Let L be an arbitrary lower word of size $`r=n(n+1)/2`$. After $`V_n`$, copy each
strict upper principal output $`U'_{ij}`$ by CNOT into $`L_{ji}`$, and apply X
to each diagonal $`L_{ii}`$. Then execute the actual $`V_n`$ inverse. The middle
word changes only L, so every principal and helper coordinate required
by that inverse is still exactly the output from (13t). The result is

```math
\begin{aligned}
|U,L,e\rangle&\longmapsto|U,L\mathbin\oplus\mathcal R(U),e\rangle,\\
\mathcal R(U)_{ji}&=[(I+U)^{-1}]_{ij}.
\end{aligned}                                        \qquad(13y)
```

This is a literal full-input query with every helper returned. Its
resources, excluding the q matrix inputs and r arbitrary output bits,
are

```math
\begin{aligned}
T_Q&=64\binom n3,\\
G_Q&\le160\binom n3+23q+4n,\\
D_T(Q)&\le128h_n,\\
W_Q&=W(n).
\end{aligned}                                        \qquad(13z)
```

Here X is expanded as $`HS^2H`$. The actual reversed query has the same
T-count and depth; if every $`S^\dagger`$ is expanded as $`S^3`$, its displayed
elementary Clifford bound is $`160\binom n3+23q+8n`$. These are literal-word
counts; no free inverse normalization is assumed.

There is also a direct controlled query. At the middle readout, replace
the q CNOTs by Toffolis controlled on h and $`U'_{ij}`$, and the n diagonal
X gates by CNOTs from h. Use q existing work bits as the private phase
helpers for those Toffolis. $`W(n)\ge q`$, and their current garbage values
are permitted by the exact private-helper echo. Each returns its current
value before $`V_n^\dagger`$ begins. All q products run in eight T layers.
Thus the controlled full-word query has

```math
\begin{aligned}
T_{Q_h}&\le64\binom n3+8q,\\
G_{Q_h}&\le160\binom n3+42q+n,\\
D_T(Q_h)&\le128h_n+8,\\
W_{Q_h}&=W(n).
\end{aligned}                                        \qquad(13aa)
```

Its actual inverse has the same displayed ledger. The control h is
used only as a native Clifford control. On the h-zero sector the entire
completed operation is identity, including all helper garbage. The
matrix inputs U return only at completion and are not native read-only.

### Cycle source and precise scope

Alternatively put the inverse-pairing phase between $`V_n`$ and its actual
inverse: $`\mathrm{CZ}(U'_{ij},L_{ji})`$ for $`i\lt j`$ and $`Z(L_{ii})`$
on the diagonal. This gives the full cycle source
$`D_\chi|U,L,e\rangle=(-1)^{\chi(U,L)}|U,L,e\rangle`$, with

```math
\begin{aligned}
T_{D_\chi}&=64\binom n3,\\
G_{D_\chi}&\le160\binom n3+25q+2n,\\
D_T(D_\chi)&\le128h_n,\qquad W_{D_\chi}=W(n).
\end{aligned}                                        \qquad(13ab)
```

Each CZ costs three elementary Cliffords, and $`Z=S^2`$ costs two.
The actual inverse changes the last Clifford term from $`2n`$ to $`6n`$
under the same H,S,CNOT convention. This agrees with the existing
phase/query duality, though the direct readout (13y) cancels unnecessary
Hadamards.

The controlled source is equally direct. Between the same computation
and its actual inverse, use $`\mathrm{CCZ}(h,U'_{ij},L_{ji})`$ on each
strict upper coefficient and $`\mathrm{CZ}(h,L_{ii})`$ on each diagonal.
The former use the phase form of the same private-helper echo, with no
target Hadamards: eight T gates and eighteen CNOTs per coefficient.
Its q helpers are borrowed from the current $`V_n`$ work and restored
before uncomputation, just as in (13aa). Thus

```math
\begin{aligned}
D_{\chi,h}|h,U,L,e\rangle
 &=(-1)^{h\chi(U,L)}|h,U,L,e\rangle,\\
T_{D_{\chi,h}}&\le64\binom n3+8q,\\
G_{D_{\chi,h}}&\le160\binom n3+40q+3n,\\
D_T(D_{\chi,h})&\le128h_n+8,\qquad W_{D_{\chi,h}}=W(n).
\end{aligned}                                        \qquad(13ac)
```

The control h occurs only as a native Clifford control, including in
the diagonal CZ realization. All its work returns on both control
sectors. Its actual inverse has the same displayed ledger, since this
word contains no S gates. No fresh lower array or controlled invocation
of $`V_n`$ is needed. The unconditional middle readout in (13ab) has the
same full-return property; phase entanglement with L or an external
reference does not change the all-input inverse identity.

Every stated circuit is exact. If a substitute V approximates its entire
garbage-output permutation (13t) in full operator norm by $`\delta`$, using
its actual inverse around any of these exact readouts gives error at most
$`2\delta`$, including the work-return error and external references.
The readout remains a full-space unitary when it borrows the computation's
dirty work, so the same two-term operator-norm estimate applies. A guarantee
only on principal coordinates would not suffice. No approximate subroutine
is used in the exact construction above.

This is a logarithmic-depth structured inverse-word query and phase
source with cubic T/G count and cubic arbitrary-work width. It is not a
returned-work in-place inversion interface, a native-read-only-U query,
or an arbitrary R98 table reader. The unchanged cycle transport $`P_x`$
still has the previously retained linear-depth schedule, so improving
the source alone does not improve the full cycle conjunction circuit.
The same conjunction already has the logarithmic-depth read-only dirty
baseline with smaller polynomial count/work, and the ordinary completed
MCX baseline is smaller still. No improvement over those comparisons or
the general R91/R96 frontier is claimed.

## 6. Fair comparison and the remaining depth question

For the **same strictly read-only-address interface**, the retained
[dirty conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
already has, with $`\alpha=\log_2 3`$,

```math
T,w=O((d+1)^\alpha),\qquad
G=O((d+1)^\alpha\log(d+2)\log\log(d+4)),\qquad
D_T=O(\log(d+2)).                                       \qquad(14)
```

The cited chapter defines $`\chi(t)=\log_2(t+2)`$; its $`O(\chi(d))`$
depth is the last bound above. Its counter, router, and output echoes
return all dirty work. It is not an ordinary clean AND tree with its
internal nodes relabeled as dirty.

For the weaker requirement that address controls return at completion,
the imported
[two-dirty-helper MCX](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates)
has $`T,G=O(d)`$ and $`D_T=O(\log(d+1))`$. That paragraph does not
establish that the original controls occur only as native Clifford
controls throughout, so it must not be substituted for (14) when comparing
the stricter R98 interface. Both comparisons show that the cycle word is
not a better conjunction circuit.

For $`Q=2^d`$, the cycle word's polynomial work eventually fits below
square-root Q, but producing one high-degree function does not give a
reader for every Q-row table. Its source phase has a circular-composition
ANF: sum over all nonempty subsets of the d vertices, multiplying the
path coefficients on the forward cyclic gaps between consecutive chosen
vertices. To derive this, expand $`\delta(u)u^{-1}`$. A length-d
contribution partitions one cycle into positive paths; exactly one selected
path contains the distinguished edge, fixing the starting factor in delta.
Each nonempty vertex subset consequently appears once. There are
$`2^d-1`$ such monomials. The factorization is a polynomial-size way of
evaluating this source; expanding every monomial in parallel is not a
square-root-Q resource construction. Section 5d now supplies a full-input
native block construction for this structured source, with arbitrary
scratch and complete return at cubic width. It still supplies neither
an arbitrary-table encoding nor a faster complete conjunction transport.

There is a precise limitation on accelerating this transport by ordinary
reversible arithmetic. Its inverse coordinate is

```math
[P_x^{-1}(z)]_{i,d}
 =z_{i,d}
  +\sum_{s=1}^{d-1}\left(\prod_{j=0}^{s-1}x_{i+j}\right)z_{i+s,d-s}
  +\prod_{j=0}^{d-1}x_{i+j}.                              \qquad(15)
```

This polynomial has degree exactly d. A classical reversible word with
arbitrary affine intervals and k disjoint Toffoli layers has coordinate
degree at most $`2^k`$. Its actual reverse has the same number of layers,
so every such realization of $`P_x`$ needs
$`k\ge\lceil\log_2 d\rceil`$. A quadratic forward coordinate map
does not eliminate its inverse cost. Likewise F has a coordinate of
degree at least d: the XOR of its top factor bits restricts to the
degree-d conjunction when nonedge coefficients are zero. This forces
the same lower bound for a factorization implementation in that model.

These are **Toffoli-layer** bounds, not unrestricted native T-depth
bounds. Nonmonomial intermediate Clifford operations can leave the
coordinate-degree argument's scope. The constructive resumption question
is whether an interference-based controlled source and transport can
beat the retained logarithmic additive depth while returning every
arbitrary wire, or whether a different noncommutative representation can
carry arbitrary tables with the full R98 ledger. This chapter establishes
neither such mechanism nor a general obstruction to one.

## 7. Two further representation boundaries

For any representation rho and fixed unitary D, suppose
$`D\rho(g)D^\dagger\rho(g)^\dagger`$ is scalar for each g in a
generating set. Then $`D\rho(g)D^\dagger=\lambda(g)\rho(g)`$;
multiplying these relations and their inverses extends lambda to a
character of the generated subgroup H. Every relation gives scalar one,
so this extension is well-defined and kills $`[H,H]`$. This applies
even when D is nonmonomial. It does not bound the Boolean degree of the
coordinate encoding $`x\mapsto g_x`$; (4) exhibits that distinction.

An apparent alternative uses dual numbers over $`\mathbb F_{2^m}`$,
with unit character
$`\chi(s+t\epsilon)=\mathrm{Tr}(ct/s)`$ for nonzero s.
Variable division suggests a nonlinear small-coordinate source, but the
full unit group has size $`(2^m-1)2^m`$, not a power of two when
$`m\gt1`$. Regular unit coordinates therefore do not fill arbitrary
physical bit strings. Padding by the nonunits and using multiplication
on the whole ring fixes the zero basis vector, where a diagonal-source
commutator has phase one. Even a nonmonomial source cannot supply a
nontrivial Boolean scalar for that natural permutation: a permutation
with a fixed point has positive trace, whereas
$`DPD^\dagger=-P`$ would require trace zero. This excludes that
padded multiplication representation, not all representations of the
group. Principal units $`1+t\epsilon`$ instead have an additive group
law and linear coefficient characters.

## 8. Focused evidence

The [cycle-character checks](../../tests/test_noncommutative_cycle_character.py)
retain the noncommutative multiplication, derivation and cyclic trace,
character identity, even- and odd-length conjunctions, factorization
order and actual inverse, all-column commutator return, and one complete
native arbitrary-output query with its literal resource counts. Separate
checks cover the shared-control four-layer phase word and colored pivot
ordering, the inverse-pairing identity and Fourier word-query equivalence,
and the complete dirty controlled-source wrapper. The commuting-pivot
checks exhaust small partial path sets and group actions, compare the full
inversion with independent back substitution, and detect the missing
cubic path in simultaneous quadratic updates. They also check every
physical column of the eight-wire native inversion, both pivot orders
and their actual inverses, and the compute-pair-uncompute source phase.
The block checks retain the contiguous partial-closure formula, the
actual helper counter-action, a literal native two-sided XOR shear, and
the full finite group of two-array shears over the dual-number ring.
The recursive checks exhaust the native two- and three-vertex computation
on all helper columns and its actual inverse, check a shared product target,
and verify the parallel native count/depth schedule through nine vertices.
They also check independent back substitution and the explicit parent
helper action with live nonlinear child garbage, all columns of the direct
and controlled small word query, and the cycle-source readout pairing. The
[commutative tests](../../tests/test_nonlinear_character_transport.py)
also check affine-cube products in a genuinely multivariate nilpotent
algebra. These are bounded evidence; the proofs above and in the
commutative theorem establish the arbitrary-size statements. No finite
fixture is an asymptotic native compiler or an unrestricted lower bound.
