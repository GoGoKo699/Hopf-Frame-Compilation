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

This is a representation escape, not a query improvement. The demonstrated
output is a single conjunction. The explicit parallel native word has
$`T,G=O(d^3)`$, extra dirty width $`d^2+d`$, and
$`D_T\le64(d-1)`$. The retained strictly read-only conjunction already
has logarithmic depth with smaller polynomial count and width; the ordinary
completed-circuit MCX baseline is smaller still. Section 6 gives their
distinct contracts. No arbitrary-table reader, sublogarithmic depth law,
or R91/R96/R98 improvement follows.

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
square-root-Q resource construction. A determinant or block-doubling
algorithm would still need a full-input native realization with arbitrary
scratch and its complete return proof.

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
ordering. The [commutative tests](../../tests/test_nonlinear_character_transport.py)
also check affine-cube products in a genuinely multivariate nilpotent
algebra. These are bounded evidence; the proofs above and in the
commutative theorem establish the arbitrary-size statements. No finite
fixture is an asymptotic native compiler or an unrestricted lower bound.
