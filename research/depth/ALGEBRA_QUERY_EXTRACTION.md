# Full-input algebra queries and fast returned multiplication

[Balanced HRT source](BALANCED_ALGEBRA_PERMUTATION.md) ·
[Conditional query interface R98](QUERY_DEPTH_FRAME_TRANSFER.md) ·
[Exact fixtures](../../tests/test_algebra_query_extraction.py)

An arbitrary Boolean table can be extracted from the HRT source without
initializing any coefficient bank. Two literal commutators remove its
unknown input coefficients. The complete query needs eight paid source
occurrences. A separate returned-work compilation of ranked subset
convolution improves the source count from exponential base three to base
two times polynomial factors, retaining logarithmic depth.

These are complete upper constructions, but neither improves R91 or R98.
For an a-bit table, the resulting coefficient carrier already occupies
$`2^{a+1}-1`$ arbitrary qubits, and the source count is too large. The
initialization gap is resolved; the count, width and depth requirements of
a faster general query remain unpaid.

## 1. An arbitrary-table scalar on the entire coefficient space

Use the square-zero algebra $`R_r`$ and HRT permutation $`F_r`$ from the
[balanced proof](BALANCED_ALGEBRA_PERMUTATION.md), with $`r=a+1`$ and
$`a\ge1`$. Its coefficient input y consists of all nonempty subsets of
$`[r]`$, every bit arbitrary. Write f for its forward coordinate map and k
for its actual inverse. The product identity is

```math
f(c+d)=f(c)+f(d)+f(c)f(d).
```

It follows directly by grouping the factors $`1+c_Se_S`$: each squares to
one, all commute, and $`e_S^2=0`$ for nonempty S.

The original address x occupies separate unchanged wires. Let
$`\iota(x)=\sum_{i=1}^a x_i e_i`$ and use the Clifford translation
$`X^{\iota(x)}`$, one CNOT from each address to its singleton coefficient.
The source conjugation

```math
J_x=F_rX^{\iota(x)}F_r^\dagger
```

has, on every coefficient input y, the literal permutation action

```math
J_x:y\longmapsto y+A_x+yA_x,\qquad
A_x=f(\iota(x))=\prod_{i=1}^a(1+x_i e_i)-1.       \qquad(1)
```

No coefficient of y is initialized. The element $`A_x`$ has no occurrence
of the dummy generator $`e_r`$. In particular its top coefficient is zero.
Let $`\lambda(y)=[e_{[r]}]y`$ and let Z be the physical Z on the top
coefficient qubit. Chronologically execute $`J_x,Z,J_x^\dagger,Z`$.
This operation $`C_x`$ is exactly

```math
C_x|y\rangle=(-1)^{\lambda(yA_x)}|y\rangle.       \qquad(2)
```

Indeed the two Z phases are $`\lambda(y+A_x+yA_x)+\lambda(y)`$.
The intervening actual inverse restores y. For each fixed address x, equation (2) is a full-space
Pauli on the coefficient register, with no omitted coefficient-dependent
scalar. It is not asserted to be a Pauli on the joint coherent address
and coefficient space.

Write the prescribed Boolean table in algebraic normal form,

```math
h(x)=\beta+\sum_{\varnothing\ne S\subseteq[a]}\alpha_S
                         \prod_{i\in S}x_i.
```

The fixed coefficient vector

```math
v=\sum_{\varnothing\ne S\subseteq[a]}\alpha_S e_{[r]\setminus S}
```

uses only nonempty basis elements; all contain the dummy generator.
Its support has size at most $`2^a-1`$. For a separate arbitrary output z,
let $`L_z:y\mapsto y+zv`$. This is a Clifford word of $`|v|`$ CNOTs,
with z only a control. The chronological word

```math
C_x,\ L_z,\ C_x^\dagger,\ L_z^\dagger
```

has the exact scalar action

```math
(-1)^{z\lambda(vA_x)}I_y
 =(-1)^{z(h(x)+\beta)}I_y.                         \qquad(3)
```

A basis term of v contributes to the top product precisely when its
complementary term occurs in $`A_x`$, proving the last equality.
Hadamards on z at the two ends, followed by X on z if beta is one, give
$`z\mapsto z+h(x)`$. Every coefficient and source helper returns, and
all occupied phases are literally plus one. Equality on the whole input
space includes arbitrary external references. No supplied eigenstate,
clean coefficient, measurement or postselection appears.

## 2. Complete query ledger

Suppose the source and its actual native inverse each have simultaneous
bounds $`T_F,G_F,D_F,B_F`$, with zero clean work and B_F extra dirty wires.
Each J or its inverse uses two source occurrences and a address CNOTs.
There are four J occurrences in (3), so the completed query has

```math
\begin{aligned}
T_Q&\le8T_F,\\
G_Q&\le8G_F+4a+2|v|+10+4\beta,\\
D_T(Q)&\le8D_F,\\
B_Q&\le(2^{a+1}-1)+B_F,\qquad a_{\mathrm{clean}}=0.
\end{aligned}                                                   \qquad(4)
```

B_Q excludes the original address and its arbitrary output. The source
helper pool is reused only between completed returning source calls.
The extra ten Cliffords are four Z words, each $`S^2`$, and two output
Hadamards. The constant X is emitted as HSSH. Reversal exchanges every
source with its actual inverse and reverses all CNOTs; adjacent inverse
S powers reduce literally to the same Z and X words. Thus the actual
query inverse has the same displayed ledger. Its Boolean action is also
an involution, but that fact is not used to omit an inverse call.

The existing source gives $`T_F\lt 48\,3^r`$, $`G_F\lt 60\,3^r`$,
$`D_F\le40\lceil\log_2r\rceil`$ and $`B_F\lt 2\,3^r`$.
Sections 3–5 supply a different fully paid source with smaller asymptotic
count and width. Either complete circuit may be selected. No source
call in (4) is a free preprocessing operation.

## 3. Linear preprocessing without zero buffers

A binary XOR straight-line program has s two-input XOR nodes. Each node
may refer only to original inputs, earlier nodes or the formal zero. Its
R outputs are references of the same kinds. There is a Clifford shear

```math
(x,A,d)\longmapsto(x,A+Lx,d)
```

using R arbitrary output bits A and s arbitrary temporary bits d, at most
$`8s+2R`$ CNOTs, and no initialized bit. Its actual inverse has the same
count and workspace.

To prove this, assign each node j its own arbitrary d_j and apply CNOTs
from its two referenced values into d_j in program order. Formal-zero
references emit no gate. Call this word C. It has at most 2s CNOTs and
an exact linear action $`d\mapsto Mx+Rd`$, where R is invertible.
The word C_0 obtained by omitting all original-input CNOTs implements
$`d\mapsto Rd`$ and also has at most 2s CNOTs.

Execute C, copy the designated output references into A, and execute the
actual C inverse. Next execute C_0, copy only the designated temporary
references into A, and execute the actual C_0 inverse. Direct-input
output references are copied in the first half only; structural-zero
outputs are never copied. The two temporary contributions Rd cancel,
and A gains exactly Lx. Both halves restore d. This proves the bound and
handles all initial temporary values and their correlations.

The lemma compiles a small XOR program itself; it does not expand each
output parity into its full support. Consequently fast zeta transforms
need not incur the $`3^t`$ flattened support count. A formal zero in that
program is an absent source term, never a physical qubit promised zero.

## 4. Fast subset convolution as a native returned query

The ranked zeta and rank-convolution algorithm is due to Björklund,
Husfeldt, Kaski and Koivisto, [*Fourier Meets Möbius: Fast Subset
Convolution*](https://arxiv.org/html/cs/0611101v1), Sections 2.3–2.4
(STOC 2007). It works over arbitrary rings, including $`\mathbb F_2`$.
The contribution here is its full-input dirty-work native compilation.

For $`t\ge0`$, put

```math
m=2^t,\qquad R=\binom{t+2}{2}2^t,\qquad
s=\begin{cases}t(t+1)2^{t-1}&t\ge1,\\0&t=0.\end{cases}
```

Rows are triples $`(T,i,j)`$, where $`T\subseteq[t]`$ and
$`i,j\ge0`$, $`i+j\le t`$. Define linear forms on arbitrary coefficient
arrays a,b and an arbitrary output-phase array z by

```math
\begin{aligned}
(Ua)_{T,i,j}&=\sum_{A\subseteq T,\ |A|=i}a_A,\\
(Vb)_{T,i,j}&=\sum_{B\subseteq T,\ |B|=j}b_B,\\
(W^Tz)_{T,i,j}&=\sum_{S\supseteq T,\ |S|=i+j}z_S.
\end{aligned}
```

Then, with all sums in the binary field,

```math
\sum_{T,i,j}(Ua)_{T,i,j}(Vb)_{T,i,j}(W^Tz)_{T,i,j}
 =\sum_S z_S\sum_{A\sqcup B=S}a_Ab_B.                 \qquad(5)
```

For fixed A,B,S the left coefficient counts T between $`A\cup B`$ and
S, subject to $`|A|+|B|=|S|`$. The count is odd exactly when
$`A\cup B=S`$; the cardinality condition then forces disjointness.
This also proves (5) in characteristic two without a division step.

Each ranked subset or superset zeta transform has at most s two-input
XOR nodes: there are t directions, $`2^{t-1}`$ updates per direction,
and t+1 ranks. Initial ranked entries are direct input references when
their sizes agree, and formal zeros otherwise. The R final rows are
references to these transformed entries, including repetitions. Thus
section 3 supplies each of the three mask shifts at cost $`8s+2R`$,
with the same reusable s temporary bits.

Reserve three separate arbitrary R-bit masks A,B,C. Shift them by Ua,
Vb and $`W^Tz`$, respectively. On the masks use the literal phase

```math
P(A,B,C)=(-1)^{\sum_{\ell=1}^R A_\ell B_\ell C_\ell}.
```

The R CCZ gates act on disjoint physical triples. The retained exact
[CCZ parity word](RETURNED_DIAGONAL_ECHO_DEGREE.md#sharp-cubic-example-and-a-literal-native-ledger)
uses seven T gates, ten CNOTs and at most five T layers per triple;
all triples run in those same five layers.

Apply P at the eight translated mask vertices

```text
000, 100, 110, 010, 011, 111, 101, 001
```

and return to 000 afterwards. The edge shifts are A,B,A,C,A,B,A,C.
Use the actual inverse shift on each returning edge. This costs eight
complete P words and eight complete linear shears. The sum of all eight
cubic phases is their third finite difference,

```math
(-1)^{\sum_\ell (Ua)_\ell(Vb)_\ell(W^Tz)_\ell},
```

independent of all masks. All masks and linear-program temporaries
return. Original a,b,z wires are used only as Clifford controls in
this phase circuit. Finally conjugate by Hadamards on all m output z
wires. Equation (5) proves the full product XOR

```math
|a,b,z,d\rangle\longmapsto|a,b,z+ab,d\rangle
```

with simultaneous native bounds

```math
\begin{aligned}
T&\le56R,\qquad D_T\le40,\\
G&\le96R+64s+2m,\\
B_{\mathrm{extra}}&\le3R+s,\qquad a_{\mathrm{clean}}=0.       \qquad(6)
\end{aligned}
```

The 80R CNOTs in the eight P words and $`8(8s+2R)`$ in the linear
shears account for the Clifford bound, including the final 2m Hadamards.
The actual inverse has the same ledger: the emitted alphabet is H,
CNOT, T and T-dagger. The masks are not copied values of any source.
Source banks may overlap each other, but every output, mask and temporary
bank is disjoint from every source; all nonzero product inputs must be
present physical wires. A structurally absent coefficient is omitted
from every linear form. This includes a b0 bank whose constant is absent.

## 5. A smaller-count complete HRT source

Replace each bilinear call in the balanced source by (6). The original
Clifford copies and shared correction identity are unchanged.
For a merge with N bank pairs and suffix dimension t, there are $`4N+2`$
product instances. The first N products and the Q product have disjoint
outputs and only unchanged shared sources, so their mask circuits run
in parallel with separate masks and temporaries. The same five-stage
schedule remains valid, now with forty T layers per stage.

Since $`s\le R`$, $`m\le R`$, and $`N\ge1`$, one merge and its actual
inverse obey the conservative bounds

```math
T\le336NR,\qquad G\le978NR,\qquad
D_T\le200,\qquad B_{\mathrm{extra}}\le9NR.                \qquad(7)
```

For G, the product calls cost at most $`(4N+2)162R`$; the unchanged
copies cost $`2Nm+4m-2`$. For work, reserve the m-bit correction bank
outside $`N+1`$ independent sets of three masks and linear temporaries:
$`m+(N+1)(3R+s)\le9NR`$. Every stage returns these scratch registers.

Put $`h_r=\lceil\log_2r\rceil`$ and
$`A_r=\binom r2 2^r`$ for $`r\ge2`$. In the balanced tree, a node with
left endpoint l and left length L has $`N=2^L-1`$ and
$`t=r-l-L\le r-2`$. Thus

```math
NR\lt \binom r2 2^{r-l}.
```

At one tree level the left endpoints are distinct, so their sum is less
than A_r. There are h_r levels. Nodes at one level have disjoint
coefficient banks. All helpers return before another level begins, so
one extra pool can be reused between levels. The complete source and
its actual inverse therefore have

```math
\begin{aligned}
T_F&\lt 336h_r A_r,\qquad G_F\lt 978h_r A_r,\\
D_F&\le200h_r,\qquad B_F\lt 9A_r,\qquad a_{\mathrm{clean}}=0.
\end{aligned}                                                   \qquad(8)
```

At r=1 the source remains identity. This improves asymptotic T/G count
from $`O(3^r)`$ to $`O(r^2 2^r\log r)`$, and extra dirty width to
$`O(r^2 2^r)`$. It pays a larger constant in the same logarithmic depth.
It makes no small-instance or optimal-count claim. Using (8) in (4)
gives a full arbitrary-table query of depth $`O(\log(a+2))`$ and
count $`O(a^2 2^a\log(a+2))`$, not the square-root count, O(table-size)
Clifford budget or sublogarithmic depth required to improve R98.

## 6. Two precise limitations of further arithmetic compression

### Deferring arbitrary helper garbage does not remove inverse degree

Suppose a full-input classical reversible computation implements
$`V(x,w)=(k(x),g_x(w))`$, where k is the quadratic HRT inverse and w
is an entirely arbitrary dirty bank. Each g_x is a bijection of its
whole fiber: V is injective and both fibers have the same finite size.
Consequently

```math
V^{-1}(y,v)=(f(y),g_{f(y)}^{-1}(v)).
```

The inverse principal coordinate has degree exactly r, independently of
all final helper values. If a word has D completed reversible stages
whose inverse coordinate degrees are at most delta, composition gives
$`\deg V^{-1}\le\delta^D`$. Therefore

```math
D\ge\lceil\log_\delta r\rceil\qquad(\delta>1).       \qquad(9)
```

Arbitrary affine coordinate changes and arbitrarily many dirty wires do
not change this argument. Bilinear XOR shears and the completed
involutory generator-split multipliers have delta at most two. Thus
keeping their helper garbage live does not collapse all balanced levels
into constantly many such completed stages. This is an architectural
bound on classical completed stages, not an unrestricted native T-depth
bound. Intermediate nonmonomial interference can leave its scope; a
principal-coordinate promise only on a clean-zero sector also lacks the
full-fiber argument.

### A full-degree table prevents this algebra quotient from shrinking

Let the coefficient of the full a-variable monomial be one. On $`R_a`$
define $`\ell(e_S)=\alpha_S`$ for nonempty S and choose its value at the
unit arbitrarily. The product pairing $`(u,w)\mapsto\ell(uw)`$ is
nondegenerate. For nonzero u, choose a support subset A of minimum
cardinality and multiply by $`e_{[a]\setminus A}`$. A term indexed by B
can survive only if $`B\subseteq A`$. Proper such subsets have zero
coefficient by minimality, leaving precisely $`e_{[a]}`$ with coefficient
one. Applying ell gives one.

If an ideal I were quotiented out while preserving all generator products
and this table functional, every $`u\in I`$ would satisfy
$`\ell(uw)=0`$ for every w, forcing $`I=0`$. Hence this coefficient-algebra
quotient route retains dimension $`2^a`$ for every full-degree table.
A one-qubit-per-coefficient carrier cannot then reach R98's square-root
width allocation by quotienting alone. This does not constrain unrelated
representations, fused native circuits, or the unrestricted model.

## 7. Evidence

The focused fixtures verify scalar extraction on arbitrary coefficient
inputs, output control and inverse signs; every physical column of the
five-wire native table query and its actual adjoint; all-input dirty linear shears
at the first nontrivial rank; the subset-convolution factorization; the
full six-wire native XOR and adjoint at t=0; larger native phase words
with arbitrary masks and temporary values; and the balanced tree ledger.
The arbitrary-width statements follow from the identities above. No
complete-frame theorem or unrestricted depth lower bound is promoted.
