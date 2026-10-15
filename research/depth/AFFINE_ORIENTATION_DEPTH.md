# Logarithmic depth for affine branch orientations

[Homogeneous common-flag core](HOMOGENEOUS_MAGNITUDE_DEPTH.md) · [Retained sign and swap gauges](../endpoint/PRECISION_GAUGE_COMPILER.md#1-contract-and-effective-normalization) · [Generic frontier](UNIFORM_SUBTREE_DEPTH.md)

**Local result R83.** A complete Hopf frame whose unordered rotation
magnitudes depend only on depth still has fixed-accuracy logarithmic
T-depth when each node's branch orientation is an affine Boolean function
of its prefix. Arbitrary node signs are permitted. The leaf relabeling is
an affine triangular permutation. Its marker correction is implemented
by a short-tag-selected affine shear, using private dirty phase queries;
it is not an unpriced selected Clifford circuit.

This enlarges R80's structured family and leaves the unrestricted R79
frontier unchanged. The proof includes the full marker signs, controlled
linear-map schedule, exact tag and helper return, certified rounding,
recursive prefix dependence, and every original dirty allocation.

## 1. Family and fixed-accuracy theorem

Number logical positions from 1 to n. At each depth j there is a
nonnegative normalized pair $`(a_j,b_j)`$. At a node with prefix
$`p\in\mathbb F_2^{j-1}`$, assume

```math
(|c_{j,p}|,|s_{j,p}|)=
\begin{cases}(a_j,b_j),&r_j(p)=0,\\(b_j,a_j),&r_j(p)=1,
\end{cases}
\qquad a_j^2+b_j^2=1,
```

where the orientation is affine over $`\mathbb F_2`$:

```math
r_j(p)=\rho_j+\sum_{k\lt j} A_{jk}p_k.
```

All binary sums in this chapter are modulo two. Magnitudes may vary with
j; zeros and equal pairs are allowed. Node signs are arbitrary. The
promise is existence of these finite affine masks, not an efficient
procedure for discovering or certifying the promise. If the masks are
not supplied, Section 2 gives terminating preprocessing that preserves
the family without exact orientation, sign, or zero tests.

Fix an integer $`L\ge6`$, let $`\eta=2^{-L}`$ and $`N=2^n`$.
For every $`n\ge1`$ and every integer dirty allocation in the following
table, one external initialized qubit suffices:

| Prescribed operator W | Required arbitrary dirty width |
|---|---|
| Complete real frame with affine orientations | $`b\ge17(L+n+7)`$ |
| Prescribed complex magnitude completion $`D_\phi W_{\mathbb R}`$ with that real-frame promise | $`b\ge17(L+n+8)`$ |

There is one Clifford+T circuit V with simultaneous guarantees

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\log_2(n+2)).
```

Here $`J_1`$ initializes the sole external qubit. Constants depend only
on fixed L, not on the angles, affine masks, dimension, or allocation.
The norm includes every logical column, literal phases, work leakage and
return, arbitrary dirty inputs, and external references. Unused supplied
wires stay untouched. G counts elementary Clifford gates, including the
arbitrary Clifford interlayers allowed by the T-depth model. All coherent
queries, guards, preparations, actual inverses, and recovery use priced
native words. There is no measurement, reset, QRAM, supplied precision
state, or free frame/permutation oracle. Classical preprocessing is
terminating without an efficiency claim.

At $`b=\Theta(\sqrt N)`$ the depth is $`O_L(\log(n+2))`$ with
square-root T-count and linear Clifford count. This is a supporting-family
upper bound, not a matching count or depth lower bound for that family.
There is no all-precision assertion with absolute constants.

## 2. Finite rounding with affine orientations and arbitrary signs

Set $`\epsilon=\eta/(64n)`$. At each depth enumerate a finite angle
grid in $`[0,\pi/2]`$ of mesh at most $`\epsilon/2`$ and all affine
Boolean masks on that depth's prefix. For each grid angle $`\phi'_j`$
and candidate mask $`r'_j`$, run certified strict tests that every node's
magnitude pair differs coordinatewise by less than epsilon from

```math
\begin{cases}(\cos\phi'_j,\sin\phi'_j),&r'_j(p)=0,\\
(\sin\phi'_j,\cos\phi'_j),&r'_j(p)=1.
\end{cases}
```

Dovetail the finitely many tests. A true affine mask and a grid point
within $`\epsilon/4`$ in angle pass with a strict margin, so this
search terminates. Dyadic multiples of $`\pi`$ give finite exactly
normalized angle descriptions. Ambiguous orientations at a balanced
pair do not obstruct the search; it returns one finite affine mask.
This is a global test of each candidate mask, not independent node
rounding that could destroy affinity.

At each node select the two signs using rational approximants to its
signed coefficients with absolute error at most epsilon. Multiply the
accepted oriented magnitude pair by those selected signs to obtain
$`(c'_v,s'_v)`$. A wrong sign implies original magnitude at most epsilon,
so each signed coordinate error is at most $`3\epsilon`$. The local
real-rotation error is at most $`3\sqrt2\epsilon`$. Each logical
layer is a direct sum of node rotations, giving

```math
\|W-W'\|\le3\sqrt2 n\epsilon\lt\eta/8.
```

The finite frame $`W'`$ has exactly affine orientations and known signs.
Use the [retained sign gauge](HOMOGENEOUS_MAGNITUDE_DEPTH.md#2-separate-signs-without-a-sign-oracle)

```math
W'=D_{\rm sign}W_{\rm mag}S_{\rm sign}.
```

The output leaf sign is the product of edge signs along that path. The
root input sign is one; the marker at v has input sign equal to the
ancestor-sign product times the two edge signs at v. This identity uses
only $`c'_v=e_0(v)|c'_v|`$ and $`s'_v=e_1(v)|s'_v|`$, so either
selected sign is valid on a zero edge. It covers every marker, including
those below zero incoming branches. Each exterior sign is an exact
Boolean phase with the all-width native ledger from R80, Section 7.
The rest of the proof treats the positive finite frame $`W_{\rm mag}`$.

## 3. Full-column leaf and marker factorization

For a local h-bit positive frame, let $`F(x)=y`$ be the tree
relabeling determined recursively by

```math
y_j=x_j+r_j(y_{\lt j}).
```

It preserves prefix length and interchanges the two children exactly
where the physical orientation says to do so. With A the strictly lower
triangular coefficient matrix and rho its constant vector,

```math
M=I+A,\qquad L=M^{-1},\qquad
F(x)=Lx+v,\quad v=L\rho,\qquad
G(y)=F^{-1}(y)=My+\rho.
```

Thus the leaf permutation P, $`P|x\rangle=|F(x)\rangle`$, is affine
and unitriangular. Let $`W_0`$ be the complete homogeneous frame with
pair $`(a_j,b_j)`$ at every node of depth j.

Write a reference marker as $`u1\,0^{h-t}`$, where $`|u|=t-1`$.
Let $`p=F_{\lt t}(u)`$ be its physical prefix. The permutation P sends
the reference scaling vector at u to the oriented scaling vector at p.
This follows inductively because a swapped child has its two coefficients
and its descendant leaf sets swapped together. For the marker complement,
a child swap gives the literal minus sign:

```math
P e_u^0=(-1)^{r_t(p)}e_p.
```

Indeed without a swap the two coefficients are unchanged. With a swap,
$`-b_t|0\rangle v_0+a_t|1\rangle v_1`$ becomes
$`a_t|0\rangle v'_0-b_t|1\rangle v'_1`$, the negative of the
oriented complement. Descendant permutations already map their scaling
vectors correctly. The root scaling vector has no sign. No nonzero
coefficient is divided out, so this argument includes zero pairs and
zero incoming branches.

Define the signed marker permutation Q by

```math
Q|0^h\rangle=|0^h\rangle,\qquad
Q|p1\,0^{h-t}\rangle
=(-1)^{r_t(p)}|G_{\lt t}(p)1\,0^{h-t}\rangle.
```

The full-column factorization is exactly

```math
W_{\rm mag}=P W_0 Q.
```

It is the retained swap/sign normalization with its marker action made
explicit. The chronological order is Q, then $`W_0`$, then P.

A preserved outside prefix p changes the local constant vector to an
affine function of p. The local matrix M remains unitriangular and
independent of p. Hence the same formulas give $`F_p,G_p,Q_p`$, and
the affine-orientation family is closed under taking conditional subtrees.

## 4. A native constant-depth guard for an affine triangular map

A linear unipotent map over $`\mathbb F_2`$ is a product of two
linear involutions. Here is a constructive proof. In a binary Jordan
basis, a block J is multiplication by $`1+t`$ on
$`\mathbb F_2[t]/(t^k)`$. Define the linear substitution

```math
\mathcal A f(t)=f\!\left(\frac{t}{1+t}\right)\pmod{t^k}.
```

The inverse of $`1+t`$ is a finite polynomial in this quotient.
Substitution twice returns t, so $`\mathcal A^2=I`$, and
$`\mathcal A J\mathcal A=J^{-1}`$. Therefore
$`J=\mathcal A(\mathcal A J)`$ is a product of two involutions.
Take the direct sum over Jordan blocks and conjugate back by the binary
change of basis. All matrices are finite classical data.

For an involution K write $`K=I+B`$. Characteristic two gives
$`B^2=0`$. A basis adapted to $`\mathrm{im}\,B\subseteq\ker B`$
conjugates K to disjoint two-coordinate CNOTs and fixed coordinates.
Every invertible binary basis change has an elementary Clifford CNOT/SWAP
implementation with polynomial count. To enable K with a preserved bit
E, apply that basis change unconditionally, the E-controlled CNOT bank,
and its actual inverse. The middle bank is a common-control Toffoli bank
on disjoint remaining pairs. The retained
[exact four-T-layer schedule](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
gives constant T-depth, linear T-count, and linear Clifford count for it;
Hadamards on its target wires reduce it to the displayed CCZ schedule.
It needs no initialized control copies or helpers. The two involution
banks are sequential, still constant depth.

An affine translation is a bank of E-controlled X gates, hence Clifford.
Thus E-controlled P or G has $`D_T=O(1)`$, polynomial Clifford count,
linear T-count, and literal identity on the entire E-zero sector.
The unconditional basis changes cancel there; E is never inside those
changes. No generic gate-by-gate control of a long Clifford circuit is
assumed.

For a local map depending affinely on a preserved outside prefix p,
view the whole linear part on $`(p,X)`$ as

```math
\begin{pmatrix}I&0\\B&M\end{pmatrix}.
```

It is again unitriangular. The same guard fixes p at the complete
interface. Intermediate Clifford basis changes may use p, but their
actual inverses restore it before another macro reads it. There is no
long-prefix lookup and no additional initialized flag.

The same interface permits E to be the conjunction of any fixed number
of preserved enable bits, disjoint from the data/prefix coordinates
used by the basis changes. No clean register for their conjunction is
required. Replace each middle controlled CNOT by the exact read-only
conjunction word on those enable bits and its pair's control, with a
private target/helper pool. Its arity is constant, so its count, width,
and T-depth are constant. Distinct pairs have disjoint non-Clifford
supports, even though their public enable bits are shared. Affine
translations use the same words and then have constant T-depth, rather
than necessarily being Clifford. All private helpers return. Thus every
constant-size enable predicate has the same asymptotic guarded-map
contract and full occupied-space inactivity.

## 5. The short-tag-selected marker correction

Use the more general local notation $`F(x)=Lx+v`$,
$`G(y)=My+w`$, with $`M=L^{-1}`$ and $`w=Mv`$.
For the moment fix the outside prefix, so v,w are constants.
Let t be the last-one position of a nonzero input
$`y=p1\,0^{h-t}`$. After G the prefix is
$`z_{\lt t}=u=G_{\lt t}(p)`$ and, for $`j\ge t`$,

```math
z_j=M_{j,\lt t}p+M_{jt}+w_j.
```

Since $`p=L_{\lt t,\lt t}u+v_{\lt t}`$, define finite binary data

```math
c_{jk}(t)=[M_{j,\lt t}L_{\lt t,\lt t}]_k,
\qquad k\lt t\le j,
```

```math
d_j(t)=M_{j,\lt t}v_{\lt t}+M_{jt}+w_j+[j=t].
```

Let $`C_t`$ fix the prefix and toggle its tail by

```math
z_j\longmapsto z_j+\sum_{k\lt t}c_{jk}(t)z_k+d_j(t),
\qquad j\ge t.
```

On $`G(y)`$ this produces $`u1\,0^{h-t}`$. On every input it is
an affine involution because all controls are in the fixed prefix.
For $`t=0`$, take $`C_0`$ to be translation by $`G(0)=w`$; this
returns the root to zero. Padded tags above h may use identity.
The required marker sign can be evaluated on the reconstructed prefix:

```math
r_t(p)=v_t+L_{t,\lt t}u.
```

This is just the t-th triangular equation for F, and is another affine
function selected by t. It is zero by definition at the root tag.

### Implement the selected shear by a quadratic phase

Let T be an s-bit tag, $`s=\lceil\log_2(h+1)\rceil`$, and let
$`H_{\rm tail}(T)`$ apply Hadamards to exactly positions $`j\ge t`$,
with all positions selected at t zero and none at padded tags. The
[short-label dirty-mask echo](SHALLOW_DIRTY_QUERIES.md#3-a-hadamard-layer-selected-by-a-short-label)
implements this exact selected involution with polynomial count/width
and $`O(\log\log(h+4))`$ T-depth. Its tag is read only, and every
dirty mask and helper returns on all inputs.

For a preserved activity E, define the diagonal phase

```math
D_E(T,z)=(-1)^{E\sum_{j\ge t}z_j
 (\sum_{k\lt t}c_{jk}(t)z_k+d_j(t))}.
```

The chronological word

```math
H_{\rm tail}(T),\quad D_E(T,z),\quad H_{\rm tail}(T)^\dagger
```

implements $`C_t`$ on E one and identity on the full E-zero space.
The two tail words are unconditional in E and cancel with their actual
inverses when the central diagonal is inactive.

Each quadratic term is a Boolean phase with address
$`(T,E,z_j,z_k)`$ of at most $`s+3`$ bits. Each linear term uses
$`(T,E,z_j)`$. Its coefficient is set to zero for irrelevant tags.
Use a distinct arbitrary dirty output and private helper pool for every
phase query, followed by its exact $`Q_f,Z,Q_f^\dagger,Z`$ echo.
All X,T,E wires are read-only controls through the complete diagonal.
The shallow query's non-Clifford supports are private, so every monomial
runs in parallel. Counts of the shared Clifford gates are summed.
For h local bits the safe bounds are

```math
T,w=O(h^{5/2}),\qquad G=O(h^3),\qquad
D_T=O(\log\log(h+4)).
```

The sign $`(-1)^{r_t(p)}`$ uses short-tag-selected affine phase terms
with the same or smaller bounds. No selected Clifford circuit or clean
parity word is supplied.

With an outside prefix p of length at most n, v,w are affine in p, so
$`d_j(t,p)`$ and the sign have additional terms involving individual
prefix bits. Query each such monomial using address $`(T,E,z_j,p_i)`$,
not the whole prefix. There are $`O(h(n+1))`$ terms in all. Thus

```math
T,w=O(n(h+2)^{3/2}),\qquad G=O(n(h+2)^2),
```

with the same local-depth bound. These explicit long-prefix terms are
why repeated recursive stages need no repeated $`\log n`$ query depth.

### Compute and erase the marker tag

The [exact last-one XOR](HOMOGENEOUS_MAGNITUDE_DEPTH.md#4-exact-last-one-xor-with-arbitrary-dirty-work)
$`F_E`$ toggles arbitrary T by E times the last-one index of X, returning
all arbitrary dirty work. It has $`O(\log(h+2))`$ T-depth and safe
$`O((h+2)^4)`$ count and width. It is identity on E zero, even on an
occupied tag.

On a tag initially zero when E is one, run chronologically

```math
F_E,\quad G_{p,E},\quad C_{E,T},\quad
\mathrm{Sign}_{E,T},\quad F_E.
```

The middle factors implement the signed marker map and preserve the
original last-one index. The final XOR returns the tag exactly. The
root uses t zero and the constant correction above. Every factor is
identity on E zero, including arbitrary occupied tags and work.
This implements $`Q_p`$ with exact dirty/tag return, polynomial count
and width, and $`O(\log(h+2))`$ T-depth. It holds by linearity on all
logical, dirty, and reference inputs in its stated initialized-tag sector.

## 6. A complete early frame using one shared tag

Suppose E and an outside prefix p are preserved. E may be one bit or
the conjunction of a fixed number of separate preserved enable bits,
as in Section 4. Those bits are disjoint from X,p,T and the helper pools.
In every query below include the enable bits in its read-only address;
this adds only a constant number of address bits. The h-bit data register
is X. A disjoint tag T of $`\lceil\log_2(h+1)\rceil`$ bits is zero
when E is one and may be arbitrary when E is zero. All remaining helpers
are arbitrary dirty inputs.

First apply the exact Q word from Section 5, returning T. Then use the
same tag for the [actual common-flag word](HOMOGENEOUS_MAGNITUDE_DEPTH.md#5-one-protected-stage-and-exact-tag-return)
on all h data bits:

```math
F_E,\quad S_E(T),\quad B^\dagger,\quad F_E,\quad B.
```

Here choose native one-qubit words $`V_j`$ and define

```math
J_j=V_j ZV_j^\dagger,\qquad
\widehat U_j=J_jZ,\qquad B=\bigotimes_{j=1}^h\widehat U_j,
\qquad
S_E(T)=\bigotimes_j\widehat U_j^{E[t\le j]}.
```

The selected source uses a complete selected-Z involution echo followed
by a separate selected-J echo. Controlled J has chronology
$`V_j^\dagger,\mathrm{CZ},V_j`$; all its T gates are on private data
targets. B is unconditional. The actual $`B^\dagger`$ exposes the
marker as a last-one position, the second XOR erases T exactly, and B
restores the desired generalized homogeneous column. On E zero the
selected source and XORs are identity, so $`B^\dagger B`$ cancels on
all occupied inputs. This is an exact identity for the actual native
$`\widehat U_j`$, which need not be real. Scalar phases of $`V_j`$
cancel with its actual inverse.

Finally apply the guarded leaf map $`P_{p,E}`$. The active logical
operator is exactly $`P_pW_{\widehat U}Q_p`$. The tag and all dirty
helpers return exactly, and the inactive sector has literal phase one
on its full occupied space. No second clean signal or nested homogeneous
compiler is needed; the Q and homogeneous words reuse the same tag.

If $`V_j`$ approximates $`R(\phi_j/2)`$ within $`e_j`$, then
$`\|\widehat U_j-R(\phi_j)\|\le2e_j`$. Homogeneous layerwise
operator telescoping gives active error at most $`2\sum_j e_j`$.
The exterior exact P and Q do not enlarge it. For a requested local error
$`\delta`$, choose $`e_j\le\delta/(2h)`$. This gives the useful
standalone local contract

```math
D_T=O(\log(h+2)+\log(1/\delta)),
```

with polynomial count and dirty width in n,h, plus
$`O(h\log(h/\delta))`$ native one-qubit gates. One-qubit synthesis
is ancilla-free, so precision adds no width. Approximation changes the
logical frame only; tag, helper, and inactive-sector return are exact.

## 7. Protected recursion and the one-clean base

For current local size $`m\gt64`$ choose

```math
r=\lceil\log_2(m+1)\rceil+2,\qquad h=m-r.
```

Split the local register into X of size h and Z of size r. Compute
$`E=[Z=0^r]`$ into the sole external initialized qubit using the
[returned two-helper conjunction](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
Save E before occupying the tag. Select the required short tag inside Z
and apply Section 6's early circuit. On E one the tag returns to zero;
on E zero its completed word is identity even if Z is occupied. Thus Z
is preserved, and the actual predicate inverse returns E exactly.

The remaining r logical layers are conditional subtree frames with the
completed upper logical bits as an outside prefix. Their orientation
masks remain affine in that prefix and the remaining local prefix.
Recurse on Z with the same sole external clean qubit and returned dirty
pool. The original depth-dependent pair for each remaining depth is
unchanged. Every original layer belongs to exactly one stage.

For $`m\le64`$ use direct local layers. At target j compute its
unchanged lower-suffix-zero predicate into the clean flag E. Its
orientation r is an affine parity of the preserved outside/local prefix.
An exact query into an arbitrary dirty mask d,

```math
d\longmapsto d+Er,
```

has constant T-depth without a clean parity register: if the parity
support is nonempty, compute it by CNOTs into one participating prefix
wire, use a Toffoli controlled by E and that wire to toggle d, and undo
the CNOT network. A constant parity adds at most CNOT(E,d). The full
query restores the entire prefix. Its elementary Clifford count is
$`O(n)`$ and its T count is $`O(1)`$. The target and lower suffix are
excluded from that parity support.

Run a complete selected-Z echo for Er, then E-controlled
$`\widehat U_j=J_jZ`$, then a complete selected-X echo for Er.
Chronological order gives

```math
E=1:\quad X^r\widehat U_jZ^r,\qquad E=0:\quad I.
```

The exact identity $`XR(\phi_j)Z=R(\pi/2-\phi_j)`$ proves equal
approximation error for both orientations. Prefix bits are restored
before each payload, and the target never enters the mask address.
Every dirty mask returns. Finally uncompute E while its suffix is
unchanged. An empty suffix uses $`X_E`$ before and after the controlled
word. There are at most 64 base layers, so their depth is $`O_L(1)`$
at the precision used below, with polynomial Clifford count in n and
only a constant dirty reservation. This closes the base without an extra
initialized orientation bit.

## 8. Summable errors and simultaneous resources

Let $`m_0=n`$ and $`m_{i+1}=\lceil\log_2(m_i+1)\rceil+2`$ until
the base is at most 64. At a stage of size m, choose each assigned
one-qubit word with

```math
e_j\le\frac{\eta}{256(m+2)^2}.
```

For $`m\gt64`$ the next size is at most $`m/2`$, and logarithms
of successive sizes decrease geometrically above an absolute cutoff.
The reciprocal-size sum is bounded by an absolute constant. In
particular,

```math
\sum_i\frac{m_i}{(m_i+2)^2}\le2,\qquad
\sum_i\log_2(m_i+2)=O(\log_2(n+2)).
```

The operator error of each stage is bounded by twice the sum of its
assigned $`e_j`$. Exact P/Q factors and preserved-prefix direct sums
retain this bound. Unitary telescoping over stages, including the base,
gives native error at most $`\eta/64`$. Combining with Section 2's
finite perturbation leaves total real error below eta. All comparisons
are full operator or initialized-isometry bounds; previous actual error
is propagated unitarily, never reset or postselected away.

At stage size m, the tag XOR costs $`O(m^4)`$ gates and dirty width;
selected shear phases cost $`O(n(m+2)^{3/2})`$ T/width and
$`O(n(m+2)^2)`$ Clifford gates; guarded affine maps cost polynomial
Clifford gates in n and $`O(n)`$ T gates. Single-qubit words add
$`O(m[L+\log(m+2)])`$ gates. Thus, for fixed L, one circuit has

```math
D_T=O_L(\log_2(n+2)),\qquad w\le C(n+2)^4,
```

```math
T,G=O((n+2)^4+n[L+\log_2(n+2)]).
```

Here C is an absolute width constant. To see the count sum, the powers
of local m sum geometrically; the global-prefix contributions sum to at
most $`O(n^3)`$. The $`O(n^2)`$ Clifford basis changes at each of
at most $`O(\log(n+2))`$ stages are also below the displayed fourth
power. Consecutive completed macros and recursive stages reuse their
returned pools. The width includes every private phase output/helper,
mask output, and predicate helper. The short tag is protected logical
work, not an extra external initialized allocation.

Physical T-layer parallelism is part of this ledger: selected-J words
use distinct data targets; each phase monomial uses private query work;
last-one conjunctions have private non-Clifford supports; and controlled
linear involutions use the exact common-control CCZ schedule. Shared
Clifford controls are serialized where needed and counted. No clean
control clones or constant elementary-depth fanout are assumed.

## 9. Signs, all widths, and prescribed phases

Each exterior sign in Section 2 is an exact Boolean phase. The retained
[all-width sign-query splice](HOMOGENEOUS_MAGNITUDE_DEPTH.md#7-exact-sign-phases-and-every-dirty-width)
uses one arbitrary dirty output, no initialized phase bit, and

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O(N/b^2+\log_2(n+2))
```

at every $`b\ge17(L+n+7)`$. Its complete echo and all actual inverses
return arbitrary dirty/reference inputs. Run S, the positive circuit,
and D sequentially, reusing the pool.

Set $`b_c=\sqrt N/n`$. Above an absolute dimension cutoff the
polynomial pool $`C(n+2)^4`$ fits $`b_c`$. For $`b\ge b_c`$ use
the structured recursive circuit. Its fixed-L polynomial count is
absorbed by $`O_L(\sqrt N)`$ T and $`O_L(N)`$ Clifford gates. For
$`b\lt b_c`$, select R79 on the original target at the original L.
It fits the literal threshold, and its fixed-L additive
$`O_L(\log^2(n+2))`$ depth is absorbed by $`N/b^2\gt n^2`$.
For the finitely many dimensions below the cutoff use R79 at all widths;
its additive term is bounded by a constant multiple of $`\log(n+2)`$
on that finite set. No hidden helper constant is presumed to fit the
coefficient 17. Each branch chooses one complete circuit supplying all
counts, depth, width, and error guarantees simultaneously.

The structured branch returns the sole clean signal exactly. The R79
fallback retains its full initialized-isometry error contract. For a
prescribed complex magnitude frame, run the real theorem and the
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
each at precision $`L+1`$. Both fit $`b\ge17(L+n+8)`$. The diagonal
has the same fixed-L count, Clifford, and depth orders. A two-factor
unitary hybrid adds the two half-error bounds while reusing the same
clean signal and dirty pool; actual intermediate leakage is propagated
by the next actual unitary. Every specified leaf phase is retained,
including phases at zero-amplitude leaves. This is not arbitrary
complex-unitary synthesis.

## 10. Scope and verification boundary

The sign/swap normalization, one-qubit approximation, shallow returned-
dirty queries, two-helper predicates, common-control CCZ schedule, and
R80 common-flag word are retained ingredients. The local additions are
the explicit affine marker shear, its short-address quadratic-phase
implementation, the guarded unipotent factorization, and the recursive
long-prefix ledger. These are analytic native constructions, not a claim
that a scalable general compiler is implemented by a finite fixture.

The three affine fixtures in
[test_orientation_depth.py](../../tests/test_orientation_depth.py) check
the marker root/sign/tail formulas for fixed outside-prefix choices,
the complete shared-tag stage against independent path-product frame
columns, and a small guarded affine word with literal Clifford+T
Toffoli phases. The shared-tag check uses nonreal actual
$`(VZV^\dagger)Z`$ factors. It includes every bounded logical/dirty
input column with initialized active tags, and every occupied inactive
tag/dirty column, including padded tags. A separate counterexample checks
the changed-address failure of a simultaneous dirty triangular update.
The queries and Clifford basis changes are supplied finite matrices;
these checks do not establish large-bank parallel scheduling, recursive
resource bounds, finite cutoffs, or a scalable native emitter. The same
file's signed deterministic fixture supports R84, not a wider R83 claim.

For arbitrary orientation tables the leaf relabeling is an arbitrary
triangular permutation, and its marker tail is generally nonlinear in
the transformed prefix. The short-address monomial decomposition and
constant-depth unipotent guard then do not follow. A complete parallel
XOR query into arbitrary dirty outputs does not supply an in-place
permutation or a clean output register. Consequently this theorem makes
no improvement to the general R79 frontier and no unrestricted
permutation or depth lower-bound claim.
