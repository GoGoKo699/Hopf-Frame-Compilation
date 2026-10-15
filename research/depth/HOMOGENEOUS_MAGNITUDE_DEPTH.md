# Logarithmic depth for homogeneous magnitudes

[Generic one-clean theorem](UNIFORM_SUBTREE_DEPTH.md) · [Shallow dirty queries](SHALLOW_DIRTY_QUERIES.md) · [Sign normalization](../endpoint/PRECISION_GAUGE_COMPILER.md#1-contract-and-effective-normalization)

**Local result R80.** Complete Hopf frames whose rotation magnitudes
depend only on tree depth admit a fixed-accuracy logarithmic T-depth
compiler at square-root dirty width. Node signs may vary arbitrarily.
The construction transports every detail band by one common product of
actual native one-qubit words. A classical last-one XOR then returns the
tag exactly. Recursing on the protected logical suffix completes the
frame with one external clean qubit.

This is a structured-family result. It does not improve R79's generic
fixed-accuracy additive depth or prove a depth or count lower bound for
the structured family. The sign normalization is retained prior work;
the common native conjugator and its complete resource ledger supply the
new compiler.

## 1. Statement and input promise

Fix an integer $`L\ge6`$, put $`\eta=2^{-L}`$, and let
$`n\ge1`$, $`N=2^n`$. At each node p of depth d write the prescribed
real rotation as

```math
R(\theta_{d,p})=
\begin{pmatrix}c_{d,p}&-s_{d,p}\\s_{d,p}&c_{d,p}\end{pmatrix}.
```

Assume there are nonnegative numbers $`a_d,b_d`$ such that

```math
|c_{d,p}|=a_d,\qquad |s_{d,p}|=b_d,\qquad a_d^2+b_d^2=1
```

for every node at that depth. Signs, zeros, and magnitudes at different
depths are unrestricted. This is a promise on the prescribed input; the
compiler need not decide equality of certified real numbers. The input
uses the repository's terminating certified evaluation interface.

For either family below and every indicated integer allocation b, there
is one native Clifford+T circuit V with one external initialized qubit
and at most b arbitrary dirty qubits:

| Prescribed operator W | Required dirty width |
|---|---|
| Complete real Hopf frame with homogeneous magnitudes | $`b\ge17(L+n+7)`$ |
| Prescribed complex magnitude completion $`D_\phi W_{\mathbb R}`$ with that real-frame promise | $`b\ge17(L+n+8)`$ |

With $`J_1`$ initializing the sole external qubit, the same circuit has

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L\!\left(\sqrt N+\frac Nb\right),\qquad
G=O_L(N),\qquad
D_T=O_L\!\left(\frac N{b^2}+\log_2(n+2)\right).
```

The constants may depend on the fixed L, but not on n, b, or the angles
and leaf phases. The norm includes every logical column, literal phases,
all work leakage and return, arbitrary dirty inputs, and reference
correlations. Unused supplied wires stay untouched. G counts elementary
Clifford gates; their arbitrary interleaving is allowed in the T-depth
model and is charged in G. All queries, preparations, actual inverses,
guards, and recovery use priced native words. There is no QRAM,
measurement, reset, supplied precision state, or free coherent evaluator.
Classical preprocessing is terminating, without an efficiency claim.

At $`b=\Theta(\sqrt N)`$ this gives $`T=O_L(\sqrt N)`$,
$`G=O_L(N)`$, and $`D_T=O_L(\log(n+2))`$. Section 8 gives an
exact signed balanced corollary with absolute constants. The complex
family is the prescribed phase-dressed magnitude completion, not
arbitrary complex-unitary synthesis.

## 2. Separate signs without a sign oracle

The [retained quadrant normalization](../endpoint/PRECISION_GAUGE_COMPILER.md#1-contract-and-effective-normalization)
already separates local signs by leaf and column gauges. Its sign-only
specialization has a useful closed form. At a node v put

```math
e_0(v)=\mathrm{sgn}(c_v),\qquad
e_1(v)=\mathrm{sgn}(s_v),\qquad \mathrm{sgn}(0)=1.
```

Let $`\alpha_v`$ be the product of ancestor edge signs on the path to
v and $`d_x`$ their product along the complete path to leaf x. Define
the output diagonal $`D=\mathrm{diag}(d_x)`$ and the input
marker diagonal S by

```math
S_0=1,\qquad
S_{p1\,0^{n-|p|-1}}=\alpha_p e_0(p)e_1(p).
```

If $`W_+`$ replaces every local pair by its magnitudes, then exactly

```math
W=D W_+ S.
```

For completeness, let $`D_v`$ be the sign diagonal internal to the
subtree at v. It maps that subtree's nonnegative scaling vector to the
signed scaling vector. On its marker complement,

```math
D_v e_v^+
=-e_0(v)|s_v|v_{v0}+e_1(v)|c_v|v_{v1}
=e_0(v)e_1(v)e_v.
```

Global D restricts to $`\alpha_vD_v`$ there; its root multiplier is
one. This proves every column with its literal sign. No division is
used, so zero edges and markers below zero incoming branches are covered.
The same identity permits either assigned sign on a zero edge: it only
uses $`c_v=e_0(v)|c_v|`$ and $`s_v=e_1(v)|s_v|`$. Thus the finite
sign choices below need no special zero test.

Exact signs cannot be demanded from a general certified real input.
Instead set $`\epsilon=\eta/(64n)`$. At each depth use one
representative node and a finite angle grid in $`[0,\pi/2]`$ with
mesh at most $`\epsilon/2`$. Certified strict component tests find a
common finite angle $`\phi'_d`$ satisfying

```math
|\cos\phi'_d-a_d|\lt\epsilon,\qquad
|\sin\phi'_d-b_d|\lt\epsilon.
```

A grid point within $`\epsilon/4`$ in angle exists, and sine and
cosine are one-Lipschitz; consequently an accepting strict test
eventually terminates. Dovetail the finitely many tests. Dyadic multiples
of $`\pi`$ give finite descriptions of exactly normalized pairs,
including the endpoints. No equality or exact-zero test is performed.

At each node obtain rational approximations to $`c_v,s_v`$ with
absolute error at most epsilon, and select each sign by whether its
rational approximation is nonnegative. Set

```math
c'_v=e_0(v)\cos\phi'_d,\qquad
s'_v=e_1(v)\sin\phi'_d.
```

A wrong selected sign implies that the original coefficient has
magnitude at most epsilon. In all cases each coordinate error is at
most $`3\epsilon`$, so each local rotation error is at most
$`3\sqrt2\epsilon`$. A logical depth layer is a direct sum of
node rotations. Its error is their maximum, and unitary telescoping
across n layers gives

```math
\|W-W'\|\le3\sqrt2 n\epsilon\lt\eta/8.
```

The perturbed frame $`W'`$ has exactly homogeneous magnitudes and known
finite sign tables. Apply the exact gauge identity to $`W'`$. Thus the
remaining approximation problem is a positive homogeneous frame with
one normalized, finitely specified real rotation per depth. Rounding
each node independently would not preserve this promise.
In Sections 3–7, $`W_+`$ denotes this positive factor of the rounded
frame $`W'`$; the preceding perturbation bound accounts for its difference
from the original target.

## 3. Actual one-qubit words and selected involution banks

Number logical positions from 1 to m, most significant first. For any
one-qubit unitaries $`U_1,\ldots,U_m`$, define their generalized
homogeneous frame by the chronological layer circuit: at position j
apply $`U_j`$ exactly when all lower logical positions are zero. Its
root column and marker columns are respectively

```math
\bigotimes_{j=1}^mU_j|0\rangle,
\qquad
|p\rangle\otimes U_t|1\rangle\otimes
\bigotimes_{j>t}U_j|0\rangle,
\quad |p|=t-1.
```

Before position t, the unchanged marker bit makes the suffix predicate
false; from t onward the still unvisited suffix is zero. This proves
the displayed columns directly and proves that they form a unitary.

For a positive real target rotation $`R(\phi_j)`$, use the retained
[one-qubit approximation premise](../../docs/FAULT_TOLERANT_COMPILER.md#primary-references)
to choose an ancilla-free native word $`V_j`$ approximating
$`R(\phi_j/2)`$ to error $`e_j`$. A scalar phase on $`V_j`$ is
immaterial in the following definition; the usual literal approximation
also suffices. Set

```math
J_j=V_j Z V_j^\dagger,\qquad
\widehat U_j=J_jZ.
```

The actual $`J_j`$ is an exact Hermitian involution, the actual
$`\widehat U_j`$ is unitary, and

```math
R(\phi_j/2)ZR(-\phi_j/2)Z=R(\phi_j),\qquad
\|\widehat U_j-R(\phi_j)\|\le2e_j.
```

All scalar phases in $`V_j`$ cancel against its actual inverse.
The actual $`\widehat U_j`$ need not be real. The factorization
$`(VZV^\dagger)Z`$ is not a parametrization of arbitrary SU(2)
matrices; its restricted family suffices to approximate these real
rotations. The generalized frame geometry above needs only unitarity.

Let T now denote an s-bit tag register, H a preserved activity bit, and
$`X_1,\ldots,X_h`$ distinct data targets, disjoint from H and T.
Set $`s=\lceil\log_2(h+1)\rceil`$ and define masks

```math
f_j(H,T)=H[t\le j],\qquad 1\le j\le h.
```

Here t is the integer encoded by T. Tag zero selects every target;
padded tags above h select none. The
[complete shallow word query](SHALLOW_DIRTY_QUERIES.md#1-parallel-output-bits-with-a-read-only-address)
XORs these masks into arbitrary dirty bits $`d_j`$. Its address is
only H,T, its table has $`O(h)`$ rows and h outputs, and it preserves
the complete address. Call it $`Q_f`$. Conservative simultaneous
bounds are

```math
\begin{aligned}
T(Q_f),w(Q_f)&=O(h^{3/2}\log^3(h+2)),\\
G(Q_f)&=O(h^2+h^{3/2}\log^3(h+2)),\\
D_T(Q_f)&=O(\log\log(h+4)).
\end{aligned}
```

Outputs are included in the polynomial width reservation below.
Every output's non-Clifford supports and helper pool are private;
shared-address interactions are Clifford gates only.

For involutions $`A_j`$ on the distinct data targets, let $`C_A`$
apply $`A_j`$ controlled by $`d_j`$. Chronologically,

```math
Q_f,\ C_A,\ Q_f^\dagger,\ C_A
```

applies exactly $`A_j^{f_j}`$ and returns every dirty output and helper.
The dirty contribution squares away because $`A_j^2=I`$. This is a
full-input operator identity, including arbitrary reference systems.

Use this echo first with $`A_j=Z`$ and then, as a separate completed
echo, with $`A_j=J_j`$. The controlled-J bank has chronology

```math
V_j^\dagger,\ \mathrm{CZ}(d_j,X_j),\ V_j
```

on each target. All non-Clifford gates lie on the private data target
$`X_j`$, so the words for different j run in parallel. The result is
the exact selected product

```math
S_H(T)=\bigotimes_{j=1}^h\widehat U_j^{\,H[t\le j]}.
```

Chronological Z then J is essential: the resulting operator is JZ.
One cannot replace the two echoes by an echo of a non-involutory U.
Each completed echo is exactly identity on H zero, even with occupied
tag, masks, and helpers. Its dirty pool is reused after it returns.

Also define the unconditional product

```math
B=\bigotimes_{j=1}^h\widehat U_j.
```

It and its actual inverse have the maximum single-target T-depth and
the sum of their target counts. No controlled global B is needed.

## 4. Exact last-one XOR with arbitrary dirty work

For $`X\in\{0,1\}^h`$, let $`t(X)=0`$ on the zero string and
otherwise let t be the largest position with a one. Its one-hot tests are

```math
u_j(H,X)=H X_j\prod_{k>j}(1-X_k),\qquad 1\le j\le h.
```

Use the [read-only conjunction interface](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
to XOR these into distinct arbitrary dirty bits $`a_j`$, each with
private helper work. The public H,X wires occur only as Clifford CNOT
controls, including negative literals. Thus the h conjunctions run
their non-Clifford layers in parallel. The safe cubic bound per row in
the [physical-support proof](SHALLOW_DIRTY_QUERIES.md#physical-t-layer-parallelism)
gives a complete query $`Q_u`$ with

```math
T(Q_u),G(Q_u),w(Q_u)=O((h+2)^4),\qquad
D_T(Q_u)=O(\log_2(h+2)).
```

Let C CNOT $`a_j`$ into precisely the tag bits set in the binary
encoding of j. This costs $`O(h\log(h+2))`$ elementary Clifford
gates. The chronological word

```math
F_H:\qquad Q_u,\ C,\ Q_u^\dagger,\ C
```

is the exact complete-input XOR

```math
|H,X,T,d\rangle\longmapsto
|H,X,T\mathbin\oplus Ht(X),d\rangle.
```

Indeed the two appearances of each arbitrary initial $`a_j`$ cancel,
and the mutually exclusive tests XOR their binary labels to t. Every
helper returns. Tag T may be arbitrary; no clean work is used. This
completed word is identity on H zero and preserves H,X throughout its
complete interface.

## 5. One protected stage and exact tag return

For local size $`m\gt64`$ put

```math
r=\lceil\log_2(m+1)\rceil,\qquad h=m-r,\qquad
s=\lceil\log_2(h+1)\rceil\le r.
```

Split the local logical register into h upper bits X and r lower bits Z.
Save $`H=[Z=0^r]`$ in the sole external initialized qubit, using the
[returned two-helper conjunction](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates).
This costs $`O(r+1)`$ T and Clifford gates and $`O(\log(r+2))`$
T-depth, with two arbitrary helpers returned. Select s coordinates of
Z as tag T. The other coordinates are spectators. Save H before T is
occupied.

Execute the following chronological word and then the actual inverse
of the saved-predicate computation:

```math
F_H,\quad S_H(T),\quad B^\dagger,\quad F_H,\quad B.
```

On H one, T starts zero. A marker $`p1\,0^{h-t}`$ first receives
tag t. The selected product then creates

```math
|p\rangle\otimes\widehat U_t|1\rangle\otimes
\bigotimes_{j>t}\widehat U_j|0\rangle.
```

The actual $`B^\dagger`$ transforms this exactly into

```math
\left(\bigotimes_{j\lt t}\widehat U_j^\dagger|p_j\rangle\right)
\otimes |1\,0^{h-t}\rangle.
```

Every computational component still has last-one index t, regardless of
the rotated prefix. The second $`F_H`$ therefore clears T exactly.
B restores the desired column. For the root column t is zero, the
selected product equals B, and $`B^\dagger`$ returns X to zero before
the second XOR. This proves every active column with its literal phase.

On H zero, both completed $`F_H`$ words and $`S_H(T)`$ are identity
on the full space of their occupied inputs. The unconditional
$`B^\dagger`$ and B consequently cancel exactly. The inactive
argument assumes no zero tag or helper. Thus Z is unchanged on both
sectors, and the saved-predicate inverse returns H exactly. Linearity
gives the same conclusion on arbitrary logical, dirty, and reference
inputs. Exact return concerns the actual native words, not merely
ideal approximating rotations.

The implemented logical operator is exactly the first h layers of the
generalized m-bit frame: the h-bit generalized frame acts when Z is
zero, and the identity acts otherwise. The remaining r layers are
$`I_X`$ tensored with the corresponding homogeneous frame on Z,
because their parameters are independent of the upper prefix. Recurse
on Z, reusing the same returned external qubit and dirty pool. This
prefix independence is precisely where the structured promise enters.

For $`m\le64`$, use the direct layer circuit. For target j compute its
unchanged lower-suffix-zero predicate into H, apply controlled
$`\widehat U_j`$ by chronological controlled Z then controlled J,
and uncompute the predicate. The controlled-J word is again
$`V_j^\dagger,\mathrm{CZ}(H,X_j),V_j`$. An empty suffix uses
$`X_H`$ before and after the controlled word. Two arbitrary helpers
suffice; every control, dirty helper, and H returns exactly. There are
only an absolute constant number of base layers.

## 6. Precision, physical support, and the recursive ledger

Let $`m_0=n`$ and $`m_{i+1}=\lceil\log_2(m_i+1)\rceil`$ until
the base size is at most 64. At a nonbase stage i, its first
$`h_i=m_i-m_{i+1}`$ original rotations are assigned to that stage;
the base receives its remaining rotations. Every original depth is
assigned exactly once. Choose its one-qubit word with

```math
e_j\le\frac{\eta}{256(m_i+2)^2}.
```

For $`m_i\gt64`$ the next size is at most $`\sqrt{m_i}`$ and at
most $`m_i/2`$. The sum of reciprocals of nonbase sizes is bounded
geometrically, reading the sequence backwards; the single base term
satisfies $`m/(m+2)^2\le1/8`$. In particular the deliberately loose
uniform bound

```math
\sum_i\frac{m_i}{(m_i+2)^2}\le2
```

suffices. Telescoping the exact generalized layer circuit against the
positive real layer circuit gives the full operator bound

```math
\|W_{\widehat U}-W_+\|
\le2\sum_j e_j
\le\frac{\eta}{128}
\sum_i\frac{m_i}{(m_i+2)^2}
\le\eta/64.
```

This is not a columnwise error estimate. Each homogeneous layer has
operator error at most its single rotation error. The tag eraser is
exact for $`W_{\widehat U}`$, so no approximate spectral comparison
or source-leakage budget is hidden in this estimate. Together with
Section 2 the structured real branch has error less than eta.

The one-qubit word at stage size m has

```math
T,G,D_T=O(L+\log_2(m+2)).
```

The selected masks use private dirty query pools. The controlled-J and
B words use separate physical data targets for every non-Clifford gate.
The last-one conjunctions have private non-Clifford supports and share
only read-only Clifford controls. Serialize shared Clifford gates at
each T-layer boundary and charge their sum in G. Padding shorter native
words with identities synchronizes the T layers. This proves actual
T-depth parallelism without initialized control copies or a free fanout
operation.

There is a fixed constant C such that each stage has the safe bounds

```math
D_{T,\mathrm{stage}}=O(L+\log_2(m+2)),\qquad
w_{\mathrm{stage}}\le C(m+2)^4,
```

```math
T_{\mathrm{stage}},G_{\mathrm{stage}}
\le C(m+2)^4+Cm[L+\log_2(m+2)].
```

The width includes every dirty mask output, conjunction output, private
helper pool, and saved-predicate helper. Pools used by consecutive
completed macros are reused, and their finite sum also obeys the same
majorant. The tag is protected logical work and is not counted as an
external dirty or clean allocation. Ancilla-free one-qubit synthesis
adds no precision-dependent width. The base obeys these bounds too.

The square-root decrease above the base implies

```math
\sum_i\log_2(m_i+2)=O(\log_2(n+2)),\qquad
\sum_i(m_i+2)^4=O((n+2)^4),\qquad
\sum_i m_i=O(n).
```

The number of stages is at most a constant multiple of
$`\log_2(n+2)`$. Therefore, at fixed L, the complete recursive
magnitude compiler has

```math
D_T=O_L(\log_2(n+2)),\qquad w\le C'(n+2)^4,
```

```math
T,G=O\!\left((n+2)^4+n[L+\log_2(n+2)]\right).
```

No all-precision logarithmic-depth claim follows from this fixed-L
statement. The constants in the final theorem are allowed to depend on L.

## 7. Exact sign phases and every dirty width

Each of D and S from Section 2 is an arbitrary Boolean phase on the
original n logical bits. They require no initialized phase qubit.
Reserve one arbitrary dirty bit d for the query output and put
$`B_q=b-1`$ for its additional dirty query pool. At the real theorem's
literal threshold, $`B_q\ge4(n+2)`$. The complete
[routed query](ONE_CLEAN_DEPTH.md#routed-schedule-at-every-literal-allocation)
with one output, $`Q=N`$, and $`r=n`$ gives

```math
T_Q=O(1+\sqrt N+N/B_q),\qquad G_Q=O(N),\qquad
D_{T,Q}=O(1+N/B_q^2+n).
```

If $`B_q\ge16c_I(n+2)`$, select the
[shallow version](ONE_CLEAN_DEPTH.md#shallow-schedule-when-its-absolute-helper-constant-fits),
which replaces n by $`\log_2(n+2)`$ in depth with the same count.
If that condition fails, $`b\le C_I(n+2)`$ for an absolute constant,
and the bounded ratio $`n(n+2)^2/2^n`$ gives

```math
n=O(N/b^2).
```

Hence every eligible width supports a complete exact Boolean phase with

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O(N/b^2+\log_2(n+2)).
```

The physical phase word is chronologically
$`Q_f,Z_d,Q_f^\dagger,Z_d`$. It returns arbitrary d and all its
helpers and has phase $`(-1)^{f(x)}`$ on every input. The two signs
and the magnitude circuit run sequentially, so their dirty pools are
reused; their widths are not added.

To fit the magnitude compiler at every literal allocation, set

```math
b_c=\sqrt N/n.
```

There is an absolute dimension cutoff $`n_0`$ above which
$`C'(n+2)^4\le b_c`$. For $`n\ge n_0`$ and $`b\ge b_c`$ use
the exact finite sign gauges and the recursive magnitude circuit.
Its polynomial counts are absorbed by $`O_L(\sqrt N)`$ T and
$`O_L(N)`$ Clifford gates, since every fixed polynomial divided by
the respective exponential is uniformly bounded. Sections 2 and 6
supply the prescribed error. Section 7's sign words supply the width
term in the displayed theorem.

For $`b\lt b_c`$, select the existing
[R79 real-frame circuit](UNIFORM_SUBTREE_DEPTH.md#1-the-all-precision-theorem)
on the original target at the original L. It already fits exactly
$`17(L+n+7)`$ and has the required count and error. At fixed L its
additional depth is $`O_L(\log^2(n+2))`$, whereas
$`N/b^2\gt n^2`$; the additional term is therefore absorbed by the
width term. For the finitely many dimensions $`n\lt n_0`$, use R79
at every width. Its additive depth is bounded by an absolute multiple
of $`\log(n+2)`$ on that finite set, with constants allowed to depend
on L. These branches preserve the literal reservation without assuming
that an unknown helper constant fits its coefficient 17.

This proves the real theorem with one selected circuit supplying every
bound simultaneously. The recursive structured branch returns its
external clean qubit exactly; the fallback uses R79's full initialized-
isometry error contract. No reset or postselection joins the branches.

## 8. Exact signed balanced corollary

Suppose every node has
$`|c_v|=|s_v|=1/\sqrt2`$, with arbitrary signs. Then, for every
integer $`b\ge17(n+13)`$, one external initialized qubit suffices for
the exact initialized-isometry identity

```math
VJ_1=J_1(W\otimes I_b),
```

and simultaneous bounds with absolute constants

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O(N/b^2+\log_2(n+2)).
```

The balanced coefficients have a fixed positive gap from zero, so their
signs are determined at fixed certified resolution. No rounding is
needed. For the positive magnitude take $`\phi=\pi/4`$. Let
$`C=SH_{\rm Had}`$, so $`CZC^\dagger=Y`$, and choose the exact
native word $`V=CTC^\dagger`$. Since

```math
T=e^{i\pi/8}e^{-i\pi Z/8},\qquad
V=e^{i\pi/8}e^{-i\pi Y/8}=e^{i\pi/8}R(\pi/8),
```

the actual word $`(VZV^\dagger)Z`$ equals $`R(\pi/4)`$ exactly.
The scalar phase cancels, not merely up to global phase. Thus every
selected echo, product conjugator, tag return, and sign phase in the
structured branch is exact, with constant-length one-qubit words.

For the narrow-width and finite-dimension branches, use an explicit
exact fallback rather than an approximate generic compiler. Execute all
n homogeneous layers directly. Each unchanged suffix predicate uses
the returned two-helper conjunction, followed by the exact controlled
Z/J words and the actual predicate inverse. This uses one initialized
qubit, two arbitrary helpers, $`O(n^2)`$ T and Clifford gates, and
$`O(n\log(n+2))`$ T-depth. Append the exact D and S sign words with
the all-width ledger of Section 7, reusing the pool. For $`b\lt b_c`$,
the magnitude depth is absorbed by $`N/b^2\gt n^2`$. For the finite
set $`n\lt n_0`$ it is absorbed by $`O(\log(n+2))`$ after enlarging
an absolute constant. The polynomial counts are absorbed by the
exponential bounds. The literal threshold $`17(n+13)`$ supplies all
of these constant helper reservations and the sign-query routed branch.

This proves an exact signed balanced compiler on every stated width.
Its square-root T-count has the same order as the general-family upper
bound; no matching lower bound for the balanced subclass is asserted.

## 9. Prescribed complex phases and verification boundary

For $`D_\phi W_{\mathbb R}`$, run the real theorem at precision
$`L+1`$ and the retained
[one-clean literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at that same precision. Each fits $`b\ge17(L+n+8)`$ and has error
at most $`2^{-(L+1)}`$. At fixed L the diagonal has

```math
T_D=O_L(\sqrt N+N/b),\qquad G_D=O_L(N),\qquad
D_{T,D}=O_L(N/b^2+\log_2(n+2)).
```

Reuse the same sole external signal and dirty pool. A two-factor unitary
hybrid bounds the total initialized-isometry error by $`2^{-L}`$.
Actual intermediate leakage is propagated by the next actual unitary;
only the comparison term uses ideal initialized work. Every leaf phase,
including at a zero-amplitude leaf, is retained. Arbitrary leaf phases
do not inherit the exactness of the balanced real corollary.

The bounded checks in
[`tests/test_homogeneous_depth.py`](../../tests/test_homogeneous_depth.py)
compare independent generalized-frame columns with the complete stage
word, including nonreal one-qubit factors, every tag sector, inactive
occupied tags, arbitrary dirty masks, and literal signs. Separate native
echo checks use actual complex V and the restricted $`(VZV^\dagger)Z`$
family. Exact rational sign fixtures include zero edges and markers
below zero incoming branches. These checks support the algebraic
identities; the resource bounds follow from the physical-support and
all-width arguments above, not from finite matrix agreement.

The generic problem remains open at the same frontier: an arbitrary
node-dependent rotation table has no established common shallow product
conjugator, and its terminal parameters retain upper-prefix dependence.
Neither generic threshold-query batching nor a supplied diagonalizing
unitary is assumed by this restricted-family proof.
