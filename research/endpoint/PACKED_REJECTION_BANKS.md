# Packed rejection banks and one final recovery

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Joint native words](JOINT_NATIVE_WORDS.md) · [One-clean source](../../docs/ONE_CLEAN_COMPILER.md)

Supported block encodings can store rejected amplitudes in disjoint logical
mode banks instead of an initialized history counter. For the Hopf tree,
those banks fit in the two supplied flags. A diagonal preweight converts
the accepted product into half of the complete target, allowing one final
coherent recovery. Balancing the three equal-support stages at each depth
also gives an exact four-port compressor and a raw-to-canonical comparison
within the literal dirty width. The direct source implementation still pays fine
precision at every depth. This does not improve the endpoint bound.

## 1. Exact multiplication with disjoint rejection spaces

Let the success space S have dimension N, and tensor every space below
with the same arbitrary dirty input space. Let P_j be projectors on S.
Allocate orthogonal failure spaces F_j of dimension $`\mathrm{rank}P_j`$
and one preweight space F_0 of dimension N. Assume

```math
N+\dim F_0+\sum_j\dim F_j\le4N.
```

These are modes of n logical qubits and two flags, not separately
initialized registers. Let K_j be a complete unitary supported on
$`P_jS\oplus F_j`$, extended by identity on the other modes. Its success
compression is A_j; other blocks may act on arbitrary dirty inputs.
If J inserts the success flags, then

```math
J^\dagger K_R\cdots K_1J=A_R\cdots A_1.
```

After stage j, each rejected component is in an earlier failure bank.
Every later stage fixes that bank, so it cannot return to success.
There is no measurement, discarded amplitude, or fresh-zero assumption.
Inverses reverse the actual complete circuit.

The support assumption is essential. Identity only on flag-zero inactive
inputs is insufficient: earlier failure banks can reach flag-one inactive
inputs after routing.

## 2. Column attenuation and exact preweight

Suppose P_j are nested increasing and V_j is unitary supported in P_j.
For positive $`\alpha_j\le1`$, write

```math
A_j=(I-P_j)+\alpha_jV_jP_j=V_j\Delta_j,\qquad
\Delta_j=(I-P_j)+\alpha_jP_j.
```

Every later Delta commutes with each earlier V, because its projector
is identity on the earlier support. Therefore

```math
A_R\cdots A_1=WD,\qquad
W=V_R\cdots V_1,\qquad D=\Delta_R\cdots\Delta_1.
```

D is generally not $`\prod_j\alpha_j`$ times identity: input columns enter
the nested supports at different depths. Treating it as a scalar gives
the wrong full frame.

In their common support basis, D and its preweight are diagonal; for the
Hopf projectors this is the computational basis. If $`D\succeq I/2`$,
the contraction $`E=(2D)^{-1}`$ has entries
in [1/2,1]. Apply its one-signal dilation first, rejecting into F_0.
All subsequent stages fix F_0. The complete word U then satisfies

```math
J^\dagger UJ=WDE=W/2.
```

E need not commute with W; it is applied before W.

## 3. Physical banks and charged routing

Write the two flags as (f,e). Success is (0,0,x); F_0 is (0,1,x).
The active depth-d Hopf labels are

```math
x=p_d\,b\,0^{n-d-1},\qquad \mathrm{rank}P_d=2^{d+1}.
```

For one axis per depth, reserve the failure bank

```math
(f,e,x)=(1,b,0^{n-d-1}1p_d).
```

The sentinel makes the banks disjoint. Their total dimension is
$`\sum_d2^{d+1}=2N-2`$, within the 2N available f=1 modes.

The canonical primitive uses f as signal and has predicate e=0 together
with the logical zero suffix. On f=1, permute the other n+1 wires:
send the prefix into the low d logical positions, b into e, and old e
into the sentinel position; flip the latter bit. The frozen suffix fills
the other positions. On f=0 the map is identity. This specifies an affine
permutation on the entire basis, not only on initialized inputs.

A wire permutation uses O(n) swaps. Each f-controlled swap is two CNOTs
around one Toffoli; the sentinel flip is a CNOT. Each routing map and
inverse therefore cost O(n) T gates and no clean work. Its conjugated
primitive is supported precisely on success P_d and its assigned bank.
Other banks map outside its full-input predicate and remain fixed.
Total routing cost is O(n^2).

### Three axes for an actual complex coarse word

An exact native coarse SU(2) coin need not be real, so its residual can
require three axis rotations. Reserve the last two depths for the existing
standalone compiler. For $`0\le d\le n-3`$ and k=0,1,2, use

```math
(f,e,x)=(1,b,\mathrm{type}_k\,0^{n-d-3}1p_d),
```

where type is a distinct two-bit string. These banks use
$`3\sum_{d=0}^{n-3}2^{d+1}=3N/2-6\lt 2N`$ modes. The same wire permutation
and conditional flips place the type bits and sentinel. The three stages
at a depth share P_d, preserving nesting. Coarse words act on success
only and fix all failure banks.

## 4. The retained native source supplies near-unit acceptance

Use the raw word $`Q=A^\dagger B_gA`$ from the
[paired-source construction](../../docs/ONE_CLEAN_COMPILER.md#3-conjugating-scalar-blocks-produces-a-rotation),
without separate recovery. Fix the head sign positive. For tail means u,v,
its accepted block is $`aI+bXZ`$, with

```math
a=\frac14+\frac u2+\frac v4,\qquad
b=-\frac18-\frac u4+\frac{3v}{8}.
```

For $`a=\alpha\cos\theta`$, $`b=\alpha\sin\theta`$, choose

```math
1-u=\frac32(1-a)+b,\qquad 1-v=(1-a)-2b.
```

If $`|\theta|\le\delta\le1/100`$ and $`\alpha=1-4\delta`$, then
$`1-a\ge4\delta`$, $`|b|\le\delta`$, and both right sides lie in [0,2].
Thus both means are legal. Certified geometric rounding gives accepted
error below $`(5/2)2^{-q}`$, with exact dirty identity in that block.
Clifford conjugations give the other axes. Rejected components retain the
actual source action.

The [inactive-sector construction](../../docs/ONE_CLEAN_COMPILER.md)
is identity on every inactive input, including the rejected signal.
Adding e=0 to its predicate supplies Section 1's support assumption.
The preweight uses the [scalar source](../../docs/OPERATOR_SOURCE_COMPILER.md),
controlled by f=0 and with e as signal. Its logical table has at most N
rows, its scalar error is below $`3\,2^{-q}`$, and it fixes all f=1 banks.

### Coarse accuracy must be weighted by depth

Choose $`\delta_d=1/[100(n-d)^2]`$. Select an actual native coarse coin C
close enough to target L that $`C^\dagger L`$ has three axis angles at most
delta_d. The derivative at identity of
$`(x,y,z)\mapsto R_x(x)R_y(y)R_z(z)`$ is invertible; a fixed local inverse
bound suffices at every depth. Certified arithmetic within this chart
requires no equality oracle, including for zero residuals.
Apply the residual factors and then C, giving L exactly in the ideal word.

The [exactly returned coarse-word construction](ORIGINAL_TREE_RESIDUAL.md#3-an-exactly-logical-native-coarse-word)
permits word length $`O(\log(1/\delta_d))`$ per address. The
[literal reflection alphabet](../../docs/BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
is $`\{T^jHT^{-j}:0\le j\lt 8\}\cup\{Z\}`$. Its dirty echo
can condition every involution G on success and suffix without a clean
buffer. Query a table bit f(x) into arbitrary z, apply the
predicate-and-z-controlled G, unquery, and apply that controlled G again.
Since $`G^2=I`$, the net exponent is exactly the predicate times f(x),
on both z inputs. Each symbol position costs $`O(2^d+n^2)`$ T gates
by the retained whole-word query and controlled-Pauli construction.
All helpers return exactly, including on failure banks. Consequently

```math
\sum_d2^d\log(1/\delta_d)=O(N),\qquad
\prod_d(1-4\delta_d)^3
\ge1-12\sum_d\delta_d>1-\pi^2/50>1/2.
```

D meets the preweight condition. All additional polynomial-in-n costs
are O(N). Uniform inverse-polynomial coarse accuracy at every address
would not give the same linear total by this word-length accounting.

## 5. One final recovery controls the entire isometry

For the rounded complete word let $`B=J^\dagger UJ`$ and
$`\epsilon=\|B-W/2\|\le1/4`$. With $`R=I-2JJ^\dagger`$, use

```math
\mathcal A=-URU^\dagger RU,\qquad
J^\dagger\mathcal AJ=3B-4BB^\dagger B.
```

R tests only the two flags and is Clifford. The literal leading minus
is the Clifford scalar $`-I=(HS)^{12}`$ on a retained wire.
The inverse includes the actual action on all banks.

Write $`B=VH`$ for its polar decomposition. Then
$`\|H-I/2\|\le\epsilon`$ and $`\|V-W\|\le4\epsilon`$. For
$`p_3(x)=3x-4x^3`$,

```math
p_3(1/2+t)=1-6t^2-4t^3,\qquad
\|\mathcal AJ-JV\|^2=2\|I-p_3(H)\|
\le(12+8\epsilon)\epsilon^2.
```

This is a full-output identity, including leakage. Hence
$`\|\mathcal AJ-JW\|\le8\epsilon`$, strictly for positive epsilon.
Norms include every dirty input and
remain unchanged upon tensoring reference identities.

The bank identity holds for rounded primitives as well. All their
accepted blocks are contractions, so telescoping bounds epsilon by
the sum of their block errors, including preweight error.

## 6. Charged cost and the surviving question

There are fewer than 3n raw stages and one preweight. At common q-bit
precision, the direct construction has

```math
T=O(N+nq),\qquad \epsilon\le3(3n+1)2^{-q}.
```

This includes O(N) table rows, coarse work, routing/predicates, every
O(q) source and inverse, and all three complete recovery calls.
The last two standalone depths add O(N+q). At the endpoint choose
$`q=N+\lceil\log_2(3n+1)\rceil+8`$. Recovery error is below
$`3\eta/32`$; give each standalone depth error at most eta/4.

For $`n\ge16`$, q+1 core wires plus two outside dirty helpers fit in
$`N+n+7`$. The core supports self-borrowed whole-word queries for all
at-most-n-bit addresses. Coarse queries temporarily borrow selectors from
that returned core; the two outside helpers hold the queried table bit
and the borrowed multiple-control helper. All are arbitrary inputs and
return exactly after the coarse word. No intermediate
initialization is assumed. Finitely many smaller n use the retained
compiler with a fixed larger constant. The endpoint cost is O(nN),
worse than the existing grouped upper bound.

The useful reduction is that a clean history counter and separate group
recoveries are unnecessary for these supported stages. The unresolved
step is joint native synthesis of this packed word at O(N), including
preweight and the actual inverse. Hoisting a common source loader
conjugates the intervening masks; those transformed masks are not free.
Lower bounds for separately accurate groups do not apply to this joint
word.

## 7. Finite evidence

The [bounded diagnostics](../../compiler_robust_hopf/packed_rejection_banks.py)
and [tests](../../tests/test_packed_rejection_banks.py) check complete routing
permutations, inactive banks, multiplication, the failure of scalar-only
normalization, preweight, and actual inverse recovery.
They are finite matrix/permutation checks, not an arbitrary-angle native
emitter. Dimension, cost, and full-isometry statements above are analytic.

## 8. Balanced scalar histories with two preweight banks

There is a stronger scalar-history identity if the attenuations are chosen
together with the preweight. It uses **ideal canonical dilations** and a
different allocation from Sections 1–6: two N-dimensional preweight banks
instead of one. It does not identify the raw source's rejected dirty action
with a scalar for free.

Let $`P_1\preceq\cdots\preceq P_R`$ be nested coordinate projectors,
and let V_k be unitary supported in P_k. Identify F_k with P_k S and use
the complete canonical unitary

```math
K_k\big|_{P_kS\oplus F_k}
=\begin{pmatrix}\alpha_kV_k&-s_kI\\s_kI&\alpha_kV_k^\dagger\end{pmatrix},
\qquad s_k=\sqrt{1-\alpha_k^2}.
```

K_k is identity on all other modes. The matrix is unitary for every V_k,
including unequal complex coins. The algebra below tensors with identity
on every arbitrary dirty input and reference. In this section require

```math
3N+\sum_{k=1}^R\mathrm{rank}P_k\le4N.
```

Write B for the baseline-preweight bank and C for the cap-preweight bank,
both copies of S. Fix $`\lambda=1/64`$ and set

```math
\begin{aligned}
g_k&=\frac14+\lambda(1-2^{k-R}),&&0\le k\le R,\\
w_k^2&=g_{k-1}-g_k=\lambda2^{k-R-1},&&1\le k\le R,\\
\alpha_k^2&=\frac{g_k}{g_{k-1}}.
\end{aligned}
```

For a coordinate x, let e(x) be its first active stage, with e(x)=R+1
when it belongs to no support. Prepare the following column of a full
unitary E:

```math
E|S,x\rangle
=\sqrt{g_{e(x)-1}}\,|S,x\rangle
+\frac{\sqrt{47}}8|B,x\rangle
+\sqrt{\lambda2^{e(x)-R-1}}\,|C,x\rangle.
\qquad\text{(9)}
```

The squared amplitudes sum to one, including the no-entry case. This is
the same success preweight as Section 2, since

```math
D_x=\prod_{k=e(x)}^R\alpha_k
=\frac{1}{2\sqrt{g_{e(x)-1}}},
\qquad (2D_x)^{-1}=\sqrt{g_{e(x)-1}}.
```

In particular $`D_x\ge4/\sqrt{17}`$, and every attenuation is at least
$`\sqrt{32/33}`$. No small accepted amplitude is introduced.

A complete E is specified by two rotations. First rotate S into S and B
with cosine $`\sqrt{17}/8`$ and sine $`\sqrt{47}/8`$. Then rotate S into
S and C, addressed by x, with

```math
\cos\phi_x=\frac{\sqrt{g_{e(x)-1}}}{\sqrt{17}/8},
\qquad
\sin\phi_x=\frac{\sqrt{\lambda2^{e(x)-R-1}}}{\sqrt{17}/8}.
\qquad\text{(10)}
```

Each rotation has the usual real block on its two copies and fixes the others:

```math
\begin{pmatrix}c&-s\\s&c\end{pmatrix}.
```

 These equations specify the action
on initially occupied and unoccupied flag sectors; E is not merely an
isometry completion left to an oracle. The entry function is a finite
table determined by the support labels, requiring no angle equality test.

Put $`W_k=V_k\cdots V_1`$, with W_0=I. Immediately before stage k, the
success contribution of an active input column is
$`\sqrt{g_{k-1}}W_{k-1}|x\rangle`$. Thus its rejected contribution in
F_k is $`w_kW_{k-1}|x\rangle`$. Define the complete bank unitary

```math
\mathcal B=I_{S\oplus B\oplus C}
\oplus\bigoplus_{k=1}^R W_{k-1}^\dagger\big|_{P_k}
\oplus I_{\rm unused}.
\qquad\text{(11)}
```

The restriction exists because earlier supported unitaries preserve every
later support. After applying this **unwinding bank**, each active
F_k column has scalar amplitude w_k in the original x coordinate. Its
native implementation has not been supplied by (11).

For k=1 through R, apply an exact Hadamard between C(x) and F_k(x) for
every x in P_k, with the positive sum going to C. Denote the chronological
product by H_bal. A column entering at e starts with cap amplitude w_e.
Its first Hadamard acts on $`(w_e,w_e)`$, produces
$`(\sqrt2w_e,0)`$, and for e less than R,
$`\sqrt2w_e=w_{e+1}`$. At the last stage,
$`\sqrt2w_R=\sqrt\lambda=1/8`$. Induction therefore zeroes every
stage bank and leaves that cap amplitude.
A no-entry column already has that amplitude. Consequently

```math
H_{\rm bal}\mathcal B K_R\cdots K_1 E J
=\begin{pmatrix}
W_R/2\\(\sqrt{47}/8)I\\I/8\\0
\end{pmatrix}.
\qquad\text{(12)}
```

The last row includes every stage bank and unused mode. This is a complete
initialized-isometry equality, with literal positive coefficients. One
uniform ordered Hadamard network handles all entry stages; it does not
look up e(x) during compression. Omitting (11) generally invalidates (12),
even though the accepted block remains W_R/2. Section 5's recovery still
applies to that accepted block and the actual uncompressed rejections.

## 9. Physical allocation and exact native scalar compression

For n at least four, reserve the last **three** Hopf depths for the
standalone compiler. Use three stages at each depth
$`0\le d\le n-4`$, so $`R=3(n-3)`$. Allocate flags (f,e)=00 to success,
01 to B, 10 to C, and 11 to the stage banks. At depth d and axis j=0,1,2,
the bank coordinate corresponding to $`x=p_d b0^{n-d-1}`$ is

```math
(f,e,y)=(1,1,\mathrm{type}_j\,0^{n-d-4}1p_db).
\qquad\text{(13)}
```

The two-bit type and the sentinel make all banks disjoint. They occupy

```math
\sum_k\mathrm{rank}P_k
=3\sum_{d=0}^{n-4}2^{d+1}=3N/4-6\lt N.
```

Thus the three full copies and all stage banks fit within the two flags'
4N modes. There is no initialized history register.

Both a canonical stage and its cap-to-bank Hadamard have full-input affine
routing maps. Here is one explicit construction. Fix c=0 for a stage or
c=1 for a compressor. On e=1, permute the remaining n+1 wires as follows:
send b to logical bit zero, the d prefix bits to the next d logical bits,
old f to the sentinel position, and one frozen suffix bit to f. Put the
remaining frozen bits in the remaining logical positions. Flip f, flip
the sentinel by $`1-c`$, and insert the two type bits by fixed flips.
On e=0 act as identity. This is a controlled affine permutation on the
entire basis, and it maps

```math
(c,1,p_db0^{n-d-1})\longmapsto
(1,1,\mathrm{type}_j0^{n-d-4}1p_db).
```

For c=0, conjugating an e-signal canonical stage with predicate
$`f=0`$ and zero suffix routes exactly its rejected sector into F_k.
For c=1, conjugating H on e with predicate $`f=1`$ and zero suffix
implements the required Hadamards between C(x) and F_k(x). The whole
e=0 sector is fixed by each route. Conjugated primitives are identity
outside their specified pair of subspaces, including on all other banks.
The inverse is the inverse of the complete affine permutation.

There are O(n) controlled wire swaps and conditional fixed flips per
route, costing O(n) T gates. A predicate-controlled H is exactly native:
the fixed Clifford+T word $`V=SHTHS^\dagger`$ satisfies
$`VZV^\dagger=H`$, so conjugate a multiply controlled Z by V. The
[borrowed predicate construction](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
implements its equivalent multiply controlled X with one returned dirty
helper and O(n^2) exact Toffolis. Therefore all of H_bal, its routes, and
its inverse cost O(n^3) T gates, independently of q. No rotation angle,
clean helper, or approximation is hidden in this scalar compression.
This is a native decomposition and cost argument, not an emitted scalable
gate word.

The two preweight rotations (10) are not claimed to be exactly native.
The [borrowed-axis compiler](BORROWED_AXIS_COMPILER.md) provides their
separate full-norm approximation with O(N+q) table/source cost, plus
polynomial predicate overhead, charged once. The other flag is a sector
predicate in each rotation and need not be supplied in zero. Its arbitrary
occupied sectors are included in the full-unitary contract.

The balanced attenuation also allows linear total coarse-word length.
For $`\delta_k=(1-\alpha_k)/4`$,

```math
\delta_k=\frac{w_k^2}{4g_{k-1}(1+\alpha_k)},
\qquad
\frac8{17}w_k^2\le\delta_k\le w_k^2\le1/128.
```

Thus the ideal near-unit source chart in Section 4 is available for
residual angles at most delta_k. At each depth choose coarse accuracy for
the smallest of its three deltas. Its word length is
$`O(R-k+1)=O(n-d)`$, giving
$`\sum_d2^d(n-d)=O(N)`$. Predicate overhead is polynomial in n and
hence O(N) asymptotically. This prices the coarse words; it does not turn
their raw rejection operators into the canonical blocks of Section 8.

## 10. The exact full-unitary prefix factorization

The operation missing from (12) can be isolated on the complete space,
including inputs initially in failure banks. Let R_k be the scalar version
of K_k, obtained by replacing V_k by I. Then

```math
K_R\cdots K_1
=\mathrm{diag}\!\left(W_R,\bigoplus_{k=1}^R
 W_{k-1}\big|_{P_k}\right)
 (R_R\cdots R_1)
 \mathrm{diag}\!\left(I,\bigoplus_{k=1}^R
 W_k^\dagger\big|_{P_k}\right).
\qquad\text{(14)}
```

The display orders success before the stage banks; B, C, and unused modes
carry identities throughout. To prove it, factor each K_k as R_k between
a left boundary with W_k on success and W_{k-1} on F_k, and a right
boundary with $`W_{k-1}^\dagger`$ on success and $`W_k^\dagger`$ on F_k.
Its four active blocks are exactly
$`\alpha_kV_k,-s_kI,s_kI,\alpha_kV_k^\dagger`$; outside P_k the success
factors cancel. Adjacent success boundaries cancel in the product. A
boundary on F_k commutes with every R_j for j different from k, since R_j
fixes F_k pointwise. Move those bank factors to the two exterior
boundaries to obtain (14).

The right boundary is identity on the preweight columns (9), which have
no F_k component. The left bank boundary is precisely the inverse of
(11). Thus (14) explains both the scalar history and the required prefix
unwinding on the full unitary. It is not a synthesis shortcut: treating
W_R or either prefix bank as an already available native gate would assume
the unresolved work.

## 11. Precision scope and checked consequences

Sections 8–10 prove an ideal canonical identity and an exact native scalar
compressor. An addressed canonical stage with two-dimensional target rows
has a separate full-norm compiler at O(S+q) cost plus predicate overhead:
factor its local matrix as
$`\mathrm{diag}(V_k,I)R_k\mathrm{diag}(I,V_k^\dagger)`$
and use the borrowed-axis compiler on its bounded number of addressed
axis rotations. Together with (10), this gives the usual separately
compiled canonical-cascade baseline. It still pays fine precision for
each stage and does not provide uniform O(N+q) synthesis of (11) or (14).

The raw paired source has its actual dirty rejection orientation. The
[borrowed-axis analysis](BORROWED_AXIS_COMPILER.md) discusses how that
orientation can be represented and separately priced. Applying a coarse
coin only on success can also change the two rejection-side coordinates.
Neither a bank-coordinate choice nor those identities make the raw word
equal to the particular canonical completion used here at no cost.

There is an additional precision issue for this **serial** balanced schedule.
The earliest canonical defect satisfies $`s_1=\Theta(2^{-R/2})`$.
If a rounded accepted coin differs from $`\alpha V`$ by epsilon and has
scalar radius $`\widehat\alpha`$, then

```math
|\widehat\alpha-\alpha|\le\epsilon,\qquad
|\widehat s-s|
\le\frac{(2\alpha+\epsilon)\epsilon}{s},
\qquad
\widehat s=\sqrt{1-\widehat\alpha^2}.
```

Consequently a direct stagewise comparison of rounded canonical blocks
with the ideal balanced blocks has bound
$`O(R2^{R/2}2^{-q})`$. To make this bound at most $`2^{-N}`$, a sufficient
choice is $`q=N+R/2+O(\log R)`$. This estimate exceeds the literal source
width allocation used in Section 6 for growing n. It is not a lower
bound: compiling the ideal canonical rotations separately remains valid,
and a different joint construction might avoid this comparison. However,
the usual $`q=N+O(\log n)`$ accepted-block estimate does not establish
the same full-norm comparison or exact Hadamard cleanup. Canonicalizing
the actual rounded defect preserves that defect rather than the balanced
coefficients. The actual rounded radii can also depend on the address;
those contractions need not commute with earlier logical coins. The
scalar-history identities here require the ideal stagewise constant
alpha_k. For the physical raw word, the accepted-block telescoping bound
and recovery of Section 5 remain the applicable guarantees. No O(N)
endpoint compiler follows from these sections.

The [balanced-history tests](../../tests/test_balanced_rejection_history.py)
check exact rational attenuation and preweight identities for every entry
stage, integer bank capacity and complete routing permutations, and a
bounded unequal-target canonical cascade with full unitarity, preweight,
unwinding, and Hadamard compression. They also check that omitting
unwinding spoils (12) while leaving ordinary recovery valid.
These checks support the stated algebra; they do not emit the unresolved
native prefix bank. The depthwise schedule below removes the stated width
overrun without supplying that bank.

## 12. Depthwise balance fits the literal precision width

Three consecutive residual-axis stages at one Hopf depth have the same
support. Their three rejection amplitudes can therefore be balanced
together with the cap. This replaces three successive two-port merges
by one four-port Walsh transform. Unlike the serial schedule, its full
rounded-source comparison fits the endpoint width.

Let D be the number of retained depths, indexed by $`0\le d\lt D`$.
Keep the two full preweight banks B,C and the three stage banks per depth
from Section 9. Set

```math
c_d=\frac{2^{d-D}}8,\qquad
g_{d,a}=\frac{17}{64}-(a+1)c_d^2,\quad 0\le a\le3,
\qquad
\alpha_{d,a}^2=\frac{g_{d,a+1}}{g_{d,a}},\quad 0\le a\lt 3.
\qquad\text{(15)}
```

In particular $`g_{d,3}=g_{d+1,0}`$ and $`g_{D,0}=1/4`$.
For an input coordinate whose first active depth is e, prepare

```math
E|S,x\rangle
=\sqrt{g_{e,0}}|S,x\rangle
+\frac{\sqrt{47}}8|B,x\rangle+c_e|C,x\rangle.
\qquad\text{(16)}
```

Use e=D for no entry. The squared amplitudes sum to one in every case.
As in (10), two complete flag rotations specify E: the first has
cosine $`\sqrt{17}/8`$ and sine $`\sqrt{47}/8`$; the second has
cosine $`8\sqrt{g_{e,0}}/\sqrt{17}`$ and sine
$`8c_e/\sqrt{17}`$. Both act on their whole pair of flag sectors
and fix the other sectors. No completion oracle is required.

Use the ideal canonical stage of Section 8 with attenuation
$`\alpha_{d,a}`$ and its actual residual-axis unitary. Immediately
before that stage, the active success amplitude has scalar modulus
$`\sqrt{g_{d,a}}`$. Its rejection has modulus exactly c_d, since

```math
g_{d,a}(1-\alpha_{d,a}^2)=c_d^2.
\qquad\text{(17)}
```

The logical factor is the complete chronological prefix before the stage.
It includes any intervening success-only coarse coins. These coarse
operations preserve the nested supports and fix every failure bank, so
they change that prefix but not (17). Keep them as separate operations:
multiplying a canonical stage on success alone does not give another
canonical stage with scalar off-diagonal blocks.

Let $`\mathcal B_D`$ be the direct sum of the actual ideal prefix
inverses on the private stage banks, identity on S,B,C and unused modes.
After applying it, the cap and the three failure coordinates at the
first active depth have amplitudes $`(c_d,c_d,c_d,c_d)`$. Apply
$`H\otimes H`$ in the ordered four-port coordinates
$`(C,F_{d,0},F_{d,1},F_{d,2})`$, with the all-positive row going to C.
This leaves $`(2c_d,0,0,0)`$. Since $`2c_d=c_{d+1}`$,
the same induction works for every entry depth, including no entry.
Writing $`H_D`$ for these chronological four-port transforms and U_can
for the ideal canonical cascade, with its separate coarse operations,
gives the full initialized-isometry identity

```math
H_D\mathcal B_D U_{\rm can}EJ
=\begin{pmatrix}
W/2\\(\sqrt{47}/8)I\\I/8\\0
\end{pmatrix}.
\qquad\text{(18)}
```

W is the complete retained-depth logical word. Identities on every
arbitrary dirty input and reference are understood throughout.

### Complete native Walsh routing and coarse cost

The four-port transform is a specified full unitary, including occupied
failure inputs. It is the product of the two disjoint pairs
$`H(C,F_{d,0})`$, $`H(F_{d,1},F_{d,2})`$, followed by
$`H(C,F_{d,1})`$, $`H(F_{d,0},F_{d,2})`$. The cap-bank pair-H
operations use Section 9's full controlled-affine routes. A bank-bank
pair-H can be implemented by swapping its first bank with the active
cap coordinates, applying the routed cap-bank H for its second bank,
and undoing the swap. Use the same routed pair gate with X instead of
H for the swap. All these gates fix every other mode exactly.

Thus the complete H_D and its actual inverse cost $`O(n^3)`$ T gates,
use one returned arbitrary predicate helper, and need no clean register
or approximated rotation. The bank count is

```math
3\sum_{d=0}^{D-1}2^{d+1}=6(2^D-1).
\qquad\text{(19)}
```

Together with S,B,C it fits 4N modes whenever $`D\le n-3`$.
For the precision ledger below choose $`D=n-8`$, $`n\ge9`$;
the last eight depths use the existing standalone compiler.

Every attenuation obeys $`\alpha_{d,a}\ge\sqrt{64/65}`$.
With $`\delta_{d,a}=(1-\alpha_{d,a})/4`$,

```math
\delta_{d,a}
=\frac{c_d^2}{4g_{d,a}(1+\alpha_{d,a})}
=\Theta(c_d^2),\qquad \delta_{d,a}\lt 1/100.
```

The near-unit raw-source chart of Section 4 is therefore available.
Coarse word length at depth d is $`O(D-d+1)`$, and
$`\sum_d2^d(D-d+1)=O(2^D)`$. The coarse words, their predicates,
and the four-port routing all have total O(N) T cost.

### Rounded canonical comparison without a width overrun

Let $`s_{d,a}=\sqrt{1-\alpha_{d,a}^2}=c_d/\sqrt{g_{d,a}}`$.
For each programmed row the native paired source gives accepted error
$`\epsilon\lt (5/2)2^{-q}`$. Its exact rejection-orientation gauge from
[the borrowed-axis theorem](BORROWED_AXIS_COMPILER.md#4-one-global-canonicalization-with-its-scope)
puts its entire raw block in canonical form. The rounded radius and
defect may depend on the row; do not assume that their histories are
exactly scalar. Instead compare complete blocks before any unwinding.
The bound (8) in that theorem is
$`\epsilon+2\epsilon/s_{d,a}`$ on the full active block.

Since $`\sqrt{g_{d,a}}\lt \sqrt{17}/8`$,

```math
\sum_{d=0}^{D-1}\sum_{a=0}^2\frac1{s_{d,a}}
\lt 6\sqrt{17}(2^D-1).
```

Full-unitary telescoping, including the unchanged coarse operations,
therefore bounds the complete cascade difference by

```math
\Delta_q
\lt \left(\frac{15D}{2}+30\sqrt{17}(2^D-1)\right)2^{-q}
\lt 128\,2^{D-q}.
\qquad\text{(20)}
```

The last inequality uses $`D/2^D\le1/2`$ and
$`15/4+30\sqrt{17}\lt 128`$. There is no additional factor D:
the reciprocal defects form a geometric sum.
At $`\eta=2^{-N}`$, choose

```math
q=N+D+12=N+n+4.
```

Then $`\Delta_q\lt \eta/32`$, while the q+1 source core and two
arbitrary outside helpers use exactly $`N+n+7`$ dirty wires.
The joint orientation bank G and the two preweight rotations have
$`O(N+q+\ell)`$ T cost at $`\ell=N+c`$, for a fixed constant c
chosen for their error budgets. Their actual shared reservation is

```math
\max\{q+1,\ell+1,2k+8\}+2,\qquad k\le n+O(1).
```

For sufficiently large n this equals q+3, as required. This explicit
maximum, rather than an asymptotic width alone, covers the implementation.
Finitely many smaller n retain the existing finite-dimension fallback;
no uniformly better count bound is claimed from that fallback.

If $`\|\widehat G-G\|\le\epsilon_G`$ and
$`\|\widehat E-E\|\le\epsilon_E`$, pair $`\widehat G`$ with
its actual inverse. For the physical raw cascade U_raw,

```math
\|\widehat G U_{\rm raw}\widehat G^\dagger\widehat E
       -U_{\rm can}E\|
\le2\epsilon_G+\epsilon_E+\Delta_q.
\qquad\text{(21)}
```

This is a full-operator comparison. Exact G fixes the preweight banks;
no intermediate initialization, postselection, or dirty reset occurs.
Choosing $`\epsilon_G\le\eta/128`$ and
$`\epsilon_E\le\eta/64`$ makes (21) smaller than eta/16.
Left multiplication by the *ideal* unwinding and Walsh operators preserves
that bound. An implemented unwinding would add its own full-input error
and cost. Equations (20)–(21) remove the serial schedule's particular
precision-width obstruction; they do not provide that implementation.
The directly emitted raw cascade still costs $`O(N+nq)=O(nN)`$.

### What a bankwise unwinding must still implement

The prefix requirement can be made quantitative. Identify one stage bank
with P_d S and let W_before be its chronological ideal prefix. If a unitary
H on that bank and arbitrary dirty work makes its failure column scalar
to error epsilon, then

```math
\|c_d H(W_{\rm before}|_{P_d}\otimes I_b)-c_d I\|
\le\epsilon
\quad\Longrightarrow\quad
\|H-(W_{\rm before}^\dagger|_{P_d}\otimes I_b)\|
\le\epsilon/c_d.
\qquad\text{(22)}
```

This follows by division by c_d and right multiplication by the unitary
prefix inverse. It covers every bank and dirty input, not just one
prepared column. In the last depth $`c_{D-1}=1/16`$, so a bankwise
endpoint-accurate scalarization still needs an endpoint-accurate logical
prefix, with error at most sixteen times its allowed column error.
Coarse prefix substitution or attenuation alone does not remove that
obligation. The statement concerns bank-preserving scalarization; joint
bank mixing or a different complete compiler is outside its scope.
It is not an unrestricted T-count lower bound.

The [depthwise tests](../../tests/test_depth_balanced_rejection.py) check
exact rational weights and width inequalities, full routing and Walsh
action, and a bounded complex cascade with separate coarse interlayers.
They do not emit the unresolved prefix bank or establish a new endpoint
count bound.
