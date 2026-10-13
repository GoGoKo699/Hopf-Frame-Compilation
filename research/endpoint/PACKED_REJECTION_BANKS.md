# Packed rejection banks and one final recovery

[Endpoint contract](../../docs/OPEN_PROBLEM.md) · [Joint native words](JOINT_NATIVE_WORDS.md) · [One-clean source](../../docs/ONE_CLEAN_COMPILER.md)

Supported block encodings can store rejected amplitudes in disjoint logical
mode banks instead of an initialized history counter. For the Hopf tree,
those banks fit in the two supplied flags. A diagonal preweight converts
the accepted product into half of the complete target, allowing one final
coherent recovery. The direct source implementation still pays fine
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
