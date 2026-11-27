# Boolean XOR queries inside the complete real frame

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Frame convention](../../docs/FAULT_TOLERANT_COMPILER.md#1-target-resources-and-theorem) ·
[Exact dirty query R91](PARALLEL_BILINEAR_WORD_QUERY.md)

The deepest corner subfamily of the prescribed real frame can implement
an arbitrary Boolean XOR query through two frame calls and two CNOTs.
The second frame call is its **actual inverse**. The reduction preserves
the number of clean and dirty qubits, including arbitrary dirty inputs
and external references. It gives a representation-independent transfer
of query lower bounds; it does not itself establish a growing lower bound.

## 1. A legal complete-frame subfamily

Let $`s\ge0`$, $`n=s+1`$, and let $`f:\{0,1\}^s\to\{0,1\}`$
be any prescribed Boolean table. Set every angle above the deepest Hopf
layer to zero and set its deepest angle at address x to
$`\pi f(x)/2`$. Under the repository's literal convention,

```math
R_y(\theta)=e^{-i\theta Y},\qquad
K=R_y(\pi/2)=XZ=-iY,
```

the resulting **entire** frame is

```math
U_f=\sum_x|x\rangle\langle x|\otimes K_t^{f(x)}.       \qquad(1)
```

Both columns of the target t are prescribed, including the minus sign
in $`K|1\rangle=-|0\rangle`$. All upper singular angles are specified;
no state-preparation completion is substituted. This is a subfamily of
both the real family and the complex magnitude family with zero leaf phases.

## 2. Two calls implement a dirty XOR query

Take another arbitrary qubit u and let $`C=\mathrm{CNOT}_{t\to u}`$.
The Pauli conjugation identity

```math
CK_tC=K_tX_u
```

holds with its literal sign: C sends $`X_t`$ to $`X_tX_u`$ and fixes
$`Z_t`$. Hence, for either value of f,

```math
C(K_t^{f}\otimes I_u)C((K_t^{f})^\dagger\otimes I_u)
=I_t\otimes X_u^f.
```

Taking the address-controlled direct sum proves

```math
CU_fCU_f^\dagger
=Q_f^{(u)}\otimes I_t,\qquad
Q_f^{(u)}|x,u\rangle=|x,u\oplus f(x)\rangle.           \qquad(2)
```

Matrix products act rightmost first. The chronological word is therefore
$`U_f^\dagger,C,U_f,C`$. The t wire is returned on every input, not
initialized. Replacing the inverse by another forward call is wrong:
$`K^2=-I`$ would leave the input-dependent phase $`(-1)^{f(x)}`$.

### Simultaneous allocation

Suppose a compiler for (1) uses c initialized qubits and b arbitrary dirty
qubits, with $`b\ge1`$. Choose u from its dirty pool. In the resulting
query, x and u are the logical address and output, while the old logical
wire t joins the remaining $`b-1`$ old dirty wires. Thus the query still
has **c clean and b dirty qubits**, and exactly the same total physical
width. The frame compiler may use u internally as arbitrary dirty work;
its complete-input contract suffices at the two call boundaries.
The same pool is reused by the inverse and forward calls. There is no
additional simultaneously live reservation.

### Approximation, leakage, and the actual inverse

Write $`J=J_c`$ for clean initialization, and let A denote (1) tensored
with identity on the original dirty pool. Suppose a native circuit V obeys

```math
\|VJ-JA\|_{\rm op}\le\eta.                            \qquad(3)
```

No exact clean return is assumed. The actual circuit inverse satisfies

```math
V^\dagger J-JA^\dagger
=V^\dagger(J- VJA^\dagger)
=V^\dagger(JA-VJ)A^\dagger,
```

so its initialized-isometry error is also at most eta. C acts only on
nonclean wires and therefore preserves the initialized subspace. The
right-to-left replacement of the inverse and then the forward call gives

```math
\|CVCV^\dagger J-J(CACA^\dagger)\|_{\rm op}\le2\eta.  \qquad(4)
```

The hybrid with exact calls stays in the initialized subspace. Error
vectors, including leakage, are propagated only by unitaries, whose
operator norms equal one. Thus (4) does not apply (3) to a leaked input.
By (2), the right side's exact target is the XOR query with complete
dirty return. Tensoring an arbitrary external reference preserves (4).
Using the actual inverse also cancels any common scalar in V; no
input-dependent target phase has been discarded.

If V has simultaneous resources $`(T,G,D_T)`$, this construction has

```math
T'\le2T,\qquad
G'\le G(V)+G(V^\dagger)+2\le4G+2,\qquad
D_T'\le2D_T.                                         \qquad(5)
```

The inverse has the same T-count, T-depth and width. Its elementary
Clifford count is charged separately: in the fixed H, S, CNOT alphabet,
$`S^\dagger=S^3`$, so reversing the supplied word gives
$`G(V^\dagger)\le3G(V)`$. The two extra gates are ordinary CNOTs.
For target error
$`2^{-L}`$, it is enough to invoke the frame compiler at precision
$`L+1`$. Its literal width threshold must be checked at that precision;
the reduction does not erase this increment.

## 3. The converse on this subfamily

An exact XOR query on target t gives (1) using

```math
U_f=Q_f^{(t)}H_tQ_f^{(t)}H_t,                           \qquad(6)
```

because $`H X^f H=Z^f`$ and $`X^fZ^f=K^f`$ for f=0,1.
The chronological word is $`H_t,Q_f^{(t)},H_t,Q_f^{(t)}`$.
Two calls reuse the same dirty pool. If the exact query needs no clean
work, neither does this frame construction; one designated clean wire
may be left untouched. Both calls are forward calls, so its simultaneous
resources satisfy $`T'\le2T(Q_f)`$, $`G'\le2G(Q_f)+2`$ and
$`D_T'\le2D_T(Q_f)`$. An approximate query with the complete initialized-isometry
contract also gives error at most twice its query error, by the same
two-call hybrid.

In particular R91, with one output bit, gives this deepest corner
subfamily an exact zero-clean construction at $`b\ge C(n+1)`$, for
one absolute constant C, with

```math
T=O(\sqrt N+N/b),\qquad G=O(N),\qquad
D_T=O(N/b^2+\log_2(n+1)),\qquad N=2^n.                \qquad(7)
```

Here the query has $`2^{n-1}`$ rows. The factor of two and the two
Hadamards are included; no precision source, recovery, or uncharged
oracle is introduced. This exact upper bound applies to the subfamily
(1), not to arbitrary signed corner frames or the full frame family.

### Scalar phases and XOR queries

The corresponding full scalar phase is
$`P_f|x\rangle=(-1)^{f(x)}|x\rangle`$. A complete XOR query gives it
on an arbitrary target t through

```math
Q_f^{(t)}Z_t(Q_f^{(t)})^\dagger Z_t=P_f\otimes I_t.      \qquad(8)
```

The chronological word is $`Z_t,(Q_f^{(t)})^\dagger,Z_t,Q_f^{(t)}`$.
It uses two actual query calls and two Pauli Z operations (four S gates
in the stated elementary alphabet), with error at most twice the query
error. Thus $`T'\le2T(Q_f)`$, $`D_T'\le2D_T(Q_f)`$ and
$`G'\le G(Q_f)+G(Q_f^\dagger)+4\le4G(Q_f)+4`$; the inverse
is charged in the same native alphabet. No target state is prepared. The
physical register set is unchanged: viewed as an s-bit scalar-phase
problem, the former XOR output t is one additional arbitrary returned
dirty wire.

Conversely, form the explicit $`(s+1)`$-bit Boolean table
$`g(x,t)=f(x)t`$. Then

```math
H_tP_gH_t=Q_f^{(t)}.                                    \qquad(9)
```

This uses one scalar-phase call, two Hadamards, the same physical
registers and its same clean/dirty helper allocation, with unchanged
error. Its table has twice as many rows as f. This is an invocation of
a compiler for an arbitrary scalar table g, **not** a free controlled
version of a supplied circuit for $`P_f`$. Equations (8)–(9) connect
the full scalar-extraction question to XOR queries with these explicit
target and row-count charges; a small trace or matrix element of a
scalar sign does not satisfy their full-operator premise.

<a id="near-uniform-positive-frames"></a>

### Strictly positive deepest layers arbitrarily close to uniform

The query reduction already applies to a strictly positive two-angle
family. Fix a density window $`0\lt\rho\lt1`$ and choose
an integer $`k\ge3`$ with $`\sin(\pi/k)\le\rho`$. Put
$`\Delta=\pi/(2k)`$. Let $`W_f`$ set every upper angle to $`\pi/4`$
and the deepest angle at x to $`\pi/4+\Delta f(x)`$. The frame
$`W_0`$ for the identically zero table is the positive balanced frame.
These are complete positive frames with zero prescribed leaf phases.
The balanced upper layers send the root to $`|+^s\rangle|0\rangle`$,
so the normalized leaf densities $`2^{s+1}|\langle x,t|W_f|0^{s+1}\rangle|^2`$
at the two leaves below x are one when $`f(x)=0`$ and $`1\pm\sin(2\Delta)`$ when
$`f(x)=1`$. Thus every density belongs to $`[1-\rho,1+\rho]`$.

The common complete upper word cancels with its adjoint. Consequently

```math
W_f W_0^\dagger=\sum_x |x\rangle\langle x|\otimes R_y(\Delta f(x)),
\qquad
(W_f W_0^\dagger)^k=U_f .                                  \qquad(9a)
```

Use actual frame words $`V_f,V_0`$ and the actual inverse $`V_0^\dagger`$
in $`A=(V_fV_0^\dagger)^k`$. If each frame has full initialized-isometry
error at most $`\eta`$, A approximates $`U_f`$ with error at most $`2k\eta`$.
Use $`C A C A^\dagger`$ from (2), with CNOT to an arbitrary borrowed
target. The query error is at most $`4k\eta`$, including all work and
reference correlations. At each hybrid boundary the ideal word returns
the clean inputs, and all actual adjoints obey the same error as their
forward words. With common forward resource bounds T,G,d for the two
frame circuits, the query has

```math
T_{\rm query}\le 4kT,\qquad G_{\rm query}\le8kG+2,\qquad
D_{T,\rm query}\le4kd .                                      \qquad(9b)
```

The same physical allocation is used. As in Section 2, the query
target comes from the old dirty pool and the old last logical bit joins
the query's dirty work. One may pad either frame implementation with
untouched work within the common allocation. Invoking both frame
compilers at precision $`L+\lceil\log_2(4k)\rceil`$ suffices for
query error $`2^{-L}`$. Any literal width threshold for the chosen
frame compiler must be met at this increased precision.

**Converse for this family.** Suppose there is a complete initialized-
isometry implementation of $`Q_f^{(u)}`$ on the s-bit address and
arbitrary target u, using at most one external clean bit and b arbitrary
dirty inputs. In the frame use u from the b dirty wires and allow the
old last logical bit t to serve as an arbitrary work wire inside the
query, together with the other $`b-1`$ dirty inputs. This gives the query
its full dirty allocation, which it returns on the ideal trajectory.

Write $`B=C_u(R_y(\Delta)_t)`$ and $`E=\mathrm{CZ}_{u,t}`$.
The exact word in operator order is

```math
E Q_f^{(u)} B (Q_f^{(u)})^\dagger B^\dagger E
 = (W_fW_0^\dagger)\otimes I_u .                             \qquad(9c)
```

On fixed x and u, the middle four factors return u and rotate t by
$`\Delta[(u\oplus f(x))-u]=\Delta(1-2u)f(x)`$. The two E factors
reverse that rotation precisely at $`u=1`$, since $`ZYZ=-Y`$.
Each identity holds on arbitrary t and u, including superpositions.
Both query calls are unconditional. To emit B, use the identity

```math
B=R_y(\Delta/2)_t\,\mathrm{CNOT}_{u\to t}\,
 R_y(-\Delta/2)_t\,\mathrm{CNOT}_{u\to t}.
```

Synthesize its two literal determinant-one rotations to any fixed
requested tolerance with no extra work and take the actual inverse of
the emitted B. If the query error is $`\epsilon_Q`$, the B word error
is $`\epsilon_B`$ and the paid balanced-word error is $`\epsilon_0`$,
then composing (9c) after that balanced word gives complete $`W_f`$
error at most $`2\epsilon_Q+2\epsilon_B+\epsilon_0`$. All three
can be chosen at most one fifth of the final tolerance. The same b
dirty wires and sole initialized wire suffice on every call. Relative
to the balanced word, the extra T-count and depth are at most twice
the query's plus a fixed rotation synthesis cost, and the extra
Clifford count is bounded by four times the forward query's plus the
fixed local gates.

At fixed precision and square-root dirty width, the paid
[positive balanced compiler](POSITIVE_HAAR_FINGERPRINT.md) fits its
sufficient width $`b\ge17(n+13)`$ for all sufficiently large n. It has
$`O_L(\log\log(n+4))`$ T-depth, $`O_L(\sqrt N)`$ T-count, and
$`O_L(N)`$ Clifford count on one circuit. Equations (9a)--(9c) therefore
give a two-way sublogarithmic-depth transfer between arbitrary returned
Boolean XOR and this positive family, with the stated fixed changes in
precision and allocation constants. The converse applies to the family
with independent angles only at the deepest layer and pays the balanced
circuit once.

## 4. What a growing lower bound would now need

Fix $`\eta_0>0`$ and a dirty allocation b. If a Boolean table f requires
XOR-query T-depth at least d at error $`2\eta_0`$ with one clean and
b dirty qubits, then its legal real frame (1) requires T-depth at least
$`d/2`$ at error $`\eta_0`$ with that same allocation. If the query
lower bound imposes count budgets, it must allow the doubled T-count
and the $`G(V)+G(V^\dagger)+2\le4G(V)+2`$ Clifford budget in (5).
This statement quantifies over
every native representation with arbitrary Clifford interlayers.

Consequently a robust dimension-growing lower bound for arbitrary
dirty Boolean queries would transfer directly to the complete real
frame problem. The converse (6) compares only this deepest subfamily.
It does not reduce arbitrary frames to two Boolean queries or show that
faster queries alone settle optimal general depth.

The missing step remains a query lower-bound invariant valid under
arbitrary Clifford interlayers, approximate clean return, and arbitrary
dirty/reference inputs. Classical multiplicative depth or Boolean
degree does not supply that invariant without a proved quantum
structural lemma. The
[earlier lower-bound scout](RESEARCH_SESSION_RESULT.md#4-reserve-lane-c-two-lower-bound-falsifiers)
already exhibits full-variable derivative rank at constant T-depth, and
the [Gauss-sum encoding](GAUSS_SUM_QUERY_ENCODING.md) hides an arbitrary
table in constant-depth amplitudes without implementing its full query.

A concrete resumption question is: **what representation-independent
constraint does full dirty return impose on the address-dependent
observable $`Q_f^\dagger Z_uQ_f=(-1)^{f(x)}Z_u`$, for a bounded
number of physical T layers?** Any proposed constraint must remain
valid after exact native rewrites and must be stable at a fixed
operator error. The present reduction exposes this question without
claiming its answer.

<a id="raw-decision-count"></a>

### 4.1 Raw decisions already require square-root T-count

Even a raw decision with arbitrarily small positive bias has a
width-independent worst-case T-count bound. The source may move the
address, use arbitrarily many clean and dirty wires, and leave arbitrary
garbage. Its output test is a fixed Hermitian Pauli. The signed-Pauli
rotation normal form is the one used in the inherited
[finite-width circuit counting](../../docs/FAULT_TOLERANT_COMPILER.md#102-fixed-width-coherent-counting)
and the [returned-channel count](NATIVE_LAYER_COUNT_DISTRIBUTION.md#2-a-description-that-does-not-pay-for-dirty-width).
The argument here counts raw bias functions directly, without the return
or fixed-accuracy premises of those bounds. The LKS/GKW counting lineage
and its imported scope remain in the [source catalogue](../../docs/reference/SOURCE_CATALOGUE.md).

Let $`a\ge1`$, $`t\ge0`$, and let V be a native word with at most t
T/T-dagger gates and arbitrary Clifford interlayers. For an a-bit
computational input x, define

```math
b_V(x)=\langle x,0\cdots0|V^\dagger P_{\rm out}V|x,0\cdots0\rangle.
```

All nonaddress inputs are set to zero only to define this mathematical
witness. A promise for arbitrary dirty inputs includes that witness;
no runtime initialization of the dirty pool is being supplied. Allow
the physical width, native word and fixed output Pauli to depend on the
table. Across all such choices, the number of distinct bias functions is
at most

```math
\#\{b_V:T(V)\le t\}
\le (t+1)\,2^{\,2(t+1)^2+a(t+1)+2(t+1)}.                 \qquad(10)
```

In particular, suppose every Boolean table f has a source satisfying

```math
(-1)^{f(x)}b_V(x)>0\qquad\hbox{for every }x.              \qquad(11)
```

The positive bias may depend on x, f and a, with no uniform lower bound.
If t is a common T-count ceiling for those sources, then

```math
2^a\le \log_2(t+1)+2(t+1)^2+a(t+1)+2(t+1),
\qquad t=\Omega(2^{a/2}).                               \qquad(12)
```

#### A finite description independent of physical width

First fix the actual count k and put $`r=k+1`$. Moving the Cliffords
to one boundary expresses $`V^\dagger P_{\rm out}V`$ as one Hermitian
Pauli $`P_0`$ conjugated successively by k signed Pauli rotations
$`e^{-i\pi P_j/8}`$. The sign of $`P_j`$ absorbs a T-dagger
orientation, and the circuit's common scalar cancels. Write

```math
P_j=i^{s_j}X(\alpha_j)Z(\beta_j),\qquad
s_j\in\mathbb Z_4,\qquad 0\le j\le k.
```

The physical binary vectors may have any length. Record only

```math
\begin{aligned}
B_{ij}&=\beta_i\mathbin\cdot\alpha_j,\qquad s=(s_0,\ldots,s_k),\\
K&=\ker\!\left(u\longmapsto\sum_j u_j\alpha_j\right)
       \subseteq\mathbb F_2^r,\\
L(u)&=\left(\sum_j u_j\beta_j\right)_{\rm address}.
\end{aligned}                                                   \qquad(13)
```

These data use at most $`r^2+2r+r^2+ar`$ bits: store B, the phases,
an r-by-r binary basis description of K, and the matrix of L. Counting
all such descriptions overcounts unrealizable data, which is harmless.
For the ordered product $`P(u)=P_0^{u_0}\cdots P_k^{u_k}`$, define

```math
\phi(u)=\sum_j s_ju_j+2\sum_{i\lt j}B_{ij}u_iu_j\pmod4.
```

Literal Pauli multiplication gives

```math
P(u)=i^{\phi(u)}
 X\!\left(\sum_j u_j\alpha_j\right)
 Z\!\left(\sum_j u_j\beta_j\right),
```

and therefore

```math
\langle x,0|P(u)|x,0\rangle=
\begin{cases}
0,&u\notin K,\\
i^{\phi(u)}(-1)^{L(u)\cdot x},&u\in K.
\end{cases}                                                     \qquad(14)
```

The phases in (14) can include plus or minus i: an ordered product of
Hermitian Paulis need not itself be Hermitian. They must be retained.

The matrix $`B+B^{\mathsf T}`$ supplies every commutation sign.
Conjugating a Pauli product Q by a new rotation leaves it unchanged
when it commutes with the axis, and otherwise gives

```math
e^{i\pi P_j/8}Qe^{-i\pi P_j/8}
=\frac{Q+iP_jQ}{\sqrt2}.
```

Thus the coefficients of the formal expansion in the ordered products
P(u) are determined by B and the fixed rotation rule. Physical
dependencies among the generators do not need separate assumptions:
expand formally, then evaluate each term by (14). Different formal
products that represent the same physical Pauli are automatically
combined in that evaluation. Equations (13)--(14) consequently determine
the entire real function $`b_V(x)`$. There are at most
$`2^{2r^2+ar+2r}`$ descriptions at count k. Summing over
$`0\le k\le t`$ proves (10).

Every bias function satisfying (11) determines a unique Boolean table.
Covering all $`2^{2^a}`$ tables gives (12), without discretizing or
lower-bounding their positive biases. This proof also allows arbitrarily
many initialized computational-zero inputs and does not require any
of them to return.

#### Signed-gap sources and the count-slack boundary

The same conclusion applies to a raw
[full-space signed-gap source](PRESCRIBED_SCALAR_TABLES.md#full-space-signed-gap-source)
with a fixed Pauli signal S. In

```math
\mathcal H=\mathrm{Re}(V^\dagger S V S),\qquad
\mathrm{sgn}(\mathcal H)=D_f\otimes I,
```

include the last S as one additional generator and take the real part
of (14). The description now has $`r=k+2`$ generators and the same
asymptotic size. The sign promise makes its zero-work expectation have
strictly the sign of f. Hence these raw sources also require worst-case
$`T(V)=\Omega(2^{a/2})`$, independently of the positive gap's size.
A non-Pauli output test or signal is not a free additional primitive in
this argument; its native implementation needs its own accounting.

In particular, no compiler for every a-bit table can have the uniform
raw bound $`T(V)=O(2^{a/2}/q(a))`$ for an unbounded function q merely
by weakening its bias or gap. This leaves table-dependent savings on
easy functions possible. It does not lower-bound the T-count of an
amplified implementation by its number of source appearances times
$`T(V)`$, or exclude joint native resynthesis across those appearances.
The conclusion is about T-count; it gives no growing T-depth lower
bound at unrestricted width. The depth question in Section 4 remains
open.

**The raw count is tight even at constant depth.** The
[paid controlled Gauss construction](GAUSS_SUM_QUERY_ENCODING.md#constant-depth-raw-decisions)
realizes every a-bit table with $`T,w=O(2^{a/2})`$, $`G=O(2^a)`$
and at most 32 physical T layers. Every computational dirty input has
bias $`2^{-m}(-1)^{f(x)}`$, where $`m=O(2^{a/2})`$.
Its private helpers return, but its principal dirty work need not.
Two fixed product superpositions give opposite deterministic signal
answers at every address, so it fails the coherent sign premise.
Thus even a band-structured count of these raw bias patterns cannot
force growing depth at an unspecified square-root count constant.
A sharper count constant is a separate question. A depth invariant
must use information absent from unrestricted strict-sign patterns,
such as a uniform quantitative bias, the coherent sign condition or
complete return. In particular, the tiny-bias construction does not
exclude a bound that uses only a fixed positive basis-state bias.

## 5. Evidence boundary

[The bounded native tests](../../tests/test_boolean_query_frame_reduction.py)
enumerate every Boolean table through two address bits using exact
$`\mathbb Q(\sqrt2,i)`$ arithmetic, check every logical and dirty
column, retain literal signs and the actual inverse, and include an
entangled-reference check. A separate small native perturbation has
nonzero clean leakage and verifies the two-error hybrid numerically.
These fixtures audit the identities and contract; the arbitrary-size
statements follow from the proofs above. The
[raw-decision fixtures](../../tests/test_raw_decision_count.py) separately
check the signed Pauli cocycle with dependent axes, the compressed
description against actual native conjugation, the extra signal in a
pair source, and physically different workspace embeddings with the
same witness bias. The counting bound is the analytic finite-description
argument above; these bounded checks establish no growing-depth lower
bound.
