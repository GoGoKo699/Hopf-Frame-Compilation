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
T'\le2T,\qquad G'\le2G+2,\qquad D_T'\le2D_T.          \qquad(5)
```

The inverse has the same T-count, elementary Clifford count, T-depth and
width. The two extra gates are ordinary CNOTs. For target error
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
may be left untouched. Its simultaneous resources obey the same factors
as (5). An approximate query with the complete initialized-isometry
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
error. No target state is prepared. The
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

## 4. What a growing lower bound would now need

Fix $`\eta_0>0`$ and a dirty allocation b. If a Boolean table f requires
XOR-query T-depth at least d at error $`2\eta_0`$ with one clean and
b dirty qubits, then its legal real frame (1) requires T-depth at least
$`d/2`$ at error $`\eta_0`$ with that same allocation. If the query
lower bound imposes count budgets, it must allow the doubled T-count
and $`2G+2`$ Clifford budget in (5). This statement quantifies over
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

## 5. Evidence boundary

[The bounded native tests](../../tests/test_boolean_query_frame_reduction.py)
enumerate every Boolean table through two address bits using exact
$`\mathbb Q(\sqrt2,i)`$ arithmetic, check every logical and dirty
column, retain literal signs and the actual inverse, and include an
entangled-reference check. A separate small native perturbation has
nonzero clean leakage and verifies the two-error hybrid numerically.
These fixtures audit the identities and contract; the arbitrary-size
statements follow from the proofs above. They establish no asymptotic
query lower bound.
