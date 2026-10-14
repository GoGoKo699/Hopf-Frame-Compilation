# Joint native words and returned-dirty compression

[Endpoint](../../docs/OPEN_PROBLEM.md) · [Source-program resynthesis](SOURCE_PROGRAM_RESYNTHESIS.md) · [Paid logical encodings](PAID_LOGICAL_ENCODINGS.md)

A complete native word has precision-sharing freedom that is absent when
each recovered group is compiled separately. This chapter gives a canonical
description of that freedom and a compression theorem for arbitrary
returned dirty work. Neither result supplies an unrestricted linear-T Hopf
compiler or a superlinear lower bound. The separate
[growing-source compiler](GROWING_SOURCE_COMPILER.md) resolves the endpoint
by a different construction.

The Pauli-word and symplectic normal forms are standard. The projected-
Clifford lemma used below is explicitly
[Ma–Joven–Liu, arXiv:2609.11153v1, Theorem C.1](https://arxiv.org/html/2609.11153v1).
The returned-work statement is presented as a specialization of this
machinery, with a self-contained contract proof and no novelty claim.

## 1. Literal words and the full error condition

There are n system qubits, a initialized qubits, and b arbitrary dirty
qubits, with $`n\ge1`$. Let $`J_a`$ insert $`|0^a\rangle`$, and let V be
the desired system unitary. The circuit U must satisfy

```math
\|UJ_a-J_a(V\otimes I_b)\|\le\eta.
```

This is an initialized-isometry norm: all logical columns, clean leakage,
dirty return, external references, and literal phase are included. Arbitrary
finite Clifford work is allowed. At the Hopf endpoint,
$`N=2^n`$, $`n\ge3`$, $`a=2`$, $`b=N+n+7`$, and $`\eta=2^{-N}`$.

Pushing all Cliffords to the left gives the exact normal form

```math
U=CQ,\qquad Q=T(P_t)\cdots T(P_1),\qquad
T(P)=\frac{I+P}{2}+\omega\frac{I-P}{2},\quad \omega=e^{i\pi/4},
```

where each axis is a signed Hermitian Pauli. Inverse T gates obey
$`T(P)^\dagger=S(P)^\dagger T(P)`$, with $`S(P)=T(P)^2`$ Clifford.
This convention preserves the native scalar that a rotor convention might
discard. The normal-form lineage includes
[Gosset–Kliuchnikov–Mosca–Russo, Proposition 1](https://arxiv.org/html/1308.4134v1).

For an explicit path sum, write the binary axis columns as
$`u_j,v_j`$, the sign bit as $`\epsilon_j`$, and
$`P_j=i^{u_j\cdot v_j+2\epsilon_j}X^{u_j}Z^{v_j}`$.
Let A and B have these columns. Then

```math
Q=\sum_{x\in\mathbb F_2^t}\alpha^{t-|x|}\beta^{|x|}
 i^{q(x)}X^{Ax}Z^{Bx},\qquad
\alpha=\frac{1+\omega}{2},\quad\beta=\frac{1-\omega}{2},
```

```math
q(x)=\sum_j(u_j\cdot v_j+2\epsilon_j)x_j
 +2\sum_{j>k}(v_j\cdot u_k)x_jx_k\pmod4.
```

The phase dot products are integer sums reduced modulo four; the labels
$`Ax,Bx`$ are binary sums. The formula follows by expanding the chronological
product and using $`Z^vX^u=(-1)^{u\cdot v}X^uZ^v`$.

Put $`V_b=V\otimes I_b`$ and $`K=J_a^\dagger UJ_a`$. Since both
compared maps are isometries, their squared error Gram matrix is exactly

```math
(UJ_a-J_aV_b)^\dagger(UJ_a-J_aV_b)
=2I-V_b^\dagger K-K^\dagger V_b.
```

Consequently the full contract is equivalent to
$`2I-V_b^\dagger K-K^\dagger V_b\preceq\eta^2I`$.
Accepted-block approximation alone, at error of order eta, does not imply
this eta-isometry guarantee. The path sum is a representation, not a
procedure for choosing a short word that satisfies the inequality.

## 2. Dirty compression theorem

Let $`d_j\in\mathbb F_2^{2b}`$ be the unsigned dirty restriction of
$`P_j`$, and define its alternating Gram matrix by

```math
(\Omega_D)_{ij}=\langle d_i,d_j\rangle_{\rm symp},\qquad
2s=\mathrm{rank}\Omega_D,\qquad
\delta_0=1-\frac1{\sqrt2}.
```

**Theorem.** If $`\eta\lt\delta_0`$, there is a coherent native circuit
with at most t T gates, the same a initialized qubits, and at most s dirty
qubits, whose full initialized-isometry error is at most eta. In particular,

```math
b_{\rm required}\le\frac12\mathrm{rank}\Omega_D
\le\lfloor t/2\rfloor.
```

The circuit contains no measurements or postselection. The projections in
the proof only restrict or analyze the original all-input contract. No
bound on the new Clifford count is asserted.

### 2.1 Freeze the commuting radical and spectators

Let the dirty-axis span have dimension r. Binary symplectic elimination
gives a dirty-only Clifford D that transforms the span into s canonical
anticommuting pairs, $`r-2s`$ commuting Z directions, and spectators.
This is the usual symplectic canonical form; see also
[Fattal et al., Theorem 1 and Lemma 2](https://arxiv.org/pdf/quant-ph/0406168).
Conjugating the complete circuit by D preserves its target and error,
because its required dirty action is identity.

Retain the s active dirty qubits and call all remaining dirty wires R.
Restrict R at input to $`|0_R\rangle`$. The conjugated word fixes R and
reduces to a word $`Q_0`$ on the system, original initialized wires, and
retained dirty wires. Reduced axes $`+I`$ disappear; reduced axes $`-I`$
give the scalar omega, an actual Clifford phase since
$`(HS)^3=\omega I`$. Thus $`T(Q_0)\le t`$.

Writing $`J_R`$ for insertion of the removed zero wires, the restricted
contract is

```math
\|CJ_RQ_0J_a-J_RJ_a(V\otimes I_s)\|\le\eta. \qquad\text{(1)}
```

The insertions act on distinct registers. No special state is assumed for
any retained dirty input.

### 2.2 The projected Clifford cannot attenuate

Set $`B_0=J_R^\dagger CJ_R`$. The established projected-Clifford theorem
of Ma–Joven–Liu gives, for a nonzero block,

```math
B_0=2^{-h/2}F\Pi,\qquad h\in\mathbb Z_{\ge0},
```

with F Clifford and Pi a stabilizer-code projector. For clarity, the
quantization and phase details needed by this contract can be seen directly.
If R has c wires, expand its zero projector:

```math
B_0^\dagger B_0=2^{-c}\sum_{z\in\mathbb F_2^c}
 J_R^\dagger C^\dagger Z_R(z)CJ_R.
```

The conjugated terms form a commuting signed Pauli group. Only the subgroup
H with I/Z restrictions on R survives. Evaluation at $`R=0`$ is a group
homomorphism. If its image contains $`-I`$, the sum is zero. Otherwise the
image G is a stabilizer group with equal-size fibers, and

```math
B_0^\dagger B_0=2^{\dim H-c}
 \left(\frac1{|G|}\sum_{g\in G}g\right)=2^{-h}\Pi,
\qquad h=c-\dim H.
```

Compressing (1) on R gives

```math
\|B_0Q_0J_a-J_a(V\otimes I_s)\|\le\eta. \qquad\text{(2)}
```

If the block vanishes or $`h\ge1`$, the error is at least
$`1-\|B_0\|\ge\delta_0`$. Hence h=0, and $`B_0`$ is a partial
isometry with source Pi and range $`\Pi'=B_0B_0^\dagger`$.

Here F can retain its **literal native phase**. Choose signed stabilizer
encoders for Pi and $`\Pi'`$. In those coordinates, C maps the canonical
source code, with R zero, deterministically onto the canonical range code,
with R zero. Its logical action normalizes logical Paulis, so a signed
tableau supplies a native projective Clifford $`F_0`$. The remaining scalar
is a deterministic basis-transition amplitude of a native Clifford circuit.
It has unit modulus in $`\mathbb Z[1/\sqrt2,i]`$ and hence is
$`\omega^k`$ by the retained
[scalar-ring lemma](../../docs/FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness).
Implement it by $`((HS)^3)^k`$ on any retained wire and undo the encoders.
This yields the literal identity $`B_0=F\Pi`$. The hypothesis $`n\ge1`$
ensures such a wire exists even for a rank-one code.

### 2.3 The entire ideal output code lies in the range

For any stabilizer projector P, either
$`PJ_a=J_a`$, or

```math
\|(I-P)J_a\|\ge1/\sqrt2. \qquad\text{(3)}
```

Indeed, if containment fails, some signed generator g of P fails to act
as identity on the initialized code. An X or Y on an initialized wire
gives $`\|(I-g)J_a/2\|=1/\sqrt2`$. Otherwise its restriction is a
nonidentity signed Pauli on the uninitialized wires, with a minus-one
eigenvector, giving norm one. Since P is contained in the plus-one
eigenspace of g, these are lower bounds for (3).

Equation (2) implies $`\|(I-\Pi')J_a\|\le\eta`$. The strict threshold
therefore forces

```math
\Pi'J_a=J_a.
```

There is no claim that $`Q_0J_a`$ lies in the source code Pi.

### 2.4 A Clifford extension preserves the complete error

Let $`Y=J_a(V\otimes I_s)`$. Since Y lies in the range of $`\Pi'`$,

```math
B_0^\dagger Y=\Pi F^\dagger Y=F^\dagger Y.
```

The actual and ideal maps in each comparison below are isometries. Their
cross-Gram matrices agree:

```math
(CJ_RQ_0J_a)^\dagger J_RY
=J_a^\dagger Q_0^\dagger B_0^\dagger Y
=(FQ_0J_a)^\dagger Y.
```

Expanding their squared error Gram matrices gives the exact equality

```math
\|CJ_RQ_0J_a-J_RY\|=\|FQ_0J_a-Y\|\le\eta.
```

Thus the physical circuit $`FQ_0`$ removes R without adding clean work,
T gates, or error. Tensoring an arbitrary reference preserves these norms.
This proves the theorem.

## 3. Canonical words and the limit of this entropy argument

After eliminating the radical, the retained dirty labels span a
nondegenerate space of dimension 2s. The kernel of their list map is exactly
$`\ker\Omega_D`$. Two spanning lists with the same Gram matrix are therefore
related by a dirty symplectic transformation. A canonical list is determined
algorithmically by $`\Omega_D`$ alone. Each axis additionally needs
$`2(n+a)`$ logical/initialized label bits and one sign bit. The number of
canonical words of length at most t is consequently bounded by

```math
2^{t(t-1)/2+(2(n+a)+1)t}. \qquad\text{(4)}
```

Shorter words can be padded by identity axes. There is no dependence on the
original dirty width.

### Terminal Cliffords contribute only logical cosets

Fix a word Q. Suppose $`C_iQ`$ satisfies the complete contract for
$`V_i`$ with error $`\eta_i`$, for i=1,2. Put
$`\varepsilon=\eta_1+\eta_2`$ and $`K_C=C_2C_1^\dagger`$. Then

```math
\|K_CJ_a-J_a[(V_2V_1^\dagger)\otimes I_b]\|\le\varepsilon.
```

If $`\varepsilon\lt1/\sqrt2`$, (3) forces this Clifford to preserve the
initialized code exactly. Its induced Clifford $`K_{\rm code}`$ satisfies,
for every dirty Hermitian Pauli P,

```math
\|K_{\rm code}PK_{\rm code}^\dagger-P\|
\le2\varepsilon\lt\sqrt2.
```

Distinct signed Hermitian Paulis have norm separation at least $`\sqrt2`$.
Thus every dirty Pauli is fixed, and

```math
K_CJ_a=J_a(L\otimes I_b)
```

for a literal system Clifford L. Successful terminal Cliffords for this Q
therefore give targets covered by at most $`|\mathrm{Cliff}_n|`$ balls of
radius $`2\eta`$, centered on logical Clifford translates of one successful
target. Here the finite literal Clifford group includes its eighth-root
phases. At $`2\eta\lt1/\sqrt2`$, (4) gives the covering bound

```math
2^{t(t-1)/2+(2(n+a)+1)t+\log_2|\mathrm{Cliff}_n|},
\qquad \log_2|\mathrm{Cliff}_n|=O(n^2).
```

### Actual Hopf packing gives only a linear lower bound

For x in the radius-one-half ball in $`\mathbb R^{N-1}`$, set
$`\psi(x)=(\sqrt{1-\|x\|^2},x)`$. This map does not shrink distances,
and every such vector is an allowed real Hopf first column. A maximal
packing with pairwise distance greater than $`4\eta`$ has at least
$`(1/(8\eta))^{N-1}`$ members: its radius-$`4\eta`$ balls cover the
original ball, so comparison of volumes gives the bound. Choose one
complete Hopf frame for each first column. No radius-$`2\eta`$ operator
ball contains two of these frames. Therefore a worst-case t-T compiler obeys

```math
\frac{t(t-1)}2+(2(n+a)+1)t+\log_2|\mathrm{Cliff}_n|
\ge(N-1)\bigl(\log_2(1/\eta)-3\bigr). \qquad\text{(5)}
```

At the endpoint, this yields $`t\ge\sqrt2N-O(n+a)`$, weaker than the
existing explicit $`2N-1`$ precision witness. Ancilla-independent counting
already has prior formulations, including
[Gosset–Kothari–Wu, Theorems 4.1–4.2 and Appendix A](https://arxiv.org/html/2411.04790v3).
Equation (5) is not an improved lower-bound scaling.

## 4. Endpoint consequences and remaining requirement

The theorem concerns jointly synthesized words. Its literal bound implies
that any construction with
$`t\le2(N+n+7)+1`$ can meet the endpoint dirty allocation even if first
described with unlimited borrowed width. An unspecified $`O(N)`$ T-count
with a larger constant does not automatically meet that allocation.
For commuting dirty axes, $`\mathrm{rank}\Omega_D=0`$, all dirty
work is removable. This includes the Clifford-only case.

The alternating dirty Gram matrix has quadratic information capacity in t.
At $`t=\Theta(N)`$, counting does not prevent it from encoding
$`\Theta(N^2)`$ bits of endpoint data. That is a capacity statement, not
an encoding or existence theorem. For this normal-form route, the
constructive requirement remains a
terminating rule selecting the complete axes and final Clifford with
$`t=O(N)`$, the literal dirty reservation, and the full error inequality
in Section 1. Separate high-precision prices for individual rotations or
recovered groups do not establish or refute such a rule.

The [packed rejection-bank construction](PACKED_REJECTION_BANKS.md)
provides a concrete supported word with one final recovery; its fine
native synthesis remains the missing step. The
[exact finite witnesses](../../tests/test_joint_native_words.py) check
literal scalar phases, projected-Clifford attenuation, and error-Gram
preservation with source leakage. They are regressions for these proof
interfaces, not a substitute for the all-width argument.
