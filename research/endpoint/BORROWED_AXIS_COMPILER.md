# Borrowed axes and jointly compiled source orientations

[Packed rejection banks](PACKED_REJECTION_BANKS.md) · [Source resynthesis](SOURCE_PROGRAM_RESYNTHESIS.md) · [Borrowed rotation signal](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations) · [Endpoint](../../docs/OPEN_PROBLEM.md)

A variable-axis exponential can transfer its axis to an arbitrary borrowed
wire, apply an addressed rotation, and transfer it back. This removes the
two initialized helpers from a particular coherent source-basis construction.
For the retained paired source, all normalized rejection orientations can
then be compiled as one addressed bank with one total precision charge.
The result prices the change from raw source blocks to canonical blocks;
it does not compile their chronological product or the fine prefix bank.

All statements below retain literal phases, actual inverses, arbitrary
dirty inputs and external references. Only T-count is optimized. No source
state, evaluator, measurement, reset or postselection is supplied.

## 1. Transfer an axis through an arbitrary wire

Let $`y`$ be an unchanged k-bit address, with $`Q=2^k`$ rows. On
registers disjoint from a temporary wire u, let $`P_y,Q_y`$ be Hermitian
involutions satisfying $`P_yQ_y+Q_yP_y=0`$. Their controlled versions
must have charged full-input implementations. Define

```math
F_y=H_u C_u(P_y)H_u,
\qquad B_y=C_u(Q_y).
```

Both are Hermitian unitaries. Direct conjugation gives

```math
F_yZ_uF_y=Z_uP_y,
\qquad B_yP_yB_y=Z_uP_y,
\qquad [B_y,Z_u]=0.
```

Consequently the following is a **full-operator identity**, without an
initialized input on u:

```math
B_yF_yZ_uF_yB_y=P_y,
```

```math
B_yF_y e^{i\theta_yZ_u}F_yB_y
=e^{i\theta_yP_y}\otimes I_u.
\tag{1}
```

Use the [borrowed-signal real-rotation primitive](../../docs/ONE_CLEAN_COMPILER.md#reusing-an-occupied-port-and-the-precision-core),
with a fixed target Clifford conjugation and the appropriate angle sign,
for the middle Rz rotation. Its ordinary target is u. Its signal v is
another arbitrary wire, distinct from u and the core. At precision
$`\ell\ge5`$,
its full-operator error is below $`43\,2^{-\ell}`$. Surrounding it by
the exact unitaries in (1) preserves that bound.

The middle primitive may reuse the source core even when the preceding
F or B has entangled it with u, y, a logical target, or a reference. Its
ideal action is identity on the entire core. This is a use of its
full-input contract, not an assumption of intermediate dirty return.
Every emitted inverse is the actual circuit inverse.

### Selected Pauli axes, including their scalar phases

For a nonidentity Hermitian Pauli $`P_y`$, an anticommuting Hermitian
Pauli $`Q_y`$ can be chosen classically from a nonidentity coordinate.
Write its complete literal word as $`i^{a(y)}X^{x(y)}Z^{z(y)}`$.
The [whole-word dirty queries](../../docs/OPERATOR_SOURCE_COMPILER.md#self-borrowed-whole-word-queries)
on address $`(y,u)`$ implement the controlled X and Z words. The
conditional scalar $`i^{u a(y)}`$ is also required. The exact
[dirty clock phase](DISJOINT_BLOCK_COMPILER.md#dirty-only-routing-and-boolean-phases)
implements its two Boolean bits at $`O(Q)`$ T cost.

The clock uses two core bits in a completed subroutine. They return
before the next Pauli query and may then resume their ordinary core roles.
Its low clock bit can be an extra address bit, while its target and
selectors remain disjoint. A conservative sufficient shared-core size is
$`m\ge2k+8`$, covering the control, clock and selector margins. Thus
both controlled selected Paulis cost $`O(Q)`$ T gates and
$`O(Qm)`$ Clifford gates. Their row-dependent phases are never discarded.

Equation (1) deliberately excludes scalar axes $`P_y=\pm I`$:
there is no anticommuting involution in that case. Inactive slots can
instead use a fixed nonidentity Pauli with zero angle. The compact
special-unitary synthesis below uses traceless Pauli generators, so it
does not require an unpriced addressed scalar phase.

### Selected source-vector axes

Suppose $`P_y=A_y`$ is a unit Clifford vector in a prescribed set of
anticommuting Majoranas. An unused Majorana K anticommutes with every
$`A_y`$, so $`Q_y=K`$ suffices. If
$`A_y=L_y\Gamma_{p(y)}L_y^\dagger`$ has a charged exact source or
selected-tail implementation, implement its control as

```math
C_u(A_y)=L_y C_u(\Gamma_{p(y)})L_y^\dagger.
\tag{2}
```

The loaders are unconditional. On u=0 their actual inverse cancels
exactly; no gatewise controlled arbitrary loader is assumed. Masked
sources use the completed mask and its actual inverse outside (2).
The selected-tail identity below supplies the same form, with its known
minus sign retained as a Z on u. Each controlled axis then costs
$`O(Q+q)`$ T gates for a q-bit geometric tail.

With a shared core of size

```math
m=\max\{q+1,\ell+1,2k+8\},
\tag{3}
```

one selected source-axis exponential therefore costs
$`O(Q+q+\ell)`$ T gates and uses $`m+2`$ arbitrary dirty wires:
the core and the distinct u,v. It requires no initialized helper.

## 2. Four tail classes suffice for the paired source

Use the q-bit source, $`q\ge5`$, from
[the one-clean compiler](../../docs/ONE_CLEAN_COMPILER.md#2-an-exact-two-tail-source-and-its-dirty-masks),
on q+1 core qubits. Write its two normalized geometric tails as

```math
A=\sum_{j=0}^{q-1}\sqrt{w_j}\,\Gamma_{2j+1},
\qquad
B=\sum_{j=0}^{q-1}\sqrt{w_j}\,\Gamma_{2j+2},
```

```math
M=\Gamma_0/2+A/\sqrt2+B/2.
```

The head and two tails are mutually orthogonal Clifford vectors. Each
tail separately has the ordinary geometric weights $`w_j`$. Its
loader contains q-1 of the native Pauli splitting rotations already
used by the paired source. The unused $`\Gamma_{2q+1}`$ anticommutes
with every head and tail vector.

For the fixed even-tail mask f and any programmed mask g,

```math
N_f=\Gamma_0/2+A/\sqrt2-B/2,
\qquad
N_g=\sigma\Gamma_0/2+A_g/\sqrt2+B_g/2,
```

where $`A_g,B_g`$ reverse signs within their respective tails and
$`\sigma\in\{-1,1\}`$. At each row y, partition each tail into
the positive and negative signs of g. There are at most four nonempty
tail classes and the head: at most five source directions, independently
of q and Q. The union of classes over different rows need not be small.

For a class C in one tail, put

```math
w_C=\sum_{j\in C}w_j,
\qquad
G_C=\sum_{j\in C}\sqrt{w_j/w_C}\,\Gamma_j,
\qquad p=\min C,
```

using local tail indices in this paragraph. The geometric first
coefficient has $`w_p\ge w_C/2`$. The normalized tail beginning at p
and its sign reversal outside C therefore have precisely the bounded
overlaps required by the
[full and half moves](SOURCE_PROGRAM_RESYNTHESIS.md#explicit-full-and-half-moves).
Their nine-exponential class alignment applies inside this tail. Empty
classes and p are determined from the finite sign table; no unknown-real
zero test is needed. Tiny class weight does not introduce an inverse-weight
quantum cost.

The [selected-tail identity](SOURCE_PROGRAM_RESYNTHESIS.md#one-full-loader-selects-every-tail)
also applies separately to either tail loader U:

```math
K_p=U\Gamma_pU^\dagger
=(-\Gamma_{p-1}+A_p)/\sqrt2,
\qquad
A_p=-K_p\Gamma_{p-1}K_p \quad(p>0),
```

with $`A_0=U\Gamma_0U^\dagger`$. The indices are again local to
that tail. Selected Majoranas are charged Pauli queries. The fixed
loaders cost $`O(q)`$ T gates. Their unconditional conjugation around
controlled centers gives the controlled tails needed in (2), including
the p=0 case and all scalar signs.

Apply Section 1 to the constant number of class-axis exponentials,
using $`\Gamma_{2q+1}`$ as the common anticommuting involution.
The ideal class rotors fix Majoranas outside their class. Pack their
at most four representatives and the head into the first five Majorana
coordinates by a row-selected signed coordinate permutation. A constant
number of coordinate-plane rotations suffices. The generator
$`i\Gamma_j\Gamma_k`$ anticommutes with the constituent
$`\Gamma_j`$, so (1) prices these moves as well. It is sufficient
to approximate this constant number of moves at precision ell; an exact
selected-Clifford implementation is not required for the bound. Track
every induced coordinate sign in the compact matrix and retain the full
Spin lift. Pair the emitted approximate basis with its actual inverse.

Thus the coherent basis costs $`O(Q+q+\ell)`$ T gates, with a fixed
multiple of $`2^{-\ell}`$ full-operator error, and the dirty reservation
(3). It puts the relevant source vectors on at most three physical core
qubits. No central parity sector is assigned an eigenvalue or discarded.

## 3. Compile the complete rejection-orientation bank

For the retained raw paired-source coin, use the exact
[accepted and rejected blocks](SOURCE_REUSE_LIMITS.md#an-inverse-logical-pair-can-leave-a-coherent-dirty-polar-operator)

```math
Q_g=\begin{pmatrix}A_g&B_g\\B_g&A_g\end{pmatrix},
\qquad A_g=sI+tX Z,
\qquad B_g=XK_g.
```

Here X and Z act on the logical target, while

```math
K_g=D_g+2tD_f=ML_g,
\qquad
L_g=N_g+2tN_f-(s+t)M,
\qquad \lambda^2=1-s^2-t^2.
\tag{4}
```

The source identities give
$`\{M,L_g\}=0`$ and $`L_g^2=\lambda^2I`$. For lambda positive,
define the literal unitary

```math
\Omega_g=X M(L_g/\lambda).
\tag{5}
```

It satisfies $`B_g=\lambda\Omega_g`$,
$`\Omega_g^\dagger=-\Omega_g`$, and
$`\Omega_g^\dagger A_g\Omega_g=A_g^\dagger`$.
At lambda zero, choose $`\Omega_g=X`$, with identity core body.
This choice still gives the last identity. Choosing an arbitrary
orientation would align the zero rejection, but would not necessarily
give the same lower-right canonical block. Entirely inactive identity
rows may instead choose identity.

**Orientation-bank theorem.** For Q conserved address rows of the raw
paired source, the complete addressed unitary (5), with the stated
zero-lambda convention, has a zero-clean compiler with

```math
T=O(Q+q+\ell),
\qquad
b=\max\{q+1,\ell+1,2k+8\}+2,
\qquad
\|\widehat\Omega-\Omega\otimes I_{\rm work}\|
\le C2^{-\ell},
\tag{6}
```

where C is an absolute constant. The original q+1 source qubits are
included in b; the norm includes their intended orientation action and
return of all additional work. Q is $`2^k`$, with identity padding
allowed. The constants do not depend on a positive lambda's magnitude.

*Proof.* Section 2 puts M and $`L_g/\lambda`$ on at most five
Majoranas, hence three physical core qubits. Their product is an
explicitly computable SU(8) matrix. Each unit Clifford-vector reflection
on eight dimensions has four positive and four negative eigenvalues,
so its determinant is one; the literal product and its sign are retained.
Factor the logical X route out as an exact addressed Pauli operation.

A fixed finite Pauli-rotation schedule covers SU(8). For example,
two-level special-unitary elimination followed by determinant-one
diagonal phases gives a constant-size schedule; expanding its selectors
uses traceless, nonidentity Pauli generators. Certified finite angle-grid
selection approximates the entire compact matrix, including singular
charts and literal phase, without an exact pivot oracle. Section 1
prices every addressed rotation in that fixed schedule. It uses the
same arbitrary core even when the core contains the compact target:
each middle rotation has ideal identity action on that core before
conjugation by F and B.

Combine the approximate coherent basis, compact body and actual basis
inverse. Their errors are bounded by $`2\epsilon_B+\epsilon_V`$,
not by the number of address rows. The number of fine appearances is
an absolute constant, proving (6).

The native sign tables make s and t dyadic, so lambda squared is exactly
rational. Its zero branch is decidable. For positive lambda, certified
normalization may require additional classical precision, but the compact
target remains unitary and its quantum synthesis cost depends on the
requested matrix error. There is no quantum operation with cost
$`1/\lambda`$. All controlled source calls, selector returns and native
phases are charged in Sections 1–2. ∎

The one-tail scalar preweight has the same argument with only two source
vectors and no logical X route. Its scalar zero-defect rows can choose
identity. Fixed Clifford conjugations for the three residual axes, their
logical routing and the bank permutations are separately charged exact
operations. They add a constant number of row tables. In the packed
layout the bank type, depth and prefix form a conserved address with
$`Q=O(N)`$; the changed logical bit is excluded from that address.
The preweight and regular banks can be processed in a constant number
of such families, with exact inactive conditioning on the success flags.

For $`q=N+O(\log n)`$, $`\ell=N+O(1)`$ and $`k\le n+O(1)`$,
(6) is $`O(N)`$ T count and fits $`N+n+7`$ dirty wires for
sufficiently large n. The theorem retains the explicit maximum in (6);
it does not infer a literal workspace bound for every small n merely
from this asymptotic expression.

## 4. One global canonicalization, with its scope

Identify each supported source stage's active success space with its
private failure bank. On that pair of spaces put

```math
G_j=\mathrm{diag}(I,\Omega_j^\dagger).
```

The exact identities above give

```math
G_jQ_jG_j^\dagger
=\begin{pmatrix}
A_j&-\lambda_j I\\
\lambda_j I&A_j^\dagger
\end{pmatrix}.
\tag{7}
```

Lambda may depend on the unchanged row. The zero-lambda convention in
Section 3 also makes (7) valid there. This is the entire raw block, not
only its accepted columns.

Because the rejection banks are disjoint and every other stage fixes
each bank, a single G, identity on success and equal to
$`\Omega_j^\dagger`$ on bank j, conjugates the complete raw program
to its canonical raw program. The intervening success-only coarse words
are unchanged: G commutes with them. Their product with a raw source
block is not silently identified with another block of form (7).

The orientation-bank theorem prices G jointly. Two appearances of G
add $`O(N+q+\ell)`$ T gates; they do not add one precision charge per
stage. Exact G fixes the initialized embedding and commutes with the
success reflection. An actual full-operator approximation to G with
error $`\epsilon_G`$, paired with its actual inverse, changes this
conjugation by at most $`2\epsilon_G`$. Known success/failure
predicates provide the retained full-input inactive conditioning.

This is an additive-cost equivalence of the raw and canonical complete
programs. Substituting the existing raw circuit inside the conjugation
retains its old cost. No linear-cost implementation of the chronological
canonical interior follows from (7).

### Rounded canonical blocks versus constant attenuation

Actual dyadic rows have $`\rho_y=\sqrt{s_y^2+t_y^2}`$ and
$`\lambda_y=\sqrt{1-\rho_y^2}`$. They do not have exactly the
chosen constant attenuation alpha. If
$`\|A_y-\alpha V_y\|\le\epsilon`$, where V is unitary, then

```math
|\rho_y-\alpha|\le\epsilon,
\qquad
|\lambda_y-\sqrt{1-\alpha^2}|
\le {2\epsilon\over\sqrt{1-\alpha^2}}.
```

The full canonical-block difference is therefore at most

```math
\epsilon+{2\epsilon\over\sqrt{1-\alpha^2}}.
\tag{8}
```

For the retained polynomial attenuation
$`\alpha_d=1-4/[100(n-d)^2]`$,
$`\sqrt{1-\alpha_d^2}\ge\sqrt7/[10(n-d)]`$. Equation (8) is
at most $`9(n-d)\epsilon`$. Summing a bounded number of axes per
depth gives $`O(n^2 2^{-q})`$, rather than the accepted-block estimate
$`O(n2^{-q})`$. A full-canonical comparison at endpoint error thus
uses $`q=N+2\log_2 n+O(1)`$; its workspace must again be checked
through (6). This polynomial schedule fits asymptotically. Exponentially
smaller balanced defects cannot inherit the same precision choice
without another analysis of (8).

In the ideal constant-attenuation model, after the complete orientation
and logical-prefix alignment, the scalar history vectors can be normalized
jointly by the existing
[zero-clean real reflection bank](CONJUGATE_SIGNAL_RECOVERY.md#zero-clean-real-reflection-and-state-extension-banks).
For the actual rounded rows, the positive contractions can depend on the
prefix and need not commute with earlier logical operations. Equation (7)
does not restore an exact factorization $`W_jD_j`$ or make those actual
histories scalar. Passing to the ideal model requires the comparison just
priced in (8). The scalar operation itself does not price its preceding
target-dependent prefix unwinding. The unresolved task is joint native
synthesis across the chronological stages, including the fine transported
prefixes of that ideal model, or an alternative canonical interior which
avoids them. No unrestricted endpoint bound changes here.

## 5. Verification

The [exact borrowed-axis tests](../../tests/test_borrowed_pauli_rotations.py)
retain every signal and core input column, complex Hermitian Pauli phases,
actual inverses and coherent row selection. They verify the algebra behind
(1). The [paired-source orientation tests](../../tests/test_paired_source_orientation.py)
check the complete small raw word, its rejection algebra, normalized bank
gauge and class dimension. Neither suite is an arbitrary-angle native
emitter or the asymptotic orientation compiler. The resource conclusion
uses the charged native primitives and constant-dimensional reduction
proved above. Finite checks are not used as a substitute for that reduction
or for the unresolved prefix compiler.
