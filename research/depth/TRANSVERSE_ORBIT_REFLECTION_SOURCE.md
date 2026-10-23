# Non-diagonal orbit reflections and a paid transverse trial

[Paid orbit coordinates](AFFINE_ORBIT_REFLECTION_SOURCE.md) ·
[Squared-product interface](COUPLED_ORDER_GRAM_SOURCE.md#squared-product-interface) ·
[Literal controlled relative phase](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#literal-controlled-t-wall)

The affine orbit normalizer permits genuinely non-diagonal reflection
trials without an initialized bank. A broad controlled-bank dressing of
a Hadamard anchor cannot have the required half-projector normalization,
even when its bank branches are arbitrary unitaries. A transverse anchor
removes that obstruction and gives an exact square-root interface. The
explicit higher-band family below pays that interface natively, moves a
bank bit quantum mechanically, and breaks the earlier blind-tail symmetry.
It nevertheless has a uniform wrong-sector witness for every p. These
are restrictions on the displayed interfaces and family, not native-depth
lower bounds.

## 1. Broad bank dressing has an exact normalization law

Use the complete orbit coordinates from the affine proof. The orbit bit
is l=z_0; the remaining principal coordinates are
w=(z without z_0,t). Their names refer to the coordinates after applying
the paid involution C. In these coordinates the first reflection is
$`\bar S=X_l`$ and the target is

```math
\bar Q_x=C_xQ_xC_x=Z_{z_0}Q_x,
```

which acts trivially on l. All w and all additional implementation work
are arbitrary, including reference entanglement. For any two unitaries
V_0,V_1 on the remaining space, which may also include the address, put

```math
V=\begin{pmatrix}V_0&0\\0&V_1\end{pmatrix},\qquad
A=V_0^\dagger V_1,\qquad
R_\theta=V^\dagger(\cos\theta Z_l+\sin\theta X_l)V.
                                                               \qquad(1)
```

This is a Hermitian involution, with c=cos theta and s=sin theta:

```math
R_\theta=\begin{pmatrix}cI&sA\\sA^\dagger&-cI\end{pmatrix}.
```

Direct multiplication, without a
commutation assumption on A and the target, gives

```math
\begin{aligned}
\mathrm{Re}[(R_\theta X_l)^2]
 &=I_l\otimes\bigl(s^2\mathrm{Re}(A^2)-c^2I\bigr),\\
G_\theta:=\frac{I+\mathrm{Re}[(R_\theta X_l)^2]}4
 &=I_l\otimes\frac{s^2(I+\mathrm{Re}(A^2))}4
 =I_l\otimes\frac{s^2(\mathrm{Re}A)^2}2.
                                                               \qquad(2)
\end{aligned}
```

The off-diagonal blocks of the square are opposite Hermitian operators
and cancel on taking the real part; unitarity gives
$`I+\mathrm{Re}(A^2)=2(\mathrm{Re}A)^2`$. In particular
$`0\le G_\theta\le s^2I/2`$. On any positive target vector,
the desired Gram is 1/2, so

```math
\left\|G_\theta-\frac{I+\bar Q}4\right\|
\ge\frac{\cos^2\theta}{2}.                              \qquad(3)
```

For the native Hadamard anchor theta=pi/4 this floor is 1/4. It holds
for every bank width, every number of bands inside V_0,V_1, and even
branch unitaries that move the address. No branch return or source-depth
assumption is needed for this algebraic screen. Dressing by a unitary
commuting with X_l instead gives the unchanged Hadamard Gram I/4.
These statements do not cover arbitrary orbit-bit mixing around the
anchor. General theta in (1) is an algebraic parameter, not a supplied
native fine-angle gate.

## 2. The transverse anchor reduces to a complete square root

For theta=pi/2, choose V_0=I,V_1=A and use the literal X anchor:

```math
\bar R=\mathrm C_l(A)^\dagger X_l\mathrm C_l(A)
 =\begin{pmatrix}0&A\\A^\dagger&0\end{pmatrix},\qquad
(\bar R X_l)^2=\begin{pmatrix}A^2&0\\0&(A^\dagger)^2\end{pmatrix}.
                                                               \qquad(4)
```

Thus the first missing operation is a unitary A whose **complete** square
is the remaining charge. Exact equality to $`I_l\otimes\bar Q`$
holds precisely when $`A^2=\bar Q`$. Approximate direct error epsilon
holds precisely when $`\|A^2-\bar Q\|\le\varepsilon`$, since taking
the adjoint gives the same norm for the other block. A native circuit
for A alone does not automatically pay its coherent control; the
following trial explicitly does so. Inserting a supplied square root
of the lookup would be circular.

## 3. A broad non-diagonal relative-phase family

Retain p at least one, m=2^p, n=pm, r=(m-1)p^2 and
$`s\le(m-1)p^3`$. The paid diagonal relative phase is

```math
D_x|z,t\rangle=\omega^{|t+Pz+B_xz|-|t+Pz|}|z,t\rangle,
\qquad \omega=e^{i\pi/4},\qquad D_0=I.
                                                               \qquad(5)
```

It is independent of l. Write b=t_0 and a=z_p, the first binary
coordinate of the second field block. Let B be the fixed CNOT wall
adding every w coordinate other than a,b into a. There are 2n-3
CNOTs. Define the chronological mixer B followed by H_b, so
$`F=H_bB`$ and $`G_F=2n-2`$. Set

```math
E_x=F^\dagger D_xF,\qquad
A_x=E_x^\dagger S_bE_xD_x,\qquad
R_x=C_x\bar R_xC_x.                                    \qquad(6)
```

The S in S_b is the native single-qubit phase gate, not the affine
reflection. Equations (4)--(6) are the actual parameterized family;
their full squared-product action is

```math
W_x=(R_xS_x)^2
 =C_x\begin{pmatrix}A_x^2&0\\0&(A_x^\dagger)^2\end{pmatrix}C_x.
                                                               \qquad(7)
```

At x=0, E=D=I, so A=S_b and (7) gives the entire correct zero-address
charge. The mixer has only one Hadamard, but its control parity reaches
every other orbit coordinate. At x=1 the resulting A does not commute
with the former blind tail X; this is not a word in that old commutant.
It is also not a fixed-address Clifford trial: Section 4 certifies a
Pauli whose transported axis differs between invariant bank fibers.

### Literal control, native schedule, inverses and full return

Denote the already-paid controlled relative-phase bounds by

```math
\begin{aligned}
T_D^{\rm c}&\le16r+140n,&d_D^{\rm c}&\le156,\\
G_D^{\rm c}&\le2G_B+178n,&h_D^{\rm c}&\le\max(r,4n),
\qquad G_B\le2n+16r+2s.
\end{aligned}                                                   \qquad(8)
```

They come from two returned quadratic queries and two private controlled-T
walls, not a controlled-query oracle. Conjugating by F pays controlled E
with two ordinary F calls. Controlled S_b is the exact three-T,
two-CNOT controlled-S word; its actual inverse has the same bounds.
Consequently V=$`\mathrm C_l(A)`$ and its actual inverse have

```math
\begin{aligned}
T_V&\le3T_D^{\rm c}+3,&d_V&\le3d_D^{\rm c}+3,\\
G_V&\le3G_D^{\rm c}+4G_F+2,&h_V&\le h_D^{\rm c}.
\end{aligned}                                                   \qquad(9)
```

The three controlled D calls are those in D,E,E-dagger. Logical l is
arbitrary and coherent, not a clean signal. The controlled-S T targets
can share logical wires because this constant-size word is executed
serially. Within each controlled D wall all physical T targets are its
private helpers, as certified in its proof home. Pools are reused only
after complete return.

Let T_C,G_C,d_C,h_C be the affine proof's normalizer bounds, with
d_C at most 40. The native direct word uses two C calls, four V or
actual-V-inverse calls, and four literal X_l=HSSH words: cancel only
the adjacent *complete* C involutions between (4) and (7). Thus

```math
\begin{aligned}
T_W&\le2T_C+4T_V,&d_W&\le2d_C+4d_V\le1964,\\
G_W&\le2G_C+4G_V+16,&h_W&\le\max(h_C,h_V).
\end{aligned}                                                   \qquad(10)
```

No native fine-angle approximation is present: all words and returns
are exact. Actual inverses reverse the literal emitted words; neither
A nor its native implementation is declared Hermitian. Arbitrary
address, bank, helper and reference inputs are included. Count and
dirty width are O(mp^3), Clifford count O(mp^3), before the fixed
complete-indicator envelope. That envelope adds
$`4c_n+2n+2m`$ Cliffords and gives dirty width $`2n-m+h_W`$,
where c_n=n^2+3n. Hence the resulting attempted indicator fits the
polynomial intermediate budgets, with constant depth, but not the
required accuracy proved false next. The Gram below is the exact
order-interface diagnostic; no unpriced controlled-R source emitter
is claimed or needed for this direct zero-clean trial.

## 4. An all-size invariant fiber gives the wrong complete charge

Take the field address x=1. In the standard binary basis,
$`B_1z=z\mathbin{\gg}p`$. Set all z coordinates other than l to zero,
leave l arbitrary, and set $`t=e_{n-p}+b e_0`$, with b=t_0 variable.
The tail and b are distinct, even
at p=1. This two-dimensional w fiber, for either l, is invariant under
D,E,A. On it,
D=I. Inside E, B changes z_p from zero to one, while H_b mixes only
b. The exponent in (5) is then 1-2b, because B_1e_p=e_0;
the other weight terms cancel. Therefore on this entire fiber,

```math
E=e^{i\pi X_b/4},\qquad
A=E^\dagger S_bE,\qquad A^2=(A^\dagger)^2=-Y_b.
                                                               \qquad(11)
```

On the corresponding tail-zero fiber, E=I and A=S_b. Hence flipping
the tail does not commute with A, for every p. More strongly,
$`A^\dagger X_bA=-Y_b`$ on the tail-zero fiber but equals $`-Z_b`$
on the tail-one fiber. Both fibers fix every other bank bit, and A
preserves them. No single Pauli can have those two different b axes,
so A at x=1 is genuinely non-Clifford. This conclusion is about its
complete operator, not merely the presence of T gates in an emitter.

By contrast v_1 has the field unit in every block, so the remaining
charge is $`\bar Q=-Z_b`$. Both orbit-bit blocks in (7) consequently
have exact restricted direct defect sqrt(2). The restricted order Gram
and target are respectively

```math
G=\frac{I-Y_b}{4},\qquad G_* =\frac{I-Z_b}{4},\qquad
\|G-G_*\|=\frac{\sqrt2}{4}.                            \qquad(12)
```

In particular their diagonal entries differ by 1/4 on b=0 or b=1;
G is not the scalar I/4 operator. This yields an all-size global lower
bound sqrt(2)/4 on the trial's Gram error, and sqrt(2) on its direct
error, for every p at least one. It is a sector-action failure despite
correct base normalization and a paid sublogarithmic schedule.

## 5. Evidence, scope and reopening condition

The normalization law (2)--(3), square-root equivalence (4), native
ledger (8)--(10), and invariant-fiber obstruction (11)--(12) are analytic.
The exact guards in
[`test_transverse_orbit_reflection_source.py`](../../tests/test_transverse_orbit_reflection_source.py)
check noncommuting block examples of (2), the literal p=1 A word on
coherent address and arbitrary private-helper columns, its actual inverse,
the full p=1 transverse square, and both small-size witness fibers.
Finite checks do not supply the asymptotic proof or a general source.

Pause this displayed completion, not the complete-frame objective. A
different A must provide its full square on both lookup sectors; replacing
the mixer or adding an extra band alone does not establish that identity.
The Hadamard-dressing cap suggests using a transverse or genuinely
multi-axis orbit-bit interface, not optimizing recovery constants.
Approximate scalar zero reflection remains a paid shallow high-degree
operation, so neither Boolean degree nor this wrong-sector example is
a general fixed-accuracy depth invariant. The complete-frame resource
frontier is unchanged. These are internal analytic proofs and exact
finite checks, not external mathematical review.
