# Two complete key-mixing tests for observable transport

[Observable contract and paid queries](OBSERVABLE_TRANSPORT_SOURCE.md) ·
[Output aliases](OUTPUT_ALIASED_CHARACTER_SOURCE.md) ·
[Quadratic XOR ledger](FULL_BANK_CHARACTER_TRANSPORT.md#5-the-quadratic-inverse-has-a-paid-xor-query-not-a-paid-in-place-map)

Two actual changes to the previously retained grammar can be assessed
without treating a dirty Fourier label as a prepared state. An outer
opposite-direction query fails a top-B readout even with arbitrary
unitary fibers in the preceding source. Clifford fibers also give a
constant defect for the top-A readout. A different candidate moves the
Fourier key through the address insertion itself. Its relative Clifford
curvature is either non-Pauli or has only affine address labels.

These are complete observable defects for the specified words and
wirings. They are not lower bounds on unrestricted native depth, on
words with further separated key mixers, or on words with several
noncommuting address insertions.

## 1. Coordinates, read-only address and physical aliases

Retain $`r=b+1`$, $`m=2^b`$, $`w=2^r-1`$ and the nonconstant
coefficient banks A and B. Write t for the coordinate vector selecting
the top coefficient, and write

```math
Q=Q_k(A;B),\qquad R=Q_k(B;A),\qquad
M=H_AQH_A,\qquad q_v(a)=v\cdot k(a).
```

Both Q and R are their own exact full-space inverse. Every paid private
query wire returns, including correlations with the two banks and
external references. Private helper wires H may be included with A in
a fiber description; R acts identically on H.

The address x is read-only in all results below. The sources do not
depend on x; x enters through the explicitly stated insertion. Define
the character label c(x) on A by

```math
Z_A^{c(x)}|a\rangle=(-1)^{\lambda(aA_x)}|a\rangle,
\qquad A_x=\prod_{i=1}^b(1+x_ie_i)-1.
```

There are two different output assignments, which must not be mixed.
For the **top-B assignment**, put the m minus one nonconstant output
characters in their A complement coordinates and put the constant
output $`Y_\varnothing`$ in the top coordinate of B. The required
observable is

```math
T_x^{(B)}=Z_B(t)Z_A^{c(x)}\otimes I_H.
```

Every remaining A and B coordinate is arbitrary private work. For the
**top-A assignment**, use the full output alias in A, including its
top coordinate, and make the entire B bank private. The required
observable is

```math
T_x^{(A)}=Z_A^{t+c(x)}\otimes I_{B,H}.
```

Each assignment uses $`3m-2`$ private bank coordinates, before the
query-private pool. Neither supplies a clean coordinate. The complete
operator, rather than its restriction to any one helper input, is the
object being compared.

## 2. One final opposite-direction query: the top-B screen

Let K preserve the B Fourier key, with arbitrary unitary fibers on A,H:

```math
K=\sum_v C_v\otimes|v\rangle\langle v|_B.
```

Let the unitary insertion P_x preserve that same key, and set

```math
S=RK,\qquad J_x=SP_xS^\dagger=R J_{0,x}R,
\qquad J_{0,x}=KP_xK^\dagger.
```

No Clifford or affine-label premise is needed in this section. The
fiber action of P_x may be arbitrary. Since R commutes with the
computational probe $`Z_B(t)`$, conjugating the complete comparison by
R gives

```math
R J_x^\dagger Z_B(t)J_x R
 =J_{0,x}^\dagger Z_B(t)J_{0,x}.
```

This candidate has exactly the B-Fourier key shift t. In contrast,

```math
R T_x^{(B)}R
 =Z_B(t)Z_A^{c(x)}D_{q_{c(x)}}(B).
```

If $`c(x)\ne0`$, the coefficient of this target at that same key
shift is zero:

```math
2^{-w}\sum_\beta(-1)^{c(x)\cdot k(\beta)}=0.
```

This uses the full-input bijectivity of k; it does not average over a
prepared key. The candidate and target are therefore orthogonal in
normalized Hilbert--Schmidt inner product. They are unitaries, so their
squared normalized Hilbert--Schmidt distance is two. Operator norm
dominates that normalized norm, yielding

```math
\bigl\|J^\dagger Z_B(t)J-T^{(B)}\bigr\|_{\rm op}
\ge\sqrt2.                                                   \qquad(1)
```

For the character family, any nonzero address supplies $`c(x)\ne0`$.
The argument allows quantum helper fibers and any number of letters
inside K. It covers the actual source $`S=RMQ`$. It does not cover
an additional mixer before R, a key-moving insertion, an alternative
probe, or address motion.

## 3. The top-A screen with Clifford fibers

Now require every C_v to be a Clifford on A,H and P_x to have a Pauli
fiber whose binary label is affine in x. It may have a key-dependent
scalar phase; this does not affect the proof. Then

```math
J_{0,x}=\sum_v J_{x,v}\otimes|v\rangle\langle v|,
```

where every J_{x,v} is a Pauli with label affine in x. Conjugating the
probe and target by R gives, respectively,

```math
Z_A^tD_{q_t}(B),\qquad
Z_A^{t+c(x)}D_{q_{t+c(x)}}(B).
```

Let $`\widehat g_h(s)=2^{-w}\sum_\beta(-1)^{h\cdot k(\beta)+s\cdot\beta}`$.
The relative candidate
and target blocks from v to v+s are

```math
\widehat g_t(s)J_{x,v+s}^\dagger Z_A^tJ_{x,v},
\qquad
\widehat g_{t+c(x)}(s)Z_A^{t+c(x)}.                           \qquad(2)
```

The two Walsh spectra in (2) have the same support and magnitude:

```math
\widehat g_h(s)=0\quad(s\cdot t=0),\qquad
|\widehat g_h(s)|=2^{-(w-1)/2}\quad(s\cdot t=1),              \qquad(3)
```

for h=t and h=t+c(x). To verify this rather than assume a flatness
promise, split $`\beta=a+e_rz`$, where a is in the augmentation
ideal of $`R_b`$ and z is in its full coefficient space. The quadratic
inverse has the literal decomposition

```math
k(\beta)=k(a)+e_r(1+a)z.
```

For $`g=1+A_x`$, the relevant quadratic form is
$`\lambda_b(g(1+a)z)`$. Multiplication by g is invertible because
$`g^2=1`$. The Frobenius pairing $`\lambda_b(az)`$ is
nondegenerate on the full algebra. Its restriction has rank
$`2(m-1)=w-1`$, and radical precisely the top-coordinate vector t.
The quadratic form takes value one at t. The standard finite Fourier
sum for a quadratic form therefore has support $`s\cdot t=1`$
and magnitude in (3). Equivalently, sum first over z to obtain one
linear constraint for each a; the remaining solution is unique modulo
that one-dimensional radical.

For fixed v,s, the binary label of the first Pauli in (2) is affine in
x. The $`2^b`$ labels $`t+c(x)`$ are affinely independent: append a
constant coordinate and their matrix is the invertible subset-monomial
evaluation matrix. An affine image of $`\mathbb F_2^b`$ has dimension
at most b, and hence can contain at most b+1 of those labels. Thus each
v,s block can have the correct binary label at no more than b+1
addresses. A mismatched Pauli label has zero trace overlap, regardless
of literal phases or helper width.

By (3), every nonzero block has the same weight. Average the normalized
Hilbert--Schmidt overlap over all addresses, all v and all supported s.
Its real part is at most $`(b+1)/2^b`$. Consequently

```math
\bigl\|J^\dagger Z_A^tJ-T^{(A)}\bigr\|_{\rm op}
\ge\sqrt{2\left(1-\frac{b+1}{2^b}\right)}\qquad(b\ge2).    \qquad(4)
```

The first value is $`1/\sqrt2`$, and the floor tends to
$`\sqrt2`$. Equation (4) concerns exactly the top-A assignment and
the stated Clifford-fiber insertion grammar. Section 2's stronger
quantum-fiber premise cannot simply be substituted into this proof.

## 4. A moving Fourier key and its exact Clifford curvature

The insertion can instead move the preserved Fourier key. This is
outside both the original fixed-key insertion theorem and Section 2.
Take the actual paid source

```math
S=MQ,
\qquad S_v=B_vA_v,\qquad
A_v=D_{q_v},\qquad B_v=H_AD_{q_v}H_A.
```

Here the B_v notation denotes a Clifford, not the bank B. In this
section write $`\mathcal B_v`$ for the polar matrix of q_v. Let

```math
P_x=P_A(x)\otimes Z_B^{u(x)}X_B^{z(x)},
\qquad J_x=SP_xS^\dagger,                                   \qquad(5)
```

where the labels of P_A, u and z are affine in the read-only address.
The A Pauli may act on independent dirty spectators as well. Any
literal scalar phase in (5) is permitted. Such phases and the X_B
factor affect only scalar block phases. The $`Z_B^{u(x)}`$ factor
translates the B-Fourier key by u(x).

First set P_A=I and z=0, and abbreviate u=u(x). The transported top-B
probe still has total key shift t. Its exact relative block is

```math
\begin{aligned}
\mathcal R_{v,t}(u)
 &=S_{v+t}S_{v+t+u}^\dagger S_{v+u}S_v^\dagger\\
 &=B_v\bigl(B_tA_uB_tA_u\bigr)B_v.                         \qquad(6)
\end{aligned}
```

The second equality follows from $`A_{v+u}=A_vA_u`$ and
$`B_{v+u}=B_vB_u`$. It retains all phases: these are equalities of
actual quadratic diagonal gates and their Hadamard conjugates, not only
symplectic matrices. Inserting P_A in (5) puts Pauli factors inside
the first line of (6). It does not change that line's symplectic action.

Put $`T=\mathcal B_t`$ and $`U=\mathcal B_u`$. On a Pauli label
$`(a,z)`$, B_t and A_u have symplectic matrices

```math
\begin{pmatrix}I&T\\0&I\end{pmatrix},\qquad
\begin{pmatrix}I&0\\U&I\end{pmatrix}.
```

Their commutator in (6) has matrix

```math
\begin{pmatrix}
 I+TU+(TU)^2&TUT\\
 UTU&I+UT
\end{pmatrix}.                                             \qquad(7)
```

Every polar matrix has zero top row and column. The top polar matrix
T is complement matching on the proper nonempty subsets, so its
radical is exactly the top coordinate. The lower-right block of (7)
shows that identity symplectic action requires UT=0. The stated radical
and zero top column then force U=0. Conversely U=0 makes (7) identity.
Finally, $`\mathcal B_u=0`$ holds precisely when u is supported on
singleton coordinates: every coefficient of size at least two occurs
as a disjoint-pair polar entry.

This gives a complete dichotomy for (5).

* If u(x) has a nonsingleton component at some address, every relative
  block there is a non-Pauli Clifford, even with the additional A Pauli
  and X_B factor. Its distance from every Pauli is at least
  $`\sqrt{2-\sqrt2}`$. Indeed the normalized trace overlap of a
  non-Pauli Clifford with a Pauli has magnitude at most
  $`1/\sqrt2`$: their normalized Choi vectors are distinct stabilizer
  states, whose nonunit overlap is at most $`1/\sqrt2`$. Normalized
  Hilbert--Schmidt distance then gives the displayed operator floor.
* If every u(x) is singleton-supported, A_u and B_u are Paulis. The
  exact key-transition factor
  $`S_{v+u}S_v^\dagger=B_vB_uA_uB_v`$ is a Pauli with label linear
  in u at fixed v. Multiplying by $`S_vP_A(x)S_v^\dagger`$ leaves
  an affine x label. The relative readout block is the product of the
  transition block at v+t with that at v, up to scalar phase. Its label
  is therefore affine in x. For $`b\ge2`$, the nonaffine character
  label c(x) gives the retained $`\sqrt2`$ floor.

Thus all insertions (5) obey, in the top-B assignment,

```math
\bigl\|J^\dagger Z_B(t)J-T^{(B)}\bigr\|_{\rm op}
\ge\sqrt{2-\sqrt2}\qquad(b\ge2).                           \qquad(8)
```

This proof does not prohibit a longer source with a different curvature,
multiple noncommuting address insertions, or a source whose fibers cease
to be Clifford. It identifies the missing operation in this first
moving-key word: a key-independent Pauli curvature with the complete
nonlinear character label.

## 5. Paid words, approximation and return

Use the common exact Q ledger from the linked query theorem. The source
$`S_3=RMQ`$ in Sections 2--3 and its conjugator with an insertion
costing $`G_P`$ Cliffords satisfy

```math
\begin{aligned}
T_{S_3}&\le3T_Q,&G_{S_3}&\le3G_Q+2w,&d_{S_3}&\le120,\\
T_{J_3}&\le6T_Q,&G_{J_3}&\le6G_Q+4w+G_P,&d_{J_3}&\le240.
\end{aligned}
```

The moving-key source $`S_2=MQ`$ and (5) satisfy

```math
\begin{aligned}
T_{S_2}&\le2T_Q,&G_{S_2}&\le2G_Q+2w,&d_{S_2}&\le80,\\
T_{J_2}&\le4T_Q,&G_{J_2}&\le4G_Q+4w+G_P,&d_{J_2}&\le160.
\end{aligned}
```

For singleton A translation $`G_P=b`$. General affine Pauli labels
are paid by their actual CNOT/CZ/Pauli word; a dense specification on
these two banks costs $`O(bw)`$ elementary Cliffords. Their controlled
phases are not supplied as arbitrary free address-dependent unitaries.

The same Q-private pool is reused sequentially. After either physical
alias assignment, the private width is at most $`3m-2+h_Q`$, plus
any explicitly added spectator wires. Every actual inverse is the
literal reversed word and has the same bounds. Q helpers return
exactly; A and B need not return in J. The completed conjugated probe
is the actual operator tested above. A failed observable identity is
not repaired by the fact that this operator is unitary.

For any general source K in Sections 2--3, the added common source
cost is one R call, and the conjugator costs
$`2(T_K+T_Q)`$, $`2(G_K+G_Q)+G_P`$, and
$`2(d_K+d_Q)`$. A common padded private pool of size
$`\max\{h_K,h_Q\}`$ suffices: the exact R word returns it on
every arbitrary input, including any current quantum fiber state left
by K. The corresponding physical alias width is
$`3m-2+\max\{h_K,h_Q\}`$. The fiber representation is not a
compiler for separate C_v branches.

If an emitted source is within eta of one of the exact screened
sources, its conjugator differs by at most $`2\eta`$ and its
transported observable by at most $`4\eta`$. Subtract $`4\eta`$
from the appropriate floor. If the completed conjugator itself is
approximated within eta, subtract $`2\eta`$. These are complete
operator estimates and include arbitrary references. They do not
assume intermediate helper return for approximate source calls.

## 6. A further mixer is outside the screens and still fails one test

The source

```math
S_4=QRMQ,\qquad J_4=S_4X_A^{\iota(x)}S_4^\dagger
```

is not of Section 2's final-R form. Indeed
$`RS_4=(RQR)(MQ)`$: MQ preserves the B Fourier key, whereas RQR
exchanges $`X_A^t`$ and $`X_B^t`$, using
$`k(a+t)=k(a)+t`$. Hence RS_4 does not preserve that key.
The conjugator satisfies $`J_4=QJ_3Q`$.
Its exact kernel is obtained without a clean label from the
Section 4 observable-source fiber formula. For $`J_0=MQX_A^uQM`$,
write

```math
J_{0,v}=(-1)^{\gamma_v(u)}X_A^{a_v(u)}Z_A^{z_v(u)},
\quad a_v=(I+\mathcal B_v^2)u,\quad z_v=\mathcal B_vu,
\quad\gamma_v=q_v(u)+q_v(\mathcal B_vu).
```

With $`d=2^w`$, the exact $`J_3=RJ_0R`$ kernel on a basis input is

```math
J_3|a,\beta\rangle
 =\frac1d\sum_{v,\beta'}
 (-1)^{v\cdot(\beta+\beta')+\gamma_v(u)
                 +z_v(u)\cdot(a+k(\beta))}
 |a+k(\beta)+a_v(u)+k(\beta'),\beta'\rangle.                 \qquad(9)
```

Applying Q before and after (9) gives the full J_4 kernel. There is no
trace extraction or scalarization in (9).

A single exact small input already stops this particular word. At
$`r=3,b=2,x=e_1,a=\beta=0`$, both target signs are +1, whereas

```math
\langle J_4^\dagger Z_A^tJ_4\rangle=0,
\qquad \langle J_4^\dagger Z_B(t)J_4\rangle=5/16.
```

The [two single-column fixtures](../../tests/test_key_mixing_observable_screens.py)
evaluate these values with integer Walsh sums, with norm exactly one;
no floating-point threshold is used. The resulting
basis-vector defects are at least $`\sqrt2`$ and
$`\sqrt{11/8}`$, respectively. This is a finite falsifier for S_4,
not an asymptotic exclusion of further mixers.

Its source uses four Q calls and 2w Hadamards. Its conjugator has
$`T\le8T_Q`$, $`G\le8G_Q+4w+b`$, $`d\le320`$, and the same
private pool and alias width as above. All Q helpers return. The extra
query does not cancel the key-dependent A translations and phases on
the displayed legal input.

## 7. What a conjugated translation additionally commits to

There is a useful exact distinction between an arbitrary observable
source and the repeatedly tested translation conjugator. Suppose
$`J_x=SX^{\iota(x)}S^\dagger`$ with x-independent S exactly meets
$`J_x^\dagger Z_{Y_0}J_x=Z_{Y_x}\otimes I_H`$ for every x.
The representation law $`J_uJ_x=J_{u+x}`$ then forces

```math
J_u^\dagger Z_{Y_x}J_u=Z_{Y_{x+u}}\otimes I_H
\quad\text{for every }u,x.                                \qquad(10)
```

Thus every fixed-u branch normalizes the full output diagonal algebra
and induces the entire output-wire translation permutation. Its action
within each computational output sector may still be an arbitrary
quantum helper unitary, with phases retained. This is an implication
of the full observable family and group law, not an added source
promise and not a reason to remove helper fibers.

Several noncommuting address insertions, intermediate address motion,
or a direct whole observable word need not have this representation
law. They remain concrete ways to leave the tested conjugator grammar.
Their complete character, native resource and actual-return identities
are still required.

## 8. Two noncommuting insertions and a protected remaining shift

A first two-insertion word can also be assessed as a complete operator.
Define

```math
L_A(u)=MQX_A^uQM,\qquad
L_B(v)=NRX_B^vRN,\qquad N=H_BRH_B,
\qquad J_x=L_A(u(x))L_B(v(x)).
```

The factors need not commute, so this J_x need not obey the translation
representation law in Section 7. Nevertheless, for affine u(x), its
top-B readout still has a protected shift. The exact Pauli-fiber formula
for L_A gives

```math
L_A(u)^\dagger Z_B(t)L_A(u)
 =\sum_v P_{v,t}(u)\otimes|v+t\rangle\langle v|_B,
```

with A Z label $`\mathcal B_tu`$ in every Pauli P. After also
Fourier-transforming A, this entire operator has the single A-key
shift $`\mathcal B_tu`$. The other labels and their literal phases
do not change that shift.

L_B preserves the A Fourier key. Conjugation by any key-preserving
unitary preserves an operator's specified key shift: its block from a
to a+s is merely multiplied on the left and right by the two unitary
fibers. This remains true when those fibers depend arbitrarily on x
and act quantum mechanically on all other banks and helpers. Therefore
the completed top-B readout of J_x has A-key shift
$`\mathcal B_tu(x)`$, while the required $`T_x^{(B)}`$ has shift
c(x). Nonaffinity of c gives an address where these shifts differ,
and the two unitary operators are then trace-orthogonal. Hence

```math
\bigl\|J^\dagger Z_B(t)J-T^{(B)}\bigr\|_{\rm op}
\ge\sqrt2\qquad(b\ge2).                                  \qquad(11)
```

This screen permits any A-key-preserving second correction, not only
L_B. It does not permit a correction that changes that remaining key
algebra. The displayed two-factor word has
$`T\le8T_Q`$, $`G\le8G_Q+8w+G_{P,u}+G_{P,v}`$,
$`d\le320`$, and the same sequentially reused query-private pool.
Its inverse is the actual reverse. For the top-A assignment,
$`L_A(u)^\dagger Z_A^tL_A(u)=(-1)^{u_{\rm top}}Z_A^t`$.
Thus L_A commutes with that probe for singleton u; a general affine u
retains the displayed address sign. In either case the remaining L_B
readout has A-key shift t, while the target has shift t+c(x). The
literal sign cannot repair that mismatch, so this assignment does not
supply an alternative successful alias.

A correction that changes the remaining key can be tested next, for
example

```math
J_x=L_A(u(x))R L_A(v(x))R.                                  \qquad(12)
```

This is outside the key-preserving-correction premise of (11). With
$`u=v=\iota(x)`$, it is an actual two-insertion commutator. Its
common ledger is $`T\le10T_Q`$,
$`G\le10G_Q+8w+2b`$, $`d\le400`$, with the same private pool
and physical alias allocation. At the single legal input
$`r=3,b=2,x=e_1,a=\beta=0`$, exact integer Walsh sums give norm
one and both top-probe expectations zero. Each corresponding target
sign is +1, giving a $`\sqrt2`$ basis-vector defect. This finite
test stops (12) with these parameters and probes; it does not exclude
other multiple-insertion words.

The remaining constructive obligation is a different complete
key-mixing correction or a whole-word interference identity. No
unobserved bank is treated as a clean label in these calculations, and
no depth bound is inferred for mechanisms outside the stated forms.
