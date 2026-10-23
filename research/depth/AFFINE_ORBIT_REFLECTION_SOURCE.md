# Affine kernel orbits and the missing reflection phase

[Canonical pencil](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md) ·
[Squared-product interface](COUPLED_ORDER_GRAM_SOURCE.md#squared-product-interface) ·
[Complete scalar zero test](EXPANDER_ZERO_REFLECTION.md)

The corank-one pencil supplies an exact, fixed-point-free native
reflection preserving the desired lookup charge. Its full orbit
normalizer is also paid at constant T-depth. This gives a concrete
partial reflection pair, not the missing complete indicator. A local
second reflection has maximal Gram error at every nonzero address.
All diagonal Boolean completions can be described exactly, and their
minimum phase degree is p+2. That degree statement is not a native
T-depth lower bound: a separate shallow scalar-zero trial below has
growing phase degree but still selects the wrong sector.

## 1. A fixed-point-free charge-preserving native reflection

Use p at least one, m=2^p, n=pm and the standard dual carrier banks z,t.
Let e=e_0 denote the first binary coordinate, not the first field block.
The [paid kernel query](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md#1-registers-and-the-actual-phase)
has $`K_x=P+B_x`$, where P deletes e and
$`B_x=\bar N_x^{\mathsf T}-I`$. In particular

```math
K_xe=0,\qquad v_x^{\mathsf T}K_x=0,\qquad v_x^{\mathsf T}e=1.
```

Write $`q_x(z,t)=v_x^{\mathsf T}(z+t)`$ and
$`Q_x|z,t\rangle=(-1)^{q_x(z,t)}|z,t\rangle`$. Define the actual
permutation reflection

```math
S_x|z,t\rangle=|z+e,t+K_xz+e\rangle.
                                                               \qquad(1)
```

Its chronological native word is $`Q_{B_x},Q_P,X_{z_0},X_{t_0}`$,
with every X emitted as HSSH. Applying (1) twice adds K_xe=0, so
S is a Hermitian involution on the complete space. It has no fixed
basis column because z_0 flips. Moreover q changes by
$`v_x^{\mathsf T}e+v_x^{\mathsf T}K_xz+v_x^{\mathsf T}e=0`$.
Thus S preserves both lookup sectors exactly, including coherent
address inputs. This is stronger than passing only the x=0 sector.

With $`r=(m-1)p^2`$, $`s\le(m-1)p^3`$ and
$`G_B\le2n+16r+2s`$, safe common forward/actual-reverse bounds are

```math
T_S\le8r,\quad d_S\le8,\quad
G_S\le G_B+n+7,\quad h_S\le r.                            \qquad(2)
```

There is no clean input and no precision error. The query's private
helpers return on arbitrary inputs. The actual reverse reverses the
native word; Hermiticity is a proved complete action, not an assumed
permission to omit that reverse.

## 2. The orbit coordinates have a paid cubic normalizer

Put $`a_x(z)=K_xz+e`$ and define

```math
C_x|z,t\rangle=|z,t+z_0a_x(z)\rangle.
                                                               \qquad(3)
```

C is an involution because x,z do not change. Since a is independent
of z_0, direct substitution gives

```math
C_xS_xC_x=X_{z_0},\qquad C_xQ_xC_x=Z_{z_0}Q_x.            \qquad(4)
```

Thus the invariant orbit coordinate is
$`u=t+z_0a_x(z)`$, together with x and all of z except z_0. The
remaining lookup charge in these coordinates is
$`v_x^{\mathsf T}(z+z_0e+u)`$; it is not a supplied local bit.

For completeness, (3) is an actual constant-depth native family.
Its terms are one CNOT from z_0 to t_0, n-1 Toffolis for z_0Pz,
and at most s three-control NOTs for the binary monomials of z_0B_xz.
The latter can be emitted with two arbitrary private bits d,g. For
logical controls a,b,c and target y, emit chronologically

```math
H_y,\quad T_{ab;d},\quad \mathrm{CCZ}(d,c,y),\quad
T_{ab;d},\quad \mathrm{CCZ}(d,c,y),\quad H_y.             \qquad(5)
```

Here $`T_{ab;d}`$ is the eight-T private-phase Toffoli, and the CCZ
is its eight-T phase word with the two target Hadamards removed;
both use g as their phase helper. The middle phase is
$`(-1)^{(d+ab)cy+dcy}=(-1)^{abcy}`$, with d,g returned.
Hence (5) is exactly C3X, with 32 T gates and 38 Cliffords.
Give each binary monomial its own d,g. The phase words have T gates
only on the private g wires; the Toffoli targets d are also private.
The outer H gates can be shared once per target y. Thus all such
monomials use 32 T layers even when their logical controls or target
bits coincide. The n-1 controlled-P terms use a separate eight-layer
phase wall. Helpers are reused only after complete return. Therefore

```math
\begin{aligned}
T_C&\le8(n-1)+32s,&d_C&\le40,\\
G_C&\le10(n-1)+38s+1,&h_C&\le\max(n-1,2s).
\end{aligned}                                                   \qquad(6)
```

Both (2) and (6) fit the polynomial intermediate source budgets.
All controls, banks and private helpers may be entangled with a reference.

<a id="diagonal-orbit-cochain"></a>

## 3. Exact diagonal completions and their sharp phase degree

Let R be a diagonal Boolean reflection with phase
$`R_x|z,t\rangle=(-1)^{f_x(z,t)}|z,t\rangle`$. Then

```math
(R_xS_x)^2|z,t\rangle
=(-1)^{f_x(z,t)+f_x(S_x(z,t))}|z,t\rangle.                \qquad(7)
```

The exact completion condition is the orbit difference equation
$`f+f\circ S=q`$. It has all the solutions

```math
f_x(z,t)=z_0q_x(z,t)+h_x\bigl(z+z_0e,\ t+z_0a_x(z)\bigr),
                                                               \qquad(8)
```

where h is arbitrary on the section z_0=0. Indeed the second term
is invariant under S, while q is invariant and z_0 flips. Conversely
subtracting z_0q from any solution leaves an invariant function, which
is determined by its value on that section. This proves existence,
but the displayed z_0q term already contains the unpaid lookup phase.

There is no choice of h making the global Boolean phase have bounded
degree in p. Precisely, the minimum algebraic-normal-form degree in
all address and carrier bits among solutions of (8) is **p+2**.
For the lower bound, restrict z to $`z=\lambda e`$. On this subspace
S is the fixed affine translation
$`(\lambda,t)\mapsto(\lambda+1,t+e)`$. A finite difference along
a nonzero fixed direction reduces Boolean degree by at least one.
The required difference is

```math
q_x(\lambda e,t)=\lambda+v_x^{\mathsf T}t,
```

which has degree p+1. Indeed the first binary coordinate in the last
field block of v_x is $`x^{m-1}=1[x\ne0]`$, whose degree is p;
its product with the corresponding independent t bit cannot cancel.
The upper bound is attained by f=z_0q, since every coordinate of v_x
is a function of p address bits. This establishes the exact minimum.

For exact Boolean R, (7) and Q have signs only. Any failed column
therefore gives direct operator error 2 and Gram error 1/2 for the
order-interference source. Consequently phase degree below p+2 cannot
even give Gram error strictly below 1/2 in this exact diagonal class.
This is **not** a bound on general native T-depth, nor a restriction
on the degree of an ideal phase approximated by a shallow unitary.

## 4. A paid local second reflection and a complete p=1 repair

The actual local trial

```math
R_{\rm loc}=Z_{z_0}\mathrm{CZ}(z_0,t_0),\qquad
f_{\rm loc}=z_0+z_0t_0
```

costs five elementary Cliffords. Let $`k_0=(K_xz)_0`$. Substitution
in (7), with all arithmetic binary, gives the complete square

```math
(R_{\rm loc}S_x)^2|z,t\rangle
=(-1)^{z_0+t_0+(z_0+1)k_0}|z,t\rangle.                  \qquad(9)
```

At x=0, k_0=0, so the entire correct base charge is obtained. For every
nonzero x, take z=0 and $`t=e_{n-p}`$. Then (9) has sign +1, whereas
Q has sign -1 because the last field block of v_x is the field unit.
The Gram gap is exactly 1/2 and the direct error is exactly 2, for
every p at least one. The column z=0, $`t=e+e_{n-p}`$ instead has
the opposite error, rejecting a correct-charge input. Both reflections
preserve Q; sector preservation
alone does not supply the required commutation/anticommutation pattern.
The direct trial has

```math
T_W\le16r,\quad d_W\le16,\quad h_W\le r,\quad
G_W\le2G_B+2n+24.                                      \qquad(10)
```

The complete-indicator envelope adds $`4c_n+2n+2m`$ Cliffords and
uses dirty width $`2n-m+r`$, where $`c_n=n^2+3n`$; it does not fix
the sector error. A literal two-clean order source is paid as well:
control S by $`C_x\mathrm{CX}(a,z_0)C_x`$, and control the local R
by CZ(a,z_0) followed by the eight-T private CCZ(a,z_0,t_0).
Before the fixed physical-basis envelope, the joint source has
$`T_U\le2T_C+16`$, $`G_U\le2G_C+45`$, $`d_U\le96`$ and
$`h_U\le\max(h_C,1)`$. Its complete accepted Gram is the one
certified by the squared-product interface, with the error above.

For p=1 there is an actual native correction, not only an abstract
solution: use

```math
R_1=Z_{z_0}\mathrm{CZ}(z_0,t_0)
\mathrm{CCZ}(x,z_0,z_1)\mathrm{CCZ}(x,z_0,t_1).           \qquad(11)
```

Its phase is z_0q, so its squared product with S is the full Q for
both addresses and every bank column. The two private-phase CCZs can
share their logical controls but use two distinct phase helpers:
$`T_R=16`$, $`G_R=21`$, $`d_R\le8`$, $`h_R=2`$.
At this size S has eight T gates and 19 Cliffords. Thus the direct
four-reflection word has 48 T gates, 80 Cliffords, depth at most 32
and two arbitrary returned helpers, before the fixed indicator envelope.
The [complete-indicator envelope](COUPLED_ORDER_GRAM_SOURCE.md#squared-product-interface)
adds 48 Cliffords at n=m=2. The resulting exact zero-clean p=1 indicator
therefore has T-count 48, Clifford count at most 128, T-depth at most 32
and four dirty bits outside its address and output. The p=1 indicator
already has a Clifford implementation; this is an identity fixture,
not a resource improvement.
This finite completed instance does not supply a uniform p-family of R.

<a id="scalar-zero-orbit-trial"></a>

## 5. Scalar zero gives a shallow large-degree phase, but a thin charge

Let $`Z_0(t)=(-1)^{[t=0]}`$ be the complete zero reflection, and
let $`Z_0(z_0,t)=(-1)^{[z_0=0,t=0]}`$. The exact ideal word

```math
R_{\rm zero}=C_x Z_0(t)Z_0(z_0,t)C_x
```

has phase $`z_0[u=0]`$, with the invariant u from (3). It is a
Hermitian involution, and the full squared product is exactly

```math
(R_{\rm zero}S_x)^2=C_x Z_0(t)C_x.                    \qquad(12)
```

This is a thin zero-test charge, not the lookup half-space. For every
address, z=0,t=e gives sign +1 in (12) and sign -1 under Q. Hence
the ideal direct error is 2 and its order-source Gram gap is 1/2.

Nevertheless this R has a paid shallow approximate native family.
Implement the two complete zero reflections to error epsilon each,
using their actual native words, and enclose them by the exact C.
Then $`\|\widehat R_{\rm zero}-R_{\rm zero}\|\le2\varepsilon`$;
the direct word using two such R calls and two exact S calls is within
4epsilon of (12), so its lookup error is at least $`2-4\varepsilon`$.
If the zero-test ledgers at k bits are $`T_0(k),G_0(k),d_0(k),h_0(k)`$
at the requested epsilon, the common R ledger is

```math
\begin{aligned}
T_R&\le2T_C+T_0(n)+T_0(n+1),\\
G_R&\le2G_C+G_0(n)+G_0(n+1),\\
d_R&\le80+d_0(n)+d_0(n+1),\\
h_R&\le\max\{h_C,h_0(n),h_0(n+1)\}.
\end{aligned}                                                   \qquad(13)
```

The two zero tests are sequential and use the same arbitrary private
pool. Uniform full-operator bounds justify telescoping despite approximate
return, including reference systems. Actual inverses reverse the
implementations; the approximate words are not declared exactly
Hermitian, and the preceding order-source Gram statement concerns the
ideal reflections. The direct
ledger is twice (13) plus twice (2). It has constant fixed-accuracy
T-depth and the permitted polynomial count/width overhead. Thus the
failure is the full sector action, not a degree-based depth argument.

## 6. Scope and first unpaid operation

The new paid component is S and its orbit normalizer C. Equations
(7)--(8) isolate the missing second reflection without a postselection
or clean-state promise. The local and scalar-zero completions fail
analytically in the stated all-size scopes. A complete p=1 repair is
checked natively. The exact guards in
[`test_affine_orbit_reflection_source.py`](../../tests/test_affine_orbit_reflection_source.py)
also check arbitrary-helper return, actual reversal, field-coordinate
orbits and the sharp phase-degree claim at small sizes.

The complete-frame frontier is unchanged. Continue only with a second
reflection whose complete orbit difference is q, or a different interface.
Merely inserting the known indicator into f=z_0q is a circular completion,
not a depth improvement. Choosing a high-degree phase does not itself
solve or obstruct its shallow approximate synthesis. These are internal
analytic proofs and exact finite checks, not external mathematical review.
