# Complete observable transport and the shared-key source limit

The first unpaid upper operation can be stated as transport of one
observable. The general two-call principle is already retained in
[joint routing, Sections 1--3](JOINT_ROUTING_INDICATOR_INTERFACES.md),
including the source-only observable promise, and
[weaker indicator interfaces](WEAKER_INDICATOR_INTERFACES.md) already
permits arbitrary address and helper gauges. Here it specializes the
current nonlinear-character source: no in-place coefficient
permutation or prescribed helper fibers are needed. This reformulation
does not by itself make the native target asymptotically easier; the
retained selected-swap and indicator reductions are equivalent up to
constant source calls and linear overhead.

The paid quadratic inverse-XOR and its partial Fourier transform do not
meet this interface through the retained translation-conjugator grammar.
In fact, arbitrarily long words that preserve their common Fourier key
have Clifford fibers, and their relative Pauli labels remain affine in
the address. The resulting obstruction is to that grammar, not to
unrestricted native depth. The general resource frontier is unchanged.

## 1. Specializing the retained observable interface

Let x be a b-bit address, let $`m=2^b`$, and let Y be the m-bit output
bank, with one wire $`Y_x`$ for each address. All address, output and
private-work inputs are arbitrary, including their external references.
For the ledger in the current coordinates, define

```math
\Omega_b=\sum_x|x\rangle\langle x|\otimes Z_{Y_x}\otimes I_H.
```

Suppose a terminating native compiler supplies a zero-clean unitary J
on the complete space, with common forward/actual-adjoint bounds
$`T_J,G_J,d_J,h_J`$, such that

```math
\bigl\|J^\dagger Z_{Y_\varnothing}J-\Omega_b\bigr\|_{\mathrm{op}}
\le\delta.                                                     \qquad(1)
```

The actual unitary J may move x, Y and H. There is no separate promise
about its basis columns or about helper return after one J call. The
identity operator on H in the target of (1) is essential.

Execute the chronological word

```math
H_Y,\qquad J,\qquad Z_{Y_\varnothing},\qquad J^\dagger,\qquad H_Y.
```

Since $`H_Y\Omega_bH_Y=I_b\otimes I_H`$, its emitted unitary E obeys

```math
\begin{aligned}
\|E-I_b\otimes I_H\|_{\mathrm{op}}&\le\delta,\\
T_E&\le2T_J,\\
G_E&\le2G_J+2m+2,\\
D_T(E)&\le2d_J,\\
h_E&\le h_J,\qquad a_{\mathrm{clean}}=0.
\end{aligned}                                                  \qquad(2)
```

The displayed Clifford allowance includes the $`2m`$ output Hadamards
and the one probe Z. Reversing this actual word has the same bounds.
Equation (1), after unitary conjugation, is already the complete error
and return estimate; no additional leakage estimate is required. It
holds on arbitrary entangled inputs because it is a full-operator norm
bound. Exact equality in (1) gives exact return of every private wire.

If the source is instead specified by
$`\|J-J_0\|_{\mathrm{op}}\le\eta`$ for an exact solution $`J_0`$,
unitary telescoping gives $`\delta\le2\eta`$. Approximation of one
prepared state or of one compressed block does not imply (1).

## 2. The zeta character fits the same interface

The square-zero character from
[output-aliased transport](OUTPUT_ALIASED_CHARACTER_SOURCE.md) uses
$`r=b+1`$, $`w=2^r-1=2m-1`$, and

```math
\begin{aligned}
R_r&=\mathbb F_2[e_1,\ldots,e_r]/(e_1^2,\ldots,e_r^2),\\
\iota(x)&=\sum_{i=1}^b x_i e_i,\\
A_x&=\prod_{i=1}^b(1+x_i e_i)-1,\\
\chi(x,y)&=[e_{[r]}](yA_x).
\end{aligned}
```

After the binary subset-zeta transform on Y, write
$`z_S(Y)=\bigoplus_{u\subseteq S}Y_u`$. Assign the actual output wire
holding $`z_S(Y)`$ to $`y_{[r]\setminus S}`$ for each nonempty
$`S\subseteq[b]`$. Assign the untouched wire
$`Y_\varnothing=z_\varnothing(Y)`$ to the top carrier coordinate
$`y_{[r]}`$. The remaining $`m-1`$ carrier coordinates are arbitrary
private inputs. No coordinate is initialized or copied.

Then

```math
\chi(x,y)
=\sum_{\varnothing\ne S\subseteq[b]}
  z_S(Y)\prod_{i\in S}x_i
=Y_x+Y_\varnothing.
```

Consequently the sufficient character condition is

```math
\bigl\|J^\dagger Z_{[r]}J-Z_{[r]}D_\chi\bigr\|_{\mathrm{op}}
\le\delta,\qquad
D_\chi|x,y,h\rangle=(-1)^{\chi(x,y)}|x,y,h\rangle.              \qquad(3)
```

After the stated wire assignment, $`Z_{[r]}D_\chi`$ is exactly the
addressed observable $`(-1)^{Y_x}`$ in zeta coordinates. Thus the
chronological word is

```math
H_Y,\quad\mathrm{zeta}_Y,\quad
J,\quad Z_{[r]},\quad J^\dagger,\quad
\mathrm{zeta}_Y^\dagger,\quad H_Y.
```

There is no outer probe Z and no separate constant-character correction.
The ledger becomes

```math
\begin{aligned}
T_E&\le2T_J,&D_T(E)&\le2d_J,\\
G_E&\le2G_J+bm+2m+2,&h_E&\le h_J,\\
\varepsilon_E&\le\delta.
\end{aligned}                                                \qquad(4)
```

Here $`h_J`$ counts every private wire outside x and Y, including any
unassigned carrier coordinates. The two actual zeta words cost bm
CNOTs. This formulation does not silently supply an extra carrier bank.

For each fixed x, $`Z_{[r]}D_\chi`$ is a Pauli. Its binary Z label is
not affine in x when $`b\ge2`$: the coefficient of
$`y_{[r]\setminus\{1,2\}}`$ is $`x_1x_2`$. This elementary label
fact, rather than the degree of an unimplemented permutation, is what
the following source test uses.

## 3. A shared Fourier key stops every word in the selected grammar

Use arbitrary banks A and B. In a fixed Fourier basis of B, suppose an
exact source has the full-space form

```math
S=\sum_v C_v\otimes|v\rangle\langle v|_B,                    \qquad(5)
```

where each $`C_v`$ is an actual Clifford on A. A may include arbitrary
dirty work; there is no prepared-state premise. No uniform circuit
description for the separate $`C_v`$ is assumed or supplied for free.
Equation (5) is used only to analyze an already paid actual word.

Let $`P_x`$ be a Pauli on A whose binary label $`\ell(x)`$ is affine
in the read-only address x, and define

```math
J_x=S P_x S^\dagger.
```

An address-controlled diagonal Pauli on the Fourier key may also
multiply $`P_x`$; it changes fiber phases but not the following labels.
In each fiber,

```math
J_{x,v}=C_vP_xC_v^\dagger,\qquad
\ell_v(x)=F_v\ell(x),                                        \qquad(6)
```

where $`F_v`$ is the binary symplectic action of $`C_v`$. Thus
$`J_{x,v}`$ is a Pauli with label affine in x. Its literal phase has
not been removed.

Take the computational-basis probe $`Z_B(t)`$, which shifts the Fourier
key by t. The transported probe has the exact block form

```math
J_x^\dagger Z_B(t)J_x
=\sum_v J_{x,v+t}^\dagger J_{x,v}
       \otimes|v+t\rangle\langle v|.                         \qquad(7)
```

Each relative block is a Pauli, possibly with phase $`\pm i`$, and its
binary label is

```math
\ell_{v+t}(x)+\ell_v(x),                                     \qquad(8)
```

which remains affine in x. This holds for every number of letters in S.

Suppose the required observable is
$`Z_B(t)Z_A^{c(x)}`$, where the Pauli label $`c(x)`$ is nonaffine; in
particular it may contain the character coordinate $`x_1x_2`$ above.
Fix any Fourier key v. Equations (6)--(8) imply that some address x has
a relative-block label different from the label of $`Z_A^{c(x)}`$.
Distinct Pauli labels are trace-orthogonal even with arbitrary literal
phases. If d is the A-space dimension, this gives

```math
\frac1d\mathrm{tr}\bigl[(R-Z_A^{c(x)})^\dagger
                         (R-Z_A^{c(x)})\bigr]=2,
\qquad
\|R-Z_A^{c(x)}\|_{\mathrm{op}}\ge\sqrt2.                     \qquad(9)
```

The full operators in (7) and its target have the same key shift, so
their difference has norm equal to the maximum norm of its blocks.
Taking also the maximum over the read-only address proves

```math
\bigl\|J^\dagger Z_B(t)J-Z_B(t)Z_A^{c(x)}\bigr\|_{\mathrm{op}}
\ge\sqrt2.                                                   \qquad(10)
```

The same argument allows a fixed Pauli factor on A in the probe: its
fixed label is merely added to (8). Clifford changes of output
coordinates do not make a nonaffine target label affine.

This is a robust defect of the completed observable, not only an exact
failure. If an emitted source $`\widetilde S`$ is within eta of an
ideal S of the form (5), the conjugator differs by at most $`2\eta`$,
and its transported probe differs by at most $`4\eta`$. Therefore the
right side of (10) can be replaced by $`\sqrt2-4\eta`$. If the
completed conjugator itself is within eta of the grammatical one, the
corresponding bound is $`\sqrt2-2\eta`$.

### Why the paid quadratic kernels satisfy the premise

Let k be the quadratic inverse coefficient map in the
[paid full-bank query](FULL_BANK_CHARACTER_TRANSPORT.md), and put

```math
Q=Q_k(A;B),\qquad M=H_AQH_A.
```

In the B-Fourier fiber v,

```math
Q_v=D_{q_v},\qquad M_v=H_AD_{q_v}H_A,\qquad
q_v(a)=v\mathbin\cdot k(a).
```

The Boolean polynomial $`q_v`$ is quadratic, so both displayed fiber
operators are Cliffords. Every fixed word in Q, M, and fixed Cliffords
on A preserves the same Fourier key and has a Clifford fiber. The
actual paid Q returns its private work exactly, so that work does not
invalidate this full-space statement. All word lengths, including
growing lengths, are covered by (10) when the address enters through
the stated affine-Pauli conjugator and the readout has the stated
relative-block target.

The theorem does not cover a word that mixes this Fourier key, a
coherent translation of the key itself, non-Clifford fiber operations,
or intermediate address motion. It does not assert that all circuits
with a given native T-depth have the form (5). Increasing the word
length within this grammar cannot fix its label problem; changing the
grammar may.

## 4. The first reverse overlap is paid and still has affine labels

Use operator order, with the rightmost factor acting first, and take

```math
S=MQ=H_AQH_AQ.
```

Retain all linear terms of $`q_v`$, and define its polar matrix by

```math
q_v(a+u)+q_v(a)+q_v(u)=a\mathbin\cdot B_vu.
```

The matrix $`B_v`$ is symmetric and alternating, and depends linearly
on v. In the Fourier-key fiber,

```math
S_v=H_AD_{q_v}H_AD_{q_v},
```

and the actual translation conjugator has the exact Pauli fiber

```math
\begin{aligned}
J_{u,v}=S_vX_A^uS_v^\dagger
 &=(-1)^{\gamma_v(u)}X_A^{a_v(u)}Z_A^{z_v(u)},\\
\gamma_v(u)&=q_v(u)+q_v(B_vu),\\
a_v(u)&=(I+B_v^2)u,\qquad z_v(u)=B_vu.
\end{aligned}                                                  \qquad(11)
```

For completeness, apply
$`D_qX^uD_q=(-1)^{q(u)}X^uZ^{Bu}`$ twice, with the intervening
Hadamards. The possible reordering signs vanish because
$`u\cdot B_vu=0`$ and
$`(B_vu)\cdot B_v^2u=0`$. Thus (11) retains the literal phase.

For the top-B probe with shift t, the relative block is exactly
$`J_{u,v+t}^\dagger J_{u,v}`$. Its labels are

```math
\begin{aligned}
\Delta z&=B_tu,\\
\Delta a&=(B_vB_t+B_tB_v+B_t^2)u.
\end{aligned}                                                  \qquad(12)
```

Both remain linear in u. For $`u=\iota(x)`$, they cannot equal the
required nonaffine character labels when $`b\ge2`$. This remains true
even if a proposed cancellation removes the X residue in (12). The
scalar phase in (11) cannot change a binary Pauli label.

The other ordering offers no shortcut for the original translation:
M commutes with all X translations, so
$`(QM)X^u(QM)^\dagger=QX^uQ`$. On the two computational banks this
last word acts as

```math
|a,z\rangle\longmapsto|a+u,z+k(u)+au\rangle,
```

using $`k(a+u)+k(a)=k(u)+au`$. With $`u=\iota(x)`$, top readout gives
the bilinear term $`[e_{[r]}](au)`$, not the required
$`[e_{[r]}](aA_x)`$.

These are actual paid words. With the common Q ledger from the linked
query theorem, S on its two declared arbitrary w-bit banks satisfies

```math
T_S\le2T_Q,\qquad G_S\le2G_Q+2w,\qquad
d_S\le80,\qquad h_S\le h_Q.
```

For $`J_x=SX_A^{\iota(x)}S^\dagger`$, the b address-controlled
translations are b CNOTs, giving

```math
T_J\le4T_Q,\qquad G_J\le4G_Q+4w+b,\qquad
d_J\le160,\qquad h_J\le h_Q.                                 \qquad(13)
```

In these two displays A and B are declared principal inputs. If B is
instead an added private bank, add w to the private width. If this
candidate is embedded in the character readout by placing the
$`m-1`$ nonconstant output characters in A and the constant output
$`Y_\varnothing`$ at the top of B, the two banks require
$`3m-2+h_Q`$ private wires in total. Every Q work wire returns exactly.
Equation (13) pays the failed candidate, not a successful indicator.

## 5. Complete normalization and the unchanged frontier

For every fixed requested delta, a qualifying observable source must
have simultaneous bounds, uniform in b,

```math
T_J,h_J=O_\delta(2^b(b+2)^c),\qquad
G_J=O_\delta(4^b(b+2)^{c_G}),\qquad
d_J(b)=o(\log(b+2)),                                         \qquad(14)
```

for fixed resource exponents. Equation (2) or (4) then supplies the
complete dirty indicator required by
[polynomial-overhead root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md).
The native source, its actual adjoint, every helper wire and the full
Clifford count must satisfy these bounds on the same circuit. A cheap
description of separate address branches is insufficient.

The retained fusion word uses eight indicator/actual-adjoint calls.
Substituting (2) gives sixteen J/actual-adjoint calls and

```math
\begin{aligned}
\varepsilon_{I_p}&\le8\delta,\qquad\delta=\tau/8,\\
D_T(I_p)&\le16d_J(b_0+1)+O_\tau(\log\log(p+4)).
\end{aligned}                                               \qquad(15)
```

where $`b_0+1=p-O_\tau(\log(p+2))`$ is the fused source address
length in the retained notation. The resulting normalized T-count and
dirty width are $`O_\tau(2^p)`$, and the Clifford count is
$`O_\tau(4^p)`$. Return and reference error are included in
$`8\delta`$; the exact tail adds no approximation error.

For the retained fixed-precision frame transfer, a word query uses at
most four normalized indicators. Its error is therefore at most
$`32\delta`$. With at most $`K_*`$ such word-query occurrences, choose

```math
\delta_L=\frac{2^{-L}}{128K_*}.
```

Their total added error is at most $`2^{-L}/4`$, within the retained
frame allocation. If J is instead approximated in unitary norm to an
exact ideal, its tolerance must account for $`\delta\le2\eta`$.
No accuracy improves automatically when the source width grows.

The first unpaid operation is a native J satisfying (1), or (3), with
the simultaneous sublogarithmic-depth budgets (14). Specializing the
retained observable interface removes unnecessary promises about this
candidate's complete source permutation and helper fibers. It is a
synthesis freedom, not a new asymptotically easier target. The shared-key
theorem explains why the currently paid quadratic/Fourier words still
do not supply it. Neither the reformulation nor this scoped obstruction
improves the unconditional general upper bound or proves a
dimension-growing unrestricted depth lower bound.

## 6. A genuine key change exposes the unpaid dirty correction

The opposite-direction two-shear word escapes the shared-key premise.
Take

```math
S=Q_k(B;A)Q_k(A;B),\qquad J_u=SX_A^uS^\dagger.
```

Use beta for the arbitrary computational value of B, to distinguish it
from the address length b. Direct substitution of the complete XOR
maps, with no initialized bank, gives

```math
\begin{aligned}
\alpha&=a+k(\beta),\qquad d=k(u)+\alpha u,\\
J_u:(a,\beta)&\longmapsto
\bigl(a+u+k(d)+\beta d,\ \beta+d\bigr).
\end{aligned}                                                  \qquad(16)
```

To verify (16), the inverse source first produces
$`(\alpha,\beta+k(\alpha))`$. The translation changes alpha to
$`\alpha+u`$. The forward source then uses
$`k(\alpha+u)+k(\alpha)=k(u)+\alpha u`$ and
$`k(\beta+d)+k(\beta)=k(d)+\beta d`$. These identities retain
every arbitrary-bank term.

Neither natural top-coordinate probe gives the desired character.
For the top-B probe, set $`a=0`$, $`u=e_i`$, and
$`\beta=e_{[r]\setminus\{i\}}`$. Then the top-B change
$`[e_{[r]}]d`$ is one, whereas the required
$`[e_{[r]}](aA_x)`$ is zero.

For the top-A probe, take $`b\ge2`$, $`u=e_i`$, $`a=e_r`$, and
$`\beta=e_{[b]\setminus\{i\}}`$. Here
$`k(\beta)=\beta`$ and
$`d=u(1+a+\beta)`$. Every monomial in d contains i, so no two can
contribute a disjoint product and $`k(d)=d`$. The top-A change is

```math
[e_{[r]}]\bigl(u+k(d)+\beta d\bigr)
=[e_{[r]}](au+a\beta u)=1,
```

while $`A_x=f(e_i)=e_i`$ and $`[e_{[r]}](aA_x)=0`$. Each example
is a literal sign disagreement on an allowed basis input, so its
observable defect is exactly two. The same conclusion holds with
arbitrary external references because the claimed full-operator
identity already fails on this input.

All private Q work returns exactly. The conjugated Z readout is
diagonal on both computational banks and returns those banks exactly
as well. Thus the failure is a wrong character caused by an arbitrary
coherent key, rather than a return failure. On two declared principal
banks, this source has $`T_S\le2T_Q`$, $`G_S\le2G_Q`$ and
$`d_S\le80`$; its translation conjugator has
$`T_J\le4T_Q`$, $`G_J\le4G_Q+b`$ and $`d_J\le160`$, with the
same reused $`h_Q`$ private pool. Added principal banks must still be
counted when embedded as private work.

This test concerns the actual word (16) and these two probes. It is
not an obstruction to other key-changing words or other readouts. Its
first unpaid correction is the cancellation of the displayed
$`k(\beta)u`$ and $`\beta d`$ terms without a clean key.

### An apparent key exchange can preserve the other key

For comparison, the two-kernel word $`S=MH_BMH_B`$ does not escape
Section 3 by merely exchanging the B-Fourier basis. In the A-Fourier
fiber a, its exact action is

```math
M_a=X_B^{k(a)},\qquad
(H_BMH_B)_a=Z_B^{k(a)},\qquad
S_a=X_B^{k(a)}Z_B^{k(a)}.
```

It therefore has Pauli fibers over the unchanged A-Fourier key.
Likewise $`M Q_k(B;A)`$ has the exact A-key fiber
$`X_B^{k(a)}D_{q_a}`$, where $`q_a(\beta)=a\cdot k(\beta)`$
is quadratic. These are Clifford fibers over the other key. Both
letters commute with every A translation, so an A insertion is trivial;
an affine B-Pauli insertion is covered by the swapped version of
Section 3. The opposite order $`Q_k(B;A)M`$ also cancels M around
a B translation. A visible Hadamard or changed query direction is
therefore insufficient to certify that the actual word leaves the
shared-key representation.

## 7. Next mechanism

The next constructive question is whether an actual full-input word
that changes the common key algebra can satisfy the observable identity
without computing every forward coefficient. A candidate must be
tested as a complete operator, including its actual inverse, literal
phases and arbitrary helper inputs. Merely renaming a register as a
clean Fourier label, or choosing a separate Clifford for each address,
does not pay that operation.

The same-key grammar has been assessed and can be left behind. The
first genuine key-changing word is also assessed: it pays the native
operation and return, but leaves an explicit dirty-key correction.
The observable interface remains available to other words with genuine
key mixing, address motion or non-Clifford fiber dynamics. Their native
resource and complete character identities remain open.
