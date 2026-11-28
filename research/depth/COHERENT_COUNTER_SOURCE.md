# Coherent counters: an exact odd-key code and the missing native schedule

[Complete source target](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Clean-source promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Paid affine phase batches](PARALLEL_AFFINE_DETUNING.md) ·
[Distributed counter and common shift](#distributed-cyclic-counter) ·
[Monomial clock repair](#monomial-clock-repair-bank) ·
[Summable approximate repair](#distance-weighted-clock-repair) ·
[Local companion bound](#local-companion-depth-count) ·
[MOD3 row and finite-angle target](#private-gauge-mod3)

A counter gives an exact complete indicator identity with logarithmically
many clean wires. A precise phase-gradient state is unnecessary: one
clean bit selects all odd Fourier keys of a modular register, and every
one of those keys gives the same indicator. All other bits of that
register may be arbitrary dirty inputs and return exactly.

This supplies a literal amplitude-$`1/\sqrt2`$ interface and identifies
its native operation. It does **not** supply a shallow source. The joint
modular accumulation and the full-bank schedule remain unpaid at the
required sublogarithmic depth. Isolated shallow rotations, a supplied
catalyst, or logical commutation do not pay either obligation.

The revision in Sections 9--10 separates this counter from a different
dirty-key mechanism. It proves a clean-dimension/error tradeoff for
one-pass row filters and pays an exact full-bank echo, whose emitted
native depth remains linear. Section 11 pays parallel local clock
visits and reduces exact row comparison to a joint repair. That
repair still costs logarithmic native T-depth in the available
realization. These results leave the general frontier unchanged.

## 1. Registers and the complete modular operation

For $`p\ge1`$, put

```math
m=2^p,\qquad k=\lceil\log_2(p+2)\rceil,\qquad Q=2^k>p+1.
```

The address x and m-bit bank Y are arbitrary. C and G are k-bit
registers, interpreted as integers modulo Q. C starts in zero. Only
the most significant bit g of G starts clean; the other k-1 bits are
arbitrary and may be entangled with x,Y and an external reference.
For a fixed bank index u define

```math
w_u(x,Y)=|x\mathbin\oplus u|+1-Y_u\in\{0,\ldots,p+1\},
\qquad
W_u|x,Y,j,z\rangle
 =|x,Y,j,z+jw_u(x,Y)\bmod Q\rangle.                 \qquad(1)
```

The operation W is defined on the **entire** C,G space, not only on
a prepared phase eigenstate. Any native implementation must also
return its private work on arbitrary inputs. Its actual inverse
subtracts the same product. Neither arithmetic encoding nor return
is supplied by writing (1).

Let

```math
H_C=H^{\otimes k},\qquad
A_u=H_CW_uH_C,\qquad
R_0=I_C-2|0^k\rangle\langle0^k|,\qquad
E_u=A_u^\dagger R_0A_u.                            \qquad(2)
```

All identities below include the literal scalar phase.

## 2. One clean bit selects every usable Fourier key

Define the Fourier eigenstates of modular addition by

```math
|\chi_t\rangle_G
 =Q^{-1/2}\sum_{z=0}^{Q-1}e^{-2\pi itz/Q}|z\rangle,
\qquad
\mathrm{ADD}_s|\chi_t\rangle
 =e^{2\pi its/Q}|\chi_t\rangle.                   \qquad(3)
```

Addition by Q/2 is exactly $`X_g`$. Thus

```math
| -\rangle_g\otimes\mathcal H_{G\setminus g}
 =\mathrm{span}\{|\chi_t\rangle:t\ \mathrm{odd}\}. \qquad(4)
```

This is equality of subspaces. No Fourier key is initialized, learned
or measured. On key t, W acts on C by the diagonal phase
$`e^{2\pi itjw_u/Q}`$. Hence

```math
\langle0_C|A_u|0_C\rangle
 =\frac1Q\sum_{j=0}^{Q-1}e^{2\pi itjw_u/Q}
 =[tw_u=0\pmod Q].                               \qquad(5)
```

Every odd t is invertible modulo Q. Since $`0\le w_u\lt Q`$, the
last expression equals $`[w_u=0]`$ on every sector in (4).
If w is zero, A is the identity. Otherwise A sends zero to a vector
orthogonal to zero. The reflection in (2) therefore gives

```math
E_u\bigl(|0_C\rangle\otimes|\xi\rangle\bigr)
 =|0_C\rangle\otimes(-1)^{[x=u]Y_u}|\xi\rangle   \qquad(6)
```

for every vector in the full odd-key code, including arbitrary
lower-G inputs and references. The Boolean phase is independent of t,
so coherences between keys return as well. Equation (6) is a complete
initialized-isometry identity, stronger than a scalar matrix element.

## 3. Full indicator and exact flat compression

Prepare g with $`P=ZH`$, which maps zero to minus, and use its
actual inverse to unprepare. For any fixed order of the bank indices,
set

```math
V=P^\dagger\Bigl(\prod_u E_u\Bigr)P,
\qquad \Omega|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle.
```

Let J insert $`C=0^k,g=0`$, leaving every other wire arbitrary.
Every ideal E returns the code, so (6) gives

```math
VJ=J(\Omega\otimes I_{\mathrm{private}}).          \qquad(7)
```

The lower G bits and every native helper return exactly. Applying H
to the same clean bit g after this completed query gives a source
$`U=H_gV`$ with

```math
K=J^\dagger UJ=\Omega/\sqrt2,\qquad K^\dagger K=I/2.
                                                               \qquad(8)
```

Thus the k+1 clean wires suffice for both the returned query and the
exactly flat compression; no additional amplitude signal is needed.
In the source U, g finishes in plus as expected, while C and all
private work return. Conjugating the Y bank by Hadamards converts the
phase indicator to the XOR indicator, adding 2m Cliffords.
The final H is only a paid normalization tail of the already returned
query. It does not make the source two-band, and no three-call recovery
is needed to recover a query that (7) already supplies.

## 4. Paid reflection, conditional ledger and actual inverses

The zero reflection requires no additional clean wire. Let C_k be the
retained [returned dirty zero-test](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence),
which toggles an arbitrary d by $`[C=0]`$. Then

```math
R_0\otimes I_{d,w_C}=C_k^\dagger Z_dC_kZ_d.        \qquad(9)
```

Use the actual C_k inverse. If its common forward/inverse bounds are
$`T_C,G_C,d_C,h_C`$, safe reflection bounds are

```math
T_R=2T_C,\qquad G_R=2G_C+12,\qquad d_R=2d_C,
\qquad h_R=1+h_C.                                \qquad(10)
```

The retained conjunction has polynomial k counts and dirty width and
$`d_C=O(\log(k+3))`$. The additive twelve covers actual inverse Z
words in the positive-S alphabet. The toggle and its helpers return on
all inputs.

Suppose each W has a zero-clean complete native implementation with
common forward/actual-inverse bounds $`T_W,G_W,d_W,h_W`$, uniform
in u. Every occurrence of an inverse in (2) is the actual native
inverse, not a new approximation. A row has bounds

```math
T_E\le2T_W+T_R,\quad
G_E\le2G_W+G_R+4k,\quad d_E\le2d_W+d_R.          \qquad(11)
```

The explicitly sequential source in (8) consequently has

```math
\begin{aligned}
T_U&\le m(2T_W+T_R),\\
G_U&\le m(2G_W+G_R+4k)+15,\\
d_U&\le m(2d_W+d_R),\\
k_{\mathrm{clean}}&=k+1,\\
h_{\mathrm{dirty}}&\le k-1+h_W+h_R.
\end{aligned}                                                   \qquad(12)
```

Separate helper pools give the safe width sum. The Clifford constant
covers P, its inverse, and the final H in both directions, emitting
$`Z=S^2`$ and charging the actual inverse S gates. These count bounds
fit the polynomial-overhead envelope if the W primitives do; the
displayed depth plainly does not qualify.

For actual primitive errors $`\eta_W,\eta_R`$ in full operator norm,
telescoping with the actual adjoints gives

```math
\delta_V\le m(2\eta_W+\eta_R),\qquad
\|\widetilde VJ-J\Omega\|\le\delta_V,
\qquad
\|J^\dagger\widetilde UJ-\Omega/\sqrt2\|
 \le\delta_V.                                    \qquad(13)
```

If the actual emitted word preserves the g code exactly, the last
bound improves to $`\delta_V/\sqrt2`$; this extra invariant is not
assumed in (13). The full
initialized-isometry error, including every return error and arbitrary
references, is at most $`\delta_V`$. Approximate W words do not
automatically retain (8)'s exact flatness. Exact native W words have
zero primitive error and all the exact return statements above.

## 5. What remains unpaid in joint arithmetic

Equation (1) can be written as p+1 controlled additions of C into G,
one for each literal in w. Each addition has finite exact native
implementations. For example, a controlled addition can be expanded
into controlled powers-of-two increments; descending target-bit order
expresses each increment using returned dirty conjunction queries.
This gives polynomial counts and width with no clean arithmetic
scratch. Its serial nonlinear depth does not qualify.

A reversible encoder of the ordinary Hamming weight also retains its
old nonlinear cost. The bit of weight $`2^r\le p`$ in
$`\sum_i x_i`$ is the Boolean elementary symmetric polynomial
of degree $`2^r`$. An XOR/Toffoli encoder with t Toffoli rounds has
degree at most $`2^t`$, even with arbitrary dirty work. Thus an
encoder explicitly producing all those bits needs logarithmically many
such rounds. This concerns that arithmetic representation and its
gate-by-gate schedule; it is not an unrestricted native T-depth bound.
Whole-word quantum resynthesis may leave this representation.

Giving all m rows their own initialized counters uses m k clean
counter wires, and serial reuse gives (12). The later
[shifted-sector promotion](PARALLEL_CLEAN_PROMOTION.md#random-shift-clean-promotion)
pays the complete bank schedule from a native-shareable, exactly
silent false-row source. The exact counter identity has the required
sector action for one row. To use this route, pay the row's modular
accumulation and actual reverse at sublogarithmic T-depth with
private T targets and read-only address boundaries. A generic
approximation of (1) does not preserve exact silence on false rows.

Nor may all W calls simply be merged before one readout. Since

```math
\sum_u w_u=mp/2+m-|Y|,                            \qquad(14)
```

their product on one counter/accumulator is independent of x. The
interleaved reflections are essential to the proven source identity.

The next question is a returned, shareable modular accumulation for
one exact-silent row, or a different complete counter echo. Either
needs a sublogarithmic native schedule with the stated simultaneous
count, precision and return bounds.

### Exact cancellation with a uniform phase seed

A uniform seed on k bits cannot replace modular accumulation by
arbitrarily low-degree binary phases with exact false-row cancellation.
Suppose each $`q_r:\mathbb F_2^p\to\mathbb F_2`$ has degree at most
$`d\ge1`$, with $`q_r(0)=0`$, and

```math
2^{-k}\sum_{r\in\mathbb F_2^k}(-1)^{q_r(x)}=[x=0].
```

Then $`k\ge\lceil p/d\rceil`$. Put $`s=\lceil p/d\rceil`$. The
Hamming weight of each $`q_r`$ is divisible by $`2^{s-1}`$. Indeed,
expand the parity of its monomial values as the integer sum over
nonempty selections of monomials. A selection of j monomials has
coefficient $`(-2)^{j-1}`$ and its product uses at most $`dj`$ variables.
Summing that term over x gives divisibility by
$`2^{j-1+\max(0,p-dj)}`$, which is a multiple of $`2^{s-1}`$.
Exact cancellation means half the seeds give one at each nonzero x,
so the total weight is $`(2^p-1)2^{k-1}`$. The odd factor gives the
claimed bound. It is sharp. Partition x into s blocks of at most d
bits and take $`q_r`$ to be the parity of the block ORs selected by r.
At any nonzero x some block OR is one, so exactly half the seeds give
each sign. The sharpness claim concerns the phase family.

For eighth-root phases $`\omega^{q_r(x)}`$, $`\omega=e^{i\pi/4}`$,
assume each of the three binary digits of $`q_r`$ has degree at most d.
The indicator of each phase value has degree at most $`3d`$, so its
weight is divisible by $`2^{\lceil p/(3d)\rceil-1}`$. Sum the exact
phase-average identity over x and expand in the independent basis
$`1,\omega,\omega^2,\omega^3`$, using $`\omega^4=-1`$. The coefficient
of 1 is $`2^k`$ and has the same divisibility. Thus
$`k\ge\lceil p/(3d)\rceil-1`$ in this variant.

Both bounds concern exact cancellation by literal low-degree phases
with uniform weights. In particular $`k,d=p^{o(1)}`$ cannot give this
exact-silent comparison. Approximate cancellation or native
interference outside this phase average calls for a separate analysis.

### Monomial SELECT with accepted work motion

Uniform SELECT can also move all the branch work. Consider k selector
bits prepared and recombined by Hadamards, and r additional initialized
bits tested for zero at the output. At each fixed address x the branch
indexed by the selector t is monomial in the computational basis. It
may implement any address-dependent permutation of the additional
clean and arbitrary dirty wires, with literal eighth-root phases.
The selector is unchanged inside SELECT. Let $`K_x`$ denote the
accepted compression, including both initialized tests.

Suppose $`K_0\ne0`$ and $`K_x=0`$ for every nonzero x.
If, at fixed computational work input and selector, output bits have
Boolean degree at most $`D_o`$ in the p address bits and the phase
exponent has integer multilinear degree at most $`D_\phi`$
modulo eight, then

```math
p\le rD_o+4k\max\{D_o,D_\phi\}.
```

Here all degrees may depend on the circuit family. In the native
monomial grammar of NOT, CNOT and Toffoli with standard diagonal
Clifford+T phases, d disjoint Toffoli rounds give the sharper
$`p\le(r+4k)2^d`$. Even false norm strictly below
$`2^{-4k}`$ forces exact false silence in this SELECT grammar.
Thus subpolynomially many initialized bits still need logarithmic
Toffoli depth in the branches. These statements allow branch work
motion and an arbitrary nonzero marked block.

Fix a computational dirty input with nonzero marked output and a
Walsh output character detecting that output. Multiply its matrix
element by the uniform selector denominator and the Walsh-vector
normalization denominator. The numerator is a sum of $`2^k`$
eighth roots or zero, with each summand the branch phase times that
output character and the r accepted-output literals. Exact false
silence makes the resulting function $`c[x=0]`$ for
$`c\ne0`$. The Walsh character is an added phase with exponent
four times a Boolean polynomial of degree at most $`D_o`$,
regarded modulo eight.

Use the prime $`\lambda=1-\omega`$ of $`\mathbb Z[\omega]`$.
Its norm is two, and 2 is a unit times $`\lambda^4`$
since the reduction of the eighth cyclotomic polynomial modulo 2
has one ramified prime. Equivalently, the four factors
$`1-\omega^j`$ for odd j multiply to 2 and each ratio
$`(1-\omega^j)/(1-\omega)`$ is an algebraic integer of norm one.
All four conjugates of c have modulus at most $`2^k`$,
so $`v_\lambda(c)\le4k`$.

Expand a phase polynomial into the product of its factors
$`1+(\omega^{a_S}-1)\prod_{i\in S}x_i`$.
Each nonconstant factor brings at most $`D_\phi`$ variables
and at least one $`\lambda`$. Include the Walsh phase with
degree $`D_o`$. For one accepted-output literal the integer
parity expansion of a Boolean degree-$`D_o`$ function
has terms from j monomials using at most $`jD_o`$ variables
with coefficient divisible by $`2^{j-1}`$. Hence
the r literals allow up to $`rD_o`$ variables without
a valuation cost. Every term that reaches the full address monomial
must otherwise pay at least one $`\lambda`$ per
$`\max(D_o,D_\phi)`$ variables. Sums preserve this
divisibility, and the full address coefficient of $`c[x=0]`$
is $`(-1)^pc`$. This proves the inequality.

For the native version, every computational label has degree at
most $`2^d`$. Instead of expanding the total phase exponent,
expand each physical phase as $`1+(\omega^j-1)b`$, where
b is the current Boolean label. A term from u monomials has
valuation at least $`1+4(u-1)\ge u`$. For CZ, use
$`1-2b_1 b_2`$ and expand the two labels, with the same
bound on the sum of the two monomial counts. This pays at
least one prime per $`2^d`$ variables outside the free
output-literal allowance. For the approximate claim, any nonzero
computational entry has a numerator with integer norm of magnitude
at least one and three other conjugates of size at most $`2^k`$.
After division by the selector denominator, that entry has modulus
at least $`2^{-4k}`$. A uniformly smaller block norm
therefore gives zero. Internal nonmonomial branch interference and
a joint readout without independent row acceptance remain open.

## 6. The whole product is a Hamming-shell observable

The E operators do commute, on the entire space and even outside the
odd-key code. To see the exact content of that commutation, fix x,Y
and a Fourier key t of G. Write $`|\psi_j\rangle=H_C|\chi_j\rangle`$
for an orthonormal basis of C. Then

```math
E_u=I-2|\psi_{tw_u}\rangle\langle\psi_{tw_u}|,
\qquad
\Bigl(\prod_u E_u\Bigr)|\psi_j\rangle
 =(-1)^{N_j}|\psi_j\rangle,
\qquad
N_j=\sum_u[j=tw_u\pmod Q]\pmod2.                 \qquad(14a)
```

Equal labels give equal projectors and unequal labels give orthogonal
ones. This proves commutation without moving any native gates.
On an odd key let $`a=t^{-1}j\pmod Q`$ use its representative
between zero and Q-1. Splitting according to Y gives

```math
N_j=\binom p{a-1}
 +\sum_uY_u\bigl([|x\mathbin\oplus u|=a]
                 +[|x\mathbin\oplus u|=a-1]\bigr)\pmod2,
                                                               \qquad(14b)
```

where out-of-range shells and binomial coefficients are zero.
The zero shell is $`Y_x`$. At Y=0 and a=1 the phase is instead
minus one, so an arbitrary counter cannot replace the zero counter.

Fixed Hamming-shell transforms of Y are binary linear maps and admit
Clifford XOR queries. Loading every shell into arbitrary dirty banks
does not supply the selected coordinate (a,x), or cancel those banks'
unknown offsets after readout. A Krawtchouk or Walsh representation
must pay that remaining operator. Formula (14b) is a concrete joint
target for a different native mechanism, not a parallel schedule.

## 7. Optional component: a fixed-accuracy gradient can also be prepared

This component is unnecessary for (4)--(8), but separates preparation
from the online arithmetic in other phase-gradient candidates. The
ordinary key-one state factorizes, up to bit order, as

```math
|\chi_1\rangle
 =\bigotimes_{r=1}^k
 \frac{|0\rangle+e^{-2\pi i/2^r}|1\rangle}{\sqrt2}.
```

Replace factors with $`r>R`$ by plus. Without any phase alignment,
the state error is at most $`\sqrt2\pi2^{-R}`$ by telescoping.
Taking $`R=\min\{k,\lceil\log_2(4\pi/\eta)\rceil\}`$ and
synthesizing the remaining finitely many single-qubit states to total
error at most eta/2 gives a preparation F with

```math
\|F|0^k\rangle-e^{i\varphi}|\chi_1\rangle\|\le\eta,
\qquad T_F,d_F=O_\eta(1),\quad G_F=k+O_\eta(1).   \qquad(15)
```

Only the k prepared wires are clean; no other helper is required. The
phase phi is fixed and input-independent. It is not discarded during
a controlled call: preparation is followed by the actual inverse,
which cancels it. Common inverse costs obey the same bounds.

If an entire exact online word L satisfies
$`L(|\chi_1\rangle\otimes I)=|\chi_1\rangle\otimes B`$,
then, for J_0 inserting zero,

```math
\|F^\dagger LFJ_0-J_0B\|\le2\eta.                \qquad(16)
```

Insert the ideal eigenstate twice and use unitarity to prove (16).
Intermediate return is unnecessary when the whole identity is known.
For words consisting of controlled modular additions and operations on
other registers, every ideal addition already factors out the same
eigenstate, so the whole identity holds with B equal to the intended
phase word. All private-return requirements belong in B. A full
operator implementation error of L adds once to (16); bounds on its
individual gates must still be composed. Exactly dyadic shifts modulo
Q have no angle-rounding error.

## 8. Scope of the recent compilation components

Vasconcelos, [arXiv:2609.34659v1](https://arxiv.org/html/2609.34659v1),
Theorem 8, gives a coherent approximation to minus a single-qubit
rotation, with O(log(1/epsilon)) clean wires and doubly logarithmic
bounded-width depth. Its error includes leakage. Theorem 11's parallel
compiler allocates a separate catalyst per simultaneous phase gate;
it does not supply one logarithmic shared register for a full bank.

Wu et al., [arXiv:2609.36574v1](https://arxiv.org/html/2609.36574v1),
charge the weighted-sum encoder in parallel phase kickback. Their cheap
encoder requires disjoint binary supports of the phase weights. The
counter's fine weights repeat p+1 times, so that stated grouping rule
needs at least p+1 sequential groups. Their numerical adder count also
allows measurement and feedforward. These are limits of the supplied
schedule, not of every unitary shared-catalyst mechanism.

Neither paper pays the missing full-bank operation in Section 5. The
general upper and lower frontier remains unchanged.

## 9. One-pass row filters: clean dimension, error and dirty mixing

A smaller coherent row filter has two distinct width questions. Let
$`A(z)=U_n^{z_n}\cdots U_1^{z_1}`$ read n mismatch literals once,
and let J insert k clean wires, with $`r=2^k`$. If the remaining dirty
space has dimension D, J has rank D in an rD-dimensional space. Allow
arbitrary mixing of these clean and dirty wires. Suppose

```math
A(0)=I,\qquad \|J^\dagger A(z)J\|\le\rho
\quad(z\ne0).
```

The reflected row $`E=A^\dagger(I-2JJ^\dagger)A`$ obeys the exact
complete error identity

```math
E(0)J=-J,\qquad
\|(E(z)-I)J\|=2\|J^\dagger A(z)J\|\le2\rho.       \qquad(17)
```

This includes every dirty-return error and arbitrary references.
Put $`P_j=U_j\cdots U_1`$ and $`\Phi_j=P_j^\dagger J`$, including
$`P_0=I`$. For $`j\lt\ell`$, the cross block
$`\Phi_\ell^\dagger\Phi_j`$ is the compression of the nonempty
interval input $`U_\ell\cdots U_{j+1}`$. Its norm is at most rho.
For $`F=\sum_{j=0}^n\Phi_j\Phi_j^\dagger`$, therefore,

```math
\mathrm{tr}F=(n+1)D,\quad
\mathrm{tr}F^2\le(n+1)D(1+n\rho^2),\quad
\mathrm{rank}F\le rD.
```

The trace/rank inequality proves

```math
r\ge\frac{n+1}{1+n\rho^2}.                       \qquad(18)
```

In particular, exact orthogonal row filtering requires $`r\ge n+1`$,
even with arbitrary dirty mixing. If a serial full-bank schedule budgets
$`\rho=O(\tau/2^p)`$ with $`n=p+1`$, (18) retains logarithmic clean
count. This is a consequence of that sufficient row-error budget; it
does not prove such a budget necessary for every coherent whole-bank word.

A stronger fixed-error bound applies when each dirty key is invariant
and the initialized register in one key fiber has dimension r. For a
unit vector psi, the inverse-prefix vectors $`P_j^\dagger\psi`$ have
pairwise overlap at most rho. Packing their separated unit vectors in
$`\mathbb C^r`$ gives $`r=\Omega_\rho(\log(n+1))`$ for fixed
$`\rho\lt1`$. More generally, a one-pass fixed preparation/readout
with literal compressed amplitudes $`\pm1/\sqrt2`$ and error $`\epsilon\le1/8`$
has separation $`\Delta=\sqrt2-2\epsilon`$ and
$`n+1\le(1+2/\Delta)^{2r}`$. No commutativity is needed: multiplying
two inverse-prefix vectors by $`P_\ell`$ identifies their difference
with the corresponding interval input minus the zero input. General
zero transitions can be absorbed into a fixed interaction frame.

Thus changing only the row transitions from commuting to noncommuting
does not enter the [paid promotion window](PARALLEL_CLEAN_PROMOTION.md#51-paying-any-sufficiently-small-packing-permutation)
while the fixed-fiber premise remains. This is a one-pass filter limit,
not an unrestricted native-depth bound. The next construction shows why
its fixed-fiber premise cannot be dropped.

## 10. A complete dirty Majorana echo pays return but retains linear depth

Take $`M=4\lceil(p+2)/4\rceil`$ and $`b=M/2`$, so M exceeds p+1.
On b arbitrary qubits use the standard Jordan--Wigner Majoranas
$`\gamma_0,\ldots,\gamma_{M-1}`$, with $`\gamma_0=X_0`$.
Choose an actual Clifford D whose adjoint is their unsigned M-cycle.
The even cycle has determinant minus one: its implementation includes
a Pauli correction, not just Majorana plane rotations.

Here is a paid construction. The cycle is the product of the matchings
$`j\mapsto-j`$ and $`j\mapsto1-j`$ modulo M, with b-1 and b pairs.
Each matching uses independent commuting $`S_P^{\pm1}`$ axes
$`P=i\gamma_a\gamma_b`$. Their product gives the cycle with signs.
Conjugation by a Majorana monomial indexed by S gives sign vector
$`v_j=|S|+[j\in S]`$ over F_2. This map is invertible for even M,
so an actual Pauli repairs every required sign, including the odd
orientation. The scalar of the emitted D is retained.

This also gives a constructive counterexample to the fixed-fiber width
claim. In k independent blocks initialize only each first qubit in plus,
leaving all other inputs arbitrary. Let J insert those k bits, and put
$`P_0=JJ^\dagger`$, $`A_w=(D^w)^{\otimes k}`$ and
$`P_w=A_w^\dagger P_0A_w`$. At every $`1\le w\le p+1`$, each
transported first Majorana anticommutes with its original one. Expanding
the projector gives

```math
J^\dagger P_wJ=2^{-k}I,\qquad
\|((I-2P_w)-I)J\|=2^{1-k/2}.                     \qquad(19)
```

Hence fixed row accuracy needs only a fixed clean bank when dirty keys
really mix. Together with the ordinary exact cyclic counter of dimension
less than $`2(n+1)`$, (18)--(19) give the component tradeoff
$`r=\Theta(\min\{n+1,\rho^{-2}\})`$, up to power-of-two rounding,
for $`0\lt\rho\lt1`$. This tradeoff makes no shallow native-depth
promise.

For the complete indicator, even this row error and initialization can
be removed. Use a single entirely arbitrary block per row and define

```math
A_u=D^{w_u},\qquad R_u=A_u^\dagger\gamma_0A_u,\qquad
Q_u=\gamma_0R_u\gamma_0R_u.
```

The Majorana anticommutation relations give the full-space identity

```math
Q_u=\begin{cases} I,&w_u=0,\\-I,&w_u\ne0.\end{cases}
\qquad
\prod_uQ_u=(-1)^{m-Y_x}I=(-1)^{Y_x}I.             \qquad(20)
```

Here m is even. Every private input and reference returns literally;
there is no clean counter, phase-gradient state or approximate row budget.
This multiple-pass Pauli readout does not assume the orthogonal
one-pass inner filter of (18).

The native bill is also explicit. Expose each matching's Pauli axes by
an actual private Clifford and its actual inverse. On each exposed target,
the chronological word
$`T_t,\mathrm{CNOT}_{c,t},T_t^\dagger,\mathrm{CNOT}_{c,t}`$
equals $`\omega^{-c}\mathrm{CS}(c,t)`$,
where $`\omega=e^{i\pi/4}`$. Thus a controlled matching has two T
bands, all targets private. The controlled Pauli repair is Clifford.
The control-only scalar cancels exactly in $`A_u^\dagger\gamma_0A_u`$;
it is not a free phase convention for an isolated controlled call.
Let common forward/actual-inverse bounds be

```math
T_{cD}\le2(M-1),\quad d_{cD}\le4,\quad G_{cD}=O(b^2),
\qquad G_{\rm pow}=O(b^2).
```

The last quantity pays a fixed power of D. Write $`w_u`$ as its constant
$`|u|+1`$ plus p+1 signed positive-control bits. This uses controlled
D or its actual inverse without changing any address bit at a native
boundary. Four A calls and four literal Majorana Paulis implement Q.
All m rows run in parallel on disjoint private blocks, sharing only
read-only address/bank controls. Therefore

```math
\begin{aligned}
T&\le8m(M-1)(p+1),& d&\le16(p+1),\\
b_{\rm dirty}&\le mb,&
G&\le m\{4[(p+1)G_{cD}+G_{\rm pow}]+32\}.
\end{aligned}                                                   \qquad(21)
```

These bounds hold for the actual inverse. The XOR indicator adds 2m
Hadamards. One clean signal and a final H would supply literal flat
compression at amplitude $`1/\sqrt2`$, but merely normalize this
already complete linear-depth query. The known
[bank-motion construction](BANK_MOTION_WALSH_INDICATOR.md) already has
the better bound $`4p`$. Thus (20) resolves return and key dependence,
while retaining the old depth scale.

### The joint-cycle shortcut has a scoped native cost

A separate implementation of $`A(x)=D^{|x|}`$ cannot become shallow
while retaining the actual
[Majorana-preserving Clifford-boundary grammar](WHOLE_WORD_MAJORANA_LIMITS.md#2-actual-clifford-boundaries-give-address-characters).
Each entry of its native adjoint matrix after d T bands then has at most
$`2^d`$ Walsh characters. Since M is divisible by four, the cycle has a
normalized complex eigenvector v with eigenvalue i. Thus
$`v^\dagger O_A(x)v=i^{|x|}`$, whose $`2^p`$ Walsh coefficients all
have magnitude $`2^{-p/2}`$.

If the emitted U approximates A tensored with identity on any helpers
on the entire space to error
$`\epsilon\lt1/2`$, up to one common scalar, use the Majorana norm identity
$`\|\sum_jz_j\gamma_j\|=\|z\|_2`$ only for real vectors z.
The conjugation bound gives $`\|(O_U-O_A)z\|_2\le2\epsilon\|z\|_2`$
for every such z. Since $`O_U-O_A`$ is a real matrix, its complex
spectral norm equals its real spectral norm. Thus the same bound applies
to the complex cycle eigenvector v. The functional
$`v^\dagger O_Uv`$ uses at most $`M^2 2^d`$ Walsh characters.
Parseval on the missing flat coefficients yields

```math
d\ge p-2\log_2M+\log_2(1-4\epsilon^2).           \qquad(22)
```

Arbitrary dirty helpers may be appended after the original b wires
in the fixed Jordan--Wigner ordering. The original M Majoranas are then
unchanged. Extend v by zero on the new Majorana coordinates, so only
the original M by M block enters. An odd Pin action may put parity
signs on later Majoranas; these do not alter that original block.
The actual Clifford boundaries, including endpoints, must
preserve that frame and fix the computational address. Equation (22)
bounds the separately required A, not a resynthesized whole Q, the
completed full bank, or unrestricted Clifford+T circuits.

The revision is therefore concrete. Dirty-key mixing can beat the
fixed-fiber clean requirement and can give exact complete return. The
present native realization still retains linear depth. Improving its
initialization is no longer the next question; a qualifying whole echo
must change its native boundary mechanism or supply a different complete
source with the simultaneous budgets.

<a id="distributed-cyclic-counter"></a>

## 11. Parallel counter visits and the missing common shift

The row sum can be moved from a shared arithmetic target to a shared
read-only control. This gives a sublogarithmic raw counter visit using
only logarithmically many clean inputs. Two completed flags then
suffice for exact comparison if the clocks admit a suitable common
shift. The shift has an exact logarithmic-depth implementation. No
sublogarithmic implementation is supplied here.

### 11.1 A paid flag for an arbitrary dirty key sum

For row $`u`$, put $`e=x\oplus u`$, $`k=\lceil\log_2(p+1)\rceil`$
and $`Q=2^k>p`$. Use a clean $`k`$-bit register C and completely
arbitrary $`k`$-bit clocks $`G_1,\ldots,G_p`$. The complete counter
visit is

```math
S_e|j,z_1,\ldots,z_p\rangle
   =|j,z_1+je_1,\ldots,z_p+je_p\rangle ,
```

with additions modulo Q. Use the
[private Remaud–Vandaele adder](../../docs/PARALLEL_DIRTY_LOOKUP.md#exact-modular-adders).
An m-bit call has $`T,G=O(m\log(m+2))`$ and
$`D_T=O(\log^2(m+2))`$, with no additional work. The first
operand can change within the call and returns at its end.

For one clock allocate k arbitrary k-bit addend registers
$`a_0,\ldots,a_{k-1}`$. Let $`L_{\mathrm{sum}}`$ add their sum into
the clock. Compute that sum into a root addend by a balanced tree of
private adders, add the root into the clock and reverse the tree.
Every addend returns. Disjoint additions within each level give
$`D_T=O(\log^3(k+2))`$ and $`T,G=O(k^2\log(k+2))`$.

Let J add $`2^rj_re_i`$ into each addend $`a_r`$ in parallel. On its
suffix of $`k-r`$ bits, apply the
[read-only signed-increment echo](../../docs/PARALLEL_DIRTY_LOOKUP.md#a-read-only-controlled-increment-using-two-additions)
with two fast private adders. Its literal is $`j_re_i`$. Toggle the
low bit of its own arbitrary helper word by a returned private-phase
Toffoli with controls $`j_r,x_i`$. If the row literal is complemented,
also add $`j_r`$ by CNOT into that helper bit. No public control is a
target and every other gate acts on private work. Thus J returns
all increment helpers, with $`D_T=O(\log^2(k+2))`$ and
$`T,G=O(k^2\log(k+2))`$. Its combined helper width is $`O(k^2)`$.

Execute chronologically $`L_{\mathrm{sum}}^\dagger,J,L_{\mathrm{sum}},J^\dagger`$.
The two initial addend sums cancel, leaving $`\sum_r2^rj_re_i=je_i`$
in the clock and returning the addends. The p clocks have disjoint
private pools and share only read-only j and x, so the complete visit
and its actual reverse both have
$`T,G=O(pk^2\log(k+2))`$, dirty width $`O(pk^2)`$ and
$`D_T=O(\log^3(k+2))=o(\log(p+2))`$.

Write $`\theta=(\theta_1,\ldots,\theta_p)`$ for the dirty Fourier
keys, with the Fourier convention of (3), and put

```math
w_e(\theta)=\sum_i e_i\theta_i\bmod Q,\qquad
P_e(\theta)=[w_e(\theta)=0].
```

Here $`P_e`$ is a projector on the entire clock space.
For $`A_e=H_C^{\otimes k}S_eH_C^{\otimes k}`$, the clean compression
$`J_C^\dagger A_eJ_C`$ equals $`P_e`$. The key-sector compression
is the uniform sum of the Q phases
$`\exp(2\pi i j w_e(\theta)/Q)`$. Thus $`A_e`$ sends $`C=0^k`$
exactly to zero or into its orthogonal complement in each key sector.
Conjugate the returned XOR at $`C=0^k`$ by $`A_e`$ and its actual reverse.
The resulting flag $`B_e`$ toggles any arbitrary target by $`P_e`$
on the entire $`C=0^k`$ code. It returns C and the arithmetic helpers
with the same resource orders as $`S_e`$. All clock superpositions
and references are allowed.

### 11.2 Exact comparison from two flags

Suppose a unitary D acts only on the clocks and possibly additional
arbitrary private work. The required condition (23) is

```math
 P_e D P_e=0\quad\text{for every }e\ne0.
```

One sufficient condition is that every nonzero Fourier transition
of D have displacement $`(r,\ldots,r)`$ for an odd r modulo Q.
Then it changes $`w_e`$ by $`r|e|`$, which is nonzero modulo Q for
$`0\lt|e|\le p`$. No choice or initialization of r is needed.

The two projections $`P_e`$ and $`D^\dagger P_eD`$ are orthogonal for
$`e\ne0`$, since their product equals zero by (23). They are both
identity for $`e=0`$. In particular they commute on every address
sector. Use $`B_e`$ and $`D^\dagger B_eD`$ as the two flag words on
arbitrary private bits a and b. Apply a target Toffoli at all four
corners of these two flag toggles, reversing the two flags in the
last two edges. In their common projector sectors the four target
toggles have parity
```math
 at+(a+\beta)t+(a+\beta)(t+\gamma)+a(t+\gamma)=\beta\gamma
 \quad\text{over }\mathbb F_2 ,
```

where t is the initial value of b and $`\beta,\gamma`$ are the
projector eigenvalues. The flag predicates act on the clocks and
private work only. The full result is therefore the completed XOR
of $`[e=0]`$. The clean counter, clocks, flags and all helpers return
exactly, including for arbitrary reference extensions. The reverse
uses the actual native reversals for each flag and D.

The word uses four completed flags, two of each type, four private
Toffolis and four appearances of D or its actual reverse. It reuses
the same k clean counter bits throughout. With a native read-only
shift satisfying (23), depth $`o(\log(p+2))`$ and polynomial counts
and width for D and its reverse would give a qualifying row source.
Use disjoint row banks and the
[shifted-sector row promotion](PARALLEL_CLEAN_PROMOTION.md#random-shift-clean-promotion),
then polynomial-overhead fusion and the
[approximate query transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md).
An approximation to a D satisfying (23) preserves the clean counter
sector exactly, but the four-call hybrid only bounds its comparison
error by four times the D error per row. That bound alone does not
justify a fixed-error parallel full bank.

### 11.3 An exact shift and the cost of making it joint

The diagonal phase

```math
D(z_1,\ldots,z_p)=(-1)^{\mathrm{MSB}(\sum_i z_i\bmod Q)}
```

depends only on the modular sum and is antiperiodic under increasing
one clock by $`Q/2`$. Its Fourier support consists of common odd
displacements. A paid variant uses the
[masked sum accumulator](../../docs/DIRTY_SUM_COMPRESSION.md#1-the-sum-interface-and-imported-primitives)
with one arbitrary k-bit accumulator a. In the chronological
accumulate–$`Z_{\mathrm{MSB}(a)}`$–actual-reverse word, its phase is

```math
 (-1)^{\mathrm{MSB}(a+\sum_i z_i+g(h)\bmod Q)}.
```

The work-dependent offset $`a+g(h)`$ is fixed in each preserved
computational work sector and leaves the same Fourier support.
Every work bit returns. The exact construction has polynomial p
counts and dirty width, zero clean inputs and $`O(k)=O(\log(p+2))`$
native T-depth. A literal classical encoder of the most significant
sum bit cannot give sublogarithmic Toffoli depth: restricting each
summand to its low bit gives an elementary symmetric polynomial of
degree $`Q/2\le p`$. This restriction is about producing that bit
with classical reversible gates.

The ideal common shift by one factors into p independent prescribed
phase gradients. These have nonnative roots when $`Q>8`$. A shared
global shift could cost only constantly many error contributions,
but sharing the physical clocks among all rows also merges their
counter adds into sums of the row control words:
$`\sum_{u:u_i\ne x_i}j_u`$ into clock i. Commutation of the completed
adds does not give a parallel native schedule for those shared
targets. A separate joint construction is needed in either case.

### 11.4 Exact independent shift words require many clock bits

Suppose $`D`$ is a tensor product of B native blocks on disjoint
nonempty clock groups with arbitrary private work. Test it on the
Fourier input tuple zero. For one block, the support of its Fourier
output tuples $`r`$ is invariant under $`r\mapsto(1+8v)r`$
simultaneously on its coordinates. Indeed the native entries and the
normalized Fourier-zero input lie in $`\mathbb Q(\zeta_8)`$.
The cyclotomic automorphism fixing $`\zeta_8`$ carries each nonzero
output matrix element on private work to the matrix element with
that transformed tuple. For $`Q\lt8`$ one can use the identity
automorphism.

Select one coordinate per block and one supported output tuple per
block. A selected residue zero already violates (23) for the
corresponding one-coordinate mismatch. Otherwise each selected
residue has a 2-adic valuation in $`\{0,\ldots,k-1\}`$. If $`B>7k`$,
eight selected residues share a valuation t. Their odd normalized
residues have nine partial sums modulo eight. Two coincide, giving
a nonempty subset with normalized sum zero modulo eight. When
$`k-t>3`$, multiply the tuple of any one of these blocks by
$`1+8v`$ to cancel the remaining multiple of eight modulo
$`2^{k-t}`$. Its selected odd residue is a unit, so such v exists.
When $`k-t\le3`$, the sum already vanishes modulo $`2^{k-t}`$.

The supported combined Fourier output thus lies in $`P_e`$ for e
supported on that nonempty subset, as does the all-zero input.
The corresponding private-work operator is a tensor product of
nonzero matrices, so $`P_eDP_e\ne0`$. Consequently (23) requires
$`B\le7k`$. This includes nonmonomial local shift words with arbitrary
private work. It does not bound a joint native circuit whose clocks
interact.

At $`p=Q-1`$ condition (23) imposes a stronger joint target. Every
supported output tuple from the Fourier-zero input must have no
nonempty subset summing to zero modulo Q. Its Q partial sums in
any order then exhaust $`\mathbb Z_Q`$. Swapping two adjacent
entries changes only one partial sum, which must remain the same.
Hence all entries equal one r, whose additive order is Q and thus
whose parity is odd. On the Fourier-zero input, D must therefore
have output supported on the common odd-key subspace

```math
 \mathrm{span}\{|\chi_r\rangle^{\otimes p}:r\text{ odd}\}
 \otimes\mathcal H_{\rm private}.
```

For $`Q>8`$, even two nonempty independent exact blocks cannot
meet this condition: after fixing a supported tuple of the first
block, all supported tuples of the second must have the same
fixed odd label. Its Galois orbit contains another label. This
necessity concerns one physical input column of the joint repair,
including all its arbitrary work. It does not imply a general
native T-depth bound.

<a id="odd-parity-row-kernel"></a>

### 11.5 An odd-parity row kernel suffices for the complete indicator

Let $`h:\mathbb F_2^p\to\mathbb F_2`$ be a fixed row predicate, and let
$`H_{x,u}=h(x+u)`$. Over $`\mathbb F_2`$, convolution by h squares to
$`(\sum_x h(x))I`$. Indeed the terms for distinct summands cancel in
pairs, and each identical pair has group sum zero. Thus an odd truth
parity makes H a binary involution. Apply H to the arbitrary bank
with $`O(4^p)`$ CNOTs by binary elimination. Apply the exact parallel
row phases $`(-1)^{h(x+u)(HY)_u}`$, then actually reverse the CNOT
word. Their total phase is $`(-1)^{(H^2Y)_x}=(-1)^{Y_x}`$, with
every bank and work input returned. The row flags can use dirty
comparison targets and an exact controlled-Z echo. A row compiler with
gatewise read-only sharing, no clean inputs, polynomial row costs and
sublogarithmic native depth therefore suffices through
[full-bank fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md). The
[exact sector cover](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
also permits a fixed number of clean inputs per row, with constant
overhead. At growing clean count the cover depth must be paid. Both
approximate row synthesis and approximate promotion need a joint error
bound because several row predicates may be true at the same address.

For a modulus $`Q=2^k\le p`$, every predicate
$`h(e)=[\sum_i a_i e_i\equiv t\pmod Q]`$ has even truth parity. The
coefficients of $`\prod_i(1+z^{a_i})`$ modulo $`z^Q-1`$ give these
truth parities over $`\mathbb F_2`$. Each factor vanishes at one,
and $`z^Q-1=(z+1)^Q`$, so the product vanishes in the quotient.
Small dyadic moduli therefore cannot be used individually in this
convolution recovery. In contrast, at least one residue predicate
of $`|e|\bmod 3`$ has odd truth parity for every p. Evaluate
$`(1+z)^p`$ at either primitive cube root in $`\mathbb F_4`$ to see
that the residue coefficients cannot all vanish.

The exact modulo-three row thus gives another sufficient operation,
with no equality predicate needed in that reduction. The standard
[quantum fanout construction](https://arxiv.org/abs/quant-ph/9903046)
uses initialized fanout copies and the nonnative cyclic eigenbasis.
A read-only, few-clean native modulo-three row with the required
exact return and depth remains unpaid.

<a id="monomial-clock-repair-bank"></a>

## 12. Monomial clock repair and dense-bank error

Let e be a p-bit mismatch, set $`Q=2^k>p`$, and let $`P_e`$ be
the distributed-counter projector onto
$`\sum_i e_i\theta_i=0\bmod Q`$ in the clocks' Fourier keys. Equivalently,
it averages translation $`z\mapsto z+e`$ on the computational
clock labels. For nonzero e, every orbit has length Q. Let the repair
$`D_e`$ be monomial on the clocks and arbitrary private work, with phase
$`\omega^{q_e(z,h)}`$ on input $`z,h`$ for $`\omega=\zeta_8`$. For each
fixed j, assume the three Boolean digits of $`q_e(je,0)`$ have degree at
most s in e. The j dependence is unrestricted. The repair may depend
on the row and the entire prescribed table.

### Exact silence requires high address-phase degree

Set $`\lambda=1-\omega`$ and

```math
F(e)=\sum_{j=0}^{Q-1}\omega^{q_e(je,0)}.
```

At e zero its value is $`Q\omega^{q_0(0,0)}`$, of
$`\lambda`$-valuation $`4k`$. Two is totally ramified in $`\mathbb Z[\zeta_8]`$, since
$`\lambda`$ has norm two and $`2`$ is associated to $`\lambda^4`$.

For a phase digit of weight $`2^b`$, expand its binary algebraic
normal form. Every nonconstant monomial M contributes the factor
$`1+(\omega^{2^b}-1)M(e)`$, with any constant phase retained. A
coefficient of a degree-a monomial in a product of these factors
is divisible by $`\lambda^{\lceil a/s\rceil}`$, since at least
that many nonconstant factors must cover those a variables.
Use the Boolean relations $`e_i^2=e_i`$. For s zero, all nonconstant coefficients vanish.

Thus F modulo $`\lambda^{4k+1}`$ is a nonzero multilinear polynomial
of degree at most $`4ks`$. A nonzero multilinear polynomial of
degree d over any coefficient ring is nonzero on at least
$`2^{p-d}`$ Boolean inputs. Choose a nonzero coefficient of highest
degree and take the finite difference in its variables. The
same nonzero coefficient is obtained on every assignment of
the complement. Each complementary fiber therefore contains a
nonzero evaluation. In particular at least
$`2^{p-4ks}-1`$ false mismatch values have nonzero F.

On a false mismatch the normalized uniform vector on
$`\{je:0\le j\lt Q\}`$, tensored with the private zero column, belongs
to $`P_e`$. The sum of the output-orbit entries of $`P_eD_eP_e`$
on that column is $`F(e)/Q`$, even if the permutation of $`D_e`$
moves all the private work. Thus exact false silence implies
$`p\le4ks`$.

### Nonzero false compression is quantized

Use normalized uniform clock-orbit vectors and computational
private columns as an orthonormal basis of $`P_e`$. Each entry
of $`P_e D_e P_e`$ is $`1/Q`$ times a sum S of at most Q eighth
roots. A nonzero S satisfies $`|S|\ge1/Q`$. Its absolute field
norm is $`|S|^2|\sigma(S)|^2`$, a nonzero integer, and the other
complex embedding obeys $`|\sigma(S)|\le Q`$ by the same sum.
Consequently

```math
\nu_e=\|P_eD_eP_e\|\ne0\quad\Longrightarrow\quad
\nu_e\ge Q^{-2}.
```

### The literal bank forces logarithmically many nonlinear rounds

For a false row put $`P=P_e`$ and $`Q'=D_e^\dagger P_eD_e`$.
Use the two completed flag XORs $`B_P,B_{Q'}`$ on private flags
a,b, respectively. With the output in its active Fourier value,
each Toffoli corner is the phase $`C(a,b)=(-1)^{ab}`$. In the
chronological word

```math
 C,\ B_P,\ C,\ B_{Q'},\ C,\ B_P^\dagger,\ C,\ B_{Q'}^\dagger,
```

the block returning the two flags from zero to zero is

```math
 E=\sum_{\alpha,\beta\in\{0,1\}}
      (-1)^{\alpha\beta} Q'_\beta P_\alpha Q'_\beta P_\alpha ,
```

with $`P_1=P,P_0=I-P`$ and similarly for Q'. Only the third
corner contributes a sign on paths returning both flags to zero.
For $`\psi\in\mathrm{ran}(P)`$ the formula gives
$`E\psi=(I-Q'-P Q')\psi`$. Choose a unit vector there maximizing
$`\langle\psi,Q'\psi\rangle=\|Q'P\|^2=\nu_e^2`$. The actual
row output has real overlap $`1-2\nu_e^2`$ with the input, whose
ideal false-row action is identity. This calculation permits
noncommuting P and Q'.

In the literal bank of disjoint row clocks and flags, activate
one false row by its output Fourier value. Keep all other rows
inactive and put both arbitrary flags of each inactive row in
plus states. Its corners and flags then act as identity on this
chosen input. A uniform literal-phase full bank error $`\tau\lt1`$ therefore forces
$`\nu_e\le\tau/2`$ in each false row. All the corresponding
overlap factors are positive.

Activate the M bad false rows at one address and use their product
of maximizing inputs, with other rows inactive as above. The overlap
with the ideal output is at most $`(1-2Q^{-4})^M`$. The error promise forces
this overlap to be at least $`1-\tau^2/2`$, hence
$`M\le C_\tau Q^4`$. A common translated template gives at least
$`2^{p-4ks}-1`$ bad rows at each address. Even for row-dependent
templates with the same degree bound, averaging their bad
addresses gives some address with at least this many bad rows.
Consequently $`4ks\ge p-4k-O_\tau(1)`$.

For a hardwired table, an unwanted active false row can instead
be disabled on a private unit vector in $`\ker P\cap\ker Q'`$.
Each rank fraction is $`1/Q`$, so this intersection is nonempty
for $`Q>2`$. No flag moves and the four identical phases cancel.
The matching row is exact for every D because both projections
are identity. Hence the same product estimate applies to bad
false rows in any fixed table's support.

For any fixed enabled set of density c, average the bad
mismatch counts over all addresses for only those rows. Some
address then has at least $`c(2^{p-4ks}-1)`$ bad enabled rows.
The argument also applies when the clean counter is reused sequentially.
Every flag returns its initialized counter, so the fixed-address row
actions are still a product on the disjoint clocks and flags.

For a classical reversible implementation with t Toffoli rounds
and bounded-degree diagonal phase exponents, the fixed-j input
restriction $`z_i=je_i`$ is Boolean-linear in e. Each wire's
degree grows by at most two per Toffoli round. Computing the
three bits of the accumulated exponent modulo eight multiplies
degree by at most a constant. Thus $`s=O(2^t)`$ and the literal
bank needs $`t\ge\log_2p-\log_2k-O_\tau(1)`$. In particular this
is logarithmic at $`k=\lceil\log_2(p+1)\rceil`$. This is a bound
on that monomial reversible-round implementation. A general
native word may resynthesize its monomial action through
nonmonomial intermediate operations. Shared row clocks, a
different extraction word, and nonmonomial repairs also fall
outside the bank obstruction.

<a id="distance-weighted-clock-repair"></a>

## 13. A distance-weighted approximate repair suffices

Use the distributed clocks of Section 11 with dyadic order
$`2p\le Q=2^k\lt4(p+1)`$. Give each clock i an arbitrary local
work space and a native companion $`D_i`$ acting only on that
clock and its work. Its Fourier grade $`D_{i,r}`$ is the complete
operator block that changes the clock Fourier key by r modulo Q,
with no restriction on transitions of its private work. Write
$`d_Q`$ for shortest distance on the cyclic group. Consider the
weighted operator condition

```math
 D_i=\sum_r D_{i,r},\qquad
 \sum_r \|D_{i,r}\|s^{d_Q(r,1)}\le B .
```

B is a fixed constant. No input of the companion is initialized.
It suffices to have this condition at $`s=C\sqrt p`$ for a
sufficiently large C depending on B and the desired fixed indicator
tolerance. If each companion and its actual reverse have polynomial
counts and width and $`o(\log(p+2))`$ native T-depth, this gives
the normalized complete indicator and improves fixed-accuracy frame
depth through the
[complete approximate transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md).
The weighted norm is a premise, with no native construction here.

### 13.1 Grade compression decays with Hamming distance

Put $`D=\bigotimes_i D_i`$ on the distributed clocks inside one
row. For a weight-w mismatch e, only grade choices on the active
coordinates S satisfying $`\sum_{i\in S}r_i=0`$ can connect the
two $`P_e`$ spaces. Since $`w\le Q/2`$, the cyclic triangle
inequality gives

```math
 \sum_{i\in S}d_Q(r_i,1)\ge
 d_Q\left(\sum_{i\in S}r_i,w\right)=w .
```

The inactive companions remain unitaries of norm one. On the
active factors expand the weighted grade sums. The entire false
compression therefore satisfies
$`\nu_e=\|P_eDP_e\|\le (B/s)^w`$. This bound uses all private
work and coherent Fourier inputs.

### 13.2 Three visits cancel the linear false error

Let $`P=P_e`$, $`Q'=D^\dagger PD`$ and let $`U(P,Q')`$ be the
actual four-corner XOR using the two completed flags of Section 11.2.
The clean counter returns exactly at each flag boundary and the
word is the exact matching XOR at e zero. At a false row
$`\|PQ'\|=\nu_e`$. For $`\nu_e\le1/2`$, let $`Q_0`$ be the
orthogonal projector onto the range of $`(I-P)Q'`$.
The projector difference $`E=Q'-Q_0`$ has norm at most $`\nu_e`$,
and its two diagonal blocks relative to P have norm at most
$`\nu_e^2`$. Indeed the map $`(I-P)|_{\mathrm{ran}\,Q'}`$ has
squared singular values in $`[1-\nu_e^2,1]`$ and
$`\|P Q'P\|\le\nu_e^2`$.

Each $`Q'`$-flag is affine in its projector. Expand the
fixed-length four-corner word around the proof-only $`Q_0`$.
The constant term is I because $`PQ_0=0`$. The linear term
from the off-diagonal blocks of E is off-diagonal relative to P
because every other factor commutes with P. The remaining linear
term and all higher terms have norm $`O(\nu_e^2)`$.

On one additional arbitrary qubit use the native order-three
Clifford $`G=i(SH)^2`$ with eigenphases $`\pm2\pi/3`$.
Its literal scalar i matters when controlled. A controlled
$`i`$ is S on the control, controlled S has a fixed native
phase word, and controlled H is a Clifford-conjugated
controlled Z between native target rotations. These gates pay
controlled G at constant cost and depth. Compute the P flag
into a new arbitrary bit, apply G controlled by that bit,
reverse the flag and subtract G controlled by the original bit.
This gives a relative conditional phase $`\exp(\pm 2\pi iP/3)`$
in each fixed sector of the bit and the G eigenkey. Call its
emitted word V. The four-corner word ignores both new wires.

Run the three conjugates $`V^j U(P,Q') V^{-j}`$ for
$`j=0,1,2`$ using actual inverses. Their linear off-diagonal
terms cancel since the three third roots sum to zero. All
other error has norm $`O(\nu_e^2)`$ at each false row. At the
match the exact XOR commutes with V and its cube is the same
XOR. The same fixed number of flags, companion visits and their
actual reverses pays the new word. The only initialized bits
are the returned k-bit counter. On superpositions of addresses,
dirty work and reference inputs all estimates are operator bounds.

### 13.3 Pay the bank and the borrowed-code promotion

At each address the total row error is bounded by

```math
 C_1\sum_{w=1}^{p}\binom pw(B/s)^{2w}
   = C_1[(1+B^2/s^2)^p-1].
```

For any fixed target error this can be made sufficiently small
by choosing C large. The
[one-selected-copy promotion](PARALLEL_CLEAN_PROMOTION.md#summable-row-error-promotion)
preserves this sum with a constant factor and returns all dirty
work within the same full-operator budget. Within a row the p
companions have disjoint supports. Across rows the original
address is used only at the read-only counter boundaries, so
both counter and companion T layers align. The promotion uses
$`k=O(\log(p+2))`$ and adds $`O(\log\log(p+4))`$ depth
and polynomial address overhead. The
[indicator fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
absorbs this overhead at fixed error, and brings the complete
source to the original count and width orders.

The weighted condition concerns the actual full operator of a
native companion. Independent binary angle synthesis with a
uniform operator-error estimate gives no comparable suppression
of the distant grades. An exactly false-silent joint companion
from Section 11 is a separate sufficient route.

### 13.4 A linear code permits a constant weight

Fix a binary linear map $`L:\mathbb F_2^p\to\mathbb F_2^m`$ with
$`m=O(p)`$ and $`|L(e)|\ge\vartheta m`$ for every nonzero input e,
for a fixed $`0\lt\vartheta\lt1/2`$. Such a family follows by choosing
independent random linear forms. For each nonzero mismatch its image
has the binomial weight law, and increasing the fixed ratio m/p
makes the union bound over $`2^p-1`$ nonzero inputs strictly less
than one. Exhaustive verification and selection give a terminating
uniform compiler. Repeating the code a fixed number of times
increases m without changing the relative distance.

Use the coordinate forms of $`L(x\oplus u)`$ as the mismatch
controls and choose a power of two $`2m\le Q\lt4m`$. The
comparison is still gatewise read-only on the address. To add the
product of one selector digit with a coordinate form into a local
target, XOR that form into an arbitrary private bit by CNOTs and X,
apply a Toffoli from that bit and the digit, reverse the parity
calculation, and repeat the Toffoli. The two occurrences cancel the
unknown private bit. Distinct private helpers and the shared-control
cubic schedule keep the substituted Toffoli batches at constant
depth, with $`O(mp\log(m+2))`$ extra Cliffords per row and
polynomial T-count and returned dirty width. The counter retains
$`O(\log(p+2))`$ initialized bits and
$`O((\log\log(p+4))^3)`$ T-depth.

The grade-compression and quadratic-echo arguments now apply with
encoded weight $`|L(x\oplus u)|`$. For any fixed $`s\gt B`$, their
total false-row error is at most

```math
C_1(2^p-1)(B/s)^{2\vartheta m}.
```

Choose the constant repetition factor to make this decay
exponentially in p. Bounded smaller arities use the exact dirty tree.
The same one-selected-copy promotion and normalization therefore
suffice at **any fixed $`s\gt B`$**. This removes a growing weight
parameter from the companion contract. The weighted spectral
condition still requires suppression at distances of order Q.

<a id="local-companion-depth-count"></a>

## 14. A depth–count constraint for a local weighted companion

Let $`Q=2^k\ge16`$. For a single companion D with t T gates
and d physical T bands, assume the weighted condition of Section 13
with $`s\gt B\ge1`$ on its entire clock and arbitrary private work. Then

```math
\sqrt{dt\,k\ln4}+4k
\ \ge\ \frac{Q}{4}\log_2 s-\log_2 B .
```

Set $`B_r=\langle\chi_r|D|\chi_0\rangle`$ on the private work.
The replacement of D by its grade-one component has norm error at
most $`B/s\lt1`$. Since D is unitary, the minimum singular value
of $`B_1`$ is at least $`1-B/s`$. Let M be the Pauli transfer
matrix of the map $`H\mapsto B_1 H B_1^\dagger`$. If the work
has w wires, M is invertible of dimension $`r=4^w`$.

Work in $`K=\mathbb Q(\zeta_Q)`$ with the valuation at its
unique dyadic prime normalized by $`v(\sqrt2)=1`$.
M factors into the input injection with $`|\chi_0\rangle`$,
the full native Pauli transfer matrix, and the output compression
with $`|\chi_1\rangle`$. The Pauli coordinates of the injection
are normalized coefficients of the stabilizer projector and have
denominator at most Q. Extraction is by clock Pauli expectations
in $`|\chi_1\rangle`$, sums of Q roots-of-unity terms divided by Q.
Thus every coefficient of either boundary has valuation at least
$`-2k`$. The
[R103 band-minor bound](PRECISION_COUNT_DEPTH_TRADEOFF.md#2-valuations-of-minors-of-a-physical-t-band)
applies over this local ring, since its one-band elimination still
uses powers of $`\sqrt2`$ and units. The full Pauli space has
dimension $`r4^k`$. Cauchy–Binet with the two rectangular
boundaries gives

```math
v(\det M)\ge-r\bigl(\sqrt{dt\,k\ln4}+4k\bigr).
```

Every field embedding $`\sigma_u:\zeta_Q\mapsto\zeta_Q^u`$,
for u odd, commutes with the adjoint operation. In the conjugated
Pauli basis, each Y acquires at most a sign. The transfer matrix
$`\sigma_u(M)`$ represents the compression of $`\sigma_u(D)`$
between the Fourier keys zero and u in that basis. The conjugate
of D is unitary, so each such determinant has absolute value at
most one. For $`u\equiv1\pmod8`$ the embedding fixes the native
matrix of D. Its compression then has norm at most
$`\|D_u\|\le B s^{-d_Q(u,1)}`$, so
$`|\det\sigma_u(M)|\le B^{2r}s^{-2r d_Q(u,1)}`$.
There are Q/8 such embeddings, with summed cyclic distance
$`Q^2/32`$. The product of the determinants therefore has magnitude
at most $`B^{rQ/4}s^{-rQ^2/16}`$.

The norm of $`\det M`$ is a nonzero rational with only dyadic
denominators. Total ramification and degree Q/2 give its dyadic
valuation as $`(Q/4)v(\det M)`$. Its magnitude is at least two
to that power. Comparing both bounds proves the inequality.

In the constant-weight coded regime $`Q=\Theta(p)`$ with fixed
$`s\gt B`$, it requires $`dt=\Omega(p^2/\log(p+2))`$.
Sublogarithmic d therefore requires
$`t=\omega(p^2/\log^2(p+2))`$ per companion. At
$`s=\Theta(\sqrt p)`$ it instead requires
$`dt=\Omega(p^2\log(p+2))`$ and sublogarithmic depth needs
$`t=\omega(p^2)`$. The estimate permits larger polynomial
and singly exponential budgets. It does not constrain a joint
repair unless that repair supplies the local weighted companion.

**Monomial companions.** At fixed $`s\gt B`$, let X denote the
cyclic translation of the clock. The weighted assumption gives

```math
\|XDX^\dagger-D\|
\le \frac{2\pi B}{Q}
     \sup_{\ell\ge0}(1+\ell)s^{-\ell}.
```

Any computational-basis monomial native unitary has eighth roots
as its nonzero entries. Indeed each field embedding of the native
ring extends to a unitary of the same monomial support, so all
conjugates of each nonzero entry have modulus one. Its field norm
is one. The only possible denominator lies at the unique dyadic
prime, where the norm gives valuation zero. The entry is thus an
algebraic integer. All its powers have bounded conjugates and hence
bounded integer coefficients in a fixed integral basis, so two
powers coincide. The finite roots of unity in $`\mathbb Q(\zeta_8)`$
form a cyclic group containing the eighth roots, and the field
degree four excludes a larger order divisible by eight. Thus two
distinct such monomials have distance at least
$`2\sin(\pi/8)`$: a differing permutation gives distance at
least $`\sqrt2`$ on some column, and a differing phase gives the
eighth-root spacing. Conjugation by the cyclic permutation X
preserves this class of monomial matrices. For large Q the bound
forces equality $`XDX^\dagger=D`$, so D has only grade zero,
with weighted norm s. This contradicts $`s\gt B`$ at any native
depth. Nonmonomial companions and joint repairs remain possible.

<a id="private-gauge-mod3"></a>

## 15. Private phase gauges pay a complete approximate MOD3 row

Sublogarithmic MOD3 rows would make the
[odd-parity row kernel](#115-an-odd-parity-row-kernel-suffices-for-the-complete-indicator)
a complete bank construction if their errors were summable across rows.
The following zero-clean row uses no ternary-valid dirty promise.
Its worst-case accuracy bill identifies a smaller sufficient precision
problem.

Let b be a stable conjunction of a fixed number of read-only literals.
On an arbitrary r-bit register, let X increment modulo $`q=2^r`$ and
$`D|j\rangle=e^{2\pi i m j/q}|j\rangle`$ for a compiled integer m.
The exact lifted word satisfies

```math
  (X^b)^\dagger D X^b D^\dagger
       =e^{2\pi i m b/q} I .
```

This uses integer wraparound, and returns the whole gauge for every
input. Choose q and m so that $`2\pi m/q`$ approximates the desired
angle to error epsilon, with $`r=O(1+\log(1/\epsilon))`$. The
[read-only increment](../../docs/PARALLEL_DIRTY_LOOKUP.md#a-logarithmic-depth-increment-using-two-dirty-bits)
uses $`O(\log(r+2))`$ native depth with its paid private controls.
Replace its literal toggle by a returned fixed-arity predicate toggle
on the same private target. To complement the gauge conditionally,
fan out an arbitrary private bit by CNOT, toggle that bit by b,
repeat fanout and actually reverse the toggle. The old bit cancels
and returns. All non-Clifford supports remain private.

The operator D is a product of r single-qubit Z rotations up to
one common phase, which cancels in the actual inverse. Approximate
each rotation to $`O(\epsilon/r)`$ using the
[zero-clean rotation word](PROGRAMMED_ROTATION_DEPTH.md#4-cover-both-failure-subspaces-without-initialized-work).
Allocate disjoint dirty pools, so all rotations align in native depth
$`O(\log(r+2))`$. The actual commutator approximates the literal
scalar phase with full error $`O(\epsilon)`$, and is identity on
b zero with its actual inverse. The T, Clifford and dirty costs
are $`2^{O(r+1)}`$. The whole word preserves the public controls
and returns dirty work to the stated tolerance.

Let G cycle the two-qubit labels $`0\to1\to2\to0`$ and fix 3.
Write $`H_{12}`$ for the usual Hadamard on labels 1 and 2, and
$`R_{01}(\theta)`$ for the real rotation sending label 0 to
$`\cos\theta\,|0\rangle+\sin\theta\,|1\rangle`$, each fixing the
other labels. Set $`\theta=\arccos(1/\sqrt3)`$ and $`\omega=e^{2\pi i/3}`$. Then

```math
\begin{aligned}
 F&=H_{12}R_{01}(\theta)\,\mathrm{diag}(1,-1,i,1)H_{12},\\
 F^\dagger G F&=\mathrm{diag}(1,\omega^{-1},\omega,1).
\end{aligned}
```

Indeed, the first three columns on the cycle are
$`(1,1,1)/\sqrt3`$, $`(1,\omega,\omega^2)/\sqrt3`$ and its conjugate.
The fourth column is the fixed label. An affine relabeling reduces
$`H_{12}`$ to a [native controlled Hadamard](../../docs/CONDITIONAL_GEOMETRIC_SOURCE.md).
The two label phases use CZ and controlled S with complemented controls
as needed. Controlled S is exact since $`2xy=x+y-(x\mathbin\oplus y)`$.
For the rotation let a indicate high bit zero. A Clifford on the low
bit exposes its Y eigenvalue $`(-1)^b`$. The required phase exponent in
units of theta is $`a(2b-1)=b-(a\mathbin\oplus b)`$. Apply
$`U_\theta=\mathrm{diag}(1,e^{i\theta})`$ on b, XOR a into b, apply
the word's actual adjoint and undo the XOR and basis change. This
costs two visits to a private word. A common phase cancels, and full
word error epsilon contributes at most $`2\epsilon`$ even with private
work or intermediate label motion.

For each mismatch literal $`e_i`$ and each of the two
nontrivial target eigenlabels, use a private gauge phase of
angle $`-2\pi/3`$ on label 1 or $`+2\pi/3`$ on label 2. Its predicate is the
conjunction of $`e_i`$ and the two computational label bits.
In that phase stage all target eigenlabel and address bits
are read-only. The $`2p`$ gauge words share no non-Clifford
targets. Request individual errors $`O(\delta/p)`$ here and
$`O(\delta)`$ for the private basis word and its actual inverse.
The basis changes and the aligned phase stage consequently
approximate $`G^{\sum_i e_i}`$ on the entire target and dirty
space at error delta. For
$`L_\delta=\Theta(1+\log(p/\delta))`$, they use
$`p\,2^{O(L_\delta)}`$ gates and dirty work, with
$`O(\log(L_\delta+2))`$ native depth.

Initialize the two target bits to the label zero, run the power,
toggle an arbitrary output on the desired residue label and
actually reverse the power. The ideal comparison returns the two
bits and all work. The full isometry error is at most twice the
power error. The [bounded-clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)
makes a zero-clean full XOR row with constant factors.
Use the promotion with all original address bits read-only.
An extra Fourier-output control yields the row phase, with
the same resource orders. Thus inverse-polynomial full row
error already has doubly logarithmic native depth.

### A finite-angle precision target suffices

The literal odd-kernel convolution still sums errors over
$`2^p`$ independent rows. Uniform row error
$`\delta=O(2^{-p})`$ suffices for fixed bank error.
Then $`L_\delta=\Theta(p)`$, recovering logarithmic depth
by the gauge word above.

For a complete word of independent row phases this error
sum is intrinsic up to constants. Fix the computational
address and suppose each row u acts on disjoint private
work, is exactly inactive at zero Fourier control, and
has active block $`U_u`$ with target scalar $`\sigma_u`$.
Every subset of Fourier controls is available. If the
bank error is $`\tau\lt1`$, then every subset product of
$`\sigma_u U_u`$ is within tau of identity. Each row's
relative eigenangles lie in
$`[-\alpha,\alpha]`$ with $`\alpha=2\arcsin(\tau/2)`$.
Choose the largest positive angle in each row having one
and omit the others. Every partial sum lies in that
interval modulo $`2\pi`$. The increments are at most
alpha, so these partial sums cannot jump to a different
translate of the interval. Hence the positive maxima
sum to at most alpha, and the negative extrema obey the
same bound. Thus $`\sum_u\|\sigma_u U_u-I\|\le2\alpha`$.
This necessity uses exactly inactive row factors and
separate work. It allows nonuniform precision and leaves
joint native readout outside its premise.

There is a weaker necessary condition with shared work and
arbitrary noncommutation between active rows. Keep exact
control sectors and fixed address, and write
$`V_u=\sigma_u U_u`$ on the common work space.
Average the ordered active words over all subsets. The result
is the ordered product A of $`(I+V_u)/2`$. Each subset word
is unitary and within tau of I, so its Hermitian part is at
least $`(1-\tau^2/2)I`$. Thus the minimum singular value
of A is at least $`1-\tau^2/2`$. If the work dimension is D,
multiplicativity of the determinant gives

```math
 \sum_u\|V_u-I\|_{2,\mathrm{tr}}^2
     \le -8\log(1-\tau^2/2),
 \qquad \|M\|_{2,\mathrm{tr}}^2=D^{-1}\mathrm{Tr}(M^\dagger M).
```

Indeed, an eigenvalue $`e^{i\theta}`$ contributes
$`|\cos(\theta/2)|`$ to the determinant modulus of its
row average, and minus its logarithm is at least
$`|e^{i\theta}-1|^2/8`$. For translated copies the mean
squared row error in normalized trace is exponentially
small even in this shared model. Large operator errors
on mutually disjoint small subspaces can meet this
condition. A comparison only on a leaking initialized
code, or a different joint gate word, need not satisfy
these product estimates.

A sufficient replacement uses just two one-qubit phase words.
For the third-root word Q, require
$`Q\approx e^{i\phi}\mathrm{diag}(1,\omega)\otimes I`$ on the public
bit c and its arbitrary private work. A sufficient gatewise convention
uses c only as a CNOT control. In target label $`2u+v`$, the eigenphase
for mismatch bit e is $`\omega^{e(u-v)}`$. The bit identity

```math
 e(u-v)\equiv(e\mathbin\oplus u)-(e\mathbin\oplus v)-u+v\pmod 3
```

implements it with Q, its actual adjoint, its actual adjoint and Q
on the four indicated parities, respectively. To substitute an
affine parity for c, replace each CNOT from c by CNOTs from its
constituents to the same private target and an X there for a
constant term. No parity bit is initialized, all T targets stay
private, and the four common phases cancel. The private basis word
$`U_\theta`$ at $`\theta=\arccos(1/\sqrt3)`$ may carry a common phase
and may move its private qubit. Each word needs zero clean inputs,
full operator error at most $`2^{-L}`$ up to its common phase,
$`2^{O(L)}`$ counts and dirty width, and $`o(\log(L+2))`$ native depth.
Using this pair in the same row
construction at $`L=\Theta(p)`$ gives a full indicator
with a common singly exponential cost exponent.
[Shared-layer fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#shared-layer-indicator-fusion)
would normalize it and improve general fixed-accuracy
frame depth. A joint bank readout that avoids exponentially
small independent row error is an alternative to this
fixed-angle precision improvement.
