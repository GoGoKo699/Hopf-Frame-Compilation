# Coherent counters: an exact odd-key code and the missing native schedule

[Complete source target](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Clean-source promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) ·
[Paid affine phase batches](PARALLEL_AFFINE_DETUNING.md)

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

There is a second independent scheduling obligation. Giving all m rows
their own initialized counters uses m k clean counter wires. Reusing
one counter gives (12). A shallow row primitive alone supplies neither
a legal parallel full-bank schedule nor a native parallel clean-removal
word. The shared odd-key bank removes gradient preparation and its
precision cost, but does not initialize these m counters.

Nor may all W calls simply be merged before one readout. Since

```math
\sum_u w_u=mp/2+m-|Y|,                            \qquad(14)
```

their product on one counter/accumulator is independent of x. The
interleaved reflections are essential to the proven source identity.

The concrete next question is a returned joint modular accumulation or
a whole-bank counter echo of sublogarithmic native depth, on arbitrary
banks and with only O(log p) clean wires. Such a word must pay its own
simultaneous count, depth, precision and return budgets. The present
identity does not exclude one, and does not construct one.

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
