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

The revision in Sections 9--10 separates this counter from a different
dirty-key mechanism. It proves a clean-dimension/error tradeoff for
one-pass row filters and pays an exact full-bank echo, whose emitted
native depth remains linear. Neither result changes the general frontier.

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
the chronological word $`T_t,\mathrm{CNOT}_{c,t},T_t^\dagger,
\mathrm{CNOT}_{c,t}`$ equals $`\omega^{-c}\mathrm{CS}(c,t)`$,
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
