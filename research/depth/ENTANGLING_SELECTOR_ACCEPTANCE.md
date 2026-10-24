# Entangling selector controls and complete acceptance screens

[Working checkpoint](../../WORKSPACE.md) ·
[Co-moving stages](CO_MOVING_LOOKUP_STAGE.md) ·
[Earlier local-selector obstruction](CO_MOVING_SELECTOR_READOUT.md) ·
[Complete recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Fixed physical-row entangling selector operations have cheap, exact native
controls. This removes a synthesis uncertainty, but not the acceptance
problem. A common two-dimensional selector module excludes arbitrary words
using linear selector permutations and global Hadamards. Partial Fourier
mixing leaves that module, yet the specified four-path source still has a
two-bit Gram floor of one sixth. Adding an affine selector flip gives an
exact finite algebraic separator and an all-size separation at address zero;
the latter has a direct central-word proof, but that word reads the fixed
physical bank bit on the full space. Neither supplies an all-address lookup
construction. A square-stage commutator gives a complete two-bit
calibration and identifies precisely why its affine-Pauli extension stops.
Non-Clifford selector controls are also paid below. A uniform-selector
witness rejects the common-Fourier dressing family at every size. A separate
rank argument excludes any one-clean source with at most two completed
co-moving stage calls and common signal/selector interlayers, without
assuming a particular path polynomial. A crossing-rank bound also rejects
every three-stage cell-factorized CCZ conjugation family, including mixed
Fourier frames. A nonlocal wrapped candidate escapes that rank bound but
still has an exact three-bit common-vector obstruction.

None of these statements changes R96 or proves a general depth lower bound.

## 1. Registers and the unpaid acceptance identity

There are p address bits, $`m=2^p`$ arbitrary bank bits Y and m arbitrary
selector bits s. Retain the complete co-moving stage

```math
R_C=i^{n_C}L_C^h,\qquad
h=\sum_uY_us_u\pmod2,\qquad
n_C=\sum_uY_us_u\sum_{c\in S_C}(x_c\oplus u_c).
```

C is a disjoint CNOT matching, $`S_C`$ its control coordinates, and L
relabels address, bank and selectors together. Its actual inverse and
returned implementation helpers are paid in the linked stage proof.
All states below may be coherent and entangled with references. Define

```math
\Omega|x,Y,s\rangle=(-1)^{Y_x}|x,Y,s\rangle,\qquad
\Pi_\pm=(I\pm\Omega)/2,\qquad H_Y=\Pi_+/2.
```

The desired bounded-clean accepted block K must satisfy
$`\|K^\dagger K-H_Y\|_{\rm op}\le\delta`$ on the complete space.
A useful first requirement is
$`K\Pi_-=0`$; positive-sector coverage and normalized strength remain
separate obligations. All impossibility statements specify their actual
operation grammar. Small-size certificates do not exclude a family patched
at those sizes.

## 2. A controlled entangling Clifford with four native T layers

For each physical row-label coordinate j, let

```math
E_j=\prod_{u:u_j=0}\mathrm{CX}(s_u;s_{u\oplus e_j}),\qquad
M_p=\prod_{j=0}^{p-1}E_j.
```

On binary selector columns these are tensor factors of
$`\begin{pmatrix}1&0\\1&1\end{pmatrix}`$. The factors commute,
$`M_p^2=I`$, and $`(M_ps)_v=\bigoplus_{u\subseteq v}s_u`$.
The direct word has pm/2 CNOTs, zero T gates and no helpers.

Put q=m/2 and split the rows by the highest label bit. With
$`A=M_{p-1}`$ (A=I for p=1), use binary matrices

```math
M_p=\begin{pmatrix}A&0\\A&A\end{pmatrix},\quad
J=\begin{pmatrix}I&0\\I&I\end{pmatrix},\quad
Q=\begin{pmatrix}I&I\\0&I\end{pmatrix},\quad
P=Q\begin{pmatrix}I&0\\0&A\end{pmatrix}Q
 =\begin{pmatrix}I&A+I\\0&A\end{pmatrix}.
```

Direct multiplication gives $`P^2=I`$ and $`M_p=PJP`$.
The explicit Q,A,Q chronology for P uses
$`G_P=m+(p-1)m/4`$ CNOTs. No generic change-of-basis cost is hidden.

To control M by an arbitrary signal a, emit P, the batch
$`\prod_{i=0}^{q-1}\mathrm{CCX}(a,s_i;s_{q+i})`$, and the actual
P-adjoint. The [shared-control bilinear batch](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#2-pay-for-a-root-controlled-complete-bank-transition)
has common forward/actual-inverse bounds
$`T\le6q+(q\bmod2)`$, $`G\le12q+3`$ and depth four, without helpers.
The common control receives no T gate for even q; common S powers are
reduced modulo four in both orientations. Thus

```math
T_{\mathrm cM}\le3m+(q\bmod2),\qquad
G_{\mathrm cM}\le(p+15)m/2+3,\qquad D_T(\mathrm cM)\le4.       \qquad(1)
```

There is zero clean work, helper work or error. The control returns;
selectors undergo exactly the requested controlled permutation. Literal
phases are +1. The actual native reverse has the same majorants.

For $`H_s=H^{\otimes m}`$, Fourier summation gives
$`H_sU_AH_s=U_{A^{-\mathsf T}}`$. Therefore
$`N_p=H_sM_pH_s=U_{M_p^{\mathsf T}}`$ is another CNOT permutation,
not a nonmonomial mixer. Transposing the explicit CNOT construction gives
the same control ledger (1), without adding Hadamard walls.

In contrast, let $`H_j=\prod_{u:u_j=0}H_{s_u}`$ and
$`C_j=H_jM_pH_j`$. Fixed coordinate relabeling reduces it to the same
block calculation. Its controlled implementation is $`H_j,\mathrm cM,H_j`$,
so

```math
T_{\mathrm cC_j}\le3m+(q\bmod2),\qquad
G_{\mathrm cC_j}\le(p+17)m/2+3,\qquad D_T(\mathrm cC_j)\le4. \qquad(2)
```

For the highest-coordinate split, write input selectors as (a,b) and
outputs as (v,w), with q bits in each half. The complete kernel is

```math
\langle v,w|C_j|a,b\rangle
 =2^{-q}(-1)^{a\cdot Aw+a\cdot b+v\cdot w+v\cdot Ab}.       \qquad(3)
```

Indeed the first H wall introduces t, M sends (t,b) to (At,A(t+b)),
and fixing w forces $`t=Aw+b`$. The two Fourier phases give (3).
This is a dense real Clifford involution. In particular,

```math
|g_j\rangle:=C_j|0^m\rangle
 =2^{-m/2}\sum_s(-1)^{q_j(s)}|s\rangle,\qquad
q_j(s)=\sum_{u:u_j=0}s_us_{u\oplus e_j}\pmod2.              \qquad(4)
```

It is not a supplied state: (4) describes an allowed input's image under
the paid unitary. Arbitrary selector inputs remain in the contract.

More generally, a binary linear involution B has $`(B+I)^2=0`$.
Put $`r=\mathrm{rank}(B+I)`$, choose paired vectors with $`(B+I)v_i=u_i`$, and
complete the u vectors to a basis of the kernel. This gives r disjoint
CNOT blocks and fixed coordinates. A literal CNOT basis change of cost
$`L_P`$ therefore controls B with $`T\le6r+(r\bmod2)`$,
$`G\le2L_P+12r+3`$ and depth four. Gaussian elimination gives the
safe bound $`L_P\le d(d-1)+3(d-1)`$ on d wires. This is not a
theorem about arbitrary controlled Cliffords or free phase corrections.

## 3. An all-size common module, including global Hadamards

At the marked address/bank input $`x=0,Y=e_0`$, every linear matching
fixes x and Y, its phase exponent vanishes, and

```math
R_C|0,e_0,s\rangle=|0,e_0,C^{s_0}s\rangle.                 \qquad(5)
```

This is a bijective phase-free selector permutation: C fixes row zero,
so $`s_0`$ is unchanged. Consequently it fixes both $`|0^m\rangle`$
and $`|+^m\rangle`$. At the unmarked input $`x=0,Y=0`$, every
stage is identity on the entire selector space.

Let $`\mathcal S=\mathrm{span}\{|0^m\rangle,|+^m\rangle\}`$.
Every selector CNOT network fixes this space pointwise, and global H
preserves it by swapping its two displayed, nonorthogonal vectors.
Thus every word in linear co-moving stages, fixed selector CNOT networks,
global H, their actual adjoints and signal-controlled versions acts
identically on this module at the two bank inputs. Arbitrary signal-only
interlayers, initialized signals and common acceptance rules retain that
equivalence. Tensoring with independent returned implementation helpers
and references does not change it.

The target Grams on these copies of the module are zero and I/2.
If their common accepted Gram is G, the triangle inequality gives

```math
\max\{\|G\|,\|G-I/2\|\}\ge1/4.                           \qquad(6)
```

This is independent of word length, number of matchings and signal-angle
precision. It includes the zeta and transposed-zeta proposal and its
global-H modification. It excludes neither partial H walls, affine
selector flips, selector phases, address-dependent interlayers nor
unrestricted native circuits. The proof uses the empty and singleton
banks; it does not assume that bank-weight-dependent signal operations
are present.

## 4. Breaking that module still does not give four-path acceptance

The graph state (4) supplies an analytic direction check. Under (5), its
amplitude changes by

```math
(-1)^{s_0[q_j(Cs)+q_j(s)]}.                                 \qquad(7)
```

If $`Ce_j=e_j`$, the matching of rows separated by $`e_j`$ is
preserved, so the state is fixed. Choosing a CNOT's target direction
therefore does not remove this particular witness. Choosing its control
direction changes the state; it is the relevant next test, not a source
proof. The same graph state is fixed by global H, by the two-vertex
graph-state identity on each disjoint pair.

For p=2 fix $`C=(0\to1)`$, $`D=(1\to0)`$ and the moved-direction
choice

```math
A_0=R_C,\quad A_1=C_0R_CC_0,\qquad
B_0=R_D\ \hbox{or}\ H_sR_DH_s,\quad B_1=C_1B_0C_1.
```

A one-signal source with three signal interlayers and the two diagonal
branch multiplexors has accepted block

```math
K=\sum_{a,b=0}^1c_{ba}B_bA_a,\qquad
c_{ba}=\langle0|S_2|b\rangle\langle b|S_1|a\rangle
        \langle a|S_0|0\rangle.                             \qquad(8)
```

Each multiplexor is implemented by controlled selector conjugation around
one uncontrolled stage; no controlled R is assumed. Exact fixed signal
interlayers give depth at most $`2(116)+4(4)=248`$ plus their own
T-depth. Approximate signal synthesis must be charged separately.

The following stronger test allows *arbitrary complex* coefficients c,
not merely those realizable in (8). Restrict all four branches to the
marked input $`x=0,Y=e_0`$, retaining the whole 16-dimensional selector
space. Order them as $`(B_0A_0,B_0A_1,B_1A_0,B_1A_1)`$ and put
$`G_{ij}=\mathrm{tr}(W_i^\dagger W_j)`$. Direct exact multiplication
gives, respectively,

```math
G=\begin{pmatrix}
16&9&9&7\\9&16&6&9\\9&6&16&9\\7&9&9&16
\end{pmatrix},\qquad
G=\begin{pmatrix}
16&9&10&13/2\\9&16&6&10\\10&6&16&9\\13/2&10&9&16
\end{pmatrix}.                                               \qquad(9)
```

In both cases $`G-8\mathbf1\mathbf1^{\mathsf T}`$ is Hermitian,
has diagonal eight and off-diagonal absolute row sum at most five.
It is positive definite. Hence the average marked acceptance is at least
$`|\sum c_{ba}|^2/2`$. On the empty bank every branch is identity,
so its acceptance is $`|\sum c_{ba}|^2`$. Complete Gram error epsilon
would imply

```math
\epsilon\ge\tfrac12|\sum c_{ba}|^2
 \ge\tfrac12(\tfrac12-\epsilon),\qquad \epsilon\ge1/6.       \qquad(10)
```

This is a complete two-bit counterexample, including coherent selectors,
for both displayed four-path sources and all their scalar coefficients.
It is not an all-p lower bound or an exclusion of other accepted words.
Removing the previous common invariant was necessary but insufficient.

## 5. Separation after an affine selector flip

### Complete two-bit certificate

An independent exact check considers the enlarged p=2 alphabet

```math
\{R_C,R_D,M_2,N_2,H_s,F\},\qquad F=X^{\otimes4},\qquad
S=R_C+R_D+M_2+N_2+H_s+F.                                   \qquad(11)
```

F takes $`|0^m\rangle`$ outside the module used in (6). More strongly, the complete positive and
negative charge representations of this alphabet have no nonzero common
intertwiner at p=2. This finite statement retains all 1024 address, bank
and selector basis states, not a singleton promise.

The reproducible certificate partitions (x,Y) into its simultaneous
C,D relabeling orbits, keeping all 16 selector columns within each block.
There are 20 such blocks and ten in each charge. Every generator preserves
each block. Compute the characteristic polynomial of S in each block over
$`\mathbb F_{65537}`$, using $`i\mapsto256`$ and the exact rational
Hadamard entries. All 100 opposite-charge polynomial pairs have greatest
common divisor one. The native stage matrices are independently compared
with the existing literal complete-stage fixture.

Why this proves the finite claim: these characteristic polynomials are
monic over $`\mathbb Z[i,1/2]`$. A coprime reduction modulo this prime
makes their resultant nonzero, so the original opposite-charge spectra
are disjoint over the complex numbers. The Sylvester equation
$`S_-T=TS_+`$ therefore has only T=0, even when S is not normal. An
intertwiner of every generator would also intertwine their sum, and
decomposes into these block pairs, so it must vanish.

The test uses exact modular arithmetic, not numerical eigenvalue
separation. Its characteristic-polynomial calculation was independently
audited by Faddeev--LeVerrier and Cayley--Hamilton checks. It does not
supply a normalized polynomial, a physical spectral filter, a stable
gap, or any all-size bound. In particular S is non-Hermitian; disjoint
finite spectra do not certify a uniform singular-value gap or a low-depth
implementation of its charge projectors.

### All-size separation within address zero

The following invariant-image proof retains its original three-image
hypothesis. The [direct central word](#fixed-row-central-word) below
strengthens the conclusion to a single nonempty matching and identifies
exactly what the separator reads away from address zero.

There is also an analytic result for every p at least two, including every
bank background. Restrict both charge representations to the invariant
address-zero space, so their charges are determined by $`Y_0`$. Use the
alphabet $`\{R_C:C\in\mathcal C\}\cup\{M_p,N_p,H_s,F\}`$, now with
$`F=X^{\otimes m}`$. Suppose the group generated by the linear matchings
has at least three distinct images of the top row $`t=(1,\ldots,1)`$.
The two opposite CNOTs already satisfy this condition for every p at least
two. Under these hypotheses the two address-zero charge representations
have no nonzero intertwiner.

First the binary zeta identities give

```math
M_pFM_p=X_{s_0},\qquad N_pFN_p=X_{s_t},\qquad
H_sX_{s_t}H_s=Z_{s_t}.                                    \qquad(11a)
```

Indeed M sends the all-ones selector column to the row-zero unit vector,
while its transpose sends it to the top-row unit vector. At x=0 the
row-zero contribution to n vanishes, every matching fixes row zero, and
flipping $`s_0`$ changes h by exactly $`Y_0`$. It leaves n unchanged and
commutes with the joint relabeling. Thus, with the actual stage adjoint,

```math
E_C:=R_CX_{s_0}R_C^\dagger X_{s_0}
       =L_C^{Y_0}\qquad\hbox{on the whole }x=0\hbox{ space}. \qquad(11b)
```

The phases cancel literally, including on nonempty backgrounds. Here L
still relabels the bank and selectors together; it is not a selector-only
substitution. Equation (11b) is identity on the positive sector and L on
the negative sector.

Let T intertwine every generator from the positive address-zero sector
to the negative address-zero sector. Since the generators are unitary,
T also intertwines their adjoints and all words. Its image W is therefore
fixed pointwise by every joint $`L_C`$, and hence by every generated
$`L_g`$. It is invariant under the selector Paulis $`X_{s_t},Z_{s_t}`$
from (11a). Choose g,h with $`t,g(t),h(t)`$ all distinct. For w in W,
both w and $`X_{s_t}w`$ are fixed by $`L_g`$, which gives
$`X_{s_{g(t)}}w=X_{s_t}w`$. The analogous Z relation uses h. Hence

```math
X_{s_t}X_{s_{g(t)}}w=w,\qquad
Z_{s_t}Z_{s_{h(t)}}w=w.                                    \qquad(11c)
```

These two Pauli operators anticommute: their row supports meet only at t.
They have no common +1 vector, so W=0 and T=0. Taking adjoints excludes
the reverse intertwiner as well. This proves separation for arbitrary
banks and coherent selectors within x=0, not merely for the earlier
empty/singleton witness.

<a id="fixed-row-central-word"></a>

### A direct central word reads the fixed physical row

Use the group commutator convention $`[A,B]=ABA^\dagger B^\dagger`$.
Let C be any nonempty matching and choose a moved row t. Then
$`t,Ct\ne0`$ and $`Ct\ne t`$. On the **complete** space, define

```math
E=[R_C,X_{s_0}],\qquad Q=[E,X_{s_t}],\qquad W=[Q,Z_{s_t}].
```

The exact identities are

```math
\begin{aligned}
E&=i^{Y_0d_C(x,0)(2s_0-1)}L_C^{Y_0},\\
Q&=(X_{s_{Ct}}X_{s_t})^{Y_0},\\
W&=(-1)^{Y_0}I=Z_{Y_0}.                              \qquad(11d)
\end{aligned}
```

Here the exponent of i is a diagonal operator. Flipping selector zero
changes h by $`Y_0`$ and changes n by
$`Y_0d_C(x,0)(1-2s_0)`$. Row zero, its bank bit and every control
coordinate are fixed by C. These facts give the first line with the
actual stage adjoint, including its nonzero-address phase. That phase
commutes with the moved-row X and Z operators. Conjugation by L takes
$`X_{s_t}`$ to $`X_{s_{Ct}}`$, giving the second line. Its two-row
Pauli has exactly one overlap with $`Z_{s_t}`$, giving the last line.
Every basis label returns with precisely the displayed phase, including
arbitrary coherent inputs entangled with a reference.

The expanded word uses eight actual stage/adjoint calls, eight
$`X_{s_0}`$, four $`X_{s_t}`$ and two $`Z_{s_t}`$ factors. In the
positive-S native alphabet, $`X=HS^2H`$ and $`Z=S^2`$ give common
forward/actual-reverse bounds

```math
T_W\le8T_R,\qquad G_W\le8G_R+52,\qquad D_T(W)\le8d_R\le928.
```

The stage's private dirty pool is reused and returns exactly; there is
no clean input. These are bounds for the displayed word, not optimal
costs: its complete action is just the two-Clifford gate $`Z_{Y_0}`$.

On x=0, W equals Omega. Taking t to be the all-ones row, every nonempty
CNOT matching moves it, and (11a) supplies both of its selector Paulis
from M, N, H and F. Thus the original alphabet with **one** nonempty
matching already has no opposite-charge intertwiner within x=0: an
intertwiner would have to intertwine +I and -I. No third row image or
spectral-gap argument is required. This strengthens the earlier
separation theorem, while explaining its limited constructive content.
For general x the same word reads $`Y_0`$, not $`Y_x`$.

<a id="address-transport-cost"></a>

### Extending the fixed-row readout requires addressed transport

Let $`(T_xY)_u=Y_{u\oplus x}`$ and define the coherent one-bank
translation $`V_Y|x,Y\rangle=|x,T_xY\rangle`$. It leaves x unchanged,
is an involution, and has the complete identity

```math
V_Y^\dagger Z_{Y_0}V_Y=\Omega.                         \qquad(11e)
```

Consequently a shallow implementation of this transport would itself
solve the phase lookup. Two global bank-Hadamard walls then give the
complete dirty indicator. An emitted zero-clean unitary with full return
error at most epsilon from $`V_Y\otimes I`$ gives error at most
$`2\epsilon`$ in (11e), by the two-factor norm estimate; its actual
native inverse is required. This is a sufficient reduction to lookup,
not an equivalence claiming that lookup cheaply implements transport.

The straightforward exact transport has p sequential batches. Batch j
swaps each pair of bank wires $`(u,u\oplus e_j)`$ with $`u_j=0`$,
controlled by $`x_j`$. Each batch has m/2 disjoint Fredkins sharing
one control. Two outer CNOTs per pair plus the
[four-layer bilinear batch](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#2-pay-for-a-root-controlled-complete-bank-transition)
give, for $`p\ge2`$, common forward/actual-inverse bounds

```math
T_{V_Y}\le3mp,\qquad G_{V_Y}\le(7m+3)p,\qquad
D_T(V_Y)\le4p.                                        \qquad(11f)
```

There are no helpers, clean inputs or errors. Distinct coordinate
batches share physical bank wires, so commutation of their logical
permutations does not make their native T layers disjoint. Equation
(11f) is an upper bound for this schedule, not a transport lower bound;
its linear depth does not improve the retained logarithmic route.

The same issue occurs if zeta is moved with the address. On selectors,
write $`M_x=T_xM_pT_x`$. Since global F commutes with every row
translation, (11a) gives

```math
M_xFM_x=T_xX_{s_0}T_x=X_{s_x}.                         \qquad(11g)
```

The right side is exactly the requested complete indicator on the
selector bank. Address-relative zeta conjugation cannot be assigned a
free cost to extend the address-zero argument. Nor can a free
$`X_{s_x}`$ be inserted into an address-relative stage echo. The first
unpaid operation in these repairs is addressed transport or selection,
which the fixed-row central word does not provide.

### Scope of the remaining all-address question

The address-zero result does not exclude intertwiners involving nonzero
address orbits, even when one endpoint has address zero. It supplies no
all-address gap, normalized accepted block or complete lookup construction.
The direct word (11d) settles its fixed-fiber readout but does not alter
this boundary; (11e)-(11g) locate the missing address dependence.
If only the same two address coordinates are used at p greater than two,
the entire alphabet commutes with every unused address X. With one bank
marking x but not its spectator-bit flip, this again gives a one-quarter
Gram floor. A scalable proposal must expose every address coordinate,
for example through the already paid full matching family, and prove its
acceptance anew. The finite separator is a seed, not that construction.

## 6. Complete two-bit readout and the boundary of square-stage words

There is a useful positive calibration. For a single-CNOT matching with
control i, let $`Q=\prod_u\mathrm{CZ}(Y_u,s_u)`$ and
$`A_i=QR_C^2`$. The invariance of n and h gives

```math
R_C^2=(-1)^{n_C},\qquad
A_i=Z_s^{a_i},\qquad (a_i)_u=Y_u(1\oplus x_i\oplus u_i).
```

Put $`B_j=H_sA_jH_s=X_s^{a_j}`$. The complete commutator is

```math
A_iB_jA_iB_j
 =(-1)^{\sum_uY_u(1\oplus x_i\oplus u_i)
                         (1\oplus x_j\oplus u_j)}I_s.       \qquad(12)
```

For p=2 and distinct i,j, (12) is exactly Omega, returning every selector.
For larger p it reads a codimension-two subcube parity, not a single row.

One need not emit eight expensive co-moving stages. Directly implement
$`A_i`$ by the shared-control batch $`\prod_u\mathrm{CCZ}(x_i,Y_u,s_u)`$
and CZ on the m/2 rows with $`u_i=0`$. The former has T-count 6m,
four T layers and Clifford count at most 10m+3; the row mask adds 3m/2
Cliffords. Hadamard conjugation adds 2m for B. Thus the word (12) and its
actual reverse have

```math
T\le24m,\qquad G\le50m+12,\qquad D_T\le16,                 \qquad(13)
```

with no initialized or private helpers and m arbitrary returned selectors.
For p=2, two bank-Hadamard walls give the complete indicator with
$`G\le52m+12`$ and the other bounds unchanged. This is a small-size
calibration, not an improvement to the growing-address frontier.

There is a precise limit to extending it with fixed selector Cliffords.
Fix a computational bank Y. All operations in this grammar preserve the
address and bank fibers. Consider any word of canonical work Paulis
$`X^{a(x)}Z^{b(x)}`$ with binary-affine labels, and fixed work Cliffords
acting on selectors and additional work, independent of x; no
address-dependent scalar phases are supplied.
Every such word has normal form

```math
U_Y(x)=\lambda i^{q_Y(x)}X^{a_Y(x)}Z^{b_Y(x)}C,\qquad
q_Y(x)=q_0+\sum_iq_ix_i+2\sum_{i\lt j}q_{ij}x_ix_j\pmod4.  \qquad(14)
```

Pauli multiplication adds twice a binary bilinear form. H adds
$`2a_jb_j`$ and exchanges labels; S adds $`a_j`$ and updates b;
CNOT acts linearly on the labels. The integer lift of a binary affine
form modulo four has only linear and even quadratic terms. These facts
prove (14) by induction, with literal constant Clifford phases in lambda.

If every work input returns with a scalar sign, C is a Pauli up to scalar,
the residual labels are constant, and q has constant parity. The returned
Boolean sign consequently has degree at most two. A singleton bank's
lookup sign has address degree p. Arbitrarily long words in this grammar,
including fixed entangling zeta interlayers, therefore cannot give exact
lookup at p at least three. Unsquared R stages, address motion, non-Clifford
selector operations and accepted signal blocks are outside this statement.

The same grammar has complete lookup error at least $`1/\sqrt2`$ for
p at least three. Otherwise, for every Hermitian work Pauli P,
$`\|UPU^\dagger-P\|\lt\sqrt2`$. Clifford Pauli discreteness forces
equality for every P, hence U is scalar. Two parities of q would give
scalar ratios $`\pm i`$, whose best simultaneous distance to the desired
signs is $`\sqrt{2-\sqrt2}>1/\sqrt2`$. Constant parity reduces to the
quadratic sign just proved; a wrong sign cannot be hidden by one common
scalar phase. This is a grammar-specific error floor, not a depth bound.

Literal phases matter. A freely assigned Hermitian normalization
$`P(a,b)=i^{ab}X^aZ^b`$ is not a free member of the above alphabet.
Indeed $`X^aZ^bP(a,b)=(-i)^{ab}`$, and its finite difference in
$`b=x_2\oplus x_3`$ gives $`(-1)^{x_1x_2x_3}`$. Assigning that phase
without synthesis would invalidate the quadratic closure.

## 7. Consequence for the next constructive question

The native control of the entangling selector involutions is paid. The
complete accepted identity remains unpaid. The CNOT/global-H grammar
cannot be repaired by longer signal processing; the displayed partial-H
four-path source is also rejected. The affine-flip enlarged algebra escapes
the common-intertwiner obstruction at two bits and, for every p at least two, between
the two address-zero sectors. The direct address-zero extraction is settled,
but its full-space action is only $`Z_{Y_0}`$. Other address orbits and a
uniform all-address source or gap remain open. The square-stage calibration cannot be extrapolated through
fixed affine-Pauli words.

Sections 8--12 pay nonlinear selector controls and screen the resulting
three-stage sources. Mixed Fourier frames alone do not repair the
cell-factorized CCZ family. Nonlocal mixing can escape that rank budget,
but the specified transpose-zeta repair still fails at three bits.

The next admissible attack must specify a genuinely different accepted
operator with complete marked cancellation and whole positive-sector
coverage, using all address coordinates. A uniformly gapped intermediate
Gram would be useful only after its normalization and complete native
recovery are paid. Neither finite algebraic distinguishability nor removal
of a known witness supplies this. Preserve arbitrary banks, selectors,
references, actual inverses and off-diagonal terms throughout.

At every fixed requested tolerance, a qualifying sublogarithmic complete
source with fixed polynomial overhead would pass through the retained
fusion and frame-transfer interfaces. Even a constant-depth source leaves
their logarithmic-logarithmic contribution; it would not by itself prove
optimal constant-depth frames.

Exact fixtures: [native selector controls](../../tests/test_entangling_selector_readout.py),
[marked acceptance](../../tests/test_entangling_marked_acceptance.py),
[finite sector separator](../../tests/test_co_moving_sector_separator.py), and
[square-stage calibration](../../tests/test_square_stage_readout.py). The
[fixed-row and transport fixtures](../../tests/test_fixed_row_transport_readout.py)
retain the full-address phase in (11d), test its disagreement with lookup,
and corroborate the translation and moved-zeta identities.
Finite fixtures corroborate the stated identities and finite certificates;
they do not replace the all-size arguments or establish an asymptotic source.

<a id="nonclifford-selector-control"></a>

## 8. A paid non-Clifford selector control and its eight-path source

Let P be a product of CCZ gates on q disjoint selector triples
$`(a_j,b_j,c_j)`$. For an arbitrary signal z, borrow q arbitrary bits
$`r_j`$ and emit, in chronological order, C,B,C,B, where

```math
C=\prod_j\mathrm{CCZ}(r_j,b_j,c_j),\qquad
B=\prod_j\mathrm{CCX}(z,a_j;r_j).
```

The two C phases add modulo two to

```math
\sum_j[r_j+(r_j\oplus za_j)]b_jc_j=z\sum_ja_jb_jc_j.
```

The second B returns every borrowed bit.
Thus this is the complete controlled-P operator, with literal phase,
arbitrary signal/selector/work inputs and arbitrary references.

The [shared-root bilinear batch](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#2-pay-for-a-root-controlled-complete-bank-transition)
pays B with $`T\le6q+(q\bmod2)`$, $`G\le12q+3`$ and four
T layers. Each C uses q disjoint copies of the four-layer singleton
CCZ schedule: $`T=7q`$, $`G=10q`$. Consequently

```math
T_{\mathrm cP}\le26q+2(q\bmod2),\qquad
G_{\mathrm cP}\le44q+6,\qquad D_T(\mathrm cP)\le16.       \qquad(15)
```

There are q dirty helpers, zero clean helpers and zero error. Actual
reversal has the same majorants after reducing adjacent positive S
powers; inverse Cliffords are not assigned a free cost.
For any subset A of the selector wires, put $`J=H_APH_A`$.
Its controlled word is $`H_A,\mathrm cP,H_A`$, adding $`2|A|`$
Clifford gates and no T layer. This includes genuinely non-Clifford
entangling selector operations, not only the earlier Clifford mixers.

For a concrete family at $`p\ge2`$, partition row labels into four-row
cells along label coordinates 0 and 1. In every cell, $`P_{j+1}`$ applies CCZ
to the three corners other than corner j, for $`j=0,1,2`$; here
$`q=m/4`$. Choose three matchings whose control coordinates jointly
cover the address. With one clean signal, define the multiplexors

```math
\mathcal B_j=|0\rangle\langle0|\otimes R_{C_j}
 +|1\rangle\langle1|\otimes J_jR_{C_j}J_j,
\qquad
U=S_3\mathcal B_3S_2\mathcal B_2S_1\mathcal B_1S_0,       \qquad(16)
```

where each S_j is an actual signal-only unitary. A multiplexor uses
controlled J_j, an **uncontrolled** R_Cj, and controlled J_j; no
controlled whole-stage synthesis is assumed. The accepted block is

```math
K=\sum_{b_1,b_2,b_3\in\{0,1\}}c_{b_3b_2b_1}
 A_{3,b_3}A_{2,b_2}A_{1,b_1},\qquad
A_{j,b}=J_j^bR_{C_j}J_j^b,                              \qquad(17)
```

with coefficients fixed by the four signal matrix entries along that
path:

```math
c_{b_3b_2b_1}=(S_3)_{0b_3}(S_2)_{b_3b_2}(S_1)_{b_2b_1}(S_0)_{b_1\,0}.
```

Arbitrary eight coefficients are only a screening relaxation,
not a supplied native linear combination. The three stages and six
selector controls have depth at most $`3(116)+6(16)=444`$ before
the charged signal synthesis. Before those signal costs, their simultaneous
bounds are $`T,b_{\rm dirty}=O(pm)`$ and $`G=O(m^2+pm)`$,
as allowed by polynomial-overhead fusion; the returned private pool can
be reused.
All signal gates, including any requested approximation, retain their
actual costs and full-word error. This pays the candidate circuit family,
not its desired acceptance identity.

<a id="uniform-selector-ray"></a>

## 9. A common Fourier frame still has a complete acceptance floor

Take the origin-singleton bank $`Y=e_0`$ at address x=0. On this
marked fiber, $`n_C=0`$, $`h=s_0`$, and the completed stage is

```math
R_C|s\rangle=|C^{s_0}s\rangle.                         \qquad(18)
```

C fixes physical row zero. The right side is therefore a permutation
within each s_0 sector, so every stage and its actual adjoint fix
$`|+\rangle^{\otimes m}`$. On the empty bank every R_C is identity.
If $`J_j=H_sP_jH_s`$, then every J_j also fixes this uniform state:
$`P_j|0^m\rangle=|0^m\rangle`$. Hence all eight branches in (17)
fix the same marked vector, while each empty-bank branch is identity.
For $`s=\sum_b c_b`$, both inputs have acceptance $`a=|s|^2`$.
Their required Gram values are respectively zero and one half, giving

```math
\|K^\dagger K-\Pi_+/2\|\ge
\max\{a,|a-1/2|\}\ge1/4.                              \qquad(19)
```

This holds for every number of such stages and every coefficient choice,
including every actual signal circuit. More generally it holds whenever
all common selector operations preserve the uniform ray with the same
induced signal action on these two fibers. Arbitrary X-basis-diagonal
mixers $`H_sDH_s`$, phase-free selector permutations and global selector
X flips obey this premise. Arbitrary signal interlayers, more signals,
full-coordinate matching coverage and more precise angles do not repair
this grammar. A computational-diagonal phase changing the uniform ray,
or a standalone Hadamard changing its Fourier frame, need not obey it.

If a physical source's complete initialized isometry is within epsilon
of one in this grammar, its accepted Gram differs by at most
$`2\epsilon`$, by contraction and unitary norm bounds. The floor is
then at least $`1/4-2\epsilon`$. This is a source-grammar limit, not
an obstruction to general non-Clifford selector mixing or native depth.

<a id="coherent-stage-rank-screen"></a>

## 10. A rank budget for arbitrary interleaved complete stage calls

The same origin-singleton fiber gives a stronger short-source screen
without a common selector eigenvector or a path-polynomial assumption.
Let $`D=2^m`$ be the selector Hilbert dimension. If a matching C contains
q disjoint CNOTs on the address, its row permutation fixes $`m/2^q`$
rows and swaps

```math
t_C=(m-m/2^q)/2
```

pairs of nonzero rows. On the $`s_0=1`$ selector sector, the trace
of that permutation is $`2^{m-1-t_C}`$: each swapped pair must have
equal bits to contribute a fixed computational basis column. The stage
is identity on $`s_0=0`$ and an involution on the other sector. Thus

```math
r_C:=\mathrm{rank}(I-R_C)
=\frac D4(1-2^{-t_C})\lt D/4.                          \qquad(20)
```

The identity matching has r_C=0. Actual adjoints have the same restriction
and rank. These statements concern completed stages; their internal
native address motion and returned implementation helpers remain paid.

Now use k initialized signal qubits, an arbitrary private dirty register
of Hilbert dimension W, and h completed stage calls. At their boundaries,
assume every interlayer preserves the marked-origin and empty-origin
bank/address fibers and acts identically on the remaining registers after
identifying those fibers. It may mix signals, selectors and private dirty
work arbitrarily, including with non-Clifford operations. There is no
read-only restriction within a completed stage. On the empty fiber each
query is identity; on the marked fiber its defect rank is
$`2^kWr_{C_j}`$. Telescoping the two complete source words gives

```math
\mathrm{rank}(U_{\rm mark}-U_{\rm empty})
\le2^kW\sum_{j=1}^h r_{C_j}.                           \qquad(21)
```

Multiplication by the common intervening unitaries preserves rank, and
rank is subadditive, which proves the estimate. The initialized code has
dimension DW. Therefore, if

```math
2^k\sum_jr_{C_j}\lt D,                                \qquad(22)
```

there is a normalized code vector on which the two complete source
outputs agree. Any common acceptance map gives identical accepted norms
on this vector. The target Grams are zero and I/2, so their complete
Gram error is again at least one quarter. No scalar empty-bank block,
branch conjugation, accepted unitary fiber or exact clean return has
been assumed. Ordinary signal-controlled stage calls have no larger
defect rank. A multiplexor of different stages must instead use its
actual direct-sum defect rank, or a justified upper bound for it.

For one clean signal, (20)--(22) exclude **every source with at most two
stage calls**, at every $`p\ge2`$, even with arbitrary common entangling
interlayers. At p=2 each nonempty matching has $`D=16,r_C=2`$, so
three calls are also excluded: $`2\cdot3\cdot2=12\lt16`$.
For p=3, a single-CNOT matching has $`D=256,r_C=48`$; three calls
give 288, and this rank screen is inconclusive. Equality in (22) with
the strict sign removed likewise gives no conclusion.

For the particular branch architecture (17), an alternative proof takes
the intersection of the fixed spaces of all branch-conjugated stages.
Its dimension is at least $`D-\sum_{j,b}r_{C_j}`$. This also excludes
arbitrary complex path coefficients under that positive-dimension
condition, without requiring their realization by signal gates.
At p=2 with three binary stages, this common space has dimension at
least four. The all-size two-stage result strengthens the earlier
specified four-path screen; the earlier explicit calculation and its
literal phase evidence remain valid.

The perturbation floor $`1/4-2\epsilon`$ from Section 9 also applies
here. Interlayers that read or move the bank/address differently on the
two fibers are outside the theorem. So are sufficiently many completed
stage calls. This is a quantitative limit of the stated source-access
architecture, not an unrestricted native T-depth bound or a general
all-address separation theorem.

<a id="selector-crossing-rank"></a>

## 11. Mixed Fourier frames still owe enough selector crossing rank

Changing the Fourier frame can break the uniform ray in Section 9.
There is nevertheless a larger common fixed subspace to check. On the
marked origin fiber let $`\mathcal H_0=\{s_0=0\}`$, of dimension
D/2, with projectors P_0 and P_1 onto the s_0=0 and s_0=1 sectors.
Every original R_C fixes this entire subspace pointwise.

Consider the binary branch architecture
$`A_{j,0}=R_{C_j}`$, $`A_{j,1}=V_j^\dagger R_{C_j}V_j`$,
where V_j is any selector-only unitary. Define its crossing rank

```math
\rho_j=\mathrm{rank}(P_1V_jP_0).
```

Inside $`\mathcal H_0`$, intersect the kernels of these crossing
maps. The resulting subspace has dimension at least

```math
\dim\mathcal F\ge D/2-\sum_j\rho_j.                    \qquad(23)
```

For every v in this intersection, both v and V_jv lie in
$`\mathcal H_0`$. Thus every A_{j,b} fixes v, using the actual
V_j adjoint. If the right side of (23) is positive, all path products
fix a common marked vector and are identity on the empty bank. The
one-quarter Gram floor (19) follows for arbitrary path coefficients.
This proof concerns the conjugated-branch architecture; arbitrary extra
interlayers require their own analysis. Satisfying the opposite rank
inequality is necessary to escape this proof, not sufficient for acceptance.

There is also a useful bound when a dressing is dense but preserves
$`\mathcal H_0`$. Such a dressing's branch fixes that subspace.
For each remaining dressing use the defect rank r_C from (20).
The possibly larger common-fixed subspace
$`\mathcal F_{\rm fix}=\mathcal H_0\cap\bigcap_j\mathrm{Fix}(V_j^\dagger R_{C_j}V_j)`$ then satisfies

```math
\dim\mathcal F_{\rm fix}\ge
D/2-\sum_{j:\,V_j\mathcal H_0\ne\mathcal H_0}r_{C_j}.
```

At most two such remaining dressings leave a nonzero common subspace,
because each r_C is strictly below D/4. Hence **any number of stages
with at most two dressings that move this subspace fails**. In particular,
a three-stage repair with one computational-diagonal selector dressing
still has the quarter floor, even when its other two dressings are
arbitrary non-Clifford unitaries.

More strongly, suppose every V_j is a product over disjoint four-row
cells, with each cell factor unitarily conjugate to CCZ on three wires
and identity on the fourth. The conjugators within a cell may be
arbitrary; Pauli conjugations and partial Fourier walls are included. Only
the cell containing physical row zero crosses P_0 and P_1. Its factor
has rank two away from identity, so

```math
\rho_j\le2(D/16)=D/8.                                 \qquad(24)
```

The other cell factors are unitary and multiply this rank by their
dimension. Equivalently one may use disjoint three-wire factors, each
of defect rank one. Equations (23)--(24) leave dimension at least D/8
for **every three-stage source in this cell-factorized CCZ family**,
at every $`p\ge2`$. Its Gram error is at least one quarter. This
includes frames that all move $`\mathcal H_0`$ and destroy the
uniform-selector witness. Removing those simpler witnesses does not
remove this rank deficit.

Independent Pauli multipliers need not have rank-two cell defects and
are not included in this statement. Different cell partitions at the
three stages are allowed. More generally, crossing rank is subadditive:
insert $`I=P_0+P_1`$ between A and B to obtain
$`\rho(AB)\le\rho(A)+\rho(B)`$. Thus at most three such
cell-factorized CCZ layers in total across all dressings still give the
quarter floor, even with intervening walls that preserve
$`\mathcal H_0`$. Four layers merely make this estimate inconclusive.

The cell factorization is essential. Conjugating a whole growing CCZ
layer by a selector Clifford that spreads the physical-zero observable
across many cells can change (24), while keeping its controlled T-depth
constant. For example the already paid transpose-zeta involution N_p
has such support; $`V=N_pJN_p^\dagger`$ is controlled with two
ordinary N_p walls and the controlled J from Section 8. The walls cost
pm CNOTs per control and no T layers. This is an admissible escape to
test, not a source or a claim that its crossing rank is sufficient.
A nonlocal conjugation of only a single CCZ still has total defect rank
D/8 and does not escape (24). Changing several cells' joint action,
adding enough stages, or changing bank/address dynamics needs a new
complete acceptance calculation. These remain source-architecture
statements and do not bound unrestricted native T-depth.

<a id="nonlocal-selector-scout"></a>

## 12. A nonlocal repair escapes the rank budget but still fails at three bits

The distinction between escaping a screen and obtaining acceptance can
be checked exactly. At p=3 use matchings $`0\to1`$, $`1\to2`$
and $`2\to0`$, so every address coordinate is exposed. Repeat the
following frames in every four-row cell. Masks use the four selector
wires numbered by their cell corners, with corner zero least significant.
Here $`J_j=Q_jH_{A_j}P_jH_{A_j}Q_j^\dagger`$ and
$`Q_j=X^{u_j}Z^{v_j}`$; each P_j omits the indicated corner.

| Stage | Omitted corner | H mask | X mask | Z mask |
| --- | --- | --- | --- | --- |
| 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 15 | 2 | 1 |
| 3 | 3 | 3 | 4 | 5 |

All three frames move $`\mathcal H_0`$, and the first breaks the
uniform-selector ray. Before nonlocal wrapping, the selector vector
$`(|s=0\rangle+|s=2\rangle)/\sqrt2`$ nevertheless witnesses
the common fixed space. In its origin cell, the first CCZ is killed by
zero selector wires 2 and 3; in the second and third frames, H on the
selector-one $`|+\rangle`$ supplies a zero control. The specified
Pauli conjugations preserve these zero controls. Other cells cannot
change s_0. This local explanation works at every p at least two and
for any matching choices, in agreement with Section 11.

Now wrap the entire growing layer, $`V_j=N_pJ_jN_p^\dagger`$,
and use branches $`R_{C_j}`$ and $`V_j^\dagger R_{C_j}V_j`$.
At p=3, the exact crossing ranks are

```math
(\rho_1,\rho_2,\rho_3)=(48,48,48),\qquad
48>D/8=32,\qquad \sum_j\rho_j=144>D/2=128.
```

The cell-factorized rank obstruction is genuinely escaped. Yet all six
branches and their actual adjoints fix the different normalized vector

```math
v=(|s=68\rangle+|s=70\rangle)/\sqrt2.
```

This has selector wire 1 in $`|+\rangle`$, wires 2 and 6 in
$`|1\rangle`$, and every other wire zero. The complete origin
stage calculation retains literal phases and both borrowed-bit inputs.
Every empty-bank branch remains identity, so (19) still gives a quarter
Gram floor for this specified three-bit source at every signal angle.
This finite certificate does **not** exclude all nonlocal CCZ words or
a family patched at this size.

The next construction must change more than these frames or angles and
first establish complete marked cancellation, then normalized coverage
of the whole positive sector. Enough crossing rank is only a necessary
screen in the conjugated-branch architecture. More stages, different
joint selector operations or bank/address dynamics require a new exact
acceptance calculation and the same simultaneous resource/return ledger.
No general native-depth lower bound or improvement to R96 follows here.

Exact fixtures: [nonlinear native controls](../../tests/test_nonlinear_selector_control.py)
check phases, actual positive-S inverses, dirty return, counts and layers;
[coherent source ranks](../../tests/test_selector_source_rank.py) check
the origin/empty restrictions and an initialized-code witness with mixed
signal/private-work interlayers;
[mixed and nonlocal frames](../../tests/test_mixed_selector_acceptance.py)
check exact cell ranks, native cell words, the two displayed witnesses,
and the wrapped crossing ranks. These corroborate the all-size proofs
and the explicitly finite counterexample; they do not supply a source.
