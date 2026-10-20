# Three- and four-layer rigidity with a preserved address

[Two-layer proof and Clifford normal form](READ_ONLY_TWO_LAYER_RIGIDITY.md) ·
[Exact dirty queries](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Conditional frame transfer](QUERY_DEPTH_FRAME_TRANSFER.md)

A zero-clean read-only XOR query with a nonaffine output requires at
least four physical T layers at a fixed positive approximation radius,
regardless of dirty width. A literal Boolean scalar phase in the same
three-layer class has degree at most two. These are native statements
about address-preserving Clifford interlayers and work-only T gates;
they do not cover arbitrary circuits that return the address only at
the end. Section 6 establishes a separate one-clean conclusion:
three layers force Boolean XOR-query degree at most four at the same
approximation radius, and literal Boolean scalar degree at most five.
Section 8 proves that four layers with zero clean qubits force XOR-query
degree at most two at the smaller radius $`\sqrt2/12`$, and literal
Boolean scalar degree at most three. Its degree-two XOR bound is sharp.
Section 8.7 extends that four-layer conclusion to every component of a
word-valued XOR query, keeping all output wires arbitrary. This needs
an additional native block argument; preparing the other outputs would
not be a zero-clean proof.

## 1. Statement and exact four-band form

Let x be an a-bit address and let W consist entirely of arbitrary work,
including every output and dirty helper. Each completed Clifford
interlayer preserves the computational address exactly. A physical T
layer consists of T, T-dagger and identity on distinct work wires.
Empty layers are permitted. Set

```math
\rho=\frac{\sqrt2}{6}\approx0.235702.
```

**Theorem (zero clean, three T layers).** Suppose

```math
V=C_3D_3C_2D_2C_1D_1C_0,\qquad
\|V-Q_f\otimes I_{\rm dirty}\|_{\rm op}\le\eta\lt\rho,
\qquad Q_f|x,t\rangle=|x,t\oplus f(x)\rangle.                 \qquad(1)
```

Then f is affine over $`\mathbb F_2`$. The conclusion applies to
each component of a word-valued XOR query. If the target instead is
$`(-1)^{f(x)}I_W`$, then $`\deg_{\mathbb F_2}f\le2`$.
The operator norm includes every dirty input and arbitrary external
references. The theorem is independent of T-count, elementary Clifford
count and finite dirty width. A common global phase in the approximation
can be optimized without changing the conclusion.

Use the [read-only Clifford normal form](READ_ONLY_TWO_LAYER_RIGIDITY.md#1-address-preserving-clifford-interlayers)
at all four interlayers:

```math
C_j(x)=\lambda_j(x)P_j(x)C_j(0),\qquad
P_j(x)=X^{u_j(x)}Z^{v_j(x)},\qquad P_j(0)=I.                 \qquad(2)
```

The labels are binary linear. The complete scalar is retained:

```math
\lambda(x)=\prod_{j=0}^3\lambda_j(x)
 =i^{\sum_i\alpha_i x_i}
   (-1)^{\sum_{i\lt j}\beta_{ij}x_ix_j},\qquad\lambda(0)=1.
                                                                    \qquad(3)
```

Thus $`\lambda(x+y)/(\lambda(x)\lambda(y))`$ is a bilinear
sign, including the diagonal signs coming from odd $`\alpha_i`$.
All additions of addresses and binary labels below are over
$`\mathbb F_2`$.

Write $`U_x=V(x)`$, $`W_3=C_3(0)D_3`$, and define

```math
A_x=W_3P_2(x)W_3^\dagger,
\qquad
B_x=W_3C_2(0)D_2P_1(x)D_2^\dagger C_2(0)^\dagger W_3^\dagger,
\qquad
R_x=C_0(0)^\dagger P_0(x)C_0(0).
                                                                    \qquad(4)
```

Direct multiplication gives the phase-retaining identity

```math
U_x=\lambda(x)P_3(x)A_xB_xU_0R_x.                            \qquad(5)
```

Both A and B are exact projective representations of the additive
address group: each is a fixed unitary conjugate of a Pauli family
with linear labels. Their multiplication factors are bilinear signs.
In particular, this remains true when B and the actual zero-address
unitary $`U_0`$ are non-Clifford. No Clifford-hierarchy closure
assumption enters (4) or (5).

## 2. A robust Pauli-label cocycle

For target fibres $`F_x=X_t^{f(x)}\otimes I`$, solve (5) using
the actual inverse $`U_0^\dagger`$:

```math
A_xB_x
 =\lambda(x)^{-1}P_3(x)^\dagger U_xR_x^\dagger U_0^\dagger.
```

Replacing both target approximations gives a specified phase times
a Pauli,

```math
M_x=\lambda(x)^{-1}P_3(x)^\dagger F_xR_x^\dagger F_0^\dagger,
\qquad \|A_xB_x-M_x\|\le2\eta.                             \qquad(6)
```

Its binary label m satisfies $`m(0)=0`$ and

```math
m(x)=\ell(x)+e_{X_t}\bigl(f(x)+f(0)\bigr),                 \qquad(7)
```

where $`\ell`$ is linear. For multiple outputs, replace the
last term by the sum of their X components.

A physical T layer conjugates a Pauli to a Clifford. This assertion
does not concern products of arbitrary hierarchy elements: on one
wire, a Pauli commuting with Z stays Pauli, while conjugating X or Y
by T or T-dagger gives a Clifford with the symplectic action of S.
Tensoring the distinct wires and conjugating by $`C_3(0)`$
therefore shows that every $`A_y`$ is exactly Clifford.

Let $`S_A(y)`$ denote that symplectic action. The projective
multiplication laws of A and B give a bilinear sign $`\zeta(x,y)`$
such that the exact products $`N_x=A_xB_x`$ obey

```math
N_{x+y}=\zeta(x,y)A_yN_xA_y^\dagger N_y.                    \qquad(8)
```

Indeed, if $`A_xA_y=a(x,y)A_{x+y}`$ and
$`B_xB_y=b(x,y)B_{x+y}`$, then
$`\zeta(x,y)=a(y,x)b(x,y)`$. Comparing three instances of
(6) with (8) yields

```math
\|M_{x+y}-\zeta(x,y)A_yM_xA_y^\dagger M_y\|\le6\eta.
                                                                    \qquad(9)
```

Both operators in (9) are phases times Paulis. Distinct projective
Pauli labels have operator-norm distance at least $`\sqrt2`$,
even after optimizing the scalar phase. Since $`6\eta\lt\sqrt2`$,
their labels must agree:

```math
m(x+y)=S_A(y)m(x)+m(y).                                     \qquad(10)
```

The estimates use full work unitaries. They neither compress the
actual inverse to an occupied sector nor assume exact return of an
approximate implementation.

## 3. Independent physical T axes force additivity

An affine symplectic matrix in (10) alone would give only a degree-two
bound. The stronger native fact is the independence of the physical T
axes in $`D_3`$.

Change the Pauli-label coordinates by the fixed Clifford
$`C_3(0)^\dagger`$. Write a label as $`(p,q)`$, with p its
X part, and let $`u_k(y)`$ be the X component of $`P_2(y)`$
on an active T or T-dagger wire k. Each $`u_k`$ is linear.
In these coordinates,

```math
S_A(y)(p,q)=(p,q+u(y)\mathbin{\odot}p),                     \qquad(11)
```

with zero entries in u on inactive wires. The signs of the T gates
do not change this binary symplectic action.

Let $`p_k(x)`$ be the corresponding X component of m(x);
linearity of this function is not presumed. Set y=x in (10), using
$`m(0)=0`$, and then compare (10) with its x,y reversal. The
independent Z components in (11) give, separately on each active wire,

```math
u_k(x)p_k(x)=0,\qquad
u_k(y)p_k(x)=u_k(x)p_k(y).                                  \qquad(12)
```

If $`u_k`$ is nonzero, choose $`y_0`$ with
$`u_k(y_0)=1`$. The first identity gives $`p_k(y_0)=0`$;
the second then gives $`p_k(x)=0`$ for every x. If
$`u_k=0`$, that wire contributes nothing. Consequently

```math
S_A(y)m(x)=m(x)\quad\text{for every }x,y,
\qquad m(x+y)=m(x)+m(y).                                   \qquad(13)
```

Thus m is linear, and (7) proves that f is affine. This coordinate
argument is essential. For example, the generic symplectic family
$`S(x)=\mathrm{CX}_{1\to3}^{x_1}\mathrm{CX}_{2\to3}^{x_2}`$
has the nonlinear cocycle

```math
m(x)=x_1e_{X_2}+x_2e_{X_1}+x_1x_2e_{X_3},
\qquad m(x+y)=S(y)m(x)+m(y).                               \qquad(14)
```

It has overlapping CNOT image directions rather than the distinct
single-wire directions in (11). Equation (14) is a falsifier of
an argument that would deduce additivity from an arbitrary affine
symplectic action. It is not a three-T-layer query construction.

## 4. Literal scalar phases have degree at most two

Set $`\sigma_x=(-1)^{f(x)}`$ and $`F_x=\sigma_x I`$ in
(6). Then

```math
M_x=\mu_xQ_x,\qquad
\mu_x=\sigma_x\sigma_0\lambda(x)^{-1},\qquad
Q_x=P_3(x)^\dagger R_x^\dagger.                            \qquad(15)
```

The Pauli family Q has linear labels and a bilinear sign cocycle:

```math
Q_xQ_y=q(x,y)Q_{x+y},\qquad q(x,y)\in\{1,-1\}.             \qquad(16)
```

This follows directly by commuting the two linear-label endpoint
families. Fixed Clifford conjugation in R preserves its original
Pauli multiplication factors; no unknown phase is suppressed.

Equation (10), now with additive m, implies that conjugation by
$`A_y`$ fixes the label of $`Q_x`$. Therefore

```math
A_yQ_xA_y^\dagger=\chi(x,y)Q_x,\qquad
\chi(x,y)\in\{1,-1\}.                                   \qquad(17)
```

The sign $`\chi`$ is a character in each variable: in x this
uses the projective law (16), and in y it uses the projective law
of A. Hence $`\chi`$ is bilinear. Substituting (15)--(17)
into (9) gives

```math
|\mu_{x+y}-\zeta(x,y)\chi(x,y)q(x,y)\mu_x\mu_y|
 \le6\eta.                                                \qquad(18)
```

The ratio of the two unit-modulus scalars is a sign: the unknown
$`\sigma`$ factors are signs, while the polarization of lambda
in (3) is a bilinear sign. Opposite signs have distance 2.
Thus (18), already with $`6\eta\lt\sqrt2`$, forces exact
equality. It follows that

```math
(-1)^{f(x+y)+f(x)+f(y)+f(0)}
 =\zeta(x,y)\chi(x,y)q(x,y)
   \frac{\lambda(x+y)}{\lambda(x)\lambda(y)}.              \qquad(19)
```

The right side is bilinear. A Boolean function whose polarization
at zero is bilinear has degree at most two, proving the scalar claim.
This controls literal phases, including address S and CZ phases.
A common target phase cancels between $`F_x`$ and
$`F_0^\dagger`$ in (6), so the same proof permits optimizing
that single global phase.

## 5. Scope, resources and evidence

The result rules out at most three physical T layers for every
nonaffine zero-clean XOR query in the stated interface, even with
unbounded finite T-count, elementary Clifford count and simultaneous
dirty allocation. Affine XOR queries and quadratic scalar phases
already have Clifford implementations. Precision is explicit in (1);
the proof supplies an operator-norm floor $`\rho`$, independent
of address length and dirty width.

There is no clean qubit hidden in the coordinate changes or Pauli
comparisons; these are proof operations, not compiler resources.
The full-work argument does not directly apply to one initialized
qubit; section 6 handles its intermediate code excursions separately.
Neither theorem allows T gates on address wires after transforming
the address basis. The repository's
[unrestricted two-layer dirty Toffoli](../../docs/PARALLEL_DIRTY_LOOKUP.md#12-the-two-bit-indicator-has-optimal-t-depth-two)
therefore remains compatible with this result. Preserving the address
only at the end does not supply (2) at intermediate boundaries.

This is a constant-layer obstruction for a useful native query
architecture. It is neither a dimension-growing lower bound nor an
optimal-general-depth theorem for the complete Hopf frame. The
unrestricted one-clean problem and the genuinely sublogarithmic query
frontier remain open.

The bounded fixtures in
[`tests/test_read_only_three_layer_rigidity.py`](../../tests/test_read_only_three_layer_rigidity.py)
check the literal native four-band identity with non-Pauli $`U_0`$,
the actual inverse, address S/CZ phases, the six-error cocycle identity,
physical-axis additivity, and the nonlinear generic cocycle (14).
They are finite falsification tests; the width-independent statements
are proved above.

## 6. One clean qubit: degree at most four

Let $`J=|0\rangle_c\otimes I`$ initialize exactly one clean work
qubit, and let $`Z=Z_c`$, $`\mathcal K=\mathrm{im}J`$.
The output t and every dirty helper remain arbitrary. All four
Clifford interlayers preserve the address, and every physical T gate
acts on work, as in (1).

**Theorem (one clean, three T layers).** If

```math
\|VJ-J(Q_f\otimes I_{\rm dirty})\|_{\rm op}
\le\eta\lt\rho=\frac{\sqrt2}{6},
```

then $`\deg_{\mathbb F_2}f\le4`$. The conclusion is independent of
the number of dirty qubits, T-count and elementary Clifford count.
It includes every logical input, all dirty inputs, external references,
literal target phases and possible leakage in the initial error
hypothesis. It asserts a necessary condition, not attainability for
every degree-four table.

The proof obtains degree at most two on every coset of a linear
subspace of address codimension at most two. Its additional work is to
control intermediate clean-code excursions without assuming that an
individual dressed Pauli family preserves the clean sector.

### 6.1 Exact terminal return and normalization on each coset

The [three-layer clean-preservation gap](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#43-one-clean-qubit-exact-preservation-of-the-terminal-sectors)
applies to every fixed-address fibre $`U_x`$. Such a fibre is still
a work Clifford+T circuit of at most three physical T layers. Its
initialized leakage is at most eta, so

```math
\eta\lt\rho\lt\frac{\sin(\pi/8)}{\sqrt2}
\quad\Longrightarrow\quad [U_x,Z]=0
\quad\text{for every }x.                                  \qquad(20)
```

This imported repository lemma concerns terminal clean sectors only.
It gives no intermediate return and no restriction on the unitary in
the uninitialized sector.

In the four-band form (5), abbreviate $`P_x=P_3(x)`$. Let
$`s_P,s_R`$ be the binary linear clean-flip syndromes of the endpoint
Paulis, and set

```math
H=\ker s_P\cap\ker s_R,\qquad \mathrm{codim}H\le2.
                                                               \qquad(21)
```

Fix any address $`x_0`$. For $`h\in H`$, define

```math
\widehat A_h=P_{x_0}A_hP_{x_0}^\dagger,\qquad
\widehat B_h=P_{x_0}A_{x_0}B_hA_{x_0}^\dagger P_{x_0}^\dagger,
\qquad L_h=\widehat A_h\widehat B_h.
```

Projective linearity gives the full-space identities

```math
U_{x_0+h}\doteq P_hL_hU_{x_0}R_h,\qquad
L_h\doteq P_h^\dagger U_{x_0+h}R_h^\dagger U_{x_0}^\dagger.
                                                               \qquad(22)
```

Here and below $`\doteq`$ is equality up to a scalar of modulus one;
it is used only to obtain a necessary invariant and never changes the
literal approximation target. Both h-families remain projective Pauli
representations, so each member squares to a scalar.
Each $`\widehat A_h`$ is Clifford, with the same symplectic action
$`S(h)=I+K(h)`$ as $`A_h`$. In particular K is binary linear
in h and has the independent-axis form (11). Pauli conjugation by
$`P_{x_0}`$ changes signs, not this action.

Write $`F_x=X_t^{f(x)}\otimes I`$. The actual inverse satisfies

```math
U_x^\dagger J-JF_x^\dagger
=U_x^\dagger(JF_x-U_xJ)F_x^\dagger.
```

Replacing the rightmost inverse in (22) and then the forward call,
with the code-preserving $`R_h^\dagger`$ between them, proves

```math
\|L_hJ-\mu_hN_hJ\|\le\delta,\qquad \delta=2\eta,             \qquad(23)
```

for unit scalars $`\mu_h`$ and physical Paulis $`N_h`$ that
preserve $`\mathcal K`$. Their physical binary labels can be chosen as

```math
m(h)=r(h)+e\,g(h),\qquad
g(h)=f(x_0+h)+f(x_0),\qquad e=\mathrm{lab}X_t,         \qquad(24)
```

where r is linear. The clean-Z label of the endpoint Paulis is
included in r; it has not been dropped by compressing to logical
labels. Equation (20), the endpoint syndromes and (22) also imply
$`[L_h,Z]=0`$ exactly. All norm comparisons apply to initialized
inputs; leakage vectors in the two-call replacement are propagated by
unitaries.

For the rest of the proof suppress the hats on A and B. Their
projective laws imply the exact cocycle

```math
L_{h+k}\doteq A_kL_hA_k^\dagger L_k.                       \qquad(25)
```

### 6.2 The actual intermediate codes form a Pauli orbit

Since L preserves $`\mathcal K`$ and A and B are involutions up
to scalars, $`B_h\mathcal K=A_h\mathcal K`$. Equality of these
half-dimensional subspaces gives the exact equality of stabilizers

```math
Z_h:=A_hZA_h^\dagger=B_hZB_h^\dagger.                      \qquad(26)
```

Every $`Z_h`$ is a signed physical Pauli, because A is Clifford.
Both representations transport the orbit by translation:

```math
A_hZ_kA_h^\dagger=B_hZ_kB_h^\dagger=Z_{h+k}.
```

Consequently

```math
[L_h,Z_k]=0\quad\text{for every }h,k.                      \qquad(27)
```

The orbit labels have the affine form

```math
z_h=z+\ell(h),\qquad z=\mathrm{lab}Z,\qquad
\ell(h)=K(h)z.                                           \qquad(28)
```

No rank bound on ell is assumed. In particular, the collection of
intermediate codes may grow with dirty width.

### 6.3 An anticommuting orbit forces a full Pauli family

Suppose some $`Z_k`$ anticommutes with Z. Order the clean sectors
as 0 and 1. A Pauli of this form and a unitary commuting with Z have
blocks

```math
Z_k=\begin{pmatrix}0&R^\dagger\\R&0\end{pmatrix},
\qquad
L_h=\begin{pmatrix}C_h&0\\0&D_h\end{pmatrix},
```

where R is a scalar times a Pauli on the nonclean work. Equation
(27) forces $`D_h=RC_hR^\dagger`$. The occupied block in (23)
is within delta of a scalar times a logical Pauli n. Its other block
is therefore within the same delta of $`RnR^\dagger=\pm n`$.
Thus L is within delta, in full operator norm, of a physical Pauli
family $`\widetilde N_h`$. This step loses no additional error.

Three comparisons in (25) now give the full Pauli-label cocycle at
$`3\delta=6\eta\lt\sqrt2`$. Section 3 applies to these labels
even though they were not known in advance: the independent physical
T axes force them to be additive. Their target-X components agree
with (24), since the possible alternative extension changes only the
clean-Z label. Hence g is affine in this case.

### 6.4 Commuting orbit codes: compare on their intersections

It remains to treat the case in which every $`Z_k`$ commutes
with Z. They then commute pairwise: conjugating $`[Z_h,Z_k]`$
by $`A_h^\dagger`$ reduces it to $`[Z,Z_{h+k}]`$.
Let $`\mathcal K_k=A_k^\dagger\mathcal K=A_k\mathcal K`$,
the positive sector of $`Z_k`$.

Solving (25) for L_h and using (23) on the actual inverse of L_k
shows that its action on $`\mathcal K_k`$ is within $`2\delta`$
of a scalar times the Pauli

```math
Q_{h,k}=A_k^\dagger N_{h+k}N_k^\dagger A_k.                \qquad(29)
```

This estimate is valid on the whole of $`\mathcal K_k`$:
$`A_k`$ first maps that input to $`\mathcal K`$, and both
ideal replacement Paulis preserve $`\mathcal K`$. It does not
apply (23) to an arbitrary excursion outside its domain.

If $`Z_k\ne-Z`$, the commuting code intersection
$`\mathcal K\cap\mathcal K_k`$ is nonzero. Compare (23) and
(29) on that intersection. The Pauli isometries differ by at most
$`3\delta=6\eta`$. For a nonzero stabilizer code, two Pauli
restrictions whose quotient label is outside its stabilizer span
have distance at least $`\sqrt2`$, even after optimizing scalar
phase. Indeed the normalized trace of that quotient on the code is
zero, so the normalized Hilbert--Schmidt distance is $`\sqrt2`$.
This includes Pauli quotients that take the code to an orthogonal
syndrome.

Consequently, with
$`E_k=\mathrm{span}_{\mathbb F_2}\{z,z_k\}`$,
the direct label comparison gives
$`m(h)+S(k)[m(h+k)+m(k)]\in E_k`$. The involution S(k)
swaps z and $`z_k`$ and therefore preserves E_k; applying it gives

```math
m(h+k)+m(k)+S(k)m(h)\in E_k
\quad\text{for every }h,\quad\text{whenever }Z_k\ne-Z.
                                                               \qquad(30)
```

For identical positive sectors the stabilizer span has dimension one;
for independent commuting stabilizers it has dimension two. The
opposite-sector case is excluded explicitly because its intersection
is zero.

### 6.5 A bounded set of exceptional directions cannot hide degree three

Assume for contradiction that $`d=\deg g\ge3`$. Substituting
(24) into (30), using $`S(k)=I+K(k)`$, gives

```math
K(k)r(h)+K(k)e\,g(h)
 + e\bigl(g(h+k)+g(h)+g(k)\bigr)\in E_k.                  \qquad(31)
```

Fix a k for which (30) applies. The first term is linear in h.
Translation lowers the degree of
$`g(h+k)+g(h)`$ to at most $`d-1`$. Since a degree-d
coefficient of g is nonzero, comparison of degree-d coefficients
in the quotient vector space modulo E_k forces

```math
K(k)e\in E_k.
```

If additionally $`e\notin E_k`$, apply a linear functional
vanishing on E_k and taking e to one. Equation (31) then implies
that $`h\mapsto g(h+k)+g(h)`$ is affine.

There is a subspace G of H of codimension at most two on which both
requirements hold. Choose a binary linear functional lambda on
physical Pauli labels with
$`\lambda(z)=0`$ and $`\lambda(e)=1`$; the clean Z and
output X are independent. On $`\ker\ell`$, equation (28)
gives $`Z_k=(-1)^{\sigma(k)}Z`$. The sign sigma is linear
there because A is a projective representation. Extend it to any
linear functional $`\bar\sigma`$ on H and set

```math
G=\ker(\lambda\circ\ell)\cap\ker\bar\sigma,\qquad
\mathrm{codim}_H G\le2.                            \qquad(32)
```

For $`k\in G`$, lambda vanishes on E_k, so $`e\notin E_k`$.
If ell(k) is zero, sigma(k) is also zero and $`Z_k=Z`$;
otherwise the stabilizer labels are distinct and the two commuting
positive sectors have nonzero intersection. Thus every G-direction
derivative of g is affine.

Choose linear coordinates with G as the first set of coordinate
directions. A Boolean monomial of degree at least three involving
one of those coordinates would give a derivative of degree at least
two in that coordinate, contradicting the preceding conclusion.
Any remaining monomial uses only quotient coordinates, of which
there are at most two. Hence $`\deg g\le2`$, a contradiction.
This proves degree at most two on every fixed coset of H.

Finally choose at most two quotient coordinates for H in address
space. In the unique multilinear expansion, every coset restriction
has degree at most two in the H-coordinates, so every coefficient of
a higher H-monomial vanishes on all quotient assignments and is zero.
A surviving monomial contains at most two H-coordinates and at most
two quotient coordinates. Therefore $`\deg f\le4`$.

The two codimension bounds in this proof have different roles:
H comes from endpoint clean flips, while G controls exceptional
intermediate-code comparisons. Neither restricts the size of the
dirty register or assumes return inside a T-layer sequence. The
argument has not been extended to four or more physical T layers,
and does not yield a dimension-growing unrestricted depth bound.

### 6.6 Scalar corollary and bounded checks

The necessary degree bound in this section also holds for a specified
phase-decorated query $`e^{i\theta(x)}X_t^{f(x)}`$, with arbitrary
fibre phases. The proof still uses its literal target and actual
inverse; the phases enter $`\mu_h`$ in (23), while (24) and the
subsequent Pauli-label comparisons are unchanged.

Consequently, a literal scalar query $`P_g=(-1)^{g(x,t)}I`$
in the same one-clean, three-layer interface has degree at most five.
Treat any one address bit t as an arbitrary output and apply two
Hadamards:

```math
H_tP_gH_t
=(-1)^{g(x,0)}X_t^{\,g(x,0)+g(x,1)}.                     \qquad(33)
```

All remaining address wires still satisfy the required interlayer
condition. The clean and dirty allocation, physical width, T-count,
T-depth and error are unchanged; the elementary Clifford count
increases by two. There is no initialized output and no controlled
oracle assumption. The phase-decorated query bound shows that the
derivative of g in every individual coordinate has degree at most
four. Hence $`\deg g\le5`$. This scalar bound is only a necessary
condition and is not asserted to be sharp.

The bounded checks in
[`test_read_only_three_layer_clean_rigidity.py`](../../tests/test_read_only_three_layer_clean_rigidity.py)
exercise the nonzero stabilizer-intersection gap, the opposite-sector
exception, the anticommuting-orbit extension and the codimension-two
derivative argument. Section 7 supplies native examples that also test
the distinction between zero clean and one clean qubit.

## 7. Native layer fixtures and the clean-qubit distinction

The zero-clean bound is sharp for a two-bit AND query. With read-only
address bits a,b and arbitrary target t, the following exact native
circuits have these resources:

| Query | External clean | Additional arbitrary dirty | T-count | T-depth | Elementary Cliffords |
|---|---:|---:|---:|---:|---:|
| $`t\mapsto t\oplus ab`$ | 0 | one bit d | 8 | 4 | 18 |
| $`t\mapsto t\oplus ab`$ | one bit c | one bit d | 8 | 3 | 28 |

The target t is arbitrary in both rows and is not included in the
additional dirty count. There are no other workspace reservations.
Every T gate acts on work. Every Clifford interlayer preserves a,b
pointwise; indeed every elementary address participation is as a CNOT
control. Counts include both endpoint Hadamards and every actual basis
inverse. These are explicit upper counts without requiring cancellation
between successive Clifford words.

The three-layer lower theorem makes four layers optimal for this
**zero-clean read-only** query, also at its stated positive error radius.
The [exact one-clean two-layer phase theorem](READ_ONLY_TWO_LAYER_RIGIDITY.md#5-exact-literal-xor-queries-with-one-clean-qubit-are-affine)
also makes three layers optimal for the **exact one-clean read-only**
AND query. That theorem supplies no positive approximation radius;
the robust one-clean two-layer degree and rank conditions still permit
AND. Neither optimality statement gives a dimension-growing depth
lower bound.

### The common literal phase identity

For bits a,b,t,d, the integer identity

```math
\sum_{s\in\{0,1\}^3}(-1)^{s_1+s_2+s_3}
 (d\oplus s_1a\oplus s_2b\oplus s_3t)
 =-4abt(-1)^d\equiv4abt\pmod8
```

implements the literal CCZ phase using eight T or inverse-T gates.
Its phase is independent of the arbitrary dirty d. A Hadamard on t
before and after the complete diagonal word converts it to the XOR
query. The basis changes below are actual CNOT circuits, not parallel
copies of unavailable clean parity values.

### Four layers with no clean workspace

For each $`(s_1,s_2)\in\{0,1\}^2`$, compute

```math
(d,t)\longmapsto
(d\oplus s_1a\oplus s_2b,
 d\oplus t\oplus s_1a\oplus s_2b).
```

Chronologically, apply CNOT a-to-d if $`s_1=1`$, CNOT b-to-d if
$`s_2=1`$, and CNOT d-to-t. On these two distinct physical work
wires apply $`T^\sigma`$ to d and $`T^{-\sigma}`$ to t,
where $`\sigma=(-1)^{s_1+s_2}`$, then reverse that exact CNOT
word. These are the two terms with the specified $`s_1,s_2`$.
The four pairs give four physical T layers and eight T gates. The CNOT
count is $`2\sum_{s_1,s_2}(1+s_1+s_2)=16`$; adding the two
Hadamards gives 18 elementary Cliffords. The common identity proves
the full operator equality $`\mathrm{CCX}_{a,b;t}\otimes I_d`$.

### Three layers with one clean workspace bit

Initialize c to zero and replace each phase parity by

```math
p_s=d\oplus s_1(a\oplus c)\oplus s_2b\oplus s_3t.
```

The complete physical diagonal phase is now
$`(-1)^{(a\oplus c)bt}`$. The three listed groups have independent
work parity rows, and the following explicit CNOT computations realize
them. Gates and parity entries within each cell are chronological or
ordered as displayed. A notation such as d-to-t denotes one CNOT.

| Layer | Actual CNOT computation | Active parity wires | T signs |
|---|---|---|---|
| 1 | d-to-t, d-to-c, a-to-c | d; $`d\oplus t`$; $`d\oplus c\oplus a`$ on d,t,c | $`+,-,-`$ |
| 2 | a-to-c, c-to-t, d-to-t, d-to-c, b-to-d, b-to-c | $`d\oplus b`$; $`d\oplus t\oplus c\oplus a`$; $`d\oplus c\oplus a\oplus b`$ on d,t,c | $`-,+,+`$ |
| 3 | t-to-d, b-to-d, d-to-c, a-to-c | $`d\oplus t\oplus b`$; $`d\oplus t\oplus c\oplus a\oplus b`$ on d,c | $`+,-`$ |

In each row, apply the displayed physical T or inverse-T gates and
then reverse the **actual** CNOT computation before starting the next
row. Thus c,d,t have their original bit values at each completed row;
no temporary parity is assumed to be initialized. The eight terms are
exactly the common phase identity with a replaced by $`a\oplus c`$.
There are $`2(3+6+4)=26`$ CNOTs and two endpoint Hadamards.

On the entire five-qubit Hilbert space the resulting operator is

```math
|a,b,t,d,c\rangle\longmapsto
|a,b,t\oplus b(a\oplus c),d,c\rangle.
```

Consequently c initialized to zero returns exactly, d returns on every
occupied input, and the required Toffoli acts on arbitrary t. Literal
full-space matrix equality proves the claims also for superposed
addresses and work entangled with arbitrary external references.
There is zero approximation error and zero clean leakage.

The [native fixture tests](../../tests/test_read_only_query_layer_fixtures.py)
check both complete matrices with exact arithmetic, including the
one-clean circuit's c=1 extension, and audit physical T layers and all
elementary gate counts. The analytic identities supply the work-return
proof; finite fixtures do not establish a general lower bound.

### A two-layer candidate retains a non-Clifford address phase

There is a cheaper exact **relative-phase** construction, which fails
the literal query contract. With c the sole clean bit and no additional
dirty helper, define the fixed work unitary

```math
D=H_t\,\mathrm{CX}_{t,c}(T_t^\dagger T_c)
            \mathrm{CX}_{t,c}H_t
 =e^{i\pi X_t/8}e^{-i\pi Z_cX_t/8}.
```

The two Pauli rotation axes commute and are independent; the displayed
Clifford basis word makes D one physical work-only T layer. The native
T and inverse-T global phases cancel. On c=0, D is I; on c=1, it is
$`e^{i\pi X_t/4}`$. Consequently

```math
DZ_tD^\dagger=|0\rangle\langle0|_c\otimes Z_t
              +|1\rangle\langle1|_c\otimes Y_t.
```

Apply, chronologically,

```math
\mathrm{CX}_{a,c},\quad D^\dagger,\quad
\mathrm{CZ}_{b,t},\quad D,\quad
\mathrm{CZ}_{b,t},\quad\mathrm{CX}_{a,c}.
```

For c initialized to zero its exact action is

```math
|0\rangle_c|\psi\rangle_t\mapsto
 |0\rangle_c(-i)^{ab}X_t^{ab}|\psi\rangle_t.
```

It uses four T gates in two physical T layers and 16 elementary
Cliffords when each CZ is expanded as H, CNOT, H. It returns the clean
bit exactly and works on arbitrary t with references. Its unwanted
address phase is not leakage or an ignorable common scalar: the
address branch $`a=b=1`$ differs from the literal query by $`-i`$.

The missing correction is $`i^{ab}`$, a controlled-S phase on the
address. It is non-Clifford and cannot be absorbed into a Clifford
interlayer. A native correction and its combined physical layer schedule
would have to be supplied before treating this word as a literal query.
The separate exact two-layer phase theorem proves that no resynthesis
within this one-clean read-only interface can complete a literal AND
query in two T layers. It does not establish an approximation floor for
that initialized problem. The exact fixture retains the candidate's
phase and verifies the actual inverse and the full c=1 extension as well.

## 8. Four layers with zero clean qubits: degree at most two

The four-layer case needs an additional invariant because its
normalized middle is no longer a product of only two dressed Pauli
representations. Peeling its outer one-T band produces a Clifford
target family; an auxiliary Pauli cocycle then handles the case in
which the target is invisible to that family's symplectic action.

**Theorem (zero clean, four T layers).** Let every Clifford interlayer
preserve the computational address, every T gate act on work, and all
work inputs be arbitrary. If

```math
V=C_4D_4C_3D_3C_2D_2C_1D_1C_0,\qquad
\|V-Q_f\otimes I_{\rm dirty}\|_{\rm op}
\le\eta\lt\rho_4=\frac{\sqrt2}{12},
```

then $`\deg_{\mathbb F_2}f\le2`$. The error and conclusion hold on
the complete operator, including arbitrary dirty inputs and external
references. T-count, elementary Clifford count and finite dirty width
are unrestricted. The degree bound is sharp: section 7 supplies an
exact degree-two query using four layers. This does not say that every
quadratic query has such a compiler, or give a growing depth bound.
There is no one-clean four-layer conclusion here.
The word-valued extension at the same radius is proved in Section 8.7.

### 8.1 Peeling the outer band and quantizing a Clifford cocycle

Use (2) at all five interlayers. Put
$`W_4=C_4(0)D_4`$,
$`R_x=C_0(0)^\dagger P_0(x)C_0(0)`$ and
$`\lambda(x)=\prod_{j=0}^4\lambda_j(x)`$. The exact factorization is

```math
U_x=\lambda(x)P_4(x)W_4P_3(x)A_xB_xW_4^\dagger U_0R_x,
```

where

```math
A_x=C_3(0)D_3P_2(x)D_3^\dagger C_3(0)^\dagger,
```

```math
B_x=C_3(0)D_3C_2(0)D_2P_1(x)D_2^\dagger C_2(0)^\dagger
D_3^\dagger C_3(0)^\dagger.
```

A and B are exact projective Pauli representations, and A is a
one-T Clifford family. For $`F_x=X_t^{f(x)}\otimes I`$, define

```math
M_x=P_3(x)^\dagger W_4^\dagger
 [\lambda(x)^{-1}P_4(x)^\dagger F_xR_x^\dagger F_0^\dagger]W_4.
                                                               \qquad(34)
```

Solving the factorization and replacing $`U_x`$ and its actual
zero-address inverse gives
$`\|A_xB_x-M_x\|\le2\eta`$. Each M is a scalar times a Clifford,
since the bracket in (34) is a Pauli and the only non-Clifford
conjugation in $`W_4`$ is one physical T layer. The exact
representation laws, followed by three such replacements, give

```math
\|M_{x+y}-\zeta(x,y)A_yM_xA_y^\dagger M_y\|\le6\eta,
                                                               \qquad(35)
```

where zeta is a unit scalar fixed by those laws.

Cliffords with different symplectic actions have projective operator
distance at least $`1/\sqrt2`$: their conjugations of some Pauli
have different projective labels, separated by $`\sqrt2`$, while
conjugation changes by at most twice the original operator distance.
Thus $`6\eta\lt1/\sqrt2`$ forces the symplectic cocycle exactly.
After that, the quotient of the two sides in (35) is a scalar times
a Pauli. A nonzero Pauli label would again have distance at least
$`\sqrt2`$ from any scalar identity. Equation (35) therefore also
forces the corresponding Pauli-displacement cocycle exactly. This
does not assert equality of the remaining scalar phases.

### 8.2 A target visible in the outer shear has degree at most two

Put $`g(x)=f(x)+f(0)`$. The bracket Pauli in (34) has label
linear(x) plus $`e_{X_t}g(x)`$. After the fixed Clifford
$`C_4(0)^\dagger`$, its X components on active $`D_4`$ wires
determine the symplectic action

```math
S_M(x)=I+H(x),\qquad H(x)=L(x)+E g(x).                    \qquad(36)
```

Here L is linear, E is fixed, and the entire outer shear space has
pairwise zero matrix products. The left factor $`P_3^\dagger`$
has no symplectic effect. Write $`\mathcal C_y`$ for conjugation
of matrices by $`S_A(y)`$. The exact symplectic cocycle becomes

```math
H(x+y)+H(y)=\mathcal C_yH(x),
```

```math
E[g(x+y)+g(x)+g(y)]
=(\mathcal C_y-I)L(x)+(\mathcal C_yE-E)g(x).              \qquad(37)
```

If $`E\ne0`$ and $`\deg g\ge2`$, comparison of the
highest-degree x coefficients gives $`\mathcal C_yE=E`$ for
every y. Equation (37) then makes every derivative of g affine.
Hence $`\deg g\le2`$. Affine g already satisfies the conclusion.

### 8.3 Vanishing outer shear and an auxiliary Pauli cocycle

Suppose now that $`E=0`$. The target Pauli after
$`C_4(0)^\dagger`$ commutes with every active $`D_4`$ axis.
Consequently $`H=W_4^\dagger X_tW_4`$ is a fixed nonidentity
physical Pauli, with label v. Choose a linear-label
Pauli family $`T_x`$, independent of g, so that projectively

```math
M_x\doteq Q_xH^{g(x)},\qquad
Q_x=P_xN_x,\qquad P_x=P_3(x)^\dagger,\qquad
N_x=D_4^\dagger T_xD_4.                                 \qquad(38)
```

The label of $`T_x`$ includes the fixed $`C_4(0)^\dagger`$
coordinate change and the linear endpoint labels. Scalars depending
on x or g are retained in the actual target and may be ignored only
in this necessary projective invariant. N is an exact one-T Clifford
Pauli representation; Q need not be a representation. Write

```math
p(x)=\mathrm{lab}P_x,\qquad
S_N(x)=S_Q(x)=I+L(x),\qquad L(x)L(y)=0.                  \qquad(39)
```

Both p and L are linear. The symplectic cocycle, now independent of
g, implies

```math
[S_A(y),S_N(x)]=0\quad\text{for every }x,y.               \qquad(40)
```

Fix y and abbreviate $`A=A_y`$. The operator

```math
C_x=N_x^\dagger A N_x A^\dagger
```

is a scalar times a Pauli by (40). The family
$`N'_x=A N_xA^\dagger=N_xC_x`$ is another exact projective
representation with the same multiplication factors as N. Comparing
their laws gives

```math
c(x+z)=S_N(z)^{-1}c(x)+c(z),\qquad
c(x)=\mathrm{lab}C_x,\quad c(0)=0.                       \qquad(41)
```

Section 3's independent physical-axis argument applies to (41).
The axes here belong to $`D_4`$, and their binary action is
$`(a,b)\mapsto(a,b+u(z)\odot a)`$. The diagonal and symmetric
cocycle equations force each active product
$`u_i(z)c_{X,i}(x)`$ to vanish. Therefore c is additive, hence
linear in x. This step does not assume a quadratic phase function for
arbitrary Clifford tableaus.

### 8.4 The remaining Pauli defect is quadratic

For the fixed y, put $`K_y=I+S_A(y)`$. The commutator of Q
with A has Pauli label

```math
k_y(x)=\mathrm{lab}(Q_x^\dagger A Q_xA^\dagger)
      =S_N(x)K_y p(x)+c(x),                             \qquad(42)
```

because its full operator factors as

```math
Q_x^\dagger A Q_xA^\dagger
=N_x^\dagger[P_x^\dagger A P_xA^\dagger]N_x\, C_x.
```

Thus k has degree at most two in x. Separately,

```math
Q_xQ_y\doteq P_{L(x)p(y)}Q_{x+y},\qquad
\mathrm{lab}(Q_{x+y}^\dagger Q_xQ_y)=L(x)p(y),             \qquad(43)
```

where $`P_l`$ denotes a Pauli of label l. The second identity
uses $`L(x+y)L(x)=0`$. The full baseline defect is the Pauli

```math
D_{x,y}=Q_{x+y}^\dagger A_yQ_xA_y^\dagger Q_y,
```

```math
d(x,y)=\mathrm{lab}D_{x,y}
      =L(x)p(y)+S_N(y)^{-1}k_y(x).                       \qquad(44)
```

For each fixed y, d has degree at most two in x, independently of
width. Moreover, its degree-two coefficient vectors lie in
$`\mathrm{im}K_y`$. Its only possible quadratic contribution is

```math
\bigl[S_N(y)^{-1}L(x)K_y p(x)\bigr]^{[2]}.               \qquad(45)
```

Equation (40) makes $`K_y`$ commute with every $`L(x)`$
and with $`S_N(y)`$, so the expression in (45) is in that image.

### 8.5 Independent inner axes exclude cubic degree

Substitute (38) into the exact Pauli-displacement cocycle from (35).
It gives

```math
v[g(x+y)+g(x)+g(y)]
=d(x,y)+[S_N(y)^{-1}S_A(y)v+v]g(x).                     \qquad(46)
```

If $`\deg g\ge3`$, highest-degree x coefficients force the
bracket multiplying g to vanish for every y, since d is at most
quadratic and translation lowers degree. Then every derivative has
degree at most two, so $`\deg g\le3`$.

It remains to exclude degree exactly three. Since $`E=0`$,
the target label v has no X component on an active outer axis;
hence $`S_N(x)v=v`$ for every x. Together with (40), the
vanishing bracket also gives $`S_A(y)v=v`$. Taking the
degree-two x parts of (46) and using (45) yields

```math
v\,(\Delta_y g)^{[2]}
\in\mathrm{im}(I+S_A(y))
\quad\text{coefficientwise for every }y.                \qquad(47)
```

In the fixed physical-axis coordinates of the one-T family A,

```math
S_A(y)(a,b)=(a,b+u(y)\odot a),\qquad
\mathrm{im}(I+S_A(y))
=\mathrm{span}\{z_i:u_i(y)=1\},                         \qquad(48)
```

with independent axes $`z_i`$ and linear functions $`u_i`$.
If v is outside the span of all axes with nonzero $`u_i`$,
(47) forces all derivatives to be affine. Otherwise choose a nonzero
axis coefficient of v, say its i-th coefficient. For every
$`y\in\ker u_i`$, v lies outside the image in (48), again making
$`\Delta_y g`$ affine. This is a codimension-one subspace of
directions. In coordinates adapted to it, a cubic monomial involving
a subspace coordinate has a quadratic derivative, while a monomial
using only the one quotient coordinate cannot be cubic. Thus cubic
degree is impossible, completing the degree-two proof.

### 8.6 Scalar corollary and scope

The projective proof also permits a specified phase-decorated target
$`e^{i\theta(x)}X_t^{f(x)}`$: its literal phases enter the scalar
in (38), without changing the Pauli-label comparisons. Consequently
a literal zero-clean scalar phase $`P_g=(-1)^{g(x,t)}I`$ at
four layers has Boolean degree at most three. Treat one address bit
t as an arbitrary output and use

```math
H_tP_gH_t=(-1)^{g(x,0)}X_t^{g(x,0)+g(x,1)}.
```

The two Hadamards add two elementary Cliffords; T-count, T-depth,
physical width, clean/dirty helper allocation and error are unchanged.
There is no prepared output or free controlled-query call. Each
coordinate derivative of g has degree at most two, giving the stated
scalar bound. No sharpness or universal construction is asserted for
this scalar corollary in the original argument. The subsequent
[phase-relocation construction](READ_ONLY_PHASE_RELOCATION.md) attains
scalar CCZ with four layers and zero clean work, making the scalar
bound sharp; it does not synthesize every cubic phase at the same
count or width.

The proof exploits native one-T axes twice, including when the target
is invisible to the outer symplectic action. It therefore resolves
that degeneracy rather than treating a visible shear as a general
invariant. It still requires address-preserving Clifford boundaries
and zero initialized work. It does not supply a one-clean four-layer
bound, a recurrence for arbitrary T-depth, or an unrestricted
complete-frame lower bound.

The [four-layer evidence tests](../../tests/test_read_only_four_layer_rigidity.py)
check the literal five-band normalization and actual inverse on all
address fibres of a native eight-T, four-layer word with 101 elementary
Cliffords, zero clean qubits and two arbitrary work qubits. A second
native fixture has commuting inner and outer shears with nonzero shear
products; it verifies the auxiliary additive cocycle, a genuinely
quadratic defect and the image constraint. Those matrix checks use
floating arithmetic with explicit residual tolerances. A third test
exhausts all three-variable Boolean tables and codimension-one
direction spaces for the derivative step. These finite checks support
the analytic, width-independent theorem; they do not replace it.

### 8.7 Word-valued queries: the outer shear is additive

The same four-layer hypotheses and radius also imply

```math
\deg f_j\le2\quad\text{for every component of }f,
\qquad Q_f|x,t\rangle=|x,t\mathbin\oplus f(x)\rangle.
                                                               \qquad(49)
```

All output bits t are arbitrary work. No other output is prepared,
discarded or absorbed into a clean allocation. A specified scalar
multiplying the target is also permitted: the following necessary
projective invariants do not alter the literal error promise. The
single-output proof above does not itself imply (49), because its
visible-shear case used one fixed target direction.

**An additive outer-shear lemma.** The normalization (34)--(35) works
unchanged for $`F_x=X^{f(x)}`$. Let $`z_1,\ldots,z_n`$ be the
independent active outer T-axis labels in its fixed Clifford coordinates,
and let Omega be the symplectic matrix. Write

```math
H(x)=S_M(x)-I=\sum_i h_i(x)z_i z_i^{\mathsf T}\Omega,
\qquad H(0)=0.
```

The robust Clifford comparison already gives exactly

```math
H(x+y)=S_A(y)H(x)S_A(y)^{-1}+H(y),\qquad S_A(y)=I+K_y,
                                                               \qquad(50)
```

where K is linear and $`K_yK_z=0`$. We prove that h is linear.
Put $`W=\mathrm{span}\{h(x):x\}`$ and delete coordinates
identically zero on W. The zero W case is immediate. Independence of
the z_i makes the image of H(x) exactly their selected coordinate span.
Thus $`L=\sum_x\mathrm{im}H(x)=\mathrm{span}\{z_i\}`$ is
preserved by S_A(y). Its restriction O_y in the z basis and the induced
automorphism phi_y of W obey

```math
\begin{aligned}
O_y\,\mathrm{diag}(w)\,O_y^{\mathsf T}
 &=\mathrm{diag}(\phi_yw)\quad(w\in W),\\
h(x+y)&=\phi_yh(x)+h(y),\qquad O_{y+z}=O_y+O_z+I.
\end{aligned}                                                  \qquad(51)
```

Both O and phi are representations of the additive address group.
Group the coordinates by their nonzero coordinate functionals
$`\chi_i\in W^*`$, obtaining independent blocks $`L_\chi`$.
The congruence in (51) forces O_y to map $`L_\chi`$ onto
$`L_{\chi\circ\phi_y^{-1}}`$. To see this without presuming
an individual-axis permutation, choose a basis of W and formal
variable lists s,t. Over their rational-function field,
$`D(t)=\mathrm{diag}(\chi_i(t))`$ is invertible, and (51) gives

```math
O_yD(s)D(t)^{-1}O_y^{-1}
 =D(\phi_y s)D(\phi_y t)^{-1}.                          \qquad(52)
```

Distinct nonzero binary linear forms have distinct ratios
$`\chi(s)/\chi(t)`$: equality of cross products makes them
proportional, hence equal over the binary field. The eigenspaces in
(52) are therefore exactly the stated blocks. This field extension is
only an algebraic proof device; it supplies no physical gate or state.

Every block orbit has size at most two. Otherwise the elementary
abelian address group has y,z taking one block to four distinct blocks
under 0,y,z,y+z. For nonzero v in the first block, the identity
$`O_{y+z}v=O_yv+O_zv+v`$ would relate four nonzero vectors in
independent blocks. Thus each block is fixed or belongs to a two-cycle
swapped by a nonzero linear character u(y).

On a two-cycle chi,psi, put $`a(x)=\chi(h(x))`$,
$`b(x)=\psi(h(x))`$ and $`v(x)=a(x)+b(x)`$. Summing the
two cocycle equations makes v linear; the first equation is

```math
a(x+y)=a(x)+a(y)+u(y)v(x).
```

Setting y=x gives $`u(x)v(x)=0`$ for every x. A nonzero linear
u forces the linear v to vanish: choose x_0 with u(x_0)=1, then
compare x_0 and x_0+x for $`x\in\ker u`$. Hence a,b are
additive. Fixed-block coordinates are additive directly. All h
coordinates are therefore linear. Since the h(x) span W, (51) also
implies phi_y is the identity and

```math
[S_A(y),S_M(x)]=0\quad\text{for every }x,y.               \qquad(53)
```

**Remove the linear visible part.** Put $`g(x)=f(x)+f(0)`$.
The outer coefficient word has the form
$`h(x)=\ell(x)+E g(x)`$, with ell linear and E extracting the
active outer X coordinates of the target Pauli. Since h is linear,
choose a linear lift g_0 with $`E g_0(x)=E g(x)`$. The residual
word $`g+g_0`$ lies in the kernel of E. Its outer conjugate is
therefore a Pauli with label v(x), fixed by every outer shear.

Fold g_0 into the linear-label baseline T_x in (38). Then

```math
M_x\doteq Q_xP_{v(x)},\qquad Q_x=P_xN_x,
\qquad S_N(x)v(z)=v(z)\quad\text{for every }x,z.          \qquad(54)
```

Equations (39)--(45) apply to this baseline: (53) supplies the required
symplectic commutation, and d(x,y), for fixed y, has degree at most two
in x with every degree-two coefficient in $`\mathrm{im}K_y`$.
These are factorizations for a necessary invariant, not additional
synthesized bands or an altered literal target. The exact Pauli-label
comparison in (35) becomes

```math
v(x+y)+v(x)+v(y)=d(x,y)+K_yv(x).                         \qquad(55)
```

**The residual vector is quadratic.** Let J be the span of the
independent inner native axes with nonzero characters in (48).
Every image of K_y lies in J, and every K_y kills J. Modulo J, the
quadratic coefficients of d disappear. Equation (55) then makes every
derivative of v modulo J affine, so v modulo J has degree at most two.
Every higher-degree coefficient of v lies in J and is killed by K_y.
The right side of (55) consequently has degree at most two in x,
proving $`\deg v\le3`$.

Write the homogeneous cubic part in the independent inner-axis basis
as $`v^{[3]}(x)=\sum_i z'_i c_i(x)`$. Taking degree-two parts
of (55), including the term K_y applied to the quadratic part of v,
gives coefficientwise

```math
(\Delta_yv^{[3]})^{[2]}\in\mathrm{im}K_y
 =\mathrm{span}\{z'_i:u_i(y)=1\}.                       \qquad(56)
```

For $`y\in\ker u_i`$, this makes every derivative of c_i
affine. Section 8.5's codimension-one argument excludes each cubic
c_i. Thus v has degree at most two. The target Pauli labels are
independent, so the residual output word maps injectively to v;
adding back its linear lift proves (49) for every component.

All quantization used only the original $`6\eta\lt1/\sqrt2`$
comparison, so the extension loses no error. It retains zero clean
work, arbitrary finite dirty width, all logical columns, references,
actual inverse and complete return. It remains a four-layer theorem
for read-only Clifford boundaries and work-only T gates, not an
unrestricted depth bound or a general layer recurrence.

The same [four-layer fixture file](../../tests/test_read_only_four_layer_rigidity.py)
now runs two retained dirty AND words together for outputs ab and ac.
It uses three shared address bits, two arbitrary outputs and two arbitrary
dirty helpers, with 16 T gates, four layers and 36 elementary Cliffords.
Exact checks retain all 128 physical columns and the actual inverse;
both outputs and helpers are arbitrary throughout. This checks the
simultaneous native contract, not the universal block lemma above.
