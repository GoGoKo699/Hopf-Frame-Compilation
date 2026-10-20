# An enlarged actual center for one all-active pair without middle singletons

[Original physical pairs](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) ·
[All-active recognition](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md) ·
[Shared physical-clean feature proof](OCCUPIED_BOOLEAN_PAIR_CENTER.md#8-a-proved-bounded-span-reopening-with-shared-feature-coordinates)

The later [uniform all-active closure](ALL_ACTIVE_SINGLETON_CLOSURE.md)
also permits arbitrary original middle singletons and their unbounded
coupling to both silent lists, with sign/XOR bound 37748865/37748864.
It sharpens the general constant below while preserving this enlarged
center proof. The one-sided 2428/2427 and rank-two six/five bounds
remain sharper where applicable; general original pairs remain open.

Products of actual silent axes can provide transported central Paulis
without providing new native phase axes. In the one-pair case with no
middle singletons, Boolean endpoint commutation makes those products
large enough to capture every growing silent direction except one at
each endpoint. The entire original axis span then has only four
directions modulo this enlarged center. This closes a complete
larger-support subclass without bounding its old coupled-silent span.

## 1. Original-list hypotheses and complete conclusion

Let

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,\qquad
P_A=K_AZ_cK_A^\dagger,\qquad P_C=K_CZ_cK_C^\dagger
```

be an actual one-clean three-native-band circuit. Each E is an
independent commuting signed-Pauli T/T-dagger list. Assume exact clean
return and occupied action

```math
G=e^{i\theta}\mathrm{diag}((-1)^{g(x)})
```

on every logical and dirty input, with one retained common scalar and
arbitrary references. The unused clean block is unrestricted and is
retained throughout. Suppose the original middle list consists of
**exactly one genuine physical pair B,D and no singleton axes**, and
its difference $`r=b+d`$ is the sum of all actual clean-active axes
at each endpoint. These support conditions refer to the original lists,
not to a chosen basis of their spans.

There is no assumption that either endpoint silent list commutes with
the pair. Both active ranks and both actual silent spans may grow.
Define

```math
D_{\rm pair}=576\cdot2^{62}+2252.
```

Then

```math
\deg g\le D_{\rm pair},\qquad
\deg f\le D_{\rm pair}-1
\quad\text{for a complete Boolean XOR query }Q_f.               \qquad(1)
```

Two separate subclasses retain sharper bounds. If either silent list
does not couple to the pair, the
[all-active recognition theorem](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md#1-statement-in-the-original-physical-lists)
gives sign/XOR bounds 2428/2427. An actual full physical shared flip,
including the zero-cross-pairing case below, gives five/four. The known
rank-two theorem remains sharper in its own overlap. Equation (1) is
uniform in width and concerns this source class; it is not an
unrestricted depth bound.

## 2. Literal orientation and the common active span

Put $`J=E_AP_AE_A^\dagger`$. Write the actual pair relation as
$`JBJ=sD`$, with s a sign, and let epsilon,delta be its original T
orientations. Use the literal repair

```math
F=T_B^\epsilon T_D^{s\epsilon},\qquad
C=T_D^{\delta-s\epsilon},\qquad E_B=CF.
```

C is Clifford and F commutes J exactly. Absorb C into the output
encoder and conjugate its actual output list, as in the
[original physical insertion identity](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md#1-keep-the-actual-middle-axes-and-their-signs).
Use corrected C-side notation from now on. The correction fixes B,D
and r, retains the actual native lists and all scalar phases, and gives

```math
E_CP_CE_C^\dagger=J,
\qquad V=K_C^\dagger E_C^\dagger F E_AK_A.          \qquad(2)
```

Intrinsic independent-native displacement images imply
$`L_A=L_C=L`$. The complete Pauli supports of the equal clean images
are $`p_A+L`$ and $`p_C+L`$, so $`p_A+p_C\in L`$. Independence
and the all-active difference hypothesis give

```math
\langle b,a\rangle=1\quad(a\in A_{\rm active}),\qquad
\langle b,l\rangle=\langle p_A,l\rangle\quad(l\in L).
                                                               \qquad(3)
```

Both silent spans commute L and P_A. For the output span this follows
from its commutation with L and P_C together with
$`p_A+p_C\in L`$. Consequently both silent spans commute r, and
their coupling to D equals their coupling to B. The genuine central
subspace supplied by the occupied Boolean theorem is

```math
R=L\cap p_A^\perp.
```

It has codimension one in L and centralizes all three actual lists
and both clean Paulis. Its genuine transport follows from the
[Boolean endpoint commutation lemma](OCCUPIED_BOOLEAN_PAIR_CENTER.md#1-actual-source-and-the-endpoint-commutation-lemma),
not merely from an intrinsic change of coordinates.

## 3. Products of silent axes force a rank-one cross matrix

Let $`u_1,\ldots,u_m`$ and $`v_1,\ldots,v_n`$ be the actual
corrected A-silent and C-silent labels. Define

```math
\alpha_i=\langle u_i,b\rangle,\qquad
\gamma_j=\langle v_j,b\rangle,\qquad
K_{ij}=\langle u_i,v_j\rangle.                    \qquad(4)
```

If alpha or gamma is zero, the sharper zero-coupling theorem already
applies. Assume both are nonzero. For any coefficient vector x with
$`\alpha^{\mathsf T}x=0`$, the literal product of the corresponding
commuting A-silent Paulis has label $`u=\sum_i x_i u_i`$. It commutes
P_A and every actual A axis. It commutes B,D as well. Because there
are **no middle singletons**, it commutes the entire middle list.
The Boolean endpoint lemma therefore makes it commute every actual
C axis. Thus

```math
\ker\alpha^{\mathsf T}\subseteq\ker K^{\mathsf T}.
```

Apply the same argument to the actual reversed circuit and products
of C-silent axes. This gives

```math
\ker\gamma^{\mathsf T}\subseteq\ker K.
```

Every column of K is therefore a multiple of alpha and every row a
multiple of gamma. Since both are nonzero,

```math
K=\lambda\alpha\gamma^{\mathsf T},\qquad
\lambda\in\mathbb F_2.                           \qquad(5)
```

These products are actual Hermitian Paulis and actual transport
identities. They are not a replacement native phase list, and no T
phase about a product is inserted or deleted. The argument keeps the
independent physical lists from (2).

## 4. The physical shared-flip alternatives

Let $`U_A,U_C`$ be the two silent spans and let

```math
M=L+\langle b\rangle+U_A+U_C                      \qquad(6)
```

be the entire corrected physical phase-axis span. If
$`p_A\notin M`$, symplectic nondegeneracy supplies a label h
orthogonal to M with $`\langle h,p_A\rangle=1`$. It also pairs
to one with p_C. The associated signed Hermitian Pauli H commutes
every actual phase axis and flips both endpoint clean Paulis.

The actual clean-controlled Pauli endpoint dressings fix every
clean-zero input and turn the two transported flips into X_c.
Together with exact clean return this gives a full-space
$`I_c\otimes G`$ extension, with the same three native bands after
the explicit orientation repair. The zero-clean theorem gives sign/XOR
degree five/four.

If $`\lambda=0`$, this alternative is forced. Otherwise a proposed
expression for p_A in (6), paired with any actual A-active axis, would
have coefficient one on b by (3). It would read

```math
p_A=b+l+u+v,\qquad l\in L,\quad u\in U_A,\quad v\in U_C.
```

Pair it with an A-silent axis whose alpha entry is one. The left side
is zero; the terms l,u vanish, the term v vanishes because K=0, and
the term b is one. This contradiction proves $`p_A\notin M`$.
For lambda one, the same shared-flip conclusion is still available
whenever p_A lies outside M. The following center argument does not
require that additional membership decision.

## 5. An enlarged actual center with a four-dimensional quotient

Assume $`\lambda=1`$. Define the physical silent kernel spaces

```math
U_A^0=U_A\cap b^\perp,\qquad
U_C^0=U_C\cap b^\perp,
\qquad Z=R+U_A^0+U_C^0.                           \qquad(7)
```

The original silent lists are independent, so these are precisely
the physical images of the coefficient kernels used in Section 3.
Each member of $`U_A^0`$ commutes both clean Paulis and every actual
phase axis by the endpoint lemma. The reversed statement holds for
$`U_C^0`$. Therefore Z is isotropic, centralizes all of M, and has
an actual full-physical-space transport between the endpoint encoders.
No bound on its dimension is assumed. It can absorb an arbitrarily
large part of either actual silent span.

Choose one actual A-active label a and one coupled actual silent label
u and v from the respective lists. Since L/R and each silent/kernel
quotient have dimension one, a,b,u,v span M modulo Z. Their binary
pairing matrix is

```math
\begin{pmatrix}
0&1&0&0\\
1&0&1&1\\
0&1&0&1\\
0&1&1&0
\end{pmatrix}.                                   \qquad(8)
```

It is nonsingular over $`\mathbb F_2`$. Since Z centralizes M,
the four representatives are independent even modulo Z. Hence

```math
\dim(M/Z)=4.                                    \qquad(9)
```

Also $`p_A\notin Z`$: every member of Z commutes a, whereas P_A
anticommutes it. Thus projecting Z onto the positive P_A code is
injective. The analogous statement holds at the output.

This is a bound on the entire physical span **modulo a newly proved
transport center**. It is not a claim that the earlier additional
coupled-silent span t is bounded. That span may grow relative to the
smaller old center R.

## 6. Keep the physical clean coordinate through the feature proof

Use the retained exact active localization with

```math
S_0=L+\langle b\rangle\subseteq M,\qquad
C_0\in\mathrm{Alg}(S_0),\qquad
G_0=C_0^\dagger E_{C,\rm active}^\dagger F E_{A,\rm active},
\qquad [G_0,P_A]=0.
```

C_0 is the actual signed Clifford sending P_A to P_C. It fixes Z
literally, because Z centralizes S_0. It also preserves M: conjugation
by an operator in $`\mathrm{Alg}(M)`$ preserves that Pauli algebra,
and the Clifford C_0 sends each Pauli to an actual Pauli. In particular
the transformed C-silent axes remain actual independent commuting axes
inside M. With $`K_C'=C_0^\dagger K_C`$, the exact physical word is

```math
V=(K_C')^\dagger(E_{C,s}')^\dagger G_0E_{A,s}K_A,\qquad
E_{C,s}'=C_0^\dagger E_{C,s}C_0.                  \qquad(10)
```

Every internal factor preserves P_A. Z commutes all of them and has
injective logical projection. Apply the proved global Boolean
Pauli-transport normalization to this enlarged Z: extract its quadratic
Boolean Clifford, make fixed affine changes, and use Hadamards only
on factors where the residual target is already identity. Z becomes
a shared computational-Z subgroup. This is the same normalization
as in the [occupied-center proof](OCCUPIED_BOOLEAN_PAIR_CENTER.md#4-computational-center-and-a-finite-feature-description),
now with the larger actual subgroup whose transport was proved above.
No coordinates are chosen separately on its eigensectors.

Include the clean-preserving lift of this output normalization in K_C'.
Conjugate the entire physical M by that output encoder and take its
logical computational X masks, omitting only the clean coordinate.
The linear X-mask map annihilates Z. By (9), its image has rank
$`k\le4`$. One fixed affine change confines it to k logical core
bits. Let $`\mathcal M`$ be the computational-Z subgroup on the
remaining logical bits. Every actual source axis in M commutes
$`\mathcal M`$ on the full physical space; C_0 does too.

Let E_* be the actual transformed C-silent subband whose axes couple
to $`S_0\cap p_A^\perp`$. Cancel only the passive phases in the
dressing, with their actual inverses. The full-logical-space factors
are then

```math
\begin{aligned}
F_{\rm out}&=(K_C')^\dagger E_*^\dagger G_0E_*K_C',\\
G&=F_{\rm out}W_3,\qquad
W_3=(K_C')^\dagger(E_{C,s}')^\dagger E_{A,s}K_A.
\end{aligned}                                                    \qquad(11)
```

These equations are restricted only after the complete factors have
been shown to preserve the clean code. One common Clifford coordinate
for P_A makes the restricted endpoints actual logical Cliffords.
Each silent list remains independent after restriction, witnessed by
an actual anticommuting active clean axis. The middle list in W_3 is
empty; retaining the three-band bound is conservative. No restriction
to a Z eigensector or dependent-list replacement is used.

F_out commutes $`\mathcal M`$. Since the residual G is diagonal,
$`W_3=F_{\rm out}^\dagger G`$ fixes $`\mathcal M`$ pointwise on
the whole logical register. The
[full-space partial-diagonal theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants)
therefore supplies its finite matrix-feature ledger in these same
computational coordinates.

## 7. Reuse the complete matrix ledger with k at most four

The [four-phase source count](OCCUPIED_BOOLEAN_PAIR_CENTER.md#83-count-exact-features-before-restricting-the-clean-bit)
can now be used directly. Keep the physical clean bit as one extra
core coordinate, so there are q=k+1 physical core bits. The three
actual active bands and E_* all have axes in M and commute
$`\mathcal M`$. Each band's core projections commute, so one fixed
q-bit Clifford diagonalizes them. Their exact eigenphase exponents
have three Boolean digits of degrees at most one, two and four.
E_* and its actual inverse share one feature list. The four distinct
phase bands therefore require at most $`12\,2^{k+1}`$ features.
The actual Clifford C_0 fixes $`\mathcal M`$ and contributes at most
$`4(k+1)+3`$ more. Thus F_out uses at most

```math
A_k=12\,2^{k+1}+4(k+1)+3
```

degree-at-most-four features. Multiplying the full physical matrices
and taking the positive clean block adds no features. No individual
active band or unused clean block is assumed Boolean.

For W_3 the retained partial-diagonal constants are

```math
\begin{gathered}
q_k=2k+4,\qquad N_k=2^{2k^2+5k+10},\qquad
B_k=9\,2^k+2k+3,\\
D_k=48q_kN_k+4B_k+k.
\end{gathered}
```

Its full matrix uses at most B_k degree-at-most-four features and
$`3N_k`$ degree-at-most-$`16q_k`$ features. Use Booleanity only
after multiplying (11). The resulting sign is a Boolean function of
these features, the A_k insertion features and at most k core input
bits. Extending that finite sign function off its attained feature
vectors and composing its algebraic normal form gives

```math
\deg g\le D_k+4A_k\le D_4+4A_4
 =576\cdot2^{62}+2252.                           \qquad(12)
```

The global normalization transports this bound back to the original
input bits: only affine changes affect the residual target's arguments,
and the extracted phase is quadratic. This proves the lambda-one case
and, together with the sharper alternatives, the entire conclusion (1).

The substitution in (12) is **k at most four**, established from the
new quotient (9). It is not a substitution into
$`K=3s+2+t`$, and it does not assert a bound on t. The native
independence and full-space hypotheses needed for W_3 were proved
before invoking its feature theorem.

## 8. Complete contract, inverse, rounding and the remaining question

The original implementation and its actual reversed inverse retain
their original T counts, at most three native T bands, the sole clean
wire and all original dirty allocation. Every orientation correction,
C_0, coordinate word and inverse is an actual Clifford operation with
its literal phase retained; emitting the fixed coordinate work costs
$`O((Q+1)^2)`$ elementary Cliffords on Q logical plus one clean wire.
The full shared-flip dressings cost $`O(Q+1)`$ and allocate no helper.
The sector-free feature description supplies no new synthesis upper
bound and declares no product-axis phase or basis change free.

All identities are exact on their stated complete physical or logical
spaces. The unused clean sector is retained until a justified complete
restriction; arbitrary logical, dirty and reference inputs and the
target's entire return contract are preserved. No accumulated error or
return allowance is omitted, because the reductions add zero error.
The actual inverse reverses every orientation and coordinate word,
with the same resource majorants.

For a Boolean XOR target, output Hadamards change only the actual
endpoint Cliffords. They leave both clean Paulis and all original-list
hypotheses unchanged. The sign polynomial $`y f(x)`$ has degree one
larger than nonzero f, giving the second bound in (1).

Apply the retained
[approximation-to-exactness theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
to the original three-native-band source first. Its conservative floor
$`2^{-530690}`$ then excludes higher-degree targets in this class,
allowing one common scalar. It is not an arbitrary fixed-tolerance
obstruction and is not inferred from a sector encoder or a dependent
phase representation.

The no-singleton premise is used precisely in Section 3. With additional
original middle singletons, a product from the pair-coupling kernel
need not commute those axes. The Boolean endpoint lemma then does not
automatically export its commutation to the opposite silent list.
Finding a sufficiently large actual shared center or a different native
matching mechanism in that case remains the next question. General
original pairs and unrestricted depth remain open, and the general
frame resource frontier is unchanged.

The [positive-singleton extension](ORIGINAL_PAIR_SINGLETON_COUPLING.md)
retains all original positive middle axes at bounded actual singleton-to-silent
coupling rank h. Its enlarged transported center has quotient dimension
at most $`2h+4`$ and recovers this theorem's constant when h is zero.
Unrestricted h and general pair differences remain open.
