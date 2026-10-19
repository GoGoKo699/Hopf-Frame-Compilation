# A common center forced by the occupied Boolean phase

[Original physical pair insertions](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md) ·
[Finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md) ·
[No-pair native orbits](GROWING_DISPLACEMENT_NATIVE_ORBITS.md)

The remaining three-layer question concerns **original physical middle
pairs** with growing endpoint clean-active ranks. The retained no-pair
theorem already handles growing ranks and displacement when those pairs
are absent. The fixed-active-rank theorem already permits arbitrary
middle pairs. Neither theorem proves a bound uniform in growing active
ranks with a fixed positive pair count.

This note establishes a constraint from the literal occupied Boolean
condition: large subspaces of the active endpoint spans must centralize
every actual native band. For m original physical pairs, the complete
clean transition has at most 3m+2 Pauli directions modulo that genuine
center, and its clean-preserving residual has at most 3m+1. Its central
dependence has finitely many bounded-degree Boolean features at fixed m.
This is a structural reduction, not a fixed-pair degree theorem. The
silent bands and their actual source dependencies remain in the word.

## 1. Actual source and the endpoint commutation lemma

Use the literal physical expression

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,
\quad P_A=K_AZ_cK_A^\dagger,
\quad P_C=K_CZ_cK_C^\dagger,
\quad J=E_AP_AE_A^\dagger,
\quad S=E_CP_CE_C^\dagger=E_BJE_B^\dagger.          \qquad(1)
```

Each E is an actual independent commuting signed-Pauli T/T-dagger
band. The equality in (1) follows from complete clean return. The
positive clean block, on the entire logical and dirty register, is

```math
G=e^{i\theta}\mathrm{diag}((-1)^{g(x)}).
```

Every phase and endpoint is retained. The unused clean block is not
assumed Boolean. Let L_A,L_C denote the spans of the clean-active
endpoint axes, and let L_B be the span of the entire original middle
list. Throughout Sections 1--5 **both endpoint active sets are
nonempty**. This holds automatically if any original physical middle
pair exists.

**Endpoint commutation lemma.** If a signed physical Pauli Q commutes
with P_A and every original A and B axis, then it also commutes with
P_C and every original C axis. The time-reversed assertion holds as
well. Consequently

```math
V(K_A^\dagger QK_A)V^\dagger=K_C^\dagger QK_C       \qquad(2)
```

is an actual full-physical-space Pauli identity, including its sign.

To prove it, Q first commutes J and S. For the non-Pauli one-band
image S, the independent native displacement formula gives

```math
N_Sq=\sum_{i\in C_{\rm act}}\langle c_i,q\rangle c_i.
```

Thus Q commutes each C-active axis. Its remaining commutation sign
with S is exactly its sign with P_C, so it commutes P_C too. This is
the actual Pauli centralizer calculation in
[clean-code rigidity, Section 1.2](CLEAN_CODE_CLIFFORD_RIGIDITY.md#12-the-one--and-two-axis-cases-with-literal-phases).

The C-silent sublist remains independent on the positive P_C code.
Indeed an actual C-active axis commutes every axis of that sublist
and anticommutes P_C. Hence P_C is not in its physical span; quotienting
by the clean Pauli introduces no relation. Let q_A,q_C be the actual
logical restrictions of the endpoint Paulis in (2), and let C_s be
the projected C-silent native band. Restricting the actual word gives
the full-logical-space identity

```math
Gq_AG^\dagger=C_s^\dagger q_CC_s.                  \qquad(3)
```

No assertion about the unused block was used to obtain (3).

The left side is a Pauli times a diagonal Boolean derivative. The
right side is a Hermitian Clifford. Its displacement image is therefore
computational Z. The independent-native intrinsic-axis argument makes
every C_s axis active on q_C a computational Z axis. The exact
one-band conjugation formula then expresses the Boolean derivative
as a scalar, a computational Z Pauli, and a product of odd S powers
on those active axes. The two Pauli cores have the same X mask.
Squaring the literal Boolean derivative gives identity, so the sum of
the active Z labels is zero. Independence of the projected list forces
the active list to be empty. This is the Boolean-reality argument of
[computational derivatives, Section 7](COMPUTATIONAL_PHASE_DERIVATIVES.md#7-boolean-reality-gives-a-shared-computational-center),
applied only after the legitimate full logical reduction (3).
Q therefore commutes the C-silent axes as well, proving the lemma.
Apply the same proof to V-dagger for time reversal.

## 2. A genuine center at bounded original pair count

The actual native clean comparison matches the middle list into signed
singletons and m disjoint pairs \{b_j,d_j\}. Their binary differences
r_j=b_j+d_j are independent, lie in L_A and L_C, and commute with both
clean Paulis and all original bands. These are the original physical
pairs; no later logical orbit permutation is being counted here.

Define

```math
R_A=L_A\cap p_A^\perp\cap L_B^\perp,
\quad R_C=L_C\cap p_C^\perp\cap L_B^\perp,
\quad R=R_A+R_C.                                   \qquad(4)
```

Every element of R_A commutes the entire A band because the original
A span is isotropic. Section 1 then exports its commutation to the
entire C band and P_C. The time-reversed argument does the same for
R_C. Consequently R is isotropic, commutes with both clean Paulis,
and centralizes **every actual A, B and C axis**. Also
r_j belongs to R_A intersect R_C.

There is a uniform codimension estimate. The image of N_A on L_B
is exactly the span of the m independent r_j: singletons have zero
displacement, and each member of the jth pair has displacement r_j.
Independence of the A-active axes identifies this rank with the
pairing rank between L_A and L_B. Therefore

```math
\mathrm{codim}_{L_A}R_A\le m+1,
\qquad \mathrm{codim}_{L_C}R_C\le m+1.        \qquad(5)
```

The extra one in each bound is the clean-Pauli commutation constraint.
Let

```math
S_0=L_A+L_C+\mathrm{span}\{b_j,d_j:1\le j\le m\}.
```

Since d_j=b_j+r_j with r_j in R, (5) gives

```math
R\subseteq\mathrm{rad}(S_0),\qquad
\dim(S_0/R)\le3m+2.                                \qquad(6)
```

The dimension of R can grow without bound. Equation (6) bounds the
quotient, not the total Pauli support dimension.

## 3. The complete clean transition and its retained logical bands

Use the actual localization and signed Clifford extension of
[finite insertion, Sections 8.1--8.2](FINITE_INSERTION_THREE_LAYER_FEATURES.md#81-localize-the-entire-clean-transition).
Every negative singleton belongs to L_A+L_C. Every positive singleton
commutes S_0 and P_A. Let E_B,act contain all pairs and all negative
singletons, and E_B,+ contain the positive singletons. Then

```math
F_0=E_{C,\rm act}^\dagger E_{B,\rm act}E_{A,\rm act}
 \in\mathrm{Alg}(S_0),\qquad
F_0P_AF_0^\dagger=P_C.
```

The actual signed action on S_0-perp and P_A extends to a Clifford C_0
which fixes S_0-perp pointwise and sends P_A to P_C. Thus

```math
C_0\in\mathrm{Alg}(S_0),\qquad
G_0=C_0^\dagger F_0\in\mathrm{Alg}(S_0\cap p_A^\perp).
                                                               \qquad(7)
```

Both C_0 and G_0 commute R, and G_0 commutes every original positive
middle axis. Since R is contained in p_A-perp and an actual A-active
axis in S_0 anticommutes P_A, (6) improves to

```math
\dim\bigl((S_0\cap p_A^\perp)/R\bigr)\le3m+1.     \qquad(8)
```

With K_C'=C_0-dagger K_C and
E_C,sil'=C_0-dagger E_C,sil C_0, the exact full-space word is

```math
V=(K_C')^\dagger(E_{C,\rm sil}')^\dagger
     G_0 E_{B,+}E_{A,\rm sil}K_A.                  \qquad(9)
```

Every internal factor preserves P_A, and both endpoints encode that
same clean Pauli. All three phase lists in (9) remain independent
after restriction to the positive code, for the following actual
physical reasons:

- An A-active axis commutes the A-silent list and anticommutes P_A.
- The same axis commutes every B-positive axis: binary fixation by J
  and independence of the A-active axes force their pairings to vanish.
- Conjugating an original C-active axis by C_0-dagger gives a Pauli
  commuting the transformed C-silent list and anticommuting P_A.

Each witness excludes P_A from the corresponding span. No clean
quotient relation is introduced. Likewise R projects injectively:
P_A cannot belong to R, since R commutes every A-active axis.
The logical image of (8) can only have smaller dimension.

This supplies three independent actual logical native lists and a
retained G_0 with bounded **noncentral** support at fixed m. It does
not make G_0 a native band or bound its total support.

## 4. Computational center and a finite feature description

For every signed Q in R, (2) gives an actual endpoint Pauli transport.
The Boolean derivative along its computational X mask is affine.
Consequently one actual quadratic Boolean phase matches all those
derivatives simultaneously. Extract that Z/CZ Clifford, then use the
fixed Clifford alignment in
[computational derivatives, Section 7](COMPUTATIONAL_PHASE_DERIVATIVES.md#7-boolean-reality-gives-a-shared-computational-center).
It makes R a shared computational Z subgroup at both logical endpoints.
The residual target remains Boolean diagonal. Hadamards in that
alignment are applied only to factors on which the residual target
is already identity. No separate coordinate circuit is selected in
different central sectors.

The central variables z are now genuine computational parities of the
aligned input register. Transporting a target-degree assertion back
uses its identity on the Hadamard-rotated factors; it does not regard
Hadamards as affine changes of the original computational bits.
The following feature
bound applies to the entire retained clean-transition matrix.
Adjoin P_A to S_0 when choosing fixed canonical coordinates. By (6),

```math
\dim\bigl((S_0+\langle p_A\rangle)/R\bigr)\le3m+3.
```

Choose a complement to R containing p_A, which is possible because
p_A is not in R, and use fixed Pauli coordinates on that complement.
It is enough to use q at most 3m+3 non-R qubits, allowing bounded
additional commuting directions as qubits in this conservative count.
P_A is a fixed core Pauli, so restricting G_0 to its positive block
requires no central-sector-dependent choice of code coordinates.

Each active band has commuting core Pauli projections. Diagonalize
them by one fixed core Clifford for that band. At any fixed core basis
string y its literal eigenphase has the form

```math
\omega^{a_y(z)},\qquad
a_y(z)=\sum_i\epsilon_i\ell_{i,y}(z)\pmod8,\quad
\epsilon_i\in\{1,-1\},                            \qquad(10)
```

where each ell is an affine Boolean parity, including the original
signed-axis constant. The three binary digits of an integer weighted
sum of affine bits modulo eight have Boolean degrees at most 1,2,4.
For example, expanding the sum of its integer coefficients into
weighted bit copies reduces the digits to the elementary-symmetric
carry formulas of orders 1,2,4; constant copies and negative weights
are handled modulo eight. Thus a band matrix depends on at most
3 times 2^q features of degree at most four, independently of its
number of original axes and the dimension of R.

The actual Clifford C_0 fixes R pointwise. In these coordinates it
has a fixed core Clifford action, affine central controls on ordered
core Pauli labels, and a central diagonal-Clifford scalar. The usual
signed Clifford normal form therefore requires at most 4q+3 further
features of degree at most four, retaining its literal scalar. This
is the same matrix feature argument used in
[finite insertion, Section 6](FINITE_INSERTION_THREE_LAYER_FEATURES.md#6-exact-matrix-features-and-boolean-phase-extraction).

It follows that the exact G_0 matrix, and its actual inverse matrix,
are functions of at most

```math
N(m)=9\,2^{3m+3}+12m+15                            \qquad(11)
```

Boolean features of central degree at most four. The bound is
deliberately coarse. It is a coordinate description of the actual
operator, not a proposed free controlled-gate implementation.

## 5. Why the finite-insertion degree theorem still does not apply

A bounded number of bounded-degree central controls does not bound
Pauli coefficient norm distortion. For r independent central pairs,
put

```math
q(z)=\sum_{j=1}^r z_{2j-1}z_{2j}\pmod2,\qquad
F_q=X_t^{q(z)}.
```

This uses one quadratic control feature and one noncentral Pauli
direction over the center. Nevertheless

```math
F_q Z_t F_q^\dagger
   =Z_t\prod_{j=1}^r\mathrm{CZ}_{2j-1,2j}     \qquad(12)
```

has Clifford displacement rank 2r, 4^r nonzero Pauli coefficients,
and Pauli coefficient l1 norm 2^r. Hence the uniform norm estimate
for an insertion of bounded **total** Pauli support cannot be replaced
by the bound (8) or (11). This is a counterexample to that proposed
inference, not a family asserted to arise from a fixed-pair three-band
clean source. The complete original native consistency conditions
have not been discarded in the actual question.

There is a second distinction. Fixing central features can make
independent original axes project to arbitrarily many equal core
axes. It need not leave an independent native band. Replacing the
features by formal independent bits gives a word whose known Boolean
action holds only on the nonlinear graph of the true feature map;
its action off that graph has not been proved Boolean. Neither move
supplies the hypotheses of the full-space finite-insertion theorem.

The first remaining inference is therefore a **global native matching
or endpoint-activation theorem with these actual central controls**.
It must use the retained A-silent, B-positive and C-silent lists and
the complete Boolean identity in (9), rather than deleting G_0,
extending selected fibers, or treating its feature table as arbitrary.
The case m=1 already has residual noncentral dimension at most four,
but no uniform degree bound is concluded here. Growing original pair
count and higher-layer asymptotics remain open as well.

The scalar on a physical pair sector alone does not supply the missing
obstruction. In the retained k divisible by four family, take
$`B=X^{\otimes k},D=Y^{\otimes k},R=BD=Z^{\otimes k}`$.
Add an arbitrary spectator t, an actual
A-silent axis Z_t, and an actual C-silent axis R Z_t. Each enlarged
endpoint list remains independent: its new axis has a nonzero t
component, absent from the active endpoint span. The actual clean
comparison and the original pair count m=1 are unchanged. On R=-1,
the pair factor is omega I, whereas the two extra endpoint factors give

```math
T_{Z_t}T_{-Z_t}^\dagger=\omega^{-1}S_{Z_t}.
```

They cancel that pair-sector scalar exactly. On R=+1 their ratio is
identity. The resulting occupied family still has distinct 1 and i
phases, so it is not a Boolean counterexample. It does show, using
actual independent native source lists, why the pair's odd phase
cannot be isolated from the retained silent endpoint factors.

## 6. Necessary hypotheses, exact checks, and resources

The opposite endpoint must have a nonempty active set in Section 1.
Here is an actual two-wire guard, with wire 0 clean:

```math
K_C=\mathrm{CNOT}_{1\to0},\quad K_A=S_0K_C,
\quad E_A=E_B=I,\quad E_C=T_0T_1.
```

Both clean Paulis are Z_0 Z_1, and
K_C-dagger E_C-dagger K_A is literally identity on every clean-zero
input. But Q=X_0 X_1 commutes P_A and the empty A/B lists and
anticommutes both C axes. Here the C-active list is empty; the
C-silent axes acquire one relation on the clean code. The already
closed Pauli-interior branch must be treated separately.

The Boolean premise is also necessary. With K_A=K_C=H_0,
E_A=T_0, E_B=I and E_C=T_0 T_1, both interior clean images are
non-Pauli and the clean wire returns, but the occupied operator is
T_1-dagger. The Pauli Q=X_1 commutes P_A and all A/B axes while
anticommuting the silent C axis Z_1.

The [exact fixtures](../../tests/test_occupied_boolean_pair_center.py)
check these two guards, an actual Boolean source with non-Pauli
endpoint images and nontrivial affine transported derivative, and
the native implementation and complete inverse of the quadratic
central-control example. They also check the actual one-pair silent-axis
phase compensation at k=4, without calling its occupied action Boolean.
Full-space identities are checked on all
physical columns, and occupied claims on every initialized-isometry
column, including arbitrary dirty and output inputs. Reference
coherence follows from the operator identities; no sampled-input
criterion is substituted. Existing original-pair fixtures retain
the non-Boolean support-growth family and are not relabeled as
Boolean examples.

The reduction preserves the original native word, its literal phases
and occupied action. Fixed coordinate Cliffords, the quadratic
correction and C_0 cost O(Q^2) elementary Cliffords each if physically
compiled on Q wires. Their actual inverse words carry the same
separate charges. No new initialized helper or free phase operation
is introduced. Equations (8) and (11) provide neither an upper
compiler nor an unrestricted lower degree theorem. The complete-frame
resource frontier is unchanged.
