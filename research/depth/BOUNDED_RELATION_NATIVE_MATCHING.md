# Bounded relations in three commuting phase bands

[Pauli-boundary theorem](PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md) ·
[Partial-diagonal features](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md) ·
[Native axis matching](ONE_BAND_AXIS_MATCHING.md)

This theorem bounds Boolean degree for an exact full-space three-band
expression when each commuting Pauli-axis list has bounded binary
relation dimension. A dependent list is an algebraic phase factor;
it is not being charged as one independent-axis native physical band.
Section 8 applies the theorem to an actual one-clean native circuit
with no paired middle orbit and bounded endpoint clean-active ranks.
The argument retains the occupied scalar and every actual inverse.

## 1. Statement and deliberately loose constants

Let
```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A
```
be a literal full-space operator, projectively a Boolean diagonal sign.
The endpoint K operators are actual Cliffords. Each E is a product of
literal T or T-dagger phases about commuting signed Hermitian Paulis.
For its binary axis matrix A_j, let
$`r_j=\dim\ker A_j`$, for j=A,B,C.

An identity axis contributes its literal scalar. Repeated binary axes
combine to an even Clifford power or one remaining odd T power.
Absorb the resulting commuting Clifford into a neighboring endpoint
and conjugate the neighboring list, retaining the scalar. These operations
do not increase any relation bound. Hence the proof uses distinct
nonzero axes.

Put
```math
\begin{gathered}
K=32r_B^2,\qquad h=4r_B+7K,\\
\rho=r_B+h,\qquad w=h+4\rho^2,\\
k_*=r_A+w+r_C,\\
D(r_A,r_B,r_C)=48\,2^{k_*}+16k_*+84.
\end{gathered}                                      \qquad(1)
```

**Bounded-relation degree theorem.** If
$`V=e^{i\theta}\mathrm{diag}((-1)^{g(x)})`$
on the entire logical space, then
```math
\deg_{\mathbb F_2}g\le D(r_A,r_B,r_C).               \qquad(2)
```

The bound is independent of the list lengths and total logical width.
It is intentionally much weaker than the retained bound 2048 when
all three lists have at most one relation. The new conclusion is
finiteness for every fixed triple of relation dimensions.

## 2. A rank profile recovers all but boundedly many axes jointly

Let A be a distinct nonzero commuting axis list of relation dimension r.
For a Pauli probe v, let I be its active axis subset and put
$`R_I=\ker A_I`$. The displacement of the Hermitian Clifford
$`EPE^\dagger`$ is
```math
N_v=A_IA_I^{\mathsf T}\Omega,
```
where Omega is the nonsingular ambient symplectic matrix. Signs and
T-dagger orientations do not change this binary matrix. Since
$`\mathrm{im}A_I^{\mathsf T}=R_I^\perp`$ in the ordinary
coordinate dot product,
```math
\begin{aligned}
w(v):=\mathrm{rank}N_v
 &=|I|-\dim R_I-\dim(R_I\cap R_I^\perp),\\
w(v)&=\mathrm{wt}_A(v)-\delta(v),
\qquad 0\le\delta(v)\le2r .
\end{aligned}                                      \qquad(3)
```

If J is Clifford and $`EJE^\dagger`$ is Clifford, its binary action
F preserves w. Indeed the exact equality
$`EJ=(EJE^\dagger)E`$ identifies Pauli images up to Clifford
conjugation, which preserves their Pauli coefficient $`\ell_1`$ norm.
Clifford flatness gives that norm as $`2^{w(v)/2}`$.
Thus
```math
|\mathrm{wt}_A(Fv)-\mathrm{wt}_A(v)|\le2r.
```

The characters $`(-1)^{\langle a,v\rangle}`$ for distinct labels
a are orthogonal on the full Pauli-label space. Expanding the active
weight as half the constant minus half the character sum therefore gives
```math
\mathbb E_v
 [\mathrm{wt}_A(Fv)-\mathrm{wt}_A(v)]^2
 =\frac{|A\mathbin{\triangle}F(A)|}{4}.
```
The distinction between F and its inverse in the character pullback
does not change this symmetric-difference cardinality. Consequently
```math
|A\mathbin{\triangle}F(A)|\le16r^2.                 \qquad(4)
```

For a whole finite group G preserving w, define the proxy set B by
strict-majority orbit membership:
```math
B=\{b:\Pr_{F\in G}[b\in F(A)]>1/2\}.               \qquad(5)
```
It is G-invariant. Each label in $`A\triangle B`$ contributes at
least one half to the average mismatch, so
```math
|A\mathbin{\triangle}B|\le32r^2.                    \qquad(6)
```
Any two labels of B occur together in some F(A): their two membership
events each have probability greater than one half. They therefore
commute. B is an isotropic set, and its relation dimension is at most
$`r+|A\triangle B|`$.

This is a joint statement about the entire group. An arbitrary profile
symmetry need not be an actual Clifford passage through E. Only the
forward implication used above is asserted.

## 3. A bounded invariant exceptional middle algebra

For each computational input Z label z define
```math
R_z=K_AZ(z)K_A^\dagger,\qquad
J_z=E_AR_zE_A^\dagger.
```
The full-space diagonal target implies
```math
E_BJ_zE_B^\dagger=E_CK_CZ(z)K_C^\dagger E_C^\dagger.
                                                        \qquad(7)
```
Both sides are Hermitian Cliffords. The exact J family is a commuting
representation, and its binary actions have the form
```math
F_z=I+N_z,\qquad
N_zv=\sum_i\lambda_i(z)\langle a_i,v\rangle a_i,
\qquad N_zN_u=0.                                    \qquad(8)
```
Peel globally silent A axes to the input as a diagonal C3 factor.
Every remaining activation character lambda_i is nonzero.

Apply Section 2 to the middle list and this finite group. Its proxy B
has relation dimension at most $`r_B+K`$.
An orbit of a proxy label b is an affine set $`b+W_b`$, where
$`W_b=\{N_zb:z\}`$. If its size is $`2^d`$, then b is outside
W_b unless d=0: every N_z kills W_b. Thus the orbit span has
dimension d+1 and its internal relation dimension is
$`2^d-d-1`$. For d at least two,
$`2^d\le4(2^d-d-1)`$.
Relations supported on distinct orbits are independent.
All proxy orbits of size at least four therefore contain at most
$`4(r_B+K)`$ labels altogether.

Remove those large orbits and every singleton or pair orbit meeting
$`B\setminus A_B`$, where A_B is the actual middle list.
The remaining proxy labels are actual axes and form an invariant
ordinary list, with only singleton or pair orbits.
All leftover actual axes form a bad list of size at most
```math
|A_B\setminus B|+4(r_B+K)+2|B\setminus A_B|
 \le h.                                             \qquad(9)
```

Write the actual middle phase as $`E_{\rm bad}E_{\rm ord}`$.
The ordinary signed axis permutation implies
```math
D_{{\rm ord},z}
 =E_{\rm ord}J_zE_{\rm ord}^\dagger J_z^\dagger
```
is a commuting Clifford: each native exponent difference is even.
It commutes E_bad. Since the total passage (7) is Clifford,
$`E_{\rm bad}J_zE_{\rm bad}^\dagger`$ is Clifford too.

Let W_bad be the span of the bad axes and define the intrinsic Pauli
kernel
```math
\mathcal K=
 \{p:E_{\rm bad}P E_{\rm bad}^\dagger
          \text{ is a Pauli up to its literal phase}\}.
```
This is a linear label space, contains $`W_{\rm bad}^\perp`$,
and is invariant under every F_z. The last assertion follows by
conjugating a Pauli image with the actual Clifford
$`E_{\rm bad}J_zE_{\rm bad}^\dagger`$.

There is a commuting-axis Clifford C_bad in the algebra of W_bad
matching E_bad's exact Pauli action on all of K. In fixed diagonal
coordinates for W_bad, the translations in K have affine Pauli
derivatives. Their compatible signed derivatives integrate into an
S/CZ phase polynomial. This is the same
[diagonal Clifford extraction](PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md#4-match-the-endpoint-cliffords-with-exact-pauli-signs)
used in the retained proof; its possible S phases and common scalar
are not discarded. Matching the derivatives leaves
```math
E_{\rm eff}=C_{\rm bad}^\dagger E_{\rm bad}
 \in\mathrm{Alg}(W),\qquad
W=\mathcal K^\perp\subseteq W_{\rm bad}.              \qquad(10)
```
W is isotropic, invariant under the complete J family, and has
dimension at most h. E_eff remains a commuting C3 phase.

Absorb C_bad into the C encoder, conjugating its axis list by an actual
Clifford. This preserves its relation dimension. The middle expression
now consists of the ordinary axes and the bounded invariant exceptional
phase E_eff. No assertion about the possibly large orbits of the
discarded raw axes is needed.

## 4. All but boundedly many ordinary pairs peel exactly

Each ordinary pair has a nonzero linear swap character lambda and a
common linear sign character chi. The exact commuting involution law
proves these signed-character statements. Its Hermitian axis product
R has binary label r, the difference of the paired labels.
R is fixed literally by every J_z. Since
$`r=N_zb`$ for some z, r lies in the A span and E_A commutes R.
Thus every R_z commutes R: r is an input computational parity in the
A coordinates. This remains true with dependent A axes.

Normalize the ordinary native orientations by commuting Clifford
corrections absorbed into the C encoder, retaining all scalar phases.
The exceptional defect is Clifford and lies in Alg(W), so it vanishes
on the symplectic quotient $`W^\perp/W`$.
All ordinary pair differences commute W because the original middle
axes commute. On that quotient, the only potentially quadratic contribution to
the displacement identity (7) is
```math
\sum_j\chi_j(z)\lambda_j(z)\,
              \bar r_j\langle\bar r_j,\cdot\rangle. \qquad(11)
```
The endpoint displacement maps and singleton contributions are linear
in z.

The relation dimension of the projected difference list is at most
$`r_B+\dim W\le\rho`$: disjoint pair incidence vectors are
independent before mapping to the original middle list, and quotienting
by W adds at most dim W relations.
Consider the tensor map $`e_j\mapsto\bar r_j\bar r_j^{\mathsf T}`$.
Its kernel has dimension at most rho, since tensors of an independent
label basis are independent. If a kernel vector has support size t,
let R be the matrix of those t labels, and let
$`\sigma=t-\mathrm{rank}R\le\rho`$.
Sylvester's inequality gives
```math
0=\mathrm{rank}(RR^{\mathsf T})
 \ge t-2\sigma,
```
so t is at most $`2\rho`$.
The union of supports of a kernel basis therefore has size at most
$`2\rho^2`$. This includes zero projected differences.

Outside that exceptional pair set, a dual functional isolates each
tensor in (11). It forces chi lambda to be linear, hence
$`\chi=0`$ or $`\chi=\lambda`$.
Orient these pairs equally or oppositely, respectively. Their entire
literal phase commutes every J_z and peels as a diagonal factor.
The identities $`T_{-P}=\omega T_P^\dagger`$ retain the signed
native scalars; in the opposite-orientation pairing the two scalar
adjustments cancel exactly.

Adjoin the axes of the at-most-$`2\rho^2`$ unpeeled pairs to W.
The resulting invariant isotropic middle space W' has
```math
\dim W'\le h+4\rho^2=w.                              \qquad(12)
```
The exceptional middle phase is still commuting C3, now on W'.

For completeness, the peeled pair exponents have the retained forms
$`u+2(1-u)q`$ or $`-u+2uq`$ modulo eight, where u is affine
and q is quadratic Boolean. The conjugated pair axis is a Hermitian
Clifford and is diagonal on the relevant Pauli parity sector; integrating
its Pauli derivatives gives the quadratic q on that sector. These are
full logical-space diagonal factors, not extensions from a non-Pauli
initialized code. Their ordinary integer degree is at most five.

## 5. Relations also bound singleton coupling and hidden C directions

We use one elementary product-space lemma. Suppose binary subspaces
C,U,R of a coordinate space obey
$`C*U\subseteq R`$ under coordinatewise multiplication, and every
coordinate occurs in the support of U. Then
```math
\dim C\le\dim R.                                    \qquad(13)
```
Pass to any infinite extension field of characteristic two. A generic
u in the extended U has no zero coordinate: it avoids finitely many
proper coordinate hyperplanes. Multiplication by u is invertible, and
sends the extended C into the extended R. Extension preserves all
original binary ranks, proving (13).

For an ordinary singleton b, (8) gives
```math
(\langle a_i,b\rangle)_i*(\lambda_i(z))_i
       \in\ker A_A .
```
The activation subspace has full coordinate support after silent A
peeling. Applying (13) to the span of all singleton pairing vectors
shows their cross-pairing rank with A is at most r_A.
Both lists are internally isotropic, so their combined span S_0 has
symplectic rank at most $`2r_A`$.
Adjoin W' to obtain S. Then
```math
\mathrm{rank}_{\rm symp}S\le2(r_A+w).          \qquad(14)
```

The residual source U, consisting of A, the singleton middle band,
and the exceptional middle C3 phase, belongs to Alg(S).
Its comparison image $`UR_zU^\dagger`$ is Clifford.
The Clifford $`(UR_zU^\dagger)R_z^\dagger`$ belongs to Alg(S),
so it fixes every Pauli in $`S^\perp`$ literally. Its symplectic
displacement is therefore supported in S.

Let c_i be the globally active C axes and let U_C be their activation
subspace. For every $`v\in S^\perp`$,
```math
(\langle c_i,v\rangle)_i*U_C\subseteq\ker A_C .
```
Apply (13) again. The rank of the map
$`v\mapsto(\langle c_i,v\rangle)_i`$ on $`S^\perp`$ is at
most r_C. By symplectic duality this is exactly
$`\dim((L_C+S)/S)`$. Thus at most r_C additional C directions
are needed. Peel silent C axes to the output and put
$`\widehat S=S+L_C`$. Then
```math
\mathrm{rank}_{\rm symp}\widehat S
 \le2(r_A+w+r_C)=2k_* .                             \qquad(15)
```

All four remaining commuting C3 factors belong to this one Pauli
algebra: active A, singleton B, exceptional middle phase, and active C.
They need not commute with one another. The number of quantum
coordinates in their common algebra is at most k_*; its central
dimension may still grow with width.

## 6. Restore the phase using finitely many low-degree features

Write the diagonal residual as
$`V_0=K_C^\dagger W_0K_A`$, where W_0 is the ordered product
of those four factors, with the C factor adjointed.
It need not be real, because the peeled phases have not yet been restored.
On $`K_A^\dagger\widehat S^\perp K_A`$, V_0 has a literal Pauli
action. Extract a diagonal Clifford C with the same signed action and put
```math
J=K_CCK_A^\dagger,\qquad
V_1=C^\dagger V_0=K_A^\dagger J^\dagger W_0K_A.
```
Then J belongs to Alg($`\widehat S`$), and V_1 is diagonal with
Pauli support in
$`K_A^\dagger(\widehat S\cap K_A\langle Z\rangle K_A^\dagger)K_A`$.

Use canonical central Z coordinates and k quantum coordinates for
Alg($`\widehat S`$), with $`k\le k_*`$.
The core projections of the final isotropic support commute. A fixed
core Clifford makes them Z-type; each eigenvalue bit has at most one
central parity offset. This does not assume that the canonical Clifford
preserves all input computational Paulis.

Those target parity eigenbits are affine functions of the original
computational input bits. More explicitly, prescribe their central
eigenbits, set complementary central coordinates to zero, set each
selected core coordinate to its target eigenbit plus its central parity
offset, and set other core coordinates to zero. Pauli-support uniqueness
makes this a section evaluating exactly the same diagonal phase.
Its substitutions are affine in the original computational input bits.

Each of the four remaining C3 factors has commuting axes in the
canonical algebra. A fixed core Clifford diagonalizes its core
projections. For each fixed core string, its central exponent modulo
eight is weighted cubic: linear terms have arbitrary coefficients,
quadratic terms even coefficients, and cubic terms coefficients
divisible by four. All three binary digits have Boolean degree at most
four by the binary carry formulas. The four factors therefore require
at most $`12\,2^k`$ features of degree at most four.

The blocks of J have the exact form
```math
J(z)=\gamma i^{a(z)}X^{u(z)}Z^{v(z)}J_* ,
```
using ordered X then Z representatives. The labels u,v are affine;
after removing their controlled Pauli Cliffords, the remaining scalar
is a central diagonal Clifford phase. Its digits have degree at most
four, conservatively including a change to Hermitian representatives.
Together with the central offsets and core eigenbits selecting the
entry, these need at most $`4k+3`$ further features of degree at
most four. The fixed common scalar gamma is retained.

All peeled pair exponents have ordinary degree at most five modulo
eight. The silent-band exponents are weighted cubic, so their sum with
the pair exponents has the same ordinary-degree upper bound five.
Its three binary digits have degree at most twenty: the largest carry
has degree four in the individual monomial bits. The three digits of
the extracted diagonal Clifford C have degree at most four.

The original Boolean target guarantees that, after every exact scalar
and peeled phase is restored, its sign is a Boolean function of these
features on their attainable assignments. Extend that function
arbitrarily elsewhere. Composition bounds its degree by the sum of
feature degrees:
```math
\begin{aligned}
\deg g
 &\le4(12\,2^k+4k+3)+3\cdot20+3\cdot4\\
 &=48\,2^k+16k+84
 \le D(r_A,r_B,r_C).
\end{aligned}                                      \qquad(16)
```
No entry of an arbitrary intermediate core matrix was assumed Boolean
or an eighth root. The exact target and the restored diagonal phases
supply that final phase condition.

## 7. What this theorem does and does not compile

The proof is full-space and permits arbitrary endpoint Clifford mixing.
It does not convert an arbitrary dependent phase list into one physical
native band. Identity axes, repeated axes, diagonal Clifford extraction,
and the effective exceptional C3 phase are exact algebraic operations;
their actual synthesis costs must be charged in any implementation.

Likewise, a bounded-support operation inserted between otherwise native
bands cannot simply be deleted before applying a partial-diagonal
theorem. Conjugating it through a silent native band can introduce
arbitrarily many independent axis constraints. The present argument
uses commuting phase lists and their bounded relation dimensions to
control that operation. It makes no such claim for a general bounded
noncommuting insertion.

Taking adjoints throughout gives the actual inverse, with conjugated
scalars and reversed T orientations. No separately selected occupied
extension or inverse is substituted.

## 8. One clean qubit, no middle pairs, and bounded active ranks

Consider an actual one-clean native three-band circuit with complete
clean and dirty return and occupied Boolean sign action.
Let
```math
\begin{gathered}
P_A=K_AZ_cK_A^\dagger,\qquad
P_C=K_CZ_cK_C^\dagger,\\
J=E_AP_AE_A^\dagger,\qquad
S=E_CP_CE_C^\dagger=E_BJE_B^\dagger .
\end{gathered}
```
Suppose both clean images are non-Pauli, their intrinsic active ranks
are k_A,k_C, and native matching for this clean probe has no paired
middle orbit. Set $`d=k_A+k_C`$ and
```math
s(d)=d+3\binom d2+7\binom d3.
```
Then the occupied sign and XOR query obey
```math
\deg g\le D(1,1+s(d),1),\qquad
\deg f\le D(1,1+s(d),1)-1.                         \qquad(17)
```

Here is the actual reduction. Split A and C into active and silent
subbands for their clean Paulis. Since every middle orbit is a singleton,
native independence forces each B axis to commute every active A and C
axis. Let E_{B,-} be the axes negated by J and E_{B,+} the others.
The negative signs are exactly anticommutation with P_A. The clean
displacement identity is
```math
N_C=N_A+N_{B,-}.
```
Native independence implies
$`L_{B,-}\subseteq L_A+L_C`$ and
$`L_C\subseteq L_A+L_{B,-}`$.
Consequently
```math
L_0=L_A+L_{B,-}=L_A+L_C
```
is isotropic and has dimension at most d. Every positive singleton
commutes L_0 and P_A.

The literal middle active word
```math
F_0=E_{C,\rm active}^\dagger
        E_{B,-}E_{A,\rm active}
```
is a commuting C3 phase in Alg(L_0), and
$`F_0P_AF_0^\dagger=P_C`$.
Extract a commuting-axis Clifford C_0 in Alg(L_0) with the same
signed action on P_A. This is the one-translation instance of the
diagonal Clifford extraction above. The residual
$`G_0=C_0^\dagger F_0`$ remains commuting C3 and commutes P_A.
Pauli-support uniqueness puts it in
$`\mathrm{Alg}(L_0\cap P_A^\perp)`$.

Absorb C_0 into the C endpoint and conjugate its silent band.
The complete actual word is now
```math
V=(K_C')^\dagger(E_{C,\rm silent}')^\dagger
             G_0E_{B,+}E_{A,\rm silent}K_A ,
```
where all displayed factors preserve the same Pauli P_A and both
endpoints map Z_c to P_A. This identity holds on both clean sectors,
not just the occupied code.

Restrict to the positive P_A code using fixed Clifford coordinates.
The three silent/positive native sublists project to full logical
commuting lists with at most one relation each. The residual G_0
projects to a commuting C3 phase on at most d coordinates and commutes
the middle projected list. Modulo a retained diagonal Clifford, a
weighted cubic phase on d coordinates uses at most s(d) odd parity-T
factors: one per linear monomial, three for a quadratic monomial
using $`2xy=x+y-(x\oplus y)`$, and seven for a cubic monomial
using the usual parity expansion of $`4xyz`$.
The absorbed Clifford commutes the middle factors, and its scalar is
retained. Adding at most s(d) axes to a relation-at-most-one list
raises its relation dimension by at most s(d). Identities and duplicates
are handled as in Section 1.

The occupied full logical operator therefore satisfies Sections 1--6
with relation bounds $`(1,1+s(d),1)`$, proving the sign bound.
For an XOR query, conjugating its target by Hadamards gives
$`(-1)^{z f(x)}`$ without changing either clean image or the middle
pair condition. This proves the second bound in (17).
The [R74 approximation-to-exactness consequence](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
gives the same obstruction below occupied operator error
$`2^{-530690}`$, with an optimized common phase allowed.

The dependent middle phase expression is not claimed to be a single
paid native band after projection. Its justified role is as an input
to the bounded-relation theorem. The later
[finite-insertion theorem](FINITE_INSERTION_THREE_LAYER_FEATURES.md#8-one-clean-qubit-with-bounded-active-ranks-including-middle-pairs)
closes fixed active ranks even with arbitrary middle pairs, retaining
the residual as a bounded Pauli-supported operation rather than one
commuting list. A rank-independent bound remains open. The
[easy-transition family](CLEAN_TRANSITION_PARAMETERS.md)
shows that uniformly compressing the residual to bounded relation
dimension cannot supply that bound. The complete-frame resource
frontier is unchanged.

The [exact dependent-axis fixtures](../../tests/test_bounded_relation_native_matching.py)
check the rank profile with signed native phases, the nontrivial
strict-majority proxy for a relation-four list, and the distinction
between a norm-profile symmetry and an actual Clifford passage.
The width-independent theorem follows from the argument above.
