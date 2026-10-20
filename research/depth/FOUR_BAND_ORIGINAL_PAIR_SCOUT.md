# Four-band scout: actual pair transport and an added-band coupling certificate

[Three-band original-pair collapse](BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md) ·
[Actual native matching](ONE_BAND_AXIS_MATCHING.md) ·
[Complete two-band quadratic source](QUADRATIC_NATIVE_QUERY.md)

This is a bounded continuation of the original physical-pair scout.
It does not claim a four-layer Boolean-degree theorem. One additional
actual native band removes the automatic intermediate Clifford matching
used at three layers. Under explicitly stated physical-center hypotheses,
complete Booleanity still gives a constraint on the additional band's
couplings. A literal four-band identity saturates that constraint and
shows why the original-representation pair collapse itself cannot extend
unchanged to four bands. A second, nontrivial Boolean example has disjoint
original middle pairs, but its complete endpoint removal still gives two
native bands. The final sections prove cubic Boolean degree throughout
the stated fixed-endpoint, one-pair orbit class. They do not turn its
commuting phase lists into a free native layer or cover arbitrary
four-band return.

## 1. The first inference that fails with an additional native band

Use the actual physical word

```math
V=K_C^\dagger E_C^\dagger E_D E_B E_AK_A,
\qquad
P_A=K_AZ_cK_A^\dagger,\quad P_C=K_CZ_cK_C^\dagger,
\quad J=E_AP_AE_A^\dagger.                         \qquad(1)
```

Each E is an independent commuting signed Pauli T/T-dagger list, with
its literal scalar phases. There is one initialized clean qubit c;
every other logical and dirty input is arbitrary. The occupied action
is a complete Boolean diagonal unitary G up to one retained common
phase, with exact clean and dirty return and arbitrary references.
The unused clean block is unrestricted.

Clean return gives the physical identity

```math
E_D(E_BJE_B^\dagger)E_D^\dagger
 =E_CP_CE_C^\dagger.                              \qquad(2)
```

Both J and the right side are Hermitian Cliffords. Equation (2) alone
does not make $`E_BJE_B^\dagger`$ Clifford. Thus it does not permit
the three-band native matching step at B. In particular, a later
logical computational-probe pair orbit is not an original physical
B-pair under J.

The subclass below therefore assumes the actual signed permutation
property of the original B axes under J. It is an additional premise,
not a consequence silently borrowed from the three-band theorem.

## 2. Explicit physical hypotheses and literal pair repair

Assume J sends each original B axis to a signed member of that same
independent list. Since $`J^2=I`$, these are actual singletons or pairs.
For each of the m genuine pairs write

```math
JB_jJ=s_jD_j,\qquad JD_jJ=s_jB_j,\qquad
\Gamma_j=s_jB_jD_j,
\qquad r_j=b_j+d_j.                                \qquad(3)
```

Here D_j is the partner axis inside the B list; it is distinct from
the name E_D of the additional native band. The actual pair differences
are nonzero and independent. Apply the retained literal orientation
repair

```math
E_B=C_{\rm rep}F E_-E_+,
\qquad F=\prod_jF_j,
\qquad
F_j=T_{B_j}^{\epsilon_j}T_{D_j}^{s_j\epsilon_j}.     \qquad(4)
```

The negative and positive singleton sublists keep their actual axes
and orientations. Absorb the actual Clifford $`C_{\rm rep}`$ to
the left by conjugating E_D and E_C by it and replacing K_C by
$`C_{\rm rep}^\dagger K_C`$. This is an exact full-word identity;
it preserves the sizes and independence of the three affected native
lists. From now on all D and C notation refers to that repaired frame.
The complete source is

```math
V=K_C^\dagger E_C^\dagger E_D F E_-E_+E_AK_A,
\qquad
W=E_-JE_-^\dagger,
\quad S=E_CP_CE_C^\dagger=E_DWE_D^\dagger.          \qquad(5)
```

Let L_A and L_C be the spans of the actual clean-active endpoint axes
in this frame. The additional physical-center premises are

```math
R_B:=\mathrm{span}\{r_j\}\subseteq L_A\cap L_C,
\qquad [E_D\text{ axis},\Gamma_j]=0
\quad\text{for every actual D axis and every }j.    \qquad(6)
```

The input inclusion follows already from the genuine J-pair relation;
the output inclusion and additional-band fixation are imposed. They
must be checked on the actual physical lists of a proposed source.
They are not statements about a selected occupied extension.

Each pair parity then commutes every band in (5), both clean Paulis,
and every singleton. For the output clean Pauli, use fixation of the
parity by E_D and commutation with W to get commutation with S; its
membership in L_C makes it commute with E_C, giving commutation with
P_C. Its exact endpoint transport is a Pauli intertwining through G.
The retained shared-center normalization therefore applies: one real
quadratic aligning phase and fixed Clifford coordinates make all the
pair parities the same computational Z parities at both endpoints.
This normalization uses only those literal Pauli intertwinings and
complete Booleanity, not a count of native bands. All lists and signed
relations are carried through the same fixed changes.

The pair formula remains

```math
F_j=\omega^{-\epsilon_j\eta_j}
 \left(\Pi_{j,+}S_{B_j}^{\epsilon_j}
       +\Pi_{j,-}\omega^{\epsilon_j}I\right),
\qquad
\eta_j=(1-s_j)/2,
\quad \Pi_{j,\pm}=(I\pm\Gamma_j)/2.                \qquad(7)
```

No scalar branch is deleted from the original source. Every assignment
of pair signs is nonempty on the positive clean code: an actual active
axis commutes with all pair differences and anticommutes with the clean
Pauli, so the clean label is independent of the independent pair labels.

## 3. The actual inherited scalar-sector transport

Take a signed Hermitian Pauli P_t with

```math
t\in L_A\cap p_A^\perp.                           \qquad(8)
```

It commutes E_A, every singleton, J, W and all pair parities. The
singleton assertion uses the actual fixed label under J and independence
of the A-active axes. Let t_A be the legitimate clean-code restriction
of $`K_A^\dagger P_tK_A`$.

Transport P_t through the entire physical word before restricting.
On the joint scalar pair sector, the result is exactly

```math
(Gt_AG^\dagger)\Pi_{\boldsymbol -}
 =R_{t,0}\Pi_{\boldsymbol -},
\qquad
R_t=K_C^\dagger E_C^\dagger E_D P_tE_D^\dagger E_CK_C.
                                                               \qquad(9)
```

Here $`R_{t,0}`$ is the actual positive-clean restriction of R_t.
It exists: P_t commutes W, so $`E_DP_tE_D^\dagger`$ commutes S;
conjugating through E_C and K_C makes R_t commute with the input
clean Pauli. It also preserves every pair sector. Equation (9) holds
on every remaining logical column, including arbitrary dirty inputs
and reference states.

The right side is a physical two-band Pauli transport followed by
its legitimate code restriction. The two individual bands need not
preserve an intervening Pauli code. Thus (9) must not be relabeled as
two independent native logical bands without an additional proof.
Nor does P_t yet have to commute with P_C; a logical Pauli t_C has
not been introduced at this stage.

## 4. Directions untouched by the added band still obey Boolean rigidity

Now additionally require

```math
t\in L_D^\perp,
\qquad L_D=\mathrm{span}\{\text{actual D-axis labels}\}.          \qquad(10)
```

Then E_D commutes P_t literally. Equations (5) and (8) give
$`[P_t,S]=0`$. The independent C-active displacement formula implies
that P_t commutes every C-active axis. Removing their rotations from S
then gives $`[P_t,P_C]=0`$. Thus t_C, the clean restriction of
$`K_C^\dagger P_tK_C`$, is an actual logical Pauli.

The output-silent list remains independent after the clean and pair
quotients. Indeed

```math
\mathrm{span}(C_{\rm silent})
 \cap\bigl(\langle p_C\rangle+R_B\bigr)=\{0\}.      \qquad(11)
```

Pairing a proposed relation with an actual C-active axis excludes its
p_C coefficient. Then $`R_B\subseteq L_C`$ and full C-list
independence exclude the remaining silent-versus-active relation.
This is the same explicitly safe quotient as in the three-band proof.

Equation (9) is now the one-band output-silent comparison. In the one
fixed computational coordinate system, its joint minus sector is a
complete Boolean diagonal operator on the whole remaining register.
The one-band Boolean-reality lemma therefore forces every projected
output-silent axis to commute t_C. Individual physical commutation
signs survive the nonzero invariant code, so every actual C-silent
axis commutes P_t on the full physical space. Hence both E_D and E_C
commute P_t, and

```math
(Gt_AG^\dagger)\Pi_{\boldsymbol -}
 =t_C\Pi_{\boldsymbol -}.                          \qquad(12)
```

This makes t_A and t_C have the same full computational X mask u.
For the complete Boolean target, its derivative
$`(-1)^{g(x)+g(x+u)}`$ is invariant under translation by u.
Consequently $`Gt_AG^\dagger`$ has one global commutation sign
with t_C; (12) makes that sign positive.

If $`\langle b_j,t\rangle=1`$, choose the nonempty sector where
only pair j has plus sign. Formula (7) sends P_t to a Hermitian
Pauli anticommuting with P_t, up to the retained quarter-turn sign.
Conjugation by E_D and E_C preserves that anticommutation because
both commute P_t. No assertion that the conjugated B_j is Pauli is
needed. The corresponding block of $`Gt_AG^\dagger`$ therefore
anticommutes with t_C, contradicting its global positive sign. Thus

```math
\langle b_j,t\rangle=0
\quad\text{for every }j
\text{ and every }t\in L_A\cap p_A^\perp\cap L_D^\perp.
                                                               \qquad(13)
```

This is the portion of paired-sector Boolean rigidity that survives
the extra native band under the explicit physical premises (6).

## 5. A radical certificate for the actual original pairs

Use the actual independent A-active axes as the coordinate basis of
L_A. Let k be their number, let e be the all-ones vector, and let M
be the matrix whose row for an actual D axis d is

```math
M_{d,i}=\langle d,a_i\rangle,
\qquad
\mathcal W=\mathrm{span}\{e,\text{rows of }M\}
 \subseteq\mathbb F_2^k.                           \qquad(14)
```

The dot product below is the ordinary coordinate dot product, not
the physical Pauli symplectic form. For pair j define

```math
v_{j,i}=\langle b_j,a_i\rangle,
\qquad r_j=\sum_i v_{j,i}a_i.                     \qquad(15)
```

The second equality is the original one-band displacement formula
$`r_j=(\Phi_J-I)b_j`$. Equation (13) says that v_j annihilates
the simultaneous kernel of e and M, hence $`v_j\in\mathcal W`$.
The parity r_j is clean-even and fixed by every D axis, so

```math
e\mathbin\cdot v_j=0,\qquad Mv_j=0.
```

Therefore every actual original pair difference satisfies

```math
v_j\in\mathrm{rad}(\mathcal W)
 :=\mathcal W\cap\mathcal W^\perp.                \qquad(16)
```

Their independence gives the explicit coupling certificate

```math
m\le\dim\mathrm{rad}(\mathcal W)
 \le\min\{\dim\mathcal W,k-\dim\mathcal W\}
 \le1+\mathrm{rank}M.                              \qquad(17)
```

In particular, when M is zero, every pair difference is the full
A-active sum and at most one pair survives; its existence forces k
even. A proper-support pair requires nonzero coupling from the added
actual native band. Equation (17) permits arbitrary coupling rank; it
does not turn one additional native layer into bounded rank.

There is also a physical second-middle interpretation, with its
provenance now established. The assumed J-permutation at B makes W
in (5) a Hermitian Clifford. Both W and S are Clifford, so actual
one-band matching classifies the D axes into W-singletons and pairs.
Let R_D span those original physical D-pair differences, and let
L_- be the original B-negative-singleton span. The singleton axes
commute the A-active list, and their binary clean displacement obeys

```math
N_W=N_A+N_-,\qquad
N_A d=N_Wd+N_-d\in R_D+L_-.
```

Here $`N_Aq=\sum_i\langle a_i,q\rangle a_i`$ retains the
actual active list. Equation (13) consequently gives

```math
R_B\subseteq\langle a_{\rm total}\rangle+R_D+L_-.
                                                               \qquad(18)
```

The original B-list independence gives $`R_B\cap L_-=0`$.
If m_D is the number of these genuine D pairs, then
$`m\le1+m_D`$. This counts actual physical pairs under the
specified clean comparison; it does not count later logical probe
orbits. No symmetric reverse theorem is asserted without checking
its own endpoint-span and center hypotheses.

## 6. Two native bands already defeat the unqualified reality inference

The complete quadratic-query construction at n=2 supplies an existing
full-space counterguard. Its phase form is

```math
U=CCZ_{x_1,x_2,z}\otimes I_d,
```

with arbitrary returned dirty d and two independent native lists. Their
literal parity axes and odd coefficients are

| Band | Actual parity axes | T orientations |
| --- | --- | --- |
| A | $`d+x_1,d+x_2,d+z,d`$ | $`+,+,+,-`$ |
| B | $`d+z+x_1,d+z+x_2,d+x_1+x_2,d+z+x_1+x_2`$ | $`-,-,-,+`$ |

All additions in these axis labels are binary. The retained native
identity verifies all physical columns and the actual inverse. Take
the complete Boolean diagonal target to be U and both comparison
Paulis to be X_z. Then

```math
UX_zU^\dagger=U^\dagger X_zU=X_zCZ_{x_1,x_2}\ne X_z.             \qquad(19)
```

There is one activated axis in the first list and three in the second.
Their sum is zero across the two bands, despite independence within
each list. The Boolean derivative is a nontrivial quadratic sign.
Thus the one-band conclusion that every activated axis must disappear
does not extend to two bands, even with full operator equality, exact
dirty return and an actual Clifford derivative. This is a limitation
of that inference, not an unrestricted depth lower bound.

## 7. A literal four-band Boolean identity saturates the certificate

For any $`m\ge1`$, take $`k=2m+1`$ physical wires. Let one wire
be initialized, choose a fixed Clifford K sending its Z to
$`P=X_1\cdots X_k`$, and set

```math
\begin{gathered}
K_A=K_C=K,\qquad
E_A=E_C=\prod_{i=1}^kT_{Z_i},\\
B_j=X_{2j-1}X_{2j},\qquad D_j=Y_{2j-1}Y_{2j},\\
E_B=\prod_{j=1}^m T_{B_j}T_{D_j},\qquad E_D=E_B^\dagger.
\end{gathered}                                                   \qquad(20)
```

Each displayed band is an actual independent native list. Inside each
two-wire block, B_j and D_j are commuting independent Paulis; disjoint
blocks give independence of the full B and D lists. The literal clean
image is

```math
J=\bigotimes_{i=1}^k\frac{X_i+Y_i}{\sqrt2},
\qquad JB_jJ=D_j,
\quad JD_jJ=B_j,
\quad \Gamma_j=B_jD_j=-Z_{2j-1}Z_{2j}.             \qquad(21)
```

These are genuine original physical B pairs, with no singleton or
orientation repair. Their differences lie in both endpoint active
spans, and every band fixes each pair parity. Nevertheless

```math
V=K^\dagger E_A^\dagger E_B^\dagger E_BE_AK=I       \qquad(22)
```

on the entire physical register, including both clean sectors. Every
logical, dirty and reference input returns exactly, with literal phase
one and zero error. The actual reverse word is also exactly identity.

The matrix M has two copies of each disjoint weight-two row v_j.
Its rank is m. The all-ones e has odd weight and lies outside their
span; it is orthogonal to each v_j. Hence

```math
\mathrm{rad}(\mathcal W)=\mathrm{span}\{v_1,\ldots,v_m\},
\qquad \dim\mathrm{rad}(\mathcal W)=m.             \qquad(23)
```

Thus (17) is attained, while every pair has proper support and the
original physical pair count grows. Already m=1 gives a three-wire
four-band source with one proper-support pair and endpoint active
rank three.

The displayed word has $`2k+4m=8m+2`$ T gates and at most four
native bands, one clean qubit and $`k-1`$ arbitrary logical wires.
A fixed two-wire Clifford exposes each pair of native axes; a fanout
Clifford exposes P. Consequently the displayed Clifford cost and its
actual inverse cost are $`O(k+m)`$, with no added helper. These are
ledgers of the exhibited word. Its optimized implementation is the
empty identity circuit, with depth zero.

Accordingly, this family disproves extension of the three-band
**original-representation collapse** to four-band expressions, even
under the preserved-pair-code premises. It says nothing against a
minimal-depth lower theorem or a different reduction of the same
operator. A viable next step must track the actual phase or identify
a removable interaction; pair count by itself cannot be its hardness
parameter. The coupling certificate provides a way to locate that
interaction without replacing it by an unpaid native band.

The arguments above are analytic. The existing two-band quadratic
fixture supports (19). The [smallest additional fixture](../../tests/test_four_band_original_pair_scout.py)
checks the m=1 three-wire word in (20): all eight physical columns,
the actual inverse, native-list sizes $`[3,2,2,3]`$, proper support,
the signed J-pair identities and exact pair-center fixation. It is a
representation guard, not a degree enumeration or a minimal-depth claim.

## 8. Remove common actual pairs before studying the remaining interaction

There is a paid reduction that removes the inverse-canceling mechanism
of Section 7, including its signed and noninverse variants. Suppose an
actual D axis has the binary label of one member of an original B pair.
Because every B axis commutes the negative B singleton band, W acts on
that B pair exactly as J does. Actual D matching under W therefore
forces the other member of the same pair into the D list as well.

Both common axes commute every axis in both middle lists. Their four
odd native rotations may consequently be moved together without changing
the word. For one common axis P, with its D representative $`\sigma P`$,
the literal product is

```math
T_P^\epsilon T_{\sigma P}^\delta
 =\omega^{\delta(1-\sigma)/2}
   T_P^{\epsilon+\sigma\delta},
\qquad \epsilon,\delta,\sigma\in\{1,-1\}.         \qquad(24)
```

The exponent on the right is even. Applying (24) to both pair axes
gives an actual Clifford Q, including its scalar. It commutes all the
remaining middle axes. Move Q to the left and set

```math
E_C'=Q^\dagger E_CQ,\qquad K_C'=Q^\dagger K_C.
```

This is an exact full-operator identity after deleting the common pair
from both middle lists. Four T gates disappear; no new T gate or helper
is introduced. Every remaining list is an actual independent sublist,
and the conjugated output list stays independent. Q fixes every remaining
B-pair parity, so both endpoint-span and added-band-center hypotheses
remain valid. Literal error and every return condition are unchanged.
The actual inverse is the reversed transformed word with all scalars
conjugated.

For h common pairs, the transformation removes $`4h`$ T gates and
uses at most four native bands. The joint Clifford correction has an
$`O(q^2+h)`$ elementary Clifford implementation on q physical wires,
with no extra workspace: expose its independent commuting Pauli axes,
apply the even powers, and undo the fixed coordinates, retaining its
Clifford scalar. The same order applies to the actual inverse. If B
contains only such common pairs, its whole band disappears and the
complete source has an actual three-native-band representation. The
retained three-band theorem then applies to that new full word.

After removing every common pair, a residual B pair has neither axis
in D. In that case $`E_DE_B`$ cannot itself be Clifford. If it were,
the one-band matching theorem applied to
$`E_D=(E_DE_B)E_B^\dagger`$, with the other Clifford equal to I,
would force the two native axis sets to coincide. This excludes a
literal Clifford absorption of the remaining middle product. It does
not exclude a different three-band representation of the whole source.

Here is a concrete next question within this scout's actual-code class.
Assume B has one genuine pair and no negative singleton, and remove all
common pairs as above. Then W=J. If D had no genuine pair, every actual
D axis would be fixed by J at the binary level, giving M=0; Section 5
would force the B pair to have all-active support. Thus a remaining
proper-support pair needs at least one genuine original D pair.

In the smallest such case, D has exactly one pair, with difference r_D.
Its singleton rows in M are zero, and its two pair rows coincide. The
certificate gives

```math
r_B=r_D\quad\text{or}\quad
r_B=a_{\rm total}+r_D                              \qquad(25)
```

for a proper-support B pair. The complement alternative requires even
input active rank, since both pair differences are clean-even. These
are actual physical differences under the two specified clean probes.
Their equality does not imply equality of their actual native axes.

The next bounded question was whether the complete Boolean identities
permit this disjoint-axis, single-pair interaction in either case of
(25), or force another exact removal. Its scalar-sector transport is
the actual restricted two-band operator (9); it cannot be replaced by
a Pauli or a vanishing derivative in advance. Section 6 already shows
that two native bands can produce a real quadratic derivative. The
new obligation is to reconcile that derivative with the other genuine
pair sector using the full physical source. A generic feature count,
or the representation pair count alone, does not resolve this question.

## 9. Disjoint actual pairs can complete a nontrivial Boolean source

The equal-difference alternative in (25) does occur with complete Boolean
action. Take physical wires $`x_1,x_2,z,d,w`$, initialize x_1 to zero,
and return the arbitrary dirty wire d. The logical inputs are x_2,z,w.
Use the literal n=2 construction of
[the quadratic native query](QUADRATIC_NATIVE_QUERY.md#2-a-complete-native-word-for-an-even-arity-symmetric-quadratic),
without its two output Hadamards. Its two actual native lists are

```math
\begin{array}{c|c|c}
 &\text{positive Hermitian Z-parity axes}&\text{T exponents}\\ \hline
B&d+x_1,\ d+x_2,\ d+z,\ d &(+1,+1,+1,-1)\\
D&d+z+x_1,\ d+z+x_2,\ d+x_1+x_2,\ d+z+x_1+x_2
 &(-1,-1,-1,+1).
\end{array}                                                     \qquad(26)
```

All sums in (26) are binary parities. Each four-axis list is independent,
and their actual axis sets are disjoint. Both lists are exposed by the
explicit invertible CNOT words in the cited construction. Their complete
product, including both values of d and literal scalar phase, is

```math
U:=E_DE_B=CCZ_{x_1,x_2,z}\otimes I_{d,w}.          \qquad(27)
```

Set the two endpoint lists and endpoint coordinates to

```math
a_1=Y_1X_2,\quad a_2=X_1Y_2,\quad a_3=Z_w,
\qquad E_A=E_C=T_{a_1}T_{a_2}^{\dagger}T_{a_3},
\qquad K_A=K_C=K=H_w\,CNOT_{w\to1},
\qquad P_A=P_C=Z_1X_w.                            \qquad(28)
```

The chronological K word is CNOT followed by H. The endpoint axes
commute, are independent, and all anticommute with P_A. They also have
a literal Clifford encoder. On wires 1,2, let Q be the chronological
word $`X_2,\ CNOT_{1\to2},\ H_1,\ CNOT_{1\to2}`$; it sends
$`Z_1,Z_2`$ to $`X_1X_2,Y_1Y_2`$ with positive signs. Follow it by
$`S_1,X_2`$ to obtain R. Then R sends the same Z axes to a_1,a_2
with positive signs and fixes Z_w. Each X is the literal HSSH word.
Thus the chronological word $`R^\dagger,(T_1,T_2^\dagger,T_w),R`$
exposes (28) as one actual three-target native band; no scalar is
discarded.

Direct Pauli multiplication gives

```math
J=E_AP_AE_A^\dagger
 =\frac{Z_1+Z_2+X_1X_2+Y_1Y_2}{2}
       \otimes\frac{X_w+Y_w}{\sqrt2}
 =(SWAP_{12}CZ_{12})\otimes\frac{X_w+Y_w}{\sqrt2}. \qquad(29)
```

Consequently J positively swaps the first two B axes in (26), positively
swaps the first two D axes, and positively fixes each remaining axis.
The T-dagger orientation of a singleton does not make its J sign
negative. There are no negative J singletons, so W=J. Both original
physical pairs have

```math
\Gamma_B=\Gamma_D=Z_1Z_2=a_1a_2,
\qquad r_B=r_D=a_1+a_2\ne a_1+a_2+a_3.           \qquad(30)
```

Every band fixes this actual parity, and its difference belongs to both
endpoint active spans. In the three active coordinates the D coupling
rows are $`110,110,000,000`$; the one-dimensional radical certificate
is attained. Thus all conditional physical hypotheses of Sections 2--5
hold, with a proper-support pair and no shared actual middle axis.

The full Boolean completion is not an occupied-sector phase guess.
Put $`E_{\rm core}=T_{a_1}T_{a_2}^\dagger`$. On the positive
Z_1Z_2 sector, a_2=a_1 and $`E_{\rm core}=I`$. On the negative
sector, a_2=-a_1 and

```math
E_{\rm core}=\omega^{-1}S_{a_1}
 =\frac{I-i a_1}{\sqrt2},
\qquad U=I.                                      \qquad(31)
```

Therefore E_core commutes U on both sectors; T_Zw also commutes U.
The complete physical word is exactly

```math
V=K^\dagger E_A^\dagger U E_AK
 =K^\dagger U K,
\qquad
V|x_1,x_2,z,d,w\rangle
 =(-1)^{(x_1\mathbin\oplus w)x_2z}
       |x_1,x_2,z,d,w\rangle.                    \qquad(32)
```

On x_1=0 this is the nontrivial Boolean phase CCZ(w,x_2,z), with exact
clean and dirty return and every reference extension. Equation (32)
also specifies the entire unused clean block. The error is zero and
the actual reversed native word is its literal adjoint.

The displayed source uses one clean and one dirty helper, fourteen T
gates and native-list sizes $`[3,4,4,3]`$. All coordinate changes and
their inverses are the fixed elementary Clifford words above. Commuting
the complete U through E_A removes both endpoint bands, leaving the
same full operator with eight T gates and two four-target bands. This
is a whole-word endpoint removal, not the common-axis repair of Section
8. No claim that the displayed four-band word has minimal depth follows.

The [exact fixture](../../tests/test_four_band_original_pair_scout.py)
checks all 32 physical columns, both actual reversed words, the literal
endpoint encoder signs, both positive J-pair swaps, their positive
singletons, center fixation, list independence and disjointness, and
both displayed native ledgers. Complete Booleanity therefore does not
force the two middle lists to share an actual pair, even in this smallest
equal-difference case. A different removal can still retain the old cost.

## 10. A proved stopping point for the fixed-endpoint diagonal route

For the precise endpoints (28), the removal in (32) is forced whenever
the complete middle operator U is computationally diagonal. This holds
even if its displayed middle bands are not individually diagonal, and
it does not require Booleanity on the unused clean block in advance.

More generally, let y collect any additional untouched computational
wires, including arbitrary dirty wires. Let U be any diagonal unitary
on $`x_1,x_2,w,y`$. Suppose the occupied action of
$`K^\dagger E_A^\dagger U E_AK`$ is a complete Boolean diagonal
operator up to one common retained phase, with exact return. Then U
is independent of w and has equal diagonal entries on x_1x_2=01
and 10, for each y. Hence U commutes E_A, and the full operator equals
$`K^\dagger U K`$. If U is supplied by two actual native bands,
this is a paid two-band realization of the same complete source.

To prove the statement, put L=E_AK. For each input basis column with
x_1=0, exact occupied diagonality makes $`L|0,x_2,w,y\rangle`$ an
eigenvector of U with its required literal Boolean eigenvalue and common
phase. The support of this vector is obtained without an approximation.
The CNOT sets x_1=w, H_w gives both physical w values, and T_Zw changes
only their nonzero phases. E_core is identity on 00 and 11; on the
ordered pair 01,10 its matrix is

```math
\frac1{\sqrt2}\begin{pmatrix}1&-1\\1&1\end{pmatrix}. \qquad(33)
```

Thus the four occupied input labels $`(x_2,w)`$ have the following
physical supports, with every listed amplitude nonzero:

| Occupied input (x_2,w) | Physical x_1x_2 support | Physical w support |
| --- | --- | --- |
| 00 | 00 | 0,1 |
| 11 | 11 | 0,1 |
| 01 | 01,10 | 0,1 |
| 10 | 01,10 | 0,1 |

A diagonal U must have one eigenvalue on each row's entire support.
These supports cover every physical basis state for each y. They force
w independence and equality between 01 and 10. Such a U commutes both
T_Zw and the two-by-two block (33), proving the full-space removal.
Its phases are also the same common phase times Boolean signs throughout
the physical space; no unused-block promise was silently inserted.

This is a limitation of these fixed endpoints with a diagonal complete
middle, not a restriction on the unrestricted native model. Changing the
two diagonal phase-polynomial lists cannot escape it. The next bounded
question is whether one equal-difference, disjoint original pair on each
middle band can give complete Boolean action with these endpoints when
the complete middle U is genuinely nondiagonal, or whether that same
physical pair structure forces another exact removal. The first new
obligation is the off-diagonal two-band transport itself; a reality or
Pauli conclusion cannot be borrowed from the three-band proof. This
scout does not infer an unrestricted four-layer degree bound from this
route failure. The following sections resolve the nondiagonal existence
question and identify its remaining actual-list obligation.

## 11. A nondiagonal middle defeats endpoint commutation, but retains two bands

There is also a literal nondiagonal completion with the same endpoints
(28). Set $`p=X_1X_2`$, $`r=Z_1Z_2`$ and $`h=Z_d`$. Products in
this section are literal signed Hermitian Pauli products; in particular
$`pr=-Y_1Y_2`$. Replace the two middle lists by

```math
\begin{array}{c|c|c}
 &\text{actual signed axes}&\text{T exponents}\\ \hline
B&hp,\ hpr,\ hZ_z,\ h &(+1,-1,+1,-1)\\
D&hZ_zp,\ hZ_zpr,\ hr,\ hZ_zr &(-1,+1,+1,-1).
\end{array}                                                     \qquad(34)
```

These are again disjoint independent four-axis lists. To expose them,
start with the diagonal parity bases $`(d+x_1,d+x_1+x_2,d+z,d)`$
and $`(d+z+x_1,d+z+x_1+x_2,d+x_2,d+z+x_2)`$, with the coefficients
in (34). Both four-by-four binary maps are invertible, so their literal
CNOT elimination words give native encoders. Conjugate each phase band
by Q, whose chronological word is $`CNOT_{1\to2},H_1`$:
$`Q^\dagger Z_1Q=p`$, $`Q^\dagger Z_2Q=r`$. The signed axes
in (34), including the negative YY representatives, follow exactly.

Before conjugating by Q, these are just the eight commuting axes of
(26) repartitioned between two independent lists. The old d+x_2 axis
moves from B to D, and the old d+x_1+x_2 axis moves from D to B, with
their original coefficients. Consequently their literal complete product
is

```math
U=E_DE_B=Q^\dagger CCZ_{1,2,z}Q
       =CSWAP_{z;1,2}\otimes I_{d,w}.             \qquad(35)
```

The actual clean probe satisfies $`JpJ=-pr`$ and $`JprJ=-p`$.
Thus both first-two-axis orbits in (34) are genuine negative signed
pairs, with $`\Gamma_B=\Gamma_D=-r`$. Every remaining singleton
has positive J sign. Their shared binary pair difference is still the
proper endpoint support a_1+a_2, and every native axis preserves the
actual pair center. No common-axis repair from Section 8 applies.

On r=+1, SWAP and E_core are both identity. On r=-1, SWAP is p and
E_core is the rotation (33). Therefore, as complete operators,

```math
E_A^\dagger SWAP_{12}E_A=Z_1CZ_{12},
\qquad
E_A^\dagger U E_A=(-1)^{x_1(1+x_2)z},
\qquad
V=(-1)^{(x_1\mathbin\oplus w)(1+x_2)z}.          \qquad(36)
```

Here Boolean sums in exponents are binary. This time E_A does not
commute U. Nevertheless the middle expression in (36) is
$`X_2CCZ_{1,2,z}X_2`$. Substituting its retained literal two-band
dirty CCZ source gives the full operator V with two bands and eight
T gates. The four-band source again has fourteen T gates and native
sizes $`[3,4,4,3]`$, one clean and one arbitrary returned dirty helper.
All coordinate words and their inverses are paid fixed Cliffords; every
error is zero, every reference is preserved, and (36) specifies both
clean blocks. The actual reversed native words retain their adjoint
phases.

The additional [32-column fixture](../../tests/test_four_band_original_pair_scout.py)
checks the signed encoders, the negative J-pair signs and centers,
positive singletons, actual middle controlled SWAP, noncommutation with
E_A, full Boolean action, both native ledgers and both actual inverses.
Thus lack of common axes and lack of endpoint commutation do not imply
an increased optimal cost. This example requires a different exact
Clifford-coordinate reduction.

## 12. Exact structure of the fixed-endpoint, equal-difference class

Now retain exactly the endpoints (28) and the complete physical form
$`V=K^\dagger E_A^\dagger E_DE_BE_AK`$. Assume that J actually
permutes each of the two native lists up to signs, with exactly one
genuine pair in each list, both of binary difference r=Z_1Z_2, and
only singleton orbits otherwise. This is an explicit actual-list
hypothesis, not a consequence of arbitrary four-band clean return.
The center and native independence premises are retained. Let y denote
the remaining logical and arbitrary dirty wires.

For every such middle axis, its displacement under J is either zero
or r. From (29), the X_1 and X_2 displacement coordinates are both
the sum of its X_1 and X_2 coordinates, and its Z_w displacement
coordinate is its X_w coordinate. Both must vanish. Thus every
actual middle axis commutes with r and Z_w, and so does U=E_DE_B.

Assume exact occupied Boolean diagonality with the full return
contract. Fix r and physical w. The occupied vectors used in Section
10 give a complete basis for the remaining pair wire and y in that
block. Since U cannot mix physical w values, the two nonzero w
amplitudes in each occupied vector force both w blocks to have the
same eigenvalue on every one of these basis vectors. Hence U is
independent of w. On r=+1 it is diagonal in Z_1 and computational y;
on r=-1 it is diagonal in p=X_1X_2 and computational y. This conclusion
uses all occupied columns, not a promise about the unused block.

Equivalently there are Boolean functions A,B,C,D of y and the one
retained scalar alpha such that the complete middle operator is

```math
U=\alpha(-1)^{A(y)}r^{B(y)}CZ_{12}^{C(y)}SWAP_{12}^{D(y)}.
                                                               \qquad(37)
```

This is a structural equality, not a free implementation of arbitrary
controlled functions. Its factors commute. The physical operator is
a common scalar times a signed permutation; the output coefficients
are independent of every returned dirty wire. The actual source has
occupied Boolean exponent

```math
f(x_2,w,y)=A(y)+B(y)(x_2+w)+D(y)w
                         +(C(y)+D(y))x_2w.       \qquad(38)
```

Indeed E_A fixes r and CZ and takes SWAP to Z_1CZ by (36). Conjugating
the resulting diagonal phase by K gives (38), with its actual common
phase. This also determines the complete unused clean block by replacing
w in (38) with x_1+w.

The actual two-band factorization imposes two further constraints.
First, every B singleton axis commutes with every D singleton axis
as physical Paulis. To prove this, restrict the actual preserved r
code. The two axes of each native pair become equal up to a retained
sign, so their two odd rotations multiply to a scalar or a Clifford.
The remaining singleton list stays independent in the quotient: a
relation summing to r would, together with the original pair, contradict
the original physical list's independence. In either sector write

```math
U_r=F_D D_s B_s F_B,                              \qquad(39)
```

where the F factors are the literal pair Cliffords, including scalars,
and D_s,B_s are actual independent native singleton lists on the code.
This is a legitimate full-space two-band factorization on that code,
diagonal in the Pauli basis just identified. Apply
[two-band diagonal rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md#1-model-and-full-space-diagonal-theorem)
with encoders $`U_A=B_sF_B`$ and $`U_B=D_s^\dagger F_D^\dagger`$.
Its common-commuting-family conclusion compares the actual singleton
axes, without conjugating them by the outer F factors. Their quotient
commutation is their physical commutation, proving the claim.

Second, the binary sums of the two singleton lists agree. Let sigma
fix i and send sqrt(2) to minus sqrt(2). For any actual signed Hermitian
Pauli P and either odd orientation,

```math
\sigma(T_P^{\epsilon})=P T_P^{\epsilon}.
```

Let P_B,P_D be the literal products of all native axes in their
respective lists. They commute with their own bands. Equation (37)
gives $`\sigma(U)=\kappa U`$ for the scalar
$`\kappa=\sigma(\alpha)/\alpha`$. Substituting the actual native
factorization and canceling it yields

```math
P_D P_B=\kappa I.
```

Thus their labels agree and kappa is +1 or -1. Each pair contributes
the same binary difference r, so those contributions cancel, leaving

```math
\sum_{b\in B_s}b=\sum_{d\in D_s}d.             \qquad(40)
```

No assumption that alpha is a root of unity is needed. Equation (40)
is a real cross-band relation; it does not spoil independence within
either original native list. Nor does a common commuting singleton
family become one native band for free. The remaining noncommuting
interaction is confined to the four paired rotations and their
couplings to that commuting singleton family.

## 13. An intermediate cross-sector criterion, resolved in Section 14

The legitimate sector factorization (39) also gives a small conditional
degree statement. A literal Boolean sign with two native bands and no
initialized code work has degree at most three. On r=+1, (37) has phase
$`A(y)+qC(y)`$ in the Z_1 bit q; on r=-1 it has phase
$`A(y)+B(y)+qD(y)`$ in the p bit q. Therefore

```math
\deg A,\deg B\le3,\qquad \deg C,\deg D\le2,
\qquad \deg f\le4.                              \qquad(41)
```

This is a theorem only for the fixed endpoints and actual orbit
hypotheses of Section 12. It is not an unrestricted four-layer bound.
The code restrictions in its proof do not supply a free physical
replacement circuit on all center sectors.

Write B_3 for the homogeneous cubic part of B, and C_2,D_2 for the
homogeneous quadratic parts of C,D. Since the two endpoint logical
bits are independent, the complete quartic part in (38) is exactly

```math
f_4=(x_2+w)B_3+x_2w(C_2+D_2).                    \qquad(42)
```

A reduction to cubic occupied degree in this class is therefore
equivalent to the additional native implication

```math
B_3=0,\qquad C_2=D_2.                            \qquad(43)
```

The pair-center equality, physical singleton commutation and sum
relation (40) alone did not prove (43). At this intermediate checkpoint,
the bounded research question was to derive (43) from the actual
two-band lists, or exhibit a literal list violating it. Section 14
resolves it by using the native axes active on the returned address
algebra. Its independently controlled SWAP and CZ components must
still be synthesized and charged; equation (37) alone does not supply
them.

There is a concrete reason a fixed Clifford on the two pair wires
cannot simply diagonalize both components at once. SWAP has the singlet
as its unique negative eigenstate, while CZ has |11> as its unique
negative eigenstate. Their binary stabilizer spaces intersect in rank
one. Two distinct computational basis states instead have the same
rank-two stabilizer-label space. Clifford conjugation preserves this
intersection rank, so no such fixed two-wire Clifford exists. This
obstruction says nothing against a different complete native realization,
Cliffords involving other wires, or reductions such as (36).

Likewise one cannot fix additional y coordinates to constants and then
import a zero-clean query lower bound without accounting for those new
initializations and the induced native-axis relations. The missing step
at this checkpoint was the actual cross-sector relation (43), not
Boolean reality by itself. The following argument supplies that step
without making any additional initialization.

## 14. Native cubic closure of the fixed-endpoint one-pair class

Under exactly the actual-list and complete-operator hypotheses of
Section 12, every axis of B commutes with every axis of D. Consequently
the complete middle U is C3. More specifically, in (37),

```math
\deg A\le3,\qquad \deg B\le2,\qquad
\deg C,\deg D\le1,\qquad \deg(CD)\le1.          \qquad(44)
```

The nonconstant linear parts of C and D must therefore be equal,
or one must vanish. In particular, independently varying linear
controls on CZ and SWAP are excluded in this class. The occupied
Boolean exponent f in (38) has degree at most three, resolving (43).
Its complete physical extension has the same degree bound.

### 14.1. Address-active axes commute with both actual lists

Use one fixed CNOT coordinate change on the pair wires so that the
address bit s is x_1+x_2 and the remaining target bit is q=x_1.
Let M be the binary Pauli space generated by r=Z_s, Z_w and all
computational Z_y. It has codimension one in a maximal isotropic
space. The structural identity (37) fixes every element of M, even
though the individual native bands need not preserve M. Therefore,
for every actual Hermitian Pauli P with label in M,

```math
E_BPE_B^\dagger=E_D^\dagger P E_D.               \qquad(45)
```

These are identical Hermitian Clifford images. The intrinsic-axis
lemma for an actual independent native list identifies their
symplectic displacement image with the span of that list's axes
anticommuting with P. Hence this active span is the same in the
two lists. Every axis in either list commutes with it, because every
native list is itself commuting. Summing over P in M proves that
every axis active on any address probe commutes with every axis in
both complete lists.

Call the remaining axes M-silent. Their labels belong to M-perp.
On a fixed computational address a=(s,w,y), an M-silent axis is a
signed one-qubit Pauli on q, multiplied by a nonzero address sign.
Within each native band, all its nonidentity silent target factors
have one Pauli type: different types would anticommute, since their
address parts contain only Z operators.

### 14.2. A silent genuine pair excludes different target types

Suppose the two bands have different nonidentity silent target types
P and Q. They anticommute. The physical singleton-commutation result
of Section 12 means that at least one of these axes belongs to a
genuine pair. Its partner is also M-silent, since their difference
is r in M.

In the present fixed coordinates, the pair part of J is Z_q on s=0
and X_q on s=1; its w factor remains (X_w+Y_w)/sqrt(2). For a silent
axis, the target types X or Z have binary displacement r under J,
while I or Y have displacement zero. Address Z factors may alter a
literal sign but not that displacement. Thus a silent genuine pair
has target type X or Z. A further silent singleton with target type
Y would anticommute with it. A further X or Z orbit would be another
genuine pair, contrary to the one-pair hypothesis. The only silent
nonidentity target factors in that band are therefore its two pair
axes.

On a fixed r sector their two odd rotations combine literally into
a scalar or an S_P power. Exactly one of the two r sectors has an
S_P or S_P-dagger factor: the two possible combined T exponents are
zero and plus or minus two. On that sector this band's target action,
up to its retained scalar, is $`e^{-i\epsilon\pi P/4}`$ for
$`\epsilon\in\{1,-1\}`$. The other band's silent action has the
form $`e^{-i\theta Q}`$ up to its retained scalar.

Factor all M-active rotations into E_act. Section 14.1 shows that
they commute with every remaining factor, so

```math
U=E_{\rm act}E_{D,\rm silent}E_{B,\rm silent}.   \qquad(46)
```

The two silent factors fix M, and U fixes M. Thus E_act fixes M too.
Moreover E_act commutes with every silent axis. In any address block
with the two nonparallel target types P,Q, its one-target block must
therefore be scalar. This is a consequence of the actual factorization,
not an assumption that the active rotations preserve the address
separately.

The classified complete U in that block is a scalar times a Pauli:
I or Z_q on s=0, and I or X_q on s=1. But a quarter rotation about
P times any rotation about the anticommuting Q cannot be a scalar
times a Pauli. For one order, writing c=cos(theta) and t=sin(theta),

```math
e^{-i\epsilon\pi P/4}e^{-i\theta Q}
 =\frac1{\sqrt2}
       (cI-i\epsilon cP-itQ-\epsilon tPQ).       \qquad(47)
```

At least one of c,t is nonzero, and each supplies two nonzero
coefficients on distinct Pauli labels. Reversing the order changes
only the sign of the PQ coefficient. No retained scalar can turn
either product into a single Pauli. This contradiction occurs on
the pair's active r sector, which is a genuine nonempty physical
sector covered by (37).

Hence the two bands cannot have nonparallel silent target types.
Every silent cross-band pair of axes commutes; every active axis
already commutes with both full lists. All actual B and D axes
therefore commute globally. Their product is a commuting Pauli
T-phase polynomial and belongs to C3. This conclusion does not
assert that their union is an independent native list, or replace
two physical bands by one for free.

### 14.3. The complete Boolean coefficients are cubic

Write the literal middle operation as

```math
U=\alpha D_g\,SWAP_{12}^{D(y)},
\qquad D_g=(-1)^{g(x_1,x_2,y)},
\qquad g=A+B(x_1+x_2)+Cx_1x_2.                   \qquad(48)
```

The phase g is symmetric in the pair wires. Since U is C3, its
conjugate of Z_1 is a diagonal Clifford. That conjugate has Boolean
exponent

```math
x_1+(x_1+x_2)D(y).
```

A Boolean diagonal Clifford has degree at most two. Independence
of the pair coordinates therefore forces D to be affine.

For any y coordinate, the permutation of $`UX_{y_j}U^\dagger`$
is the fixed Clifford $`SWAP^{\partial_jD}X_{y_j}`$, since the
derivative of the affine D is constant. Its literal Boolean phase
is $`\partial_jg`$: swapping the pair does not change g. Dividing
out that Clifford permutation leaves a diagonal Clifford, so every
$`\partial_jg`$ has degree at most two. Separating the terms with
zero, one and two pair coordinates gives

```math
\deg A\le3,\qquad \deg B\le2,\qquad \deg C\le1.
```

There is also the stated obstruction to independent controls. The
permutation of $`UX_1U^\dagger`$ flips the pair bits by
$`(1+D,D)`$, an affine Clifford because D is affine. Its literal
Boolean phase is

```math
\Delta_{(1+D,D)}g
 =B+Cx_2+CD(x_1+x_2).                            \qquad(49)
```

Removing that affine Clifford again leaves a diagonal Clifford.
The first two terms of (49) already have degree at most two;
independence of x_1,x_2 forces CD to be affine. For affine C,D,
this is equivalent to their linear parts being equal or one being
zero. This proves all of (44).

Substituting (44) in (38) gives occupied degree at most three.
The full physical Boolean exponent is

```math
A+B(x_2+x_1+w)+D(x_1+w)+(C+D)x_2(x_1+w),          \qquad(50)
```

which also has degree at most three. The conclusion retains the
actual common phase, every clean and dirty column, arbitrary
references and the actual reversed native inverse. It is a proof
about the original operator, with no extra helper, approximation
error or resource substitution. In particular, it is not a same-budget
construction obtained by freely synthesizing the coefficient functions.

### 14.4. The next bounded original-pair question

The fixed-endpoint, one-pair common-orbit class is now closed at cubic
degree. A useful continuation must leave those hypotheses. At general
four-band return, the first unsupported inference is already (2):
the intermediate image $`E_BJE_B^\dagger`$ need not be Clifford.
Thus return alone does not license original B-pair matching, much less
matching both actual middle lists under the same J. This is not a
depth obstruction in the unrestricted model.

Retain one genuine original B pair and its explicitly preserved physical
center, but allow one actual negative B singleton. Its Clifford image
W in (5) need not equal J, and the additional native list is matched
under W when that matching is justified, not under J. The next bounded
question is whether this smallest change permits complete Boolean
transport outside the common-orbit class, or forces an actual reduction
back to it. The real two-band scalar transport (9), literal repairs and
both pair-sign sectors must be carried through that question. Replacing
W by J or calling its projected axes independent in advance would be
the first unsupported step.

This source-specific cubic theorem does not change the complete
approximate dirty-indicator frontier, the retained unrestricted
three-layer results, or any inherited resource, error and return
recurrence.
