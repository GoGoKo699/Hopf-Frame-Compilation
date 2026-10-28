# Scratch: a weak complete-probe decoder has a Reed--Muller Pauli domain

Independent audit of Root's continuation argument. Not part of the frozen
PR269 tree. The argument is valid with the literal phases and constants below.

## 1. Contract and conclusion

Let x range over F2^p, m=2^p, and n=m+q. The target contains an m-qubit bank
indexed by x and q additional arbitrary qubits. Let Q_x be a family of
Clifford unitaries on all n qubits, whose binary symplectic matrix F_x has
entrywise Boolean degree at most d. Additivity, involutivity and an identity
helper action for Q are not needed. The inverse action has the same degree,
because F_x^(-1)=J F_x^T J for the fixed binary symplectic form J.

Let E be an arbitrary unitary, and define

    R=E^dagger Z_0 E,       R_x=Q_x R Q_x^dagger.

Assume only the complete probe equation

    max_x ||E R_x E^dagger-Z_x||_op <= epsilon.             (1)

The identity tensor factor on the q helper coordinates in Z_x is implicit.
The comparison includes all dirty and reference inputs. It does not require
E Q_x E^dagger to be Clifford or to equal the whole wire translation.

Let H_E be the binary Pauli domain: P belongs to it when E P E^dagger is
a signed Hermitian Pauli. Let

    M_d(p)=sum_(j=0)^(min(d,p)) binomial(p,j).

If epsilon < 2^(-d-1), then every image-domain Pauli has an output X mask
in the binary Reed--Muller code RM(d,p). Consequently

    dim H_E <= m+M_d(p)+2q,
    symplectic_rank(H_E) <= 2M_d(p)+2q.                    (2)

For a native Clifford+T word for E, if A is the span of its absorbed T axes,

    dim A >= m-M_d(p),       T(E) >= m-M_d(p).             (3)

These are domain, axis-rank and count restrictions, not growing depth bounds.
For an affine field shear, d=1 and the error threshold is one quarter.

## 2. Orientation and the approximate signed conjugation

Fix a Hermitian Pauli P in H_E and write S=E P E^dagger. Let a(x) be the
X-mask bit of S on target bank wire x. Thus

    Z_x S Z_x=(-1)^(a(x)) S.

Put B_x=E R_x E^dagger. Both B_x and Z_x are Hermitian unitaries.
Equation (1) and unitary telescoping imply

    ||B_x S B_x-(-1)^(a(x)) S|| <= 2 epsilon.

Conjugate this inequality by E^dagger and then Q_x^dagger. With the
Hermitian Pauli P_x=Q_x^dagger P Q_x, it becomes

    ||R P_x R-(-1)^(a(x)) P_x|| <= 2 epsilon.              (4)

This confirms the orientation: it is the inverse-query Pauli family P_x,
not Q_x P Q_x^dagger, in the signed-conjugation identity. The symplectic
inverse formula above preserves the original degree hypothesis.

## 3. A scalar cube forces every high derivative to vanish

Take any affine (d+1)-cube of addresses, with independent directions, and
choose a fixed order of its 2^(d+1) vertices. The binary label of P_x is
F_x^(-1) v, of degree at most d. Its sum over that cube is therefore zero.
The ordered Pauli product is consequently

    Omega=product_x P_x=lambda I,

where lambda is one of 1,-1,i,-i. The Paulis need not commute and lambda
need not be positive. No sign-degree hypothesis on the Clifford family is
used: only the label sum is needed.

Conjugating the entire ordered product by R gives Omega exactly. Multiplying
the approximate identities (4) in the same order costs at most
2^(d+2) epsilon. If the cube parity of a is odd, the requested product is
-Omega, at operator distance two. Hence

    odd cube parity => 2 <= 2^(d+2) epsilon
                    => epsilon >= 2^(-d-1).               (5)

Below the threshold, every (d+1)-fold affine derivative of a is zero, which
is the Boolean degree-at-most-d condition. If d>=p, RM(d,p) is the whole
space and the same conclusion is automatic.

## 4. Domain dimension and its symplectic rank

The image of H_E is a linear Pauli-label subspace, and E induces an
injective symplectic map on it. Equation (5) places that image inside

    V_allowed={(a,b,c): a in RM(d,p),
                         b in F2^m, c in F2^(2q)}.

This space has dimension M_d(p)+m+2q. Its radical consists of bank Z masks
b orthogonal to every Reed--Muller X mask; the helper component is zero.
The radical dimension is m-M_d(p), so

    symplectic_rank(V_allowed)=2M_d(p)+2q.

Restriction to a subspace cannot increase the rank of the alternating form,
proving both statements in (2). Finally A^perp is contained in H_E because
every Pauli commuting with all absorbed rotation axes remains a Pauli after
E. Thus 2n-dim A <= m+M_d(p)+2q, proving (3).

## 5. Independent four-wire cells fail even this weaker contract

The retained four-wire six-T cell D has Pauli domain

    span{X0,X2,X3,X1 Z0},

of dimension four and symplectic rank two. Let m=4r=2^p and consider

    E=L (D^tensor r tensor I_q) K,

with arbitrary fixed joint input and output Cliffords K,L, allowed to mix
every cell and helper. Its domain contains a subspace of symplectic rank
2r+2q. Clifford walls preserve that rank. The weaker complete-probe contract
with a degree-d Clifford query is therefore impossible below 2^(-d-1) whenever

    2r+2q > 2M_d(p)+2q,
    equivalently m > 4M_d(p).                             (6)

For the field shear, M_1(p)=p+1. Condition (6) holds for every p>=5.
Thus independent copies of the four-wire cell, even with arbitrary joint
Clifford walls and any number of identity spectator helpers, cannot satisfy
the complete-probe contract below error one quarter from five address bits
onward. This is stronger in contract scope than the full-intertwiner screen;
it does not require the whole query endpoint to be a wire permutation.

The conclusion does not exclude coupling the non-Clifford bands themselves,
or appending non-Clifford helper operations that change the Pauli-domain rank.
It does not exclude every two- or three-band decoder and gives no general
native frame-depth lower bound. A successful weak-probe candidate may exploit
its large target commutant, but its exact Pauli domain still obeys (2).

## Audit checks

The companion scratch fixture checks the field cube label sums and literal
phases; an explicit noncommuting cube with product -I; Reed--Muller dimensions
and symplectic radicals; and tensor-cell symplectic rank under joint Clifford
changes. It is a bounded exact supplement to the proof, not an asymptotic
search or native construction. No current PR file was modified.

## 6. Exact covariance gives a stronger dimension bound

Assume now epsilon=0 and that Q is an additive representation, so
Q_y Q_x is Q_(x+y) (a scalar projective factor would also cancel under
conjugation). For a bank Z mask v define

    Z^v=product_x Z_x^(v_x),
    K_Z={v: E^dagger (Z^v tensor I_q) E is a Pauli}.

Exact orbit covariance gives Q_y R_x Q_y^dagger=R_(x+y). Because the
R_x are exact pullbacks of the commuting Z_x, it follows that

    Q_y E^dagger Z^v E Q_y^dagger=E^dagger Z^(Pi_y v) E.    (7)

Thus K_Z is a translation-invariant binary subspace. Its pullback Pauli-label
map phi is linear and injective: a nonzero Z word has trace zero and cannot
be pulled back to a scalar Pauli. Taking labels in (7) gives

    phi(Pi_y v)=F_y phi(v).

For any d+1 independent directions, additivity and the degree-d hypothesis
give product_j(F_yj-I)=sum_vertices F_x=0. Injectivity of phi transfers this
identity to the corresponding translation differences on K_Z. Therefore

    K_Z subset RM(d,p),       dim K_Z <= M_d(p).            (8)

Project the image-label space of H_E onto its bank X mask and all helper
labels. Its kernel is exactly K_Z, interpreted as pure bank Z labels.
The image dimension is at most M_d(p)+2q by the earlier probe theorem.
Consequently exact additive queries satisfy

    dim H_E <= 2M_d(p)+2q,
    dim A >= 2m-2M_d(p).                                  (9)

This exact argument does not assert approximate invariance of H_E. Replacing
(7) by a product of approximate probe identities costs up to wt(v)*epsilon.
Sections 8--9 retain that weight cost explicitly and use a short codeword
to obtain a uniform bound for independent cells.

## 7. The field's common fixed kernel sharpens the exact result further

For a full-rank field shear, d=1 and K_Z consists of affine masks by (8).
If p>=2, every affine mask v has a nonzero translating stabilizer: its one
linear slope has a nonzero kernel. For such a nonzero y,

    Pi_y v=v => F_y phi(v)=phi(v).

Every nonzero field query has the same fixed-label kernel. Hence every
F_x fixes phi(v). Injectivity and (7) then imply Pi_x v=v for all x.
Thus v is constant. This proves the sharper statement

    K_Z subset span{the all-ones bank mask},
    dim K_Z <= 1.                                        (10)

The field degree may be enlarged and arbitrary identity source helpers may
be present; only the common nonzero fixed kernel and affine label action
are used. More explicitly, a degree-K source on 2K core qubits plus h
identity source helpers has n=2K+h=m+q. Here q counts the physical coordinates
outside the target bank; it need not equal h or be query-fixed. The weak
contract allows any coherent query action on those remaining coordinates,
provided the probe comparison (1) holds on their entire space. No target
helper identity assumption is introduced by the kernel projection.
Combining (10) with the X-mask projection bound yields

    dim H_E <= p+2+2q,
    dim A >= 2m-p-2,       T(E) >= 2m-p-2,                 (11)

for every exact weak complete-probe decoder with p>=2. The full query
endpoint still need not be Clifford or a wire translation.

For independent four-wire D cells, dim H_E is at least m+2q. Since
m>p+2 for every p>=3, that entire class fails the EXACT weak-probe contract
already at p=3, even with arbitrary joint Clifford walls and arbitrarily
many identity spectator helpers. The p=2 finite decoder is consistent:
its domain dimension four reaches the bound p+2=4.

The exact p=2 native full intertwiner independently saturates every step:
its image Pauli domain has dimension four; its bank X-mask projection is
all of RM(1,2), of dimension three; its pure bank Z kernel is exactly
{I,Z_all}, of dimension one. Tensoring identity spectators adds 2q to the
domain and projection but does not enlarge that pure-bank kernel. Thus the
dimension count is compatible with the retained positive finite example.

The exact covariance argument alone gives no uniform positive error floor.
The following weight-dependent argument and code bound supply one without
asserting approximate invariance of the domain.

## 8. A robust low-weight pure-Z kernel lemma

Keep the literal additive full-rank field query, p>=2, but allow a nonzero
probe error epsilon. Define K_Z by exact Pauli pullback as before; it remains
a linear space, even though it need not be translation-invariant. If a
nonzero v in K_Z has weight w, then

    epsilon < 1/(3w) => v is the constant all-ones mask.    (12)

To prove this, put P=E^dagger Z^v E, an exact Hermitian Pauli. Choose a
fixed ordering of the w occupied addresses and form C_v as the product
of their R_x in that order. The probe equation gives

    ||C_v-P|| <= w epsilon.

Exact query covariance translates that ordered product. Its target factors
Z_(x+y) commute, so their product is Z^(Pi_y v) regardless of the translated
order. A second telescoping gives, for every y,

    ||Q_y P Q_y^dagger-E^dagger Z^(Pi_y v)E||
         <= 2w epsilon.                                  (13)

No reordering of the approximate R_x is used. At y=0 the comparison in
(13) is exact, so that vertex has zero error.

The source Paulis P_y=Q_y P Q_y^dagger commute pairwise. In field coordinates
their labels are v_P+N_y v_P, where each <v_P,N_y v_P>=0. The derivative
labels lie in the isotropic core space span{Z_u,X_v}, with no helper
component, so <N_y v_P,N_z v_P>=0 as well. Consequently every four-address
plane product P_0 P_y P_z P_(y+z) is a Hermitian scalar, +I or -I.
Its label vanishes by the affine query action.

Compare this source product with its four target Z words using (13).
The zero vertex is exact and the other three contribute at most 2w epsilon
each, for a total of 6w epsilon. If the second translation difference
(I+Pi_y)(I+Pi_z)v is nonzero, the target product is a conjugate of a
nontrivial Z word. Its distance from either +I or -I is exactly two.
Thus epsilon<1/(3w) forces all second differences of v to vanish: v is
affine.

At p>=2, an affine v has a nonzero translating stabilizer y. Equation (13)
then compares Q_y P Q_y^dagger directly with P. These source Paulis commute,
so their distance is either zero or two. Since 2w epsilon<2, P is fixed
by this nonzero query. The common nonzero fixed kernel of the field pencil
then fixes P literally under every query. Apply (13) again: unequal target
Z words also have distance two, so every translation fixes v. This proves
(12).

In particular, epsilon<1/(3m) makes every K_Z mask constant. This threshold
is also below one quarter for p>=2, so the earlier X-mask restriction applies.
The stronger dimension and axis bounds (11) therefore remain valid at

    epsilon < 1/(3m).                                    (14)

This is a width-dependent robust theorem for arbitrary E. It is distinct
from the uniform error floor for the special independent-cell class below.

## 9. Uniform one-twelfth error floor for independent cells

Let p>=3 and E=L(D^tensor r tensor I_q)K, with m=4r=2^p as in Section 5.
Here I_q describes the spectator factor in E; it does not require the
query to be identity on those physical coordinates. The field query may
act on a larger core, and the decoded endpoint may move every helper.
Only the complete-probe equation and the source field hypotheses are used.

Suppose for contradiction that epsilon<1/12. The earlier one-quarter
theorem bounds the image-domain X-mask projection by p+1 dimensions.
Including all helper labels makes that projection at most p+1+2q
dimensional. Since the independent cells have domain dimension at least
m+2q, its pure-Z kernel satisfies

    dim K_Z >= m-p-1,
    [F2^m:K_Z] <= 2^(p+1)=2m.                            (15)

This binary linear code has a nonzero word of weight at most four. Indeed,
if its minimum distance were at least five, the Hamming balls of radius
two centered at its codewords would be disjoint. Counting their elements
would imply

    1+m+binomial(m,2) <= [F2^m:K_Z] <= 2m.

But the left side minus 2m is (m-1)(m-2)/2, positive for m>2. Thus there
is a nonzero v in K_Z with weight w<=4. Since m>=8, it is not the all-ones
mask. On the other hand epsilon<1/12<=1/(3w), so (12) forces it to be
that constant mask, a contradiction.

Therefore the independent four-wire-cell class obeys the complete-probe
error floor

    epsilon >= 1/12                 for every p>=3.       (16)

The earlier symplectic-rank argument improves this to epsilon>=1/4 for
p>=5. Both statements permit arbitrary joint fixed Clifford walls and any
number of identity spectator factors in E. Neither excludes coupling the
non-Clifford bands, changing the spectator factor with a non-Clifford word,
or using a different source family. No growing native-depth lower bound or
improved frame compiler follows.
