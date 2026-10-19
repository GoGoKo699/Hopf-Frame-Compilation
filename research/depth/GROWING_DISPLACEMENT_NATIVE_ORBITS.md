# Growing middle displacement reduces to native pair permutations

[External-center correction](EXTERNAL_CENTER_MIDDLE_CORRECTION.md) ·
[Actual-source support](MIDDLE_DEFECT_SUPPORT_CLOSURE.md#8-the-original-clean-transition-removes-the-mixed-case) ·
[Partial-diagonal features](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md)

The full normalized cubic identity closes the growing-total-displacement
case of the retained actual no-middle-pair branch. Its coordinatewise
algebra forces every growing displacement direction, except one common
transvection, to be a permutation of a pair of original middle axes.
Those pair phases are retained as an actual native middle band. A
three-band version of the previous bounded-displacement factorization
then gives a uniform Boolean-degree bound.

The result uses the full literal Boolean source to obtain its tensor
identity. It is not a statement about arbitrary Gram-compatible label
models, and it does not construct a growing-degree Boolean family.
The remaining original middle-pair cases and the higher-layer asymptotic
problem are separate. The complete-frame resource frontier is unchanged.

The two uses of pairs are distinct. The excluded original middle pairs
refer to the physical clean-transition reduction. The pairs classified
here are orbits of the full logical computational-probe action on the
positive middle axes. Such logical pair orbits can occur within the
original no-middle-pair branch.

## 1. Actual-source hypotheses and conclusion

Retain the actual no-middle-pair logical normal form, its original
one-clean physical source, and all paid normalizations:

```math
V=K_C^\dagger E_C^\dagger H^{b(z)}E_BE_AK_A
 =e^{i\theta}D_g.                                      \qquad(1)
```

Each E is an original independent commuting signed Pauli T/T-dagger
list. The word acts on the entire logical register, with arbitrary dirty
inputs and references. The scalar, original central cubic signs,
quadratic normalizations, clean return and actual inverse are retained.

Consider the remaining branch in which both interior clean images are
non-Pauli, the shared-clean-flip criterion does not already close the
problem, and the quadratic control has polar rank at least four with
some nonconstant probe flag. Let beta be its polar symbol and S_beta
its support. The actual-source theorem gives S_beta as the original
clean-active intersection, embedded faithfully in the positive middle
span and in the genuine shared center R.

Use the original independent positive-middle axes as coordinates on
the binary vector space with its ordinary dot product. This metric is
distinct from physical Pauli symplectic pairing. The existing matching
and common-preimage results give

```math
\begin{gathered}
O_u=I+N_u,\qquad N_u^T=N_u,\qquad N_uN_{u'}=0,\\
U=\sum_u\mathrm{im}N_u,\qquad K=U^\perp,
\qquad U\cap S_\beta=0,\\
c(O_u)=\beta\wedge p_u,\qquad N_uv=p_u
\quad\hbox{for one fixed }v,                         \qquad(2)
\end{gathered}
```

where the full alternating tensor is

```math
c(O)(x,y,z)
 =\sum_i(Ox)_i(Oy)_i(Oz)_i+\sum_i x_i y_i z_i.          \qquad(3)
```

In addition, genuine-center fixation gives

```math
S_\beta\subseteq K,\qquad
\beta(r,x)=0\quad(r\in U).                            \qquad(4)
```

Indeed every O_u fixes S_beta pointwise, and symmetry makes its image
orthogonal to that support. Equation (4) is stronger than the disjointness
in (2); it must not be inferred from disjointness alone. The restriction
of beta to K retains its full rank by the existing actual-source
disjointness argument.

**Native-orbit conclusion.** There is a partition of the original middle
coordinates into disjoint pairs, a set I of singleton coordinates, and
other fixed singletons. Put p equal to the indicator vector of I, and
let r_j be the indicator of the jth pair. Then

```math
N_u=\alpha_u pp^T+\sum_j\epsilon_{u,j}r_jr_j^T,
\qquad p_u=\alpha_up,                               \qquad(5)
```

with binary linear characters alpha and epsilon. Some alpha_u is one,
the size of I is even and at least six, and

```math
\mathrm{rank}\beta=|I|-2.                     \qquad(6)
```

Every r_j term is the actual permutation of its two original native
axes. The only exceptional action is the common transvection on I.
No new orthogonal native basis has been chosen.

**Degree conclusion.** With the constants of the partial-diagonal
finite-feature theorem, the high-rank branch obeys

```math
\deg g\le C_*:=D_3+576.                              \qquad(7)
```

Combining the retained rank-zero, rank-two, constant-flag and shared-flip
closures gives, throughout this actual no-middle-pair reduction, the
uniform bound

```math
\deg g\le
\max\{C_*,D_{\mathrm{ins}}(0,0,0,3)\}.               \qquad(8)
```

This is independent of total displacement dimension, clean-active ranks,
quadratic support and dirty width. For a literal XOR query, first check
this reduction for its actual output-Hadamard sign representation; the
predicate degree is then at most one less than (8).

## 2. The common fixed space is one dimension from a coordinate algebra

Write x*y for coordinatewise multiplication in the original native
coordinate space. Restrict (3) to x,y in K. They are fixed by every O_u,
and every p_u is orthogonal to K. Consequently

```math
N_u(x*y)=\beta(x,y)p_u,
\qquad x*y+\beta(x,y)v\in K.                        \qquad(9)
```

Some flag is nonzero, so v is not in K. Define

```math
L=K+\langle v\rangle.                              \qquad(10)
```

Products of two K vectors lie in L. To prove that L is an algebra,
put T(z)=v*z modulo L for z in K. Expand the associativity identity
(x*y)*z=x*(y*z), using (9) for the intermediate K terms. Modulo L it
gives

```math
\beta(x,y)T(z)=\beta(y,z)T(x)
\qquad(x,y,z\in K).                               \qquad(11)
```

For any z, choose y in its beta-orthogonal hyperplane but outside the
radical of beta restricted to K. Such a y exists because that rank is
at least four. Choose x with beta(x,y)=1. Equation (11) gives T(z)=0.
Thus v*K is contained in L. Since v*v=v, L is closed under every
coordinatewise product.

It is a unital algebra. Every orthogonal binary matrix fixes the
characteristic all-ones vector: its i-th column has dot square one,
so the column sum is one. This gives O-transpose times that vector
equal to itself, and orthogonality gives the same identity for O.
Thus every O_u fixes that vector, which
therefore belongs to K. Also U is totally isotropic by the pairwise
zero-product identities in (2), so

```math
L^\perp\subseteq U\subseteq K\subseteq L.            \qquad(12)
```

A unital subalgebra of a binary coordinate algebra consists precisely
of the vectors constant on blocks of a coordinate partition. One proof
defines two coordinates to be equivalent when all algebra elements
agree there. The algebra separates distinct blocks; products of
separating elements and their complements supply each block indicator.

The coisotropy (12) allows only blocks of size one or two. Within a
block, every even-weight vector lies in L-perp. A block of size at
least three has an even-weight vector which is not constant on that
block, contradicting (12). This partition uses actual original axes.

Since K is a hyperplane in L, write K as the kernel of one nonzero
functional lambda on L, normalized by lambda(v)=1. Every pair-block
indicator lies in L-perp and hence in K, so lambda involves only
singleton coordinates. Let I be its support. The all-ones vector lies
in K, so I has even size. On K, (9) gives

```math
\beta(x,y)=\lambda(x*y)=\sum_{i\in I}x_i y_i.        \qquad(13)
```

The form on the even-parity subspace of an even set I has rank |I|-2.
The other coordinates contribute only radical directions. Full-rank
restriction of beta to K now proves (6), including the size condition.
Moreover

```math
U=\mathrm{span}\{p,r_1,\ldots,r_s\}.           \qquad(14)
```

## 3. The full tensor forbids coupling between the native pairs

Each symmetric N_u with image in U has a unique symmetric matrix
representation in the independent basis (14). Its output on either
coordinate of pair j is the same linear functional, denoted L_j.
Let q_j(x)=r_j dot x be that pair's parity functional.

Evaluate (3) with its first argument r_j. Equation (4) and total
isotropy make the right side zero. The left side is, directly on the
two coordinates of that pair,

```math
c(O_u)(r_j,y,z)
 =q_j(y)L_j(z)+L_j(y)q_j(z)=0.                     \qquad(15)
```

Therefore L_j is either zero or q_j. In the symmetric matrix of N_u,
every row corresponding to a pair has only its diagonal entry.
Symmetry also removes its coupling to the remaining p direction.
This proves the first equality of (5).

The vector v belongs to L. Its dot product with each r_j is zero,
and its dot product with p is lambda(v)=1. Thus N_uv=alpha_u p,
which proves the flag equality in (5). Linearity in u follows from
the existing native identities N_{u+u'}=N_u+N_{u'}.

The support of beta is also entirely within I. Choose u with alpha_u=1.
All disjoint pair permutations preserve the native cubic tensor, so
(3) and (5) give

```math
\beta\wedge p=p\wedge e_2(I),
\qquad \beta=e_2(I)+p\wedge a                     \qquad(16)
```

for some vector a. If a had a nonzero coordinate outside I, evaluating
the corresponding polar map on that coordinate would put p in S_beta,
contradicting p in U and (2). Thus a is supported on I. Similarly,
applying the polar map to p gives (1+p dot a)p, so p dot a=1.
Its kernel on I is the span of p and a, consistent with (6).

In particular one can now choose the common correcting vector to be
any single original coordinate e_i with i in I:

```math
N_ue_i=p_u\quad\hbox{for every }u.                 \qquad(17)
```

The corresponding correction is one actual middle Pauli. Its controlled
Boolean phase is still not a free native rotation; the degree proof
below does not synthesize that correction as a new middle band.

## 4. Retain the permutation phases inside a three-band residual

Let E_exc be the original middle sublist on I and let E_perm contain
every other original middle axis. These are literal signed native
sublists with their original T orientations. Their product is E_B.
Their binary actions are invariant, but no phase of E_perm is dropped.
Signs on its permuted axes and its Clifford probe defects stay in the
word.

On the exceptional native-coordinate space, the total displacement
dimension is one. Its common fixed space is

```math
K_{\rm exc}=p^\perp\subseteq\mathbb F_2^I.          \qquad(18)
```

Apply the paid normalization from Section 7 of the external-center
proof to this fixed family. The adaptation is as follows.

First peel the actual endpoint axes with zero input or output
computational X mask as diagonal native factors D_A,D_C. The residual
is complex diagonal; it is not assumed Boolean. Native independence
shows that every remaining endpoint axis commutes K_exc. Both middle
sublists and H^b commute K_exc as well. The actual input/output Pauli
intertwining therefore admits one literal diagonal S/CZ/Z aligner L_0.
The aligned complex diagonal residual centralizes the common family.

Use fixed affine computational coordinates and diagonal decorations T
to normalize that family jointly with the genuine center. On its
selected X modes the whole residual is exactly identity, so applying
Hadamards to those modes leaves the entire residual unchanged and makes
K_exc computational. Keep D_A,D_C,L_0 conjugated only by T. This is
the same whole-word identity-mode argument as the previous proof; no
separate factor is conjugated back through those Hadamards.

The common matrix-entry mask space for the restricted actions on I has
dimension at most one. Every remaining endpoint axis active against
E_exc has its input X mask in that same space, by native independence.
An axis active against H is already active against E_exc: compare its
native coefficient in the full physical identity

```math
(\Phi_{J_u}-I)h=N_ue_i=p_u                         \qquad(19)
```

and its output-family counterpart. Here N_u acts on the actual middle
axis e_i, not on an external h. Every remaining endpoint activation
mask is nonzero, so the comparison is valid.

Split both remaining endpoint lists into active and passive sublists
relative to E_exc and H. Define, in these actual normalized frames,

```math
\begin{aligned}
A_{\rm in}&=K_A^\dagger A_{\rm act}K_A,&
C_{\rm out}&=K_C^\dagger C_{\rm act}K_C,\\
D_{\rm in}&=K_A^\dagger(H^bE_{\rm exc})K_A,&
W_3&=K_C^\dagger C_{\rm pass}^\dagger
                E_{\rm perm}A_{\rm pass}K_A.
\end{aligned}                                                        \qquad(20)
```

The exact factor order is

```math
T V_0T^\dagger
 =C_{\rm out}^\dagger W_3D_{\rm in}A_{\rm in}.      \qquad(21)
```

Only the passive endpoint factors move past H^b E_exc. The entire
E_perm remains as the middle band of W_3. Accordingly W_3 has three
actual independent native bands, rather than the two bands of the
previous all-middle factorization. This is exactly what the existing
partial-diagonal theorem permits.

The endpoint active masks need at most one direction. The quotient
of the exceptional middle by K_exc needs one more, and H needs at
most one further direction; its quadratic controls remain genuine
computational parities. Choose the computational-Z subgroup M_0
annihilating these at-most-three masks. Then A_in,C_out,D_in fix M_0
pointwise. Diagonality of the complete residual in (21) forces W_3
to fix M_0 pointwise on the full register as well.

The partial-diagonal three-band theorem supplies its full matrix
features with bound D_3. There are six other commuting-C3 factors:
A_in,C_out,D_in and the restored D_A,D_C,L_0. Each is internally a
commuting Pauli phase and fixes M_0; the six factors need not mutually
commute. For each factor, the retained feature calculation costs at
most three times eight Boolean features of degree at most four.
Multiplying the matrices in their actual order therefore adds at most

```math
6\cdot3\cdot4\cdot2^3=576                         \qquad(22)
```

to the Boolean degree bound. This proves (7). Restoring the affine T
coordinates cannot increase the Boolean feature degrees. The retained
central cubic and real quadratic signs do not increase this bound.

All normalization and alignment words use fixed Clifford operations,
with O(Q squared) gate cost on Q logical wires, no helper, and their
actual reversed inverse words. The argument is a full-space degree
factorization, not a new native compiler or a free conditional decoder.

## 5. Scope of the resulting obstruction

The mathematical uncertainty changed here is whether growing total
middle displacement can hide arbitrarily many non-permutation native
directions in this actual Boolean no-middle-pair branch. It cannot:
the global phase tensor, genuine-center fixation and source independence
leave one exceptional direction. Total displacement itself may still
grow through genuine native pair permutations.

The obstruction applies to every actual source satisfying the stated
normal-form hypotheses. It does not merely forbid one chosen native
basis for a corrected residual. Conversely, it does not cover the
remaining original middle-pair branch, all one-clean three-layer
circuits without that reduction, or higher-layer asymptotics.
Those are the next lower-bound questions. R96 and the complete-frame
resource frontier are unchanged by this component theorem.

The [bounded exact fixtures](../../tests/test_growing_displacement_native_orbits.py)
check a six-axis transvection together with a disjoint native pair,
the full trilinear identity, the common-preimage product algebra,
literal signed native defects and their actual inverses. A cross-coupled
displacement satisfies the binary symmetry and nilpotence conditions
but fails the full cubic tensor. These are algebra and native-defect
checks, not a completed one-clean Boolean source. The uniform theorem
is the analytic argument above.
