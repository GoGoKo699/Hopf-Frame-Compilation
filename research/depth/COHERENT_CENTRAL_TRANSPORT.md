# Coherent central transport for completed parity shears

[Coordinate-weighted pulses](COORDINATE_WEIGHTED_PULSE_LIMITS.md) ·
[Clifford-interleaved lookup](CLIFFORD_INTERLEAVED_LOOKUP.md) ·
[Dirty cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md) ·
[Complete derivative contract](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence)

This note concerns **monomial completed parity shears**, not the fractional
quarter-angle pulses. Arbitrarily many fixed physical Clifford mask frames are
permitted between the completed shears. The conclusion counts those specified
completed shears; it does not lower-bound unrestricted native T-depth or exclude
joint synthesis in a different representation. The low-degree factor is paid below;
the coherent central transport remains a sufficient source interface, not an
implemented lookup source. The complete-indicator and general frame frontiers
remain unchanged.

## 1. The exact binary frame group

Fix p at least two. Work over the binary field on the m-dimensional coefficient
space, with $`m=2^p`$ and coordinates indexed by subsets S of the
p labels. Define

```math
(N_i a)_S=\begin{cases}a_{S\cup\{i\}},&i\notin S,\\0,&i\in S,\end{cases}
\qquad N_u=\sum_i u_iN_i,\qquad N_S=\prod_{i\in S}N_i,
\quad N_\varnothing=I.
```

The $`N_i`$ commute, square to zero, and their distinct square-free monomials are
linearly independent: the root row of the monomial $`N_S`$ is the coordinate row
at S. In particular $`N_u^2=0`$. Let $`P_E`$ and $`P_O`$ be the projections
onto even and odd subset labels and put

```math
E_u=P_E N_uP_O,\qquad O_u=P_O N_uP_E,\qquad
 e_u=I+E_u,\qquad o_u=I+O_u.
```

These are literal phase-free reversible binary maps. Each is an involution;
same-parity maps add their masks. Also

```math
e_uo_u=I+N_u,\qquad
[e_u,o_v]=I+N_uN_v.                                      \qquad(1)
```

Here the commutator is $`e_uo_ve_uo_v`$. The first identity follows from
$`N_u^2=0`$. For the second, $`E_uO_v+O_vE_u=N_uN_v`$,
while $`E_uO_vE_u`$, $`O_vE_uO_v`$ and $`E_uO_vE_uO_v`$ vanish. Every polynomial
of even total degree in the $`N_i`$ commutes with both parity projections and with
all $`e_u`$ and $`o_v`$, so the commutator in (1) is central in this frame group.

For an alternating binary $`p\times p`$ matrix K, meaning a symmetric
matrix with zero diagonal, define

```math
F_K=\prod_{i\lt j}(I+K_{ij}N_iN_j),\qquad
(u\wedge v)_{ij}=u_iv_j+v_iu_j.                          \qquad(2)
```

The factors commute and are involutions, hence $`F_KF_L=F_{K+L}`$ and
$`F_K^{-1}=F_K`$. The coefficient of $`N_S`$ in $`F_K`$ is zero for odd
cardinality S and is the principal Pfaffian $`\mathrm{Pf}(K_S)`$ for even cardinality,
with the empty Pfaffian equal to one. This assertion follows by expanding the
product: a surviving choice of edges must be a perfect matching of S.

In particular,

```math
F_{u\wedge v}=I+N_uN_v.                                 \qquad(3)
```

An elementary proof also covers characteristic two explicitly. In the
coefficient on a set of 2k distinct labels, expand every matching factor
$`u_iv_j+v_iu_j`$. For each partition of those labels into k positions
carrying u and k carrying v, the same monomial occurs $`k!`$ times, once for
each bijection between the two parts. For k at least two this multiplicity
is even, so the coefficient vanishes. At k equal to one it is precisely
the displayed alternating-matrix entry. Thus no characteristic-zero
exponential formula or division by factorials is being used.

Every word of fixed mask frames has the normal form

```math
M=F_K e_a o_b.                                           \qquad(4)
```

To obtain it, commute every o factor past every e factor using (1)--(3),
collect the central terms, and add masks of equal parity. For p at least two
the parameters are unique: the degree-one coefficients in root and singleton
rows identify a and b, and the remaining degree-two root coefficients identify K.
Only existence is needed below. The high-degree Pfaffian terms remain in $`F_K`$;
centrality does not mean that its action has degree at most two in the $`N_i`$.

## 2. The complete derivative and its right-commutant freedom

Allow an arbitrary phase multiplying the monomial source:

```math
V|x,c\rangle=e^{i\phi(x,c)}|x,M(x)c\rangle,
\qquad M(x)=F_{K(x)}e_{a(x)}o_{b(x)}.                    \qquad(5)
```

The phase may include the physical S and address-to-bank CZ phases of the
Clifford mask frames. All address inputs and coefficient inputs remain
arbitrary. Because the probe is diagonal, these phases cancel in its actual
conjugation:

```math
VZ_{A_\varnothing}V^\dagger|x,c\rangle
 =(-1)^{(M(x)^{-1}c)_\varnothing}|x,c\rangle.
```

The root row of $`o_b`$ is the root row of the identity, while
$`M^{-1}=o_be_aF_K`$. Consequently

```math
\mathrm{row}_{\varnothing} M^{-1}
 =\mathrm{row}_{\varnothing}F_K(I+N_a).             \qquad(6)
```

The degree-one coefficients of this row are $`a_i`$ and its degree-two
coefficients are $`K_{ij}`$. Define the alternating clique matrix

```math
\mathcal K(x)_{ij}=x_ix_j\quad(i\ne j),\qquad
\mathcal K(x)_{ii}=0.
```

The complete coefficient-lookup derivative

```math
Q|x,c\rangle=(-1)^{\sum_Sc_S\prod_{i\in S}x_i}|x,c\rangle
```

therefore holds if and only if

```math
a(x)=x,\qquad K(x)=\mathcal K(x).                        \qquad(7)
```

Necessity follows already from singleton and pair bank inputs. For
sufficiency, the Pfaffian of the clique matrix on a set of 2k active labels
is the parity of its $`(2k-1)!!`$ perfect matchings, hence one. Thus $`F_K`$ contains
exactly every even monomial in the active $`N_i`$. Multiplication by $`I+N_x`$
adds every odd monomial once: an odd set has an odd number of possible
removed labels. This proves (6)--(7) on every coefficient input.

The mask $`b(x)`$ is unconstrained. It is the exact right-commutant freedom
visible in this group, since $`o_b`$ fixes the root coordinate and commutes
with the probe. Requiring the entire zeta permutation would unnecessarily
remove this allowed freedom.

## 3. Rank bound with fixed mask-frame interlayers

Consider a source built from s completed parity shears $`e_{v_j(x)}`$ or
$`o_{v_j(x)}`$, separated by arbitrary fixed mask-frame Cliffords. A fixed frame
may contain any number of the physical monomial Clifford pulses whose
binary bank action is generated by constant-mask $`e_u`$ and $`o_u`$. Their S/Z
phases are retained; only their bank permutation is independent of x.
The masks $`v_j`$ may be arbitrary Boolean functions of the whole address in
this algebraic statement. Implementing those functions is **not** supplied
for free by the result.

If this source supplies the complete derivative, then

```math
s\ge\lfloor p/2\rfloor.                                 \qquad(8)
```

Normalize at the level of binary bank maps; all diagonal phases cancel
in the derivative. Replace every $`v_j(x)`$ by $`v_j(x)+v_j(0)`$, absorbing its
constant same-parity shear into an adjacent fixed frame. The normalized
masks vanish at x equal to zero. Move all fixed frames to the right; what
remains in the binary bank action is a product of s fixed-frame conjugates
of the normalized shears followed by one fixed frame C. At x equal to zero
the product of normalized factors has identity binary bank action. By (7),
the binary action of C has $`a=0`$ and $`K=0`$, hence it is $`o_b`$ for
some fixed b. Its entire physical monomial unitary commutes with the
root Z probe: $`o_b`$ fixes that bit and all its diagonal phases commute
with it. Right multiplication by its actual inverse therefore removes
C without changing the complete derivative.

For a fixed frame with binary action $`F_Le_ao_b`$, centrality and (1) give

```math
C e_v C^{-1}=F_{v\wedge b}e_v,\qquad
C o_v C^{-1}=F_{a\wedge v}o_v,                            \qquad(9)
```

as equalities of binary bank maps. Physical conjugation may also produce
phases, which cancel in the derivative as in (5). Thus the j-th normalized
factor contributes a central term $`v_j(x)\wedge c_j`$ with a fixed vector $`c_j`$.
Reordering all e factors before all o factors adds only terms
$`v_i(x)\wedge v_j(x)`$. For some set of pairs determined by their ordering,

```math
K(x)=\sum_j v_j(x)\wedge c_j
     +\sum_{(i,j)\in\mathcal J}v_i(x)\wedge v_j(x).
```

The image of every summand lies in the span of its two displayed vectors.
Therefore

```math
\mathrm{im}K(x)\subseteq
\mathrm{span}\{v_1(x),\ldots,v_s(x),c_1,\ldots,c_s\},
\qquad \mathrm{rank}K(x)\le2s.                    \qquad(10)
```

At the all-one address, the target matrix is $`J+I`$. Its rank over the binary
field is p when p is even and p-1 when p is odd. Indeed $`(J+I)y=0`$
says y is the constant vector whose entry is the sum of its entries;
that nonzero constant vector is in the kernel exactly when p is odd.
Combining (7) and (10) proves (8).

If (8) fails, the complete derivative error is exactly two: some singleton
or pair coefficient input has the wrong sign. If a native implementation
approximates the whole stated source, tensored with identity on arbitrary
work, within epsilon in operator norm, its derivative error is at least
$`2-2\epsilon`$. The latter includes return errors and arbitrary references.

## 4. A paid low-degree factor and the remaining source interface

For the coordinate mask $`v(x)=x`$, a completed parity shear is a commuting
CCX family on $`r=pm/4`$ edges. Hadamards on its tail bank turn it into a
cubic CCZ family. The existing dirty cubic compiler pays one such half
with at most 8 T bands, $`8r`$ T gates and r returned dirty helpers; the
Hadamard walls are Clifford. This is a paid component, not an assumption
that arbitrary masks $`v(x)`$ are cheaply computable.

In particular the complete phase-free map

```math
\mathsf J|x,c,w\rangle=|x,(I+N_x)c,w\rangle,
\qquad I+N_x=e_xo_x,
```

has common forward and actual-inverse bounds

```math
\begin{aligned}
D_T(\mathsf J)&\le16,&T(\mathsf J)&\le4pm,\\
G_{\rm all}(\mathsf J)&\le8pm+2m,
&w_{\rm dirty}(\mathsf J)&\le pm/4,\\
a_{\rm clean}&=0,&\epsilon&=0.
\end{aligned}                                                   \qquad(11)
```

Each of the two parity classes has r edges. The cubic-only eight-step
Gray word costs eight T gates and eight CNOTs per edge. Its two Hadamard
walls use m gates, so both classes cost 16r T gates, 16r CNOTs and 2m
Hadamards. The same r private dirty targets are reused after the first
class returns them. Equation (1) proves the complete binary action;
the Gray-word identity gives its literal phase-free action on every
logical and dirty input. Reversing that native word gives the actual
inverse with the same counts, without a fresh initialization.

The missing full even transport has the precise sufficient interface

```math
\mathsf F|x,c,w\rangle=|x,F_{\mathcal K(x)}c,w\rangle,
\qquad
F_{\mathcal K(x)}=\sum_{|S|\ \mathrm{even}}
       \left(\prod_{i\in S}x_i\right)N_S.
                                                               \qquad(12)
```

Suppose a native unitary $`\widetilde{\mathsf F}`$ approximates this
entire map within operator error delta, with arbitrary dirty work
returned in that same bound. Suppose its forward word and actual inverse
have common bounds d, t, g and w for native T-depth, T-count, all native
gates and extra dirty width, with zero clean inputs. Then

```math
\begin{aligned}
\widetilde V&=\mathsf J\widetilde{\mathsf F},&
V_0&=\mathsf J\mathsf F,\\
\|\widetilde V-V_0\|&\le\delta,&
\|\widetilde VZ_{A_\varnothing}\widetilde V^\dagger-Q\|
 &\le2\delta,\\
D_T(\widetilde V)&\le d+16,&
T(\widetilde V)&\le t+4pm,\\
G_{\rm all}(\widetilde V)&\le g+8pm+2m.
\end{aligned}                                                   \qquad(13)
```

Identity factors on all helpers are implicit in the norm comparison and
the target. The proof is multiplication by the exact unitary $`\mathsf J`$
and the standard two-term conjugation bound. Equations (6)--(7) give
$`V_0Z_{A_\varnothing}V_0^\dagger=Q`$. With reusable returned helpers,
the completion needs at most $`\max(w,pm/4)`$ extra dirty wires. If the
candidate's private pool is kept separate, the bound is $`w+pm/4`$.
The estimate holds coherently with every external reference, and the
actual inverse is the reversed composition
$`\widetilde{\mathsf F}^{\dagger}\mathsf J^{\dagger}`$ with the same
common bounds. It does not assume that an approximate source has an
exactly clean intermediate output.

Thus paying the full map (12) would be sufficient for the complete
coefficient-lookup derivative, with the stated additional costs. This
is not asserted to be equivalent to an unrestricted derivative source:
the derivative leaves the b(x) right-commutant freedom in (7), and a
different source can have further freedom outside this frame group.
The fixed-address Clifford $`F_{\mathcal K(\mathbf1)}`$ supplies no
coherent implementation of (12).

Equation (8) explains why repeatedly inserting these completed monomial
halves into the same fixed mask-frame group does not give a bounded-size
source. It remains a **specified-source count**. It does not establish a
native T-depth lower bound: the word might admit a joint synthesis, and
fractional nonlinear pulses, arbitrary Clifford bank transformations outside
this mask-frame group, and address motion lie outside the proof.

The [exact fixtures](../../tests/test_clifford_interleaved_lookup.py)
corroborate the binary commutator, Pfaffian, derivative and rank identities.
The all-size claims and native ledger follow from the displayed algebra
and the returned-dirty Gray word, rather than from the finite fixtures.

<a id="partial-character-diagonalizer"></a>

## 5. A partial character transform with arbitrary residual fibers

The full returned map (12) is stronger than necessary. An additional arbitrary
bank lets a source move every coefficient sector while a final probe echo
returns that bank. This gives a different sufficient interface: diagonalize
the central characters, with no prescribed coordinates on their residual fibers.
The character-dependent middle is a paid cubic phase; the fixed transform
itself remains unpaid at sublogarithmic native depth.

### Embed every coefficient input in a formal unit

Work in the square-zero algebra on the p labels and an additional generator
$`e_0`$. Complement the coefficient labels by writing

```math
\widehat c=\sum_{S\subseteq[p]}c_{[p]\setminus S}e_S,
\qquad u=1+d+e_0\widehat c.
```

All m coefficients of c and all $`m-1`$ nonconstant coefficients of d are
arbitrary physical inputs. The displayed 1 is a formal algebraic constant,
with no physical wire or initialized-input assumption. Multiplication by
$`e_i`$ on $`\widehat c`$ is precisely the lowering map $`N_i`$ on c.
For alternating K, put

```math
g_K=\prod_{i\lt j}(1+K_{ij}e_ie_j),\qquad
\mathcal A_K:\quad
(\widehat c,d)\longmapsto
\bigl(g_K\widehat c,\ g_K(1+d)-1\bigr).                 \qquad(14)
```

This is a phase-free Clifford involution on all $`2m-1`$ coefficient wires.
Its c action is $`F_K`$. Its affine d action is retained, rather than claimed
to return inside the source. The group law is
$`\mathcal A_K\mathcal A_L=\mathcal A_{K+L}`$.

Let $`q=\binom p2`$ and let s be the q pair coefficients $`d_{\{i,j\}}`$.
Under (14), $`s\mapsto s+K`$. Define $`r=g_su`$; its pair coefficients
are zero. Omitting those formal zero coordinates, the map

```math
\mathcal S:\quad u\longmapsto
\bigl(s,\ r\text{ with its pair coordinates omitted}\bigr)
```

is a permutation of all $`2m-1`$ physical bits. Its inverse is $`u=g_sr`$.
There are $`2m-1-q`$ arbitrary residual bits; no extra zero-filled bank is
introduced. Since $`g_{s+K}g_K=g_s`$, one has

```math
\mathcal S\mathcal A_K\mathcal S^\dagger
 =X_s^K\otimes I_{\rm residual},\qquad
B=H_s\mathcal S,\qquad
B\mathcal A_KB^\dagger=Z_s^K\otimes I_{\rm residual}.  \qquad(15)
```

Thus a common character transform exists. Arbitrary unitaries inside the
joint character fibers give other valid transforms; (15) does not prescribe
their residual coordinates. Another sufficient choice is
$`B=H_{\rm pairs}F_{p+1}^\dagger`$ using the
[HRT source](FULL_BANK_CHARACTER_TRANSPORT.md): multiplication by $`g_K`$
adds K only to its additive pair coordinates. Its known compiler retains
logarithmic depth. Requiring that entire HRT coordinate change is unnecessary
for the interface (15).

### A paid coherent middle and the complete source ledger

The fixed B acts on bank and work, independently of the original address x.
For the clique matrix $`\mathcal K(x)`$, its coherent diagonal middle is

```math
D=\prod_{i\lt j}\mathrm{CCZ}(x_i,x_j,s_{ij}),\qquad
\mathcal A(x)=B^\dagger D B.                            \qquad(16)
```

The cubic-only dirty compiler pays D with at most eight T bands,
$`8q`$ T gates, $`8q`$ CNOTs and q returned arbitrary helpers. Its total
native gate count is therefore $`16q`$, with zero error and no clean input.
All shared address and character controls are included.

Let $`B_0`$ be a specified ideal transform satisfying (15), including identity
on any private work in that conjugation identity. B itself need not return
that work between calls. Suppose a literal native $`\widetilde B`$ has
full operator error at most eta relative to $`B_0`$, common forward and
actual-inverse bounds $`d_B,t_B,g_B`$, and $`h_B`$ private dirty wires.
Here $`g_B`$ counts every elementary native gate, including T gates.
Combine (16) with the paid $`\mathsf J`$ from (11):

```math
\widetilde V=\mathsf J\widetilde B^\dagger D\widetilde B.
```

In the ideal source, the c action is $`(I+N_x)F_{\mathcal K(x)}`$ and is
independent of d. Consequently the root-probe conjugation is exactly Q,
tensored with identity on d and every private helper. The actual adjoint
of the source cancels the arbitrary intermediate d action. Equations
(6)--(7) and unitary telescoping give the complete bounds

```math
\begin{aligned}
D_T(\widetilde V)&\le2d_B+24,&
T(\widetilde V)&\le2t_B+8q+4pm,\\
G_{\rm all}(\widetilde V)&\le2g_B+16q+8pm+2m,&
w_{\rm extra}&\le m-1+h_B+\max(q,pm/4),\\
\|\widetilde VZ_{A_\varnothing}\widetilde V^\dagger
        -Q\otimes I_{\rm work}\|&\le4\eta.
\end{aligned}                                                   \qquad(17)
```

The conservative width keeps B's private pool separate and reuses the
returned cubic-middle and J pools. Error (17) includes every address, bank,
work and external-reference input, including all final return errors.
The actual reversed source has the same common resource bounds. The
[two-call derivative readout](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence)
then gives the complete indicator with twice the source depth and T-count,
the same derivative error, and its already paid Clifford walls.

An exact simultaneous identity on pair generators proves (15) for every K.
For an approximate transform, separate generator errors must be accumulated
over the up to q factors or replaced by a direct uniform bound on the whole
coherent middle. A per-generator tolerance is not automatically the complete
source tolerance in (17).

### The forward diagonalizer requires growing hierarchy level

The partial interface does not admit a bounded Clifford-hierarchy shortcut.
If an exact B satisfying (15) belongs to hierarchy level $`\mathcal C_k`$,
then

```math
k\ge\lfloor p/2\rfloor+2.                              \qquad(18)
```

This concerns B, not its inverse. The explicit HRT choice above has
$`B^\dagger=F_{p+1}H_{\rm pairs}\in\mathcal C_3`$; the hierarchy is not
closed under taking adjoints. Neither hierarchy statement pays native depth.

For the proof, let $`\mathcal A_{ij}`$ denote (14) with only pair ij present.
Starting from any physical Pauli R, conjugation by B puts
$`BRB^\dagger`$ in $`\mathcal C_{k-1}`$. Each group commutator with a
Pauli character image lowers the hierarchy level by one; after $`k-1`$
commutators the result is scalar. Conjugating back, the corresponding nested
commutators of R with $`\mathcal A_{ij}`$ must therefore be scalar.

For a Clifford with binary symplectic matrix $`S_{ij}`$, the Pauli label
of a commutator is multiplied by $`S_{ij}-I`$. On c, this increment is
$`N_iN_j`$ on X labels. Choose R to flip $`c_{[p]}`$ in the original
coefficient ordering. Products for any sequence of disjoint pairs remain
nonzero through $`\lfloor p/2\rfloor`$ factors: they move that label to
the coordinate with those labels removed. The affine d translation and
identity private work do not cancel this c block. This proves (18).

There is a uniform third-level error floor as well. For p at least four,
any $`B\in\mathcal C_3`$ and Hermitian Pauli targets for two disjoint pairs
obey

```math
\|B\mathcal A_{ij}B^\dagger-P_{ij}\|\le\epsilon_{ij},\quad
\|B\mathcal A_{kl}B^\dagger-P_{kl}\|\le\epsilon_{kl}
\quad\Longrightarrow\quad
4\epsilon_{ij}+2\epsilon_{kl}\ge2.                     \qquad(19)
```

Indeed, put $`\widetilde A_{ab}=B^\dagger P_{ab}B`$. The double commutator
$`[[R,\widetilde A_{ij}],\widetilde A_{kl}]`$ is scalar $`\pm I`$,
because $`BRB^\dagger`$ is Clifford. For the c witness above, its ideal
counterpart is a nonidentity Hermitian Pauli, at distance two from either scalar.
Two telescoping estimates bound their distance by the left side of (19).
Thus the uniform pair-image error is at least one third. Arbitrary private
work and reference inputs are included.

These are limits on a common character diagonalizer, not growing lower
bounds on unrestricted native T-depth or complete-indicator error. A circuit
outside $`\mathcal C_3`$ may still have constant native depth; joint synthesis
of $`B^\dagger D B`$ need not implement B separately. The remaining question
is whether nonmonomial native dynamics can pay this partial character
transform, or its whole coherent conjugation, within the complete
sublogarithmic-depth and polynomial-overhead budgets.

The [complete observable-basis fixtures](../../tests/test_complete_observable_bases.py)
check the free pair action, arbitrary residual sectors and disjoint-pair
commutators. The all-size identities, hierarchy argument and native/error
ledger are the analytic claims above.

<a id="clique-polarization-integration"></a>

## 6. Cheap polarization and its sharp coherent-query cost

Write $`F(x)=F_{\mathcal K(x)}`$ and $`J(x)=I+N_x`$ for their binary
bank actions. The quadratic clique map obeys

```math
\mathcal K(x+y)+\mathcal K(x)+\mathcal K(y)=x\wedge y,
\qquad
F(x+y)F(x)F(y)=I+N_xN_y=[e_x,o_y].                 \qquad(20)
```

Thus its polarization is cheap even though F itself remains unpaid at
sublogarithmic native depth. The same cocycle appears in
$`J(x)J(y)J(x+y)=I+N_xN_y`$. Their product removes it:

```math
\rho_x=J(x)F(x)=\prod_i(I+x_iN_i),\qquad
\rho_x\rho_y=\rho_{x+y}.                              \qquad(21)
```

These identities follow from (1)--(3), or by expanding the even and odd
subset products. Equation (21) is the complete zeta action, so replacing F
by the full zeta is not a cheaper construction.

For $`p\ge2`$, four paid parity halves compile the polarization on arbitrary
coefficient inputs, with arbitrary coherent x,y, zero clean work and

```math
D_T\le32,\qquad T\le8pm,\qquad
G_{\rm all}\le16pm+4m,\qquad w_{\rm dirty}\le pm/4.     \qquad(22)
```

The same returned pool is reused after each half; reversing the four actual
native words gives the same ledger. At p at most one the polarization is
identity. Cheap polarization does not by itself integrate the quadratic map.

### Fixed-address polarization loses the quadratic refinement

On central character s, write
$`\chi_s(K)=(-1)^{\sum_{i\lt j}s_{ij}K_{ij}}`$ and let S be the symmetric
zero-diagonal matrix with upper triangle s. Then

```math
\chi_s(\mathcal K(x))=(-1)^{q_s(x)},\qquad
q_s(x)=\sum_{i\lt j}s_{ij}x_ix_j,\qquad
\chi_s(x\wedge y)=(-1)^{x^{\mathsf T}Sy}.              \qquad(23)
```

All characters occur: in the complemented algebra ordering of Section 5,
the orbit of the unit 1 under $`F_K`$ is free, since the pair coefficients
of $`g_K`$ are K. Its Hilbert span is the regular representation of the
central group. This is a subspace for a lower-bound witness, not a supplied
initialized resource for a construction.

If x has at least three active coordinates, take S to be a triangle on
three of them. Then $`Sx=0`$ but $`q_s(x)=1`$. The trivial character and this
character therefore agree on every query $`F_{x\wedge y}`$, for every y,
while the required phase has opposite signs. Arbitrary coherent y queries
and character-independent interlayers cannot distinguish them, regardless of
query count. The resulting direct-unitary error is at least square root two;
for the half-projector Gram it is at least one quarter. Varying the first
argument or mixing central character sectors leaves this restricted screen.

### A precise quantum-query model

Allow a query with coherent controls u,v whose character fiber is

```math
Q_s|u,v,w\rangle=\chi_s(u\wedge v)|u,v,w\rangle.        \qquad(24)
```

Interlayers, preparation and acceptance may depend on the original address x
and act arbitrarily on quantum work, but are independent of s. Clean or dirty
work and arbitrary accepted output fibers are allowed. No physical cost for
these interlayers is asserted. The oracle is (24), not an arbitrary
character-dependent gate. A fixed unconditional central Clifford $`F_{K_0}`$
may also be included; all such scalar characters commute and accumulate
into one offset K_0 independent of x.

After t queries, every resulting matrix or accepted block has expansion

```math
A_s=\chi_s(K_0)\sum_{\mathrm{rank}R\le2t}\chi_s(R)A_R. \qquad(25)
```

Each query contributes only rank-at-most-two character labels. Matrix
multiplication convolves labels, rank is subadditive, and character-independent
interlayers or compression preserve their support. This includes coherent
Fourier sums over query labels; it does not assume a classical query schedule.

For a full unitary or initialized isometry approximating
$`\chi_s(\mathcal K(x))I`$, if
$`\mathrm{rank}(\mathcal K(x)+K_0)>2t`$, the target Fourier coefficient in
(25) vanishes. Uniformly average normalized squared Hilbert--Schmidt error
over s. Both norm terms are one and the cross term is zero, giving average
two. Some character sector therefore has operator/isometry error at least
$`\sqrt2`$.

For an accepted block with target Gram
$`(I+\chi_s(\mathcal K(x))I)/4`$, the offset cancels from $`A_s^\dagger A_s`$,
whose Fourier support has rank at most 4t. If
$`\mathrm{rank}\mathcal K(x)>4t`$, its target Fourier coefficient is I/4
while the actual one is zero. The normalized trace and
$`|\mathrm{Tr}(E)/d|\le\|E\|`$ give complete Gram error at least 1/4.
No condition on the accepted polar action was used.

Put $`h=\lfloor p/2\rfloor`$. The worst-case clique rank is 2h. A fixed
offset cannot reduce the worst-case direct rank: for even p, the Pfaffian
of $`K_0+\mathcal K(x)`$ is a multilinear polynomial whose top monomial
$`\prod_i x_i`$ has coefficient one. Indeed the all-clique perfect matchings
are odd in number, while using a K_0 edge omits its two vertex variables.
The polynomial is nonzero on some Boolean input, giving rank p. For odd p,
apply the same argument on any p-1 principal coordinates. Consequently

```math
t\ge h\quad\text{for direct error below }\sqrt2,
\qquad
t\ge\lceil h/2\rceil\quad\text{for Gram error below }1/4. \qquad(26)
```

### Both query bounds are attained

An explicit fixed-mask decomposition avoids adaptive rank factorization.
For $`1\le j\le h`$, let

```math
a_j=e_{2j-1}+\sum_{k>2j}e_k,\qquad
b_j=e_{2j}+\sum_{k>2j}e_k,\qquad D_x=\mathrm{diag}(x),
\qquad
\mathcal K(x)=\sum_{j=1}^{h}(D_xa_j)\wedge(D_xb_j).     \qquad(27)
```

Each term supplies all edges from its pair to later coordinates, and no
edge within that later tail. Their sum is exactly the clique matrix.
The product of h polarization queries therefore implements F(x) exactly.

For the Gram, divide these h terms into two lists of length at most
$`t=\lceil h/2\rceil`$, defining $`K_1(x)+K_2(x)=\mathcal K(x)`$. Pad a
shorter list with zero queries. One signal selects the two lists coherently;
t calls implement $`\mathrm{SELECT}_s=\mathrm{diag}(\chi_s(K_1),\chi_s(K_2))`$.
With convention $`R_y(\theta)=\exp(-i\theta Y)`$, use

```math
U_s=R_y(-\pi/8)\,\mathrm{SELECT}_s\,R_y(3\pi/8),
\qquad
\langle0|U_s|0\rangle
 =\frac{\chi_s(K_1)+\chi_s(K_2)}{2\sqrt2}.             \qquad(28)
```

Its Gram is exactly $`(I+\chi_s(\mathcal K(x))I)/4`$. The accepted sign
can depend on the character; recovery is allowed to retain that polar
freedom. Only one acceptance signal is initialized. The unequal boundary
angles are essential for the factor one half in the projector Gram.

These boundaries have exact native words. Put C=SH. Then
$`CTC^\dagger=e^{i\pi/8}R_y(\pi/8)`$ and $`HZ=R_y(\pi/4)`$.
Use $`(HZ)CTC^\dagger`$ for PREP and $`CT^\dagger C^\dagger`$ for POST;
their opposite scalar phases cancel literally in the complete word.
They add two T-type gates, with their actual inverses. The oracle upper
bound does not supply a free native branch-controlled query or label
workspace; those operations still require their own full allocation proof.

Equations (26)--(28) determine this coherent-query model exactly. They do
not bound unrestricted native T-depth. There is also a useful extension
which permits some physical mixing of central characters.

### Orbit phases still leave a rank-distance boundary

Use the regular-orbit basis $`|g_K\rangle`$. A polarization query translates
K by $`u\wedge v`$, of rank at most two. Permit additional interlayers
of the form $`\sum_K|g_K\rangle\langle g_K|\otimes U_K`$, with arbitrary
unitaries on work, preparation and acceptance on that work. These can depend
on K, unlike the preceding character-independent interlayers. After t queries,
an input at K reaches only $`K+K_0+R`$ with $`\mathrm{rank}R\le2t`$.
Fixed unconditional central offsets remain allowed and accumulate into K_0.
The desired direct output is orthogonal to that support whenever
$`\mathrm{rank}(\mathcal K(x)+K_0)>2t`$, again giving error square root two.
Two accepted paths from input coordinates differing by rank greater than 4t
have disjoint final-coordinate support. Their Gram matrix element is zero,
whereas the target between $`|g_0\rangle`$ and $`|g_{\mathcal K(x)}\rangle`$
is 1/4. Thus (26) remains sharp in this larger model as well.

For example, at p=4 a single physical Z on the top coefficient restricts to

```math
|g_K\rangle\longmapsto
(-1)^{K_{12}K_{34}+K_{13}K_{24}+K_{14}K_{23}}|g_K\rangle. \qquad(29)
```

Its central-character matrix has every entry of magnitude 1/8, since the
three complementary pairs give three independent quadratic Fourier sums.
It therefore mixes even character differences of rank four at zero T-depth,
leaving the first model's premise. Nevertheless it is an orbit phase and
cannot overcome the extended rank-distance bound. This orbit calculation
does not supply a source on all arbitrary coefficient inputs. A successful
changed word must leave the relevant model through other bank dynamics or
joint native synthesis, rather than add these phases alone.

The Gram target here is the half-projector of the central F action, not the
complete lookup projector. The weaker root-probe interface also need not
realize F as a separate word; its separate readout remains essential.
The [polarization fixtures](../../tests/test_central_polarization.py) check
the algebra, ranks, shifted cliques, orbit phase and literal one-signal boundaries.

<a id="weak-pair-gauge-source"></a>

## 7. A weak pair gauge removes the full-basis correction

The complete carrier action needs less than the full affine action (14).
Retain the dummy-generator embedding of Section 5, with $`p\ge2`$,
$`m=2^p`$, $`q=\binom p2`$, and all $`2m-1`$ nonconstant inputs arbitrary:

```math
y=d+e_0\widehat c,\qquad s=(d_{ij})_{1\le i\lt j\le p},
\qquad g_s=\prod_{i\lt j}(1+s_{ij}e_ie_j).
```

Let C be the product of the local $`C_{\{i,j\}}`$ matchings from the
[local HRT compiler](COHERENT_FEEDBACK_SOURCES.md#1-an-exact-local-change-of-hrt-coordinates),
over pairs excluding the dummy generator. Its exact action is

```math
C:y\longmapsto g_s y,\qquad C^2=I.                    \qquad(30)
```

There is no constant coefficient in y. Each matching therefore preserves
every selected pair coordinate; all its targets have degree at least three.
Consequently their controls stay fixed, the completed matchings commute,
and their product is multiplication by $`g_s`$. Since $`g_s^2=1`$, C is an
involution. This algebra does not parallelize the overlapping native gates.

For an alternating K, put $`k_2=\sum_{i\lt j}K_{ij}e_ie_j`$ and let
$`Q_K:y\mapsto y+k_2`$ toggle its pair coordinates. Rightmost-first
composition gives the full physical action

```math
C Q_K C(y)=g_{s+K}(g_s y+k_2)
          =g_K y+g_{s+K}k_2.                           \qquad(31)
```

The additive term contains no $`e_0`$. Hence
$`\widehat c\mapsto g_K\widehat c`$ exactly, independently of every
arbitrary d input. The other coefficients may move. With
$`K=\mathcal K(x)`$, the source $`V=\mathsf J C Q_K C`$ therefore has the
required complete root-probe derivative by (6)--(7). The actual source
inverse cancels the arbitrary intermediate work motion in that derivative.
No d input is initialized or discarded.

### Conditional native and error ledger

The pair toggle $`Q_{\mathcal K(x)}`$ is the shared-address batch
$`\prod_{i\lt j}\mathrm{CCX}(x_i,x_j;d_{ij})`$. H on its q targets turns it
into q cubic phases. The eight-band dirty compiler gives
$`T\le8q`$, $`D_T\le8`$, $`G_{\rm all}\le18q`$ and q returned private
dirty helpers, including both H walls. No clean input is used.

Suppose an actual native $`\widetilde C`$ approximates $`C\otimes I`$
in full operator norm by eta, with common forward/actual-inverse bounds
$`t_C,d_C,g_C,h_C`$. Use
$`\widetilde V=\mathsf J\widetilde C Q_K\widetilde C^\dagger`$. Then

```math
\begin{aligned}
T(\widetilde V)&\le2t_C+8q+4pm,&
D_T(\widetilde V)&\le2d_C+24,\\
G_{\rm all}(\widetilde V)&\le2g_C+18q+8pm+2m,&
w_{\rm extra}&\le m-1+h_C+\max(q,pm/4),\\
\|\widetilde V-V\|&\le2\eta,&
\|\widetilde VZ_{A_\varnothing}\widetilde V^\dagger-Q\otimes I\|
 &\le4\eta.
\end{aligned}                                                   \qquad(32)
```

Here $`w_{\rm extra}`$ excludes x and the principal c bank. It includes
all $`m-1`$ arbitrary d wires and keeps C's private pool separate from the
returned Q/J pool. Every actual inverse has the same common bounds.
Unitary telescoping includes intermediate leakage and all references.
The usual two-call probe readout doubles the source T-count and T-depth
and gives the complete indicator; it does not supply C. A direct
implementation of the entire conjugation with only its
probe guarantee can be weaker than this sufficient full-C approximation.

### Which simultaneous-Hadamard correction can be omitted

For the selected pair family, let $`E=F_{p+1}H_{\rm pairs}F_{p+1}^\dagger`$
and let B toggle each $`y_{ij}`$ by $`y_i y_j`$. Regrouping the individual
encoded-Hadamard palindromes into $`C B H_{\rm pairs}B C`$ does not give E.
The [intersecting-family audit](COHERENT_FEEDBACK_SOURCES.md#pair-star-hadamard-fusion)
proves the first disjoint-pair witnesses.
In particular, the six pairs on four labels give squared zero-column error
7/8: the true output contains a quartic Pfaffian that the regrouped word lacks.

The full simultaneous encoding uses
$`T(y)=g_s(1+y+s_2)-1`$, with $`s_2=\sum_{i\lt j}s_{ij}e_ie_j`$.
It preserves all pair and singleton coordinates and has actual inverse
$`T^{-1}(y)=g_s(1+y)-1+s_2`$; it need not be an involution. For a selected
pair S, direct substitution gives
$`T X_S T^{-1}(y)=(1+e_S)(1+y)-1`$, the action of $`FX_SF^\dagger`$.
B commutes with $`X_S`$, while its Z image is $`Z_SD_{q_S}`$, unchanged
by T. Thus TB has F's selected-pair Pauli images, proving
$`E=T B H_{\rm pairs}B T^{-1}`$ with that actual inverse.
Relative to C, T's additional translation is $`g_s(1+s_2)-1`$.
It has no dummy-generator component. Thus repairing that full encoding is
unnecessary for (31); its failed full-Hadamard column does not refute the
weaker root source. Exact star batches do have a paid 24-band encoded-Hadamard
schedule, but parallel execution of different stars is not supplied.

The alternative E sandwich also has a paid middle. The inverse HRT pair
coordinate is $`k_{ij}(y)=y_{ij}+y_i y_j`$, so

```math
\mathcal A_{\mathcal K(x)}=E\widetilde D(x)E,\qquad
\widetilde D(x)=(-1)^{\sum_{i\lt j}x_ix_j(y_{ij}+y_i y_j)}. \qquad(33)
```

Allocate q arbitrary outer bits $`h_{ij}`$ and a separate reusable q-bit
Gray pool. Set
$`U=\prod_{i\lt j}\mathrm{CCX}(x_i,x_j;h_{ij})`$ and
$`D=\prod_{i\lt j}\mathrm{CCZ}(h_{ij},y_i,y_j)\mathrm{CZ}(h_{ij},y_{ij})`$.
Chronological $`U,D,U^\dagger,D^\dagger`$ gives the phase in (33), since
the two D exponents differ by $`x_ix_j(y_i y_j+y_{ij})`$. Every helper
returns. The cubic-only compiler prices U at 8q T and 18q total gates,
and D at 8q T and 19q total gates, including $`\mathrm{CZ}=H\mathrm{CX}H`$.
Thus the complete middle and its actual inverse have

```math
T\le32q,\qquad D_T\le32,\qquad G_{\rm all}\le74q,
\qquad w\le2q,\qquad a_{\rm clean}=0,\qquad\epsilon=0. \qquad(34)
```

All shared controls and literal phases are included. This pays the middle
of this particular sandwich; it does not pay E or its joint use.

### A cheap truncation and a nonmonomial trial still fail

The involution $`C_0(y)=(1+s_2)y`$ preserves all pair controls. Multiplication
by the homogeneous degree-two form $`s_2`$ flips
$`\lfloor|S|/2\rfloor\bmod2`$ on coefficient labels and squares to zero.
Splitting into the two directed XOR shears therefore implements C_0 in
sixteen T bands. With $`R=q(2^{p-1}-1)`$ raw monomial occurrences, valid
bounds are $`T\le8R`$, $`G_{\rm all}\le16R+4m`$ and at most R returned
dirty helpers, with zero clean inputs. Nevertheless,

```math
C_0Q_KC_0(y)=(1+k_2+s_2k_2)y+k_2+s_2k_2.             \qquad(35)
```

The carrier retains the unwanted $`s_2k_2`$ term. At p at least four, use
$`s_2=e_1e_2`$, $`k_2=e_3e_4`$ and
$`\widehat c=e_{[p]\setminus\{1,2,3,4\}}`$: its root sign differs from
the target. Even at $`s_2=0`$, a clique on the first four active labels
has a nonzero quartic Pfaffian omitted by $`1+k_2`$. Both are exact error-two
probe witnesses, unaffected by improved gate precision.

A genuinely nonmonomial variation also has a direct failure. At p=4,
restrict to the invariant sector
$`y=s_2+e_0(1+c_2+z e_{1234})`$, with pair banks ordered
$`(12,13,14,23,24,34)`$. Here
$`C_0:(s,c,z)\mapsto(s,c+s,z+\sum_{j=1}^6s_jc_{7-j})`$.
For any six-bit physical Hadamard mask a, put
$`E_a=C_0H_s^aC_0`$ and $`V_a=E_aX_s^{\mathbf1}E_a`$.
Fourier summation in the $`X_z=\pm1`$ sectors gives

```math
V_aZ_zV_a^\dagger
 =X_s^aX_c^aZ_s^aX_z^{|a|\bmod2}Z_z.                 \qquad(36)
```

On the zero vector this is $`|a,a,|a|\bmod2\rangle`$ with positive amplitude,
whereas the central clique probe gives minus the zero vector. The empty
mask has error two and all 63 nonempty masks have error at least square root
two. This is a specified p=4 source-family failure, not a restriction on
arbitrary nonmonomial dynamics.

The higher matching terms can be displayed without treating them as free.
Let $`q(v)`$ be the quadratic correction in the HRT inverse, and let
$`\gamma_j(s_2)`$ sum all products of j pairwise-disjoint selected edges.
Then $`q^{\circ j}(s_2)=\gamma_{2^j}(s_2)`$ and

```math
g_s=\prod_{j\ge0}\bigl(1+q^{\circ j}(s_2)\bigr),       \qquad(37)
```

with iterates zero beyond the available degree. For the first identity,
partitioning $`2^{j+1}`$ distinct edges into two unordered equal groups has
odd multiplicity. In the product, the multinomial coefficient for the
binary decomposition of any matching size is odd. This proves (37).
Its minimum degrees double, but computing these iterates or controlling
their products is not a free operation or a native-depth lower bound.

The first unpaid operation remains the whole $`C Q_{\mathcal K(x)}C`$
conjugation, or a different native word with the same complete root image.
The [pair-gauge fixtures](../../tests/test_pair_gauge_fusion.py) check the
full arbitrary-bank identity, probe, middle, truncation and changed-word
witnesses. No sublogarithmic source or new general frame bound follows.

<a id="reduced-carrier-gauge"></a>

## 8. Only the pair register is needed, and the probe leaves more freedom

The higher auxiliary coefficients in (30) are unnecessary for the source
interface. Return to the original downward coefficient ordering of Section 1.
Use an independent q-bit pair register s and the m-bit carrier c, both entirely
arbitrary, and define the phase-free permutation

```math
U_F|s,c\rangle=|s,F_s c\rangle,
\qquad F_s=\prod_{i\lt j}(I+s_{ij}N_iN_j).
```

The central group law gives $`U_F^2=I`$. If Q_K adds K to s, then

```math
U_F Q_K U_F|s,c\rangle=|s+K,F_Kc\rangle.                \qquad(38)
```

This proves the complete carrier action without supplying any higher d
coordinates or an initialized unit coefficient. The q pair inputs replace
the $`m-1`$ additional inputs in (32). Only the root row of the linear
completion is needed, so replace the full J by the smaller involution

```math
L_x:\ c_\varnothing\longmapsto
c_\varnothing+\sum_i x_i c_{\{i\}},
\qquad c_S\longmapsto c_S\quad(S\ne\varnothing).
```

With $`K=\mathcal K(x)`$, use $`V=U_FQ_KU_F L_x`$, with L on the right.
The root row of its inverse is the complete lookup row
$`\mathrm{row}_\varnothing L_xF_K=\sum_{S\subseteq\mathrm{supp}(x)}e_S^{\mathsf T}`$
by (6)--(7). The pair register moves inside V and returns under its actual
inverse in $`VZ_{A_\varnothing}V^\dagger`$. No promise on s, c or a
reference is used.

### Paid factors and the reduced conditional ledger

For one pair i,j, its factor in U_F is a shared-control batch

```math
\prod_{S:\,S\cap\{i,j\}=\varnothing}
 \mathrm{CCX}(s_{ij},c_{S\cup\{i,j\}};c_S).
```

The $`r=m/4`$ carrier pairs are disjoint. Hadamards on their targets and
the [four-band shared-control phase word](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
therefore give, including the target Hadamards,

```math
D_T\le4,\qquad T\le6r+(r\bmod2),\qquad
G_{\rm Cliff}\le12r+3,\qquad
G_{\rm all}\le18r+(r\bmod2)+3.                         \qquad(39)
```

There are no private helpers or clean inputs. The actual native inverse has
the same majorants; replace inverse S by three S gates and reduce each
adjacent root S power modulo four, as allowed by the literal identity $`S^4=I`$.
The literal phase is one on every basis column. Serial composition
of all q factors gives q times these bounds. Their completed permutations
commute, but the overlapping carrier wires do not give a parallel native
schedule. In particular this explicit baseline has depth at most 4q, not
sublogarithmic depth.

L is a shared-target batch of p Toffolis. Two Hadamards on the root carrier
turn it into p CCZ gates sharing that root, with otherwise disjoint control
pairs $`(x_i,c_{\{i\}})`$. The same four-band word gives common forward
and actual-inverse bounds

```math
D_T(L)\le4,\qquad T(L)\le6p+(p\bmod2),\qquad
G_{\rm Cliff}(L)\le10p+5,\qquad
G_{\rm all}(L)\le16p+(p\bmod2)+5,
```

with no private helpers or clean inputs and literal error zero. Thus the
smaller source removes both the higher auxiliary bank and the full J word;
(11) remains a useful paid full linear map in its other applications.

If an actual native $`\widetilde U_F`$ approximates $`U_F\otimes I`$
within eta in full operator norm, with common forward/actual-inverse bounds
$`t_F,d_F,g_F,h_F`$, use
$`\widetilde V=\widetilde U_FQ_K\widetilde U_F^\dagger L_x`$.
The paid Q from (32) and the smaller L give

```math
\begin{aligned}
T(\widetilde V)&\le2t_F+8q+6p+(p\bmod2),&
D_T(\widetilde V)&\le2d_F+12,\\
G_{\rm all}(\widetilde V)&\le2g_F+18q+16p+(p\bmod2)+5,&
w_{\rm extra}&\le q+\max(h_F,q),\\
\|\widetilde V-V\|&\le2\eta,&
\|\widetilde VZ_{A_\varnothing}\widetilde V^\dagger-Q\otimes I\|
 &\le4\eta.
\end{aligned}                                                   \qquad(40)
```

Here Q in the last line is the complete lookup sign, not the pair toggle.
The helper count excludes x and c and includes the q arbitrary pair bits.
The source-private pool and the q Gray bits are reused: Q is exactly identity
on its whole private pool, including any components left by an approximate
source call, so the same full-operator telescoping argument remains valid.
This lowers the sufficient source's auxiliary allocation and the paid
completion cost; it supplies no better asymptotic depth. Two actual source
calls still pay the complete indicator readout.

### Exact freedom for carrier-linear gauges

Requiring (38) for every pair translation would overprescribe the weaker
probe. Let $`H_s\in\mathrm{GL}(m,\mathbb F_2)`$ be any family, define
$`U_H|s,c\rangle=|s,H_sc\rangle`$, and put

```math
\ell_x=e_\varnothing^{\mathsf T}+\sum_i x_i e_{\{i\}}^{\mathsf T},
\qquad v_x=\sum_{S\subseteq\mathrm{supp}(x)}e_S^{\mathsf T},
\qquad V_H=U_HQ_{\mathcal K(x)}U_H^\dagger L_x.
```

The position of L matters for a general H. Its actual inverse sends the carrier
to $`L_xH_{s+\mathcal K(x)}H_s^{-1}c`$. Thus the complete root derivative
is correct if and only if

```math
\ell_x H_{s+\mathcal K(x)}=v_xH_s
\quad\text{for every }s,x.                              \qquad(41)
```

Arbitrary monomial phases on U_H do not change this criterion: they cancel
in conjugation of the diagonal probe. No phase-dependent error is discarded.
For a phase-free U_H, approximation by a full native operator gives the same
resource and telescoping bounds as (40), with the H costs substituted.

Equation (41) has a useful exact parametrization before invertibility is
imposed. Choose arbitrary row-vector-valued functions a(s) and b_i(s), and set

```math
f_s(x)=a(s+\mathcal K(x))+\sum_i x_i b_i(s+\mathcal K(x)),
\qquad
\mathrm{row}_S H_s=\sum_{T\subseteq S}f_s(\mathbf1_T).   \qquad(42)
```

All sums are binary. Boolean interpolation gives $`v_xH_s=f_s(x)`$.
Since the zero and singleton cliques vanish, the empty and singleton rows
of H_s are exactly a(s) and b_i(s); hence the other side of (41) is also
f_s(x). Conversely, every solution of (41) has this form by taking its first
p+1 rows. Every H_s must separately be invertible for U_H to be a physical
permutation. That condition is not supplied by interpolation.

The freedom is strictly larger than $`H_s=F_sH_0`$. At p=4, write
$`t=\{1,2,3,4\}`$, let e_S denote a coefficient column, and define

```math
u_s=s_{12}s_{13}e_\varnothing+s_{13}e_{12}
                   +s_{12}e_{13}+e_{123},\qquad
H_s=F_s+u_se_t^{\mathsf T}.                             \qquad(43)
```

The column u obeys (42) with $`a(s)=s_{12}s_{13}`$ and every b_i zero,
since $`f_s(x)=(s_{12}+x_1x_2)(s_{13}+x_1x_3)`$.
Moreover $`e_t^{\mathsf T}F_s^{-1}u_s=0`$: F_s fixes the top row and u has
zero top coefficient. Therefore $`I+F_s^{-1}u_se_t^{\mathsf T}`$ is an
invertible rank-one shear, proving every H_s invertible. Yet for s zero
and K supported on pair 12,
$`u_K=e_{13}+e_{123}`$ differs from $`F_Ku_0=e_3+e_{123}`$.
Thus $`H_K\ne F_KH_0`$ and full carrier transport (38) fails, while the
complete probe remains exact. This example demonstrates allowed freedom;
it is not a faster native compiler.

There is also a precise algebraic constraint on this gauge class. If every
entry of a and b_i has Boolean degree at most d in the q pair variables,
f_s has address degree at most $`2d+1`$. When $`2d+1\lt p`$, (42) gives
an identically zero top row, contradicting invertibility. Consequently

```math
d\ge\lfloor p/2\rfloor.                               \qquad(44)
```

This is a constraint on an address-independent carrier-linear gauge's low
rows. It is neither a native-depth lower bound nor an obstruction to a joint
source that mixes the address or uses nonmonomial dynamics. The reduced
[fixtures](../../tests/test_reduced_pair_gauge.py) check complete arbitrary
pair/carrier action, literal paid factors, the interpolation criterion and
its strictly weaker example. The general frame frontier remains unchanged.
