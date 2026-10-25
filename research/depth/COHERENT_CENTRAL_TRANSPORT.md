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
