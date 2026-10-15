# Four T layers and middle crossing rank

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Three-layer obstruction](CLEAN_AWARE_DEPTH_OBSTRUCTION.md) · [Research index](../README.md)

Four-layer approximation at fixed clean budget has a quantitative structural
condition: the minimum symplectic pairing rank of the two active middle
T layers grows at least as the square root of the accuracy exponent for
the complete-frame witness below. This result, R74, gives an explicit
initialized-sector gap at each rank. The separate [operator-rigidity
theorem](FOUR_LAYER_OPERATOR_RIGIDITY.md), R75, excludes four layers
altogether in the zero-clean model at integer precision L at least three;
a five-layer lower bound with positive clean work remains open.

## 1. Circuit class and theorem

Write a native circuit with at most four T layers, padding empty layers
if necessary, as

```math
V=C_4D_4C_3D_3C_2D_2C_1D_1C_0.
```

Each D is a physical layer of T or T-dagger gates on distinct wires;
the C factors are arbitrary Clifford unitaries. The active axes of
$`D_3`$ are $`Z_j`$, $`j\in J_3`$, and those of
$`E=C_2D_2C_2^\dagger`$ are the signed Hermitian Paulis
$`P_k=C_2Z_kC_2^\dagger`$, $`k\in J_2`$. Let L and M be their
respective spans in the binary symplectic Pauli space, discarding phases
only for this linear-algebra definition. Both spans are isotropic.
Define the **middle crossing rank** of this decomposition by

```math
r=\mathrm{rank}_{\mathbb F_2}
\bigl(\langle Z_j,P_k\rangle_{\rm sp}\bigr)_{j\in J_3,k\in J_2}.
```

The pairing is zero for commuting axes and one for anticommuting axes.
Its rank is unchanged by a basis change within either span or a common
Clifford conjugation. Axis signs and T versus T-dagger do not affect r.
Different four-layer decompositions of the same circuit unitary can
have different ranks; we do not assume uniqueness.

**Restricted obstruction theorem.** For every fixed pair of integers
$`a,R\ge0`$, there is $`\varepsilon_{a,R}>0`$ such that the complete
real frame

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3),\qquad n\ge1,
```

has initialized-isometry distance at least $`\varepsilon_{a,R}`$
from every circuit **admitting** a four-layer decomposition with
$`r\le R`$, at most a clean qubits, and arbitrary dirty width. The
threshold is uniform in n and dirty width. Section 5.1 supplies an
explicit value and a quantitative rank--precision consequence.
Zero leaf phases give the same witness in the complex magnitude family.

Consequently every sequence of at-most-four-layer approximations to
these frames, with fixed clean budget and error tending to zero, has
**minimum crossing rank tending to infinity**, where the minimum is
over all four-layer decompositions. Both middle layers must then have
diverging numbers of active T gates. Section 5 first proves the qualitative claim using R73's imported
polynomial-limit premise and the normalization bound $`\beta\le r`$.
Section 5.1 gives an independent elementary quantitative proof, using
finite sums of affine phase states. No supplied states, new initialized
history, or additional circuit layers enter either analysis.

## 2. Peeling the outer layers

For a q-qubit circuit, absorb the outer Clifford factors into the
Hermitian witness Paulis P',Q'. A full-space transfer entry becomes

```math
x=2^{-q}\mathrm{Tr}
\bigl(A D_3C_2D_2B D_2^\dagger C_2^\dagger D_3^\dagger\bigr),
```

```math
A=C_3^\dagger(D_4^\dagger P'D_4)C_3,\qquad
B=C_1(D_1Q'D_1^\dagger)C_1^\dagger.
```

A single physical T layer conjugates a Pauli to a Hermitian Clifford,
so A and B are Hermitian Clifford unitaries. Put
$`B_0=C_2BC_2^\dagger`$ and $`E=C_2D_2C_2^\dagger`$. Then

```math
x=2^{-q}\mathrm{Tr}(A D_3E B_0E^\dagger D_3^\dagger).
```

Only the relative position of the two commuting axis spans is needed
to reduce the remaining middle Clifford transformation.

## 3. Simultaneous symplectic normal form

**Lemma.** Two isotropic subspaces L,M whose pairing has rank r admit
symplectic coordinates in which L is Z-type and contains
$`Z_1,\ldots,Z_r`$, while M is spanned by
$`X_1,\ldots,X_r`$ and a Z-type subspace on coordinates beyond r.
The residual Z subspaces of L and M may intersect.

Choose bases whose pairing matrix is $`\mathrm{diag}(I_r,0)`$
by invertible row and column operations. The first r pairs span a
nondegenerate symplectic subspace H: same-span pairings vanish, and the
cross-pairing matrix there is the identity. Every residual vector lies
in $`H^\perp`$. Their combined span is isotropic, since residual
cross-pairings also vanish. Choose a basis of that combined span,
extend it to an isotropic basis of $`H^\perp`$, and complete to a
symplectic basis. Sending the first pairs to $`Z_i,X_i`$ and the
combined residual span to later Z coordinates proves the lemma.
This construction includes a nonzero intersection $`L\cap M`$;
it never treats the two residual bases as jointly independent.

The symplectic change has a Clifford representative G
([Dehaene–De Moor, Sections II and IV](https://arxiv.org/pdf/quant-ph/0304125),
F44). Conjugation retains the actual Hermitian Pauli signs. In these
coordinates $`GD_3G^\dagger`$ is a diagonal parity-phase word F,
and

```math
GEG^\dagger=H_r F_0H_r,\qquad H_r=H^{\otimes r}\otimes I,
```

for another diagonal parity-phase word $`F_0`$. Each original axis
of E is a product of the first r X coordinates and later Z coordinates;
applying $`H_r`$ therefore makes it diagonal. Every original axis
and its rotation sign remains in the phase word. A basis change of the
span does not replace that word by rotations on only its new generators.

With $`A'=GAG^\dagger`$ and
$`B''=H_rGB_0G^\dagger H_r`$, both still Hermitian Clifford,

```math
x=2^{-q}\mathrm{Tr}
\bigl(A'F H_r F_0B''F_0^\dagger H_rF^\dagger\bigr).
```

These are equalities analyzing the given circuit. The diagonal words
are not assigned an uncharged implementation. Common scalar phases
of the words cancel in conjugation.

## 4. Affine supports and the normalization bound

Use normalized vectorization
$`|K\rangle\!\rangle=2^{-q/2}\sum_{x,y}K_{xy}|x,y\rangle`$.
The preceding identity is

```math
x=\langle\!\langle A'|
(F\otimes\overline F)(H_r\otimes H_r)
(F_0\otimes\overline{F_0})|B''\rangle\!\rangle.
```

The boundary states are stabilizer Choi states. By F44's affine
stabilizer form, their nonzero amplitudes have magnitudes
$`2^{-r_A/2}`$, $`2^{-r_B/2}`$ on affine supports of dimensions
$`r_A,r_B`$ in $`\mathbb F_2^{2q}`$. Relative phases are linear
modulo four and quadratic signs. Hermiticity restricts their common
phases to eighth roots, by the transposed-entry argument in
[R73, Section 3](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#3-three-layers-as-a-stabilizer-overlap-with-one-t-layer).

For doubled output and input indices z,w, the nonzero kernel entries
of $`H_r\otimes H_r`$ satisfy $`z_j=w_j`$ on all untouched
coordinates. Each row and column has $`2^{2r}`$ nonzero entries,
each of magnitude $`2^{-r}`$, with a quadratic sign. Intersect this
linear relation with the two affine boundary supports. If the
intersection is empty, x=0. Otherwise it is affine of dimension k.
Parameterizing it by $`u\in\mathbb F_2^k`$ gives

```math
\boxed{\quad x=2^\beta\mathbb E_u\zeta_8^{p(u)},\qquad
\beta=k-\frac{r_A+r_B}{2}-r.\quad}
```

Here p is a polynomial modulo eight of weighted cubic form

```math
p=c+\sum_i\ell_i u_i+
2\sum_{i\lt j}q_{ij}u_iu_j+
4\sum_{i\lt j\lt k}t_{ijk}u_iu_ju_k\pmod8.
```

Diagonal Pauli rotations contribute signed affine parities. Their
integer parity expansions truncate at degree three modulo eight,
with the displayed divisibilities. Boundary phases and the Hadamard
quadratic sign retain these divisibilities under affine restriction.
The possibly sixteenth-root scalar of each diagonal rotation word
cancels between the conjugate vectorization halves. The remaining
boundary scalar is an eighth root and enters the constant c. There
is no unrestricted continuous phase multiplying the expectation.

Projection of the compatible-pair space onto either boundary support
has fibers of dimension at most 2r. Therefore

```math
k\le r_A+2r,\qquad k\le r_B+2r,
```

and hence

```math
\boxed{\qquad\beta\le r-\frac{|r_A-r_B|}{2}\le r.\qquad}
```

The doubled transform contributes fiber dimension 2r and amplitude
$`2^{-r}`$. These two factors are both essential. In particular r=0
recovers a prefactor at most one even when the original middle
Clifford has large Hadamard rank.

## 5. Bounded-rank closure and initialized work

Let $`\mathcal S_{4,R}`$ be the set of full-space transfer entries
over all widths for circuits admitting a four-layer decomposition
with $`r\le R`$, and let
$`\mathcal K_{4,R}=\overline{\mathcal S_{4,R}}\subseteq[-1,1]`$.
We claim

```math
\mathcal K_{4,R}\subseteq\mathbb Z[1/\sqrt2].
```

For a sequence converging to a nonzero x, the bound
$`|x_j|\le2^{\beta_j}`$ eventually bounds $`\beta_j`$ below.
Section 4 bounds it above by R. The exponents are half-integers, so
pass to a subsequence with constant $`\beta`$, then to one on which
all eight residue probabilities of p converge. By
[R73, Section 4](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#4-polynomial-limits-and-the-closed-coefficient-set),
these residue events have Boolean degree at most 21, and F43's
polynomial-limit lemma makes their limiting probabilities dyadic.
Thus $`x\in\mathbb Z[1/\sqrt2,i]`$; reality puts it in
$`\mathbb Z[1/\sqrt2]`$. The zero limit belongs to that ring too.

For fixed a, the compact Minkowski sum of $`2^a`$ copies of
$`\mathcal K_{4,R}`$ is still contained in this ring. Its rational
elements are dyadic, so it excludes 1/3. Define

```math
\varepsilon_{a,R}=\tfrac12\mathrm{dist}\!\left(
\frac13,\underbrace{\mathcal K_{4,R}+\cdots+\mathcal K_{4,R}}
_{2^a\ \text{terms}}\right)>0.
```

The clean-projector expansion in R73 places the initialized transfer
entry in that sum, with no additional prefactor. Every summand uses
the same circuit decomposition and rank bound; only the Pauli pair
changes. An initialized-isometry error $`\eta`$ perturbs the
transfer witness by at most $`2\eta`$, including leakage. Padding
idle clean wires handles circuits using fewer than a clean qubits.
The complete frame in Section 1 has the target entry 1/3, proving the
theorem without replacing an initialized-isometry requirement by a
full-unitary requirement. Arbitrary dirty inputs and external
references retain the contract in R73.

If the minimum ranks of an approximating sequence failed to diverge,
an infinite subsequence would admit rank at most some fixed R,
contradicting this positive gap. Since
$`r\le\min\{|J_2|,|J_3|\}`$, both active middle layer sizes diverge
as well. The following proof strengthens this divergence to an explicit
rank--precision bound; it does not force more layers at positive clean
budget.

## 5.1. Quantitative rigidity on the initialized sector

The bounded-rank obstruction also has an elementary quantitative proof.
It uses the initialized observables themselves and does not require the
polynomial-limit premise of Section 5. Fix a clean budget a, let
$`\Pi=|0^a\rangle\langle0^a|\otimes I`$, and let P be a
Hermitian Pauli on the remaining logical and dirty wires. Normalize

```math
E=\Pi P,\qquad O=2^{a/2}E,\qquad
\|O\|_{2,\mathrm{tr}}=1,\qquad
x_E(V)=\langle O,VOV^\dagger\rangle_{2,\mathrm{tr}}.
```

The normalized Hilbert--Schmidt inner product uses the full physical
space. The scalar $`x_E`$ is real and at most one. In terms of ordinary
transfer entries it is a sum of $`2^{2a}`$ terms with prefactor
$`2^{-a}`$; it is different from Section 5's one-sided initialized
transfer witness.

**Initialized identity-isolation theorem.** Put

```math
K=2^{a+1},\qquad A=2(K^2-1),\qquad
S_a=\frac{A^K-1}{A-1},
```

```math
\mathcal T_a(r)=S_a\bigl[2Kr^2+(5K+8)r+7K+5a+5\bigr].
```

For every native four-layer decomposition of crossing rank r,

```math
\boxed{\qquad
x_E(V)=1\quad\hbox{or}\quad
1-x_E(V)\ge2^{-[\mathcal T_a(r)+1]}.
\qquad}
```

The constants are conservative; their useful feature is
$`\mathcal T_a(r)=O_a(r^2)`$. For one clean qubit they are explicitly

```math
\mathcal T_1(r)=27931(8r^2+28r+38).
```

### Conditional boundary states

Expand $`O=2^{-a/2}\sum_s PZ^s`$, where s ranges over the a clean
Z generators. Peeling the outer T layers as in Section 2 turns each
summand into a Hermitian Clifford. Thus each Choi boundary is a sum of
$`2^a`$ flat affine stabilizer states, with coefficient $`2^{-a/2}`$.
Section 3 uses the same middle normal form for every summand. On doubled
Choi coordinates it has diagonal weighted-cubic phases and a Hadamard
transform on $`m=2r`$ active coordinates.

The Choi discrepancy is consequently

```math
f=\sum_{i=1}^K s_i2^{-a/2}\xi_i,\qquad
s_i\in\{1,-1\},\qquad
\|f\|^2=2(1-x_E(V)).
```

A left state $`\xi_i`$ is flat on an affine support with a weighted
phase; a right state has the same form followed by the common diagonal
phases and active Hadamards. All are normalized. Their marginal on the
untouched coordinates is uniform on an affine space $`T_i`$, of
dimension $`h_i`$. Each nonempty active fiber is a normalized m-qubit
state. The middle operations preserve these marginals.

For nonempty intersections define

```math
d_{ij}=\frac{h_i+h_j}{2}-\dim(T_i\cap T_j),
```

and set $`d_{ij}=\infty`$ for disjoint supports. Cauchy--Schwarz inside
the active fibers gives $`|\langle\xi_i,\xi_j\rangle|\le2^{-d_{ij}}`$.

### A quantitative gap for one cluster

Consider projected supports connected by edges $`d_{ij}\le t`$.
Their affine hull H has codimension at most
$`\ell=2(K-1)t`$ above every member. Indeed, root a spanning tree at
that member. Adding a child support increases the hull dimension by at
most $`h_{\rm child}-\dim(T_{\rm child}\cap T_{\rm parent})\le2t`$.
A common intersection of all supports is unnecessary.

Parameterize H by $`u\in\mathbb F_2^h`$ and put
$`e_i=\mathrm{codim}_H T_i\le\ell`$. A cluster coordinate at
active output z is $`2^{-h/2}\lambda_z(u)`$, where

```math
\lambda_z(u)=\sum_i s_i2^{(e_i-a)/2}
\mathbf1_{T_i}(u)\,v_i(u,z).
```

Each $`v_i`$ is a normalized conditional m-qubit state. A left term is
flat on an affine active fiber of dimension $`b_i\le m`$. A right term
has the form

```math
v_i(u,z)=2^{-(b_i+m)/2}
\sum_{w\in\mathbb F_2^{b_i}}\zeta_8^{q_i(u,w)},
```

with output phases included in the weighted-cubic polynomial $`q_i`$.
Affine sections and their phases can be extended from $`T_i`$ to H
by an affine retraction.

Suppose the cluster is nonzero, and choose $`u_0,z_0`$ with a nonzero
coordinate. Fix all projected-support membership bits to their values
at $`u_0`$; this costs Boolean degree at most $`\sum_i e_i\le K\ell`$.
For a left term, fix its active-fiber membership at $`z_0`$ using at
most m affine constraints, and its phase residue using degree at most
seven. For a right term, expand $`q_i(u,w)`$ in w. Its constant,
linear-w, and quadratic-w coefficients have residue-event degrees at
most seven, three, and one respectively; cubic-w coefficients are
constant. Thus the per-term degree cost is at most

```math
d_m=7+3m+\binom m2\ge m+7.
```

These digit bounds follow from the weighted divisibilities: the three
binary digits of a weighted-cubic polynomial modulo eight have Boolean
degrees at most one, two, and four. Fixing one residue therefore costs
degree at most seven. The lower coefficient classes have digit degrees
one and two, or just one, respectively.
To check the digit bound directly, express the integer polynomial as a
sum of Boolean monomial tokens, including coefficient multiplicities.
Its binary digit j is the coefficient of $`z^{2^j}`$ in the product
of $`1+t z`$ over these tokens t, over the binary field. Quadratic
tokens occur in even multiples and first contribute at $`z^2`$; cubic
tokens occur in multiples of four and first contribute at $`z^4`$.
Every resulting monomial in that coefficient has degree at most
$`2^j`$. Constant shifts and the lower coefficient classes retain the
stated bounds.

The combined event contains $`u_0`$, preserves the chosen nonzero
coordinate, and has Boolean degree at most $`K\ell+Kd_m`$. Its measure
is at least $`2^{-K\ell-Kd_m}`$ by the minimum-weight bound for a
nonzero Boolean polynomial.

There is also a uniform arithmetic bound on that coordinate. With a
common denominator $`(\sqrt2)^{a+2m}`$,

```math
y=|\lambda_{z_0}(u_0)|^2
=\frac{u+v\sqrt2}{2^{a+2m}},\qquad u,v\in\mathbb Z.
```

The embedding $`\sigma(i)=i`$, $`\sigma(\sqrt2)=-\sqrt2`$ takes
each conditional flat phase state to another normalized state and
changes the active Hadamards only by a common sign. Therefore

```math
y'=\sigma(y)=|\sigma(\lambda_{z_0}(u_0))|^2
\le K^2 2^{\ell-a}.
```

The nonzero algebraic norm $`yy'`$ is at least
$`2^{-2(a+2m)}`$. Hence the chosen squared coordinate is at least
$`2^{-a-4m-\ell}/K^2`$. Averaging over the event proves that every
nonzero cluster has squared norm at least

```math
\delta(t)=\frac{2^{-(K+1)\ell-Kd_m-a-4m}}{K^2}
=\frac{2^{-At-Kd_m-a-4m}}{K^2}.
```

### Separating clusters without cancellation

Set

```math
C=Kd_m+a+4m+4\log_2K+1,\qquad
t_0=0,\qquad t_{j+1}=At_j+C\quad(0\le j\lt K).
```

Form the graph with edges $`d_{ij}\le t_j`$. Its component partition
can strictly merge at most $`K-1`$ times, so some two consecutive scales
among $`t_0,\ldots,t_K`$ have the same partition. Use those components
as clusters. Distinct clusters then have $`d_{ij}\gt t_{j+1}`$ for
every cross pair. Their total interference in $`\|f\|^2`$ is at most

```math
K^2 2^{-t_{j+1}}\le\delta(t_j)/2.
```

If f is nonzero, some cluster is nonzero. It follows that

```math
\|f\|^2\ge\delta(t_j)/2\ge2^{-t_{j+1}}\ge2^{-t_K}.
```

Since $`m=2r`$, the recurrence gives $`t_K=\mathcal T_a(r)`$,
proving the theorem.

### Finite-boundary form

The same proof applies to any sum of K normalized conditional boundary
states of the displayed form, with coefficients $`\pm2^{-b/2}`$,
where $`K\ge2`$ and the nonnegative integer b are independent.
Neither $`K=2^{a+1}`$ nor unit norm of an entire boundary sum is used
in the cluster argument. With

```math
A=2(K^2-1),\qquad d_m=7+3m+\binom m2,\qquad
C=Kd_m+b+4m+4\log_2K+1,
```

every nonzero such sum f satisfies

```math
\|f\|^2\ge2^{-T},\qquad
T=C\frac{A^K-1}{A-1}.
```

The Choi states of non-Hermitian Clifford products are also allowed
when their common phases are eighth roots. Fixed ring-valued Clifford
representatives have this property. Arbitrarily chosen continuous scalar
phases are not allowed. The [one-clean midpoint theorem](ONE_CLEAN_MIDPOINT_RIGIDITY.md#3-a-clifford-midpoint-with-unrestricted-crossing-rank)
uses this form at K=8, b=4, m=0 for an eight-term compressed discrepancy.

### Initialized return and a rank--precision bound

Let J append the a clean zeros and suppose
$`\|VJ-JU\|\le\eta`$, with U preserving the data Pauli P.
For the isometries $`F=VJ`$, $`G=JU`$ on a d-dimensional data-plus-dirty
space,

```math
\|FPF^\dagger-GPG^\dagger\|_F
\le2\|F-G\|_F\le2\eta\sqrt d.
```

Normalizing by the full dimension $`2^a d`$ gives

```math
1-x_E(V)=\tfrac12\|VOV^\dagger-O\|_{2,\mathrm{tr}}^2
\le2\eta^2.
```

Thus $`\eta\lt2^{-[\mathcal T_a(r)+2]/2}`$ forces exact preservation
of $`\Pi P`$. Taking P=I first gives exact terminal preservation of
the occupied clean sector. For the complete-frame witness, the
spectator and dirty X,Z observables and the logical Y then force the
unitary on that sector to be $`e^{i\phi}R_y(\alpha)\otimes I`$.
Its initialized transfer entries belong to $`\mathbb Z[1/\sqrt2]`$,
so R73's elementary unit-circle classification gives
$`\alpha\in(\pi/8)\mathbb Z`$. No conclusion about the unused clean
sectors or intermediate return is needed.

The witness therefore has the explicit distance lower bound

```math
\min\left\{
2^{-[\mathcal T_a(r)+2]/2},\quad
c_0\right\},\qquad
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}.
```

In particular, a four-layer approximation with error $`2^{-L}\lt c_0`$
must satisfy

```math
\boxed{\qquad \mathcal T_a(r)\ge2L-2.\qquad}
```

For every fixed clean budget this is a necessary
$`r=\Omega_a(\sqrt L)`$ condition, valid for the minimum over all native
four-layer representations. The same bound applies to each middle
layer's number of active T gates. Idle padding handles fewer than a
clean qubits. This does not bound unrestricted four-layer width or
prove a five-layer obstruction for a positive clean budget.

For three layers, peeling the two outer layers leaves a diagonal middle
word, so the same argument applies with m=0. It supplies an elementary
explicit fixed-clean obstruction with
$`\mathcal T_a(0)=S_a(7K+5a+5)`$. R73's sharper explicit no-clean
constant remains preferable in that case.

## 6. Growing normalization and the remaining scalar question

The bound $`\beta\le r`$ can be attained by actual native boundary
data. In the peeled two-middle-layer expression, take
$`A=B=C_2=H^{\otimes q}`$ and
$`D_2=D_3=T^{\otimes q}`$. The crossing rank is q. In Section 4's
form the boundary supports are both full, $`r_A=r_B=2q`$, and
the compatibility space has dimension $`k=4q`$, so $`\beta=r=q`$.
A one-qubit trace calculation, followed by tensor multiplication, gives

```math
x_q=\left(\frac{1+2\sqrt2}{4}\right)^q.
```

These boundaries are allowed by the outer-layer peeling: conjugating
$`TXT^\dagger=(X+Y)/\sqrt2`$ by an X-axis quarter-turn Clifford
gives $`(X+Z)/\sqrt2=H`$. The example certifies growth of the
prefactor within the native constraints. Its coefficient tends to
zero; it proves neither a nondyadic limit nor that this decomposition
minimizes crossing rank.

With the same unreduced boundaries but identity middle layers, a raw
Hadamard-kernel expansion would also have prefactor $`2^q`$, with
bias $`2^{-q}`$ and coefficient one. Its active crossing rank is
zero, and Section 3 removes that artifact. Raw Hadamard rank is thus
not a substitute for the active-axis invariant.

There is also a simple **abstract** obstruction to dropping the
bounded-prefactor hypothesis. For independent uniform bits, let

```math
f_k(a,x,y)=a\sum_{j=1}^k x_jy_j\pmod2,\qquad
\mathbb E(-1)^{f_k}=\frac{1+2^{-k}}2.
```

The identity follows by first conditioning on a, then averaging the
independent pairs. Sum m disjoint copies to obtain $`F_{k,m}`$ and
set $`p_{k,m}=4F_{k,m}\pmod8`$. This is a weighted cubic, with

```math
\mathbb E\zeta_8^{p_{k,m}}=2^{-m}(1+2^{-k})^m.
```

Choose $`m_k=\lfloor2^k\log(4/3)\rfloor`$, where log is natural,
and $`\beta_k=m_k-2`$. Then

```math
2^{\beta_k}\mathbb E\zeta_8^{p_{k,m_k}}
=\tfrac14(1+2^{-k})^{m_k}\longrightarrow\frac13.
```

Every finite value is dyadic and lies in $`[1/4,1/3]`$, since
$`\log(1+t)\le t`$; the exponent tends to infinity. The unscaled
phase probabilities tend to 1/2 and 1/2, consistently with F43.
Thus bounded phase degree, a real final value of magnitude at most
one, and an unrestricted dyadic prefactor cannot by themselves imply
ring closure. **No native four-layer realization of this abstract
family is claimed.** Its support and unitarity constraints would
require a separate argument.

An unrestricted scalar extension must control the scaled bias as
crossing rank grows, with native support and unitarity retained.
[R75](FOUR_LAYER_OPERATOR_RIGIDITY.md) instead uses operator norm and
exact Pauli commutants to prove a zero-clean five-layer obstruction. It
also reduces any near-one scalar pair to discrete phases in a common
active-axis algebra, but does not decide whether native rare-sector
phases realize a scalar counterexample. Positive-clean unrestricted
five-layer rigidity remains open. Optimal depth remains an open problem.
