# Two-layer rigidity with a preserved address

[Exact dirty queries](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Query/frame reduction](BOOLEAN_QUERY_FRAME_REDUCTION.md) ·
[Conditional frame transfer](QUERY_DEPTH_FRAME_TRANSFER.md)

An exact zero-clean dirty query with a nonaffine output cannot use two physical T
layers if the address is preserved at every Clifford-layer boundary and
the T gates act only on work. This restriction includes the native
read-only address interface of R91 and conditional R98. It does not
include every circuit whose address returns only at the end.
The theorem below gives a width-independent positive approximation
floor with zero initialized work. Its proof separates the address
dependence of Clifford interlayers from their fixed work operations.
Section 4 gives a separate initialized-work invariant: with one clean
qubit, two layers force degree at most two and quadratic polar rank at
most two, also at a fixed positive approximation radius.

## 1. Address-preserving Clifford interlayers

Let x be an a-bit computational address and let W contain w arbitrary
work qubits, including the query output and all dirty helpers. A Clifford
C that preserves x exactly has blocks

```math
C=\sum_x|x\rangle\langle x|\otimes C(x).
```

There are an x-independent work Clifford $`C(0)`$, binary-linear
labels $`u(x),v(x)\in\mathbb F_2^w`$, and a diagonal Clifford
phase $`\lambda(x)`$, normalized at zero, such that

```math
C(x)=\lambda(x)P(x)C(0),\qquad P(0)=I,\quad\lambda(0)=1.
                                                               \qquad(1)
```

Here P can be written as an ordered product of fixed Hermitian work
Paulis, one controlled by each address bit. Equivalently, including its
known phase in lambda, it is $`X^{u(x)}Z^{v(x)}`$. In this latter
convention the complete scalar is

```math
\lambda(x)=i^{\sum_i\alpha_i x_i}
 (-1)^{\sum_{i\lt j}\beta_{ij}x_ix_j},
\qquad \alpha_i\in\mathbb Z_4,\quad\beta_{ij}\in\mathbb F_2.
                                                               \qquad(2)
```

No address phase is discarded. In particular, address S and CZ gates
are included in (2). Each controlled-Pauli word and the diagonal
address phase are genuine Clifford operations.

To prove the factorization, first use exact preservation of the
stabilizers $`Z_{x_i}`$ to identify the induced work Clifford on the
zero-address code. After removing that fixed work Clifford, conjugation
of any work Pauli is itself times a Z-type address Pauli. The resulting
linear signs are implemented by an ordered list of address-controlled
work Paulis, using nondegeneracy of the work symplectic form. Removing
that list leaves an operator acting as a scalar on every work block.
It is still Clifford and therefore is a diagonal address Clifford.
Its standard generators S and CZ give (2). Reordering the controlled
work Paulis changes only a quadratic address sign, which stays in (2).

This proof permits arbitrary Clifford circuits inside an interlayer;
only its completed action must preserve x. The more restrictive
elementary interface in which x occurs only as a Clifford control
automatically satisfies the hypothesis.

## 2. Robust zero-clean obstruction

Set

```math
\kappa=\frac{\sqrt{2-\sqrt2}}2\approx0.382683.
```

**Theorem (zero-clean, two T layers).** Let V preserve the address
and have the native form

```math
V=C_2D_2C_1D_1C_0,                                  \qquad(3)
```

where each C is an address-preserving Clifford interlayer and each D
is one physical layer of T, T-dagger and identity gates on distinct
work wires. Empty layers are allowed. There are no initialized work
qubits. If, for a Boolean table f,

```math
\|V-Q_f\otimes I_{\rm dirty}\|_{\rm op}\lt \kappa,
\qquad Q_f|x,t\rangle=|x,t\oplus f(x)\rangle,
                                                               \qquad(4)
```

then f is affine over $`\mathbb F_2`$. The same conclusion holds
componentwise for multi-output XOR queries. If instead V approximates
the literal scalar $`(-1)^{f(x)}I_W`$ within this radius, then f
has algebraic degree at most two. The norm includes every occupied
work input and arbitrary external references, at any finite dirty
width. Optimizing a common global phase does not change either
conclusion.

### A flat Pauli-expansion gap

For a physical T layer D and a work Pauli P, let k be the number
of its active wires on which P anticommutes with Z. Then
$`DPD^\dagger`$ has exactly $`2^k`$ distinct Pauli terms of equal
coefficient magnitude $`2^{-k/2}`$. In normalized Hilbert--Schmidt
norm, its distance to any phase times a Pauli is therefore at least

```math
\sqrt{2-2^{1-k/2}}\ \ge\ \sqrt{2-\sqrt2}
\quad\text{if }k\ge1.                              \qquad(5)
```

Normalized Hilbert--Schmidt norm is at most operator norm. Thus
distance below the right side of (5) forces k=0, in which case D
commutes with P exactly. A scalar phase in P does not affect this
argument. Clifford conjugation preserves both norms and the Pauli
comparison class.

### Reduction to a Clifford address family

Apply (1) to all three interlayers. Write $`U_x=V(x)`$,
$`U_0=V(0)`$, and define

```math
M(x)=C_2(0)D_2P_1(x)D_2^\dagger C_2(0)^\dagger,
\qquad R(x)=C_0(0)^\dagger P_0(x)C_0(0).
```

The complete, phase-retaining identity is

```math
U_x=\lambda(x)P_2(x)M(x)U_0R(x),                    \qquad(6)
```

with $`\lambda=\lambda_2\lambda_1\lambda_0`$. All factors
other than M and the fixed $`U_0`$ are known projective Paulis
with linear address labels. Importantly, this identity does not assume
that the physical $`U_0`$ is exactly a Pauli or even Clifford.

Let the target fibre be $`F_x=X_t^{f(x)}\otimes I`$ and let
$`\eta=\|V-Q_f\otimes I\|\lt \kappa`$. Every fibre has error
at most eta, including $`\|U_0-F_0\|\le\eta`$. Solving (6)
for M and replacing both $`U_x`$ and the actual $`U_0^\dagger`$
by their target Paulis shows that M(x) is within $`2\eta`$ of
a phase times a Pauli. The estimate uses only unitary multiplication.
By (5), $`D_2`$ commutes with $`P_1(x)`$ for every x. Hence

```math
M(x)=C_2(0)P_1(x)C_2(0)^\dagger.                    \qquad(7)
```

Replace only $`U_0`$ in (6) by $`F_0`$ and call the resulting
family $`\widetilde U_x`$. Equations (1), (2) and (7) show that
$`\sum_x|x\rangle\langle x|\otimes\widetilde U_x`$ is a
Clifford on the full address/work space. Every fibre is a phase times
a Pauli, its Pauli label is affine in x, and

```math
\|\widetilde U_x-F_x\|\le2\eta\lt \sqrt2.             \qquad(8)
```

Distinct projective Pauli labels have operator-norm distance at least
$`\sqrt2`$, even after optimizing their scalar phases. Thus (8)
forces the fibre label to equal that of $`F_x`$ at every x. Its
target X component is f(x), proving that f is affine. The same argument
applies to every output coordinate of a word query.

For the scalar target, put $`F_x=(-1)^{f(x)}I`$ throughout.
Equation (8) forces every work Pauli label of the constructed Clifford
to vanish. It is therefore a diagonal address Clifford. Since
$`\widetilde U_0=F_0`$, its phases relative to the zero-address
value belong to $`\{1,i,-1,-i\}`$. A phase $`\pm i`$ has
distance $`\sqrt2`$ from either sign, so (8) forces
$`\widetilde U_x=(-1)^{f(x)}I`$ exactly. Formula (2) then gives
degree at most two. A common scalar phase is removed in both comparisons
and in the normalization at x=0, without changing (5).

## 3. Meaning and limits

For a nonaffine Boolean table, any zero-clean circuit in this native
read-only class requires at least three T layers at error below kappa.
At the repository's target precision $`L\ge6`$, this threshold is
more than sufficient. The conclusion is a constant-layer obstruction,
not a dimension-growing lower bound and not optimality of R91.

The [earlier unrestricted full-input result](../../docs/PARALLEL_DIRTY_LOOKUP.md#11-a-full-input-lower-bound)
already excludes one T layer for a nonaffine permutation. The new
two-layer statement needs the preserved-address interfaces of (3).
There is an explicit unrestricted two-layer dirty Toffoli in the
[same proof home](../../docs/PARALLEL_DIRTY_LOOKUP.md#12-the-two-bit-indicator-has-optimal-t-depth-two):
its two parity bases place physical T gates on transformed address
wires. It therefore lies outside (3) and prevents promoting this
read-only theorem to an unrestricted three-layer bound.

Initializing work changes the argument: the zero-address target no
longer determines a full work unitary, and occupied-sector compression
cannot be substituted in (5). Section 4 instead uses a projective
Pauli-representation invariant. With one clean qubit it excludes cubic
queries and quadratic queries of polar rank greater than two at two
layers, including approximate return. It does not exclude every
nonaffine one-clean query. Likewise, preserving the address only at the
end does not provide (1) at each T-layer boundary.
The query-to-frame lower-bound transfer accepts unrestricted native
queries; this scoped theorem does not, by itself, supply its missing
dimension-growing premise. Conditional R98 requires the read-only
interface but a genuinely sublogarithmic construction, which a
three-layer exclusion does not rule out.

## 4. Initialized work: a robust degree and polar-rank bound

Let $`J_c=|0^c\rangle\otimes I`$ initialize c clean work qubits;
all remaining work, including the output t, is arbitrary. The full
address circuit still satisfies (3), with every Clifford interlayer
preserving the address and every physical T gate acting on work. Set

```math
\rho=\frac{\sqrt2}{6}\approx0.235702.
```

**Theorem (two layers with initialized work).** If

```math
\|VJ_c-J_c(Q_f\otimes I_{\rm dirty})\|_{\rm op}
\le\eta\lt \rho,                                          \qquad(9)
```

then the following necessary conditions hold:

| Initialized work | Necessary condition on f |
|---|---|
| c=0 | f is affine; section 2 gives the stronger radius kappa. |
| c=1 | f has degree at most two, and its quadratic polar form has rank at most two. |
| $`c\ge1`$ | $`\deg_{\mathbb F_2}f\le2c`$. |

The same conclusions hold coordinatewise for an XOR word query. There
is no bound on the number of dirty qubits. The hypothesis permits
nonzero clean leakage; (9) is an initialized-isometry bound, not a
promise of exact work return. Tensoring external references preserves
every norm used below. The argument establishes necessary conditions,
not constructions attaining them.

### Three Pauli-control bands

The factorization (6) has the abstract form

```math
U_x\doteq P(x)A(x)WR(x),\qquad
A(x)=V_*S(x)V_*^\dagger,                              \qquad(10)
```

where P, S and R are Pauli families with binary-linear labels,
$`P(0)=S(0)=R(0)=I`$, and $`V_*,W`$ are fixed work unitaries.
The symbol $`\doteq`$ means equality up to a scalar of modulus one;
it is used only for the invariant, never to weaken the target contract
(9). In (6), take $`P=P_2`$, $`S=P_1`$,
$`V_*=C_2(0)D_2`$, and $`W=U_0`$. The argument below actually
allows arbitrary fixed $`V_*,W`$ in (10). Thus it needs no closure
claim about the Clifford hierarchy.

For a physical Pauli P, its clean syndrome is the c-bit string recording
which initialized Z stabilizers it anticommutes with. Let
$`s_P,s_R:\mathbb F_2^a\to\mathbb F_2^c`$ be the linear
syndrome maps of the endpoint families, and define

```math
H=\ker s_P\cap\ker s_R,\qquad \mathrm{codim}H\le2c.
                                                               \qquad(11)
```

For $`h\in H`$, both endpoint Paulis preserve the initialized
subspace. Write their restrictions there as $`p(h),r(h)`$.
They are logical Paulis with linear labels.

Fix any address $`x_0`$ and put

```math
B_{x_0}(h)=P(x_0)A(h)P(x_0)^\dagger.
```

Projective linearity of each Pauli family, including the conjugated
family A, gives the full-space identity

```math
U_{x_0+h}\doteq P(h)B_{x_0}(h)U_{x_0}R(h),
\qquad
B_{x_0}(h)\doteq
P(h)^\dagger U_{x_0+h}R(h)^\dagger U_{x_0}^\dagger.
                                                               \qquad(12)
```

For each $`x_0`$, $`h\mapsto B_{x_0}(h)`$ is an exact
projective representation of the additive group H.

### The actual inverse and leakage-safe comparison

Write $`F_x=X_t^{f(x)}\otimes I`$ on the nonclean work. Each
address fibre obeys $`\|U_xJ_c-J_cF_x\|\le\eta`$. Its actual
inverse obeys the same bound, because

```math
U_x^\dagger J_c-J_cF_x^\dagger
=U_x^\dagger(J_cF_x-U_xJ_c)F_x^\dagger.               \qquad(13)
```

In (12), first replace the rightmost inverse on $`J_c`$, then the
forward call. The intervening $`R(h)^\dagger`$ preserves the
initialized subspace, so both uses of (9) have initialized inputs.
Leakage errors are propagated only by unitaries. For an appropriate
scalar $`\mu_{x_0,h}`$ this proves

```math
\|B_{x_0}(h)J_c-\mu_{x_0,h}J_cL_{x_0}(h)\|\le2\eta,
\qquad
L_{x_0}(h)=p(h)^\dagger F_{x_0+h}r(h)^\dagger F_{x_0}^\dagger.
                                                               \qquad(14)
```

Every L is a logical Pauli. Its binary Pauli label is

```math
\mathrm{lab}L_{x_0}(h)
=\mathrm{lab}p(h)+\mathrm{lab}r(h)
 + e_{X_t}\bigl(f(x_0+h)+f(x_0)\bigr).                \qquad(15)
```

For distinct Pauli labels, the minimum operator distance over scalar
phases is $`\sqrt2`$: the quotient of their Hermitian representatives
has two opposite eigenvalues, either $`\pm1`$ or $`\pm i`$.
Apply (14) twice to $`B(h)B(k)J_c`$ and once to
$`B(h+k)J_c`$. The exact projective multiplication law then puts
the corresponding ideal Pauli product within $`6\eta\lt \sqrt2`$
of a scalar times $`L(h+k)`$. Their labels must coincide. Thus
$`h\mapsto\mathrm{lab}L_{x_0}(h)`$ is linear. By (15),

```math
h\longmapsto f(x_0+h)+f(x_0)
\quad\hbox{is linear on H}.                           \qquad(16)
```

This comparison applies (14) only to logical Pauli images of the
initialized subspace. It does not require B to return the clean qubits
exactly, nor apply an initialized estimate to a leaked vector.

### Why one clean gives degree two rather than three

Affineness on codimension-two cosets alone would permit degree three.
The additional constraint is that the slopes in (16) depend on the
coset representative only through its **left** syndrome $`s_P(x_0)`$.
Indeed, if $`s_P(x_1)=s_P(x_0)`$, then

```math
D=P(x_1)P(x_0)^\dagger
```

is a code-preserving Pauli, and
$`B_{x_1}(h)=D B_{x_0}(h)D^\dagger`$ up to a scalar.
Its logical restriction conjugates every logical Pauli without changing
that Pauli's label. Comparing the two uses of (14) gives a distance at
most $`4\eta\lt \sqrt2`$ between the candidate labels. Hence those
labels, and by (15) their slopes, coincide.

Choose linear coordinates $`(z,h)`$ on address space, with h spanning
H and z the quotient of dimension $`k\le2c`$. The slope conclusion
has the precise form

```math
f(z,h)=a(z)+\sum_j h_jm_j(s_P(z)).                    \qquad(17)
```

Here a is an arbitrary Boolean function on k bits, and each m is an
arbitrary Boolean function of at most c independent syndrome bits.
Consequently
$`\deg f\le\max(k,c+1)\le2c`$ for $`c\ge1`$. With no
clean qubits, H is the entire address space and (16) gives affineness.

For $`c=1`$, if $`s_P=0`$ then k is at most one and the slopes
are constant, so f is affine. Otherwise take $`z_1=s_P(x)`$ and,
if necessary, a second quotient coordinate $`z_2`$. Equation (17)
becomes

```math
f(z,h)=a(z_1,z_2)+\ell(h)+z_1m(h),
```

with linear $`\ell,m`$. Its quadratic part is therefore the
product of two linear forms
$`z_1(\alpha z_2+m(h))`$, up to absorbed linear terms. Its polar
bilinear form has rank at most two. For example, both

```math
f=x_1x_2x_3,\qquad
f=x_1x_2+x_3x_4
```

are excluded at two layers with one clean qubit and error below rho:
the first has degree three, and the second has polar rank four.

The bound $`\deg f\le2c`$ can also be read as a clean-width
requirement within this two-layer interface. It supplies no
dimension-growing T-depth bound at fixed c: the argument controls
exactly three Pauli-control bands and has not been extended to an
arbitrary number of T layers. Moreover, it does not apply when address
wires participate in physical T layers or are restored only at the
end. The unrestricted two-layer example in section 3 remains a decisive
test against such a scope enlargement.

### Bounded verification

The analytic proof is supplemented by
[`test_read_only_clean_rigidity.py`](../../tests/test_read_only_clean_rigidity.py).
It checks the projective Pauli gap, the stated polar-rank witnesses, and
an exact native three-band example with one initialized qubit and every
dirty input. The example uses two arbitrary fixed conjugating bands;
it is **not** asserted to have two physical T layers. Its literal word
has 14 T gates, 26 elementary Clifford gates, T-depth at most 10, two
address bits, one arbitrary output, one dirty helper and one clean
helper. It implements $`f(a,b)=ab`$ and demonstrates why the
three-band argument with general fixed unitaries cannot replace its
one-clean quadratic conclusion by affineness. A second native fixture
adds a small nonzero clean-leakage perturbation and checks (12)--(14)
using actual circuit inverses. These finite checks supply no growing
lower bound or general synthesis theorem.

[`test_read_only_two_layer_rigidity.py`](../../tests/test_read_only_two_layer_rigidity.py)
independently checks the literal Clifford normal form, including address
S and CZ phases, the identity (6) with the actual non-Pauli zero-address
unitary, the physical T-layer flat Pauli expansion, and the unrestricted
two-layer dirty Toffoli on all sixteen occupied inputs. The latter test
also records its physical T gates on address wires, the precise reason
it lies outside this chapter's hypotheses.
