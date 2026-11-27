# Four abstract Pauli representations can return any Boolean phase

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Five-layer native reduction](READ_ONLY_FIVE_LAYER_REDUCTION.md) ·
[Earlier partial-normalizer counterexample](ONE_BAND_AXIS_MATCHING.md#5-an-abstract-normalizer-can-have-arbitrarily-high-tableau-degree)

The number of factors in a scalar-returning product of fixed-unitary
conjugates of linear-label Pauli representations does **not** bound the
Boolean degree of the returned phase. Four factors suffice for every
Boolean function vanishing at zero, in arbitrary dimension. The factors
below are honest representations with fixed literal Pauli gauges. There
is no address-dependent scalar hidden in their definition.

The conjugating sources have an explicit, growing native cost. Therefore
this is a counterexample to an abstract representation-product invariant,
not a constant-depth query construction or a counterexample to the
retained native five-layer conjecture. It identifies information that a
successful degree lower bound must retain: the source circuits and their
common placement in the native layer sequence.

For [two projective address visits](#two-visit-marked-width), exact
signed conjugation can require exponential work, with a matching
[paid representation](#paid-two-visit-exact-table). Completely returned
words in this shape also have a [width-independent fixed-error screen](#two-visit-returned-word-packing)
from the overlap of their two target eigenprobe sectors.

## 1. A family of Pauli generators carrying an arbitrary table

Let a be positive, let $`Q=2^a`$, and identify the bank labels with
the integers from zero to Q minus one, with bitwise addition. Let
$`f:\mathbb F_2^a\to\mathbb F_2`$ satisfy f(0)=0. Use Q arbitrary
bank qubits, one at each label u, and define literal Hermitian Paulis

```math
\gamma_u=X_u\prod_{v\lt u}Z_v^{f(u+v)}.
```

Each has square I, and direct Pauli multiplication gives

```math
\gamma_u\gamma_v=(-1)^{f(u+v)}\gamma_v\gamma_u,
\qquad \gamma_0=X_0.                                      \qquad(1)
```

The assumption f(0)=0 is exactly the alternating condition on the
commutation table. No assumption on the degree, sparsity or symmetry
beyond the displayed dependence on u+v is needed.

For a bank basis string s, write its coordinates as s_u. In increasing
order of the generator indices,

```math
\gamma_0^{s_0}\gamma_1^{s_1}\cdots
\gamma_{Q-1}^{s_{Q-1}}|0^Q\rangle=|s\rangle.              \qquad(2)
```

The rightmost generators act first; each Z factor therefore sees a zero
bit. Equation (2) is an identity used to define operators, not an
initialized work assumption.

## 2. An exact translation action with its literal signs

Let $`\pi_x(u)=u+x`$. Define the quadratic bank phase

```math
q_x(s)=\sum_{\substack{u\lt v\\\pi_x(u)\gt \pi_x(v)}}
 f(u+v)s_us_v\pmod2
```

and the monomial Clifford

```math
C_x|s\rangle=(-1)^{q_x(s)}|\pi_xs\rangle,
\qquad (\pi_xs)_{u+x}=s_u.                               \qquad(3)
```

Reordering the permuted generator word in (2), using (1), gives exactly
the sign in (3). Translation preserves every commutator in (1).
Consequently substitution of translated generators is an automorphism,
including all literal signs. It follows that

```math
C_x\gamma_u C_x^\dagger=\gamma_{u+x},\qquad
C_xC_y=C_{x+y},\qquad C_0=I,\qquad C_x^\dagger=C_x.       \qquad(4)
```

These identities can also be checked directly on (3): the signs count
the generator pairs whose relative order is reversed, and two successive
reorderings count the same pair twice when it is reversed twice. The
vacuum is fixed, so there is no residual projective scalar.

Each C_x is a wire permutation after the CZ gates specified by q_x.
Its dependence on x can have high degree; being Clifford-valued at every
x does not make this a linear-label Pauli family.

## 3. Four honest Pauli images return the prescribed scalar

Introduce a further register z of a arbitrary work qubits and the fixed
unitary source

```math
U=\sum_z|z\rangle\!\langle z|\otimes C_z,
\qquad U^\dagger=U,\quad U^2=I.
```

Using the exact representation law in (4), on every z and every bank
state,

```math
A_x=U(X_z^x\otimes I)U^\dagger=X_z^x\otimes C_x,
```

```math
U_\gamma=(I\otimes\gamma_0)U,\qquad
B_x=(I\otimes\gamma_0)A_x(I\otimes\gamma_0)
   =U_\gamma(X_z^x\otimes I)U_\gamma^\dagger.             \qquad(5)
```

Thus both A and B are fixed-unitary conjugates of the same canonical
linear Pauli representation. They obey the **exact** laws
$`A_xA_y=A_{x+y}`$ and $`B_xB_y=B_{x+y}`$; every member is an
involution. From (4),

```math
A_xB_x=I_z\otimes\gamma_x\gamma_0,
```

and hence

```math
A_xB_xA_xB_x=(-1)^{f(x)}I_{z,\mathrm{bank}}.             \qquad(6)
```

All work returns, with the literal scalar shown, for arbitrary work
inputs and arbitrary external references. No measurement, postselection,
leakage or approximation is involved. In particular, choosing the
degree-a monomial f gives unbounded Boolean degree at the fixed factor
count four. Allowing projective representations, approximate return, or
more than four factors cannot restore a degree bound based on that factor
count alone, since the exact honest example is already included.

This does not contradict the
[returned diagonal-echo theorem](RETURNED_DIAGONAL_ECHO_DEGREE.md).
The source U in (5) is a controlled, generally non-diagonal Clifford
action. It is not the diagonal source in that theorem.

## 4. A fully charged native source and returned circuit

The source U is explicitly realizable without any initialized qubit or
extra helper. Since the C_x form an exact commuting representation,
implement it by the controlled operators C_ej, with work bit z_j as
control, for the a coordinate vectors e_j. For each coordinate:

1. Apply the CZ gates in q_ej, controlled by z_j. These are CCZ gates.
2. Apply the bank translation by e_j, controlled by z_j. These are
   Q/2 Fredkin gates, on disjoint bank pairs.

The first step precedes the second, as required by (3). Across all
coordinates put

```math
s=\frac{aQ}{2},\qquad m=\frac{Q\,\mathrm{wt}(f)}2,
```

where wt(f) counts the nonzero entries of the table. There are s Fredkin
gates and m CCZ gates. For m, each unordered pair u<v is reversed by
exactly the basis translation corresponding to its highest differing
bit. Thus the inversion sets partition the unordered pairs. Each
nonzero difference with f equal to one occurs in Q/2 pairs.

Use the retained literal CCZ expansion with seven T/T-dagger gates,
ten CNOTs and five T layers. A Fredkin uses two CNOTs and a Toffoli;
the latter is H, CCZ, H, for a total of fourteen elementary Cliffords.
A conservative serial source ledger is therefore

```math
T(U)=7(s+m),\qquad G(U)=14s+10m,\qquad
D_T(U)\le5(s+m).                                        \qquad(7)
```

All of these gates act on the a+Q arbitrary work qubits. The original
address is absent from U. There are zero clean qubits and zero additional
helpers. No simultaneous parallelism of the controlled gates is assumed.

For a coherent original address x, each canonical word X_z^x in (5)
uses a CNOTs. Expanding all four factors in (6) uses eight source or
actual inverse calls, four such CNOT words, and four gamma_0 gates.
Using X=HSSH in the fixed H,S,CNOT Clifford alphabet gives

```math
\begin{aligned}
T&=56(s+m), & G&=112s+80m+4a+16,\\
D_T&\le40(s+m), & b&=a+Q,\qquad c_{\rm clean}=0.
\end{aligned}                                          \qquad(8)
```

Each inverse in this ledger is the reversed native inverse circuit;
the abstract involution identity does not make a call free. Every
physical T is on work, and completed Clifford interlayers preserve the
original address pointwise. Equation (6) proves zero error, zero
leakage and complete workspace return for this expanded circuit.

Algebraically, the X_z factors cancel and (6) becomes a word in four
C_x factors and four gamma_0 factors. This is not permission to charge
C_x as a free address-dependent Clifford: its coherent control must be
constructed. In particular, directly expanding controlled Clifford
gates with original address controls can place T on those controls and
leave the retained read-only interface. The source ledger (8) stays
within that interface.

The construction is not competitive with the retained query baseline.
Apply the exact one-output [R91 query](PARALLEL_BILINEAR_WORD_QUERY.md)
twice around two Z gates on its arbitrary output t:

```math
Z_tQ_fZ_tQ_f=(-1)^{f(x)}I.
```

The two Z gates cost four S gates. With a dirty allocation of order Q,
the R91 hypotheses hold for all sufficiently large a, and its output
and helpers fit in the a+Q work budget. This baseline has
$`T=O(\sqrt Q)`$, G=O(Q), and
$`D_T=O(\log_2(a+2))`$, with exact return, zero clean qubits and
both queries charged. Thus (8) supplies evidence about the failed
invariant, not a resource saving. The displayed source synthesis is not
claimed optimal. Section 5 supplies a mechanism-specific source lower
bound; no lower bound on the final scalar is asserted.

## 5. A source-cost barrier survives Pauli-carrier compression

The source cost in (7) is not solely a defect of its gate-by-gate
implementation. A rank obstruction applies to every **Clifford
translation source with Pauli carriers** of the following form. Allow w
arbitrary carrier qubits, Hermitian Pauli operators gamma_u satisfying
(1), and a literal Clifford representation C_x satisfying (4). They need
not use the bank embedding in Section 1. Retain the coherent source U in
(5), on a arbitrary selector qubits and these w carriers. Extra arbitrary
helpers are allowed if the full source returns them as identity. There
are zero initialized helpers. The Clifford hypothesis applies to C_x on
the entire carrier space, not only to its action on the listed gamma_u.

Let F be the binary Q by Q matrix, and P_x the permutation matrix,
defined by

```math
F_{uv}=f(u+v),\qquad P_x e_u=e_{u+x}.
```

For every such source and every nonzero x,

```math
\begin{aligned}
w&\ge\frac{\mathrm{rank}_{\mathbb F_2}F}{2},\\
T(U)&\ge\mathrm{rank}_{\mathbb F_2}\bigl(F(P_x+I)\bigr).
\end{aligned}
                                                               \qquad(9)
```

These are lower bounds for this specified source mechanism. They are not
lower bounds on synthesizing the final scalar, or on general returned
queries. In particular, they do not charge independent source calls
after a compiler has replaced the entire word by a different mechanism.

### 5.1 Binary commutation rank and displacement

Write the binary Pauli labels as columns of a matrix Gamma with 2w
rows, and let J be the nondegenerate binary symplectic form. Equation
(1), including its signs, gives

```math
F=\Gamma^{\mathsf T}J\Gamma.
```

This proves the width bound. If S_x is the symplectic action of C_x,
the translation identity gives

```math
\begin{aligned}
S_x\Gamma&=\Gamma P_x,\\
F(P_x+I)&=\Gamma^{\mathsf T}J(S_x-I)\Gamma.
\end{aligned}
```

Consequently

```math
\mathrm{rank}(S_x-I)
 \ge\mathrm{rank}\bigl(F(P_x+I)\bigr).              \qquad(10)
```

Literal Pauli signs and Clifford phases do not alter these binary
equalities. They remain part of the source and returned-word identities;
no phase is discarded from the implemented operator.

### 5.2 A physical T gate can at most double Pauli support

A Clifford V with symplectic action S has exactly
$`2^r`$ nonzero Pauli coefficients, all of magnitude
$`2^{-r/2}`$, where $`r=\mathrm{rank}(S-I)`$.
Here the Pauli basis is orthonormal in normalized Hilbert--Schmidt norm.
For completeness, Paulis fixed by S impose affine sign constraints on
the support of V. Its support is therefore contained in one coset of
$`(\ker(S-I))^\perp=\mathrm{im}(S-I)`$.
The relation between V P and its Clifford image times V translates
coefficients by every vector in this image and preserves their
magnitudes. A nonzero coefficient thus fills that entire coset.
Normalization gives the stated equal magnitudes.

In contrast, conjugating a Pauli by a circuit with t physical T or
T-dagger gates produces at most $`2^t`$ Pauli terms. Clifford gates
permute terms; each T conjugation fixes a commuting term or splits an
anticommuting term into two. This is a count statement, independent of
how many wires or T layers the circuit uses.

The specified full source obeys

```math
U(X_z^x\otimes I)U^\dagger=X_z^x\otimes C_x.
```

The right side is Clifford, with displacement rank
$`\mathrm{rank}(S_x-I)`$. Its exact Pauli support and the
doubling bound prove $`T(U)\ge\mathrm{rank}(S_x-I)`$,
and hence (9). Appending returned arbitrary helpers tensors both sides
with identity and changes neither rank nor support. An initialized-only
source approximation is not covered by this argument.

The bound is robust when a native source is supplied in full operator
norm. If a unitary circuit V, with t T gates and the same arbitrary
workspace contract, satisfies

```math
\|V-U\otimes I_{\rm helpers}\|_{\rm op}\le\delta\lt\frac12,
```

then its conjugated Pauli is within $`2\delta`$ of the displayed
Clifford. Projection onto at most $`2^t`$ Pauli terms leaves squared
normalized Hilbert--Schmidt error at least
$`1-2^{t-r}`$ whenever t is less than r. Thus

```math
\begin{aligned}
t&\ge r+\log_2(1-4\delta^2),\\
r&=\mathrm{rank}(S_x-I).
\end{aligned}                                          \qquad(11)
```

In particular, the integer lower bound $`t\ge r`$ is unchanged for
$`\delta\lt1/(2\sqrt2)`$. This hypothesis concerns each approximated
source itself; a small error of the completed scalar word does not imply
it. When eight such source or actual-inverse calls approximate the eight
specified exact calls, their errors add to at most $`8\delta`$ by
telescoping. No reverse implication is used.

With h simultaneous arbitrary returned helpers, the source has a+w+h
physical wires. A T layer contains at most that many T gates, giving the
elementary count-to-width consequence

```math
D_T(U)\ge
\left\lceil
\frac{\mathrm{rank}(F(P_x+I))}{a+w+h}
\right\rceil.
```

The same consequence holds for the small full-source approximation
regime above. It is only a count-to-width bound. In particular, (12)
already requires width of order Q, so this argument supplies no
dimension-growing depth lower bound at large width. No independent
elementary Clifford-count lower bound is claimed.

### 5.3 Full-degree tables cannot enter the R91 count envelope this way

Over the binary field, the translation matrix satisfies

```math
F^2=(\mathrm{wt}(f)\bmod2)I.
```

For an off-diagonal entry the summands pair under translation by its
nonzero difference; the diagonal entry is the parity of the table.
If wt(f) is odd, F is invertible. Every nonzero translation P_x is Q/2
disjoint transpositions, so (9) becomes

```math
w\ge Q/2,\qquad T(U)\ge Q/2.                            \qquad(12)
```

Odd weight is exactly the condition that the coefficient of the
degree-a monomial in the algebraic normal form of f is one. Thus (12)
applies to every full-degree table, not only to a sparse example.

More generally, if f has Boolean degree d at least one, restrict F to
the subgroup supported on a degree-d monomial. The restricted table has
odd weight, so that principal submatrix has full rank $`2^d`$.
For a coordinate e_i in that monomial, the derivative
$`\Delta_{e_i}f`$ has degree d minus one. The same argument applies
to its translation matrix $`F(P_{e_i}+I)`$. Therefore

```math
w\ge2^{d-1},\qquad T(U)\ge2^{d-1}.                      \qquad(13)
```

For the original uncompressed Q-qubit bank there is an additional
implementation-specific bound $`T(U)\ge Q`$, even when f is zero.
For nonzero x the Clifford C_x is a quadratic diagonal phase followed
by a fixed-point-free bank permutation. Its symplectic displacement has
both diagonal blocks P_x-I, so its rank is at least Q. Since C_x is an
involution, its displacement squares to zero in dimension 2Q, giving
rank at most Q. Its rank is consequently exactly Q. This illustrates
why an expensive source can persist even when the final scalar is
trivial.

Multiplying U by the fixed Pauli gamma_0, or taking its actual inverse,
does not change its minimum T-count: the inverse transformation is Clifford
multiplication or native circuit reversal, respectively. The same statement
holds for full-operator source approximation.

If the eight source calls in (8) are implemented separately, (12)
forces at least $`4Q`$ T gates for any full-degree table, including
compressed carriers and the small source-approximation regime above.
Increasing dirty width cannot reduce that count. R91 at dirty width
of order Q implements the final scalar with $`O(\sqrt Q)`$ T
gates, G=O(Q), and logarithmic address depth. The source-call mechanism
therefore cannot provide the general simultaneous count envelope needed
by R98, regardless of a faster schedule for those same sources.

This does not forbid cancellation or a different native synthesis of
the **whole** scalar word. Such a synthesis is exactly the necessary
constructive change: it must avoid separately paying for the full
Clifford translation source. It must also retain literal return, price
every control and inverse, and meet the T-count, Clifford-count, depth
and simultaneous-width bounds together. No lower bound here excludes
non-Clifford carriers, different source structures, or general queries.

## 6. Whole-word fusion removes carrier motion but leaves the table bit

The source lower bound in Section 5 cannot be transferred to an arbitrary
resynthesis of the completed word. Indeed, the selector Pauli factors
cancel exactly, and (6) reduces to

```math
\bigl(C_x\gamma_0C_x^\dagger\gamma_0\bigr)^2
 =(\gamma_x\gamma_0)^2=(-1)^{f(x)}I.
```

Since $`\gamma_0=X_0`$, only the Z_0 coefficient of gamma_x
affects this commutator. That coefficient is exactly f(x), including
x equal to zero because f(0)=0. The selected X factor and every other
Z factor commute with gamma_0 and disappear from the completed phase.
Thus fusion can legitimately eliminate the unused carrier motion; it
leaves the task of loading f(x) into the surviving probe coefficient.
Treating that coefficient as already available would omit the query.

For example, splitting the address as x=(u,v) expresses the table as

```math
f(x)=e_u^{\mathsf T}D e_v.
```

The retained native way to supply those one-hot coordinates is precisely
the [R91 dirty-indicator construction](PARALLEL_BILINEAR_WORD_QUERY.md).
Its exact returned scalar implementation, with both queries charged,
has $`O(\sqrt Q)`$ T gates, G=O(Q), and
$`O(\log_2(a+2))`$ T-depth at dirty width of order Q. The source
and whole-word count bounds are therefore consistent. The proposed
algebraic fusion alone does not replace the logarithmic address step.

A constructive improvement needs an explicit native probe-transport or
table-bit loading identity that avoids that step while retaining the
simultaneous count, elementary Clifford and dirty-width envelope. No such
replacement was obtained in this attempt. This diagnoses the candidate;
it is not a lower bound against a different fused synthesis, and it adds
no circuit, uncharged source, approximation or workspace assumption.

## 7. Consequence for the general depth question

A scalar obtained by commuting a returned XOR query with a fixed probe
can indeed be written as a bounded product of Pauli representations.
Discarding the native conjugators at that point loses necessary
information: (6) realizes every truth table within precisely that
abstract class. A successful invariant must use more than factor count,
even with literal gauges and exact full return.

The [native-defect theorem](READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
now uses the shared prefix structure to bound five-layer read-only XOR
degree. A depth-dependent recurrence remains the next invariant question;
the independently supplied, growing-cost sources in (5) cannot replace
that native structure. This counterexample does not give a five-layer realization,
does not improve R96, and does not change the scoped R97 obstruction.

The [bounded checks](../../tests/test_scalar_representation_product.py)
verify the Pauli commutators, literal translation action and group law,
canonical-source identity, returned scalar on every work column, actual
inverse source word, and exact source gate counts in small instances.
The arbitrary-width conclusion is the analytic construction above.

<a id="two-visit-marked-width"></a>

## 8. Two projective visits need exponential work for exact false silence

The four-factor scalar identity above is distinct from a two-factor observable conjugation. Let $`x\in H=\mathbb F_2^r`$, let B and C be arbitrary projective unitary representations on q work qubits, and put $`S_x=B_xC_x`$, normalized at zero. Their fixed work conjugators may contain arbitrary noncommuting dynamics. Let P and Q be Hermitian observables with $`P^2=\Pi_P`$, $`Q^2=\Pi_Q`$, and both projectors of rank $`2^{q-k}`$. Allow different supports. Suppose

```math
\|S_xPS_x^\dagger-Q\|\le\delta\quad(x\ne0),\qquad
\|P+Q\|\le\eta\lt2.
```

At $`\delta=0`$, writing $`m=2^r`$ and $`\beta=2^{-k}(1-\eta/2)`$, one necessarily has

```math
q\ge(m-1)\beta/4.
```

The same width order holds for $`\delta\le(1-\eta/2)/(4\pi m)`$. This is a representation bound, independent of the actual T-depth and arbitrary dirty width. A read-only row probe with exactly two interior affine Pauli visits and two affine Pauli endpoints has this form on the subspace fixing the endpoints' two clean-X syndromes. Its dimension satisfies $`r\ge p-2k`$. Absorb the endpoint target signs, which are one character, by multiplying P and Q by Z on one analytical borrowed qubit and multiplying B by the conjugating X character there. This keeps the support fraction and changes q by one. Exact clean-sector preservation and the stated row conjugation give the displayed contract directly. When false return is approximate, the observable error delta must also include any false-sector leakage.

For $`g\ne0`$ write $`P_g=B_gPB_g^\dagger`$ and $`Q_g=B_gQB_g^\dagger`$. The projective laws give the exact conjugation identity

```math
B_g^\dagger S_xB_g\ \propto\ S_{g+x}S_g^\dagger.
```

Both sides are proportional to $`B_{g+x}C_xB_g^\dagger`$. With exact false silence, conjugation by S_x fixes every $`Q_g`$ for $`g\ne0,x`$, and sends $`Q_x`$ to $`P_x`$. For approximate false silence these comparisons have errors at most $`2\delta`$ and $`\delta`$, respectively.

Set $`U_g=\exp(i\pi Q_g/2)`$ and $`V_g=\exp(i\pi P_g/2)`$. Encode a subset W of nonzero labels by the normalized Choi vector of the ordered product $`U_W=\prod_{g\in W}U_g`$. No such encoding or measurement is supplied to the compiler. On its $`2q`$ analytical qubits measure the effect $`(I-\mathrm{Re}\mathcal T_x)/2`$, where $`\mathcal T_x(K)=S_xKS_x^\dagger`$ is a unitary on operator coordinates. If x is absent the detection probability is zero. If x is present, only U_x changes, to V_x, and all other factors cancel in the trace. With normalized trace $`\tau`$, the detection probability is

```math
\frac{1-\mathrm{Re}\tau(U_x^\dagger V_x)}2
=\frac14\bigl(\|\Pi_P-\Pi_Q\|_{2,\tau}^2+\|P-Q\|_{2,\tau}^2\bigr)
\ge\beta.
```

Indeed $`\mathrm{Re}\tau(QP)\le -2^{-k}+\eta 2^{-k}`$ by $`Q(P+Q)`$, and $`\tau(\Pi_P\Pi_Q)\le2^{-k}`$. The squared-distance expression also applies to marked comparisons specified by Frobenius contrast instead of near negation. Under independent uniform subset bits, the detection event conveys at least beta/2 bits of mutual information about membership of x. Summing the single-bit mutual informations is bounded by the Holevo information $`2q`$, by independence, the conditional mutual-information chain rule and data processing. This proves the stated width bound.

For nonzero delta, each wrong conjugation changes its exponential by at most $`\pi\delta`$ by Duhamel's formula. Its difference is supported on the union of two rank-$`2^{q-k}`$ projectors, so the normalized trace-norm error per factor is at most $`2\pi\delta 2^{-k}`$. The product hybrid keeps the mark's exact one-factor replacement and bounds detection error by $`\pi m\delta 2^{-k}`$. The displayed delta threshold bounds absent detection above by beta/4 and present detection from below by 3beta/4. A detection event thus has posterior probability at least 3/4 of membership and occurs with probability at least 3beta/8. The same information bound gives the asserted constant-fraction width. At k=0 with an exact flip, the $`2^{m-1}`$ subset products are exactly Hilbert–Schmidt orthogonal, as their conjugation signs differ at every membership discrepancy. This sharpens the zero-clean constant to $`q\ge(m-1)/2`$ in that case.

<a id="two-visit-translation-rank"></a>

### Exact signed tables have a derivative-rank bound

Suppose instead $`S_xPS_x^\dagger=(-1)^{f(x)}P`$ exactly, with f(0)=0. Then for every y the same group identity gives

```math
S_xP_yS_x^\dagger=(-1)^{f(x+y)+f(y)}P_y.
```

Let R be the binary rank of the displayed derivative table. Select R independent columns. Products of the actual conjugations by S_x can flip any one selected observable and fix the others, by solving the linear system for that sign pattern. The preceding subset encoding on these columns now gives $`q\ge R/(4\cdot2^k)`$, or $`q\ge R/2`$ when k=0 and the products are orthogonal. The proof needs exact signed conjugation. A general compressed signed gap before returned amplification does not imply it.

The translation matrix with entries $`f(x+y)`$ is convolution by f in the binary group algebra $`\mathbb F_2[t_1,\ldots,t_r]/(t_i^2)`$, where $`t_i`$ is one plus the ith group generator. Its rank drops by at most one on subtracting row zero from all rows. For odd total truth parity the constant algebra coefficient is one, so f is a unit and R=m−1. More generally choose a nonzero algebra monomial $`t_S`$ of smallest degree d. Multiply inputs using only complementary variables and project the products onto monomials containing S. In the complement this is multiplication by a unit, because the coefficient of $`t_S`$ has constant term one. Thus $`R\ge2^{r-d}-1`$.

For uniformly random truth bits, the algebra coefficients are independently uniform since the change of coordinates is invertible. The probability that degrees zero through two all vanish is $`2^{-1-r-\binom r2}`$. A union bound over at most $`2^{O(p(k+1))}`$ affine endpoint cosets of codimension at most 2k therefore gives $`d\le2`$ on every available coset with probability tending to one, for $`k=p^{o(1)}`$. A constant sign normalization affects only the degree-r algebra monomial. Thus typical exact signed tables have the same exponential width order even after an adaptive endpoint restriction. This condition concerns the two-visit representation of an exact observable, rather than the four-factor scalar return of (6).

<a id="paid-two-visit-exact-table"></a>

## 9. A paid two-visit exact table at exponential width

Let x=(u,v) with $`|u|=k\le r`$, and put $`m'=2^{r-k}`$. Give the compiler k clean bits A and an arbitrary bank with an even independent Majorana Pauli list containing $`\Gamma_j`$ for all $`j\in\mathbb F_2^{r-k}`$ and spare axes. Its width is $`O(m')`$. Define B_x to translate A by u and map each indexed $`\Gamma_j`$ to $`\Gamma_{j+v}`$, fixing spare axes. These signed Clifford translations form a projective representation, since their prescribed action on the full Pauli frame determines each product up to scalar. Let W on A=a be the ordered Pauli $`W_a=X^{s(a)}Z^{t(a)}`$ satisfying

```math
W_a\Gamma_jW_a^\dagger=(-1)^{f(a,j)}\Gamma_j.
```

Any signature is possible because pairing with the full independent even frame is invertible. On A=0, $`B_x W B_x^\dagger W^\dagger`$ preserves A and conjugates $`\Gamma_0`$ with sign $`(-1)^{f(u,v)+f(0,0)}`$ on every arbitrary bank column. A fixed anticommuting Pauli corrects the constant sign if needed.

For a literal two-visit implementation with affine Pauli address controls, use r further arbitrary bits z and the fixed work word $`U=\sum_z|z\rangle\langle z|\otimes B_z`$. Then $`A_x=U X_z^x U^\dagger`$ acts on the bank by B_x times a monomial unitary on z, because $`B_{z+x}B_z^\dagger\propto B_x`$. The z coefficient can depend on x,z but acts trivially on the bank. It cancels in $`A_x W A_x^\dagger W^\dagger`$, where W is independent of z. The two actual address representations are A_x and $`W A_x^\dagger W^\dagger`$. They are fixed-unitary conjugates of the same affine Pauli representation, and include every incoming z and reference input.

To compile U, factor B_z into the affine A translation and r−k addressed Majorana coordinate translations. Each such translation consists of $`m'/2`$ disjoint Majorana swaps. On pair (a,b), $`\exp(-i\pi P_{ab}/4)`$ with $`P_{ab}=i\Gamma_a\Gamma_b`$ swaps the two axes with one minus sign and fixes the others. A fixed Pauli corrects the unwanted signs, and its controlled version is Clifford. The other controlled factors are exact:

```math
\mathrm{ctrl}_z\exp(-i\pi P_{ab}/4)
=T(P_{ab})T(Z_zP_{ab})^\dagger,
\qquad T(P)=e^{i\pi/8}e^{-i\pi P/8}.
```

The two T scalars cancel. The axes for a matching commute and are independent. A fixed work Clifford takes them to distinct Z targets, at $`O((m'+2)^2)`$ elementary Cliffords including its actual reverse. One private T band implements the first factors, and CNOTs from the selected z bit to the same targets realize the second private T band. Each controlled translation thus uses m' T gates, two T bands and no private helper. The affine A translation is Clifford. This proves common bounds for U and its actual inverse.

To implement W, use the [exact k-bit dirty indicator](NONUNIFORM_DIRTY_INDICATOR.md) toggling one bit of an arbitrary $`2^k`$-bit bank. Before and after that toggle, apply a fixed Clifford matrix shear adding the binary parity of all X masks whose corresponding bank bits are one. Actually reverse the indicator. The two old bank contributions cancel on every column, giving the X part of W with literal phase one. The Z part uses separately the same echo with controlled-Z shears, so its unknown offsets and phases also cancel. Emit the two factors in their prescribed order. This uses $`O(2^k)`$ T gates and dirty work, $`O(2^km')`$ elementary Cliffords and $`O(\log(k+2))`$ T-depth. For k=0 the fixed W is directly a Pauli.

The comparison and its actual inverse have common bounds

```math
\begin{aligned}
T&=O((r-k)m'+2^k),
&G&=O((r-k+1)(m'+2)^2+2^km'+r),\\
D_T&=O(r-k+\log(k+2)),
&b_{\rm dirty}&=O(m'+2^k+r),\quad a_{\rm clean}=k.
\end{aligned}
```

Every original address bit stays read-only at the completed Clifford boundaries. The inverse uses a constant-factor positive-S expansion. A fixed Clifford identifies the target Z with $`\Gamma_0`$. For returned XOR into an arbitrary separate target, conjugate a CZ between these two targets by the completed comparison, using its actual inverse, and multiply by the original CZ. Hadamards on the new target turn the resulting exact phase into XOR. All bank inputs, z, A and private helpers return on the entire initialized code, with no cost beyond constant factors in the ledger. For $`k\le r/2`$ the specified Majorana bank gives width $`\Theta(m')`$, matching the order of the restricted lower bound for typical exact tables before any endpoint loss. The width and native depth do not improve the general frame frontier.

<a id="two-visit-returned-word-packing"></a>

## 10. Fixed-error packing for completely returned two-visit words

Suppose an actual read-only comparison word has work fibers

```math
V_x=\lambda_x P_x A' B_x C_x A R_x,
```

where A,A' are arbitrary fixed work unitaries, B,C projective representations of the additive address group, and P,R affine Pauli endpoints with binary-linear labels after absorbing their zero-address factors. The phase $`\lambda_x`$ may be arbitrary. All fixed words can have arbitrary native depth and dirty width. In particular, a gatewise read-only word with at most three native T bands has this shape and two internal visits. The two internal visit representatives can be any fixed-unitary conjugates of affine Pauli families.

Let J append k clean ancillas and leave target t and every dirty wire arbitrary. Consider either a returned Boolean XOR word $`U_x=X_t^{f(x)}`$ or the deepest-corner frame $`U_x=(X_t Z_t)^{f(x)}`$, with the identity on other work and arbitrary refs. Assume the full initialized-isometry promise

```math
\|V_xJ-JU_x\|_{\mathrm{op}}\le\epsilon\qquad\text{for all }x.
```

The proof also works with unit-modulus target-only eigenphases whose ratio is $`(-1)^{f(x)}`$, on a fixed Pauli eigenbasis. It concerns complete returned words, not a raw signed-gap decision with arbitrary garbage.

Restrict to the linear subspace H on which the two endpoint Paulis preserve the clean code. Their clean-X labels impose at most 2k equations, so $`r=\dim H\ge p-2k`$. Adding a fixed affine address shift moves a chosen coset to zero if needed. Let

```math
\mathcal L=\{h\in H:\ x\mapsto f(x+h)+f(x)\text{ is affine on }H\}.
```

It is a subspace, since the derivative in a sum of directions is the sum of one derivative and the other shifted derivative. Then necessarily

```math
2^{-2k-1}\le |\mathcal L|2^{-r}
 +(1-|\mathcal L|2^{-r})(12\epsilon)^2.
```

In particular, with one actual compiler-clean wire and error at most $`2^{-6}`$, typical tables cannot have such a returned word irrespective of dirty width or fixed middle-word depth. This is not a lower bound on native words with three or more genuine internal visits, or on circuits that move the address physically.

### Lifted target parity removes scalar and projective phases

For a work operator K write $`K^{[2]}=K\otimes\overline K`$. Both B and C become ordinary representations under this lift, with the phases canceling. The lifted endpoint Paulis are ordinary representations too, and every lifted Pauli commutes with every other lifted Pauli, as scalar Pauli commutators cancel. No doubled physical word is emitted.

Let $`\sigma=X_t`$ for XOR and $`\sigma=Y_t`$ for the deepest-corner frame. Take the fixed mathematical code with two copies of the k clean zero bits and target parity $`\sigma\otimes\bar\sigma=-1`$. Let $`J'`$ be its insertion isometry. Its relative rank is $`\mu=2^{-2k-1}`$ because each eigenspace of the Pauli parity has half the dimension on the clean tensor code. Every lifted endpoint Pauli preserves this parity, even if its individual factor negates sigma, because the two signs multiply to one. Hence the endpoints preserve $`J'`$ for all x in H. The ideal lift $`U_x^{[2]}`$ restricted to this sector is precisely the scalar $`s(x)=(-1)^{f(x)}`$, for both target choices. For the XZ corner, $`X_tZ_t=-iY_t`$, so its tensor product with its complex conjugate is the lifted target parity itself.

The actual lifted initialized-isometry error is at most $`2\epsilon`$ by a two-term tensor hybrid; its restriction to the parity sector has the same bound. Define two fixed isometries

```math
E=A^{[2]}J',\qquad F=(A'^{[2]})^\dagger J'.
```

Write $`p_x,r_x`$ for the restricted lifted endpoints and set $`L_x=p_x^\dagger r_x^\dagger`$. Their representations commute, so L is ordinary. Removing the lifted endpoints and the fixed factors from the promised action proves

```math
\|C_x^{[2]}E-s(x)(B_x^{[2]})^\dagger F L_x\|\le2\epsilon.
```

The state spaces E and F need not coincide, so no additional marked-word approximation is charged. By changing the constant phase of F, assume s(0)=1, which does not alter affine-derivative directions.

### Returned word overlaps force an approximate character relation

For brevity suppress the superscript [2] on B,C. Fix h and define $`A_h=F^\dagger B_hF`$. With $`y=x+h`$, compare the two isometries in the preceding display on the pair x,y. The ordinary representation identities give

```math
\left\|E^\dagger C_h E-s(x)s(x+h)L_x^\dagger A_h L_{x+h}\right\|
\le4\epsilon.
```

The same comparison at x=0 gives $`E^\dagger C_h E`$ within $`4\epsilon`$ of $`s(h)A_h L_h`$. Since $`L_{x+h}=L_xL_h`$, subtraction and right multiplication give

```math
\left\|L_x^\dagger A_h L_x-d_h(x) A_h\right\|\le8\epsilon,
\qquad d_h(x)=s(x)s(x+h)s(h).
```

If h is outside $`\mathcal L`$, d_h is not a character of H. Choose x,y for which $`d_h(x+y)\ne d_h(x)d_h(y)`$. Compose the comparison twice at x and y and compare it to the one at x+y. Since the signs differ and conjugation is isometric, $`2\|A_h\|\le24\epsilon`$. Thus $`\|A_h\|\le12\epsilon`$ for every h outside L.

### Translated-code packing gives the fixed-error bound

For a transversal of the $`M=2^r/|\mathcal L|`$ cosets form the projector $`\Pi_x=B_xFF^\dagger B_x^\dagger`$. In total Hilbert dimension D each has rank $`R=\mu D`$. Any two distinct x,y in the transversal have a difference outside L. Their joint trace overlap is the squared Frobenius norm of $`A_{x+y}`$, at most $`R(12\epsilon)^2`$. Positivity and trace Cauchy--Schwarz applied to the sum give

```math
(MR)^2/D\le\mathrm{Tr}\Bigl(\sum_x\Pi_x\Bigr)^2
\le MR+M(M-1)R(12\epsilon)^2.
```

Division by $`M^2R`$ proves the claimed lower bound. This proof uses arbitrary fixed E,F and the actual returned target action. It makes no assertion about approximate-sign witnesses before complete return.

For a fixed large H and nonzero h the differences $`f(x+h)+f(x)`$ of a uniformly random table are independent uniform bits over the $`2^{r-1}`$ input pairs. The probability they agree with any affine value function is at most $`2^{r+1-2^{r-1}}`$. A union bound over h and all affine endpoint cosets of codimension at most 2k shows that, for fixed k, almost all tables have $`\mathcal L=\{0\}`$ on every possible retained H. At k=1 and error $`2^{-6}`$, $`(12\epsilon)^2=9/256<1/8$. The inequality requires $`[H:\mathcal L]\le8`$ because the index is a power of two. The theorem applies directly to tested deepest-corner frame columns and does not require a separate frame-to-XOR inverse call. Four physical read-only T bands generally leave three genuinely interior projective visits. Grouping two there requires a further closure or commutation argument, absent from the proof above.

<a id="biased-two-visit-approximate-row"></a>

## 11. Fixed false error allows a polynomial-width two-visit row

The exponential width in §8 needs exact or sufficiently accurate false action. For a fixed rational $`0\lt\delta\lt1/2`$, choose a power of two $`M=O((p+1)/\delta^2)`$ and a multiset of M binary characters $`\chi_j`$ on $`\mathbb F_2^p`$ satisfying

```math
\beta_x=\frac1M\sum_j\chi_j(x),\qquad \beta_0=1,\qquad
|\beta_x|\le\delta\quad(x\ne0).
```

Independent uniform characters and the independent-sign tail bound give simultaneous existence by a union bound. Finite search of the integer sums supplies the table. Use M independent Pauli Majoranas on an arbitrary $`M/2`$-qubit bank, and define

```math
R=M^{-1/2}\sum_j\Gamma_j,\qquad
R_x=B_xRB_x^\dagger=M^{-1/2}\sum_j\chi_j(x)\Gamma_j,
```

where the projective Pauli family B_x supplies the indicated signs. Its coordinate generators are the Pauli solutions to the binary dual-pairing equations of §9. Their controlled gates are Clifford and cost $`O(pM)`$ elementary Cliffords, with the address read-only. Both R and R_x are Hermitian involutions. Set $`C_x=R B_x^\dagger R`$, so $`S_x=B_xC_x=R_xR`$. Anticommutation gives the full-space identity

```math
S_x R S_x^\dagger=-R+2\beta_xR_x.
```

Choose a fixed work Clifford identifying $`\Gamma_0`$ with the target Z. Implement a fixed loader F with $`F\Gamma_0 F^\dagger=R`$. At each of $`\log_2 M`$ stages, split every current supported axis into an equal sum of itself and one unused partner. The disjoint pairs have independent commuting axes $`i\Gamma_a\Gamma_b`$ and each split uses one Pauli T rotation with its sign chosen to give the positive coefficient. Its common fixed Clifford frame costs $`O(M^2)`$ elementary gates per stage and puts all the rotations in one physical private T band. This exact F uses M−1 T gates and no helper. Use the actual $`F\Gamma_0 F^\dagger`$ to implement R. The chirality Pauli D of the even frame, the normalized product of all its generators, anticommutes with each generator. Thus the emitted comparison $`F^\dagger D S_x F`$ exactly negates $`\Gamma_0`$ at x=0 and preserves it to operator error at most $`2\delta`$ on every false x. Its actual inverse has the same observable errors. The common ledger is $`T=O(M)`$, $`G=O(M^2\log M+pM)`$, $`D_T=O(\log M)`$, $`b_{\rm dirty}=O(M)`$, and zero clean bits. Both visits, all fixed loader occurrences and their actual adjoints are charged. All arbitrary dirty and reference inputs obey the operator comparison.

This per-row error cannot be substituted into a dense independent product. Let $`n_\chi`$ count the multiplicity of character chi. Parseval gives

```math
\sum_{x\ne0}|\beta_x|^2
=2^p\sum_\chi(n_\chi/M)^2-1\ge 2^p/M-1.
```

For the literal target-Z echo the false row's full error equals $`2|\beta_x|`$. At any fixed address, every nonzero mismatch occurs in one false row. If an independent tensor product of row words has full error $`\tau\lt1`$ for every occupation subset, each row has eigenangles at most $`\alpha=2\arcsin(\tau/2)`$ in absolute value after factoring its ideal sign. Choose one largest absolute eigenangle per row. The positive chosen angles must sum to at most alpha, since a first partial sum exceeding alpha would still be below pi and would give a violating subset. The negative angles obey the same bound. Hence the sum of squared row errors is bounded by $`4\alpha^2`$, independent of the number of rows. The displayed Parseval estimate requires $`M=\Omega_\tau(2^p)`$ in that independent schedule. To use this approximate gap for dense lookup one must change the joint native schedule or the coherence of the false responses. A polynomial-width standalone probe at fixed tolerance does not supply that operation.
