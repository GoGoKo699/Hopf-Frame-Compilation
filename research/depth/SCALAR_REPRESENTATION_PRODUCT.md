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
