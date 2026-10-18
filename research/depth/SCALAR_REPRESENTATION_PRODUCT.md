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
B_x=(I\otimes\gamma_0)A_x(I\otimes\gamma_0)
   =[(I\otimes\gamma_0)U](X_z^x\otimes I)
      [(I\otimes\gamma_0)U]^\dagger.                     \qquad(5)
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
invariant, not a resource saving. No lower bound against a better
native synthesis of U or the final scalar is asserted.

## 5. Consequence for the general depth question

A scalar obtained by commuting a returned XOR query with a fixed probe
can indeed be written as a bounded product of Pauli representations.
Discarding the native conjugators at that point loses necessary
information: (6) realizes every truth table within precisely that
abstract class. A successful invariant must use more than factor count,
even with literal gauges and exact full return.

The next question is whether the **shared native prefix structure** of
the representations obtained from a depth-d circuit imposes a relation
that the independently supplied, growing-cost sources in (5) do not
satisfy. This counterexample does not give a five-layer realization,
does not improve R96, and does not change the scoped R97 obstruction.

The [bounded checks](../../tests/test_scalar_representation_product.py)
verify the Pauli commutators, literal translation action and group law,
canonical-source identity, returned scalar on every work column, actual
inverse source word, and exact source gate counts in small instances.
The arbitrary-width conclusion is the analytic construction above.
