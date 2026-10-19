# Original middle pairs as a finite physical insertion

[Occupied pair center](OCCUPIED_LOGICAL_TRANSPORT.md#4-the-exact-remaining-diagonal-pair-center-question) ·
[No-pair closure](GROWING_DISPLACEMENT_NATIVE_ORBITS.md) ·
[Full-space finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md)

The existing pair-center reduction already records the sectorwise
Clifford action and the relations created by a quotient. The reduction
here keeps the original physical axes instead. It makes the remaining
bounded-pair question a precise stability question for the completed
no-pair theorem, without claiming that stability has been proved.

An original clean-transition pair is not one of the logical pair orbits
classified in the no-pair theorem. The latter theorem does not include
an extra insertion merely because its support is finite before clean
compression. No new Boolean-degree bound or frame frontier follows here.

## 1. Keep the actual middle axes and their signs

Use the actual one-clean three-band word

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,
\qquad P_A=K_AZ_cK_A^\dagger,
\qquad J=E_AP_AE_A^\dagger.
```

Every E is an independent commuting signed Pauli T/T-dagger list.
Assume complete clean return, arbitrary dirty and reference inputs,
and an occupied Boolean sign, with its one fixed common scalar retained.
The actual clean comparison is

```math
E_BJE_B^\dagger=E_CP_CE_C^\dagger,
\qquad P_C=K_CZ_cK_C^\dagger.
```

Native matching gives signed singleton and pair orbits of J on the
original B axes. Let the m pairs have actual Hermitian axes B_j,D_j,
original orientations epsilon_j,delta_j in {1,-1}, and signs s_j in
{1,-1} such that

```math
JB_jJ=s_jD_j,\qquad JD_jJ=s_jB_j.
```

Define the following literal operators, with all factors commuting:

```math
\begin{aligned}
F_j&=T_{B_j}^{\epsilon_j}T_{D_j}^{s_j\epsilon_j},\\
C_j&=T_{D_j}^{\delta_j-s_j\epsilon_j},\\
F&=\prod_{j=1}^mF_j,\qquad C=\prod_{j=1}^mC_j.
\end{aligned}                                                    \qquad(1)
```

Each C_j is I or a signed Pauli S/S-dagger phase. In particular C is
an actual Clifford, not a Clifford specified only modulo scalar.
The identity

```math
T_{-P}^{a}=\omega^aT_P^{-a},\qquad
\omega=e^{i\pi/4},\qquad a\in\mathbb Z
```

proves JF_jJ=F_j. When s_j=-1, the two factors produce opposite
powers of omega, which cancel exactly. For s_j=1, J simply exchanges
the equal orientations. Thus

```math
[F,J]=0,\qquad E_B=CFE_0,                          \qquad(2)
```

where E_0 consists of the original singleton axes and their original
orientations. F commutes each axis of E_0. Its original 2m axes are
independent, so its Pauli-support span has dimension at most 2m.

Absorb only the actual Clifford C into the output data:

```math
K_C'=C^\dagger K_C,\qquad E_C'=C^\dagger E_CC,
\qquad P_C'=C^\dagger P_CC.
```

Then the complete operator is exactly

```math
V=(K_C')^\dagger(E_C')^\dagger F E_0E_AK_A,          \qquad(3)
```

and its clean comparison is

```math
E_0JE_0^\dagger=E_C'P_C'(E_C')^\dagger.            \qquad(4)
```

Both endpoint lists remain physically independent. Every remaining
middle axis is an original singleton under J, so (4) is an actual
no-middle-pair clean comparison. No central sector was selected and no
axis was replaced by its quotient image.

The pair differences R_j=B_jD_j remain the globally transported center
of the existing occupied-pair theorem. C fixes them pointwise. Equations
(1)--(4) supplement that theorem with a finite insertion in the original
physical coordinates; they do not replace its central normalization.

## 2. What the bounded-pair question now requires

Let V_0 denote the word (3) with F deleted. It still returns the clean
qubit, because of (4). But

```math
V=V_0Q_F,\qquad
Q_F=K_A^\dagger E_A^\dagger F E_AK_A,
\qquad [Q_F,Z_c]=0.                                \qquad(5)
```

The factor Q_F need not be diagonal or a Boolean sign. Consequently
V being an occupied Boolean sign does not prove that V_0 is one.
The no-pair degree theorem cannot be applied to V_0 with that missing
hypothesis silently supplied.

A sufficient new lemma would be the following, presently **unproved**:
for every fixed m, the occupied Boolean degree of the actual words (3)
is bounded by a function of m alone, under (4), independent endpoint
lists, and the literal native insertion (1) that commutes J and each
E_0 axis. The bound must be uniform in endpoint clean-active ranks,
dirty width, and all central controls. Such a lemma would close fixed
physical-pair count. It would still leave growing physical-pair count
and the higher-layer asymptotic problem open.

The retained full-space finite-insertion theorem does not establish this
lemma. It requires full-space logical-probe comparisons after the clean
restriction has been justified. Here only the initialized clean sector
has the specified Boolean action; the other sector is retained and may
have a different action. Commuting with a non-Pauli clean image does not
turn an occupied comparison into a full-space Clifford identity.

Nor is it valid to pull F through E_A and keep its support bound 2m.
A native band sends a Pauli to a Clifford, generally with growing Pauli
support. The example below shows this failure already for one physical
pair. It is a failure of that relocation argument, not an obstruction
to the proposed initialized-code lemma.

One finite description does survive the pullback. Expanding the original
2m commuting T phases gives at most 4^m Pauli products, with their fixed
cyclotomic coefficients. Therefore

```math
E_A^\dagger F E_A=\sum_{a\in\mathbb F_2^{2m}}c_aH_a,\qquad
H_a=E_A^\dagger P_aE_A,                             \qquad(6)
```

where each H_a is an actual Hermitian Clifford. The complete sum preserves
P_A. Its individual terms need not preserve P_A. Equation (6) is an
operator identity, not a free linear-combination circuit or a replacement
for the finite Pauli-support hypothesis of the existing theorem.

## 3. One original pair can acquire growing Pauli support

Let k be a positive multiple of four. On k physical wires put

```math
E_A=E_C=T^{\otimes k},\qquad
P=X^{\otimes k},\qquad Q=Y^{\otimes k},\qquad
F=T_PT_Q.
```

P and Q are independent commuting Hermitian Pauli axes, so F is an
actual two-axis native band. Let K be the GHZ Clifford consisting,
in time order, of H on wire zero followed by CNOT from that wire to
every other wire. Only wire zero is initialized. Then

```math
KZ_0K^\dagger=P,\qquad
J=E_APE_A^\dagger
 =\bigotimes_{i=1}^k\frac{X_i+Y_i}{\sqrt2}.
```

J exchanges P and Q with positive signs. Thus the source

```math
V_k=K^\dagger E_A^\dagger F E_AK                  \qquad(7)
```

has exactly one original middle pair, both clean-active ranks equal to
k, and exact clean return. The other k-1 wires are arbitrary data.
The family is not being asserted to have Boolean diagonal action.
In particular it cannot disprove a lemma that uses that additional
hypothesis.

Write alpha=(1+omega)/2 and beta=(1-omega)/2. Since k is a multiple
of four, PQ=Z^{\otimes k}; hence

```math
E_A^\dagger FE_A
 =\alpha^2I+\beta^2Z^{\otimes k}
  +\alpha\beta\,2^{1-k/2}
    \sum_{|S|\ {\rm even}}X_{S^c}Y_S.               \qquad(8)
```

The 2^{k-1} displayed X/Y Paulis have distinct labels and nonzero
coefficients. They do not collide with I or the all-Z term.
Their differences span the even-weight Z subspace, and one X/Y label
adds an independent direction. The complete Pauli-support span has
therefore dimension k, although F's original support span has dimension
two. This is actual native support growth, not a synthetic Gram model.

A concrete coordinate word also supplies the ledger. Before K, apply
CNOT from each wire other than one to wire one, and call the resulting
Clifford U. In operator order U=K times that CNOT word. For k divisible
by four,

```math
UZ_0U^\dagger=P,\qquad UZ_1U^\dagger=Q,
\qquad F=U(T_0T_1)U^\dagger.
```

The U word contains 2k-2 CNOTs and one H. Thus (7) and its actual
reversed inverse each use 2k+2 T/T-dagger gates, at most three T layers,
6k-2 elementary Cliffords, one clean wire, and no additional work.
All column identities hold literally, including the unused clean sector
and arbitrary references. The zero-leakage statement concerns the sole
initialized wire, not identity return of the data acted on by V_k.

## 4. Native accounting and next question

For a general original source, F has 2m T/T-dagger gates in one native
band on the existing wires. Its inverse has the same T-count and reverses
every orientation. C can be emitted in the original B coordinates using
at most m S/S-dagger gates, or at most 3m positive S gates, with both
coordinate words charged. Fixed Clifford synthesis gives O(Q^2)
elementary gates on the same Q wires.

Equations (1)--(4) preserve the original total T-count and at most three
T bands: F and E_0 share the original independent B coordinates.
Explicitly emitting the Clifford correction and its inverse adds at most
O(Q^2) elementary Cliffords to each respective original forward/inverse
ledger. There is no additional clean or dirty allocation, measurement,
postselection, approximation error, or discarded scalar. The normal form
is a degree-proof interface, not a faster compiler.

The [exact native checks](../../tests/test_original_middle_pair_insertions.py)
verify all signed swap/orientation cases on a small actual source,
including the correction, normalized clean comparison, complete columns,
and actual inverse. They also verify (8) for k=4 and its paid source.

The next lower question is whether the literal global Boolean condition
controls (3) uniformly at fixed m while retaining the finite Clifford
sum (6), rather than assuming it has finite Pauli support after clean
compression. Neither this note nor the existing sector reduction proves
that step. R96, the no-pair degree theorem, the fixed endpoint-rank theorem,
and the remaining unrestricted frontier retain their previous scopes.
