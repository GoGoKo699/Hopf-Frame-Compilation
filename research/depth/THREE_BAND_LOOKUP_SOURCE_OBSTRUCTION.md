# Three native T bands cannot supply complete lookup at growing address size

[Two-block obstruction](CLIFFORD_INTERLEAVED_LOOKUP.md) ·
[One-phase trace gap](FOUR_LAYER_OPERATOR_RIGIDITY.md#2-the-one-phase-gap-with-arbitrary-hermitian-clifford-boundaries) ·
[Three-layer exactness](MOD3_ZETA_TRANSPORT.md#three-band-transport-exactness)

For at least seven address bits, no zero-clean native source with at most
three physical T layers has complete lookup as a Pauli derivative. The
statement permits arbitrary dirty width, unrestricted Clifford interlayers,
and intermediate address/bank/work movement. Its uniform operator-error
floor is $`\tau_* =\sin(\pi/8)/(2\sqrt2)`$.

The new step is a structural theorem at five address bits: the final two
pure Pauli bands, after absorbing the original Clifford walls, must preserve
the address-Z algebra pointwise in any exact candidate. A low-rank point
reflection and the retained one-phase trace gap force this preservation;
it is not assumed. The existing degree-six theorem then closes the exact
three-band case at seven bits.

This is a constant source-depth lower bound, not a growing unrestricted
native-depth lower bound. It does not exclude four or more source layers,
initialized clean-work promises, or every direct implementation of a query
that is not presented as the stipulated Pauli derivative. R96 is unchanged.

## 1. Full-space model and claims

Let $`p\ge1`$, $`m=2^p`$, and let c contain one arbitrary coefficient
bit for every subset of the p address positions. On any number of
additional arbitrary helper qubits h, define

```math
Q|x,c,h\rangle=(-1)^{f(x,c)}|x,c,h\rangle,
\qquad f(x,c)=\sum_{S\subseteq[p]}c_S\prod_{i\in S}x_i.
                                                               \qquad(1)
```

Every input and every external reference is included. The fixed bank
Möbius Clifford converts this to ordinary addressed lookup; that change
of coordinates changes no native T depth. Let P be any fixed Hermitian
Pauli on the whole physical space. A source V need not separately return
its helpers: the target $`VPV^\dagger=Q`$ includes their identity action.

The native alphabet is H, S, CNOT and T/T-dagger. Each physical T layer
acts on distinct wires. Write $`q=p+m+|h|`$ and $`D=2^q`$.
There is no bound on q or on the number of Clifford gates.

**Three-band obstruction.** For $`p\ge7`$ and every native V of
T-depth at most three,

```math
\|VPV^\dagger-Q\|\ge\tau_*,\qquad
\tau_* =\frac{\sin(\pi/8)}{2\sqrt2}.                  \qquad(2)
```

The exact structural result used to prove (2) holds already for
$`p\ge5`$. Its threshold is distinct from the final degree threshold.
The cases p=3 through p=6 are not all resolved by this note.

## 2. Normalize the actual physical bands

Pad empty bands if necessary and write
$`V=C_3D_3C_2D_2C_1D_1C_0`$, with each D a physical T/T-dagger
band and every C Clifford. Define

```math
\begin{aligned}
E_3&=C_3D_3C_3^\dagger,\\
E_2&=(C_3C_2)D_2(C_3C_2)^\dagger,\\
K&=(C_3C_2C_1)D_1(C_0PC_0^\dagger)D_1^\dagger
             (C_3C_2C_1)^\dagger.
\end{aligned}                                                   \qquad(3)
```

The actual prefix $`U_0=(C_3C_2C_1)D_1C_0`$ is third-level and
satisfies $`K=U_0PU_0^\dagger`$ and $`V=E_3E_2U_0`$. A physical
band sends a Hermitian Pauli to a Hermitian Clifford, so K is literally
a Hermitian Clifford. The derivative is exactly

```math
VPV^\dagger=E_3E_2K E_2^\dagger E_3^\dagger.           \qquad(4)
```

Each E is a pure band of commuting Pauli rotations: its nontrivial odd
axes are independent binary Pauli labels, obtained by Clifford-conjugating
distinct physical Z axes. Their signs and the native scalar phases remain
in the word. This use of independence does not apply to arbitrary lists
of dependent phase gadgets.

For an exact candidate put

```math
G=E_3^\dagger,\qquad B=E_2,\qquad
M=GQG^\dagger=BKB^\dagger.                            \qquad(5)
```

Both G and B, and their actual inverses, send every Pauli to a Clifford.
This follows from their explicit physical-band forms. No general closure
of the third Clifford-hierarchy level under inversion is assumed.

## 3. A rank threshold from the one-phase trace gap

The retained [one-phase lemma](FOUR_LAYER_OPERATOR_RIGIDITY.md#2-the-one-phase-gap-with-arbitrary-hermitian-clifford-boundaries)
states that for Hermitian Cliffords L,H and one pure physical Pauli band B,
with $`R=BHB^\dagger`$,

```math
L=R\quad\text{or}\quad
D^{-1}\mathrm{tr}(LR)\le1-\Delta,\qquad
\Delta=\frac{1-1/\sqrt2}{2}.                          \qquad(6)
```

It allows arbitrary Hermitian Clifford boundaries, not merely one-band
Pauli images. The overlap is real. For an unequal pair,

```math
2D\Delta\le\|L-R\|_{\rm HS}^2
       \le4\,\mathrm{rank}(L-R),
\qquad
\mathrm{rank}(L-R)\ge\frac{\Delta D}{2}.              \qquad(7)
```

The upper inequality uses $`\|L-R\|\le2`$. Thus equality is forced
whenever another argument gives rank at most $`2D/m`$ and
$`m\gt4/\Delta`$. Since $`4/\Delta=16+8\sqrt2\lt32`$,
p at least five suffices. No trace gap for two arbitrary two-band Pauli
images, and no general bounded-depth rank gap, is being asserted.

## 4. Every good Pauli commutes with the dressed address projectors

Call a Pauli label good for M when its conjugate by M is Clifford, and
write $`\mathrm{Good}(M)`$ for this binary space. It is a linear space:
conjugation preserves products and Cliffords form a group. Let
$`L_B`$ be the isotropic span of the independent odd axes of B.
For every Pauli R in its full centralizer, B commutes with R and
$`MRM=B(KRK^\dagger)B^\dagger`$ is Clifford. Therefore

```math
L_B^\perp\subseteq\mathrm{Good}(M),\qquad
\mathrm{Good}(M)^\perp\subseteq L_B.                  \qquad(8)
```

Fix an address a. The complete independent coefficient bank supplies a
binary mask s whose coefficient polynomial is the point function
$`[x=a]`$. Include the constant coefficient when it is needed. Let
$`\Pi_a=|a\rangle\langle a|\otimes I`$ and define

```math
A_s=GX_c^sG^\dagger,\qquad
S_a=G(I-2\Pi_a)G^\dagger.                            \qquad(9)
```

The first operator is a Hermitian Clifford, even when G mixes address
and bank. Direct conjugation of (1) gives

```math
A_sMA_s=MS_a,\qquad [M,S_a]=[A_s,S_a]=0,\qquad
\mathrm{rank}(I-S_a)=D/m.                             \qquad(10)
```

Take a Hermitian Pauli R with a good label, and set
$`R_s=A_sRA_s`$, another Hermitian Pauli. Compare

```math
L=A_s(MRM)A_s,\qquad
F=MR_sM
 =B\bigl[K(B^\dagger R_sB)K^\dagger\bigr]B^\dagger.
                                                               \qquad(11)
```

L is a Hermitian Clifford by goodness. The bracket in F is also a
Hermitian Clifford, because the actual inverse band $`B^\dagger`$
sends Paulis to Cliffords. These are exactly the boundaries required
by (6), with no hypothesis on the middle Pauli's goodness.
Equation (10) gives the exact difference

```math
L-F=M(S_aR_sS_a-R_s)M,\qquad
\mathrm{rank}(L-F)\le2\,\mathrm{rank}(I-S_a)=2D/m.      \qquad(12)
```

For p at least five, (7) forces $`L=F`$. Hence R_s commutes with
S_a, and so does R, since A_s commutes with S_a. This works for every
address a. Taking linear combinations of the corresponding projectors
proves that every good Pauli commutes with every

```math
J_i=GZ_{x_i}G^\dagger.                                \qquad(13)
```

The earlier bank-independent-final argument used ordinary coefficient
fibers to obtain this commutation. Equations (9)–(12) instead retain the
dressed Clifford bank translations, so they permit arbitrary final-band
bank coupling and intermediate address movement.

## 5. Independent native axes force both later bands to hold the key

A Pauli commuting with J_i commutes with each nonzero Pauli term in its
expansion: conjugating by a Pauli only changes the individual signs.
Therefore every support label of J_i lies in
$`\mathrm{Good}(M)^\perp\subseteq L_B`$.

For a pure physical band G acting on a Hermitian Pauli P, let
$`L_P`$ span the axes anticommuting with P. Independence makes the
Pauli support of $`GPG^\dagger`$ exactly $`P+L_P`$, with every
coefficient nonzero. If an axis a is active, that support contains P
and P+a, and they anticommute. An isotropic space cannot contain this
pair. Applied to J_i, whose entire support lies in the isotropic space
$`L_B`$, this proves that no axis of G is active on $`Z_{x_i}`$.
Thus, literally,

```math
GZ_{x_i}G^\dagger=Z_{x_i},\qquad
Z_{x_i}\in L_B,
\qquad [G,Z_{x_i}]=[B,Z_{x_i}]=0.                     \qquad(14)
```

This is the promised p-at-least-five structural theorem. It concerns
the pure bands in (3), not each original Clifford subword before their
absorption. Every original Clifford network and every dirty input were
allowed in deriving it.

For completeness, these pointwise commutation relations give exactly
the key-preserving normal form needed below. Each axis has only I/Z
factors on the address. Its bank/work restriction is therefore a member
of a commuting Pauli family. A fixed Clifford on bank/work alone
simultaneously diagonalizes that family. Hence each of $`E_3=G^\dagger`$
and $`E_2=B`$ has a decomposition into a diagonal third-level phase
between Clifford walls fixing every address Z. Writing these decompositions as $`E_j=F_jD'_jF_j^\dagger`$,
combining adjacent walls puts the actual source into

```math
V=C'_3D'_3C'_2D'_2U',\qquad\begin{cases}
C'_3=F_3,\quad C'_2=F_3^\dagger F_2,\\
U'=F_2^\dagger U_0\in\mathcal C_3.
\end{cases}                                               \qquad(15)
```

Both displayed C walls fix the key, both D phases are diagonal
third-level, and the actual first-band prefix supplies U'. No arbitrary
Hermitian Clifford is presumed to have a one-band Pauli preimage.
This is a proof normalization, not a
free synthesis replacing the original physical circuit.

## 6. Degree closure, approximation and scope

The retained [late-key theorem](CLIFFORD_INTERLEAVED_LOOKUP.md#three-block-late-key-mixing)
applies to (15), including genuinely non-Clifford bank/work fibers. It
bounds the Boolean degree of an exactly transported Pauli-label family
by six. The coefficient lookup (1) has label coordinate
$`\prod_{i=1}^p x_i`$ on $`Z_{c_{[p]}}`$, hence degree exactly p.
For p at least seven this contradicts exact transport.

To obtain (2), suppose the complete error were strictly below tau_*.
The retained [three-band exactness theorem](MOD3_ZETA_TRANSPORT.md#three-band-transport-exactness)
then forces $`VPV^\dagger=Q`$: folding the first Pauli image leaves
four T layers, R75 forces exact commutation with every target diagonal Z,
and distinct diagonal Boolean involutions have distance two. The exact
contradiction proves the stated uniform floor.

At least four source layers are consequently necessary below this floor
for p at least seven, regardless of dirty width or Clifford cost. This
does not establish a logarithmic lower bound, classify four-layer
sources, or rule out an initialized-sector construction. The result is
also not a claim that every general Hopf-frame compiler must use this
particular lookup-source representation.

## 7. Exact finite checks

`tests/test_three_band_lookup_source_obstruction.py` checks the dressed
point-reflection identity and its rank estimate with actual address/bank
mixing, the one-phase rank separation, the independent-axis support
argument, and the absorbed-Clifford normalization. Small complete matrices
retain literal phases and arbitrary spectator columns. These fixtures
check the proof ingredients; the dimension-independent implication and
its p thresholds follow from the argument above, not a finite search.
