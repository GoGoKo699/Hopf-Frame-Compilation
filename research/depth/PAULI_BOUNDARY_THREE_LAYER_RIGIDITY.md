# Three layers with Pauli clean images have bounded Boolean degree

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Zero-clean three-layer theorem](THREE_LAYER_NATIVE_STRUCTURE.md) ·
[Projected native matching](ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation)

This is a restricted one-clean theorem. It allows arbitrary Clifford mixing
of address, output and dirty wires, but requires the actual clean image to
remain a Pauli stabilizer code at every native-band boundary. It does not
cover a non-Pauli intermediate code and does not improve the complete-frame
resource frontier.

The later [partial-diagonal extension](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#7-one-pauli-interior-clean-image-suffices)
covers either actual interior clean image being Pauli with a larger
degree bound. Its rank-one/rank-one reduction also reaches this theorem
after a literal repacking when both original interior images are
non-Pauli. The present theorem retains its stronger 2048 bound under
its own stated hypothesis.

## 1. Statement and physical hypothesis

Consider an actual Clifford+T circuit with at most three native T bands,
at most one clean qubit, and arbitrarily many dirty qubits. Suppose its
occupied action is projectively a Boolean sign operation
$`(-1)^{g(x)}`$, with complete clean and dirty return. Suppose also that
the image of its clean-zero space at every native-band boundary is the
positive eigenspace of a signed Pauli. Then

```math
\deg_{\mathbb F_2}g\le2048.                         \qquad(1)
```

Consequently a literal Boolean XOR query in this class obeys
$`\deg f\le2047`$: conjugate its target by Hadamards to obtain
$`(-1)^{zf(x)}`$. The constant is deliberately loose. The result is that
the bound is independent of arity, dirty width and Clifford cost.

The hypothesis concerns the actual physical circuit. A native band has
independent commuting physical axes. Conjugating a Pauli that anticommutes
with k active axes gives a flat expansion on $`2^k`$ distinct Paulis.
Thus a band carrying a Pauli clean image to another Pauli clean image
must have k=0: it commutes with that stabilizer. Clifford changes of code
coordinates identify each positive code with a full logical Hilbert space.
Each induced band is a product of signed T-angle phases about commuting
logical Pauli axes, with at most one binary relation, because only one
stabilizer direction was quotiented out.

A projected identity contributes its actual scalar. If two projected axes
coincide, their signed T or T-dagger factors combine to a scalar or an
ordinary Clifford rotation on that axis. Absorb that Clifford into an
adjacent logical Clifford, retaining the scalar. Conjugating the neighboring
band preserves its one-relation property. The absorbed logical Clifford
only changes adjacent axis coordinates and uses no new band. This is a
proof rearrangement, not an uncharged synthesis claim. We may therefore
assume all lists have distinct nonzero axes and any relation has length at
least three.

After these reductions the occupied diagonal operation is a full logical
operator of the form

```math
V=K_C^\dagger E_C^\dagger E_B E_AK_A.               \qquad(2)
```

Each E is one of the projected bands. Every scalar and every axis
orientation remains in (2). All subsequent logical identities hold on the
whole logical space, not merely on a nonlinear subset of computational
strings.

## 2. Ordinary middle-band pairs still peel

For an input diagonal Pauli label z put

```math
\begin{gathered}
R_z=K_A Z(z)K_A^\dagger,\qquad
J_z=E_AR_zE_A^\dagger,\\
S_z=E_BJ_zE_B^\dagger=E_CK_CZ(z)K_C^\dagger E_C^\dagger.
\end{gathered}                                      \qquad(3)
```

These are exact commuting representations; J and S are Hermitian
Cliffords. A commuting T-angle phase polynomial is C3 even when its axes
are dependent. Write the A-axis labels as a_i and put

```math
\lambda_i(z)=\langle a_i,R_z\rangle,\qquad
N_zv=\sum_i\lambda_i(z)\langle a_i,v\rangle a_i.
                                                        \qquad(4)
```

Then the symplectic action of J_z is $`I+N_z`$, and
$`N_zN_u=0`$ because the A axes commute. Peel globally silent A axes
to the input as diagonal C3 factors. Every remaining lambda_i is nonzero.

The projected matching theorem recovers every B axis individually except
possibly a unique length-four relation. In that exception it recovers the
invariant three-dimensional span W of the four relation axes. Call all
other B axes ordinary; take W=0 when the exception is absent.

The ordinary axes have only singleton or paired orbits under the J family.
An orbit of size at least four would give

```math
0=N_zN_u b=b+F_zb+F_ub+F_zF_ub
```

with four distinct ordinary axes. This would be a length-four relation,
which has already been isolated in W. For a pair $`\{b,b+r\}`$,
there is a nonzero linear swap character lambda with
$`N_zb=\lambda(z)r`$. The signed action has one common sign character
chi on the pair. Indeed an involutory signed swap has equal off-diagonal
signs, and every commuting nonswap has equal diagonal signs; the exact
group law then makes chi linear.

The difference r lies in the input maximal Pauli space
$`M_A=K_A\langle Z\rangle K_A^\dagger`$. To prove this without A-axis
independence, first choose z with lambda(z)=1 in (4), giving
$`r\in\mathrm{span}(A)`$. If B,D are the actual commuting
Hermitian representatives of the pair, their product R=BD is fixed
literally by every J_z: the two signs cancel even when the pair swaps.
Since E_A commutes R, every R_z commutes R. Thus $`r\in M_A`$.

The sign restriction from the zero-clean proof also survives. Normalize
the ordinary B rotations to positive T orientations, absorbing their
Clifford corrections into the C encoder with all phases retained. The
ordinary defect has displacement

```math
N_{B,z}v=\sum_b\chi_b(z)\langle b,v\rangle b.
```

The exceptional defect, if present, is Clifford: the total defect and
the ordinary defect are Clifford. It belongs to the Pauli algebra of W,
so its displacement image lies in W. Modulo W, (3) consequently gives

```math
F_{S_z}-I=N_z+N_{B,z}+N_{B,z}N_z\pmod W.             \qquad(5)
```

Singletons contribute zero to the last term. A pair contributes
$`\chi(z)\lambda(z)r\langle r,\cdot\rangle`$. In the exceptional
case the ordinary axes are independent modulo W, so these projected
rank-one tensors are independent. Without an exception the pair
differences have at most one relation and have no duplicates. Their
rank-one tensors are still independent: if a relation occurs, eliminate
one vector; its square tensor has a nonzero off-diagonal entry because
the relation has at least three members, whereas the other square tensors
are diagonal in that basis.

The symplectic action of the native C image is linear in z, including for
a dependent projected list. Equation (5) therefore forces each chi lambda
to be linear. Hence $`\chi=0`$ or $`\chi=\lambda`$.
Orient a pair equally in the first case and oppositely in the second.
Its full phase P then commutes every J_z literally. The native identities
$`T_{-B}=\omega T_B^\dagger`$ cancel the two scalars in the opposite
orientation case. Thus
$`K_A^\dagger E_A^\dagger P E_AK_A`$ is a diagonal factor that can
be peeled exactly.

Each peeled pair has an exponent, modulo eight, of one of the forms

```math
p=u+2(1-u)h,\qquad p=-u+2uh,                       \qquad(6)
```

where u is an affine input parity and h is a quadratic Boolean function.
This is the same parity-sector Clifford argument as in the zero-clean
proof: the conjugate of B is Hermitian Clifford and is diagonal on the
relevant parity sector. A quadratic extension from that sector remains
quadratic after returning to the original affine coordinates. Integer
lifting in (6) gives ordinary degree at most five modulo eight. No
simultaneous independence of the different pair parities is required.

## 3. The residual Pauli algebra has at most five quantum qubits

After peeling all ordinary pairs, consider an ordinary singleton b.
Equation (4) gives

```math
\sum_i\langle a_i,b\rangle\lambda_i(z)a_i=0
\quad\text{for every }z.                            \qquad(7)
```

If A is independent, b commutes every A axis. If A has one relation
with incidence vector r_A, the coefficient vector in (7) is either zero
or a multiple of r_A. Since each remaining lambda_i is nonzero, a
noncommuting singleton must anticommute with exactly the relation axes,
and those axes must all have the same activation character. Hence the
cross-pairing matrix between A and all ordinary singleton axes has rank
at most one. Both lists are internally isotropic. Their combined span
$`S_0`$ therefore has symplectic rank at most two.

Adjoin the possible exceptional middle span W, of dimension at most three,
and call the resulting space S. Adding one direction increases symplectic
rank by at most two, so

```math
\mathrm{rank}_{\mathrm{symp}}S\le8.             \qquad(8)
```

The residual source $`U=E_{B,\mathrm{remaining}}E_A`$ belongs to the
Pauli algebra of S. For each input R_z, the Clifford
$`(UR_zU^\dagger)R_z^\dagger`$ also belongs to that algebra, since
conjugation by a Pauli preserves it. It commutes every Pauli in
$`S^\perp`$. Symplecticity then gives

```math
\mathrm{im}(F_{UR_zU^\dagger}-I)\subseteq S.    \qquad(9)
```

Recover the active C axes from these images. For an independent active
subset its displacement image is the full active span. If the active
subset contains the unique even relation, its image contains every axis
outside the relation and every pairwise difference within the relation.
An axis is recovered individually whenever a probe activates it but not
another relation member. The only unrecovered direction occurs when all
relation members have one common nonzero activation character. Thus all
globally active C axes lie in S except for at most one additional
direction beta. This is the intrinsic-span argument of the
[two-layer theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md#31-the-one-relation-quotient-and-its-intrinsic-span);
it does not require S to be isotropic.

Peel globally silent C axes to the output as a diagonal C3 phase. The
remaining three bands all belong to the Pauli algebra of
$`\widehat S=S+\langle\beta\rangle`$, with

```math
\mathrm{rank}_{\mathrm{symp}}\widehat S\le10.   \qquad(10)
```

In canonical Pauli coordinates this algebra has arbitrary central Z
coordinates and at most five quantum qubits. This does not bound the
number of central coordinates or their quadratic ranks.

Write the residual occupied operation as
$`V_0=K_C^\dagger W_0K_A`$, where W_0 is the product of the three
remaining bands. V_0 is diagonal. It need not be real: all peeled pair
phases and diagonal C3 factors must still be restored.

## 4. Match the endpoint Cliffords with exact Pauli signs

We use a general diagonal extraction fact, without assuming a Boolean
residual. Suppose a diagonal unitary V maps every Pauli in a subspace A
to a Pauli under conjugation. There is a diagonal Clifford C with the
same signed action on every member of A.

To prove it, let T be the computational translation space of A and use
an invertible binary coordinate change making T the first s coordinates.
For each translation e_i in that space, diagonality and the Pauli image
give the literal ratio

```math
\frac{V(x+e_i)}{V(x)}=i^{a_i}(-1)^{\ell_i(x)},
\qquad a_i\in\mathbb Z/4\mathbb Z,                 \qquad(11)
```

with ell_i linear. Any constant sign is included in a_i. Applying the
translation twice gives
$`a_i\bmod2=\ell_i(e_i)`$. Equality of the two orders of translations
gives $`\ell_i(e_j)=\ell_j(e_i)`$. A phase polynomial with linear
coefficient $`2a_i`$ on x_i and pair coefficient
$`4\ell_i(e_j)`$ on x_i x_j consequently reproduces (11). Include
also its prescribed cross terms for j outside T, and choose the remaining
coefficients arbitrarily. It is a diagonal Clifford. Its differences
match every translation in T by the exact cocycle law, and therefore
match every signed Pauli in A. Returning through the binary coordinate
change preserves the Clifford property. This proof retains the possible
S phases and does not replace V by a real sign operation.

Apply this fact with
$`A=K_A^\dagger\widehat S^\perp K_A`$. Since W_0 commutes this
commutant in its own coordinates, V_0 has a Pauli action on A. Choose
the matching diagonal Clifford C and put

```math
J=K_CCK_A^\dagger,\qquad
V_1=C^\dagger V_0=K_A^\dagger J^\dagger W_0K_A.    \qquad(12)
```

J commutes every Pauli in $`\widehat S^\perp`$ with its exact sign.
The Pauli commutant theorem therefore puts J in the algebra of
$`\widehat S`$. V_1 is diagonal and belongs to the Pauli algebra of
$`K_A^\dagger\widehat SK_A`$, so its support lies in the isotropic
space corresponding to
$`\mathcal T=\widehat S\cap M_A`$. This is an exact full-logical-space
operator statement.

## 5. A finite set of low-degree phase features

Put the Pauli algebra of $`\widehat S`$ in canonical form: central
Z bits z, k quantum qubits with $`k\le5`$, and an identity factor on
any remaining wires. The central directions are the radical of
$`\widehat S`$. Take a basis of the central part of $`\mathcal T`$
as some of these central coordinates. Its core projections form an
isotropic space of dimension at most k. A fixed core Clifford makes
them Z-type. The eigenvalue bit of each such generator differs from its
core Z bit by one central parity; at most k such parity offsets occur.

These are choices of Pauli bases, not a nonlinear change of the target
variables. Pauli-basis uniqueness gives

```math
\mathrm{Alg}(K_A^\dagger\widehat S K_A)
\cap\mathrm{Diag}
=\mathrm{span}\{Z(t):Z(t)\in
K_A^\dagger\mathcal T K_A\}.                       \qquad(12a)
```

The phase therefore factors through the eigenvalue bits of these input
Z parities, which are affine functions of the original computational
input. Evaluate that same phase by the following linear section in the
canonical coordinates: prescribe the central part of the target parity
bits, set the other central coordinates to zero, set each selected core
coordinate to its prescribed target bit plus its central parity offset,
and set remaining core coordinates to zero. The operator belongs to the
algebra of $`\mathcal T`$, so this section gives exactly its value.
Every substitution in the section is affine. Include the at-most-k
central offsets and core eigenvalue bits among the features below. The
canonical Clifford need not preserve the rest of $`M_A`$.

Each residual band has commuting Pauli axes in this canonical algebra.
Their core projections commute, so a fixed core Clifford diagonalizes
them. At a fixed central character, the band is consequently a fixed
core-Clifford conjugate of

```math
\mathrm{diag}_u\omega^{p_j(z,u)},\qquad
\omega=e^{i\pi/4},                                  \qquad(13)
```

where p_j is a weighted cubic phase polynomial modulo eight. Its linear
coefficients are arbitrary, pair coefficients are even, and triple
coefficients are divisible by four. For each fixed core string u, every
binary digit of $`p_j(z,u)\bmod8`$ has Boolean degree at most four.
Indeed the weight-one terms have degree one, weight-two terms degree at
most two, and weight-four terms degree at most three. A binary carry into
the four's digit has total weight at most four, hence degree at most four.
This also follows from the elementary-symmetric formulas for binary
digits of a sum of bits.

The central blocks of the Clifford J have the form

```math
J(z)=\gamma\,i^{a(z)}X^{u(z)}Z^{v(z)}J_* ,           \qquad(14)
```

where J_* is a fixed k-qubit Clifford, u,v are affine bit vectors, and
the scalar is a Clifford phase on the central bits. The Pauli convention
in (14) is ordered X then Z. One proof removes the affine core Pauli
labels using controlled X and controlled Z Cliffords from the central
wires; the remaining operation is a central diagonal Clifford. Its
scalar exponent is linear plus twice quadratic modulo four. Switching
to Hermitian Pauli representatives only adds the explicit normalization
factor $`i^{u(z)\cdot v(z)}`$. Even with that convention, every scalar
digit has degree at most four. The common gamma is fixed, not an
input-dependent scalar.

The entire finite core matrix, and any specified diagonal entry after
the fixed core basis change, is determined by the following Boolean
features:

| Features | Number at k at most five | Degree bound |
| --- | ---: | ---: |
| Three digits of each band phase at every core string | $`9\,2^k\le288`$ | 4 |
| Affine Pauli labels in J | $`2k\le10`$ | 1 |
| Scalar phase digits in J, with a spare digit allowed | 3 | 4 |
| Central offsets of the final commuting core axes | $`k\le5`$ | 1 |
| Core eigenvalue bits selecting the entry | $`k\le5`$ | 1 |

There are fewer than 320 features, each of degree at most four. V_1 has
eighth-root entries up to one common scalar: the original target is a
Boolean sign, every phase peeled in Sections 2 and 3 is an eighth root,
and the extracted diagonal Clifford C also has eighth-root entries.
Equal feature vectors determine the same exact core matrix entry.
Thus each digit of its phase exponent is a Boolean function of these
features on the attainable feature vectors. Extend that function
arbitrarily to the other vectors. Its algebraic normal form has degree
at most 320 in the features, so composition gives degree at most 1280
in the original input bits. This is a finite exact algebra argument;
there is no continuity estimate or count of physical dirty wires.

## 6. Restore every peeled phase and conclude

All peeled pair phases have ordinary exponent degree at most five modulo
eight, and all peeled silent-band phases have weighted cubic exponents.
Their aggregate exponent therefore has ordinary degree at most five.
For any integer multilinear polynomial of degree d modulo eight, each
binary digit has Boolean degree at most $`4d`$: expand its coefficients
into binary weights and use that the four's digit of a sum has carry
degree at most four in the individual monomial bits. Hence the three
digits of the aggregate peeled exponent have degree at most twenty.
The three digits of the extracted diagonal Clifford C have degree at
most four, a conservative bound.

The final Boolean sign bit is a Boolean function of the core features,
these three Clifford digits, and the three aggregate peeled digits.
Its degree is therefore at most

```math
1280+3\cdot4+3\cdot20=1352\lt2048.                  \qquad(15)
```

This restores all literal phases before taking the Boolean degree and
proves (1). The projective version removes only one common scalar.

The retained
[three-layer approximation-to-exactness result](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
also applies:
below $`2^{-530690}`$ occupied operator error, allowing a free common
phase, the Boolean occupied action is exact. Thus the same degree
obstruction has that uniform error floor within the Pauli-boundary
class. The statement includes arbitrary dirty/reference entanglement.

The unrestricted one-clean three-layer degree question remains open.
A counterfamily of unbounded degree would have to use a non-Pauli clean
image at some native-band boundary. Conversely, no change of representation
that merely removes the one-relation axis ambiguity is by itself needed
for this restricted theorem: the exceptional relation-four span was
retained and paid for in the five-qubit core bound.

The [exact regression fixtures](../../tests/test_pauli_boundary_three_band.py)
exercise the hidden-relation singleton with its unused-code sign and
actual inverse, the ordered-Pauli central Clifford phase, and the quartic
binary carry of a native phase. These check the delicate identities;
the width-independent theorem is established by the argument above.
