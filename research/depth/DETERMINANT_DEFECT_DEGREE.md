# A computational degree bound from the common determinant defect

[Non-Pauli code decoding](NONPAULI_CODE_DECODING.md#5-the-zero-pair-zero-character-logical-defect-is-abelian) ·
[Partial-diagonal three-layer features](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md) ·
[Clean-code Clifford rigidity](CLEAN_CODE_CLIFFORD_RIGIDITY.md)

The common determinant cocycle does give a bridge from intrinsic code
coordinates to computational input coordinates. Its constant-phase kernel
identifies a full computational-Z subgroup that can be corrected with
one actual Clifford fixing the occupied clean sector. A bounded
codimension of that kernel consequently bounds Boolean degree. This
note does not prove that its codimension is bounded for every native
three-band representation.

## 1. The actual circuit and its common cocycle

Let V be an actual circuit on Q physical wires, with one clean qubit c,
at most three native T bands, and arbitrary dirty wires. Suppose its
occupied action is a literal or projective Boolean sign, with complete
clean and dirty return. Write

```math
V=K_C^\dagger E_C^\dagger E_BE_AK_A,
\qquad U_C=E_CK_C,\qquad U_{AB}=E_BE_AK_A.
                                                               \qquad(1)
```

Let z range over the full $`Q-1`$ dimensional computational-Z label
space of the non-clean wires, including dirty wires. Define

```math
C_z=U_CZ(z)U_C^\dagger,\qquad
A_z=U_{AB}Z(z)U_{AB}^\dagger,\qquad W_z=C_zA_z.
```

The C_z are actual Hermitian Cliffords. The two families C_z,A_z are
exact commuting representations and agree on the shared clean code
$`\Pi=(I+S)/2`$, where $`S=U_CZ_cU_C^\dagger`$.

Assume the hypotheses of
[the common determinant lemma](NONPAULI_CODE_DECODING.md#5-the-zero-pair-zero-character-logical-defect-is-abelian):
S is non-Pauli with intrinsic active rank different from two, the clean
matching has no middle pairs, and every logical probe character nu_z is
zero. These hypotheses provide one common Pauli algebra K of symplectic
rank two. In fixed intrinsic central coordinates u and its one quantum
core, they give

```math
\begin{gathered}
W_z(u)=\Pi_+(u)+\alpha_z(u)\Pi_-(u),\\
\alpha_{z+v}(u)=\alpha_z(\sigma_vu)\alpha_v(u),\\
\alpha_z(u)=\kappa_z(u)
              \frac{\delta(\sigma_zu)}{\delta(u)}.
\end{gathered}                                                    \qquad(2)
```

Here each sigma_v is the actual affine central permutation induced by
C_v. The same delta is the one-core determinant of the actual middle
band. It is a central diagonal Clifford phase up to one retained fixed
scalar. Each kappa_z is a fixed fourth root times an affine Boolean
sign. Thus every alpha_z is a central diagonal Clifford phase. The u
coordinates are not assumed to be computational inputs of V.

## 2. The constant-defect subgroup has a paid computational correction

Define

```math
H=\{z:\alpha_z(u)\text{ is independent of }u\},
\qquad r=(Q-1)-\dim H.                            \qquad(3)
```

The cocycle law makes H a binary linear subspace. In particular, the
product of two constant alpha-functions is still constant after the
central permutation in (2). For z in H, write its value as a_z. Since
$`\alpha_0=1`$ and $`z+z=0`$, the cocycle gives $`a_z^2=1`$.
Consequently there is a linear character
$`\chi:H\longrightarrow\mathbb F_2`$ such that

```math
\alpha_z=(-1)^{\chi(z)},\qquad W_z=S^{\chi(z)}
\quad(z\in H).                                   \qquad(4)
```

These are full-space identities, not occupied restrictions. Pull them
back through the actual endpoint encoder U_C. Since
$`U_{AB}=U_CV`$, equation (4) gives

```math
VZ(z)V^\dagger=Z(z)Z_c^{\chi(z)}
\quad(z\in H).                                   \qquad(5)
```

Extend chi to a linear functional on all non-clean computational-Z
labels. Choose an X-only Pauli Q_X realizing it by commutation, and put

```math
D=|0\rangle\langle0|_c\otimes I
  +|1\rangle\langle1|_c\otimes Q_X,
\qquad V'=D^\dagger V.                           \qquad(6)
```

D is the product of CNOTs from c to the support of Q_X. It fixes every
clean-zero vector literally, commutes Z_c, and obeys
$`DZ(z)D^\dagger=Z(z)Z_c^{\chi(z)}`$ on H. Thus

```math
V'PV'^\dagger=P\quad
\text{for every }P\in M=\langle Z_c,Z(H)\rangle,
\qquad\dim M=Q-r.                                \qquad(7)
```

Exact clean return implies $`[V,Z_c]=0`$ on the whole space by
unitarity, so the Z_c assertion in (7) includes the unused sector.
The occupied action of V' is exactly that of V, with the same literal
or projective scalar.

The correction has at most $`Q-1`$ CNOT gates, no T gates, and no
added wire. Its actual inverse is the same commuting CNOT word. In
particular,

```math
\begin{gathered}
D_T(V')\le D_T(V)\le3,\qquad T(V')=T(V),\\
G(V')\le G(V)+Q-1,\\
G((V')^\dagger)\le G(V^\dagger)+Q-1.
\end{gathered}                                                    \qquad(8)
```

The inverse of the original native word is charged separately. Every
identity above holds on all physical columns and external references;
the correction is not a chosen extension of an occupied Clifford.

## 3. Bounded defect codimension bounds Boolean degree

Choose fixed invertible linear coordinates on the non-clean wires in
which H is generated by the first $`Q-1-r`$ computational Z operators.
Keep c as its own central coordinate. Equation (7) is exactly the
full-space preservation premise of the
[partial-diagonal finite-feature theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants)
with r quantum coordinates. Its explicit constant is

```math
\begin{gathered}
q_r=2r+4,\qquad N_r=2^{2r^2+5r+10},\\
B_r=9\,2^r+2r+3,\qquad
D_r=48q_rN_r+4B_r+r.
\end{gathered}                                                    \qquad(9)
```

Apply the theorem to the actual three-band V'. Then set the central
clean bit to zero in its finite-feature functions. Restricting that
bit does not increase their Boolean degrees. On this sector the
remaining quantum-core matrix is the prescribed Boolean diagonal sign,
so the theorem's matrix-entry phase extraction proves

```math
\deg_{\mathbb F_2}g\le D_r.                       \qquad(10)
```

The linear change of the non-clean computational coordinates does not
increase algebraic degree. The sign is independent of dirty input bits,
but the preserved subgroup and the full operator proof include those
bits from the start.

When r=0, V' is fully diagonal on all Q physical wires. The sharper
[zero-clean three-layer theorem](THREE_LAYER_NATIVE_STRUCTURE.md)
gives its modulo-eight phase polynomial ordinary degree at most five.
Setting c=0 preserves that degree, and extracting the occupied Boolean
sign gives $`\deg g\le5`$. The unused clean block need not itself
be a Boolean sign. In particular, all interior clean images
may be non-Pauli in this argument; neither an interior Pauli hypothesis
nor a shared Pauli clean flip was assumed.

For a literal Boolean XOR query, conjugate its output by Hadamards
before forming (1)--(3). The resulting sign is
$`(-1)^{z f(x)}`$ and has degree $`1+\deg f`$ for nonzero f. Thus
the same hypotheses and defect codimension for this sign representation
give $`\deg f\le D_r-1`$, or $`\deg f\le4`$ when r=0. No claim is
made that the logical defect subgroup is unchanged by this output
basis change; it is the sign representation to which the theorem is
applied.

The retained [R74 approximation-to-exactness consequence](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
first forces exact occupied Boolean action below operator error
$`2^{-530690}`$, including a freely optimized common scalar. For a
supplied native circuit satisfying the stated code and cocycle
hypotheses, (10) therefore has the same uniform error floor. The
complete dirty/reference contract is unchanged.

## 4. One common bounded central quotient is sufficient

There is an intrinsic-coordinate sufficient condition for bounded r.
Two labels z,v belong to the same H coset if and only if alpha_z and
alpha_v differ by a constant scalar. Indeed the cocycle and v+v=0 give

```math
\alpha_v(\sigma_vu)\alpha_v(u)=1.
```

If $`\alpha_z=c\alpha_v`$, equation (2) gives
$`\alpha_{z+v}=c`$. Conversely, if $`\alpha_{z+v}=c`$, use
$`\sigma_v^2=I`$ in the same equality to obtain
$`\alpha_z=c\alpha_v`$. Hence

```math
2^r=\#\{\alpha_z\text{ modulo constant scalars}:z\}.
                                                               \qquad(11)
```

Suppose the entire alpha-family factors through one fixed
t-dimensional linear quotient of the intrinsic central coordinates.
Choose a linear section of that quotient. Each induced phase is still
a diagonal Clifford phase. Modulo a constant, such a phase on t bits
has t linear coefficients modulo four and $`t(t-1)/2`$ quadratic
coefficients modulo two. Therefore

```math
2^r\le4^t2^{t(t-1)/2},\qquad
r\le\frac{t(t+3)}2.                              \qquad(12)
```

A fixed common t thus yields the width-independent bound
$`\deg g\le D_{t(t+3)/2}`$. This conclusion uses the complete
cocycle family to produce the computational subgroup (7). It does
not require identifying those intrinsic central bits individually with
target input bits. A bounded number of controls for each separate
probe, with different controls for different probes, is not the common
quotient hypothesis in (12).

## 5. The first remaining inference

The native determinant lemma proves that each alpha_z is a low-degree
phase in intrinsic coordinates. It does not bound the number of
distinct such phases across the full logical family. Inferring a
uniform bound on r from that pointwise degree statement would be the
first unsupported step. The common quotient in Section 4 is an
additional hypothesis, not a consequence of the one-core algebra.

The [actual identity-target family](CLEAN_TRANSITION_PARAMETERS.md#5-unbounded-determinant-codimension-for-literal-occupied-identity)
answers the first question negatively: all hypotheses of Section 1
hold in a three-native-band split, but H is zero and r grows with the
number of non-clean wires. The same family needs that many common
linear central coordinates. Its first interior clean image is Pauli
and its second has active rank one; it does not settle codimension
growth under a stricter both-non-Pauli hypothesis.

The remaining quantitative question is whether the actual bands and
literal Boolean occupied target provide a different paid representation
with bounded r, or a degree argument that avoids this parameter.
All native bands, clean-silent factors, and their joint dependence
remain part of that question. Thus r is a parameter of the supplied
representation, not a measure of target hardness. This note closes the
fixed-r subcase and does not assert a general one-clean three-layer
degree theorem or a complete-frame resource improvement.
