# Native transport through a diagonal third-level source and a bounded core

[Six-layer comparison](READ_ONLY_SIX_LAYER_REDUCTION.md) ·
[Diagonal-source transport](NATIVE_CLIFFORD_TRANSPORT.md) ·
[Bounded Pauli deficiency](BOUNDED_DEFICIENCY_TRANSPORT.md) ·
[Native source structure](NATIVE_BAND_SOURCE_STRUCTURE.md)

A diagonal third-level source coupled to a fixed-size quantum core has
source-uniform native transport bounds even when its Pauli deficiency
is unbounded. The proof retains the diagonal phase expression while
normalizing the exact signed Pauli domain. It then reconstructs each
output label coordinate from at most four fixed control-sector blocks.
The literal sign follows with the modulo-four carry retained.

The actual hybrid factorization and third-level membership are explicit
source hypotheses. Neither is inferred from full workspace return.
All operators, adjoints, signs and phases below are literal; fixed
coordinate changes introduce no physical source or helper.

## 1. Statement and native contract

Let the actual shared source F of the native six-layer comparison, up to
fixed input/output Cliffords, have either form

```math
F=D(U_A\otimes I_B)\quad\text{or}\quad F=(U_A\otimes I_B)D,
```

where $`|A|=a`$, $`U_A`$ is an arbitrary fixed unitary, and D is a
weighted-cubic modulo-eight diagonal phase on A+B:

```math
D|z\rangle=\omega^{p(z)}|z\rangle,\qquad \omega=e^{i\pi/4},\qquad
p(z)=c+\sum_i a_i z_i+2\sum_{i\lt j}b_{ij}z_i z_j
       +4\sum_{i\lt j\lt k}c_{ijk}z_i z_j z_k\pmod8.
```

A fixed scalar on F can be retained in either end factor, understood
as a Clifford up to that actual scalar. Assume F
itself belongs to the third level $`\mathcal C_3`$. Let the exact native
inner Clifford family L and actual outer comparison K of the
[six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md#2-a-clifford-valued-native-pullback)
satisfy

```math
\|K_x-FL_xF^\dagger\|_{\rm op}\le2\eta.
```

Both families are Clifford up to their retained actual scalar phases.
The [native input ledger](NATIVE_CLIFFORD_TRANSPORT.md#4-native-degrees-including-modulo-eight-carries)
gives label degree two and literal fixed-probe sign degree six for L
and its actual inverse, preserved under fixed Clifford coordinates.
Its moving-probe ledger gives label degree $`d+2`$ and sign degree
$`3d+6`$ for an input label of Boolean degree d.

**Hybrid transport theorem.** For $`a\ge1`$, define

```math
N(a)=3\,2^{5a}+15a+2,\qquad
B(a)=24N(a)+8a^2+12a.
```

Then the following source-uniform bounds hold:

```math
\eta\lt \sqrt2/8\quad\Longrightarrow\quad
\deg S_K\le4B(a)+2,\qquad
\deg\kappa_K(p)\le12B(a)+6
```

for every fixed Hermitian probe p. Consequently the actual zero-clean
read-only six-layer query satisfies $`\deg f\le12B(a)+6`$, improving to
$`4B(a)+2`$ when the outer target shear is visible. These constants are
deliberately loose; neither sharpness nor attainment is claimed. The
radius is independent of a, work width, the coefficients of D, and
$`U_A`$. The degree bounds depend only on a.

For $`a=0`$, the retained diagonal theorem gives the sharper label/sign
bounds 26/78 at the same radius. For $`a=1`$, the one-core collapse in
section 8 gives those same sharper bounds. The general theorem does not
claim that F is diagonal C3 in fixed Clifford coordinates.

Use the tensor-factorizing Hermitian section
$`\tau(x,z)=i^{x\cdot z}X^xZ^z`$, where the exponent's dot product is
over the integers. Thus
$`\tau(p_A,p_B)=\tau(p_A)\otimes\tau(p_B)`$ literally. All Boolean
degrees below are algebraic degrees over $`\mathbb F_2`$.

## 2. Strip many exact spectator directions using only B coordinates

Define the exact Pauli-preserving domain, including the zero label, by

```math
E=E_F=\{p:F\tau(p)F^\dagger\text{ is a signed Hermitian Pauli}\}.
```

It is a binary linear space, since conjugation preserves Pauli products
and their literal scalar phases. F commutes literally with every
$`Z_B`$ Pauli. Let

```math
R_B=\{b:X_B^b\in E\}.
```

For $`b\in R_B`$, the output Pauli has B-X label b: this follows from its
commutation pairing with every fixed $`Z_B`$. Its active A label g(b)
is linear into $`\mathbb F_2^{2a}`$. Put $`R_0=\ker g`$, so
$`\dim(R_B/R_0)\le2a`$.

For b in $`R_0`$ the output is entirely on B. D has exactly the same signed
action on these $`X_B`$ translations in either ordering. For F=D $`U_A`$ this
is immediate because $`U_A`$ commutes with $`X_B`$. For F=$`U_A`$ D, conjugate the
B-only output back by $`U_A`$, which also fixes it literally.

Thus

```math
D X_B^b D^\dagger=(-1)^{s(b)}\tau(b,c(b))\qquad(b\in R_0).
```

The linear label map $`c:R_0\to B^*`$ is symmetric on $`R_0`$ by
commutation. Extend it to a symmetric map on B, retaining its prescribed
pairings with all of B. Choose a diagonal B-only Clifford $`G_B`$
realizing this map, and correct its basis signs by a Z-type B Pauli.
This matches D's literal action on every $`X_B^b`$ for $`b\in R_0`$.

The phase $`D'=G_B^\dagger D`$ commutes literally with those translations,
and remains weighted cubic. Because $`G_B`$ is B-only, it commutes with
$`U_A`$, so removing $`G_B`$ on the left leaves D' $`U_A`$ or $`U_A`$ D' in the same
ordering. A B-only invertible binary coordinate change makes $`R_0`$ the
translation space of some B coordinates. The phase of D' is independent
of those coordinates, so they split off as an exact identity subsystem.
No B initialization is used. Restrict the phase polynomial to zero on
the identity coordinates as an analytical description of that factor.

After this stripping, retain the names A,B,F,E. The remaining $`R_B`$ has
dimension at most $`2a`$. The kernel of the projection from E to the
A Pauli labels consists of B Paulis and has dimension
$`|B|+\dim R_B`$. Since all $`Z_B`$ lie in E,

```math
\dim(E/Z_B)\le4a.
```

Consequently the B-X projection of E has dimension at most $`4a`$. Another
B-only invertible binary change puts its support on at most $`4a`$ B wires.
Include them in a new core C, with

```math
m=|C|\le5a.
```

Let R denote the remaining B wires. E has no $`X_R`$ label and contains all
$`Z_R`$. Hence, as a binary subspace,

```math
E=E_C\oplus Z_R.
```

Indeed any $`Z_R`$ tail of an element of E may be removed by another element
of E. Up to the fixed Cliffords already extracted, F still has D U or
U D form, D is weighted cubic, and U is a fixed operator supported on C.

## 3. Match the remaining signed action with controlled Paulis

F fixes $`Z_R`$. For p in $`E_C`$ write its actual image as

```math
F\tau(p)F^\dagger=(-1)^{s(p)}\tau(g(p))Z_R^{f(p)}.
```

The core map g is symplectic and injective. If its image label were zero
for nonzero p, the image would be a signed $`Z_R`$ Pauli, and conjugating
back by F (which fixes $`Z_R`$) would equate a nonidentity core Pauli to a
$`Z_R`$ Pauli, impossible. The Pauli multiplication phases show that the
signed core action extends to a fixed Clifford $`C_0`$ on C, including when
$`E_C`$ is degenerate. After $`C_0^\dagger`$, the image is literally

```math
\tau(p)Z_R^{f(p)}.
```

For each R coordinate j, choose a core Pauli label $`v_j`$ representing
the linear functional $`p\mapsto f(p)_j`$ under the symplectic pairing
on C. Define
the actual controlled-Pauli Clifford

```math
W=\prod_j \bigl(|0\rangle\langle0|_j\otimes I_C
                 +|1\rangle\langle1|_j\otimes\tau(v_j)\bigr)
```

in one fixed order. Its conjugation on $`\tau(p)`$ gives precisely the
same $`Z_R`$ tail. Therefore $`H=W^\dagger C_0^\dagger F`$ commutes
literally with all Paulis of $`E_C`$ and all $`Z_R`$. Left multiplication
by these Cliffords leaves the exact Pauli-preserving domain unchanged.
The ordered product is essential: replacing it by the Hermitian section
of the sum of its labels can change the scalar to a non-Clifford phase.

A fixed Clifford on C now puts $`E_C`$ into canonical form: a full Pauli
space on $`S_C`$, a Z space on $`R_C`$, and no domain directions on a remaining
quantum core Q. Conjugate H by this fixed core Clifford. Include $`R_C`$ in
the control register and $`S_C`$ in the identity spectator register. Together
with the identity wires stripped in section 2, the final source has

```math
E_H=\mathcal P_S\oplus Z_R,\qquad
H=I_S\otimes\sum_z |z\rangle\langle z|_R\otimes H_z.
```

All $`H_z`$ are q-qubit unitaries, where $`q=|Q|\le a`$. To verify this
last bound, write $`E_H`$ as full Pauli space on s wires plus Z space on
r wires, with q remaining wires. Its largest isotropic dimension is
$`s+r=w-q`$. The original E contained the isotropic space $`Z_B`$ of
dimension $`w-a`$. Fixed symplectic coordinates preserve this dimension,
so $`q\le a`$. Removing an
identity tensor factor changes w and s equally.

## 4. The controlled source retains a bounded phase expression

This is the key preservation assertion, stated with an explicit feature
budget. Split a final control string z into the original wide controls
$`z_R`$ and at most m core-control bits $`z_C`$ introduced in section 3.

Before the last fixed core Clifford and block restriction, the source
matrix on C is obtained from

```math
W_{z_R}^\dagger C_0^\dagger
       [D_{z_R}U\ \text{or}\ U D_{z_R}]
```

by fixed core matrix multiplications. Here U is a fixed m-qubit matrix,
$`D_{z_R}`$ is diagonal with entries $`\omega^{p_c(z_R)}`$, for the
$`2^m`$ fixed core strings c, and each $`p_c`$ is weighted cubic modulo
eight. B-only linear substitutions, diagonal Clifford removal, and
restriction of an identity variable preserve that phase form.

Write $`v_j=(\alpha_j,\beta_j)`$ in Hermitian Pauli coordinates. The
actual ordered product $`W_{z_R}`$ has the computational-basis form

```math
W_{z_R}=i^{\psi(z_R)}X^{u(z_R)}Z^{v(z_R)},
```

where u,v are binary linear maps and

```math
\psi(z_R)=\sum_j (\alpha_j\cdot\beta_j)z_j
  +2\sum_{j\lt k}(\beta_j\cdot\alpha_k)z_jz_k\pmod4.
```

Integer dot products and the chosen order are retained. Taking its
adjoint requires no new variable, since the full matrix is already
determined by these features.

The final q-qubit $`H_z`$ is a block of fixed core Clifford conjugations
of this matrix, at the specified core-control bits $`z_C`$ and any fixed
identity-spectator basis state.

More explicitly, writing w for the old-wide controls, r for the new core
controls, and V for the final fixed core basis Clifford, its entries are

```math
H_{(r,w)}[y,y']
=\langle0_S,r,y|V^\dagger W_w^\dagger C_0^\dagger
   [D_wU\ \text{or}\ UD_w]V|0_S,r,y'\rangle.
```

The displayed $`0_S`$ selects a matrix coefficient of an already proved
identity tensor factor; it is not a physical initialization. The entire
$`H_z`$ is consequently a deterministic function of the following bits:

1. The three bits of each $`p_c(z_R)`$, for at most $`2^m`$ strings c.
2. The $`2m`$ bits $`u(z_R),v(z_R)`$.
3. The two bits of $`\psi(z_R)`$.
4. At most m core-control bits $`z_C`$ selecting the block.

There are at most

```math
3\,2^m+3m+2\le N(a)
```

such bits. If the coordinates of $`z=z(x)`$ have Boolean degree at most h,
every feature bit has degree at most $`4h`$. The
[weighted-cubic carry lemma](NATIVE_CLIFFORD_TRANSPORT.md#4-native-degrees-including-modulo-eight-carries)
gives $`h,2h,4h`$ for the three phase bits; the modulo-four phase gives
at most $`2h`$, and the linear labels and selectors give h. Arbitrary fixed U
entries are constants in the finite feature lookup. There is no source-
dependent separation or approximation in this representation.

## 5. Exact endpoint and domain compatibility

At $`\eta\lt \sqrt2/8`$, the
[third-level rounding lemma](NATIVE_CLIFFORD_TRANSPORT.md#1-a-robust-exactness-lemma-for-every-third-level-source)
makes $`F L_x F^\dagger`$ projectively equal to $`K_x`$: the actual
endpoint error is at most $`2\eta\lt \sqrt2/4`$. The lemma does not require
$`F^\dagger\in\mathcal C_3`$. Make all fixed coordinate changes above
and keep the names H,L,K. Their input degree ledgers remain unchanged.

L normalizes $`E_H`$: for $`p\in E_H`$, the identity

```math
K(H\tau(p)H^\dagger)K^\dagger
 =H(L\tau(p)L^\dagger)H^\dagger
```

makes the right-hand side Pauli. Hence $`S_L`$ maps $`E_H`$ into itself,
and invertibility makes it onto. It also preserves
$`\mathrm{rad}(E_H)=Z_R`$ and $`E_H^\perp`$.
Since H commutes $`E_H`$ literally, K and L have identical signed actions
there. Put

```math
J_x=L_x^\dagger H L_x H^\dagger.
```

J is an exact Clifford up to its actual scalar and commutes $`E_H`$. The
actual K differs from L J only by a scalar, irrelevant only now to its
Pauli conjugation action.

The Clifford L permutes computational R blocks by an affine map
$`z\mapsto\phi_x(z)`$. Each fiber is a Clifford on Q+S. Preservation
of $`E_H`$
forces that fiber to factor into Q and S Cliffords up to scalar. Thus

```math
J_{x,z}=L_{Q,x,z}^\dagger H_{\phi_x(z)}L_{Q,x,z}H_z^\dagger.
```

The S factors and all fiber scalar phases cancel in this literal formula.

For fixed z, $`\phi_x(z)`$ has degree at most six: read it from
$`L_x^\dagger Z_{R_j}L_x`$, whose labels have degree two and signs
degree six. Similarly the projective tableau of
$`M_{x,z}=L_{Q,x,z}^\dagger`$ is read from
$`L_x^\dagger\tau(p_Q)L_x`$. Such images lie in Q Paulis times $`Z_R`$
because $`E_H^\perp`$ is invariant. Restricting to the fixed input z
adds a fixed linear combination of degree-two $`Z_R`$ labels to its
degree-six sign. The local inverse tableau therefore has $`4q^2`$ label
bits of degree two and $`2q`$ sign bits of degree six. This uses the
inverse fiber at the fixed input sector z; no substitution by
$`\phi_x^{-1}`$ is needed.

By section 4, $`H_{\phi_x(z)}`$ has at most $`N(a)`$ feature bits, each
of degree at most 24. $`H_z`$ is fixed. The local scalar of M cancels,
so $`J_{x,z}`$ is a
well-defined function of those source features and its local projective
tableau. The sum of their degree bounds is at most

```math
24N(a)+8q^2+12q\le B(a).
```

Every Boolean function of this finite feature tuple therefore has degree
at most $`B(a)`$, by multilinear interpolation and substitution. In
particular, any Boolean quantity uniquely determined by a fixed number
c of blocks $`J_{x,z_i}`$, at fixed $`z_i`$ independent of x, has degree
at most $`cB(a)`$. A function specified only on the feature tuples that
occur can be extended arbitrarily to all Boolean tuples before
interpolation. No equality test over wide spectator labels occurs.

## 6. Recover a fixed probe, including its literal sign

J acts identically on S, so omit that factor. Fix a Hermitian input probe
$`\tau(p_Q)\otimes\tau(u,v)_R`$, with u,v fixed binary strings. Because
J is Clifford and commutes $`Z_R`$, its output has form

```math
(-1)^{\kappa(x)}\tau(t(x))\otimes\tau(u,v'(x)).
```

For every fixed control string z, define

```math
A_z=J_{x,z+u}\tau(p_Q)J_{x,z}^\dagger.
```

Reading the computational action of the tensor-factorizing Hermitian
section gives the exact identity

```math
A_z=\alpha(x)(-1)^{(v'(x)-v)\cdot z}\tau(t(x)),
\qquad
\alpha(x)=(-1)^{\kappa(x)}i^{u\cdot(v'(x)-v)}.
```

The last dot product is over the integers modulo four. At $`z=0`$,
$`A_0`$ determines t and $`\alpha\in\{1,i,-1,-i\}`$ from the two
blocks at 0 and u. Their bits consequently have degree at most
$`2B(a)`$. For each R coordinate j,
$`A_{e_j}A_0^\dagger=(-1)^{v'_j-v_j}I`$. This uses only the four
blocks at $`0,u,e_j,e_j+u`$. Thus every $`v'_j`$ has degree at most
$`4B(a)`$.

Write $`\alpha=i^\mu`$, $`\mu\in\mathbb Z_4`$. The scalar relation says

```math
2\kappa=\mu-u\cdot(v'-v)\pmod4.
```

The low bit of the integer sum $`u\cdot v'`$ has degree at most
$`4B(a)`$, and its high bit is a second elementary symmetric polynomial
of those bits, of degree at most $`8B(a)`$. Subtracting the fixed
$`u\cdot v`$ and the variable two-bit $`\mu`$ preserves the bound
$`8B(a)`$: the extra borrow multiplies bits of degrees at most
$`4B(a)`$ and $`2B(a)`$. Divisibility by two is guaranteed by
the exact Pauli identity. Hence

```math
\deg S_J\le4B(a),\qquad
\deg\kappa_J(p)\le8B(a).
```

This sign recovery retains the Hermitian-section carry. It does not
sum an unbounded list of output-generator signs or treat $`\alpha`$ as real.

## 7. Return to the native comparison and full query

Conjugate the fixed probe first by J and then by L. The incoming label
to L has degree $`d=4B(a)`$, and the incoming sign has degree at most
$`2d`$. The native moving-probe ledger gives output label degree $`d+2`$
and sign degree at most $`3d+6`$. Restoring the fixed output Clifford
adds a cubic sign on a degree-$`(d+2)`$ label, also bounded by $`3d+6`$.
This proves the stated label bound $`4B(a)+2`$ and fixed-probe sign bound
$`12B(a)+6`$.

The existing visible outer-shear extractor adds only a quadratic baseline;
the invisible-shear fixed-probe detector adds only a degree-six baseline.
They give the target bounds stated in section 1.

All norm estimates begin at the actual complete endpoint comparison and
the original full query norm. Exact rounding is used only with the
explicit C3 source assumption. Every logical column, literal target
phase, arbitrary dirty and external reference input, actual inverse,
leakage bound and full workspace return remains required. All original
T and elementary Clifford gates, source occurrences, controls, inverses
and simultaneous work remain charged. The coordinate changes and block
evaluations are analysis; there is no initialized sector, free source,
new control, helper or physical source synthesis.

This theorem does not establish that an arbitrary three-band prefix has
the hybrid form, nor that its core size is bounded by an overlap parameter.
It does not assert the source is semi-Clifford. Those remain separate
source-classification questions. Its error radius is uniform in work
width and in $`U_A`$, while its degree bound depends on the specified a.

## 8. A sharper one-core consequence

If the core is one qubit, third-level membership forces the whole hybrid
source to be diagonal C3 in fixed Clifford coordinates. This certifies
the retained 26/78 transport bounds without the general feature budget.

### Restricting the fixed core

Let D be a weighted-cubic modulo-eight diagonal source on A+B, and let
$`U_A`$ act only on A. Suppose either

```math
F=D(U_A\otimes I_B)\in\mathcal C_3
\quad\text{or}\quad
F=(U_A\otimes I_B)D\in\mathcal C_3.
```

For each fixed computational string z of B, denote the diagonal A block
of D by $`D_z`$. Then respectively

```math
D_zU_A\in\mathcal C_3(A)
\quad\text{or}\quad
U_AD_z\in\mathcal C_3(A).
```

Proof. F commutes literally with every Z Pauli on B. For any Hermitian
Pauli P on A, $`J=F(P\otimes I)F^\dagger`$ is a Clifford and also
commutes with all $`Z_B`$. Thus J is block diagonal in B. Each block
$`J_z`$ is a Clifford on A, up to its actual scalar: for every A Pauli Q,
the full Clifford image $`J(Q\otimes I)J^\dagger`$ is a Pauli commuting
with all $`Z_B`$, so it is a signed A Pauli tensored with a $`Z_B`$
Pauli. Restriction to z is therefore a signed A Pauli. But
$`J_z=F_z P F_z^\dagger`$. This proves the C3 claim for $`F_z`$.

This is an algebraic restriction of an identity valid on every B column;
it does not initialize B or discard the other columns. Every source
phase is retained in $`D_z`$ and $`F_z`$.

In particular, with $`z=0`$, put $`V=D_0U_A`$ in the first ordering and
$`V=U_AD_0`$ in the second. Since D and the embedded $`D_0`$ are diagonal, one
gets literal factorizations

```math
F=D'(V\otimes I_B)
\quad\text{or}\quad
F=(V\otimes I_B)D',
\qquad D'=D(D_0^\dagger\otimes I_B),
```

where $`D'`$ is still diagonal C3 and V is a C3 unitary on A. The scalar
contained in $`D_0`$ is canceled literally in $`D'`$.

### One-qubit C3 classification

Every one-qubit C3 unitary V has

```math
V=C_1T^\epsilon C_2,\qquad \epsilon\in\{0,1\},
```

up to a scalar that can be assigned to one fixed Clifford factor.

To see this, the images of X,Y,Z are traceless Hermitian Cliffords whose
Bloch vectors form an orthonormal frame. Flat Clifford Pauli support
implies that a traceless one-qubit Hermitian Clifford has either one or
two nonzero Pauli coefficients: there are only three nonidentity Paulis,
and the support size is a power of two. Its Bloch direction is therefore
either a signed coordinate vector or a signed edge vector
$`(\pm e_i\pm e_j)/\sqrt2`$.

If one direction is an edge vector, after coordinate signs and a
permutation take it to $`(e_1+e_2)/\sqrt2`$. Among the allowed directions,
the only orthogonal choices are $`\pm e_3`$ and
$`\pm(e_1-e_2)/\sqrt2`$. Thus the frame has exactly two edge vectors
and one coordinate vector, and is a fixed Clifford transform of the
frame made by T. If no direction is an edge, V is Clifford. Orientations
are handled by Clifford signs and permutations; $`T^\dagger`$ differs from T
by a Clifford $`S^\dagger`$. This proves the factorization.

### One-core collapse

**Theorem.** Suppose A consists of one qubit, D is an arbitrarily wide
diagonal C3 source, $`U_A`$ is any fixed one-qubit unitary, and either
$`F=D(U_A\otimes I)`$ or $`F=(U_A\otimes I)D`$ is C3. Then F is diagonal C3 in
fixed Clifford coordinates.

For the first ordering, the block restriction gives $`F=D'V`$ and the
classification gives $`V=C_1T^\epsilon C_2`$. The source $`W=D'C_1`$
is diagonal C3 in fixed Clifford coordinates, and
$`WT^\epsilon=FC_2^\dagger`$ is C3. For $`\epsilon=0`$ the conclusion
already holds; for $`\epsilon=1`$, apply the
[native diagonal-C3 right-extension corollary](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension).

For the second ordering, write $`F=C_1T^\epsilon C_2D'`$. The source
$`W=C_2D'`$ is diagonal C3 in fixed Clifford coordinates and
$`T^\epsilon W=C_1^\dagger F`$ is C3. Apply the
[native diagonal-C3 left-extension corollary](NATIVE_BAND_SOURCE_STRUCTURE.md#4-native-left-extension-without-an-inverse-c3-premise).
Neither argument assumes that $`F^\dagger`$ is C3.

Fixed end Cliffords on F do not change the statement. The theorem does
not require $`U_A`$ to be Clifford, C3, exactly synthesizable by a specified
small word, or separated by a fixed gap from Clifford. C3 membership of
the full source provides the needed core restriction.

The retained [diagonal-source transport theorem](NATIVE_CLIFFORD_TRANSPORT.md#5-consequence-for-the-actual-six-layer-query) therefore gives

```math
\eta\lt \sqrt2/8\quad\Longrightarrow\quad
\deg S_K\le26,\qquad\deg\kappa_K(p)\le78,
```

and the respective visible-shear/general target bounds 26/78. This
auxiliary does not assume that the arbitrary one-qubit factor itself is
third-level; its fixed block restriction supplies the needed C3 core.

## 9. Application: commuting extensions with bounded relation dimension

The hybrid hypothesis follows from a diagonal-C3 source followed by a
commuting T list with bounded binary relation dimension. The list may
be dependent; its length and the width of the diagonal phase may be
arbitrarily large.

**Relation-count theorem.** Let U be diagonal C3 in fixed Clifford
coordinates, and let E be a product of actual T/T-dagger phases about
a commuting list of signed Hermitian Pauli axes. Let r be the dimension
of the binary relation space of that list, counting repeated and zero
labels. Suppose the full source $`F=UE`$ is C3. Put

```math
M(r)=4(2^r-1).
```

There are fixed Cliffords $`C_{\rm out},C_{\rm in}`$, a wide diagonal
C3 phase D, and a unitary V supported on at most $`M(r)`$ qubits such that

```math
F=C_{\rm out}\,D(V\otimes I)\,C_{\rm in}
```

literally, including its scalar. In fact V is the product of at most
$`M(r)`$ commuting Pauli T phases in fixed coordinates and is itself
C3. Its active-qubit count can be bounded more sharply by the dimension
of the isotropic span of those axes.

Thus for $`r\ge1`$ the six-layer comparison has the source-uniform
bounds of section 1 with $`a=M(r)`$, at $`\eta\lt \sqrt2/8`$. For
$`r=0`$ there is no residual core, and the sharper 26/78 bounds apply.
In particular, relation dimension two gives at most twelve bad axes and
at most twelve core qubits. This is a conditional source certificate;
no relation bound is inferred merely from a returned query.

### Remove duplicate axes before applying paired half angles

Absorb U's fixed end Cliffords analytically, conjugating every E axis
together, so U is diagonal. For each nonzero binary label choose one
Hermitian representative and collect all rotations about it. Changing
the representative's sign changes the orientation and a fixed scalar,
which must be retained. An even total exponent is a Clifford power;
an odd exponent leaves one T phase and a Clifford power. Zero-label
rotations are scalars. Since the whole list commutes, collect these
Cliffords and the scalar into an actual fixed right factor:

```math
E=E'C_E.
```

The surviving list $`E'`$ has distinct nonzero binary axes and relation
dimension $`r'\le r`$. Removing a duplicate pair deletes two list
entries and lowers their span rank by at most one; removing a zero
label deletes one relation. The source $`UE'=FC_E^\dagger`$ remains
C3. This preprocessing is necessary because repeated axes have
nonorthogonal images and do not satisfy the paired half-angle lemma.
It is an analytical identity; the original native ledger is unchanged.

### Count the non-Pauli images by parity-check columns

Write $`a_1,\ldots,a_n`$ for the surviving axes and let

```math
\mathcal R=\{v\in\mathbb F_2^n:\textstyle\sum_i v_i a_i=0\},
\qquad\dim\mathcal R=r'.
```

Use a basis of $`\mathcal R`$ as parity-check rows, and denote its
i-th column by $`\sigma_i\in\mathbb F_2^{r'}`$. Nondegeneracy of
the ambient symplectic form gives the exact activation space

```math
\{(\langle a_i,p\rangle)_i:p\text{ a Pauli label}\}
=\mathcal R^\perp.
```

A probe can therefore isolate axis i when $`\sigma_i=0`$, and can
isolate exactly i,j when $`\sigma_i=\sigma_j`$. Put
$`H_i=U\tau(a_i)U^\dagger`$, including the actual axis signs.
These are commuting, traceless, pairwise orthogonal Hermitian Clifford
involutions. Cancel the Clifford image of the isolating probe from its
full-source image. For a single axis, the retained
[half-angle lemma](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension)
makes $`H_i`$ a Pauli. For a pair, the
[paired half-angle classification](ONE_OVERLAP_SOURCE_STRUCTURE.md#2-paired-half-angle-rigidity)
says that either both images are Paulis, or both have the same
four-element Pauli support. Opposite T orientations are absorbed into
the Hermitian involutions' signs.

Hence a nonzero column class either has only Pauli images or has at
most four non-Pauli images: mutually orthogonal images sharing one
four-dimensional real operator space cannot be more numerous. A
singleton class also obeys this bound. The zero-column class has no
bad images. The number of bad axes is consequently at most

```math
4(2^{r'}-1)\le M(r).
```

Write their commuting phase products as
$`E'=E_{\rm good}E_{\rm bad}`$. The good list can still be dependent.

### Absorb the good phases with their literal signs

The [signed diagonal Clifford extension](NATIVE_CLIFFORD_TRANSPORT.md#2-the-diagonal-source-and-its-pauli-directions)
gives a diagonal Clifford G matching U on its whole Pauli-preserving
domain. Thus $`H=G^\dagger U`$ fixes every good axis literally. If
$`R_X`$ is U's preserved translation space, H commutes with all
$`X^v`$ for $`v\in R_X`$ and with all Z Paulis. Its weighted-cubic
phase therefore descends to quotient Z observables with labels in
$`R_X^\perp`$. Restriction to a fixed linear section preserves the
weighted-cubic form modulo eight.

Every good axis has X label in $`R_X`$, so it commutes with these
quotient observables. The good axes and quotient observables together
span an isotropic space. A fixed Clifford K carries it into the standard
Z Lagrangian, giving

```math
D_0=K(H E_{\rm good})K^\dagger
```

as a diagonal C3 phase. Dependencies among good axes merely add phase
terms. Fixed Clifford signs become fixed affine offsets in their
arguments; the signs and scalar remain literal.

### Localize the bad span without losing diagonality

Let $`\mathcal B`$ be the isotropic span of the K-transformed bad
axes, of dimension m, and let $`L_Z`$ be the standard Z Lagrangian.
Put $`b=\dim(\mathcal B\cap L_Z)`$ and $`h=m-b`$. There is a
monomial Clifford J preserving $`L_Z`$ such that

```math
J\mathcal B=\mathrm{span}\{X_1,\ldots,X_h,
                                      Z_{h+1},\ldots,Z_m\}.
```

Here is a construction. Isotropy says that $`\mathcal B\cap L_Z`$
annihilates the X projection of $`\mathcal B`$. An invertible linear
coordinate change, implemented by CNOTs and swaps, puts those spaces
on $`Z_{h+1},\ldots,Z_m`$ and $`X_1,\ldots,X_h`$, respectively.
The remaining basis vectors of $`\mathcal B`$ then have the form

```math
b_i=X_i+\sum_j Q_{ji}Z_j,\qquad 1\le i\le h.
```

Commutation gives $`Q_{ij}=Q_{ji}`$ for $`i,j\le h`$. The displayed
columns extend to a symmetric matrix. Its diagonal Clifford shear,
implemented by S and CZ gates, removes every Z tail while fixing all
Z labels. The resulting J is monomial, so $`J D_0J^\dagger`$ stays
diagonal C3: only a fixed linear or affine substitution acts on its
phase. The Hermitian signs added to the bad axes are retained as their
actual rotation orientations.

Set

```math
D=JD_0J^\dagger,\qquad
V\otimes I=JK E_{\rm bad}K^\dagger J^\dagger.
```

The second operator is supported on m qubits, with $`m\le M(r)`$.
In the normalized coordinates the exact factorization is

```math
F=(G K^\dagger J^\dagger)\,D(V\otimes I)\,(JKC_E).
```

Restoring U's original end Cliffords proves the theorem. Every identity
holds on the full work space, including arbitrary dirty inputs and
external references. The statement bounds the surviving non-Pauli
images, not the raw number of exceptional axes in an original native
description; duplicate cancellations can leave that raw count
arbitrarily large. Sources, actual inverses, helpers and all physical
gates remain charged exactly as in the original circuit.

## 10. Finite verification

The [focused fixtures](../../tests/test_hybrid_diagonal_core_transport.py)
check signed localization in both source orders, the ordered controlled-
Pauli phase formula with noncommuting core axes, and four-block recovery
on every Hermitian probe of an entangling Clifford. The latter exercises
imaginary scalar phases and the integer-dot modulo-four carry. A binary
symplectic fixture checks simultaneous monomial localization with a
nontrivial Z intersection, diagonal shear entries, and tails outside the
final core. These are checks of the proof mechanisms; the arbitrary-width
degree and uniform radius claims follow from the argument above. The localization fixture
need not be C3 because sections 2--4 do not use that hypothesis.
