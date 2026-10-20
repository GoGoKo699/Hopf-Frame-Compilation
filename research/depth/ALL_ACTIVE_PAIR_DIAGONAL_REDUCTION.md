# All-active one-pair reduction to a diagonal native three-band word

[Original physical hypotheses and central insertion](ALL_ACTIVE_SINGLETON_CLOSURE.md) ·
[Full-space diagonal native-three theorem](THREE_LAYER_NATIVE_STRUCTURE.md) ·
[Rank-two endpoint case](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#5-all-paired-sources-with-two-clean-active-axes-at-both-endpoints)

The large feature constant in the first all-active closure is unnecessary.
Once the actual dressed pair axis has one computational X mask, native
independence identifies every endpoint rotation which dresses that axis.
Those rotations are already computationally diagonal and can be peeled.
The remaining axis is an actual Pauli. One fixed change of basis then
makes both the target and the pair insertion diagonal on the **entire**
clean-reduced logical space. The full-space diagonal three-band theorem
applies directly, without enumerating one-qubit matrix fibers.

## 1. Statement and scope

Keep exactly the original physical hypotheses of
[the all-active singleton closure, Section 1](ALL_ACTIVE_SINGLETON_CLOSURE.md#1-original-physical-hypotheses-and-uniform-bound):
one clean qubit, three actual native bands, one genuine original middle
pair, and pair difference equal to the sum of all actual clean-active
axes at each endpoint. The occupied action is a complete Boolean sign
operation on all logical and dirty inputs, with arbitrary references,
exact clean return, a retained common scalar, and an unrestricted unused
clean block. Original singletons and silent lists remain unrestricted.

Let k be the common active rank. It is even and at least two, as proved
by the actual active-axis pairing in
[the retained central-insertion reduction](ALL_ACTIVE_SINGLETON_CLOSURE.md#2-a-literal-central-controlled-quarter-turn).
Then

```math
\begin{array}{c|cc}
&\deg g&\deg f\text{ for a complete Boolean XOR query}\hline
k\ge4&\le5&\le4\\
k=2&\le6&\le5
\end{array}
```

Thus the uniform all-active one-pair bounds are six and five. The
k=2 row is the retained rank-two endpoint theorem; the proof below
establishes the k at least four row. It does not revise the generic
partial-diagonal feature theorem. It replaces the use of that theorem
in this particular physical source class.

## 2. Retained reduction through the actual one-mask Clifford

Use Sections 2--5 of
[the all-active singleton closure](ALL_ACTIVE_SINGLETON_CLOSURE.md).
These first perform the literal orientation repair, retain the common
physical clean coordinate, and restrict that coordinate only after all
retained internal factors preserve it. The shared computational-center
alignment keeps the occupied target Boolean diagonal. On the resulting
entire logical space the exact word is

```math
G=K_C^\dagger C^\dagger G_0 E A K_A,
```

where A,C,E are the actual independent native silent and singleton
lists. The Hermitian Pauli P_h commutes with every E axis: in the
[original-list reduction](ORIGINAL_PAIR_SINGLETON_COUPLING.md#2-literal-repair-and-the-common-active-span)
each retained positive singleton commutes B and every actual active A
axis, while h=b+beta a for one such active axis a. P_h also commutes
the genuine central Pauli group R. All three bands commute R. The
insertion is the literal full-space identity

```math
\begin{aligned}
G_0&=D_R\omega^{-\epsilon\eta}
 \left(\Pi_+S_{Q_RP_h}^{\epsilon}
                       +\Pi_-\omega^{\epsilon}I\right),\\
\Pi_\pm&=(I\pm\Gamma)/2,\qquad \omega=e^{i\pi/4}.
\end{aligned}
```

Here D_R is diagonal C3 in R coordinates, Q_R is a Hermitian diagonal
Clifford in R, and Gamma is an actual nontrivial central Pauli. In
particular, G_0 commutes P_h literally. The output image of R is a
computational-Z subgroup.

For k at least four,
[the nonaffine central quadratic and scalar-branch argument](ALL_ACTIVE_SINGLETON_CLOSURE.md#5-quadratic-signs-force-the-actual-dressed-axis-to-be-monomial)
proves that

```math
H_C=K_C^\dagger C^\dagger P_h C K_C=X(u)D
```

is an actual Hermitian Clifford with one computational X mask u;
D is an actual diagonal Clifford. This is the last result used from
the previous large-constant proof. Its generic matrix-feature step is
not used below.

Define also the actual Hermitian Clifford

```math
H_A=K_A^\dagger A^\dagger P_h A K_A.
```

Conjugating a Pauli through a commuting native band produces a Clifford:
each axis anticommuting with that Pauli contributes an S rotation, with
the native scalar retained. Since E and G_0 commute P_h, cancellation
in the displayed complete word gives

```math
G H_A G^\dagger=H_C.
```

Conjugation by the diagonal G cannot change the computational X mask
of any matrix entry. Therefore H_A also has exactly the mask u. This
obtains the input conclusion from the actual full target; no second
sector argument or extension into the unused clean block is needed.

## 3. Native independence peels all dressing of P_h

The following observation applies separately to A and C. Let

```math
B=\prod_j T_{P_j}^{\epsilon_j},\qquad
H=K^\dagger B^\dagger P_h B K,
```

with an actual independent commuting axis list. Suppose H is an actual
Clifford with only one computational X mask. Split B as B_0 B_1,
where B_1 contains exactly the axes anticommuting with P_h.

Because H has one mask, conjugation by H fixes the binary label of
every computational Z. The binary symplectic action of
B^\dagger P_h B is the product of the transvections on the B_1 axes.
For every computational Z label z, this gives

```math
\sum_{j:\langle p_j,h\rangle=1}
       \langle p_j,Kz\rangle p_j=0.
```

The p_j are independent. Hence each coefficient in this sum vanishes
for every z. Every K^\dagger P_j K from B_1 is therefore a signed
computational Pauli. Consequently

```math
\Delta_B=K^\dagger B_1 K
```

is the literal product of computational T/T-dagger phase rotations,
and is diagonal C3. The remaining B_0 commutes P_h literally. This
argument uses the actual independent axes, including their signs and
orientations; it does not replace a dependent commuting phase list by
a native list.

Apply this observation to the two actual one-mask Cliffords H_A,H_C.
Write A=A_0 A_1, C=C_0 C_1 and set

```math
\Delta_A=K_A^\dagger A_1K_A,\qquad
\Delta_C=K_C^\dagger C_1K_C.
```

Both Delta factors are diagonal C3. The exact rearrangement is

```math
\begin{aligned}
G'&=\Delta_C G\Delta_A^\dagger
   =K_C^\dagger C_0^\dagger G_0 E A_0K_A,\\
W&=K_C^\dagger C_0^\dagger E A_0K_A,\\
F&=K_C^\dagger G_0K_C,\qquad G'=FW.
\end{aligned}
```

The simplification of F uses the fact that C_0 commutes P_h and R.
The three bands of W are actual sublists of the retained native lists,
so W is an actual full-space native-three word. Define the actual
endpoint Paulis

```math
P_A=K_A^\dagger P_hK_A,\qquad
P_C=K_C^\dagger P_hK_C.
```

Every internal factor of G' and W commutes P_h. Therefore

```math
G'P_A(G')^\dagger=P_C,\qquad WP_AW^\dagger=P_C.
```

G' remains diagonal, although it need not remain a Boolean sign.
Thus P_A and P_C have the same computational X mask. It is u:
conjugating H_A by the peeled diagonal Delta_A does not change its
mask, and gives P_A. The corresponding statement holds at the output.

## 4. A fixed full-space change of basis makes the remainder diagonal

If u=0, both endpoint Paulis are diagonal. Diagonal G' then forces
P_A=P_C, with the sign included. The insertion F is already diagonal:
it is a function of the computational central R and P_C. Consequently
W=F^\dagger G' is a full-space diagonal native-three word. No additional
change of basis is needed in this case.

Suppose u is nonzero. There is an actual diagonal Clifford L with

```math
L^\dagger P_C L=P_A.
```

For completeness, write the two binary labels as (u,v_A),(u,v_C).
A diagonal Clifford implements any symmetric binary shear
(x,z) to (x,z+Bx). Choose symmetric B with Bu=v_A+v_C; such a B exists
for every nonzero u and every right-hand side. A diagonal Pauli then
corrects the signed action if necessary, since some computational Z
anticommutes with X(u). This constructs one fixed Clifford, without
conditioning a frame on a central character.

Put

```math
\widehat G=L^\dagger G',\quad
\widehat W=L^\dagger W,\quad
\widehat F=L^\dagger F L.
```

Then the actual identity is

```math
\widehat G=\widehat F\widehat W,\qquad
[\widehat G,P_A]=[\widehat W,P_A]=0.
```

The target hat G is diagonal, and the insertion hat F is the same
central-controlled quarter-turn, now about Q_R P_A. The computational
central Paulis are unchanged by the diagonal L.

Choose one fixed monomial Clifford J with

```math
J^\dagger P_AJ=X_t
```

for one actual logical coordinate t. A computational linear change
sends u to that coordinate, and a diagonal Clifford removes the
remaining Z part and sign. Conjugation by J preserves computational
diagonality. Each central Pauli in J^\dagger R J is computational Z
and commutes X_t, so it has no Z_t component.

Since J^\dagger hat G J is diagonal and commutes X_t, its phase is
independent of t. A Hadamard H_t therefore leaves this entire operator
unchanged. Define the full-space conjugates

```math
\begin{aligned}
\overline G&=H_tJ^\dagger\widehat GJH_t
            =J^\dagger\widehat GJ,\\
\overline F&=H_tJ^\dagger\widehat FJH_t,\\
\overline W&=H_tJ^\dagger\widehat WJH_t.
\end{aligned}
```

The insertion bar F is diagonal because the only noncentral Pauli
X_t has become Z_t. The operator bar G is diagonal by the preceding
independence of t. Thus

```math
\overline W=\overline F^\dagger\overline G
```

is diagonal on the entire logical space. It still has an actual
three-native-band implementation: L,J,H_t only change fixed endpoint
Cliffords of W, and the internal axes are unchanged. In particular,
there has been no restriction to an eigenvalue of P_h or X_t, no new
native relation, and no uncharged synthesis of a projected phase band.

## 5. Literal phase degree and restoration

[The full-space diagonal native-three theorem, Section 1](THREE_LAYER_NATIVE_STRUCTURE.md#1-exact-statement-and-conventions)
gives an ordinary integer multilinear phase exponent of degree at most
five, modulo eight, for bar W in these fixed computational coordinates.
A common native scalar may be retained outside that exponent. The
theorem applies also when the diagonal operator is not Boolean.

In these same computational coordinates write gamma for the affine
negative-eigenvalue bit of Gamma and q for the quadratic sign bit of
Q_R. For u nonzero, the literal exponent of bar F is

```math
p_F(z,t)=p_D(z)-\epsilon\eta
 +2\epsilon(1-\gamma(z))\bigl(q(z)\mathbin\oplus t\bigr)
 +\epsilon\gamma(z)\pmod8.
```

p_D is the retained diagonal-C3 exponent of D_R and has ordinary degree
at most three. The scalar minus branch and the prefactor are displayed
explicitly. For u=0, replace t in this formula by the affine sign bit of
the diagonal Pauli P_C; the same estimate holds.

This exponent has ordinary degree at most five. One way to check the
only nontrivial product is to write, modulo four, an affine Boolean
gamma as L_1-2L_2 with degrees at most one and two, and the quadratic
Boolean s=q XOR t as H_2-2H_4 with degrees at most two and four. Then

```math
2\gamma s
 =2L_1H_2-4(L_2H_2+L_1H_4)\pmod8
```

has degree at most five. Also 2s has degree at most four modulo eight,
and the integer lift of gamma has degree at most three modulo eight.
The estimates include affine constants and remain valid after an
invertible affine change of computational coordinates.

It follows that bar G has a phase exponent of ordinary degree at most
five in these same fixed J coordinates. Restore all endpoint phases
**in those coordinates**:

```math
J^\dagger GJ
 =(J^\dagger\Delta_C^\dagger L J)\,
       \overline G\,(J^\dagger\Delta_AJ).
```

For u=0 take J=L=I. The restored factors are diagonal C3 or diagonal
Clifford in these coordinates: Clifford conjugation preserves their
hierarchy level, while the monomial J preserves their diagonality.
Their ordinary phase degrees are at most three. The quadratic Boolean
phase from the earlier shared-center alignment may likewise be restored
here. All scalar factors remain included.

The restored target in these fixed coordinates is e^{i theta} times
the Boolean sign (-1)^{tilde g}, where tilde g is g composed with the
affine computational permutation of J. Subtract its phase at one
computational input to remove only the common scalar. The resulting
exponent is 4(tilde g(x) XOR tilde g(0)) modulo eight and has an integer
multilinear representative of degree at most five. Uniqueness of
multilinear coefficients modulo eight gives deg tilde g at most five.
Finally use affine invariance of **Boolean degree** to conclude

```math
\deg_{\mathbb F_2}g\le5\qquad(k\ge4).
```

This last step does not assert that an arbitrary ordinary degree-five
polynomial modulo eight remains ordinary degree five under affine
Boolean substitution. The native-three phase theorem and every
restoration factor were applied in one fixed coordinate system before
extracting the Boolean degree.

Hadamards on the output of a complete Boolean XOR query give the
complete sign (-1)^{z f(x)}, hence deg f at most four in this row.
Combining with the retained k=2 result gives the uniform six/five bound.

## 6. What the revision changes

The first large bound used a generic finite alphabet for a one-qubit
matrix core. In the present class that core is removable: actual native
independence makes its endpoint dressing diagonal, and the remaining
Pauli symmetry permits a fixed full-space diagonalization. This is a
structural replacement of the feature step, not an optimization of the
number of alphabet elements.

The proof preserves arbitrary dirty and reference inputs, all occupied
logical columns, exact clean return, and the unrestricted unused clean
block. It adds no initialized helper and no approximation error.
Fixed Clifford changes, if emitted, retain their usual
O((Q+1)^2) elementary-Clifford costs on the same register. The proof is
a necessary condition on an existing circuit, not a new upper compiler.

The same original-circuit approximation consequence remains available:
the initialized rounding theorem at error strictly below 2^{-530690}
first yields exact clean return and exact occupied Boolean action for
that same native circuit. If its actual original lists satisfy the
all-active one-pair hypotheses, the sharper degree bounds above then
apply. The source-class hypothesis remains essential; this does not
exclude all approximate three-band circuits or all fixed accuracies.

General original pairs, several original pairs, unrestricted one-clean
three-layer degree, and dimension-growing native-depth lower bounds
remain open. The upper resource frontier is unchanged.
