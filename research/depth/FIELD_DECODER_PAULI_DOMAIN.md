# Pauli-domain requirements for a complete field decoder

[Flat four-address decoder](FLAT_FIELD_BANK_PILOT.md) ·
[Encoded-copy contract](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Pauli-domain localization](BOUNDED_DEFICIENCY_TRANSPORT.md#2-signed-localization-of-a-possibly-degenerate-domain)

A complete decoder from the additive field query to wire translation must
move almost every physical Pauli out of the Pauli group. At three or more
address bits, its absorbed T axes span at least $`2m-2`$ binary directions,
where $`m=2^p`$ is the translated bank size. This remains necessary when the
decoded helper register has an arbitrary Clifford action.

The conclusion concerns axis rank and T count. It gives no growing T-depth
bound and does not apply to a decoder constrained only on one physical probe.
In particular it does not improve R96 or pay an all-size decoder.

<a id="robust-pauli-domain"></a>

## 1. Full-space statement, including moving helpers

Let $`p\ge3`$, $`m=2^p`$, and let $`Q_x`$ be a literal additive field
shear of degree K, tensored with identity on h arbitrary helper qubits:

```math
Q_x|u,v,w\rangle=|u,v+M_xu,w\rangle,
\qquad M_x=\sum_{i=0}^{p-1}x_iM_i.
```

Every nonzero $`M_x`$ is invertible. The source acts on $`n=2K+h`$
qubits, and we require $`n\ge m`$. Its symplectic matrices satisfy
$`F_x=I+N_x`$ and $`N_xN_y=0`$. These conditions permit any injective
address subspace of the field, not only consecutive powers of its generator.

The target has an m-qubit bank whose wires are indexed by
$`z\in\mathbb F_2^p`$, and $`q=n-m`$ other arbitrary qubits. Let

```math
T_x=W_x\otimes R_x,\qquad W_x:z\longmapsto z+x,
```

where every $`R_x`$ is a Clifford. It may move the helper register; no
identity promise or additive representation law for R is needed here.
For an arbitrary unitary E define

```math
\eta=\max_x\|EQ_xE^\dagger-T_x\|_{\rm op},\qquad
H_E=\{v:EP(v)E^\dagger\text{ is a signed Hermitian Pauli}\}.
```

The [retained domain argument](BOUNDED_DEFICIENCY_TRANSPORT.md#2-signed-localization-of-a-possibly-degenerate-domain)
makes $`H_E`$ a binary linear subspace, with an injective symplectic
label map under conjugation. Let A be the binary span of the actual axes
after all Cliffords in a native word for E have been absorbed to one side.

**Theorem.** For arbitrary Clifford helper actions,

```math
\eta\lt\frac{\sqrt2}{6}
\quad\Longrightarrow\quad
\dim H_E\le2+2q,\qquad \dim A\ge2m-2,\qquad T(E)\ge2m-2.
                                                        \qquad(1)
```

If every target query maps every Hermitian Pauli to a commuting Pauli,
the same conclusion holds under the stronger threshold

```math
\eta\lt\frac13.                                      \qquad(2)
```

Identity helper action and a residual field shear both satisfy this
condition. All comparisons concern the complete operator, with arbitrary
dirty inputs and reference entanglement. The theorem makes no assertion
that an approximately transported Pauli domain is query-invariant.

## 2. Literal source derivatives and the target gap

Use the Hermitian section $`P(a,b)=i^{a\cdot b}X^aZ^b`$, with the
dot product in the exponent taken over the integers. Field CNOT queries
map the ordered section $`X^aZ^b`$ without an extra phase. Their common
fixed-label space is

```math
F=\mathrm{span}\{Z_u,X_v\}\oplus V_h,
\qquad \ker N_x=F\quad(x\ne0),                       \qquad(3)
```

where $`V_h`$ is the full helper Pauli space. Every fixed label is fixed
literally, including its sign.

For any Hermitian Pauli P put $`P_x=Q_xPQ_x^\dagger`$. The derivative
$`P_xP`$ has label $`N_xv\in F`$, hence is fixed literally by every
query. Moreover P commutes with $`P_x`$: writing
$`v=(a_u,a_v,b_u,b_v)`$, their pairing is

```math
\langle v,N_xv\rangle
=a_u\cdot M_x^Tb_v+b_v\cdot M_xa_u=0.                 \qquad(4)
```

Thus unequal source Paulis of this form have operator distance two.

Fix $`v\in H_E`$, and write $`S=EPE^\dagger`$ and
$`U_x=EQ_xE^\dagger`$. Define

```math
A_x=U_xSU_x^\dagger S=E(P_xP)E^\dagger,\qquad
B_x=T_xST_x^\dagger S.
```

The first is fixed by every $`U_y`$. The second is a scalar times a
Pauli, since the target is Clifford. Unitary telescoping gives

```math
\|A_x-B_x\|\le2\eta,\qquad
\|T_yB_xT_y^\dagger-B_x\|
\le2\eta+2\|A_x-B_x\|\le6\eta.                      \qquad(5)
```

Two unequal Pauli operators with the same overall scalar convention
have distance at least $`\sqrt2`$. Hence (1) forces every $`B_x`$ to
be fixed literally by every target query. In the commuting-image case,
$`B_x`$ is Hermitian and commutes with its target conjugate; the gap is
two, giving (2). The latter condition holds for involutive wire
permutations and field shears separately, and for their tensor product.

## 3. The helper quotient removes the apparent escape

Let $`\mathcal U`$ be the image-label space of $`H_E`$. At the label
level, the conclusion of (5) says that every $`s\in\mathcal U`$
has zero second target differences. In particular its two bank masks
obey

```math
(\Pi_y-I)(\Pi_x-I)a=0,\qquad
(\Pi_y-I)(\Pi_x-I)b=0\quad\text{for every }x,y.
```

First translation differences are therefore constants, so

```math
a(z)=a(0)+\alpha\cdot z,\qquad
b(z)=b(0)+\beta\cdot z.                              \qquad(6)
```

This conclusion concerns only the bank; the helper masks need not be
fixed. Let $`F_R`$ be the subspace of helper Pauli labels whose operators
are **literally** fixed by every $`R_x`$. It is linear because a product
of literally fixed Pauli operators is literally fixed, including its
scalar multiplication phase. Put $`f=\dim F_R`$ and project

```math
\pi:\mathcal U\longrightarrow V_q/F_R.                \qquad(7)
```

For $`s\in\ker\pi`$, the helper Pauli is fixed. The two linear forms
in (6) have a common nonzero kernel element x because $`p\ge3`$.
Thus $`T_xST_x^\dagger=S`$ literally. Comparing with the source gives

```math
\|Q_xPQ_x^\dagger-P\|
=\|U_xSU_x^\dagger-S\|\le2\eta\lt2.
```

The source gap and (3) fix P under every query. Comparing again with
the target fixes S under every target query: $`2\eta\lt\sqrt2`$
in (1), and the stronger commuting gap is also available in (2).

The full target fixed space has dimension $`2+f`$: both bank masks
are constant, and the helper label belongs to $`F_R`$. Therefore

```math
\dim\ker(\pi|_{\mathcal U})\le2+f,\qquad
\dim\mathcal U\le(2+f)+(2q-f)=2+2q.                  \qquad(8)
```

This proves the domain bound. It also explains why simply requiring an
individual moving-helper Pauli to have a nonzero target stabilizer would
be invalid: a nonfixed residual-field Pauli has none. Only the kernel
in (7) needs such a stabilizer.

Every Pauli commuting with all absorbed axes remains a Pauli after E,
so $`A^\perp\subset H_E`$. Nondegeneracy of the n-qubit pairing gives

```math
2n-\dim A\le2+2q,\qquad \dim A\ge2(n-q)-2=2m-2.
```

The T count is at least the span rank of its individual axes, proving
the remaining conclusions of (1).

<a id="identity-and-residual-helper-actions"></a>

## 4. Identity helpers and a residual field query

For the original field degree $`K=m/2`$ and identity target helpers,
$`q=h`$. Every element of $`H_E`$ is then source-fixed and its image
is target-fixed, not merely dimension bounded. At error below one-third,

```math
\dim H_E\le2+2h,\qquad
\mathrm{rank}_{\rm symp}H_E\le2h.                   \qquad(9)
```

The two bank fixed labels, all-X and all-Z, commute because m is even.
The boundary $`p\ge3`$ is essential to this proof: at p equal to two,
the masks $`a(z)=z_0`$ and $`b(z)=z_1`$ have no common nonzero
translating stabilizer. The finite four-address decoder is consistent
with (9).

For an enlarged field degree $`K=m/2+\delta`$, with
$`\delta\ge p`$, consider the allowed target

```math
T_x=W_x\otimes Q_{\delta,x}\otimes I_h.
```

The residual full-rank field shear acts on $`2\delta`$ helper qubits.
The source and target characters agree: for nonzero x both traces are
$`2^{m/2+\delta+h}`$. This supports an abstract simultaneous
intertwiner, but does not pay its native word. In this target, (2) gives

```math
\dim H_E\le2+4\delta+2h,\qquad
\mathrm{rank}_{\rm symp}H_E\le4\delta+2h,\qquad
\dim A\ge2m-2.                                      \qquad(10)
```

For the symplectic rank statement, all affine bank masks in (6) pair
to zero at $`p\ge3`$: the summed product has degree at most two on
the p-dimensional binary cube. Only the helper coordinates contribute.
Thus moving Clifford helpers do not remove the extensive axis requirement.
Their actual motion and inverse calls must still be charged in a complete
indicator echo.

<a id="independent-flat-cells"></a>

## 5. Independent finite cells cannot give a growing decoder

The [six-T four-address decoder D](FLAT_FIELD_BANK_PILOT.md) has exact
Pauli domain

```math
H_D=\mathrm{span}\{X_0,X_2,X_3,X_1Z_0\}.
```

It has dimension four and symplectic rank two; each of these Paulis is
fixed with positive sign. The domain equals the symplectic orthogonal
complement of D's absorbed-axis span.

Let $`m=4r=2^p`$ with $`p\ge3`$. A tensor product of r such cells
and q untouched physical helpers has a Pauli domain of dimension at
least $`4r+2q`$. Arbitrary fixed input and output Clifford walls preserve
that dimension, even when they mix cells and helpers. It exceeds the
$`2+2q`$ permitted by (1). Consequently this entire independent-cell
class fails full field-to-bank intertwining below $`\sqrt2/6`$ for
arbitrary Clifford target-helper actions, and below one-third for identity
or residual-field helper actions.

Coupling the non-Clifford bands can change this domain and remains open.
The necessary axis count is compatible with the desired linear T budget;
neither this failure nor the count bound rules out constant-depth coupled
decoders. The separate [two-band screen](TWO_BAND_FIELD_DECODER_OBSTRUCTION.md)
also requires any successful two-band word to have intersecting absorbed
axis spans; it uses a different sparse-image argument.
The [focused fixtures](../../tests/test_field_decoder_pauli_domain.py)
check source phases, second differences, the helper quotient's boundary
case, and the finite domain. The all-size conclusion is the proof above.
