# Affine Pauli-orbit copy: a Clifford-dual obstruction and a dirty-mask gauge

[Encoded dirty copy](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Noncommuting copy grading](DIRTY_COPY_GRADING.md) ·
[Native copy interactions](NATIVE_COPY_INTERACTION.md) ·
[Field-shear source](FIELD_SHEAR_INTERTWINER.md)

An affine Pauli query cannot supply a complete dirty indicator through
a copy whose output-Z pullbacks are Clifford, even with an arbitrary
non-Clifford seed. For at least two address bits, its uniform copy error
is at least one quarter. This includes every native one-band copy.
A native one-bit endpoint shows why the address-size qualification matters.

This is a query-and-copy-class restriction, not a general native-depth
lower bound. In particular, the field-shear query is not Pauli. A dirty
address mask can remove its nonzero character, but a paid gauge then
turns that masked family into an ordinary Pauli translation. The gauge
need not preserve one-band depth of the copy.

All identities concern complete carrier, output and arbitrary private-helper
spaces. Algebraic products act rightmost first; gate lists are chronological.

## 1. Interface and scope

Let $`x\in\mathbb F_2^p`$, $`m=2^p`$, and let R be any unitary on a
carrier and its private helpers. Let $`Q_x`$ be Paulis whose binary
labels are affine in x. Their scalar phases are arbitrary; the family
need not be a literal representation. Set

```math
R_x=Q_xRQ_x^\dagger.
```

Let Y consist of m arbitrary output qubits. The fixed copy C may mix
the carrier, outputs and helpers, but does not act on the address. Its
uniform full-space operator error is

```math
\eta=\max_x\|C(R_x\otimes I_Y)C^\dagger
                     -R_x\otimes X_{Y_x}\|.          \qquad(1)
```

Assume that every output dual

```math
H_j=C^\dagger Z_{Y_j}C                              \qquad(2)
```

is Clifford. These duals are automatically Hermitian involutions.
The assumption holds for any native copy of T-depth at most one:
a physical T/T-dagger band maps a Pauli, under conjugation, to a tensor
product of single-qubit Clifford operators, and Clifford walls preserve
that property. Deeper copies satisfying (2) are also included.

**Theorem.** Under these hypotheses, $`p\ge2`$ implies $`\eta\ge1/4`$.
The theorem makes no claim about a non-Pauli field query or an arbitrary
two-band copy.

## 2. A central four-point identity

Fix a Hermitian Clifford involution H and define
$`H_x=Q_x^\dagger H Q_x`$. For every pair of addresses a and b,

```math
H_aH_0H_bH_{a+b}=\lambda_{a,b}I,
\qquad |\lambda_{a,b}|=1.                           \qquad(3)
```

The scalar need not be one. For example, $`H=H_u`$, with query
generators $`Q_{10}=X_u`$ and $`Q_{01}=Z_u`$, gives $`-I`$ in (3).

Absorb the affine label offset into $`H_0`$, and let q(x) denote the
remaining linear Pauli label. Let S be the binary symplectic action of
$`H_0`$. Since $`H_0^2=I`$, $`S^2=I`$. The operator
$`P_x=H_xH_0`$ is a scalar times a Pauli with label

```math
\ell(x)=(I+S)q(x).
```

This label is linear in x and satisfies $`S\ell(x)=\ell(x)`$. Therefore

```math
H_aH_0H_bH_{a+b}=P_aP_b(H_0P_{a+b}H_0)
```

has label $`\ell(a)+\ell(b)+S\ell(a+b)=0`$ and is a scalar Pauli.
This proves (3), retaining its literal phase. Shifting the address
origin gives the same statement for every affine two-dimensional plane.
Neither independence of the query generators nor commutation of the
carrier family is assumed.

## 3. Robust copy obstruction

Conjugating $`Z_{Y_j}`$ by the two unitaries compared in (1) changes
the operator by at most twice their distance. Pulling this comparison
through C gives

```math
\|R_xH_jR_x^\dagger-(-1)^{[x=j]}H_j\|\le2\eta.
```

Here and below carrier operators are tensored with output identity.
Conjugate by $`Q_x^\dagger`$ to obtain

```math
\|R H_{j,x}R^\dagger-(-1)^{[x=j]}H_{j,x}\|\le2\eta,
\qquad H_{j,x}=Q_x^\dagger H_jQ_x.                  \qquad(4)
```

For $`p\ge2`$, choose independent a and b with
$`j\in\{0,a,b,a+b\}`$. The four requested one-hot signs have product
minus one. Multiply the four comparisons (4) in the order
a, zero, b, a+b. Telescoping costs at most $`8\eta`$. Equation (3)
says conjugation by R leaves the original product $`\lambda I`$
unchanged, whereas the requested product is $`-\lambda I`$. Thus

```math
2=\|\lambda I-(-\lambda I)\|\le8\eta,
```

proving the theorem. More generally, exact copying of a sign table
requires even parity on every affine two-dimensional plane. Such
Boolean functions are affine; a singleton is not affine for $`p\ge2`$.

## 4. The one-bit endpoint has a native one-band copy

One carrier qubit u, two output qubits $`Y_0,Y_1`$, and no initialized
workspace suffice when p is one. Take

```math
Q_0=I,\quad Q_1=Z_u,\quad R_0=H_u,\quad R_1=Z_uH_uZ_u.
```

The two carrier involutions anticommute. Define the Clifford $`C_P`$
by the chronological word

```text
H(Y0), CX(Y0,u), CX(u,Y1), H(Y0).
```

Its action is $`Z_u\mapsto Z_uX_{Y_0}`$,
$`X_u\mapsto X_uX_{Y_1}`$, and hence
$`Y_u\mapsto Y_uX_{Y_0}X_{Y_1}`$.
Write $`T(P)=\exp(i\pi(I-P)/8)`$ and set $`U=T(Y_u)`$.
Then $`UZ_uU^\dagger=R_0`$ and $`UX_uU^\dagger=-R_1`$.
Consequently

```math
C=UC_PU^\dagger
 =T(Y_u)T(Y_uX_{Y_0}X_{Y_1})^\dagger C_P           \qquad(5)
```

satisfies $`CR_xC^\dagger=R_xX_{Y_x}`$ on every input. The two
T axes in (5) commute and are independent. A literal native diagonalizer
D has chronological gates

```text
CX(u,Y0), CX(Y1,Y0), H(u), S(u), H(Y0), H(Y1).
```

It sends $`Z_u\mapsto Y_u`$ and
$`Z_{Y_0}\mapsto Y_uX_{Y_0}X_{Y_1}`$. Thus C has the native word

```text
C_P, actual_inverse(D), T(u), TDG(Y0), D.
```

There are two T-type gates in one physical band, 20 total gates after
expanding the inverse S into three S gates, and no helper wires. The
actual inverse has the same ledger. In particular C is not an involution:

```math
C^2=H_{Y_0}\,\mathrm{CX}(Y_0;Y_1)\,H_{Y_0}\ne I.
```

The coherent source $`\mathcal A=QH_uQ^\dagger`$, with
$`Q=\mathrm{CZ}(x,u)`$, is Clifford. The indicator echo has chronological
word $`\mathcal A^\dagger,C^\dagger,\mathcal A,C`$: four T-type
gates in at most two bands, returning the carrier on every input.
This is a finite native endpoint, not an asymptotic compiler.

## 5. A dirty address-mask gauge

Suppose a paid coherent family obeys the literal additive relations

```math
Q_xQ_t=Q_{x+t},\qquad Q_0=I.
```

Add p arbitrary dirty label qubits t, and define

```math
Q'_x=X_t^x\otimes Q_x,\qquad
B=\sum_t|t\rangle\langle t|\otimes Q_t.
```

On the complete space,

```math
B Q'_x B^\dagger=X_t^x\otimes I.                    \qquad(6)
```

Indeed, the work factor for an input label t is
$`Q_{t+x}Q_xQ_t^\dagger=I`$. For the existing field-shear compiler,
B is the same paid coherent query with t replacing its address, at the
same eight-band depth and polynomial resource ledger. Its private
helpers return exactly; t is never initialized.

The mask makes $`\mathrm{Tr}Q'_x=0`$ for nonzero x, but (6)
removes the field factor by this paid gauge. Conjugating both the seed
and copy by B turns a masked-field pair into a Pauli-translation pair,
with unchanged error and constant additional depth. The conjugated
copy need not remain one-band. Therefore the preceding obstruction
does not automatically exclude a one-band copy written before the
gauge. No native construction is supplied by the gauge itself.

## 6. Exact finite evidence

[The exact fixtures](../../tests/test_pauli_orbit_copy.py) check the
central identity with affine offsets, dependent and noncommuting query
generators, and arbitrary query phases, including the literal $`-I`$
example and the two-bit odd-sign contradiction. They verify the native
one-bit copy, its actual inverse, and every column of its coherent
indicator echo. The dirty-mask gauge is checked with an arbitrary
returned helper in a native one-bit word and independently over every
two-bit field-shear address, mask and bank input.

These checks validate finite identities and gate accounting. They do
not replace the analytic all-size obstruction or supply the general
field seed and copy interaction.
