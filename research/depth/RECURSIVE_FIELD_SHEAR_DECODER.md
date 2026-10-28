# A recursive finite-field decoder for the full lookup representation

[Field-shear source and contract](FIELD_SHEAR_INTERTWINER.md) ·
[Returned cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md) ·
[Polynomial-overhead discovery budget](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)

The full bank-permutation representation has an explicit all-size decoder
built from finite-field coordinate changes and Fourier transforms. The
construction recursively decodes two half banks and couples their common
translation action to the swap of the halves. It specifies every basis
column, including the exceptional zero sector, without choosing dense
joint eigenbases.

This is an algebraic construction, not a shallow native compiler. Two
punctured, field-dependent coordinate permutations still need paid native
words. Compiling the recursion one level at a time also retains a depth
cost at every level. The query itself is fully paid below, even for the
dense multiplication matrices of a field tower. R96 is unchanged.

The later [projective recoupler](PROJECTIVE_FIELD_RECOUPLER.md) jointly compiles
the two conditional coordinate maps at p equal to three, with a different
allowed zero fiber. Its full eight-wire decoder costs 172 T gates in at most
83 bands and has no extra helper. The projective identity is general, but
enumerating its fibers grows exponentially in bank width; the uniform native
obligation and the separate-level depth limitation below remain.

## 1. A compatible tower and address pencil

Write $`K_1=\mathbb F_2`$. For $`j\ge2`$, define

```math
K_j=K_{j-1}[\theta_j]/(\theta_j^2+\theta_j+\delta_{j-1}),
\qquad
\delta_{j-1}=\prod_{h=2}^{j-1}\theta_h,
\qquad [K_j:\mathbb F_2]=2^{j-1}.                    \qquad(1)
```

An empty product is one. These quadratic polynomials are irreducible.
Indeed the relative trace of each $`\theta_h`$ is one. Repeatedly applying
the trace tower identity gives absolute trace one for
$`\delta_{j-1}`$. In a binary field, the image of $`z\mapsto z^2+z`$
is exactly the trace-zero hyperplane: the map has kernel of size two,
its image has half the field's size, and every image has trace zero.
Thus the displayed quadratic has no root in $`K_{j-1}`$.

For p address bits, put $`m=2^p`$, $`k=2^{p-1}`$ and take the two
k-qubit source registers U and V to be elements of $`K_p`$. Use

```math
\alpha_0=1,\qquad \alpha_i=\theta_{i+1}\quad(1\le i\lt p),
\qquad
\beta_p(x)=\sum_{i=0}^{p-1}x_i\alpha_i,
\qquad
Q^{(p)}_x|U,V\rangle=|U,V+\beta_p(x)U\rangle.       \qquad(2)
```

The address coefficients are binary-linearly independent: at each stage
the new $`\theta_j`$ lies outside the preceding field, which contains
all earlier coefficients. Hence every nonzero address multiplies by a
nonzero field element, and its matrix is invertible. This is a different
classical choice of the full-rank pencil from the consecutive-power
pencil in the [source proof](FIELD_SHEAR_INTERTWINER.md#1-a-full-rank-linear-pencil-and-one-native-shear).
Its full characters, complete-action contract and indicator wrapper are
identical. No conversion between the two choices is assumed free.

The decomposition $`x=(x',b)`$ obeys

```math
\beta_p(x',b)=\beta_{p-1}(x')+b\theta_p.              \qquad(3)
```

The field tower and all its multiplication matrices are classical circuit
descriptions. No field-element state is initialized or loaded.

### An eight-band query for any binary multiplication matrices

Let $`M_i`$ be multiplication by $`\alpha_i`$ in the chosen binary
basis. Ordinary Hadamards on V turn the coherent query into the cubic
phase

```math
(-1)^{\sum_{i=0}^{p-1}\sum_{j=1}^{k}
                   x_i U_j\ell_{ij}(V)},
\qquad
\ell_{ij}(V)=(M_i^{\mathsf T}V)_j.                  \qquad(4)
```

Each of the kp summands is a product of three binary linear forms on
disjoint principal registers. The returned parity-cube word applies
unchanged. When it calls for adding the third control to its private
arbitrary mask, add all bits in the support of $`\ell_{ij}`$ by CNOTs.
They are consecutive Clifford gates between the same physical T bands;
no clean copy of the linear form is needed. The closed walk returns the
mask and has exactly the intended cubic sign on every input.

Using one private mask per summand gives the forward and actual reversed
word bounds

```math
T=8kp,\qquad D_T\le8,\qquad h=kp,\qquad
G_{\rm all}\le8kp+8k^2p+2k,\qquad w_{\rm clean}=0.  \qquad(5)
```

The gate bound conservatively charges k CNOTs for each of the eight
parity updates per mask, as well as both Hadamard walls. Thus dense
tower matrices still give $`T,h=O(mp)`$ and
$`G_{\rm all}=O(m^2p)`$, within the retained discovery budgets. This
bound does not pay any variable field-coordinate transformation used
by the decoder below.

## 2. The intermediate representation that is actually coupled

For p=1, take $`E_1=\mathrm{CX}(V;U)`$, as in the one-address
[native word](FIELD_SHEAR_INTERTWINER.md#4-the-one-address-native-word).
For $`p\ge2`$, write $`F=K_{p-1}`$, $`K=K_p`$,
$`|F|=q=2^\ell`$, and $`\ell=2^{p-2}`$. Suppose the two half-bank
decoders $`E_{p-1}`$ are already specified. Before applying them, the
required representation on two pairs $`(u_1,v_1),(u_2,v_2)\in F^2`$
is generated by

```math
A_{x'}=Q^{(p-1)}_{x'}\otimes Q^{(p-1)}_{x'},
\qquad
S:(u_1,v_1)\leftrightarrow(u_2,v_2).                 \qquad(6)
```

The two half decoders are identical, so their tensor product commutes
with S. Conjugating (6) by them gives the physical lower-bit bank
translations and the physical upper-bit half swap, respectively.

For a binary field J, use the trace-pairing Fourier transform

```math
\mathcal H_J|z\rangle
  =|J|^{-1/2}\sum_{w\in J}(-1)^{\mathrm{Tr}_{J/\mathbb F_2}(zw)}|w\rangle.
                                                               \qquad(7)
```

It is a literal self-adjoint involution. The trace pairing is
nondegenerate, so in any binary basis it is a bank of ordinary
Hadamards and an invertible binary linear map. It is therefore
Clifford, with $`O((\log_2|J|)^2)`$ elementary gates. These transforms
act on arbitrary input registers; they do not prepare Fourier states.

Apply $`\mathcal H_F`$ to both v coordinates and denote the resulting
pairs by $`a=(u_1,w_1)`$ and $`c=(u_2,w_2)`$. Put

```math
d=a+c,\qquad q_0(a)=a_0a_1,\qquad
B(a,d)=a_0d_1+a_1d_0.
```

For $`\lambda=\beta_{p-1}(x')\in F`$, the lower translation is now
the phase

```math
(-1)^{\mathrm{Tr}_{F/\mathbb F_2}
  [\lambda(q_0(a)+q_0(c))]},
\qquad
q_0(a)+q_0(a+d)=B(a,d)+q_0(d),                         \qquad(8)
```

whereas S sends $`a\mapsto a+d`$ and fixes d. In particular, every
generator in (6) is the identity on the entire sector $`d=0`$.
This is an intermediate sector after the recursive half decoders have
been pulled back. It is not the fixed physical Bell zero sector
excluded by the [direct-extension obstruction](FIELD_SHEAR_INTERTWINER.md).

## 3. A full-coordinate recoupler, including zero

Write $`\theta=\theta_p`$ and $`\delta=\delta_{p-1}`$, so
$`\theta^2+\theta+\delta=0`$. Identify
$`U=d_0+d_1\theta`$ with $`d=(d_0,d_1)`$. Multiplication by theta
and its associated determinant are

```math
C=\begin{pmatrix}0&\delta\\1&1\end{pmatrix},
\qquad
\nu(d)=B(d,Cd)=d_0^2+d_0d_1+\delta d_1^2
       =\mathrm{Norm}_{K/F}(U).                       \qquad(9)
```

Thus $`\nu(d)\ne0`$ precisely when $`d\ne0`$, and the two columns
d and Cd then form an F-basis. For $`a=r d+sCd`$, (8) and the swap
take the particularly simple form

```math
q_0(a)+q_0(a+d)=\nu(d)s+q_0(d),\qquad
S:r\longmapsto r+1.                                  \qquad(10)
```

Here r and s are whole F registers, not single qubits. The translation
of r by the field element one is diagonalized by $`\mathcal H_F`$,
with eigenvalue $`(-1)^{\mathrm{Tr}_{F/\mathbb F_2}(\widehat r)}`$.

The following chronological word defines a unitary $`R_p`$ on all
$`4\ell=m`$ bank qubits.

1. Apply $`\mathcal H_K`$ to V. Reuse V for its Fourier coordinate and
   define the F-valued linear functions of that coordinate

   ```math
   \phi_0=\mathrm{Tr}_{K/F}(UV),\qquad
   \phi_1=\mathrm{Tr}_{K/F}(\theta UV).                \qquad(11)
   ```

2. Apply the classical, d-controlled affine permutation $`D_p`$ of V:

   ```math
   (\widehat r,s)=
   \begin{cases}
   (\phi_1,(\phi_0+q_0(d))/\nu(d)),&d\ne0,\\
   (V_0,V_1),&d=0,
   \end{cases}
   \qquad V=V_0+V_1\theta.                            \qquad(12)
   ```

3. Apply $`\mathcal H_F`$ to the first register $`\widehat r`$,
   unconditionally, and call the output r.

4. Apply the d-controlled linear permutation $`F_p`$:

   ```math
   a=
   \begin{cases}
   r d+sCd,&d\ne0,\\
   (r,s),&d=0.
   \end{cases}                                      \qquad(13)
   ```

5. Replace the retained d by $`c=a+d`$, and regard the two outputs as
   $`a=(u_1,w_1)`$ and $`c=(u_2,w_2)`$. Apply $`\mathcal H_F`$
   to $`w_1,w_2`$ to obtain $`(u_1,v_1),(u_2,v_2)`$.

Each step is a bijection or a Fourier transform on existing registers.
For nonzero U, multiplication by U is invertible. The F-linear map
$`Z\mapsto(\mathrm{Tr}_{K/F}Z,\mathrm{Tr}_{K/F}(\theta Z))`$
is also invertible: its two functionals are independent by
nondegeneracy of the relative trace pairing and independence of
$`1,\theta`$. Therefore (12) is a bijection on V for every nonzero d.
Equation (9) proves the same for (13). Their zero branches are explicit
fixed bijections on the whole remaining register, not clean-state
conventions. All maps have literal computational phase one.

### Proof of simultaneous intertwinement

After the first Fourier wall, the source generator with address
$`(x',b)`$ is the diagonal phase

```math
(-1)^{\mathrm{Tr}_{F/\mathbb F_2}
 [\beta_{p-1}(x')\phi_0+b\phi_1]}.                   \qquad(14)
```

For $`d\ne0`$, (12) makes its first exponent exactly the target
exponent in (10), while the second is the Fourier eigenvalue of
$`r\mapsto r+1`$. Thus the lower generators and S have the required
action under the word. For $`d=0`$, U is zero, so every source generator
is the identity; every target generator in the intermediate coordinates
is also the identity by (8). Consequently the unconditional Fourier
wall and the arbitrary but bijective zero branches preserve all
$`q^2`$ states of this sector.

The complete operator identity is therefore

```math
R_p Q^{(p)}_{(x',b)}R_p^\dagger=A_{x'}S^b
\quad\text{for every }(x',b),\qquad
E_p=(E_{p-1}\otimes E_{p-1})R_p,
\qquad
E_pQ^{(p)}_xE_p^\dagger=W_x.                          \qquad(15)
```

The last equation follows inductively from (6). All coordinates have
been specified, so this is a terminating algebraic rule at every size,
with no spectral-basis choice. The actual inverse is the reversed word
using the inverses of (12) and (13), the inverse coordinate addition,
the self-adjoint Fourier walls, and the two recursive actual inverses.
The identities hold on arbitrary coherent addresses and arbitrary bank
states, and hence on all external references.

The final half decoders redistribute the intermediate d=0 subspace
among physical translation sectors. Keeping d=U within $`R_p`$ therefore
does not retain the invalid fixed physical Bell-zero assignment.

## 4. What remains to be compiled and what the recursion costs

Equations (12) and (13) define finite classical permutations. Calling
them permutations does not supply efficient reversible circuits for
them. In particular, the norm division, multiplication by a variable
field element, and the whole zero branch cannot be treated as free
coordinate notation, a prepared nonzero flag, or a promised nonzero
input. A native implementation must return every arbitrary helper and
pay both forward and actual inverse words.

Suppose those two complete controlled permutations at level p have
per-direction T counts $`t_{D,p},t_{F,p}`$, T-band bounds
$`d_{D,p},d_{F,p}`$, total gates $`g_{D,p},g_{F,p}`$ and arbitrary
returned helper requirements $`h_{D,p},h_{F,p}`$. Running the two
recursive half decoders in parallel gives the conservative recurrences

```math
T_p\le2T_{p-1}+t_{D,p}+t_{F,p},
\qquad \Delta_p\le\Delta_{p-1}+d_{D,p}+d_{F,p},
```
```math
G_p\le2G_{p-1}+g_{D,p}+g_{F,p}+O(k^2),
\qquad h_p\le\max\{2h_{p-1},h_{D,p},h_{F,p}\}.        \qquad(16)
```

The base values are $`T_1=\Delta_1=h_1=0`$ and $`G_1=1`$. Helper reuse
in (16) is valid only because each supplied primitive returns its
entire helper register on arbitrary inputs. The same recurrences pay
the actual inverse. All Fourier and fixed-coordinate costs appear in
the $`O(k^2)`$ Clifford term.

For example, per-level $`T,h=O(kp^c)`$ and
$`G=O(k^2p^c)`$ for a fixed c would give decoder
$`T,h=O(mp^{c+1})`$ and $`G=O(m^2p^c)`$ by summing the binary
recursion. This is a conditional ledger, not a claim that the
controlled permutations already have those bounds.

These recurrences describe exact primitives. If the two implemented
coordinate changes at level j instead have complete operator errors
$`\epsilon_{D,j},\epsilon_{F,j}`$, a telescoping bound for the complete
decoder is

```math
\|\widetilde E_p-E_p\|
\le\sum_{j=2}^{p}2^{p-j}
                   (\epsilon_{D,j}+\epsilon_{F,j}).              \qquad(17)
```

Consequently a fixed per-cell approximation tolerance cannot be used
unchanged as p grows. Any approximate compiler must allocate accuracy
over all the recursive cells, and include the resulting synthesis and
work-return costs. No approximate native ledger is inferred here.

### Independent compilation at each level does not improve depth

Every $`R_p`$ for $`p\ge2`$ is non-Clifford. On binary Pauli labels
every pair of source deviations $`N_x=F_x-I`$ has $`N_xN_y=0`$,
as in the [field-shear boundary](FIELD_SHEAR_INTERTWINER.md#clifford-simultaneous-conjugacy-is-impossible-for-p-at-least-two).
For any nonzero lower address, the target deviation of
$`A_{x'}`$ is $`\mathrm{diag}(N,N)`$ with $`N\ne0`$, while the
half-swap deviation has block matrix

```math
N_S=\begin{pmatrix}I&I\\I&I\end{pmatrix}.
```

Their product has every block equal to N and is nonzero. Simultaneous
Clifford similarity would preserve the zero-product property, a
contradiction. In particular, each separately implemented exact
recoupler requires a non-Clifford band. The separate-stage depth sum
in (16) has at least $`p-1`$ positive contributions.

This is a limitation of that compilation ledger, not a lower bound on
the resynthesized $`E_p`$ or on unrestricted frame compilation. A depth
improvement through this construction would require a joint native
implementation across levels, rather than merely faster isolated
field-arithmetic calls. No such implementation is supplied here.

## 5. A paid multiplication-XOR component

One substantial arithmetic component can be paid without clean
workspace. For any fixed binary field J of degree k, there is an exact
native word for

```math
|U,V,Z,h\rangle\longmapsto|U,V,Z+UV,h\rangle          \qquad(18)
```

on three arbitrary k-qubit registers and arbitrary returned helpers.
All field bases are fixed classical data. In particular, a fixed
binary-linear field isomorphism converts a degree-k tower field to a
polynomial basis on each of the three existing registers and back.
These are invertible binary linear maps, with no T cost, no clean
copies and $`O(k^2)`$ Clifford gates. The Fourier proof below uses the
target's dual coordinates in that chosen basis.

In a polynomial basis $`J=\mathbb F_2[t]/(f)`$, put
$`n=2k-1`$ and $`e=\max\{1,\lceil\log_2n\rceil\}`$.
Take an auxiliary field $`A=\mathbb F_{2^e}`$ and n distinct
classically specified points $`\zeta_i\in A`$. This auxiliary field
only describes binary linear forms; it does not introduce an input
register or a prepared state. The product of the degree-below-k
polynomials U and V has degree below n, so evaluation at these points
and Lagrange interpolation reconstruct it exactly.

Let $`L_i(t)\in A[t]`$ be the interpolation basis polynomials and
$`b_1,\ldots,b_e`$ a binary basis of A. Choose a binary-linear
projection $`\pi:A\to\mathbb F_2`$ with $`\pi(1)=1`$. For each
triple $`(i,a,b)`$, the polynomial
$`L_i(t)b_ab_b`$, reduced modulo f and then projected coefficientwise
by pi, gives a fixed binary k-vector $`v_{iab}`$. The projection is
needed because individual interpolation summands have coefficients in
A. On the sum representing an actual product, the coefficients are
binary, so the projection preserves the answer exactly.

Writing $`u_{ia}(U)`$ for bit a of $`U(\zeta_i)`$ and similarly
for $`v_{ib}(V)`$, all of which are binary linear forms, gives

```math
Z\cdot(UV\bmod f)
  =\sum_{i=1}^{n}\sum_{a,b=1}^{e}
       u_{ia}(U)v_{ib}(V)\,[Z\cdot v_{iab}].          \qquad(19)
```

Thus the multiplication phase is a batch of at most
$`s=(2k-1)e^2`$ products of three linear forms on disjoint registers.
The same returned parity-cube argument as in (4) gives

```math
T\le8s=O(k\log^2(k+1)),\qquad D_T\le8,\qquad h\le s,
\qquad G_{\rm all}\le8s(k+1)+2k=O(k^2\log^2(k+1)).   \qquad(20)
```

The $`2k`$ Hadamards on Z turn the phase into (18). Arbitrary basis
changes, when required, add $`O(k^2)`$ Clifford gates to the displayed
polynomial-basis ledger. Each helper receives its own eight T phases;
all other evaluation, interpolation and reduction operations are fixed
parity updates between those phases. The actual inverse is the literal
reversed native word and has the same resource bounds.

This is an XOR oracle with both multiplicands preserved. It does not
provide an in-place overwrite by $`UV`$, an inverse, a nonzero flag,
or a register containing a known product. Substituting a borrowed
arbitrary register for a clean product register changes the action;
any use in (12) or (13) must pay that complete transformation separately.

## 6. Why a bounded number of affine field queries cannot compile one cell

There is a sharper screen for one particular implementation grammar.
Isolate the U-read-only cell

```math
K_U=F_p\,\mathcal H_F^{(\widehat r)}\,D_p.           \qquad(21)
```

For each fixed U this is an affine Clifford on the two F registers
that encode V. The fixed Fourier and coordinate walls outside the cell
do not pay its U dependence. The following statement concerns words
that keep U unchanged, whose query stages are Clifford on V and any
helpers with binary symplectic matrix entries affine in the k input
bits of U, interleaved with fixed Clifford gates on V and helpers.

For $`U=d_0+d_1\theta\ne0`$, the linear parts of (12) and (13) are

```math
A_U=\begin{pmatrix}
d_0+d_1&d_0+(\delta+1)d_1\\
d_1/\nu&(d_0+d_1)/\nu
\end{pmatrix},
\qquad
B_U=\begin{pmatrix}d_0&\delta d_1\\d_1&d_0+d_1\end{pmatrix}.
                                                               \qquad(22)
```

Use field Pauli labels $`X_a`$ and
$`Z_b|v\rangle=(-1)^{\mathrm{Tr}_{F/\mathbb F_2}(b\cdot v)}|v\rangle`$.
In these coordinates the partial Fourier transform swaps the first
field X and Z labels exactly. Ordinary binary labels differ by only
a fixed trace-dual basis change, preserving the maximum coordinate
degree. The XX block
of the cell's symplectic matrix is
$`B_U\,\mathrm{diag}(0,1)\,A_U`$. Its upper-left field entry is

```math
g(U)=\frac{\delta d_1^2}{\nu(U)}\quad(U\ne0),
\qquad g(0)=0.                                       \qquad(23)
```

The zero value follows from the specified identity zero branches,
for which the cell is just the partial Fourier transform. The affine
offset in (12) affects Pauli signs but not this symplectic entry.

Put $`q=2^\ell`$. Since $`d_1=U+U^q`$ and
$`\nu(U)=U^{q+1}`$, the function in (23), including U=0, satisfies

```math
g(U)/\delta=U^{q-1}+U^{q^2-q}.                       \qquad(24)
```

It has exact binary algebraic degree ell. Here is a direct proof of
the needed lower bound. For any nonzero $`a\in F`$, the binary
function $`\mathrm{Tr}_{F/\mathbb F_2}(a g(U)/\delta)`$ expands
into $`2\ell`$ monomials whose exponents are the distinct cyclic
rotations, in $`2\ell`$ binary positions, of a block of ell ones
followed by ell zeros. All coefficients are nonzero and every
exponent has binary weight ell. This is the unique representing
polynomial of degree below $`q^2`$ for that function on K.

A binary linear coordinate function of U is a linearized polynomial
$`\sum_j c_jU^{2^j}`$. A product of r coordinates can therefore
produce, after reduction modulo the finite-field function identity,
only nonzero exponents of binary weight at most r: addition of
powers of two and cyclic carrying do not increase that weight.
Thus no coordinate polynomial of degree below ell can represent
the preceding trace function. Conversely, a monomial whose exponent
has weight ell is a product of ell Frobenius-linear functions, giving
the matching upper bound. In particular, some binary symplectic entry
of (21) has degree ell.

A word with r affine-entry query stages has symplectic entry degree
at most r, since its matrices are multiplied with only constant
Clifford matrices between them. Equality to the cell, including
identity on arbitrary helpers, consequently requires

```math
r\ge\ell=2^{p-2}.                                  \qquad(25)
```

This statement has a constant operator-error gap within the same
grammar. Every fixed-U candidate block is Clifford, as is (21).
If its uniform operator error is below
$`\gamma=\sqrt{2-\sqrt2}`$, the retained
[Clifford trace gap](CLIFFORD_ORBIT_QUERY_RANK.md) forces its binary
symplectic action to agree at every U: a nonscalar relative Clifford
has distance at least gamma, independently of its scalar phase.
The same degree comparison proves (25). Identity helper blocks do
not remove the displayed entry.

This excludes bounded affine-query factorizations of the specified
cell. It is not a native T-depth lower bound. A stage with nonlinear
symplectic dependence on U, movement of U, a non-Clifford fixed-U block,
another choice of decoder, or joint synthesis across levels falls
outside the premises. A complete multiplication-XOR query (18) with
read-only U acts as a Clifford shear on V and Z, with symplectic
entries affine in U. Such calls are covered by the stage bound,
even though their internal native decompositions use non-Clifford
gates. Nonlinear address-dependent Pauli signs alone remain allowed
and do not defeat the degree bound. Any implementation outside this
grammar still requires a paid circuit with complete work return.

## 7. Verification and the next operation

`tests/test_recursive_field_shear_decoder.py` checks every column and
address of the recursive word for p=1, 2 and 3, including its actual
abstract inverse and every exceptional sector. It checks the cell's
symplectic coefficient and its algebraic degree through p=4. The dense
query's emitted phase core and actual inverse are exhausted over all
principal and dirty columns at p=2; at p=3, all principal columns are
checked with two opposite dirty masks, alongside the linear-form
identity. The multiplication tensor is checked at field degrees two
and four, with every principal and dirty column of the native phase
core and its actual inverse checked at degree two. The explicit
Hadamard walls turn these phase cores into the stated queries.

These are fixtures for the all-size algebraic proofs. They do not
compile (12) or (13), or certify an asymptotic native decoder bound.

The constructive task is now explicit: jointly synthesize the controlled
field coordinate changes across levels, including their complete zero
branches and arbitrary work return, or change the mechanism. Faster
separate cells retain the serial-stage cost, and a bounded number of
read-only affine queries cannot implement even one specified large cell.
The eight-band query and multiplication-XOR component are paid; the
decoder's native depth remains open.
