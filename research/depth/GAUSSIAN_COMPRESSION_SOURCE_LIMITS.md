# Gaussian compression with a small clean signal

Constant-amplitude Gaussian compression with a small clean signal still
has bounded logical Pauli-label and sign degree. This applies to the
[nilpotent-pencil source](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md) and to
direct indicator compression within the actual query grammar specified
below. It does not give an unrestricted native T-depth lower bound.

## Statement and conventions

Let the read-only address have p bits. For each computational address x,
let V_x be Clifford on an arbitrary bank B and a k-qubit signal A, with
k at least one. B includes every logical output and every dirty helper;
its dimension is unrestricted. Only A is initialized. Define

```math
\begin{aligned}
K_x&=(I_B\otimes\langle0^k|)V_x(I_B\otimes|0^k\rangle),\\
\|K_x-J_x/\sqrt2\|_{\rm op}
  &\le\varepsilon\lt\frac{\sqrt2-1}{2}.
\end{aligned}                                           \qquad(1)
```

where J_x is unitary on the complete bank B. The norm includes arbitrary
dirty inputs and external references. No intermediate return of A, and
no projection as a physical operation, is assumed.

Use Hermitian tensor Paulis
$`P(a,b)=\bigotimes_j i^{a_jb_j}X_j^{a_j}Z_j^{b_j}`$.
Suppose the physical binary symplectic entries of V_x have Boolean
address degree at most D, and the signs of its Hermitian generator
images have degree at most E. All degree statements use algebraic normal
form over $`\mathbb F_2`$.

Then $`U_x=\sqrt2K_x`$ is exactly Clifford and
$`\|U_x-J_x\|\le\sqrt2\varepsilon`$. Its logical symplectic entries
have degree at most 2kD, and its logical Hermitian generator signs have
degree at most

```math
2kD+\max(E,4D)\le E+(2k+4)D.                            \qquad(2)
```

These bounds do not depend on the number of dirty wires. The weaker
label bound $`D(k+1)^2`$ follows without using rank one and is proved
below as well.

For q literal affine-address XOR shear queries interleaved with arbitrary
fixed, address-independent Cliffords on B and A, the physical degrees
satisfy $`D\le q`$ and $`E\le3q`$. Therefore a compressed source for
the addressed inverse-column observable requires $`2kq\ge p`$.
A direct complete indicator compression, whose information appears in
Pauli signs rather than its symplectic matrix, requires the conservative
bound $`(2k+7)q\ge p`$.

## Invertibility fixes the compression amplitude

Let M_x be the k-by-k binary matrix whose j-th column is the
output-clean X label of $`V_xZ_{A_j}V_x^\dagger`$, and let
$`r_x=\mathrm{rank}\,M_x`$. Expanding the output-zero projector gives

```math
K_x^\dagger K_x
 =2^{-k}\sum_{t\in\ker M_x^{\mathsf T}} R_t.             \qquad(3)
```

Here R_t is the signed data Pauli obtained by restricting
$`V_x^\dagger Z_A(t)V_x`$ to the input-clean zero sector. Terms outside
the displayed kernel have an input-clean X factor and vanish. The
surviving R_t form a commuting Pauli group, including their literal
signs. If this group contains negative identity, its sum vanishes. If
it has a nontrivial data image, its average is a proper stabilizer
projector. Both alternatives contradict invertibility of K_x, which
follows from (1). Thus every surviving R_t is positive identity, and

```math
K_x^\dagger K_x=2^{-r_x}I.                              \qquad(4)
```

Every singular value is consequently $`2^{-r_x/2}`$. Its distance from
$`1/\sqrt2`$ is at most epsilon. The closest different value is 1/2,
at distance $`(\sqrt2-1)/2`$, so (1) forces $`r_x=1`$ at every
address. This proves exact flatness; it does not merely bound average
success probability.

## Clean-stabilizer elimination determines the logical tableau

The following argument initially needs only invertibility of K_x.
Suppress x, and choose an input data Pauli generator P. Let b be the
output-clean X label of $`VPV^\dagger`$. Then

```math
Ms=b                                                        \qquad(5)
```

has a solution. Otherwise there is a vector lambda with
$`\lambda^{\mathsf T}M=0`$ and $`\lambda\cdot b=1`$.
The Pauli $`V^\dagger Z_A(\lambda)V`$ has no input-clean X factor,
so its restriction is a signed data Pauli R. Since the output-zero
bra is fixed by $`Z_A(\lambda)`$, one has KR=K. Invertibility gives
R=I. But $`\lambda\cdot b=1`$ says that R anticommutes with P,
a contradiction.

Dress P by the input-clean stabilizer $`Z_A(s)`$ from (5).
Its image has no output-clean X factor. Its data restriction R_P obeys
$`R_PK=KP`$, so the normalized K conjugates every data Pauli to a
Pauli and is Clifford. If $`s\in\ker M`$, the same argument applied
to $`Z_A(s)`$ gives identity data restriction, including positive sign.
The resulting logical Pauli is therefore independent of the chosen
solution of (5).

For any one output-data label coordinate, write this value as
$`d+c\cdot s`$. It is determined by the k-squared entries of M,
the k entries of b, the k entries of c, and d. Extend that function
arbitrarily on invalid tuples. Its degree in these $`(k+1)^2`$
arguments is at most $`(k+1)^2`$, proving the general bound
$`D(k+1)^2`$ after substitution.

Rank one sharpens the result. Write $`M=u v^{\mathsf T}`$ with
nonzero u,v. Solvability gives $`b=\beta u`$, and solution independence
gives $`c=\gamma v^{\mathsf T}`$, for binary beta,gamma. Hence

```math
d+c\cdot s=d+\mathrm{OR}(b)\mathrm{OR}(c).                \qquad(6)
```

Each OR has degree at most kD in the address, giving 2kD. The same
degree bound holds for the inverse logical symplectic matrix: symplectic
inversion is a fixed transpose and relabeling of entries.

## The logical sign, including the Pauli multiplication phase

Let the images of P and the j-th input-clean Z generator have labels
v,w_j and signs $`\sigma_P,\sigma_j`$. Define binary selectors

```math
\begin{aligned}
\rho_i&=b_i\prod_{h\lt i}(1+b_h),\\
\xi_{ij}&=M_{ij}\prod_{r\lt j}(1+M_{ir}),\\
\eta_j&=\sum_i\rho_i\xi_{ij}.
\end{aligned}                                             \qquad(7)
```

Their address degrees are at most 2kD. If b=0, all selectors vanish.
Otherwise they choose the first nonzero row indicated by b and the
first nonzero entry in that row. Since M has rank one and b is its
nonzero image vector, exactly one selected column equals b. Thus the
dressing can be chosen as a single clean Z generator.

For globally commuting Hermitian Paulis define chi by
$`P(v)P(w)=(-1)^{\chi(v,w)}P(v+w)`$. The logical sign is exactly

```math
\sigma_{\rm logical}
 =\sigma_P+\sum_j\eta_j\bigl(\sigma_j+\chi(v,w_j)\bigr).
                                                               \qquad(8)
```

Restriction to the output-clean zero sector adds no sign once its
X label vanishes.

For completeness, chi has degree at most 4D, independently of width.
Write each one-qubit multiplication phase as $`i^{f_l}`$. Its low
bit is the local symplectic parity, of degree at most 2D. Its high bit
is a Boolean function of four label bits, hence has degree at most 4D.
The high bit of the sum of all f_l is the sum of those high bits plus
the pairwise products of their low bits. The global low bit is zero
because the two Paulis commute. Therefore the resulting sign has degree
at most 4D. Substitution in (8) proves (2).

## Physical signs for the literal Gaussian query grammar

Consider q complete affine-address binary XOR shear queries, in either
direction, with arbitrary fixed Cliffords on all non-address wires
between them. Each query is a literal phase-free linear permutation
at fixed x. Its action on the non-Hermitian convention
$`R(a,b)=X^aZ^b`$ has no phase and an affine symplectic matrix.
Consequently every intermediate label has degree at most q.

A fixed Clifford acts in this convention by
$`R(\ell)\mapsto i^{\phi(\ell)}R(S\ell)`$, where phi is a
multilinear polynomial over $`\mathbb Z_4`$ with linear terms and
even quadratic terms. A Boolean polynomial f of degree d has a
modulo-four lift $`f=f_0+2f_1`$, with
$`\deg f_0\le d`$ and $`\deg f_1\le2d`$: expand the XOR of its
algebraic-normal-form monomials and retain only single and pair terms
modulo four. Substituting current labels into phi therefore gives
degree at most 2q. Summing the phases of all fixed interlayers does
not increase that bound.

Conversion of the output to Hermitian Paulis subtracts
$`\sum_j a_j(x)b_j(x)`$ from the raw exponent. A product of two
degree-q Boolean lifts has modulo-four degree at most 3q, since the
product of their two even parts vanishes modulo four. The final
exponent is pointwise even. Uniqueness of multilinear polynomials over
$`\mathbb Z_4`$ (or Boolean-cube inversion of their values) makes
every coefficient even, so division by two gives a Boolean physical
sign polynomial of degree at most 3q. Thus

```math
D\le q,\qquad E\le3q.                                   \qquad(9)
```

This is a statement about actual query words and their signed lifts.
An arbitrary family with affine symplectic matrices need not satisfy
the sign bound: an address-dependent Pauli can conceal any prescribed
Boolean sign without changing its symplectic matrix.

Read-only global Clifford interlayers need a small qualification.
They can contribute an affine-in-address Pauli displacement in addition
to a fixed bank Clifford. Its extra raw phase is even and has degree
at most the current label degree plus one. Therefore
$`E\le\max(3q,q+1)`$; (9) still holds for q at least one.
At q=0, for example, CZ from an address wire to a bank wire has sign
degree one. A uniform coarse statement may replace q by q+1.
Fixed affine classical relabelings of the address preserve degree.
Actual address Hadamards are outside this read-only grammar.

## Observable consequences and complete return

From (1), the normalized compressed Clifford satisfies
$`\|U_x-J_x\|\le\sqrt2\varepsilon`$. Its conjugation of any unit
Pauli differs from the target by at most $`2\sqrt2\varepsilon`$,
which is below $`\sqrt2`$. Distinct Pauli labels have distance at
least $`\sqrt2`$; opposite signs of the same Pauli have distance two.
Thus both labels and signs of a prescribed target Pauli conjugation
are fixed exactly.

For $`J_x^\dagger Z_{Y_0}J_x=Z_{Y_x}\otimes I_{\rm dirty}`$,
the addressed label has Boolean degree p. Applying the logical label
bound, in the required conjugation orientation, gives $`2kD\ge p`$.

For a direct complete indicator $`J_x=X_{Y_x}\otimes I_{\rm dirty}`$,
fix any output coordinate u. Its conjugation of $`Z_{Y_u}`$ has sign
$`(-1)^{[x=u]}`$, whose exponent has degree p. Equation (2) gives
$`E+(2k+4)D\ge p`$. The information is in the sign even though the
indicator has identity symplectic matrix. An address-dependent global
phase cannot hide this conjugation sign. Substitution of (9) gives

```math
\begin{aligned}
2kq&\ge p&&\text{for inverse-column label transport},\\
(2k+7)q&\ge p&&\text{for the direct indicator}.
\end{aligned}                                           \qquad(10)
```

No postselection is required to use an available compression. The
[retained general-k recovery, Section 6.3](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#63-recovering-a-general-clean-bank-compression-with-one-extra-flag)
already pays the zero-sector phase reflection and its extra clean flag.
Writing $`\Pi=I_B\otimes|0^k\rangle\langle0^k|`$ and
$`R_\Pi=I+(i-1)\Pi`$, the flatness (4) at rank one gives the exact
operator identity

```math
V_xR_\Pi V_x^\dagger R_\Pi V_x\Pi
 =e^{3\pi i/4}U_x\Pi.                                  \qquad(11)
```

The paid realization of these reflections returns its extra flag and
zero-test pool. Equation (11) returns the original k-qubit signal
exactly as well. It has full initialized-isometry error at most
$`\sqrt2\varepsilon`$ against the prescribed J_x, including dirty
return and arbitrary references, up to the displayed fixed removable
phase. Its inverse is the actual reversed word and has the same error
on the returned sector. No new recovery recurrence or free reflection
is asserted here.

This constant-amplitude Gaussian compression grammar requires
$`q=\Omega(p/(k+1))`$, hence $`q=\Omega(p/\log\log p)`$ when
$`k=O(\log\log p)`$. It does not convert q into a lower bound on
native T-depth. Joint native resynthesis, genuine non-Clifford carrier
dynamics, and address motion
that destroys the stipulated fixed-address Gaussian structure remain
outside the proof. In particular, the result supplies no new shallow
source and no unrestricted frontier improvement.
