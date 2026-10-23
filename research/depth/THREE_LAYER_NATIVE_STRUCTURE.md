# Three native T layers: a degree bound and a paid cubic query

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Two-layer diagonal rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md) ·
[Native axis matching](ONE_BAND_AXIS_MATCHING.md) ·
[Session outcome](THREE_LAYER_NATIVE_SESSION_RESULT.md)

This proof concerns unrestricted physical Clifford+T circuits. Clifford
blocks may mix every address, output and dirty wire. No intermediate
address-preservation or third-level-source premise is imposed. The exact
degree theorem below is **zero clean**; it does not cover a circuit that
is specified only on a clean-zero sector.

## 1. Exact statement and conventions

Let V be diagonal on its entire q-qubit physical space and have native
T-depth at most three. In its original computational coordinates, its
diagonal entries have the form

```math
V|x\rangle=\omega^{p(x)}|x\rangle,\qquad
\omega=e^{i\pi/4},
```

where p is an integer multilinear polynomial modulo eight of ordinary
degree at most five. An irrelevant common native scalar can instead be
kept outside this expression. The assertion is an ordinary degree bound
in any invertible affine computational coordinates, not a diagonal-hierarchy
classification or a weighted-degree-five assertion.

Consequently, if V is projectively a Boolean sign operation
$`(-1)^{g(x)}`$, then $`\deg_{\mathbb F_2}g\le5`$.
If a three-layer circuit implements a literal Boolean XOR query
$`Q_f:|x,z,w\rangle\mapsto|x,z+f(x),w\rangle`$ with every arbitrary
dirty wire w returned, then

```math
\deg_{\mathbb F_2}f\le4.                              \qquad(1)
```

Hadamards on z turn the query into the sign $`(-1)^{zf(x)}`$ without
changing T-depth. Its Boolean degree is $`\deg f+1`$ for nonzero f.
The conclusion permits arbitrary dirty width and references, since it
is a full operator identity. It is a necessary condition, not a
converse synthesis theorem.

Use signed Hermitian Pauli axes P, binary symplectic labels p, and the
literal native convention

```math
\begin{gathered}
T_P=\frac{I+P}{2}+\omega\frac{I-P}{2},\\
S_P=T_P^2,\qquad T_{-P}=\omega T_P^\dagger.          \qquad(2)
\end{gathered}
```

A native band is a product of T or T-dagger on independent commuting
Pauli axes, in fixed Clifford coordinates. Independence is always
physical independence on the full space. Every fixed Clifford and
scalar introduced in the proof is retained in the operator identities;
the proof rearranges an existing circuit and is not a free synthesis
claim for any rearranged factor.

## 2. One-band versus two-band encoders

Split a three-layer circuit, absorbing arbitrary Clifford blocks into
axis coordinates and endpoint Cliffords, as

```math
V=K_C^\dagger E_C^\dagger E_B E_A K_A.              \qquad(3)
```

The E factors are native bands; empty bands are allowed. Let
$`R_z=K_AZ(z)K_A^\dagger`$ and
$`R'_z=K_CZ(z)K_C^\dagger`$, where z ranges over all binary q-vectors.
Their binary spans $`M_A,M_C`$ are maximal isotropic spaces. Diagonality
is exactly the simultaneous family of identities

```math
J_z=E_AR_zE_A^\dagger,\qquad
S_z=E_BJ_zE_B^\dagger=E_CR'_zE_C^\dagger.             \qquad(4)
```

Both J_z and S_z are Hermitian Cliffords, and each family is the exact
commuting representation of the input Z group. The
[one-band matching theorem](ONE_BAND_AXIS_MATCHING.md) therefore says
that J_z permutes the individual binary axes of B, up to signs. This
uses the fact that both ends of this particular one-band passage are
Clifford; it does not apply matching through an unknown non-Clifford
intermediate.

Write a_i for the A-axis labels and
$`\lambda_i(z)=\langle a_i,R_z\rangle`$. The symplectic action of J_z is

```math
F_z=I+N_z,\qquad
N_zv=\sum_i\lambda_i(z)\langle a_i,v\rangle a_i.
                                                               \qquad(5)
```

Every lambda_i is linear in z. Since the A axes commute,
$`N_zN_u=0`$ for all z,u. Globally silent A axes have
$`\lambda_i\equiv0`$ and hence lie in $`M_A`$. Their band factors
peel to the input as diagonal C3 operations:

```math
E_A=E_{Aa}E_{As},\qquad
E_AK_A=E_{Aa}K_A(K_A^\dagger E_{As}K_A).
```

The peeled factor commutes the input Z group, so removing it preserves
diagonality. It does not change J_z. From now on A denotes only the
globally active axes, and this peeled diagonal factor is restored at
the end.

## 3. Middle-band orbits are singletons or pairs

The permutation representation induced by F_z on the independent
B-axis set has no orbit of size four or larger. An elementary abelian
group acts transitively on each of its orbits. If an orbit had at least
four points, choose z,u whose translations are independent in that
orbit. Then

```math
N_zN_u b=b+F_zb+F_ub+F_zF_ub
```

would be a sum of four distinct B axes. Independence makes it nonzero,
contradicting $`N_zN_u=0`$.

For a paired orbit $`\{b,d\}`$, put $`r=b+d`$. There is a nonzero
linear swap character lambda such that

```math
N_zb=N_zd=\lambda(z)r.                              \qquad(6)
```

Equation (5) and A-axis independence show that every A axis having
nonzero pairing with b has activation character exactly lambda, and
r is the sum of those axes. Since b and d are B axes, they commute,
so $`\langle b,r\rangle=0`$. Thus the number of these A neighbors is
even. It follows that the activation character of r is zero, hence

```math
r\in M_A\cap\mathrm{span}(A).                 \qquad(7)
```

The differences r of distinct pairs are independent: each is the sum
of its own two distinct members of the independent B list. A singleton
B axis has $`N_zb=0`$ for every z, so it commutes every remaining A
axis; the globally silent A axes were already peeled.

The signed action on each paired orbit has the exact form

```math
\begin{gathered}
J_z B J_z^\dagger=(-1)^{\chi(z)}B_{\lambda(z)},\\
J_z D J_z^\dagger=(-1)^{\chi(z)}B_{1+\lambda(z)},    \qquad(8)
\end{gathered}
```

where $`B_0=B,B_1=D`$ and chi is linear. To verify this, a signed
swap that squares to identity has equal off-diagonal signs. Any
nonswap commuting with that swap must have equal diagonal signs.
The exact group law then makes both chi and lambda characters.
Singleton signs likewise give linear characters.

## 4. Native matching restricts the paired signs

First orient every B rotation as T rather than T-dagger. This is an
exact rearrangement: if $`E_B=F\widetilde E_B`$, where F is the
product of the necessary $`S_B^\dagger`$ factors, then

```math
K_C^\dagger E_C^\dagger F\widetilde E_B
=(F^\dagger K_C)^\dagger
 (F^\dagger E_CF)^\dagger\widetilde E_B.             \qquad(9)
```

F is Clifford; conjugating E_C by F preserves a native independent
axis list. Thus this normalization changes neither the number of
bands nor J_z, and it keeps all literal phases.

For the positive B orientation, the exact defect
$`E_BJ_zE_B^\dagger J_z^\dagger`$ is, up to its specified scalar,
the product of S rotations on B axes whose signs in (8) are negative.
Its symplectic action is $`I+N_{B,z}`$, with

```math
N_{B,z}v=\sum_b\chi_b(z)\langle b,v\rangle b.
```

Consequently

```math
F_{S_z}=I+N_z+N_{B,z}+N_{B,z}N_z.                   \qquad(10)
```

All terms except the last are linear functions of z. The right side
of (4), being a native one-band Pauli image, also has symplectic
action linear in z. A singleton contributes zero to the last term.
The pair $`\{b,b+r\}`$ contributes exactly

```math
\chi(z)\lambda(z)\,r\langle r,\cdot\rangle.        \qquad(11)
```

The operators $`r\langle r,\cdot\rangle`$ for the different pairs
are linearly independent: use symplectic dual vectors to their
independent r labels. Therefore every product chi lambda is linear.
Since lambda is nonzero and both are linear Boolean functions, this
forces

```math
\chi=0\quad\hbox{or}\quad\chi=\lambda.             \qquad(12)
```

Choose equal T orientations on a pair with chi=0 and opposite
orientations on a pair with chi=lambda. Further orientation changes
are absorbed by (9), so the opposite encoder remains native. Each
entire pair phase now commutes every J_z **literally**. For chi=0,
swapping B and D preserves $`T_BT_D`$. For chi=lambda, the only
nontrivial action is $`(B,D)\mapsto(-D,-B)`$, and (2) gives

```math
T_{-D}T_{-B}^\dagger
=(\omega T_D^\dagger)(\omega^{-1}T_B)=T_BT_D^\dagger.
```

The native scalars cancel in this last equality; no input-dependent
phase has been dropped.

## 5. Peeling the pairs and closing the native remainder

Write $`E_B=E_{Bf}\prod_j P_j`$, where f denotes the singleton
axes and P_j the orientation-adjusted pair phases. All these factors
commute. The operators

```math
F_j=K_A^\dagger E_A^\dagger P_j E_AK_A              \qquad(13)
```

are diagonal, because P_j commutes every J_z. Hence

```math
\begin{gathered}
V=V_0\prod_jF_j,\\
V_0=K_C^\dagger E_C^\dagger U K_A,\\
U=E_{Bf}E_A,                                       \qquad(14)
\end{gathered}
```

with the earlier peeled input-diagonal factor understood. Each
singleton B axis commutes every active A axis, so U is a C3 phase
polynomial supported on one isotropic Pauli space L. The entire
V_0 is diagonal.

The following closing argument is needed; diagonality alone does not
justify simply assuming that a two-band source is C3. For each input
Z, let $`R=K_A Z K_A^\dagger`$ and
$`R'=K_C Z K_C^\dagger`$. Then

```math
S=URU^\dagger=E_C R'E_C^\dagger.                    \qquad(15)
```

S is Clifford because U is C3. For any Pauli Q with label in
$`L^\perp`$, U commutes Q, and therefore $`SQS^\dagger=\pm Q`$.
The symplectic action F_S fixes $`L^\perp`$ pointwise. Symplecticity
then implies

```math
\mathrm{im}(F_S-I)\subseteq L.                 \qquad(16)
```

For the native C-band image on the right of (15), this image is
exactly the span of the C axes anticommuting with R'. Indeed its
symplectic displacement is the sum of the independent transvections
on those axes; physical independence makes its image their full
span. Varying Z shows that every globally active C axis lies in L.

Split $`E_C=E_{Ca}E_{Cs}`$ into globally active and silent axes.
The silent axes lie in $`M_C`$, and
$`K_C^\dagger E_{Cs}^\dagger K_C`$ is diagonal C3. Peel that factor
on the output side of (14). The remaining middle operation
$`E_{Ca}^\dagger U`$ is a C3 phase polynomial, since all its axes
belong to the same isotropic L. Multiplying it by the endpoint
Cliffords preserves C3. The remaining operation is also diagonal,
because both V_0 and the peeled factor are diagonal. Products of
diagonal C3 operations are diagonal C3, so **V_0 is diagonal C3**.
This proves the required closure without a missing source premise.

<a id="6-exact-pair-phases-and-the-common-coordinate-bound"></a>

## 6. Exact pair phases and the original-coordinate degree bound

For a paired orbit choose its signed Hermitian difference to be
$`R=BD`$; the commuting B,D make this literal equality unambiguous.
Its label r is in $`M_A\cap\mathrm{span}(A)`$, so E_A
commutes R and $`K_A^\dagger R K_A`$ is a signed input Z parity.
Let u be its eigenvalue bit. Also put

```math
H=K_A^\dagger E_A^\dagger B E_AK_A.
```

H is a Hermitian Clifford, commutes this input parity, and the
conjugated D is exactly HR. Formula (2) gives the full branchwise
pair phases:

| Pair orientations | u=0 | u=1 |
| --- | --- | --- |
| $`T_BT_D`$ | $`S_H`$ | $`\omega I`$ |
| $`T_BT_D^\dagger`$ | $`I`$ | $`\omega^{-1}S_H`$ |

On the branch containing S_H, diagonality of F_j implies diagonality
of H: $`S_H=((1+i)/2)I+((1-i)/2)H`$ has nonzero H coefficient.
The restriction of H to that parity sector is a logical Clifford.
To see this directly, make the parity one coordinate Z_1. Since H
commutes Z_1, it preserves both sectors; it conjugates any logical
Pauli to a physical Pauli commuting Z_1, which restricts to a logical
Pauli. Thus its diagonal Hermitian eigenvalue is
$`(-1)^{h(y)}`$ for a Boolean polynomial h of degree at most two in
the other coordinates. Extend h to be independent of u. The literal
phase exponents of the two rows of the table are respectively

```math
p_j=u+2(1-u)h(y),\qquad
p_j=-u+2u h(y)\pmod8.                               \qquad(17)
```

All pair-difference labels r_j are independent. Therefore one common
invertible affine change of the input coordinates makes every signed
parity bit u_j a distinct coordinate. For each j the restriction of H
above is still Clifford in these same coordinates, so h_j is quadratic
in the remaining coordinates. This is a convenient way to construct
the h_j. Returning to the original coordinates keeps each u_j affine
Boolean and each h_j quadratic Boolean. The particular divisibilities
in (17) preserve the degree-five bound under this substitution, although
an arbitrary integer degree-five polynomial need not retain its degree.

For a Boolean quadratic h, write its algebraic normal form as an XOR
of monomials of degree at most two, including a possible constant.
Its integer lift satisfies

```math
h=\sum_\alpha m_\alpha
 -2\sum_{\alpha\lt\beta}m_\alpha m_\beta\pmod4.
```

More explicitly, for an affine Boolean v and a quadratic Boolean h in
any affine coordinates, their integer lifts have the forms

```math
v=L_1-2L_2\pmod4,\qquad h=H_2-2H_4\pmod4,
```

where each subscript bounds the ordinary degree of its integer
polynomial. Constants are allowed. Therefore

```math
2vh=2L_1H_2-4(L_2H_2+L_1H_4)\pmod8
```

has degree at most five. The apparent degree-six product has coefficient
eight and vanishes. An affine Boolean u itself has integer degree at
most three modulo eight. Applying these facts to v=u or v=1-u in (17)
proves degree at most five in the original coordinates, and in every
invertible affine coordinate system. The residual V_0 and the peeled
silent-A phase are diagonal C3; their standard phase polynomials have
linear coefficients modulo eight, quadratic coefficients divisible
by two, cubic coefficients divisible by four, and no higher terms.
This weighted cubic property is preserved by affine coordinate changes.
Adding these exponents and all (17) proves the claimed degree-five
bound.

If the final diagonal is a common scalar times $`(-1)^g`$, subtract
its value at the zero string. Its exponent is $`4(g(x)+g(0))`$
modulo eight. Multilinear coefficients over $`\mathbb Z/8\mathbb Z`$
are unique by Boolean Möbius inversion. Hence every coefficient of
g above degree five vanishes over $`\mathbb F_2`$ in the original
coordinates, proving (1).

## 7. Fixed-error zero-clean consequence

Let an at-most-three-layer circuit V approximate a Boolean diagonal
sign G on the entire physical space, allowing one freely optimized
common phase. Put

```math
\tau=\frac{\sin(\pi/8)}{2\sqrt2}\approx0.135299.
```

If $`\|V-e^{i\phi}G\|\lt\tau`$, the retained
[R75 operator-rigidity theorem](FOUR_LAYER_OPERATOR_RIGIDITY.md)
forces V to preserve each individual input Z exactly. It applies to
at most four layers, so three is included. This conclusion holds
simultaneously without a union bound, and makes V diagonal.

Every diagonal entry is a modulus-one element of
$`\mathbb Z[i,1/\sqrt2]`$, hence an eighth root by the retained
[unit-circle lemma](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#42-exact-dirty-return-and-the-explicit-zero-clean-floor).
If the ratios of these entries to G are not constant, two distinct
eighth roots occur. Their smallest possible maximum distance to a
common phase is $`2\sin(\pi/16)\gt\tau`$. Thus closeness below tau
forces $`V=\gamma G`$ exactly for one common scalar gamma.
Section 6 then gives $`\deg G\le5`$.

Therefore a Boolean XOR query with $`\deg f\gt4`$ has full operator
error at least tau at zero-clean T-depth at most three, independently
of dirty width and allowing a free common phase. For integer precision
exponents $`L\ge3`$, such a query needs at least four T layers.
The same conclusion applies with arbitrary external references and
requires complete dirty return. It is a constant-layer obstruction,
not a dimension-growing lower bound.

## 8. Paying the third band: an exact cubic XOR query

There is an exact three-layer AND3 query with one clean qubit and one
arbitrary dirty helper. This completes the phase correction that could
not fit the preceding two-layer construction. Its depth is optimal at
this allocation by the [two-layer lower theorem](TWO_LAYER_DIAGONAL_RIGIDITY.md).
No optimal T-count or elementary Clifford count is asserted.

### 8.1. Literal phase cancellation

Use arbitrary bits $`t,x_1,x_2,x_3,x_4`$ and clean c=0. Let e_j be
the elementary symmetric integer polynomial of degree j in the four
x bits. Put

```math
q=e_1+e_2\pmod2,\qquad
2q=2e_1-2e_2+4e_4\pmod8.                            \qquad(18)
```

The nine-CNOT circuit C applies, chronologically, t to each x_i,
then each x_i to c, then t to c. It computes the six physical parities
$`t,t+x_1,\ldots,t+x_4,c+t+\sum_i x_i`$ over the binary field.
Define

```math
\begin{gathered}
E=C^\dagger T^{\otimes6}C,\\
A=H_tS_t^\dagger H_t,\qquad K=S_tA,\\
V=K^\dagger E^\dagger H_tEH_t.                       \qquad(19)
\end{gathered}
```

This is the [literal two-layer construction](TWO_LAYER_DIAGONAL_RIGIDITY.md#5-exact-one-clean-diagonal-non-c3-counterexample).
For completeness, let p_0 be the integer sum of the six bits
$`(0,x_1,x_2,x_3,x_4,\bigoplus_i x_i)`$. It is even, and on c=0

```math
E=\omega^{p_0}S_t^\dagger Z_t^q,\qquad
H_tEH_t=\omega^{p_0}A X_t^q.
```

The x-dependent scalar cancels in (19). Since
$`A^\dagger Z_tA=-Y_t`$, the occupied action is exactly

```math
V_0=(iZ_t)^q,\qquad p_V=2q+4tq\pmod8.              \qquad(20)
```

Define M by three chronological CNOT groups: x_j to x_4 for j=1,2,3;
x_4 to x_j for j=1,2,3; then x_j to x_4 again. This nine-CNOT
involution replaces every x_i by the parity of the other three.
Set $`C'=CM`$, so its chronological word is M followed by C, and put

```math
\begin{gathered}
E'=C'^\dagger T^{\otimes6}C',\\
D=S_t\prod_{i=1}^4 S_{x_i}\prod_{i=1}^4\mathrm{CZ}(t,x_i),\\
R=DE'.                                      \qquad(21)
\end{gathered}
```

The six band rows are t, the four t-plus-three-x parities, and
$`c+t+\sum_i x_i`$. They are physically independent because C' is
an invertible CNOT circuit; the last row's c cannot be omitted when
counting physical layers. On c=0 the sum of the five nonzero x parities
is $`4e_1+2e_2`$ modulo eight. Including t gives E' exponent
$`4e_1+2e_2+6t+4te_2`$. Adding D's exponent yields

```math
p_R=-2e_1+2e_2+4tq\pmod8.                          \qquad(22)
```

Equations (18), (20) and (22) cancel both the imaginary phase and the
unwanted t-dependent sign:

```math
\begin{gathered}
(RV)_0=(-1)^{x_1x_2x_3x_4}\otimes I_t,\\
Q=H_{x_4}RVH_{x_4},\\
Q_0|x,t\rangle=|x_1,x_2,x_3,x_4+x_1x_2x_3,t\rangle. \qquad(23)
\end{gathered}
```

Every occupied query column has literal phase one. The correction is
an actual paid native band; it is not supplied as an abstract C3 gate.

### 8.2. Complete ledger and inverse

Each band in (19) and (21) acts once on each of the six wires. C and
its actual inverse each cost nine CNOTs; C' and its inverse each cost
eighteen. Expand each CZ as H-CNOT-H and S-dagger as three S gates.
The resulting exact ledgers are:

| Word | T count | T-depth | Elementary H/S/CNOT |
| --- | ---: | ---: | ---: |
| V | 12 | 2 | 44 |
| R | 6 | 1 | 53 |
| RV | 18 | 3 | 97 |
| Q | 18 | 3 | 99 |
| Actual raw inverse of Q | 18 | 3 | 117 |
| Actual inverse after adjacent-S reduction | 18 | 3 | 109 |

The forward query has 76 CNOT, 14 H and nine S gates. Reversing the
actual word exchanges T with T-dagger and expands every inverse S.
Nine consecutive S gates from the old S-cubed segment reduce to one,
saving eight gates in the inverse. The occupied query is an involution,
but its unused c=1 block is not Hermitian: the forward word must not
be substituted for the full actual inverse.

The simultaneous allocation is four logical bits, one arbitrary dirty
helper t and one clean helper c, with no additional work. Both V and R
preserve the physical c sectors. Equation (23) proves zero leakage,
exact clean and dirty return, and error zero at every precision.
Tensoring that occupied operator identity with an external identity
covers arbitrary reference entanglement. The actual inverse has the
same occupied Boolean action and the separately charged ledger above.

The [exact native fixtures](../../tests/test_three_layer_native_mechanisms.py)
check all 32 occupied columns with literal phases, the three independent
physical bands and all 64 columns of the actual inverse. The analytic
identities, rather than numerical sampling, establish the construction.

### 8.3. A fairly simplified baseline

The exact relative Toffoli
$`J_{ab,c}=(-i)^{ab}\mathrm{CCX}_{ab,c}`$ uses four T gates,
two layers and eight elementary Cliffords. Between H_c boundaries,
its chronological diagonal word is

```text
CX(c,b), T(c), T-dagger(b), CX(c,b),
CX(c,a), CX(a,b), T(b), T-dagger(a), CX(a,b), CX(c,a).
```

Its phase identity, with parenthesized sums denoting XOR, is

```math
c-(b+c)+(a+b+c)-(a+c)=-2ab+4abc\pmod8.
```

Compute J into clean c from x_1,x_2, apply an exact dirty-helper
Toffoli from c,x_3 into x_4 using t, and apply the actual J-dagger.
The relative phase cancels because the middle gate preserves both
controls and c. The middle gate uses the already simplified
[eight-T, two-layer, eighteen-Clifford parity word](DIAGONAL_ENCODING_SESSION_RESULT.md#2-constructive-check-and-first-unpaid-cost).
The complete baseline and its actual inverse each have

```math
T=16,\qquad G=34,\qquad D_T\le6,
```

with the same one-clean/one-dirty allocation and zero error. The exact
fixtures verify its whole occupied action and actual inverse. The new
construction improves the demonstrated depth while increasing both
gate counts; six is only the displayed baseline's depth upper bound.
The two-layer lower theorem, not a claim about this baseline's optimal
synthesis, proves that the new depth three is optimal.

### 8.4. Growing arity does not remove the remaining address phase

The construction extends to k x bits for every positive k divisible
by four. C costs $`2k+1`$ CNOTs and the same complementary-parity
map M costs $`3(k-1)`$, with x_k as pivot. Thus C' costs $`5k-2`$.
For $`q=e_1+e_2\pmod2`$, (18) holds at every k. One proof expands
the XOR of the linear and quadratic monomials modulo four: the
pair-products have parity $`e_2+e_4`$, since each triple occurs six
times and each quadruple three times. The complementary parity sum is

```math
\sum_i\left(\bigoplus_j x_j+x_i\right)
 +\bigoplus_jx_j=ke_1+2e_2\pmod8.                  \qquad(24)
```

Here additions inside parentheses are binary. Expanding a j-variable
parity gives coefficient $`(-2)^{j-1}`$; each j-subset occurs
$`k-j+1`$ times, proving (24) when k is divisible by four.

If k is four modulo eight use positive T in E; if k is zero modulo
eight use T-dagger. In both cases its t exponent is six modulo eight,
so (20) is unchanged. Keep the correction E' positive, and replace
each S in D by $`S^r`$, where r=1 in the first case and r=3 in the
second. Equation (24) again gives (22). The full occupied operation is
$`(-1)^{e_4(x)}\otimes I_t`$, with

```math
T=3(k+2),\qquad D_T=3,\qquad
G=21k+8+r(k+1),
```

and actual inverse Clifford count
$`21k+8+(4-r)(k+1)`$ after adjacent-S reduction. There is one clean
c and one dirty t, no other work, zero error and complete return as
above. This is a growing-arity family of degree four, not growing
Boolean degree. An output Hadamard on x_k gives

```math
(-1)^{e_4(x_1,\ldots,x_{k-1})}
X_{x_k}^{e_3(x_1,\ldots,x_{k-1})}.
```

For k=4 the address phase is one. For larger k it remains a nontrivial
paid-operation requirement; dropping it would change the full operator.
This family supplies no faster arbitrary-table query for R98.

## 9. What fails when one clean qubit is allowed

The zero-clean proof does not automatically apply to an occupied
clean-zero sector. The [one-relation matching theorem](ONE_BAND_AXIS_MATCHING.md#6-a-projected-band-with-one-relation)
now recovers the projected axes of a band that restricts to a Pauli
stabilizer code, except for one possible three-dimensional ambiguity.
This is a bounded component result. A general non-Pauli intermediate
code need not admit that quotient at all.

Nor does the previous two-layer clean-removal construction extend to
an arbitrary one-band/two-band split. Here is a fully paid example on
clean c and arbitrary d. Put

```math
\begin{gathered}
K=\mathrm{CX}_{c,d}H_c,\\
C=\mathrm{CX}_{c,d}H_c\mathrm{CX}_{c,d},\\
E=T_cT_d,\qquad B=C(T_cT_d^\dagger)C^\dagger.
\end{gathered}
```

Let $`U_A=EK`$, $`U_B=B^\dagger EK`$ and
$`W=U_B^\dagger U_A=K^\dagger E^\dagger BEK`$. Both encoders
have the same non-Pauli clean image
$`J=(X_c+Y_c)(X_d+Y_d)/2`$. The first encoder's only binary Pauli
axes that both flip J and have Pauli preimages are Z_c and Z_d.
Neither has a Pauli preimage under U_B. For example, at the outer
band its image is a sum of four Paulis; the inner E preserves the
nonzero antisymmetric XY/YX part. The exact fixture exhausts all
sixteen binary Pauli labels, including their irrelevant overall signs,
and verifies that there is no common clean-flipping Pauli of this type.
Thus the shared Pauli required by the two-layer construction cannot be
inferred for this split.

This displayed W uses six T gates, ten elementary Cliffords and three
layers; its actual inverse has the same costs. Its fair complete
simplification is merely

```math
W=\mathrm{CX}_{c,d}(T_cT_d^\dagger)\mathrm{CX}_{c,d},
```

which uses two T gates, two Cliffords and one layer. Its c=0 action
is just $`T_d^\dagger`$. The example disproves a representation-specific
inference; it is no lower bound on the target, no difficult query and
no evidence that another split or native representation cannot remove
initialization.

The remaining lower question is a degree bound independent of arity
and dirty width for **one-clean three-layer literal queries**, or a
complete native counterfamily. A proof must control the entire logical
Pauli algebra on a possibly non-Pauli intermediate code, without assuming
the failed common flip or a source-wide C3 property. The constructions
above must remain admissible. This is a test toward a depth-dependent
invariant; neither the constant-layer bound nor the finite optimal
query changes the complete-frame resource frontier.

