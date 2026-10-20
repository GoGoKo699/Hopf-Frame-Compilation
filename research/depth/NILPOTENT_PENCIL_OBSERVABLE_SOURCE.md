# A nonhomomorphic affine-pencil source for the complete bank

[Complete observable contract](OBSERVABLE_TRANSPORT_SOURCE.md) ·
[Returned bilinear query](PARALLEL_BILINEAR_WORD_QUERY.md) ·
[Earlier inverse-map boundary](GAUSSIAN_INVERSE_EXTRACTION.md)

This is an explicit sufficient source and a paid quadratic XOR query,
not a faster in-place compiler. A fixed Vandermonde change of Pauli-label
coordinates puts all addressed output observables into one inverse
column of the affine nilpotent pencil below. Its address-dependent
family is not required to be a representation of XOR translations.
The unpaid operation is a complete shallow implementation of its
in-place map, or a different whole word transporting that one probe.

## 1. Coordinates and the exact inverse column

Let $`p\ge1`$, $`m=2^p`$, and $`F=\mathbb F_{2^p}`$. Fix a binary
basis $`\beta_0=1,\beta_1,\ldots,\beta_{p-1}`$ of F. The physical
p-bit address x denotes $`\sum_i x_i\beta_i`$. This is a fixed linear
identification, not a loader or an initialized register.

Work first in $`F^m`$, with standard coordinate vectors
$`e_0,\ldots,e_{m-1}`$. Put

```math
Le_j=e_{j+1}\ (j\lt m-1),\qquad Le_{m-1}=0,\qquad
\bar N_x=I+xL,\qquad
v_u=(1,u,u^2,\ldots,u^{m-1})^{\mathsf T}.
```

Here $`u^0=1`$ also at u=0. Nilpotence gives

```math
\bar N_x^{-1}=\sum_{j=0}^{m-1}x^jL^j,\qquad
\bar N_xv_x=e_0=v_0.                                      \qquad(1)
```

The m vectors $`v_u`$, with $`u\in F`$, form an F-basis. Their
Vandermonde determinant is the product of the nonzero differences of
the distinct field elements. Consequently

```math
\mathcal B=\{\beta_s v_u:u\in F,\ 0\le s\lt p\}
```

is a binary basis of the n=pm dimensional carrier. Assign the existing
output wire $`Y_u`$ to its s=0 coordinate. All remaining
$`(p-1)m`$ coordinates are arbitrary private inputs. In particular,
there is no initialized polynomial coefficient or supplied impulse.

Let E be the invertible binary matrix whose columns are the vectors of
$`\mathcal B`$ in the standard binary basis
$`\{\beta_s e_j\}`$. All matrices below act on binary column vectors.
Define

```math
N_x=E^{-1}\bar N_x E.
```

These are matrices on **Pauli labels** in the physical carrier
coordinates. If $`a_u`$ is the standard binary label of $`Z_{Y_u}`$,
then (1) says exactly

```math
N_x^{-1}a_0=a_x.                                         \qquad(2)
```

The field notation in (1) and the binary conjugation in (2) must not
be interchanged with a data-register transpose.

## 2. A sufficient source and its full error recurrence

For an invertible binary matrix A, use the convention
$`U_A|z\rangle=|Az\rangle`$. With
$`Z(v)|z\rangle=(-1)^{v^{\mathsf T}z}|z\rangle`$,

```math
U_A^\dagger Z(v)U_A=Z(A^{\mathsf T}v).
```

Thus the direct source

```math
J=\sum_x|x\rangle\langle x|\otimes U_{N_x^{-\mathsf T}}
```

obeys the literal complete identity

```math
J^\dagger Z_{Y_0}J
=\sum_x|x\rangle\langle x|\otimes Z_{Y_x}
 \otimes I_{(p-1)m}.                                    \qquad(3)
```

The target is the whole addressed bank. It is not one scalar inverse
entry evaluated on an initialized carrier.

It suffices to implement the opposite orientation
$`W=J^\dagger`$, whose principal coordinate map is

```math
(x,z)\longmapsto(x,N_x^{\mathsf T}z).                    \qquad(4)
```

Every coordinate of (4) has binary degree at most two. A still weaker
sufficient ideal allows arbitrary coherent private-helper fibers:

```math
W_0=\sum_{x,z}|x,N_x^{\mathsf T}z\rangle\langle x,z|
                   \otimes R_{x,z},                     \qquad(5)
```

where each $`R_{x,z}`$ is a unitary on all additional private work.
This is a complete-space condition; no helper starts in a preferred
state. Set $`J=W_0^\dagger`$. Conjugation of the carrier-only probe
by W cancels the two matching helper fibers, so (3) remains true
tensored with identity on that work. Intermediate address motion and
helper nonreturn are permitted in an implementation of (5).

Suppose an emitted V is within eta in full operator norm of a specified
$`W_0`$ in (5), and its forward and actual reverse have common ledger
$`T_V,G_V,d_V,h_V`$. Use $`J=V^\dagger`$. Unitary telescoping gives

```math
\|V Z_{Y_0}V^\dagger-\Omega_p\|\le2\eta.
```

The chronological complete indicator word is

```math
H_Y,\quad V^\dagger,\quad Z_{Y_0},\quad V,\quad H_Y.
```

Its simultaneous common ledger is

```math
\begin{aligned}
T_I&\le2T_V,&D_T(I)&\le2d_V,\\
G_I&\le2G_V+2m+2,&h_I&\le(p-1)m+h_V,\\
\varepsilon_I&\le2\eta,&a_{\rm clean}&=0.
\end{aligned}                                           \qquad(6)
```

All return and reference errors are included. If only the transported
observable is approximated with error delta, replace $`2\eta`$ by
delta in (6); a full-source approximation is sufficient, not necessary.
The sixteen-source-call fusion recurrence is then the retained one:
error at most $`16\eta`$ and depth
$`16d_V+O_\tau(\log\log(p+4))`$, at the smaller fused address size.
One may request $`\eta=\tau/16`$. Polynomial source overhead is
normalized only if the same V also has sublogarithmic depth.

## 3. The full quadratic XOR query is paid

Consider two arbitrary n-bit principal banks z,t. The available query is

```math
Q_N:(x,z,t,d)\longmapsto
       (x,z,t+N_x^{\mathsf T}z,d).                      \qquad(7)
```

It preserves x,z and returns every additional dirty helper. It is not
the in-place map (4).

For explicit common bounds, put

```math
n=mp,\qquad r=(m-1)p^2,\qquad s\le(m-1)p^3,
\qquad c_n=n^2+3n.
```

An invertible binary n-by-n matrix and its inverse can each be compiled
using at most $`c_n`$ CNOTs by Gaussian elimination, with each row swap
expanded into three CNOTs. These are fixed matrices, independent of x.

First apply $`E^{-\mathsf T}`$ to both principal banks. In the resulting
standard dual coordinates, add z into t using n CNOTs. The remaining
bilinear term is

```math
t_j\longmapsto t_j+M_x^{\mathsf T}z_{j+1},
\qquad0\le j\lt m-1,                                 \qquad(8)
```

where each block has p bits and $`M_x`$ is the binary matrix for
multiplication by x in F. Its entries are linear in the address bits.
Finally undo both fixed bank changes with $`E^{\mathsf T}`$.

For every output bit in (8) and every address bit, allocate one arbitrary
R91 helper. The associated linear form in $`z_{j+1}`$ has at most p
terms. The returned R91 phase derivative applies all these terms in
eight T layers: physical T targets are distinct private helpers,
while shared address, source and phase-output wires occur only as
Clifford controls. The total number of such helpers is r, and the total
linear-form incidence is s. Enclose the phase word by Hadamards on the
$`(m-1)p`$ affected target bits.

The two four-T relative-phase words use sixteen parity CNOTs per helper;
the linear toggles and their actual inverses use at most $`2s`$ CNOTs.
Thus (7), including its fixed changes of basis and linear term, has

```math
\begin{aligned}
T_{Q_N}&\le8r,&D_T(Q_N)&\le8,\\
G_{Q_N}&\le4c_n+3n+16r+2s,&h_{Q_N}&\le r,\\
\varepsilon_{Q_N}&=0,&a_{\rm clean}&=0.
\end{aligned}                                           \qquad(9)
```

The actual reversed native word has the same bounds. The extra n-bit
target t is a principal input of (7), not part of $`h_{Q_N}`$. When
used as an added source-building bank, its n wires must be allocated
in addition to the $`(p-1)m`$ private carrier wires and r query helpers.
Consequently one such embedded component needs
$`(2p-1)m+r`$ private wires outside the specified address and output Y.
All of these are arbitrary inputs. No field-arithmetic clean scratch,
zero bank, or private phase depending on its initial value is omitted.

The count orders in (9) are
$`T,h=O(mp^2)`$, $`G=O(m^2p^2+mp^3)`$, and constant depth.
They fit the relaxed source envelopes. They pay the XOR query only.

## 4. Why an ordinary completion retains logarithmic depth

The inverse of (4) has the probe coordinate

```math
(N_x^{-\mathsf T}z)_{Y_0}=z_{Y_x}.                      \qquad(10)
```

As a Boolean polynomial in the p address bits and the arbitrary bank,
the right side has degree exactly p+1. The forward map (4) is quadratic,
but its inverse is not uniformly bounded-degree. Therefore one must
not infer a cheap inverse from the degree of the forward formula.

Suppose an actual proposed J consists of r monomial reversible rounds
of coordinate degree at most k, interspersed with arbitrary affine
changes of all coordinates. Phases are unrestricted. Its transported
Z probe is diagonal with signs determined by a coordinate polynomial
of degree at most $`k^r`$. Any observable error below two forces those
signs to agree exactly with $`(-1)^{Y_x}`$ on every input. Hence

```math
k^r\ge p+1.                                            \qquad(11)
```

If applying (11) by reversing a circuit for W, the inverse of each
declared round must also have degree at most k. This holds for affine
rounds and the paid XOR-query involutions; it does not hold for an
arbitrary quadratic permutation, as (4) itself demonstrates. No
monomial-round bound is asserted for general native Clifford+T words.

## 5. The Vandermonde pencil is a star with one zero leaf pivot

There is a useful exact description in F-coordinates. Let B have the
columns $`v_u`$, and put $`A=B^{-1}LB`$. For $`p\ge2`$,

```math
A_{00}=0,\qquad
A_{0u}=A_{u0}=A_{uu}=u^{-1}\quad(u\ne0),
```

and all remaining entries vanish. Indeed
$`Lv_u=u^{-1}(v_u+v_0)`$ for nonzero u, while the finite-field power
sums give $`Lv_0=e_1=\sum_{u\ne0}u^{-1}v_u`$. Thus the same pencil
acts on F-coordinate labels as

```math
\begin{aligned}
(N_xz)_0&=z_0+\sum_{u\ne0}(x/u)z_u,\\
(N_xz)_u&=(1+x/u)z_u+(x/u)z_0\quad(u\ne0).
\end{aligned}                                           \qquad(12)
```

For nonzero x, exactly the leaf u=x has zero diagonal coefficient;
its output is $`z_0`$. At x=0 the pencil is identity. These are
Pauli-label coordinates: the physical data map remains the binary
transpose prescribed in Section 2. Symmetry over F does not silently
identify a chosen binary multiplication matrix with its transpose.
For p=1, the separate two-by-two matrix A has all four entries one;
the source identity (1)--(3) is unchanged.

Equation (12) is a concrete zero-pivot target for a whole quantum
readout. Ordinary leaf elimination introduces inverses of
$`1+x/u`$ and needs separate treatment of the selected zero pivot.
Those reciprocals or a coherent choice of the pivot have not been
implemented. Likewise, an address Hadamard around a phase linear in x
and bilinear in two banks is just a quadratic address-feedback shear;
grouping such classical completed shears retains (11). A useful next
word must keep genuine quantum interference across those cuts.

## 6. Two-bank Fourier dilation still owes the inverse column

At each fixed x, (7) is a Clifford shear. Its binary symplectic entries
are affine in x. Both query directions, their transposes, arbitrary
fixed Clifford changes of either bank, and arbitrary dirty spectators
are allowed in the following observation. A word containing q such
affine-symplectic query stages and fixed Clifford boundaries has
symplectic entries of Boolean degree at most q. Thus its transported
Pauli label has degree at most q.

The desired label $`a_x`$ in (2) has degree p: each physical output
coordinate is the equality function for one address. Distinct Pauli
labels have operator-norm distance at least $`\sqrt2`$, including
arbitrary literal phases. Therefore this complete read-only-address
Clifford-fiber dilation needs at least p affine query stages at
observable error below $`\sqrt2`$. This is the retained
[affine-symplectic round argument](GAUSSIAN_INVERSE_EXTRACTION.md#4-a-tight-affine-symplectic-round-obstruction)
applied to the full inverse column, not a new unrestricted depth bound.
It allows Fourier/key changes on the arbitrary banks; no common key is
required. It does not cover actual address motion, non-Clifford fiber
dynamics, or joint native resynthesis that removes these stage cuts.

For example, the paid Choi/Bell change $`C=Q_NH_z`$ only gives

```math
C^\dagger Z_t(b)C=X_z(N_xb)Z_t(b),
\qquad
CZ_z(a)C^\dagger=X_z(a)X_t(N_x^{\mathsf T}a).
```

Neither identity produces $`N_x^{-1}a_0`$ or removes the arbitrary
second-bank factor. Replacing that bank by zero changes the contract.

This differs from the earlier Gaussian scalar-entry proposal in what
is encoded: (2) is every output label in a fixed full-bank embedding,
with no table loader and no initialized impulse. It shares the same
unpaid distinction between an XOR shear and an in-place map. The first
open operation remains a complete nonmonomial transport of the single
inverse column with the qualifying native ledger. No frontier
improvement follows from these identities alone.

### A paired completion is equivalent to the weak in-place source

For clarity, put $`A_x=N_x^{\mathsf T}`$. A supplied weak source (5)
gives the complete inverse XOR query with the chronological word

```math
W_0^\dagger(t,H),\quad\mathrm{COPY}(t\to z),\quad W_0(t,H).
```

It adds $`A_x^{-1}t`$ into arbitrary z. COPY leaves t unchanged, so
the two private fibers cancel literally. Therefore the chronological
word

```math
Q_A(z;t),\quad Q_{A^{-1}}(t;z),\quad Q_A(z;t),\quad
\mathrm{SWAP}(z,t)
```

is exactly
$`(z,t,H)\mapsto(A_xz,A_x^{-1}t,H)`$. With a source V at error eta,
and the exact query (9), this paired completion has

```math
\begin{aligned}
T_{\rm pair}&\le2T_V+2T_{Q_N},&
d_{\rm pair}&\le2d_V+2d_{Q_N},\\
G_{\rm pair}&\le2G_V+2G_{Q_N}+4n,&
h_{\rm pool}&\le\max(h_V,h_{Q_N}),\\
\varepsilon_{\rm pair}&\le2\eta.
\end{aligned}                                           \qquad(13)
```

The extra principal bank t costs n arbitrary wires in addition to the
carrier and the displayed pool. Conversely, this exact paired map is
itself a weak source (5), with t as a helper fiber. No controlled source
call is needed. This equivalence does not remove the inverse operation
or make a word satisfying the single-probe promise itself a paired
completion. Section 8 gives a separate paid conversion from that promise.

## 7. Bounded evidence and next operation

The [GF(4) matrix fixture](../../tests/test_nilpotent_pencil_observable_source.py)
checks the binary Vandermonde basis, all four inverse columns, the
physical transpose and output aliases, the star formula, and the XOR
query's complete block matrix and actual involution. Matrix identities
cover every arbitrary bank value; no initialized-bank state test or
full Hilbert-space enumeration is used. The fixture does not emit the
unpaid in-place source or prove an asymptotic schedule.

The selected next question is a whole quantum readout of the star's
zero pivot, with arbitrary carrier and second-bank inputs, outside the
fixed-address Clifford-stage and classical-feedback cuts above. The
first unpaid operation is the complete inverse-column transport, not
the finite-field basis, its quadratic XOR query, actual reversal, or
the wrapper recurrence. The equivalent interfaces and their fixtures
do not improve the unrestricted frontier.

## 8. A paid converse through arbitrary-bank convolution

The weak in-place completion is also obtainable from a supplied complete
indicator with only constant-factor depth loss and polynomial overhead.
This is a resource equivalence in the relaxed envelopes, not a cheaper
mechanism for obtaining either object. Keep $`A_x=N_x^{\mathsf T}`$ and
let a supplied native indicator have complete error tau and common
forward/actual-inverse bounds $`T_I,G_I,d_I,h_I`$. All its helper inputs
are arbitrary and its ideal helper action is identity.

### A returned bilinear convolution query

For an arbitrary n-bit coefficient bank $`C=(C_0,\ldots,C_{m-1})\in F^m`$,
define the binary linear map $`P_C`$ in standard field coordinates by

```math
(P_Cy)_j=\sum_{d=0}^j C_dy_{j-d},\qquad 0\le j\lt m.
```

On two other arbitrary n-bit principal banks z,t, pay the exact query

```math
B:(C,z,t,H)\longmapsto(C,z,t+P_C^{\mathsf T}z,H).          \qquad(14)
```

Its target-Hadamard phase is $`(-1)^{z^{\mathsf T}P_Ct}`$;
the transpose is binary, with no implicit trace-pairing basis change.
To obtain a trilinear decomposition, embed F in
$`F'=\mathbb F_{2^{2p}}`$, evaluate the two degree-at-most-m-1
polynomials at $`2m-1`$ distinct fixed elements of F', multiply their
values, and interpolate. Such points exist because
$`m^2\ge2m-1`$. Take the first m product coefficients and their fixed
binary pairing with z. Schoolbook binary multiplication in F' uses
$`(2p)^2`$ scalar products. Thus the phase is a sum of at most

```math
R=(2m-1)(2p)^2\le8mp^2
```

products of three fixed binary linear forms, respectively in C,t,z.
Evaluation and interpolation here specify these forms; no evaluated
register is initialized or supplied for free.

Apply the [returned cubic-mask construction](ALGEBRA_QUERY_EXTRACTION.md#4-fast-subset-convolution-as-a-native-returned-query)
to three separate arbitrary R-bit masks. Each form shear can be emitted
directly with at most nR CNOTs. The eight mask vertices
$`000,100,110,010,011,111,101,001`$, followed by return to 000, use eight
shears and eight layers of R disjoint literal CCZ words. The
[seven-T, ten-CNOT, five-layer CCZ word](RETURNED_DIAGONAL_ECHO_DEGREE.md#sharp-cubic-example-and-a-literal-native-ledger)
therefore gives the simultaneous ledger

```math
\begin{aligned}
T_B&\le56R,&d_B&\le40,\\
G_B&\le8nR+80R+2n,&h_B&\le3R,\\
\varepsilon_B&=0,&a_{\rm clean}&=0.
\end{aligned}                                           \qquad(15)
```

The 2n term is the target Hadamard envelope. All three masks return
exactly on every input; C,z,t are principal banks excluded from
$`h_B`$. Actual reversal has the same native bounds. In particular
$`T_B,h_B=O(mp^2)`$ and $`G_B=O(m^2p^3)`$.

### The coefficient echo pays the inverse query

On C, the chronological fixed change $`E^{-1}`$, supplied indicator
on its s=0 coordinates, and change E gives a word $`\widetilde D_x`$
within tau of the complete translation
$`D_x:C\mapsto C+v_x`$. The other $`(p-1)m`$ coefficient coordinates
are arbitrary, not initialized. This word has T-count $`T_I`$,
depth $`d_I`$, Clifford count at most $`G_I+2c_n`$, and private
helper width $`h_I`$.

In standard dual coordinates, execute chronologically

```math
B,\quad\widetilde D_x,\quad B^\dagger,
\quad\widetilde D_x^\dagger.                             \qquad(16)
```

For the exact indicator, C follows $`C,C+v_x,C+v_x,C`$ and the
two target additions sum to

```math
(P_C^{\mathsf T}+P_{C+v_x}^{\mathsf T})z
 =P_{v_x}^{\mathsf T}z=\bar N_x^{-\mathsf T}z.
```

The coefficient bank and every mask/helper are restored. With the
emitted indicator, use its actual inverse in the last position;
unitary telescoping bounds the complete query and return error by
$`2\tau`$. No control of the indicator, clean impulse or replacement
inverse is assumed.

Wrap (16) by $`E^{-\mathsf T}`$ on both z,t before it and
$`E^{\mathsf T}`$ on both afterwards. The resulting ideal query is
$`Q_{A^{-1}}:(z,t)\mapsto(z,t+A_x^{-1}z)`$, because
$`A_x^{-1}=E^{\mathsf T}\bar N_x^{-\mathsf T}E^{-\mathsf T}`$.
Its common forward/actual-inverse ledger is

```math
\begin{aligned}
T_{\rm inv}&\le2T_I+2T_B,&d_{\rm inv}&\le2d_I+80,\\
G_{\rm inv}&\le2G_I+2G_B+8c_n,&
h_{\rm inv}&\le n+\max(h_I,3R),\\
\varepsilon_{\rm inv}&\le2\tau,&a_{\rm clean}&=0.
\end{aligned}                                           \qquad(17)
```

Here both n-bit query banks are principal. The additional n bits in
$`h_{\rm inv}`$ are C. The eight fixed bank changes consist of four
on the principal banks and four in the two coefficient translations;
their actual inverses are included. The mask and indicator-helper
pools can be reused because their ideal actions return the entire
pool, and the approximation is on its complete input space.

### Paired source and precise contract distinction

Insert (17) in the three-query and SWAP completion of Section 6.
The resulting source approximates
$`(z,t)\mapsto(A_xz,A_x^{-1}t)`$ with common bounds

```math
\begin{aligned}
T_{\rm src}&\le2T_I+2T_B+2T_{Q_N},&d_{\rm src}&\le2d_I+96,\\
G_{\rm src}&\le2G_I+2G_B+2G_{Q_N}+8c_n+3n,&
h_{\rm src}&\le2n+\max(h_I,3R,h_{Q_N}),\\
\varepsilon_{\rm src}&\le2\tau,&a_{\rm clean}&=0.
\end{aligned}                                           \qquad(18)
```

The width in (18) excludes the principal carrier z but includes the
second bank t, coefficient bank C, and the common helper pool. Thus
the total private width outside address and the specified output Y
is $`(p-1)m+h_{\rm src}`$. The second bank has the explicitly allowed
fiber action $`U_{A_x^{-1}}`$ in (5); C and the pool have ideal identity
action. Every error above includes their return and arbitrary external
references. Choosing $`\tau=\eta/2`$ gives source error eta. Together
with (6), this proves equivalence of the complete indicator and weak
in-place/paired source under the relaxed polynomial count/width and
quadratic-times-polynomial Clifford envelopes, with constant-factor
depth changes. The additive 96 is a paid depth bound, not a new shallow
source construction.

The literal single-probe promise (3) is nevertheless weaker than (5).
For example $`J_x=\mathrm{SWAP}(Y_0,Y_x)\otimes I_{\rm private}`$
satisfies (3). For $`p\ge2`$, its x=1 principal action has order two,
whereas $`N_1`$ and $`A_1`$ have order $`m\ge4`$:
$`(I+L)^{m/2}=I+L^{m/2}\ne I`$ and $`(I+L)^m=I`$.
Consequently this J's inverse is not itself a source of the form (5),
even with arbitrary helper fibers. This distinguishes contracts, not
their computational difficulty under paid conversions: the observable
wrapper followed by (18) is such a conversion and retains the supplied
depth. Neither the weak contract nor this converse improves the full
resource frontier by itself.

## 9. Spectral processing alone does not discriminate the zero pivot

For one star leaf, put $`a=1+x/u\in F`$, and consider the exact
field query and its phase form on two arbitrary p-bit banks:

```math
Q_a|z,t\rangle=|z,t+az\rangle,\qquad
D_a|z,t\rangle=(-1)^{\mathrm{Tr}_{F/\mathbb F_2}(azt)}
                  |z,t\rangle.
```

The phase form is the target-bank Fourier conjugate of the XOR query
in fixed trace-dual coordinates. Both are returned bilinear-query
components of the type paid in Section 3. This observation does not
grant a free controlled native implementation.

Consider an actual word whose only a-dependence is through Q_a, D_a,
their actual adjoints, or signal-controlled occurrences of these
queries. Between occurrences allow arbitrary unitaries on a separate
signal system, including logical targets, clean signals and dirty
spectators, but no change of the two query banks or of a. Signal
operations may be non-Clifford, and their number and width are
unrestricted in this test.

The common bank subspace

```math
\mathcal K=\mathrm{span}\{|0,t\rangle:t\in F\}
```

has eigenvalue +1 for every Q_a and D_a, for every a. Their actual
adjoints have the same property. Every controlled occurrence is
therefore identity on the entire signal system tensored with
$`\mathcal K`$. If the intervening signal words are
$`V_0,\ldots,V_s`$, the full word restricts literally to
$`(V_s\cdots V_0)\otimes I_{\mathcal K}`$, independently of a.
This retains coherent bank superpositions and references; it is not
a trace or average-input statement.

Consequently this grammar cannot implement a helper-independent
zero-pivot discriminator. More quantitatively, if the required signal
unitaries at a=0 and at some nonzero a are U_0 and U_1, both with
identity on the banks, a common complete error bound epsilon requires
$`2\epsilon\ge\|U_0-U_1\|_{\rm op}`$. A conditional toggle or
sign with distance two between those branches thus has error at least
one. More clean signals, precision or repeated spectral processing
do not remove this common-eigenspace defect.

This statement concerns using the leaf query as a spectral primitive
for a separate returned discriminator. It is not a lower bound for
direct transport of a probe inside the live banks, for bank-mixing
interference, for intermediate address motion, or for a whole word
whose a-dependence enters elsewhere. Such an operation must supply its
own complete observable and return identity.

## 10. A paid address-motion trial retains translation symmetry

There is a literal nonmonomial trial using the paid query and address
rotations. For this trial, assign the specified output Y to the s=0
coordinates of the query's target bank t, with the source bank z
entirely private. This makes the transported probe nontrivial; keeping
Y in z instead would leave its Z probe unchanged. For any fixed signs
$`\sigma_i\in\{-1,1\}`$, put

```math
R_X=\prod_{i=1}^p H_{X_i}T_{X_i}^{\sigma_i}H_{X_i},
\qquad J=Q_N^\dagger R_XQ_N.
```

The chronological word is the emitted query, R_X, then that query's
actual inverse. Its ideal query is an involution, so this is also the
ideal word $`Q_NR_XQ_N`$. R_X moves the address coherently and is
nonmonomial. Both J and its actual inverse have the common exact ledger

```math
\begin{aligned}
T_J&\le2T_{Q_N}+p,&d_J&\le2d_{Q_N}+1\le17,\\
G_J&\le2G_{Q_N}+2p,&
h_J&\le(2p-1)m+h_{Q_N},\\
\varepsilon_J&=0,&a_{\rm clean}&=0.
\end{aligned}                                           \qquad(19)
```

Here h_J counts every private wire outside the address and specified
output Y: the $`(p-1)m`$ remaining target coordinates, the source bank
of n=pm bits, and the reused query pool. This is a wire assignment for
the trial, so no SWAP conversion is needed. The helper pool is returned
exactly, and R_X uses only the existing address wires.

Nevertheless this complete word has the wrong symmetry. Write
$`\Delta_iN=N_{x+e_i}+N_x`$, which is independent of x because the
pencil is affine. Let $`Q_{\Delta_iN}`$ be the fixed bank shear
$`(z,t)\mapsto(z,t+(\Delta_iN)^{\mathsf T}z)`$. On the complete
space, including the returned query pool,

```math
Q_N^\dagger X_{X_i}Q_N=X_{X_i}Q_{\Delta_iN}.              \qquad(20)
```

Every right-hand side commutes with every address translation
$`X_{X_j}`$. Up to its irrelevant scalar phase, R_X is a product of
rotations about $`X_{X_i}`$, so (20) proves
$`[J,X_{X_j}]=0`$ for every j. The transported bank probe
$`O=J^\dagger Z_{Y_0}J`$ therefore has the same symmetry. In contrast,

```math
\|X_{X_j}\Omega_pX_{X_j}-\Omega_p\|_{\rm op}=2,
\qquad \|O-\Omega_p\|_{\rm op}\ge1.                     \qquad(21)
```

The equality follows by fixing an address x and assigning opposite Z
eigenvalues to the distinct output coordinates x and x+e_j. The second
claim is the triangle inequality, since
$`X_{X_j}OX_{X_j}=O`$. Thus even though (19) pays an actual coherent
address-motion word at constant depth and the required count/width
scales, it does not meet the observable contract at error below one.

More completed conjugated-rotation blocks and address-independent bank
interlayers preserve this symmetry. An unmatched query or a
noncommuting computational-address phase can break it, but this alone
does not supply a zero-pivot readout. The selected next attempt must
also supply genuine nonmonomial carrier mixing across the query stages
and be analyzed as a complete word. This bounded failed trial does not
exclude broader address-motion or non-Clifford bank interference.
