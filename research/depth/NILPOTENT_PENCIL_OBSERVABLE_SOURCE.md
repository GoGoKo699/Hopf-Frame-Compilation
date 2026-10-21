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
The [Gaussian compression revision](GAUSSIAN_COMPRESSION_SOURCE_LIMITS.md)
now bounds both logical labels and signs for the stated small-clean,
fixed-address Clifford grammar. Its affine-query stage bound does not
lower-bound unrestricted native depth or exclude non-Gaussian readout.

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

Within the query architecture, changing only one bank's basis has an
additional exact limitation. Every query, source-bank Fourier change
and address-only gate commutes with all target translations X_t(v).
Any word built from them retains the fixed target-translation character
of its probe, whereas the addressed observable changes that character
with x. The same orthogonal-eigenspace argument gives a gap of at least
$`\sqrt2`$. Changing only the target bank instead leaves the arbitrary
source sector z=0 invariant, where the query is identity and no bank
selection is performed. These statements test these actual primitives;
they do not cover general native circuits. A new query-based attempt
must mix both carrier banks across stages and prove complete readout,
not merely add address phases or a one-bank Fourier boundary.

<a id="latent-center-readout"></a>

## 11. A paid forward subgroup has the addressed observable in its center

Retain the physical coordinates of Sections 1–3, with n=pm, and enumerate
the binary basis so that e_0=a_0. To distinguish the inverse label from
the field vectors v_u of Section 1, write

```math
\begin{aligned}
 A_x&=N_x^{\mathsf T},\qquad \nu_x=N_x^{-1}a_0=a_x,\\
 w_j(x)&=A_xe_j\quad(j\ne0),\\
 \mathcal H_x&=\mathrm{span}\{w_j(x):j\ne0\}
             =\mathrm{ker}(\nu_x^{\mathsf T}).
\end{aligned}                                                       \qquad(22)
```

The last equality follows from invertibility and
$`\nu_x^{\mathsf T}A_x=e_0^{\mathsf T}`$. The specified output wires
and the other (p−1)m carrier wires are all arbitrary. On this complete
carrier use the literal real Pauli subgroup

```math
 \Gamma_x=\langle X(w_j(x))\ (j\ne0),\ Z(e_i)\ (0\le i\lt n)\rangle,
 \qquad X(a)|z\rangle=|z+a\rangle.                       \qquad(23)
```

These are already paid native generators. Since A_x is affine in x and
A_0=I, each w_j(x) is e_j plus p fixed binary masks controlled by the
address bits. Emit X on j and the corresponding address-to-carrier
CNOTs. Every completed boundary fixes x; there is no helper or clean
wire. Safe common forward/actual-reverse bounds are T=d=0 and
G≤np+8 for an X generator, and G≤6 for a Z generator. The constants
cover unfused inverse S words in the positive-S alphabet. A literal
concatenation using each of the 2n−1 generators once therefore has

```math
 T=d=0,\qquad G\le(n-1)(np+8)+6n=O(m^2p^3).             \qquad(24)
```

Its only width beyond x and specified Y is the existing (p−1)m carrier
coordinates. This is a complete family, not an extraction claim for
that concatenation. The whole coherent translation by A_xy is instead
the Section 3 query Q_N on arbitrary coefficient bank y and target z,
with the separate n-bit principal bank and ledger (9). The subgroup
corresponds algebraically to y_0=0; this does not initialize y, erase its
first coordinate, or grant a prepared coherent selection of generators.

The center is exactly the desired label. Commutation with every Z(e_i)
forces a Pauli's X label to vanish. Commutation with all translations
in H_x then forces its Z label into span(nu_x). All Z labels are
present, so for $`\Omega_x=Z(\nu_x)=Z_{Y_x}`$,

```math
 Z(\Gamma_x)=\{I,-I,\Omega_x,-\Omega_x\}.                \qquad(25)
```

On each of the two Omega_x sectors, translations by H_x and all
characters generate the full matrix algebra on that affine coset.
Thus the two sector representations are irreducible, each of dimension
$`2^{n-1}`$. Their characters differ only on scalar multiples of
Omega_x; every noncentral element has trace zero on both sectors.
The subgroup and this claim are fiberwise. No common abstract expander
or compiler for its proposed generators follows merely from the fibers.

Ordinary invariant mixing does not select the central character. Since
n≥2, some nonzero w in H_x anticommutes with a fixed-coordinate Z.
Their commutator is −I, so neither sector has an invariant vector.
A one-dimensional character sends every commutator to +1; twisting
cannot change this conclusion. More explicitly, the available Weyl
commutator is $`(-1)^{b^{\mathsf T}A_xe_j}I`$: it reads a forward-entry
scalar, not the bank-dependent inverse-column observable. An ordinary
group expander therefore approaches the zero invariant projector on
both sectors, even if its complete native implementation is granted.

Conjugation has a different invariant space. For normalized trace
$`\tau_n(O)=2^{-n}\mathrm{tr}(O)`$, Pauli orthogonality gives the exact
twirl on the complete carrier:

```math
\begin{aligned}
 \Phi_x(O)&=|\Gamma_x|^{-1}\sum_{g\in\Gamma_x}gOg^\dagger\\
           &=\tau_n(O)I+\tau_n(\Omega_xO)\Omega_x.
\end{aligned}                                                       \qquad(26)
```

For a fixed starting probe O independent of x with $`\|O\|\le1`$,
the m addressed Paulis are orthogonal. Hence

```math
 \sum_x|\tau_n(\Omega_xO)|^2
 \le\tau_n(O^\dagger O)\le1.                            \qquad(27)
```

A uniform approximation of Phi_x(O) to alpha Omega_x thus has error at
least $`|\alpha|-2^{-p/2}`$. This is a fixed-probe twirl statement;
a separately supplied x-dependent probe requires its own paid complete
word. Vectorizing conjugation gives only
$`\mathrm{span}\{\mathrm{vec}(I),\mathrm{vec}(\Omega_x)\}`$,
a rank-two invariant space in dimension $`4^n`$. This does not prepare
a code on the two arbitrary banks. Spectators and references tensor
these identities with the appropriate untouched identity action.

The first unpaid operation is full-space central-character readout with
constant amplitude and the required native return and sharing contract.
Preserving Omega merely preserves its sectors; using an invariant
projector equal to $`(I+\Omega_x)/2`$ would require every generator to
fix that entire sector pointwise. The paid subgroup instead acts
nontrivially on both. Its abelianization has binary rank 2n−1, so the
Pauli subgroup itself has no constant-size generating set. A proposed
ambient automorphism expander still needs its actual generators and
selected representation; the canonical conjugation A_x S A_x^{-1}
exposes the unpaid inverse. A coherent SELECT over the listed basis
generators likewise needs a priced index and preparation. None of
these observations rules out a different full central-character
protocol or supplies a qualifying shallow source.

<a id="paired-projector-readout"></a>

## 12. A half-rank relative word and its paid one-signal readout

The paired representation has a useful target that is not a thin
invariant-vector projector. It also makes the unpaid normalization
explicit. Let C and D be two arbitrary n-bit carriers, put
$`d=2^{n-1}`$, and let S swap their entire states. For binary labels
a,b define the Hermitian Pauli
$`P(a,b)=i^{\sum_j a_jb_j}X(a)Z(b)`$. Then

```math
\begin{aligned}
 \Pi_x&=\frac{I+\Omega_{x,C}\Omega_{x,D}}2,\\
 B_x&=\frac1{2d}\sum_{a\in H_x,\ b\in\mathbb F_2^n}
           P(a,b)\otimes P(a,b)=\Pi_x S,\\
 B_x^\dagger B_x&=B_xB_x^\dagger=\Pi_x.
\end{aligned}                                                       \qquad(28)
```

To prove (28), change binary basis so that nu_x is the first Z label.
The subgroup's Hermitian Paulis are I or Z on that coordinate and
arbitrary Paulis on the remaining n−1 coordinates. The ordinary
Pauli expansion of swap on those coordinates is
$`d^{-1}\sum_P P\otimes P`$. On the equal-central-bit sector the
first-coordinate swap is identity. Conjugating back on both carriers
preserves their full swap; Pauli signs cancel between the two copies.
Thus B_x has norm one and rank $`2^{2n-1}`$. The same formulas hold
coherently in x and tensored with every private reference identity.

There are $`2^{2n-1}=2d^2`$ summands in (28). Their coefficient
one-norm is d, not one. Ordinary uniform coherent SELECT therefore
compresses to $`B_x/d`$, not B_x. It also has a 2n−1 bit index whose
preparation and return must be paid. Writing a=A_xy with y_0=0 lets the
forward pencil supply labels; it does not pay this normalization or
initialize the coefficient bank. The full-space identity escapes the
thin-code rank argument, but by itself supplies no native advantage:
after a fixed swap its desired projector is already the lookup parity
on two banks.

### 12.1 Six actual source calls read a partially flat block

Here is a paid conditional extraction at the original amplitude
$`1/\sqrt2`$. Let U have exactly one clean signal a, with J inserting
a=0 and P=JJ-dagger. Every other register is arbitrary. In the next
formulas B can be any square partial isometry with initial projector
$`\Pi=B^\dagger B`$, including the coherent direct sum of B_x above.
Suppose its entire accepted block satisfies

```math
 K=J^\dagger UJ,\qquad
 \|K-B/\sqrt2\|\le\eta\le1/8,\qquad
 e(\eta)=3\sqrt2\eta+2\eta^2.
                                                                    \qquad(29)
```

No promise on U's other blocks is required. Use the literal words

```math
 W=US_aU^\dagger S_aU,\qquad V=W^\dagger Z_aW.
                                                                    \qquad(30)
```

The actual source adjoints, not separately approximated inverse targets,
occur in (30). For $`\delta=-1-i`$, direct block expansion gives

```math
\begin{aligned}
 WJ&=-\delta JK+UJ(iI-2iK^\dagger K),\\
 J^\dagger WJ&=[(1+2i)I-2iKK^\dagger]K.
\end{aligned}                                                       \qquad(31)
```

At eta=0 the first identity maps Pi into the clean sector as
$`e^{i\pi/4}JB`$, while the second vanishes on Pi-perp. Consequently
V has exact initialized action $`J(2\Pi-I)`$, with exact signal and
private return. Its accepted phase is retained in (31) and cancels in
the literal conjugation (30); no phase correction is assumed free.

The approximate bound remains linear, rather than using a generic
square-root leakage estimate. Since
$`\|K^\dagger K-\Pi/2\|\le\sqrt2\eta+\eta^2`$,

```math
 \|WJ\Pi-e^{i\pi/4}JB\|\le e(\eta),\qquad
 \|PWJ(I-\Pi)\|\le\sqrt5\eta.
                                                                    \qquad(32)
```

For the second bound, K is a contraction, the prefactor in the second
line of (31) has norm at most sqrt(5), and
$`\|K(I-\Pi)\|\le\eta`$. Put Q=W-dagger P W. Projection
geometry gives the exact relation

```math
 \|(Q-J\Pi J^\dagger)J\|
 =\max\{\|(I-P)WJ\Pi\|,\|PWJ(I-\Pi)\|\}
 \le e(\eta).
```

Indeed, the squared operator on the left is block diagonal in Pi and
Pi-perp, with blocks $`\Pi(I-J^\dagger QJ)\Pi`$ and
$`(I-\Pi)J^\dagger QJ(I-\Pi)`$. Therefore the complete return bound is

```math
 \|VJ-J(2\Pi-I)\|\le2e(\eta).                         \qquad(33)
```

There is a weaker exact sufficient interface: **only**
$`K^\dagger K=\Pi/2`$ is needed for exact readout. The polar partial
isometry $`B=\sqrt2 K`$ then has arbitrary final subspace and arbitrary
logical/private fibers; (31) still proves (33) with zero error. This
does not infer block closeness (29) from approximate Gram closeness.
The stated linear approximate bound uses (29) itself.

An approximate Gram promise is also sufficient, with its own weaker
error estimate. If A=K-dagger K and
$`\|A-\Pi/2\|\le\delta`$, the second line of (31) gives

```math
 \bigl(J^\dagger WJ\bigr)^\dagger\bigl(J^\dagger WJ\bigr)
   =g(A),\qquad g(t)=5t-8t^2+4t^3.
```

Since A is a positive contraction, $`g(\Pi/2)=\Pi`$ and
telescoping the three powers yields

```math
 \|g(A)-\Pi\|\le[5+8(1+1/2)+4(1+1/2+1/4)]\delta=24\delta.
```

The same projection identity proves

```math
 \|VJ-J(2\Pi-I)\|\le2\sqrt{24\delta}.                 \qquad(33a)
```

This source interface does not require a chosen accepted partial
isometry, any phase alignment, or any promise on its final subspace.

### 12.2 Complete indicator, resources and arbitrary-work return

For (28), $`2\Pi_x-I=\Omega_{x,C}\Omega_{x,D}`$. Bitwise CNOTs
from D into C conjugate this product to Omega_{x,C}. Hadamards on the
m specified Y coordinates of C convert it to the complete XOR
indicator. The other n−m coordinates of C and all n coordinates of D
are private arbitrary work, not initialized code states.

Let T_U,G_U,d_U,b_U be common bounds for the actual U and its actual
adjoint; b_U excludes x, both n-bit carriers and the clean signal.
The one-clean complete indicator V_I just constructed has

```math
\begin{aligned}
 T_{V_I}&\le6T_U,&
 G_{V_I}&\le6G_U+2n+2m+10,\\
 d_{V_I}&\le6d_U,&
 b_{V_I}&\le2n-m+b_U,\\
 \|V_IJ-JI_p\|&\le2e(\eta).
\end{aligned}                                                       \qquad(34)
```

The ten extra Cliffords inside V count two S, two actual S-daggers
emitted as three positive S each, and Z emitted as SS. For the
indicator's actual inverse emit **the same full word V_I**: (30) is
exactly Hermitian for every actual unitary U, even an approximate
source, and the outer Hadamard/CNOT conjugation preserves this fact.
This is a literal adjoint implementation with the same phase and count,
not a separately approximated inverse or a quotient by global phase.

Apply the [serial one-clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#3-covering-all-sectors-and-paying-the-complete-recurrence)
with q=m. Its full recurrence becomes

```math
\begin{aligned}
 T_{\rm out}&\le48T_U+48m+8,\\
 G_{\rm out}&\le48G_U+16n+98m+152,\\
 D_T({\rm out})&\le48d_U+32,\\
 b_{\rm out}&\le2n+b_U+2,\qquad b_{\rm clean}=0,\\
 \|\widetilde I_p-I_p\otimes I_{\rm borrowed}\|
   &\le32e(\eta)=96\sqrt2\eta+64\eta^2.
\end{aligned}                                                       \qquad(35)
```

For p≥1, m is even and the final eight T gates can be omitted. For
$`0\lt\tau\le1`$, request $`\eta=\tau/144`$ in the source and
evaluate every source resource at that tolerance. All source-signal,
second-carrier, private-work and promotion-toggle return errors are
already included in (35), on arbitrary coherent address, output, dirty
and reference inputs. The serial promotion has **no** shared-address
native-boundary premise: U may move the address internally and have
arbitrary off-code action. The stronger native-sharing hypotheses of
the separately parallel spectral wrapper are not imported here.

For the Gram-only source promise (33a), all resource lines of (35) are
unchanged and its full-operator error is instead at most
$`32\sqrt{24\delta}`$. Request
$`\delta=\tau^2/24576`$ and evaluate every source resource at that
Gram tolerance. This is not the linear block-error claim in (35).

The extraction uses one source-clean bit and removes it. An exact
projector-Gram source gives exact final return and action. A source
satisfying (29) uniformly at each fixed requested tolerance, with
T/width $`O(m\mathrm{poly}(p))`$, Clifford count
$`O(m^2\mathrm{poly}(p))`$ and depth $`o(\log p)`$, would
enter the retained polynomial-overhead fusion window. No such source
is emitted by the paired-Pauli sum or the forward pencil. Its
correct normalized lookup sector remains the first unpaid operation;
the [new normalized trials](CROSSED_CUBIC_KERNEL_SOURCE.md#normalized-forward-gram-source)
pay half-projector normalization for other sectors. The later
[signal-interleaved words](SIGNAL_INTERLEAVED_KERNEL_SOURCE.md) pay their
controlled phases but retain wrong-sector acceptance. The
[scalar-kernel source](SCALAR_KERNEL_FILTER_SOURCE.md) rejects the wrong
charge but misses most of the correct sector. The
[corank-one Fourier/Bell construction](CORANK_FOURIER_BELL_SOURCE.md)
identifies the exact desired partial isometry; shared outcomes and the
paid coherent-SWAP repair do not implement its complete normalized Gram.
Neither these scoped failures, (28), nor the ordinary SELECT loss is an
unrestricted native depth bound.

The [small full-column fixture](../../tests/test_nilpotent_pencil_observable_source.py)
checks the paired-Pauli identity and the six-call phase-S readout,
including a nontrivial private Hadamard fiber, every rejected signal
input, the literal phase, and the actual inverse. It is a finite
identity check, not a compiler for the missing asymptotic source.

<a id="four-call-gram-readout"></a>

### 12.3 Four actual calls suffice under the Gram-only promise

The Gram-only interface has a shorter recovery than (30). This does
not improve the linear-in-eta block-error estimate (33), and does not
supply a source. Retain one clean signal, let $`Z=Z_a`$, and define
the actual native word

```math
 M=ZU^\dagger ZU,\qquad F=M^2,\qquad
 G=K^\dagger K,\quad K=J^\dagger UJ.
                                                               \qquad(35a)
```

Its chronological order is U,Z,U-dagger,Z,U,Z,U-dagger,Z. It uses
four source/actual-adjoint calls and four literal Z gates. In the clean
and rejected block decomposition, $`U^\dagger ZU`$ is a Hermitian
involution with upper-left block $`E=2G-I`$. Its off-diagonal block
therefore has squared norm operator $`I-E^2`$. Multiplication gives
the exact complete accepted block

```math
 J^\dagger FJ=2E^2-I=8G^2-8G+I=:f(G).                 \qquad(35b)
```

Suppose $`\|G-\Pi/2\|\le\delta\le1`$ for a projector Pi.
No commutation of G and Pi is assumed. Put
$`S_0=I-2\Pi`$, $`\Delta=G-\Pi/2`$, and
$`D=FJ-JS_0`$. Since F and S_0 are unitary and (35b) is Hermitian,

```math
\begin{aligned}
 D^\dagger D
 &=2I-S_0f(G)-f(G)S_0\\
 &=16\Pi\Delta^2\Pi
   +16(I-\Pi)G(I-G)(I-\Pi).                         \qquad(35c)
\end{aligned}
```

The two summands have orthogonal Pi and Pi-perp supports. The first
has norm at most $`16\delta^2`$. Because G is a positive contraction,
the second is bounded above by $`16(I-\Pi)G(I-\Pi)`$, whose norm
is at most $`16\delta`$. Thus the full initialized-isometry estimate is

```math
 \|FJ+J(2\Pi-I)\|\le4\sqrt\delta.                    \qquad(35d)
```

This includes rejected-signal leakage, every private fiber and arbitrary
references. At delta zero it returns the signal exactly and implements
the **negative** of the requested reflection. At nonzero delta F need
not be Hermitian or preserve the clean sector; its actual inverse is
$`(U^\dagger ZUZ)^2`$, not a second forward F. The elementary
native example $`U=H_aT_aH_a`$ gives $`F=-iX_a`$, so omitting this
distinction would already fail on one wire.

Apply the same CNOT/Hadamard envelope as (34). Its initialized target
is $`-I_p`$, a fixed phase times the complete indicator; that phase
cancels in the serial promotion's source conjugations. No phase
correction is needed. Common forward/actual-inverse bounds are
$`T_V\le4T_U`$, $`G_V\le4G_U+2n+2m+8`$,
$`d_V\le4d_U`$, $`b_V\le2n-m+b_U`$. Each Z and its actual
adjoint are emitted as two positive S gates. Substitution in the
one-clean promotion gives

```math
\begin{aligned}
 T_{\rm out}&\le32T_U+48m+8,\\
 G_{\rm out}&\le32G_U+16n+98m+136,\\
 D_T({\rm out})&\le32d_U+32,\\
 b_{\rm out}&\le2n+b_U+2,\qquad b_{\rm clean}=0,\\
 \|\widetilde I_p-I_p\otimes I_{\rm borrowed}\|
   &\le64\sqrt\delta.                                \qquad(35e)
\end{aligned}
```

Again the final eight T gates can be omitted for even m. For requested
$`0\lt\tau\le1`$, $`\delta\le\tau^2/4096`$ suffices. All
source resources must be evaluated at that Gram precision. Actual
adjoints are retained throughout the promotion, including after leaked
components; no fresh clean signal or reset is used. The polynomial
fusion envelope is unchanged. This is a smaller conditional recovery
ledger, not a change to the unconditional complete-frame frontier.

The [exact four-call fixtures](../../tests/test_gram_four_call_readout.py)
check the full-column negative reflection with an arbitrary private
Hadamard fiber, the noncommuting-Gram identity (35c), and the literal
non-Hermitian inverse example. These finite identities do not establish
a normalized higher-band source or an asymptotic native schedule.

<a id="two-band-projector-gram-gap"></a>

## 13. One-clean projector-Gram sources need more than two bands

The Gram interface in Section 12 leaves the accepted partial isometry
free, but a sufficiently accurate **one-clean, two-band** source still
cannot supply the complete lookup projector. This conclusion permits
arbitrary intermediate address motion and arbitrary dirty width. It is
not the read-only source theorem, and it does not follow by treating a
Gram promise as a unitary-compression promise.

**Theorem.** Let U be any actual native circuit of T-depth at most two,
with one initialized clean signal a and all other inputs arbitrary.
Put $`K=J^\dagger UJ`$ and $`G=K^\dagger K`$. For any projector Pi,
if

```math
 \|G-\Pi/2\|\le\delta\le1/256,                       \qquad(36)
```

there is an actual native circuit of T-depth at most two, with the same
one clean signal, that returns that signal exactly and approximates
$`2\Pi-I`$ on the entire initialized input space with error at most
four delta. Consequently, for p at least three, no such U satisfies
(36) for the paired lookup projector of Section 12.

### 13.1 Principal angles and a literal Clifford rounding

Factor the actual two bands as $`U=U_B^\dagger U_A`$, where U_A and
U_B are one-band encoders on all physical wires. Define

```math
 A=U_AZ_aU_A^\dagger,\qquad
 B=U_BZ_aU_B^\dagger,\qquad C=AB,\qquad
 Q=(C+C^\dagger)/2.                                  \qquad(37)
```

Each one-band encoder maps a Pauli to a literal Clifford: conjugation
through a physical T/T-dagger band produces a product of local
Cliffords, and the outer native Cliffords preserve this property.
Thus A and B are Hermitian Clifford involutions and C is a Clifford.
All their actual phases are retained. Also
$`ACA=C^\dagger`$, and Q commutes with A.

Conjugate Q by U_A, and write the actual U in its original equally
sized clean blocks as

```math
 U=\begin{pmatrix}K&M\\N&L\end{pmatrix}.
```

Then $`U_A^\dagger Q U_A`$ has diagonal blocks $`2G-I`$ and $`2L^\dagger L-I`$.
The **diagonally opposite** blocks K and L have the same singular
values. Indeed, unitarity gives
$`K^\dagger K=I-N^\dagger N`$ and
$`LL^\dagger=I-NN^\dagger`$, and N is square. Therefore the spectrum
of Q is the spectrum of $`2G-I`$, repeated twice. This uses the sole
clean signal's equal half-dimensional sectors; it is not a general
growing-clean identity.

Since G is a positive contraction, (36) places its spectrum within
delta of $`\{0,1/2\}`$. For $`\delta\le1/4`$, every corresponding
eigenvalue t of Q obeys

```math
 t^2(1-t^2)\le4\delta.
```

For example, near the zero eigenvalue of G one has
$`1-t^2=4\lambda(1-\lambda)\le4\delta`$; near one half one has
$`t^2\le4\delta^2`$. Functional calculus for the unitary C gives
the exact polynomial identity

```math
 (C^4-I)^\dagger(C^4-I)=16Q^2(I-Q^2).
```

Consequently

```math
 \|C^4-I\|\le8\sqrt\delta\le1/2.                    \qquad(38)
```

Here is a phase-safe Clifford rounding, without a finite scalar-phase
assumption. Set $`D=C^4`$. For every Hermitian physical Pauli P,

```math
 \|DPD^\dagger-P\|\le2\|D-I\|\le1\lt\sqrt2.
```

D is Clifford, and distinct signed Hermitian Paulis have operator
distance at least square root of two. Thus D commutes with every Pauli
and is scalar. The identity $`ACA=C^\dagger`$ implies
$`ADA=D^\dagger`$, so that scalar is either plus or minus one.
Equation (38) excludes minus one. Therefore

```math
 C^4=I.                                               \qquad(39)
```

### 13.2 An actual two-band reflection, not a free spectral function

Equation (39) makes the spectrum of G a subset of $`\{0,1/2,1\}`$.
The bound $`\|G\|\le1/2+\delta\lt1`$ excludes one. Hence

```math
 \Pi'=2G\quad\hbox{is a projector},\qquad
 \|\Pi'-\Pi\|\le2\delta.                             \qquad(40)
```

The required physical reflection has the explicit form

```math
 R=-U_A^\dagger C^2U_A.                               \qquad(41)
```

C squared is an actual Clifford. Since
$`AC^2A=C^{-2}=C^2`$, R commutes with the original clean Z_a and
returns its clean sector exactly. Moreover,
$`C^2=2Q^2-I`$, and on that sector

```math
 J^\dagger RJ
 =I-2(2G-I)^2=4G-I=2\Pi'-I.
```

It follows on all original logical, output, private and reference
inputs that

```math
 \|RJ-J(2\Pi-I)\|\le4\delta\le1/64.                 \qquad(42)
```

This does not implement an uncharged polynomial in G: (41) is the
native implementation. Use the actual U_A and its adjoint and an
exact Clifford synthesis of C squared, retaining its scalar. The
minus sign is a literal Clifford scalar, for example $`(XZ)^2=-I`$.
There are two physical T bands, at most twice the T count of U_A,
no new workspace and an additional $`O(q^2)`$ Clifford synthesis cost
on q physical wires. The word is exactly Hermitian, so its actual
inverse can be emitted as the same word. These facts are needed only
to retain the native two-band class; no efficient classical synthesis
or resource-frontier improvement is inferred.

### 13.3 Applying the Boolean obstruction and selecting the next source

For the paired projector, $`2\Pi-I=\Omega_{x,C}\Omega_{x,D}`$.
The bitwise CNOT conjugation from Section 12 changes this to the complete
lookup phase Omega_p, using no additional clean work. This phase is
$`(-1)^{F(x,Y)}`$, where

```math
 F(x,Y)=\bigoplus_{u\in\mathbb F_2^p}Y_u[x=u],
 \qquad\deg_{\mathbb F_2}F=p+1.
```

The [robust two-band Boolean proof](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction)
applies directly to every real diagonal Boolean phase of degree greater
than three: its preserved-Z rounding, native eighth-root phases and
real-diagonal C3 conclusion do not require the phase to have the special
single-output form $`zf(x)`$. Thus its error floor
$`\sin(\pi/8)/(2\sqrt2)`$ also applies to Omega_p when p is at least
three. This is strictly larger than one over 64, contradicting (42).
No output bank is initialized or fixed to a test state in this reduction.

In particular, the Gram accuracy $`\delta=\tau^2/4096`$ requested
by Section 12.3 for $`0\lt\tau\le1`$ is already inside this excluded
two-band radius. Intermediate address motion therefore does not rescue
a one-clean two-band projector-Gram source. The next constructive
candidate at this interface must use at least three source bands, or
change its clean allocation together with a separately paid recovery
and promotion. This does not retire the projector-Gram mechanism:
three or more bands are outside (37)'s Clifford-product rounding, and
one clean signal still has the complete recurrence (35).
The conclusion is specific to this projector-Gram target and one clean
signal. It does not exclude a general unitary-compression source whose
Gram approximates $`I/2`$, nor a growing-clean projector-Gram source.

### 13.4 Higher-band discovery must match two spectral sectors

The principal-angle algebra also gives an exact diagnostic beyond two
bands. Factor any one-signal source as $`U=U_B^\dagger U_A`$ and
retain A,B,C from (37). Let
$`P_A=(I+A)/2`$ and $`\widehat\Pi=U_AJ\Pi J^\dagger U_A^\dagger`$.
Then the complete exact Gram contract is equivalent to the **joint**
conditions

```math
 \begin{aligned}
 (C+I)(C^2+I)&=0,\\
 P_A(I-C^2)/2&=\widehat\Pi.
 \end{aligned}                                          \qquad(43)
```

For necessity, the equally sized one-signal sectors and the
diagonally opposite singular-value identity in Section 13.1 give the
spectrum of Q as that of $`2G-I`$, repeated twice. If
$`G=\Pi/2`$, Q has only eigenvalues minus one and zero. Since
C is unitary and $`Q=(C+C^\dagger)/2`$, C has only eigenvalues
$`-1,i,-i`$, proving the first line. On those eigenvalues,
$`(I-C^2)/2=I+Q`$. Compression to the A-positive sector then gives
the second line. Conversely the first line implies this same polynomial
identity, and

```math
 P_A(I+Q)=2P_AP_BP_A=2U_AJGJ^\dagger U_A^\dagger,
 \qquad P_B=(I+B)/2.
```

The second line therefore gives $`G=\Pi/2`$ on the entire input
space. The one-signal rank balance is essential to the necessity
argument; (43) is not asserted for a general growing clean bank.

For an actual three-band source, U_A can be a two-band encoder and
U_B a one-band encoder. A, and hence C, need not be Clifford. Fourth
order closure $`C^4=I`$ alone still admits the unwanted Gram-one
sector (C=I); the first line of (43) removes it. Even that first line
does not identify the surviving sector with the lookup projector:
a logical input rotation changes it while preserving the spectrum.
The exact fixtures include both failures. Neither line supplies a free
native implementation of $`C^2`$ or an output endpoint. A candidate
must satisfy the sector matching and pay its emitted source, actual
inverse, precision and return together.

The exact principal-angle fixture in
[`test_nilpotent_pencil_observable_source.py`](../../tests/test_nilpotent_pencil_observable_source.py)
checks (39)--(42) for a small two-band source whose projector depends
only on a two-bit bank parity. All physical columns, scalar phases
and actual inverses are retained. This permitted low-degree example
checks the algebra rather than the excluded asymptotic source.
