# A finite-field shear with the full lookup character

[Complete indicator contract](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) ·
[Returned cubic batches](DIRTY_CUBIC_PHASE_COMPILER.md) ·
[Queried Clifford boundary](CLIFFORD_ORBIT_QUERY_RANK.md)

A single bipartite shear can have exactly the same simultaneous unitary
representation as address-controlled permutation of the whole lookup
bank. Its coherent native implementation has eight T bands and
$`O(mp^2)`$ T gates, arbitrary returned helpers and total gates. The fixed
unitary intertwining the two representations exists exactly at every
size, but its all-size native synthesis is unpaid. It must be nonmonomial and
non-Clifford for $`p\ge2`$.

This is a sufficient route for the original commuting lookup family.
It does not identify that family with the noncommuting Pauli targets
obtained after the minimal-tree pullback. It supplies neither an
all-size shallow decoder nor an improvement to R96.

## 1. A full-rank linear pencil and one native shear

Let $`p\ge1`$, $`m=2^p`$ and $`k=m/2`$. In particular $`p\le k`$.
Choose a monic irreducible polynomial f of degree k over the binary
field. In the polynomial basis of $`\mathbb F_2[t]/(f)`$, let C be
the column-oriented multiplication matrix of t. Define

```math
M_x=\sum_{i=0}^{p-1}x_iC^i,
\qquad
Q_x|u,v\rangle=|u,v+M_xu\rangle,
\qquad u,v\in\mathbb F_2^k.                         \qquad(1)
```

The address is an arbitrary p-bit register, while u and v are two
arbitrary k-qubit banks. For nonzero x, its polynomial has degree
strictly below k and is nonzero in the field. Thus $`M_x`$ is invertible.
There is no exceptional address other than zero and no omitted
zero-input sector. On every computational column, with literal phase
one,

```math
Q_xQ_y=Q_{x+y},\qquad Q_x^2=I,\qquad Q_0=I.           \qquad(2)
```

All arithmetic in (1) defines a circuit; no field-element state is
prepared or loaded. The binary polynomial f is a classical description.
For example, it may be selected by a finite exhaustive irreducibility
test. The bounds below count the emitted native circuit, not a claim
of polynomial classical search time for either f or the decoder.
Once f is given, its companion powers and their sparse occurrence list
are computed directly from its coefficients.

### A paid forward word and its actual inverse

For $`0\le i\lt p\le k`$, the first $`k-i`$ columns of $`C^i`$
are unreduced monomials and each has one nonzero entry. The remaining
i columns have at most k entries each. If s counts raw occurrences,
before cancellation between any equal terms, then

```math
s=\sum_{i=0}^{p-1}\mathrm{nnz}(C^i)
 \le kp+\frac{(k-1)p(p-1)}2=O(mp^2).                  \qquad(3)
```

For each nonzero entry $`(C^i)_{\ell j}`$, the query adds
$`x_i u_j`$ to $`v_\ell`$. All address and u controls are disjoint
from the v targets. Hadamards on the k v wires therefore turn the
entire query into one cubic sign batch

```math
(-1)^{\sum_{i,j,\ell}(C^i)_{\ell j}x_i u_jv_\ell}.
                                                               \qquad(4)
```

Use one private arbitrary helper for each raw cubic occurrence and the
closed parity-cube word from the
[returned cubic compiler](DIRTY_CUBIC_PHASE_COMPILER.md). Each helper
receives eight T or inverse-T gates and eight CNOTs; every physical T
band has one phase on each distinct helper. The cubic phase is
independent of every helper value, and the closed walk returns each
helper with literal phase one relative to the intended sign. Shared
controls and repeated targets cause no serialization of these eight
T bands.

The coherent query $`Q=\sum_x|x\rangle\langle x|\otimes Q_x`$
and its actual reversed native word obey

```math
T=8s,\qquad D_T\le8,\qquad
G_{\rm all}=16s+2k,\qquad
w=p+m+s,\qquad w_{\rm clean}=0.                    \qquad(5)
```

Here $`G_{\rm all}`$ counts both Clifford and T gates; the Clifford
count is $`8s+2k`$. In particular, T count, width and total gate count
are $`O(mp^2)`$, within the allowed polynomial discovery overhead.
The inverse is paid by reversing these gates and inverting their
phases, not by substituting an abstract assertion that Q is an
involution. No reset, measurement, QRAM or clean flag is used.

## 2. Simultaneous character equivalence to the full bank permutation

Index the m bank wires by $`a\in\mathbb F_2^p`$. Let $`W_x`$ be
the phase-free physical wire permutation taking wire a to wire
$`a+x`$. Then

```math
W_xW_y=W_{x+y},\qquad W_x^2=I,\qquad
W_xX_0W_x^\dagger=X_x.                                \qquad(6)
```

For nonzero x, the shear in (1) fixes precisely the basis strings
with $`u=0`$, leaving all $`2^k`$ v values. The wire translation
$`W_x`$ consists of k disjoint swaps; a fixed basis string is constant
on each pair, also giving $`2^k`$ choices. Consequently the complete
group characters agree:

```math
\mathrm{tr}(Q_x)=\mathrm{tr}(W_x)=
\begin{cases}
2^{2k},&x=0,\\
2^k,&x\ne0.
\end{cases}                                                   \qquad(7)
```

This compares every product of the generators because both families
obey the group law. Matching isolated generator spectra alone would
not establish simultaneous equivalence.

For character $`\chi_z(x)=(-1)^{z\cdot x}`$, the orthogonal joint
projectors are $`m^{-1}\sum_x\chi_z(x)Q_x`$ and the corresponding
sum with W. Their equal ranks are

```math
d_0=\frac{2^{2k}+(m-1)2^k}{m},\qquad
d_z=\frac{2^{2k}-2^k}{m}\quad(z\ne0).                \qquad(8)
```

These are nonnegative integers; $`p\le k`$ also makes integrality
immediate. Match orthonormal bases in each character sector. This
constructs one fixed, address-independent unitary E satisfying

```math
E Q_x E^\dagger=W_x\qquad\text{for every }x.           \qquad(9)
```

The two commuting families are real symmetric involutions, so the
bases and E can be chosen real. Equation (9) is an exact existence
and dense-construction statement on the full $`2^m`$-dimensional
bank Hilbert space. It is not a bound on a native word for E.

### What a paid decoder would buy

Suppose a native E, including any arbitrary returned helper registers,
satisfies (9) tensored with their identity. Its actual inverse and the
paid query implement $`V=E Q E^\dagger`$, and the complete indicator is

```math
V X_0 V^\dagger
 =\sum_x|x\rangle\langle x|\otimes X_x=I_p.          \qquad(10)
```

Thus every address, bank, helper and external reference input is
included. If E and its inverse each use at most t T gates, d T bands,
g total gates and h additional helper wires, then a conservative
ledger for (10), with native $`X=HSSH`$, is

```math
T\le4t+16s,\qquad D_T\le4d+16,\qquad
G_{\rm all}\le4g+32s+4k+4,\qquad
w\le p+m+s+h.                                         \qquad(11)
```

A terminating uniform decoder compiler, at each fixed tolerance,
with $`t,h=O(mp^{O(1)})`$, $`g=O(m^2p^{O(1)})`$ and
$`d=o(\log(p+2))`$ would therefore supply the missing sufficient
route through the retained polynomial-overhead fusion. The constants
and fixed polynomial exponents may depend on the requested tolerance.
For approximate intertwinement, the uniform operator bound
$`\|EQ_xE^\dagger-W_x\|\le\eta`$ gives at most $`2\eta`$ error
in (10), by the two-factor conjugation estimate. The implemented
unitary and its actual inverse are still required. A dense E, a list
of its eigenvectors or equality only on chosen input states does not
pay this operation.

## 3. Two boundaries for the fixed decoder

These obstructions concern simultaneous intertwinement (9). They do
not exclude every circuit that has the weaker probe (10), nor every
different shallow native source.

### Computational monomial decoders have a fixed error gap

The source has $`2^k`$ computational basis strings fixed by every
$`Q_x`$: exactly $`u=0`$. The target has only two strings fixed by
every $`W_x`$, since all wire translations together act transitively
on the m wire labels. They are the all-zero and all-one strings.
For $`p\ge2`$, $`2^k\gt2`$.

A monomial E maps these source fixed strings to distinct computational
strings, up to phases. At least one image b is moved by some target
$`W_x`$. But $`EQ_xE^\dagger`$ fixes b literally. The two outputs
are orthogonal, proving

```math
\max_x\|EQ_xE^\dagger-W_x\|\ge\sqrt2
\quad\text{for every computational monomial }E,
\quad p\ge2.                                             \qquad(12)
```

Arbitrary identity helpers multiply both common-fixed counts by
their dimension, so the argument survives monomial decoders mixing
those helpers with the bank. This explains why a reversible field
relabeling, including a separately handled zero branch, cannot by
itself implement (9). It is not an optimality claim for the gap.

### Clifford simultaneous conjugacy is impossible for p at least two

On binary Pauli labels the shear has matrix $`F_x=I+N_x`$.
Its X-label block sends u labels into v labels through $`M_x`$;
its Z-label block sends v labels into u labels through
$`M_x^{\mathsf T}`$. Each block annihilates the image of every
other such deviation. Therefore

```math
N_xN_y=0\qquad\text{for all }x,y.                      \qquad(13)
```

For wire translations, let $`P_x`$ act on the m X-coordinate labels.
For independent x and y,

```math
(P_x-I)(P_y-I)e_0=e_0+e_x+e_y+e_{x+y}\ne0.             \qquad(14)
```

A fixed Clifford conjugation induces simultaneous symplectic
similarity, which preserves products of deviations. Equations
(13) and (14) contradict such an E for $`p\ge2`$. Identity helper
blocks do not change the contradiction. This is an exact decoder
boundary and not an unrestricted native-depth lower bound.

It also has a uniform approximation gap. The retained
[Clifford trace-gap argument](CLIFFORD_ORBIT_QUERY_RANK.md) gives
distance at least $`\gamma=\sqrt{2-\sqrt2}`$ from the identity for
every nonscalar Clifford. The obstruction above forces at least one
relative Clifford $`W_x^\dagger E Q_xE^\dagger`$ to have nonidentity
symplectic action. Therefore every Clifford decoder, including one
with arbitrary identity helpers, has maximum intertwinement error
at least $`\gamma`$. No optimality of this constant is asserted.

## 4. The one-address native word

For $`p=k=1`$, take $`f=t+1`$. The shear is
$`Q_x=\mathrm{CX}(u;v)^x`$, and the fixed Clifford
$`E=\mathrm{CX}(v;u)`$ satisfies

```math
E\,\mathrm{CX}(u;v)\,E=\mathrm{SWAP}(u,v).            \qquad(15)
```

The coherent source is the chronological word E, a CCX controlled by
x and u with target v, then E. Expand that CCX using the retained
[two-band returned-helper primitive](WEAK_PAULI_ORBIT_COMPRESSION.md).
The resulting source uses eight T gates, two T bands, twenty Clifford
gates and one arbitrary returned helper, hence twenty-eight total gates. The
complete echo in (10) uses sixteen T gates, at most four bands and
sixty total gates. Both actual reversed words are checked literally.

## 5. An exact two-bit joint-eigenbasis fixture

For $`p=2`$, use $`f=t^2+t+1`$, giving

```math
C=\begin{pmatrix}0&1\\1&1\end{pmatrix},\qquad
(d_0,d_1,d_2,d_3)=(7,3,3,3).                            \qquad(16)
```

The finite test constructs exact real orthonormal bases from
computational permutation orbits. For an orbit represented by b,
with stabilizer H, each character z trivial on H contributes

```math
|b,z\rangle=
\frac1{\sqrt{|\mathcal O_b|}}
\sum_{x\in\mathbb F_2^p/H}\chi_z(x)|Q_xb\rangle,      \qquad(17)
```

and likewise for W. Representatives and cosets are chosen
deterministically. Orbit sizes are one, two or four, so the exact
entries lie in $`\mathbb Q(\sqrt2)`$. Sorting columns by character
gives basis matrices $`U_Q,U_W`$ and the explicit dense
$`E=U_WU_Q^\dagger`$. The test checks full orthogonality, all four
$`16\times16`$ identities in (9), and their complete X-probes,
with literal phases. These block identities also cover coherent
addresses and any external reference.

The source has four fixed orbits and three four-point orbits. The
target has two fixed orbits, three two-point orbits and two
four-point orbits. Their equal character counts conceal this
different orbit structure; converting the remaining six-dimensional
pieces requires quantum superposition. The orbit Fourier fixture
does not compile that change into native gates.

## 6. A paid four-wire decoder for two address bits

There is also a native p=2 decoder, independent of the canonical
eigenbasis choice in (17). It uses only the four bank wires. The
following construction explains its full action, including the zero
sector; it does not extend the native bound to growing p.

Apply Hadamards to v in the source. The two generator characters
become $`\chi_0=u\cdot v`$ and $`\chi_1=(Cu)\cdot v`$, where
these binary exponents determine the eigenvalue signs. For the
target, use Bell coordinates a and b for physical wire pairs
$`(0,1)`$ and $`(2,3)`$. The first wire translation acts by the
phase $`(-1)^{a_0a_1+b_0b_1}`$; the second exchanges a and b.

Put $`d=a+b`$ and, for $`d\ne0`$, write
$`a=r d+sCd`$. The symmetric invertible matrix
$`F_d=[d,Cd]`$ is I, C and $`C^2`$ for binary d=1, 2 and 3,
respectively. The target generators now act as

```math
W_1:\ (-1)^{s+d_0d_1},\qquad
W_2:\ r\longmapsto r+1\quad(d\ne0).                   \qquad(18)
```

Both target generators act trivially when $`d=0`$. Set $`F_0=I`$
and relabel the source by

```math
d=u,\qquad r=(Cu)\cdot v,\qquad s=u\cdot v+d_0d_1
\quad(d\ne0).                                           \qquad(19)
```

On $`d=0`$, use $`(r,s)=(v_1,v_0)`$, preserving all four
arbitrary input states. Apply H to r unconditionally: both generators
are the identity on this entire four-dimensional zero sector, so any
basis change there preserves their action. No nonzero flag is needed. Then
convert back using $`a=F_d(r,s)`$, $`b=a+d`$, and apply the two
Bell basis transforms. Equation (18) has the required signs
$`(-1)^{\chi_0}`$ and $`(-1)^{\chi_1}`$. Hence both generators
and every product satisfy (9), including the full zero sector.

The native word implements each $`F_d`$ with two controlled powers
of C, using the identity
$`F_d=C^{\,2d_0+(d_0\mathbin\oplus d_1)}`$. Each controlled C
or its inverse is two ordinary Toffolis. A separate Toffoli adds
$`d_0d_1`$ in (19). Expand those Toffolis using the retained four-band
native word, and implement every inverse literally. The intervening
Fourier change is just one ordinary Hadamard.

The resulting E and its actual inverse each contain 63 T or
inverse-T gates, 107 CNOTs and 23 Hadamards:
193 gates in 36 chronological T bands, with no helper wires.
Combining this decoder with (5), where $`s=5`$, gives the common
forward/actual-inverse ledgers:

| Native word | T gates | T bands at most | Total gates | Total wires | Clean wires |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fixed four-bank decoder E | 63 | 36 | 193 | 4 | 0 |
| Coherent source $`V=E Q E^\dagger`$ | 166 | 80 | 470 | 11 | 0 |
| Complete indicator $`V X_0 V^\dagger`$ | 332 | 160 | 944 | 11 | 0 |

The source and indicator use five arbitrary returned cubic masks;
their total width also includes two address wires. These are finite
concatenation bounds, not finite optimality or asymptotic claims.

The decoder also gives a concrete non-Clifford fixed seed
$`R_0=E^\dagger X_0E`$. With basis integers ordered by wire weights
1, 2, 4 and 8, the exact test obtains

```math
R_0X_0R_0|0\rangle=
\frac{|1\rangle+|5\rangle-|9\rangle+|13\rangle}{2}.
                                                               \qquad(20)
```

A Pauli maps each basis column to one basis column. Thus
$`R_0X_0R_0`$ is not Pauli, certifying that $`R_0`$ is not a
Clifford unitary. This is a finite joint source-and-decoder witness;
no all-size seed or decoder efficiency follows from it.

## 7. Verification and remaining operation

`tests/test_field_shear_intertwiner.py` checks the displayed field
examples, all nonzero pencil ranks and the sparse native ledger
through p=4. It exhausts the full characters and common-fixed sets
through p=3. For p=2 it evaluates the emitted cubic core on all 2048
address, bank and helper columns, including its actual inverse;
the two explicit Hadamard walls turn this exact diagonal operator
into (1). It also checks the p=1 native source and complete indicator
as full matrices, and both p=2 intertwiners as full operators.
The native p=2 decoder's actual inverse equals its exact adjoint;
all four source blocks and complete probes agree literally. Its
physical T bands, the non-Clifford seed witness and the complete
source/indicator ledgers are checked as well. The component operator
identities include every one of the five arbitrary cubic masks.

The first unpaid all-size operation is native synthesis of E and
its actual inverse with the width, gate and depth budgets in (11).
The source query and the representation-theoretic existence result
are complete. A fast decoder remains open; R96 and the selected
manuscript lane are unchanged.
