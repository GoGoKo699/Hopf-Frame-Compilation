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

## 7. A punctured regular coordinate, with an exact spectator count

The source representation has a useful full-space normal form. Regard u
and v as field elements and split a field coordinate into the first p
polynomial coefficients z and the remaining $`k-p`$ coefficients t.
Define a computational permutation A by

```math
A|u,z,t\rangle=
\begin{cases}
|u,u(z+t)\rangle,&u\ne0,\\
|0,z+t\rangle,&u=0.
\end{cases}
\qquad
A^\dagger Q_x A|u,z,t\rangle
 =|u,z+x[\,u\ne0\,],t\rangle.                         \qquad(21)
```

Here z and t occupy disjoint polynomial coordinates; their sum denotes
that concatenated field element. Multiplication by nonzero u is a
bijection, while the second branch preserves every zero-sector state.
Thus $`k-p`$ spectator qubits can be factored out exactly. The inverse
coordinate change divides by u only on the nonzero branch. Neither this
nonlinear permutation nor its actual inverse has been paid natively.

This is the largest possible exact identity factor, even after an
arbitrary change of quantum basis. From (8),

```math
d_0=2^{k-p}(2^k+2^p-1),\qquad
d_{z\ne0}=2^{k-p}(2^k-1),\qquad
\gcd_z d_z=2^{k-p}.                                   \qquad(22)
```

The gcd follows because $`2^k-1`$ is odd and the two parenthesized
integers differ by $`2^p`$. An identity factor on r qubits divides every
joint character multiplicity by $`2^r`$. Consequently at least $`k+p`$
active bank qubits remain. Adding arbitrary identity helpers multiplies
all multiplicities equally and does not reduce this active minimum.
In particular, a sole p-qubit translation register with all other
qubits as spectators cannot represent the source at growing p.

There is also a robust version within exact group representations.
If an alternative representation has a mismatched character multiplicity,
choose a unit vector in the larger character sector orthogonal to the
smaller one. Character orthogonality gives average squared address
discrepancy two on that vector. The identity term contributes zero, so

```math
\max_x\|B_x-W_x\|\ge\sqrt{\frac{2m}{m-1}}.             \qquad(23)
```

This applies, in particular, to an exact representation with too many
identity spectators. It does not assume that an approximate circuit with
an inexact group law has this form, and it is not a native-depth bound.

## 8. Growing addresses require redistribution of the Bell fibers

The p=2 construction uses $`d=u`$. That convention cannot persist for
$`p\ge3`$, even if the entire basis inside the zero fiber is changed.
To see the precise missing operation, Bell-transform each pair of
physical wires differing in the first address bit. Collect the two-bit
Bell labels with second address bit zero into a k-bit string a, and
those with second address bit one into b, in matching order of the
remaining address bits. Put $`d=a+b`$.

Let $`R_0`$ be the physical-bank projector corresponding to $`d=0`$
in these Bell coordinates. Equivalently, on every four-wire coset of
$`\langle e_0,e_1\rangle`$, it is the joint positive eigenspace of
$`X^{\otimes4}`$ and $`Z^{\otimes4}`$. Its dimension is $`D=2^k`$.
The residual translation subgroup
$`H=\langle e_2,\ldots,e_{p-1}\rangle`$ has size
$`h=2^{p-2}`$ and preserves this projector. Every nonidentity element
permutes the residual Bell cells in pairs. On $`R_0`$, its trace is
$`2^{k/2}`$, so its negative eigenspace has dimension
$`(2^k-2^{k/2})/2`$.

The source projector $`P_0`$ onto $`u=0`$ is fixed pointwise by every
query. For a candidate E define

```math
P=EP_0E^\dagger,\qquad
\epsilon=\|P-R_0\|,\qquad
\eta=\max_x\|EQ_xE^\dagger-W_x\|.
```

The projectors have equal rank, hence
$`\epsilon=\|(I-R_0)EP_0\|`$. A negative target eigenvector inside
$`R_0`$ has projection onto P of norm at least
$`\sqrt{1-\epsilon^2}`$. Since the conjugated source acts identically
on P, the negative eigenspace yields

```math
\eta\ge2\sqrt{1-\epsilon^2}.                          \qquad(24)
```

In particular, mapping the whole source zero fiber onto the whole target
zero fiber forces error two. A small perturbation of that assignment
cannot repair it.

The necessary redistribution is extensive in mean as well. Let
$`\Pi_H=h^{-1}\sum_{y\in H}W_y`$. On $`R_0`$ this common-positive
projector has rank $`J=(D+(h-1)2^{k/2})/h`$. Averaging squared
intertwinement errors on P gives
$`\mathrm{tr}((I-\Pi_H)P)\le D(h-1)\eta^2/(2h)`$.
Because $`R_0`$ commutes with $`\Pi_H`$, its overlap with P is at
most J plus this error term. Therefore

```math
\frac{\mathrm{tr}(P_0E^\dagger(I-R_0)EP_0)}{D}
\ge\frac{h-1}{h}
       \left(1-2^{-k/2}-\frac{\eta^2}{2}\right).        \qquad(25)
```

At p=3, an exact decoder must export at least $`6/16`$ of the
zero-fiber mass on average; at growing p this fraction tends to one.
An exact decoder has worst-case leakage one. Tensoring both projectors
with arbitrary helper identities multiplies D and J equally and leaves
(24) and (25) unchanged.

### Even the weaker lookup probe cannot preserve this fiber

The same zero-fiber convention also fails when only the complete probe
(10) is required. Set $`S_x=EQ_xE^\dagger`$ and
$`\delta=\max_x\|S_xX_0S_x^\dagger-X_x\|`$. Then

```math
\delta+2\epsilon\ge1,\qquad p\ge3.                    \qquad(26)
```

Here is a full-space witness. On wires zero through seven take
$`\psi=\Phi^+_{01}\Phi^+_{23}\Phi^-_{45}\Phi^-_{67}`$,
where $`\Phi^\pm=(|00\rangle\pm|11\rangle)/\sqrt2`$;
use positive Bell pairs on the remaining wires. The equal signs within
each quartet put $`\psi\in R_0`$. With $`a=e_0,z=e_2`$,
$`X_a\psi=X_0\psi`$ and $`X_{z+a}\psi=-X_z\psi`$.
Also $`\|(S_x-I)\psi\|\le2\epsilon`$ because every $`S_x`$
fixes P pointwise.

Put $`\kappa=\delta+2\epsilon`$,
$`v=X_0\psi,w=X_z\psi,b=S_zv,c=S_{z+a}v`$. The probe errors
give $`\|b-w\|,\|c+w\|\le\kappa`$ and
$`\|(S_a-I)v\|\le\kappa`$. The exact source group law then gives
$`\|b-c\|\le\kappa`$. All four vectors are unit, hence

```math
2\kappa^2\ge\|b-w\|^2+\|c+w\|^2
 =4-2\mathop{\rm Re}\langle w,b-c\rangle
 \ge4-2\kappa,
```

which forces $`\kappa\ge1`$. Thus exact preservation of the zero
fiber has probe error at least one. Arbitrary helpers and references can
be tensored into the witness; no basis-state promise is used. These are
restrictions on this fiber assignment, not on unrestricted decoders.

## 9. Exact probe freedom and the remaining constructive choice

The weaker probe nevertheless has useful freedom. Suppose only
$`S_xX_0S_x^\dagger=X_x\otimes I`$ is exact on every arbitrary
helper input.
The group law implies $`S_xX_yS_x^\dagger=X_{x+y}`$ for all x and y.
In the joint X eigenbasis, S therefore permutes the bank eigenvalue
strings exactly as W does, with a unitary action $`C_x(b)`$ on the
helper fiber over each string b. The group law is the cocycle equation

```math
C_{x+y}(b)=C_x(W_yb)C_y(b).                            \qquad(27)
```

For nonzero x, precisely $`2^k`$ strings are fixed by $`W_x`$.
The source character equals $`2^k`$ times the helper dimension.
Its trace is the sum of $`\mathrm{tr}(C_x(b))`$ over these fixed
strings. The real part of each trace is at most the helper dimension; equality
of their sum forces $`C_x(b)=I`$ on every fixed string. Thus the
cocycle is trivial on every stabilizer. The literal source character is
essential here; untracked address-dependent phases would not justify
this saturation.

Choose one representative $`b_0`$ per orbit and a translation
$`t_b`$ taking it to b. Define $`D_b=C_{t_b}(b_0)`$; the stabilizer
identity makes this independent of the choice of translation. Then
$`C_x(b)=D_{W_xb}D_b^\dagger`$, proving the exact characterization

```math
S_x=D(W_x\otimes I)D^\dagger,
\qquad [D,X_y\otimes I]=0\quad\text{for every }y.       \qquad(28)
```

The gauge may act nontrivially within each helper fiber; it need not be
scalar diagonal on the full space. Consequently $`D^\dagger E`$ is a
full simultaneous intertwiner.
Conversely, multiplying any such intertwiner by any D in this
commutant preserves every probe. D is an abstract gauge and has no
free native implementation; it can cancel from the desired probe,
so optimizing the probe may still be cheaper than synthesizing one
particular full intertwiner. No approximate version of (28) is asserted.

A concrete recoupling mechanism is visible in the orbit modules.
The source has $`2^k`$ fixed basis states and
$`(2^k-1)2^{k-p}`$ free translation orbits. The target also has
periodic strings, whose stabilizers are proper nonzero subgroups.
For subgroups $`H\subset K\subseteq G=\mathbb F_2^p`$ with
$`\dim(K/H)=2`$, the three intermediate subgroups obey

```math
\bigoplus_{H\subset L\subset K}\mathbb C[G/L]
 \simeq \mathbb C[G/H]\oplus2\mathbb C[G/K].           \qquad(29)
```

This follows by Fourier transforming the four-element quotient and
counting each character, then inducing to G. In that quotient it is the elementary
conversion of three two-point modules into one regular module and
two trivial modules. At p=3, the source orbit counts are sixteen fixed
and thirty free orbits; the target has two fixed, seven size-two,
fourteen size-four and twenty-three free orbits. Their character ranks
are $`(46,30,\ldots,30)`$, but their computational orbits must be
redistributed. Equation (29) identifies a local representation-theoretic
move, not a native routing algorithm or a cost bound.

The [recursive field-shear decoder](RECURSIVE_FIELD_SHEAR_DECODER.md)
now supplies a constructive all-size choice of E using finite-field
coordinate changes and Fourier transforms. Its intermediate zero fiber
is common-fixed for two previously decoded half-banks; the subsequent
half-bank decoders redistribute the physical Bell fibers. Retaining
$`d=u`$ in those intermediate coordinates therefore differs from retaining
the physical projector $`R_0`$ excluded above.

Native synthesis of the recursion's punctured affine relabeling and
field multiplication, or a cheaper joint choice using the probe freedom
(28), is the next operation. These coordinate maps have no free depth or
helper budget. Computing a canonical stabilizer and matching all
character fibers, or merely listing the moves (29), likewise leaves the
cost unpaid. Repairing only the basis inside the physical Bell zero
branch remains excluded by (26).

## 10. A diagonal seed in an invariant stabilizer basis cannot suffice

There is a separate exact restriction on the starting reflection. Let
$`Q_x\otimes I_h`$ act on the m bank qubits and h arbitrary helpers.
Choose any stabilizer basis that every query permutes up to phases;
equivalently, its Pauli Lagrangian is invariant under every query's
binary symplectic action. The basis may mix bank and helper wires.
Let R be any Hermitian involution diagonal in this basis. Its Boolean
sign function may have arbitrary degree and complexity.

For $`p\ge2`$, no fixed unitary E can satisfy

```math
E Q_xR Q_x^\dagger E^\dagger=X_x\otimes I_h
\qquad\text{for all }x\in\mathbb F_2^p.               \qquad(30)
```

Thus an arbitrarily nonlinear diagonal seed in preserved stabilizer
coordinates does not solve the complete probe. E need not be Clifford,
monomial or shallow. This statement concerns exact operators; no
constant-error extension is asserted.

### Common affine fixed points have the wrong possible cardinalities

Write $`\mathcal V`$ for the binary Pauli space of dimension
$`2(m+h)`$. For the field shear, the symplectic deviations
$`N_x=F_x-I`$ of every nonzero address have the same image W and
kernel $`W^\perp`$. Here W is the m-dimensional isotropic space
spanned by X on the v bank and Z on the u bank; helper labels lie in
the kernel. These facts only require invertibility of every nonzero
pencil matrix, so they also apply to repeated smaller-field blocks.

Let $`\mathcal L`$ be the invariant Pauli Lagrangian of the chosen
stabilizer basis, of dimension $`m+h`$, and set

```math
r=\dim(W/(W\cap\mathcal L)).
```

Computational labels for this basis identify with
$`\mathcal V/\mathcal L`$ through the symplectic pairing. Their
query permutations can be affine, even though the original field
shears are linear. The linear deviation of each generator has image
in the same r-dimensional projection of W. Consequently their common
linear fixed space has dimension t with $`t\ge m+h-pr`$.

That space also contains the projection of $`W^\perp`$. The
restriction of the symplectic pairing between $`\mathcal L`$ and W
has rank r, because $`\mathcal L^\perp=\mathcal L`$. Therefore

```math
\dim(\mathcal L\cap W^\perp)=m+h-r,
\qquad
\dim((W^\perp+\mathcal L)/\mathcal L)=h+r.
```

Combining the two lower bounds gives

```math
t\ge h+\max\{m-pr,r\}
 \ge h+\left\lceil\frac{m}{p+1}\right\rceil.          \qquad(31)
```

The affine common fixed-point equations either have no solution or
have one coset of this linear fixed space. Their number of basis
labels is therefore zero or at least
$`2^{h+\lceil m/(p+1)\rceil}`$. Existence of a common affine fixed
label must not be assumed.

### Complete lookup forces exactly two helper fibers

All the reflections $`R_x=Q_xRQ_x^\dagger`$ are diagonal in the
chosen basis. If (30) held, their joint sign patterns would be all
$`2^m`$ strings, with exactly $`2^h`$ basis columns per pattern.
The source group law permutes these patterns by the regular
translation of their m coordinates.

For any nonzero x, precisely $`2^{m/2}`$ sign patterns are fixed by
this translation. The trace of $`Q_x\otimes I_h`$ is exactly
$`2^{m/2+h}`$. Its trace comes only from these fixed-pattern fibers,
each of dimension $`2^h`$. As in (27)–(28), saturation of the sum of
their maximum possible real traces forces the query to be the
identity on every such fiber. This accounts for phases and arbitrary
helper mixing; a permutation of columns inside a fixed fiber would
not suffice.

Every computational label fixed by all query permutations must have
a constant joint sign pattern. Conversely, both constant-pattern
fibers are fixed pointwise by every query by the saturation argument.
The common permutation fixed set must therefore have exactly
$`2^{h+1}`$ columns. But for $`p\ge2`$,
$`\lceil2^p/(p+1)\rceil\ge2`$, contradicting both possible
cardinalities in (31). This proves (30) impossible.

### A small affine example with no common fixed label

The empty case in (31) occurs already at p=2. Write
$`\mathbb F_4=\mathbb F_2[t]/(t^2+t+1)`$ and use the trace-Fourier
coordinates on v, where the query is the phase
$`(-1)^{\mathrm{Tr}(\lambda uw)}`$. The following is a Clifford
stabilizer basis, indexed by $`b,c\in\mathbb F_4`$:

```math
|b,c\rangle_{\mathcal L}
 =\frac12\sum_{a\in\mathbb F_4}
      (-1)^{\mathrm{Tr}(ca)}|u=a,w=ta+b\rangle.
```

Since $`\mathrm{Tr}(za^2)=\mathrm{Tr}(z^2a)`$ in this field,
the query permutes its labels, with literal phase one, by

```math
(b,c)\longmapsto(b,c+(\lambda t)^2+\lambda b).        \qquad(32)
```

The generator $`\lambda=1`$ fixes only labels with $`b=t+1`$;
$`\lambda=t`$ fixes only those with $`b=1`$. Thus each nonzero
query has four fixed labels, but there is no common fixed label.
There are six two-point orbits and one four-point orbit. The original
computational basis instead has four common fixed labels. Neither
basis can support the two constant sign fibers required by (30).

This screen does not exclude seeds that are nondiagonal in every
query-invariant stabilizer basis, general non-stabilizer joint
coordinates, or a different source family. A synthesis using one of
those mechanisms still needs its complete native ledger.

<a id="7-verification-and-remaining-operation"></a>
<a id="10-verification-and-remaining-operation"></a>

## 11. Verification and remaining operation

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

The continuation fixtures check (21) on all p=3 columns and addresses,
including the zero branch and the one spectator qubit; the character
multiplicity gcd and all computational orbit counts agree with (22).
Sparse exact p=3 checks verify the Bell zero-fiber code, its residual
negative eigenvector, invariant rank ten, mean export $`3/8`$ and the
opposite pair signs used in (26). A separate complete operator fixture
reconstructs (28) with noncommuting one-helper blocks, checking every
probe and the helper-dependent cocycle rather than assuming scalar phases.
An explicit four-wire Clifford basis checks (32) on every query and
column, including its empty common affine fixed set, its nontrivial
linear fixed space and both orbit counts. This guards against replacing
the zero-or-coset alternative in (31) by an unjustified fixed-point
existence claim.

For the full-intertwiner route, the first unpaid all-size operation is
native synthesis of E and its actual inverse within (11). The alternative
[encoded dirty-copy interface](ENCODED_DIRTY_COPY_INTERFACE.md) can use
a suitable transported companion without synthesizing that particular
decoder; its complete copy and return operations still need payment.
The source query and representation-theoretic existence result are
complete. Neither route currently supplies a new depth bound; R96 and
the selected manuscript lane are unchanged.
