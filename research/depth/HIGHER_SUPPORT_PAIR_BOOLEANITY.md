# Higher-support physical pairs and the Boolean completion test

[Original physical-pair scout](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md) ·
[Center transport and paid resynthesis](FOUR_BAND_CENTER_TRANSPORT.md) ·
[Two-layer diagonal rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md)

The fully matched support-two mechanism has an actual paid two-band
resynthesis. This note locates a precise boundary of that reduction
and then tests the complete Boolean condition beyond it. An isolated
pair with larger active support, or with coupling to a passive endpoint
axis, cannot be compressed to one native band by Clifford dressing.
Nevertheless, a literal five-wire source with this obstruction and
genuine physical-center motion still fails every diagonal completion
by input/output Clifford changes.

The reusable test is the Pauli-preserving domain of the complete
occupied unitary. Exact clean return alone does not supply enough
Pauli transport for a diagonal target. All pairs below are original
physical pairs of independent native lists. No unrestricted four-layer
degree theorem, target-depth lower bound or resource-frontier
improvement follows.

## 1. A necessary Pauli-domain condition for diagonal completion

For an actual unitary W on n logical wires, define the binary space

```math
M_W=\{p:WP_pW^\dagger\text{ is a Pauli up to scalar}\}.
```

It is a linear subspace of the Pauli-label space, with its inherited
symplectic form. An occupied unitary must first be justified on its
actual preserved code; this definition does not license a quotient
of a nonpreserved physical space.

For arbitrary logical Cliffords L_A,L_C,

```math
G=L_C^\dagger WL_A
\quad\Longrightarrow\quad
M_G=L_A^\dagger M_WL_A.
```

The notation on the right denotes the Clifford action on Pauli labels.
If G is computationally diagonal, M_G contains every logical Z.
Therefore M_W must contain an n-dimensional isotropic subspace.
Equivalently, its largest isotropic dimension must obey

```math
\nu(M_W)
 =\dim M_W-\frac12\mathrm{rank}(\Omega|_{M_W})
 \ge n.
```

This is a necessary condition, not a converse synthesis theorem.
In particular it applies to Boolean diagonal targets, with their
literal common scalar retained.

Tensoring W with s arbitrary identity spectators adds their full
symplectic Pauli space. It increases both the required diagonal
dimension and the maximal isotropic dimension by s. A deficit in
this test therefore persists under arbitrary endpoint Clifford mixing
with those spectators, including wires designated as dirty.

### A controlled-family version

Suppose W preserves c computational control wires and has target
fibers W_z. For every individual fiber,

```math
\nu(M_W)\le c+\nu(M_{W_z}).
```

To prove this without assuming that Pauli transport fixes each fiber,
project M_W onto its control X mask, of rank k. Its kernel contains
the c-dimensional control Z space C. After deleting those Z coordinates,
write the remaining target-label space as K. A zero-X-mask Pauli has
a Pauli image in each fiber, so $`K\subseteq M_{W_z}`$ for every z.
Choose k control-Z/lift hyperbolic pairs. Adding control-Z vectors to
K removes its pairings with the lifts without changing its own form.
The remaining c minus k control directions are radical. Hence

```math
\begin{aligned}
\dim M_W&=c+\dim K+k,\\
\mathrm{rank}(\Omega|_{M_W})
 &=2k+\mathrm{rank}(\Omega|_K),\\
\nu(M_W)&=c+\nu(K)\le c+\nu(M_{W_z}).
\end{aligned}
```

No classification of cross-fiber transport or equality between K and
any full fiber domain is required. In particular, if one fiber has
trivial Pauli-preserving domain, the control Z space and this bound give

```math
\nu(M_W)=c.
```

## 2. The exact boundary of the isolated-pair resynthesis

Let

```math
E_A=\prod_jT_{a_j}^{\epsilon_j},\qquad
J=E_APE_A^\dagger
```

be an actual independent commuting native list and its clean probe.
The a_j are signed Hermitian Paulis, and each epsilon_j is plus or
minus one. A clean-active axis anticommutes P; a passive axis commutes
P. Suppose b and $`b'=JbJ`$ form a genuine commuting original pair.
Its difference has support S of even size k at least two among the
actual clean-active A axes. Let T contain the t passive A axes
anticommuting b, and put m=k+t. Every other A axis commutes b.

If the displayed second native axis is c with $`JbJ=\sigma c`$,
retain the literal representative conversion

```math
T_c^{\theta_0}
 =\omega^{\theta_0(1-\sigma)/2}T_{b'}^{\sigma\theta_0},
\qquad \omega=e^{i\pi/4}.
```

Below eta and theta are the effective orientations. Put
$`PbP=\kappa b`$ and $`a_S=\prod_{j\in S}a_j`$. The signed
pair formula is

```math
\begin{aligned}
b'&=\lambda a_Sb,\qquad r=bb'=\lambda a_S,\\
\lambda&=\kappa(-1)^{k/2}\prod_{j\in S}\epsilon_j.
\end{aligned}
```

Passive A factors cancel out of J. Native independence determines S
from its active displacement. Since k is even, r is a signed
Hermitian Pauli commuting with every A axis.

Consider just the actual pair factor

```math
U=E_A^\dagger(T_b^\eta T_{b'}^\theta)E_A.
```

For any i in S, both pair axes anticommute a_i, so

```math
Ua_iU^\dagger
 =\frac{a_i-\eta\theta r a_i
   -iE_A^\dagger(\eta b+\theta b')E_Aa_i}{2}.
```

The literal endpoint expansion is

```math
E_A^\dagger bE_A
 =2^{-m/2}\prod_{j\in S\cup T}(I+i\epsilon_j a_j)b.
```

Multiplication by $`\eta(1+\eta\theta r)`$ makes the coefficient
of a subset H of S acquire the factor

```math
1+\eta\theta\kappa(-1)^{|H|}.
```

Exactly half the subsets of S survive, with doubled magnitude;
all passive choices survive. The third summand in the Pauli image
therefore has exactly $`2^{k-1+t}`$ distinct terms, each of magnitude
$`2^{-m/2}`$ after the outer factor. Its affine support lies in
$`b+a_i+\mathrm{span}(A)`$. It is disjoint from the first two
labels a_i and r+a_i: b cannot lie in the isotropic A span while
anticommuting a_i. The first two labels are distinct because the pair
difference is nonzero. Thus the exact support cardinality is

```math
2+2^{k-1+t},
```

with the first two coefficient magnitudes equal to one half.
Physical native independence, not quotient independence, rules out
all collisions.

For m greater than two, this cardinality is not a power of two.
A Clifford has power-of-two Pauli support with equal nonzero
magnitudes: expand its stabilizer Choi state in the Bell basis,
itself a Clifford transform of the standard basis. Hence the displayed
Pauli image is not Clifford, and U is outside C3. It cannot be a
single native band with arbitrary Clifford words on either side.
Common scalar phases do not alter this conclusion.

The remaining possibility is k=2,t=0, where the
[literal support-two identity](FOUR_BAND_CENTER_TRANSPORT.md#f-a-paid-two-band-resynthesis-of-the-fully-matched-mechanism)
does resynthesize the pair. This proves a boundary of that local
operation. It does not exclude cancellation with other actual factors
in the complete source, nor impose a lower bound on its optimal depth.

## 3. A literal center-moving source with exact clean return

Use five physical wires and

```math
\begin{aligned}
E_A=E_C&=\prod_{j=1}^5T_{Z_j},\\
P&=X_1X_2X_3X_4X_5,\\
J&=E_APE_A^\dagger
  =\prod_{j=1}^5\frac{X_j+Y_j}{\sqrt2}.
\end{aligned}
```

Use the following actual middle lists, with positive Hermitian
representatives and positive T orientations:

```math
\begin{aligned}
E_B&=T_{X_1X_2X_3X_4}T_{Y_1Y_2Y_3Y_4},\\
E_D&=T_{X_1X_5}T_{Y_1Y_5}.
\end{aligned}
```

Each list is commuting and independent. J swaps each displayed pair
with positive sign, so both bands commute J exactly. The original
physical B center is

```math
r=(X_1X_2X_3X_4)(Y_1Y_2Y_3Y_4)=Z_1Z_2Z_3Z_4.
```

Both D axes anticommute r. Their independent native conjugation gives
four distinct Pauli terms, so D genuinely moves r. This B difference
uses four original active A axes and lies beyond the support-two
resynthesis.

Let K_A,K_C be actual Clifford encoders taking the sole initialized
clean Z to P. One explicit encoder applies H on wire 1 followed by
CNOTs from wire 1 to wires 2,3,4,5. The complete word is

```math
V=K_C^\dagger E_A^\dagger E_DE_BE_AK_A.
```

It returns the clean wire exactly because the middle word preserves
J. The four native band sizes are [5,2,2,5], with 14 T/T-dagger gates
and depth at most four. There are four arbitrary logical wires, one
clean wire and no additional helper. The independent lists have
Clifford exposures using $`O(q^2)`$ gates on q physical wires;
those exposures and both endpoint Cliffords are paid. The reversed
literal word is the actual full-space inverse, with the same T count
and depth and its own inverse Clifford exposures.

Arbitrary identity spectators and external references may be tensored
with this source. Clean return and center motion do not establish
Boolean occupied action. The following proof excludes it for every
logical Clifford choice of the endpoints, including spectator mixing.

## 4. Exact transport on the actual clean code

Identify the P=+ code by

```math
\begin{aligned}
|y_2,y_3,y_4,y_5\rangle_L
 =\frac{1}{\sqrt2}\bigl(&|0,y_2,y_3,y_4,y_5\rangle\\
 &+|1,1+y_2,1+y_3,1+y_4,1+y_5\rangle\bigr).
\end{aligned}
```

Logical X_j is physical X_j, and logical Z_j is physical Z_1Z_j,
for j from 2 to 5. Introduce the invertible logical coordinates

```math
f=y_2+y_3,\qquad g=y_2+y_4,\qquad
u=y_2+y_3+y_4,\qquad v=y_5.
```

All these sums are binary. This is a Clifford coordinate change for
the analysis, on the explicitly preserved code.

The D pair commutes E_A because its generator XX+YY on wires 1,5
commutes Z_1+Z_5. Its code axes are X_u and -Z_vX_u, giving the
literal operator

```math
D_L=T_{X_u}T_{-Z_vX_u}
 =\begin{cases}
\omega I,&v=0,\\
S_{X_u},&v=1.
\end{cases}
```

For $`C=E_A^\dagger E_BE_A`$, the original center r is logical Z_u.
On u=1 the pair is omega times identity. On u=0 it is S_b, with
$`b=X_1X_2X_3X_4`$. Let w be the ordinary integer sum of
y_2,y_3,y_4. Direct endpoint conjugation gives

```math
E_A^\dagger bE_A|0,y_2,y_3,y_4,y_5\rangle
 =-i^w|1,1+y_2,1+y_3,1+y_4,y_5\rangle.
```

When u=0, w is even and the complementary component has the same
real factor. The induced operator is
$`-(-1)^{w/2}X_v`$. On that sector, w/2 modulo two is the elementary
symmetric quadratic in y_2,y_3,y_4, which becomes f+g+fg. Thus, for

```math
h(f,g)=(1+f)(1+g),
```

the exact C fiber is

```math
C_h
 =\begin{cases}
S_{(-1)^hX_v},&u=0,\\
\omega I,&u=1,
\end{cases}
\qquad
C_h=T_{(-1)^hX_v}T_{(-1)^hZ_uX_v}.
```

The occupied operator is block diagonal in (f,g), with fiber
$`U_h=D_LC_h`$. The value h=1 occurs only at f=g=0. In particular,
the omega factors on inactive sectors have been retained; they are
not discarded as input-dependent phases.

## 5. Its Pauli domain is too small for any diagonal completion

For h=0, the four fiber axes are

```math
b_1=X_v,\quad b_2=Z_uX_v,\quad
d_1=X_u,\quad d_2=-Z_vX_u.
```

Their only cross-band anticommutations are the pairs (b_1,d_2) and
(b_2,d_1). These are two commuting nondegenerate Pauli planes, with
four independent labels. A fixed two-qubit Clifford identifies U_0
with a tensor product of two one-qubit words T_ZT_X, with signed
orientations retained. Such a word has trivial Pauli-preserving
domain: X has a two-term image, while Y and Z have three-term images.
The same holds for either rotation orientation. Tensoring the two
factors therefore gives

```math
M_{U_0}=\{0\}.
```

A general Pauli on the four code qubits has an (f,g) translation a,
an (f,g) Z character, and a Pauli Q on (u,v). If its image under U
is Pauli, each fiber must be proportional to one fixed two-qubit
Pauli. For every a there is an address z such that neither z nor z+a
is zero: two choices remain if a is nonzero, and three if a is zero.
That fiber gives $`U_0QU_0^\dagger`$, so Q must be identity.

If a is nonzero there is also a mixed h=0/h=1 fiber, whose matrix is
$`U_1U_0^\dagger`$ or its adjoint, up to the character phase. But

```math
\begin{aligned}
C_1C_0^\dagger
 &=\begin{cases}iX_v,&u=0,\\I,&u=1,\end{cases}\\
\mathrm{tr}(U_1U_0^\dagger)
 &=\mathrm{tr}(C_1C_0^\dagger)=2.
\end{aligned}
```

A two-qubit Pauli up to scalar has trace magnitude zero or four.
The mixed fiber is therefore impossible. Hence a=0. The surviving
Z characters do commute the entire word, proving exactly

```math
M_U=\mathrm{span}\{Z_f,Z_g\}
   =\mathrm{span}\{Z_2Z_3,Z_2Z_4\}.
```

This space is isotropic of dimension two. A diagonal completion on
the four logical code qubits requires dimension four by Section 1,
so no independent input/output logical Clifford choices can provide
one. With s arbitrary identity spectators, the maximal isotropic
dimension is 2+s while the required dimension is 4+s. Endpoint
Cliffords mixing every spectator cannot remove this deficit.

### 5.1. Commuting D-phase repairs retain the deficit

Keep the same five wires, B band and endpoint bands, and add any
positive J-singletons that can coexist in the actual D list. Their
labels are Z-type, have even weight, and commute X_1X_5. Hence they
lie in

```math
\mathrm{span}\{Z_1Z_5,Z_2Z_3,Z_2Z_4\}.
```

In the (f,g) fibers they add only a phase $`F_{f,g}(Z_v)`$.
Allowing an arbitrary such phase enlarges the native class and makes
the following obstruction stronger; it does not supply a free gate.

In U_0's two-plane basis, $`Z_v=-d_1d_2`$ is the product of the
two D Z axes. Put

```math
D_Z=T_{Z_1}T_{Z_2},\qquad B_X=T_{X_1}T_{X_2}.
```

Any h=0 to h=0 Pauli fiber has arbitrary functions of Z_1Z_2 on the
left and right of $`D_ZB_XQB_X^\dagger D_Z^\dagger`$.
If Q has a nonzero Z mask, conjugation by B_X produces several
distinct X-translation supports. Diagonal factors on either side
preserve these supports and cannot merge or cancel them. Thus a
Pauli image requires $`Q=X^x`$ up to scalar.

For nonzero x, compare the remaining diagonal phase at computational
inputs q and $`q+(1,1)`$. The arbitrary Z_1Z_2 phases cancel from
this ratio. The D_Z contribution has ratio

```math
\begin{cases}
\pm i,&x=(1,0)\text{ or }(0,1),\\
(-1)^{1+q_1+q_2},&x=(1,1).
\end{cases}
```

A Pauli diagonal character has one constant ratio plus or minus one.
Both cases are impossible, forcing Q=I. Every (f,g) translation has
an h=0 to h=0 fiber, so the entire Pauli-preserving domain is contained
in the Pauli space on (f,g). Its maximal isotropic dimension is at
most two, still below four.

Thus no positive commuting D-singleton correction on these five wires
can supply a diagonal completion. Identity spectators do not repair
the deficit. New quantum couplings to spectator axes, modifications
to B, or changes to the endpoint bands are outside this extension.

### 5.2. Every one-band output endpoint on the same five wires

Now retain the original A, B and D bands of Section 3, and allow
any actual one-native-band output endpoint C on these same five
physical wires. Exact clean return requires

```math
E_CP_CE_C^\dagger=J.
```

The middle word preserves J and maps its entire positive eigenspace
onto itself; equality of these half-dimensional clean codes is
equality of their probe involutions. The intrinsic active span of J
has rank five and is exactly the span of the five original Z axes.
Native independence therefore forces C to have exactly five active
axes in that span and no passive axes. Its binary basis matrix O
obeys $`OO^{\mathsf T}=I`$.

Every such five-by-five binary orthogonal matrix is a permutation,
or is a four-coordinate complement-parity block up to row and column
permutations. Here is a short classification. Each row has odd
weight; weight five cannot be orthogonal to another odd row. If all
five rows had weight three, their two-element zero sets would be
distinct and pairwise intersecting. A pairwise-intersecting family
of two-element sets is a star or a triangle, of size at most four
on five coordinates. Thus a weight-one row isolates a column.
In the remaining four coordinates, a further weight-one row forces
all rows to have weight one; otherwise all four rows are the
complements of distinct coordinates.

First give every axis its positive Z representative and positive
orientation, obtaining a band E_0. For a permutation E_0=E_A.
For a complement block on a four-wire set S, put

```math
\begin{aligned}
F_0&=\prod_{i\in S}S_i,\\
H_0&=E_A^\dagger E_0F_0
    =(-1)^{e_1^S+e_2^S+e_3^S}.
\end{aligned}
```

The e_j are ordinary integer elementary-symmetric sums on those four
physical bits, reduced modulo two in the sign H_0. The four
complementary parities have integer
phase exponent $`3e_1^S-4e_2^S+4e_3^S`$ modulo eight.
Subtracting A's exponent and adding F_0's exponent leaves
$`4(e_1^S+e_2^S+e_3^S)`$. This literal Boolean phase is invariant
under complementing all four bits, so it commutes P and

```math
E_0(F_0PF_0^\dagger)E_0^\dagger=J.
```

In the permutation case take F_0=H_0=I.

For general literal signed axes and orientations, write
$`E_C=\alpha L E_0`$, retaining the exact scalar alpha from signed
representative conversion. L is the product of the corresponding
$`S_{\mathrm{axis}}^\dagger`$ factors for negative effective
orientations. All these diagonal factors commute. With
$`F_C=L^\dagger F_0`$,

```math
E_A^\dagger E_CF_C=\alpha H_0.
```

Thus the actual output encoder may be chosen as F_CK_A; every other
encoder of its clean probe differs by a legitimate logical Clifford.
For the displayed encoder the literal occupied source is
$`\overline{\alpha}H_0U`$; the other paid logical endpoint choices
act on its left or right.

On the actual P code, H_0 is computationally diagonal in the logical
y coordinates and hence preserves f,g. In each fixed (f,g) sector
it is a Boolean phase on the two remaining coordinates u,v, so its
fiber is a Clifford. Any ordinary h=0 fiber of H_0U consequently
still has trivial Pauli domain. The controlled-family lemma gives
$`\nu(M_{H_0U})=2\lt4`$. No choice of this actual output endpoint
and no logical Clifford input/output gauges can provide a diagonal
completion.

This is a paid endpoint classification. A complement-parity block
has a nine-CNOT encoder: for one chosen wire in S, apply the three
CNOTs into it, the three out of it, and the three into it again.
Its reversed word supplies the inverse. The band exposure costs
eighteen CNOTs, with any coordinate permutations also paid.
F_0 costs at most four S gates; every factor in L has an ordinary
Pauli-parity Clifford exposure. The actual C list still has five T
gates, so the full source has sizes [5,2,2,5], 14 T gates and at most
four T layers, with no extra helper. Its literal reversed word
retains the common scalar, the actual inverse and exact clean
return. Identity spectators may be added and mixed by the endpoint
Cliffords; new native quantum couplings to them are outside this
five-wire classification.

### 5.3. Common quantum spectator attachments still retain the deficit

There is also a genuine quantum-coupling extension of the first trial.
Attach the same arbitrary signed spectator Pauli P_s to both actual B
axes and the same Q_s to both actual D axes, on s additional wires.
The endpoint bands act as identity on these spectators. The two
attachments may anticommute. Each original pair remains independent,
commuting and J-aligned, with its old center: the common attachment
squares to identity. Complete clean return therefore still holds on
every spectator input, without implying dirty return.

In an ordinary h=0 fiber the four actual induced axes are

```math
\begin{aligned}
B:&\quad X_vP_s,\ Z_uX_vP_s,\\
D:&\quad X_uQ_s,\ -Z_vX_uQ_s.
\end{aligned}
```

Let e be the binary commutation pairing of P_s and Q_s. Their
cross-band pairing matrix is

```math
\begin{pmatrix}e&1+e\\1+e&e\end{pmatrix}.
```

It is invertible for either e. The four labels consequently span a
nondegenerate four-dimensional Pauli space. An analytical Clifford
change of coordinates identifies the fiber with two signed factors
$`T_ZT_X`$ and identity on s remaining wires. Its maximal isotropic
Pauli-domain dimension is s. The general fiber inequality above gives

```math
\nu(M_W)\le2+s\lt4+s.
```

No logical Clifford endpoints, including arbitrary spectator mixing,
can make this complete occupied operator diagonal. This excludes a
specific common-attachment mechanism, not all spectator interactions.
The actual native sizes remain [5,2,2,5], with fourteen T gates, no
added initialized helper, and paid Pauli exposures and actual inverses
on all 5+s wires. Their elementary Clifford cost is $`O((5+s)^2)`$.
All statements concern the full spectator space, without measurement
or a supplied eigenstate. The coordinate change is an analytic test,
not an uncharged replacement of the physical source.

## 6. The first remaining obligation

This source has actual higher-support physical pairs, a moved center,
and exact clean return, and its isolated B factor cannot be made one
native band by Clifford dressing. It nevertheless has no complete
Boolean target even after every output-endpoint refactorization in
Section 5.2. The proof is
analytic; no finite search or numerical fixture is needed.

The next chosen trial changes the actual middle transport or adds a
native quantum coupling outside the common attachments above, and must supply a
logical Lagrangian in the Pauli domain. Boundary Clifford changes,
the commuting phase repairs and the five-wire output-endpoint
refactorizations above cannot do that. Any proposed Boolean source
must also retain the
[paid directional derivative](FOUR_BAND_CENTER_TRANSPORT.md#e-a-paid-derivative-for-a-general-original-physical-pair)
and its indicator center-mask restriction. These tests do not imply
a total-degree bound or an unrestricted depth lower bound. The main
approximate-indicator resource, error and return recurrences remain
unchanged.
