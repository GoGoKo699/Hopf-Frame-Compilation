# An original-pair bound with a bounded middle-singleton coupling rank

[All-active original-pair recognition](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md) ·
[Pair-only enlarged center](PAIR_ONLY_ENLARGED_CENTER.md) ·
[Physical-clean matrix features](OCCUPIED_BOOLEAN_PAIR_CENTER.md#8-a-proved-bounded-span-reopening-with-shared-feature-coordinates)

Positive original middle singletons obstruct the pair-only kernel
argument: a silent-axis product that commutes the pair need not commute
those singletons. The full Boolean endpoint lemma supplies a stronger
replacement. It identifies two entire physical Pauli spans, not just
the kernels generated within individual native lists. Their common
symplectic rank is controlled by the actual singleton-to-silent coupling
rank. The radical is then a genuine transported center.

This proves a source-class degree bound at fixed coupling rank, with
arbitrarily many actual positive singletons and growing silent spans.
It does not bound that rank by the number of original pairs or close
the unrestricted one-pair problem. No product-axis phase is substituted
for a native gate.

## 1. Original-list hypotheses and the bound

Use the actual one-clean three-band word

```math
\begin{aligned}
V&=K_C^\dagger E_C^\dagger E_BE_AK_A,\\
P_A&=K_AZ_cK_A^\dagger,\qquad
P_C=K_CZ_cK_C^\dagger.
\end{aligned}
```

Every phase list consists of independent commuting signed Hermitian
Paulis with literal T/T-dagger orientations. Assume exact clean return
and complete occupied action

```math
G=e^{i\theta}\,\mathrm{diag}((-1)^{g(x)})
```

on the entire logical and dirty register, including arbitrary reference
entanglement. The common scalar is retained; the unused clean block is
unrestricted. Suppose the original middle list has exactly one genuine
physical pair B,D, and its difference $`r=b+d`$ is the sum of all
actual clean-active axes at each endpoint. Arbitrarily many original
singleton axes are allowed. These are support conditions on the
physical lists, not on a freely chosen basis after quotienting.

Let $`U_A`$ be the original input-silent span and T the span of the
original middle singletons. Section 2 proves that all those singletons
are positive. Define the actual bilinear coupling rank

```math
h=\mathrm{rank}\bigl(\langle\cdot,\cdot\rangle:
                         T\times U_A\longrightarrow\mathbb F_2\bigr).
```

The corrected output-silent span has the same rank against T. For
nonnegative integer k set

```math
\begin{aligned}
A_k&=12\,2^{k+1}+4(k+1)+3,\\
q_k&=2k+4,\qquad N_k=2^{2k^2+5k+10},\\
B_k&=9\,2^k+2k+3,\\
D_k&=48q_kN_k+4B_k+k.
\end{aligned}
```

Then the complete occupied bounds are

```math
\begin{aligned}
\deg g&\le D_{2h+4}+4A_{2h+4},\\
\deg f&\le D_{2h+4}+4A_{2h+4}-1
       &&\text{for a Boolean XOR query }Q_f.
\end{aligned}
```

The sharper parameter $`2\rho+2`$ may replace $`2h+4`$, where the
augmented actual coupling rank $`\rho\in\{h,h+1\}`$ is defined
in Section 4. No bound on h is inferred from the single pair.

Retain sharper results in their respective subclasses. A physical
shared clean flip gives sign/XOR degree five/four. If either endpoint's
silent list commutes the entire original middle list, the
[one-sided theorem](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md#1-statement-in-the-original-physical-lists)
gives 2428/2427. With no middle singletons, h is zero and the displayed
bound is exactly the retained
[pair-only constant](PAIR_ONLY_ENLARGED_CENTER.md)
$`D_4+4A_4=576\cdot2^{62}+2252`$. The existing rank-two endpoint
bounds remain sharper where applicable.

## 2. Literal repair and the common active span

Put $`J=E_AP_AE_A^\dagger`$. For the signed original pair relation
$`JBJ=sD`$ and original orientations epsilon,delta, use

```math
F=T_B^\epsilon T_D^{s\epsilon},\qquad
C=T_D^{\delta-s\epsilon},\qquad E_B=CFE_0.
```

C is an actual Clifford, F commutes J, and $`E_0`$ retains all
original singleton axes and orientations. Absorb C into the output
encoder and conjugate the actual output phase list as in the
[original insertion identity](ORIGINAL_MIDDLE_PAIR_INSERTIONS.md#1-keep-the-actual-middle-axes-and-their-signs).
Use corrected C-side notation below. C fixes the pair difference, so
the output all-active sum remains r, and all lists stay independent.

The [all-active clean comparison](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md#2-orientation-repair-removes-every-negative-singleton)
then removes every negative singleton: the two displacement quadratic
forms give $`r=r+\sum_{u\in B_-}u`$, and original middle
independence forces the negative list to be empty. Hence

```math
E_CP_CE_C^\dagger=J,\qquad L_A=L_C=L,
\qquad p_A+p_C\in L.
```

The complete nonzero Pauli supports, not just their spans, establish
the last relation. Let $`U_C`$ be the corrected output-silent span.
The all-active condition also gives

```math
\langle b,l\rangle=\langle p_A,l\rangle
\quad(l\in L),\qquad R=L\cap p_A^\perp.
```

R is the genuine transported center of the occupied-center theorem.
Both silent spans commute L and $`P_A`$. T is isotropic and commutes
$`L,b,P_A`$: its axes are original positive singletons. These are
statements about the actual signed Pauli operators as well as labels.
The entire corrected phase-axis span is

```math
M=L+\langle b\rangle+T+U_A+U_C.
```

## 3. Boolean endpoint transport identifies entire physical spans

For any physical Pauli label q, commuting with $`P_A`$ and every
actual A and B axis is exactly membership in

```math
\langle p_A,L,b,T,U_A\rangle^\perp.
```

The [Boolean endpoint commutation lemma](OCCUPIED_BOOLEAN_PAIR_CENTER.md#1-actual-source-and-the-endpoint-commutation-lemma)
sends this domain into
$`\langle p_C,L,b,T,U_C\rangle^\perp`$. Applying the same lemma
to the actual reversed circuit gives the reverse inclusion. Taking
symplectic orthogonal complements and using $`p_A+p_C\in L`$
therefore proves

```math
\begin{aligned}
\widetilde M
 &=\langle p_A,L,b,T,U_A\rangle\\
 &=\langle p_A,L,b,T,U_C\rangle\\
 &=M+\langle p_A\rangle.
\end{aligned}
```

This uses arbitrary physical Pauli products in the endpoint lemma; it
does not assign T phases to a new spanning set.

In particular, a member of T commuting with $`U_A`$ commutes both
clean Paulis and all A and B axes. It therefore commutes $`U_C`$,
and conversely. The two coupling maps on T have the same kernel,
proving equality of their ranks h.

If $`p_A\notin M`$, symplectic nondegeneracy supplies a label eta
orthogonal to M with $`\langle\eta,p_A\rangle=1`$. Its Hermitian
Pauli also flips $`P_C`$, since $`p_A+p_C\in L`$. It commutes
every actual phase axis. The retained
[clean-controlled Pauli dressings](RANK_ONE_CLEAN_PHASE_INTERFACE.md#23-a-full-odd-coset-clean-flip-criterion-and-its-dual)
fix every clean-zero vector and give the full-space operator
$`I_c\otimes G`$, with the same three native bands and charged
Clifford dressings. The zero-clean theorem gives five/four. The rest
of the proof treats $`p_A\in M`$, so $`\widetilde M=M`$.

## 4. The exact symplectic rank with the singleton band retained

Choose an actual A-active label a and define

```math
\begin{aligned}
\beta&=\langle p_A,b\rangle,\qquad
q=p_A+b+\beta a,\\
\alpha(u)&=\langle u,b\rangle\quad(u\in U_A),\\
\gamma(v)&=\langle v,b\rangle\quad(v\in U_C),\\
u'&=u+\alpha(u)a,\qquad v'=v+\gamma(v)a.
\end{aligned}
```

Let $`U_A',U_C'`$ be the images of these linear shifts, and put
$`S_0=L+\langle b\rangle`$ and $`B_0=T+\langle q\rangle`$.
Both shifted silent spaces and $`B_0`$ commute $`S_0`$. Each is
isotropic. For example $`q`$ commutes L because the pairings of
$`p_A`$ and b with L coincide, and it commutes b because the two
beta terms cancel. A shifted u commutes b because its two alpha
terms cancel. Silent axes commute a, so the shifts preserve isotropy.

The span identity of Section 3 implies

```math
\begin{aligned}
W&=R+B_0+U_A'=R+B_0+U_C',\\
M&=S_0+W,\qquad [S_0,W]=0,\qquad S_0\cap W=R.
\end{aligned}
```

For completeness, the two displayed descriptions of W follow because
both centralize $`S_0`$ and generate the same M together with it.
Their possible difference lies in
$`\mathrm{rad}(S_0)=R`$, already included in each W. The quotient
$`S_0/R`$ is a nondegenerate symplectic plane.

Define the augmented coupling rank

```math
\rho=\mathrm{rank}\bigl(\langle\cdot,\cdot\rangle:
                         B_0\times U_A'\longrightarrow\mathbb F_2\bigr).
```

Pairing a shifted u with T gives its original singleton couplings;
pairing it with q gives $`\alpha(u)`$. Thus this rank is obtained
by adjoining the one pair-coupling row to the original singleton
coupling map, and

```math
h\le\rho\le h+1.
```

The same rank is obtained from $`U_C'`$: a $`B_0`$ element
commutes one shifted silent space exactly when it commutes W, since
$`B_0`$ and R are isotropic and centralize $`B_0`$. The two
coupling maps on $`B_0`$ therefore have the same kernel. In
particular, the same $`t\in T`$ satisfies
$`\langle t,u\rangle=\alpha(u)`$ for every u exactly when it
satisfies $`\langle t,v\rangle=\gamma(v)`$ for every v.

W is the sum of two isotropic spaces, together with the central R.
Its form rank is twice their cross-pairing rank, even if their
spanning lists overlap: the Gram matrix on their combined generators
has blocks zero, A, $`A^{\mathsf T}`$, zero. Consequently

```math
\mathrm{rank}_{\rm symp}M=2+2\rho\le2h+4.
```

Take the actual physical radical

```math
Z=\mathrm{rad}(M).
```

Every member of Z commutes every actual phase axis. Since both clean
Paulis lie in M, it commutes those too. Thus its signed Pauli operators
have full physical transport between the endpoint encoders. Z is
isotropic. Its projection to the positive clean code is injective:
$`p_A\notin Z`$, witnessed by its anticommutation with the actual
active axis a in M. The output statement follows identically. Hence

```math
\dim(M/Z)=2+2\rho\le2h+4.
```

This converts a form-rank calculation into a quotient by a newly proved
transport center. It does not bound the old coupled-silent span t or
declare that bounded pairing rank bounds a physical axis span.

## 5. The complete physical-clean feature proof

Use the actual active localization

```math
\begin{aligned}
C_0&\in\mathrm{Alg}(S_0),\qquad C_0P_AC_0^\dagger=P_C,\\
G_0&=C_0^\dagger E_{C,\rm active}^\dagger F E_{A,\rm active},\\
K_C'&=C_0^\dagger K_C,\qquad
E_{C,s}'=C_0^\dagger E_{C,s}C_0.
\end{aligned}
```

The exact physical word is

```math
V=(K_C')^\dagger(E_{C,s}')^\dagger
             G_0E_{B,+}E_{A,s}K_A.
```

Every internal factor preserves $`P_A`$. Because $`S_0\subseteq M`$,
the actual Clifford $`C_0`$ fixes Z literally and preserves M. The
conjugated output list remains independent and lies in M. The positive
middle list is retained with every original phase and orientation;
it is not replaced by the generators $`B_0`$ or shifted silent labels.

Apply the retained global Boolean Pauli-transport normalization to Z:
extract the quadratic Boolean Clifford, make fixed affine coordinate
changes, and use Hadamards only on factors where the residual target
is already identity. Z becomes computational Z at both endpoints.
The logical X-mask map of the output-conjugated physical M annihilates
Z, so its image has rank

```math
k\le2\rho+2\le2h+4.
```

One fixed affine change confines these X masks to k logical core bits.
Let $`\mathcal M`$ be the computational-Z subgroup on all remaining
logical bits. Every output-conjugated physical source axis in M
commutes $`\mathcal M`$, as does the actual $`C_0`$.

Keep the physical clean bit as one additional core coordinate. The
[enlarged-center feature proof](PAIR_ONLY_ENLARGED_CENTER.md#6-keep-the-physical-clean-coordinate-through-the-feature-proof)
now applies with the positive singleton band retained. In detail,
let $`E_*`$ contain the actual transformed C-silent axes coupling to
$`S_0\cap p_A^\perp`$. Cancel only passive phases in the dressing
against their actual inverses. The complete logical factors are

```math
\begin{aligned}
F_{\rm out}&=(K_C')^\dagger E_*^\dagger G_0E_*K_C',\\
W_3&=(K_C')^\dagger(E_{C,s}')^\dagger
                         E_{B,+}E_{A,s}K_A,\\
G&=F_{\rm out}W_3.
\end{aligned}
```

These restrictions occur only after the full factors preserve the
clean code. All three separate lists in $`W_3`$ remain independent.
For the positive middle list, any actual A-active axis commutes its
span and anticommutes $`P_A`$, ruling out a relation created by
restriction. The input-silent and transformed output-silent lists
have their respective retained active witnesses. A common fixed
clean coordinate makes both restricted endpoints actual logical
Cliffords. No dependent union is treated as a native band.

The four distinct physical phase bands used by $`F_{\rm out}`$
all lie in M; its source ledger is therefore the same $`A_k`$
degree-at-most-four features as in the retained proof. The clean
coordinate remains present during matrix multiplication. On the
remaining logical coordinates, $`F_{\rm out}`$ commutes
$`\mathcal M`$. Since the residual G is diagonal,
$`W_3=F_{\rm out}^\dagger G`$ fixes $`\mathcal M`$ pointwise on
the entire logical space. The
[partial-diagonal theorem](PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#1-statement-and-explicit-constants)
applies to this actual three-native-band word, including all positive
singleton gates. Its retained matrix-feature bounds contribute $`D_k`$
to the resulting degree bound.

Only the product is assumed Boolean. Combining the two full matrix
feature lists and at most k core input bits gives

```math
\deg g\le D_k+4A_k
       \le D_{2\rho+2}+4A_{2\rho+2}
       \le D_{2h+4}+4A_{2h+4}.
```

The shared normalization returns this bound to the original input
bits: its residual argument changes are affine, and the extracted
phase is quadratic. This is a complete sign bound, not a separate
bound on independently selected central sectors.

## 6. The remaining matching block is a symmetric coupling form

The span argument also narrows what is left when h is unrestricted.
In the symplectic quotient $`W/\mathrm{rad}(W)`$, the images of
$`B_0,U_A',U_C'`$ have dimension rho. The two silent images are
Lagrangian complements of $`B_0`$. Choose matched coupling coordinates
$`\lambda(u),\mu(v)\in\mathbb F_2^\rho`$ for original silent
labels. Their actual cross pairing has the form

```math
\langle u,v\rangle
 =\lambda(u)^{\mathsf T}H\mu(v),\qquad H=H^{\mathsf T}.
```

Indeed choose $`u_i',v_i'`$ with the same coupling functional on
$`B_0`$, and choose dual labels $`b_j\in B_0`$. Modulo the radical,

```math
v_i'=u_i'+\sum_jH_{ji}b_j.
```

Isotropy of $`U_C'`$ makes H symmetric. The shifts by a preserve
the original cross pairings because both silent spans commute a.
With no singleton span, rho is at most one and this reduces to the
rank-one pair-only formula.

These are identities among actual Pauli labels and products. The
independent native T-phase lists retain their original axes. A basis
change, a symmetric coupling form, or a product in $`B_0`$ does not
provide a free phase band or a free Clifford correction to the word.
When rho grows, the unpaid step is matching or eliminating this
interaction while retaining every original phase and a complete
shared feature ledger. No argument here bounds h from one pair, and
the form alone is not an asserted realizable Boolean counterfamily.

## 7. Complete contract, inverse and approximation

The original source and its actual inverse retain their T counts,
at most three native T bands, one clean wire and their entire original
dirty allocation. Orientation repair, $`C_0`$, the fixed global
coordinate words and their actual inverses carry their
$`O((Q+1)^2)`$ elementary-Clifford charges on Q logical wires and
the existing clean wire. The shared-flip branch uses its actual
$`O(Q+1)`$ Clifford dressings and no helper. No emitted phase list
is replaced by the shifted generators in Section 4.

Every reduction is exact. The unused clean sector is retained until
legitimate full logical restriction, and all logical, dirty and
reference inputs obey the target's complete return contract. The
actual inverse reverses the original orientations and coordinate
words, with the same resource majorants and zero added error.
The matrix-feature description supplies no new synthesis upper bound.

Output Hadamards turn a complete XOR query into the Boolean sign
$`(-1)^{y f(x)}`$ without changing its clean Paulis, physical lists,
or coupling ranks. Its degree is one greater for nonzero f, proving
the stated XOR bound.

For approximate Boolean action, apply the
[original-circuit rounding theorem](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
before any restriction or coordinate reduction. The retained
conservative floor $`2^{-530690}`$ excludes higher-degree targets
in this fixed-h class, allowing one common scalar. It is not an
arbitrary fixed-tolerance obstruction.

This extends the pair-only result to bounded actual middle coupling,
including arbitrary numbers of uncoupled positive singletons when
h is zero. The older coupled-silent span may still grow relative to
its smaller center. General original pairs, unrestricted middle
coupling and the dimension-growing depth question remain open; the
general frame resource frontier is unchanged.
