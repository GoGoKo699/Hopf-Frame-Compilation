# Three native bands with at most one-dimensional overlap

[Native band structure](NATIVE_BAND_SOURCE_STRUCTURE.md) ·
[Signed Clifford transport](NATIVE_CLIFFORD_TRANSPORT.md) ·
[Six-layer reduction](READ_ONLY_SIX_LAYER_REDUCTION.md)

An exact third-level source with three physical T bands is diagonal
third-level in fixed Clifford coordinates when its first two active-axis
spans intersect in dimension at most one. The proof closes both exceptional
one-relation cases. It retains a genuine order-three intermediate action;
no nilpotence or intermediate Clifford-membership premise is introduced.
A separate extension lemma also handles a commuting T-gadget list with one
binary dependency, including its four-axis exception.

## 1. Exact statement and native scope

In common fixed Clifford coordinates write the specified source as

```math
F=E_CE_BE_A.
```

Each E is a Clifford conjugate of one physical T/T-dagger band. Its active
binary axis list is independent and isotropic. Let A and B denote the
first two lists, with spans $`U_A,U_B`$. Axis signs, inactive wires, native
scalar phases and actual adjoints remain in the specified operators.

**One-overlap theorem.** Suppose F is C3 on the full work space and

```math
t=\dim(U_A\cap U_B)\le1.
```

Then F is diagonal C3 in fixed Clifford coordinates: it has the form
$`C_{\rm out}D C_{\rm in}`$ with fixed Cliffords on the two sides and D a
weighted cubic phase modulo eight, as defined in the
[native source theorem](NATIVE_BAND_SOURCE_STRUCTURE.md#1-common-native-coordinates-and-exact-hypotheses).
The end Cliffords need not be inverses.

The t=0 case is already proved there. For t=1, let r be the unique nonzero
coefficient relation expressing a B combination in $`U_A`$. The earlier
overlap proof identifies the activation code of Paulis centralizing A as
$`r^\perp`$. If r has weight at least three, every coordinate functional
on this code is nonzero and distinct. The proof therefore makes every B
axis commute with every C axis, and the existing outer-two-band certificate
applies. Only relation weights one and two remain; sections 3 and 4 close
them.

All conclusions here require exact C3 membership of the specified source.
No robust t=1 source certificate is asserted. Complete returned-query action
does not itself imply C3 membership or the overlap bound. A different native
resynthesis can change the axis lists and their overlap.

## 2. Paired half-angle rigidity

**Lemma.** Let $`H_1,H_2`$ be commuting Hermitian Clifford involutions with

```math
\mathrm{tr}H_1=\mathrm{tr}H_2=\mathrm{tr}(H_1H_2)=0,
```

where tr is normalized trace. If

```math
K=\exp(i\pi H_1/4)\exp(i\pi H_2/4)
 =\tfrac12(I+iH_1+iH_2-H_1H_2)
```

is projectively Clifford, then either:

1. Both $`H_1`$ and $`H_2`$ are signed Paulis; or
2. Both have the same four-element Pauli support, each coefficient is plus or minus 1/2, and $`H_1+H_2`$ is a sum of two signed commuting Paulis.

In both cases $`H_1H_2`$ is a signed Pauli. The conclusion also applies with either exponent orientation reversed, by replacing the corresponding H with -H.

**Proof.** The normalized trace of K is exactly 1/2. Flat Pauli coefficients of a Clifford imply that K has exactly four nonzero Pauli coefficients, each of magnitude 1/2. Its identity coefficient is real. Taking adjoints conjugates the coefficients in the Hermitian Pauli basis and preserves their support. Hence

```math
H_1+H_2=(K-K^\dagger)/i
```

has at most three nonzero Pauli coefficients, none on identity.

Let $`N_j`$ be the Pauli support size of $`H_j`$. Each $`N_j`$ is a power of two and its nonzero coefficients are real, of magnitude $`N_j^{-1/2}`$.

If $`N_1`$ and $`N_2`$ differ, no coefficient on their intersection cancels. The union of their supports has size at most three, so each $`N_j`$ is at most two and one equals one. That $`H_j`$ is Pauli, so its half-angle rotation is Clifford. Cancelling it from K and applying the single half-angle lemma forces the other $`H_j`$ to be Pauli too. Thus unequal support sizes are in fact impossible.

If $`N_1=N_2=N`$, orthogonality gives normalized Hilbert--Schmidt norm squared $`\|H_1+H_2\|_2^2=2`$. Each coefficient has magnitude at most $`2/\sqrt N`$, and at most three survive. Thus $`2\le12/N`$, so $`N\le6`$. Since N is a power of two, N is one, two or four.

$`N=1`$ is the Pauli case. For $`N=2`$, the common-support size must be even by orthogonality, since each common coefficient product is plus or minus 1/2. If there are no common terms, the sum has four terms, impossible. If both supports coincide, orthogonality leaves one cancelled term and one surviving term of magnitude $`\sqrt2`$. Thus $`(H_1+H_2)^2=2I`$, contradicting commutation, which gives $`2I+2H_1H_2`$. Hence $`N=2`$ is impossible.

For $`N=4`$, let m be the intersection size. Orthogonality forces m even with exactly m/2 cancellations. The sum support has size 8-m-m/2. A bound of three therefore forces m=4. There are precisely two surviving terms, each of magnitude one. Write $`H_1+H_2=P+Q`$ for signed Hermitian Paulis P,Q. If P,Q anticommuted, their sum squared would be $`2I`$, giving the same contradiction. Therefore they commute and

```math
H_1H_2=PQ
```

is Pauli. This proves every assertion. QED.

## 3. A weight-one relation

### 3.1 A support-coset consequence for a native image

Suppose $`U=E_CE_G`$, where G is an isotropic band whose axes commute with every C axis. Let b be a Pauli axis commuting with $`E_G`$. Then

```math
H_b=U b U^\dagger=E_C b E_C^\dagger.
```

If exactly m C axes anticommute with b, its Pauli support is the affine coset $`b+\mathrm{span}(C_{\mathrm{bad}})`$, of size $`2^m`$, with all coefficient magnitudes $`2^{-m/2}`$. Native axis independence prevents collisions.

Suppose a is a distinct Pauli axis commuting with b, and $`H_a=UaU^\dagger`$ is Clifford. If the paired half-angle criterion holds for $`H_a,H_b`$ and $`H_b`$ is non-Pauli, the preceding lemma forces m=2 and gives $`H_a,H_b`$ the same four-element support. Write the two $`C_{\mathrm{bad}}`$ axes as $`c_1,c_2`$. Their product $`R=H_aH_b`$ is a nonidentity Pauli (its trace is zero). Multiplication of equal coset supports puts the label of R in $`\mathrm{span}\{c_1,c_2\}`$. R commutes with $`H_b`$. Since U commutes both $`c_1,c_2`$, commutation of any Pauli in that span with $`H_b`$ is identical to its commutation with b. Each $`c_j`$ anticommutes with b, so the only nonzero commuting label in their span is $`c_1+c_2`$. Therefore

```math
\mathrm{label}(H_aH_b)=c_1+c_2,
\qquad a+b=c_1+c_2.
```

The second equality follows by conjugating the Pauli R back through U, which fixes every Pauli in $`\mathrm{span}\{c_1,c_2\}`$ literally. Actual Hermitian signs may change the signed equality, but cannot change these label equalities.

### 3.2 Closing the source case

**Theorem.** Let $`F=E_CE_BE_A`$ be C3. Suppose the one-dimensional space $`U_A\cap U_B`$ is generated by one individual B axis b. Then F is diagonal C3 in fixed Clifford coordinates.

The established overlap proof, with the unique relation supported on b, says that every other B axis commutes with every C axis. Split $`E_B=E_GE_b`$. All B axes commute, and b belongs to the isotropic span $`U_A`$, so $`E_b`$ commutes with $`E_A`$. Set $`U=E_CE_G`$. It is diagonal C3, and

```math
F=U E_A E_b.
```

Write the binary label b uniquely as $`\sum_{i\in J}a_i`$. J is nonempty because b is nonzero.

If $`|J|=1`$, the corresponding actual signed native A rotation and $`E_b`$ are rotations about the same binary axis. Their product is a scalar times a Clifford (their signed $`\pi/8`$ angles sum to an even multiple of $`\pi/8`$), commuting with the remaining A band. Thus $`F=UE_{A\setminus J}C_b`$ for an actual Clifford $`C_b`$, with the scalar retained. Right multiplication by its actual inverse shows $`UE_{A\setminus J}`$ is C3. The diagonal native right-extension corollary makes it diagonal C3, and restoring $`C_b`$ on the right preserves that source class.

Now suppose $`|J|\ge2`$. For any i in J choose a Pauli P whose commutation vector against A is the unit vector at i. Then it also anticommutes with b, since b is the sum of those A labels. The exact image of P through $`E_AE_b`$ is P times the product of the two Clifford half-angle rotations about the actual signed axes $`a_i`$ and b, with orientations determined by their native T/T-dagger choices. F is C3 and U is C3, so cancelling the Clifford $`UPU^\dagger`$ proves the paired half-angle criterion for $`H_i=Ua_iU^\dagger`$ and $`H_b=UbU^\dagger`$. Both are Hermitian Clifford involutions; their individual traces and product trace vanish because $`a_i,b`$ are distinct nonidentity commuting Paulis.

If $`H_b`$ is Pauli, its native image through $`E_C`$ has spread zero, so b commutes with every C axis. All B axes then commute with C, and the established outer-two-band certificate proves the theorem.

If $`H_b`$ is non-Pauli, the support-coset consequence applies to every i in J. It gives two fixed C axes $`c_1,c_2`$ and

```math
a_i+b=c_1+c_2\quad\text{for every }i\in J.
```

This makes all selected A labels identical, contradicting independence when $`|J|\ge2`$. The non-Pauli case is impossible, completing the theorem. QED.

This proof makes no assumption that $`F^\dagger`$ is C3 and no assertion that the inner two-band source is C3. The given unitary and actual source factorization remain fixed. It introduces no free physical coordinate transformation, source, inverse, auxiliary qubit, or prepared-state premise. The all-column native source equality and original resource accounting are unchanged.

## 4. A weight-two relation

Assume the unique B-coordinate relation has weight two and some B axis
fails to commute with C. If all B/C axes commute, the retained outer-two-band
certificate already proves the conclusion.

### 4.1. There is a single bad pair in each outer band

Relabel the two exceptional B axes as $`b_1,b_2`$, and put
$`r=b_1+b_2\in U_A`$. The activation code of Paulis centralizing A is
$`u_1=u_2`$, with all other B coordinates free. The overlap argument gives
that all other B axes commute with C. A C axis can anticommute with neither
or both of $`b_1,b_2`$: anticommutation with exactly one would produce a
nontrivial transvection permutation of an isotropic C set. If c anticommutes
with both, the activation $`u_1=u_2=1`$ sends c to $`c+r`$. Thus its orbit
is the pair $`\{c,c+r\}`$. Two distinct such pairs would sum to zero,
contradicting independence of C. Consequently there is exactly one bad C pair.

Choose representatives b and c from the bad B and C pairs. Then

```math
W=\mathrm{span}\{r,b,c\},\qquad
\langle b,c\rangle=1,\qquad
r\in\mathrm{rad}(W).
```

Every remaining B or C axis commutes with W, and all remaining B/C axes
commute with one another. Their product $`G=E_{C,\mathrm{good}}E_{B,\mathrm{good}}`$
is a commuting diagonal-C3 phase source. The outer two bands are
$`G E_{C,\mathrm{bad}}E_{B,\mathrm{bad}}`$.

### 4.2. Isolate one A coordinate of the relation

Write $`r=\sum_{j\in J}a_j`$ in the independent A list. Fix any
$`j\in J`$. Choose a Pauli p with $`\langle p,a_k\rangle=\delta_{jk}`$.
Then $`\langle p,r\rangle=1`$. Use r and p as the Z and X axes of a
control qubit. Within each bad pair choose the representative commuting with
p; this is possible and unique because its two members differ by r. The
chosen b,c anticommute and commute with both control axes, so form the two
axes of a logical target qubit. Complete to a symplectic basis by spectator
qubits. This changes representatives within W, and does not change W itself.

All A axes commute with r. Every A axis other than $`a_j`$ also commutes
with p, whereas $`a_j`$ has a control-Z component. In the two control sectors,
the actual native phases therefore obey

```math
E_{A,+}E_{A,-}^{\dagger}
=\exp(\pm i\pi\tau(a_*)/4)
```

up to its retained scalar, where $`a_j=r+a_*`$ in these coordinates and
$`a_*`$ has no control component. The case $`a_*=0`$ means a scalar relative
phase. Only this isolated native phase survives; all other commuting A
factors cancel literally.

Because G commutes with r,b,c, its sector restrictions $`G_+,G_-`$ act only
on spectators. They are commuting diagonal C3 phases in common spectator
coordinates. Their ratio is a diagonal Clifford: each native good phase
contributes either the identity or a Pauli quarter-turn. Thus

```math
D=G_+G_-^{\dagger},\qquad S_D^2=I.
```

The bad outer core restricts in each sector to an exact one-qubit Clifford
$`R_+,R_-`$, retaining scalar powers of $`\omega=e^{i\pi/4}`$. Each bad
paired band contributes a Pauli quarter-turn in exactly one sector and a
scalar in the other. Explicitly, put

```math
N_e(P)=\frac{I+P}{2}+\omega^e\frac{I-P}{2},\qquad e\in\{1,-1\}.
```

Axis signs are retained using $`N_e(-P)=\omega^eN_{-e}(P)`$. After this
identity puts a pair in the form $`N_e(P)N_f(rP)`$, its exact restrictions
are

```math
N_e(P)N_f(rP)\big|_{r=+1}=N_{e+f}(P),\qquad
N_e(P)N_f(rP)\big|_{r=-1}=\omega^fN_{e-f}(P).
```

Exactly one of $`e+f,e-f`$ is zero and the other is $`2`$ or $`-2`$.
Thus all these block operators are ordinary one-qubit Cliffords with
eighth-root scalars. The same identity makes the ratio of the two sectors
of each good phase a Clifford. Hence the binary symplectic action

```math
M=S_{R_-}^{-1}S_{R_+}
```

has order three on the target plane: it is the product, in one order or the
other, of the two distinct transvections along b and c. Indeed, if
$`B_\pm,C_\pm`$ denote the binary branch actions, then
$`B_-B_+=T_b`$, $`C_-C_+=T_c`$, and each branch action is either the
indicated transvection or the identity. Therefore

```math
M=B_-C_-C_+B_+=B_-T_cB_+\in\{T_bT_c,T_cT_b\}.
```

### 4.3. Full C3 membership removes every spectator component

The source commutes with r, so write its blocks as $`F_+,F_-`$. The image of
the control X is the block operator with off-diagonal blocks
$`V=F_+F_-^{\dagger}`$ and $`V^{\dagger}`$. Since F is C3, that operator is
Clifford. Conjugating spectator/target Paulis by it shows that V is itself
Clifford, and that V and its inverse have the same binary Pauli action:

```math
S_V^2=I.
```

On the other hand,

```math
V=R_+\bigl(G_+\exp(\pm i\pi\tau(a_*)/4)G_+^{\dagger}\bigr)
 R_-^{\dagger}D.
```

All factors except the parenthesized one are Clifford. For nonzero $`a_*`$,
the half-angle lemma therefore makes
$`H=G_+\tau(a_*)G_+^{\dagger}`$ a signed Pauli. Write its label as q+s,
where q is on the target plane and s is on spectators. After conjugating
$`S_V`$ by $`S_{R_-}^{-1}`$, its binary symplectic action is

```math
S=M T_{q+s}S_D.
```

Here M acts only on the target, $`S_D`$ only on spectators, and
$`T_v=I+v\ell_v`$, with $`\ell_v(x)=\langle v,x\rangle`$. If q is zero,
the target diagonal block of $`S^2`$ is $`M^2\ne I`$, impossible. This
also excludes the scalar case $`a_*=0`$.

For nonzero q, the one-qubit identities

```math
A=M T_q,\qquad A^2=I,\qquad A Mq=q
```

hold. The block matrix of S on target plus spectators is

```math
S=
\begin{pmatrix}
A & Mq\,\ell_s S_D\\
s\,\ell_q & T_s S_D
\end{pmatrix}.
```

Put $`\alpha=\langle s,S_Ds\rangle`$. The target diagonal block of
$`S^2=I`$ gives $`\alpha=0`$, because it is
$`I+Mq\,\alpha\ell_q`$. Its upper-right block then gives

```math
q\,\ell_s S_D+Mq\,\ell_s=0.
```

The vectors q and Mq are independent, since M has order three on a binary
plane. Hence $`\ell_s=0`$ and s is zero. Since $`G_+`$ acts only on
spectators, it cannot turn a nonidentity spectator factor into the identity.
Thus $`a_*`$ itself has no spectator factor, and $`a_j\in W`$.

Repeat for every $`j\in J`$, choosing a fresh isolating p each time. The
space W is unchanged by these choices. Its maximal isotropic dimension is
two, so independence and isotropy of A give $`|J|\le2`$. The case
$`|J|=1`$ would mean $`a_j=r`$, already excluded. Therefore A contains
exactly the pair

```math
\{a,a+r\}\subset W,
```

whose sum is r. No assertion is made that the remaining A axes commute with
all of W.

### 4.4. Retain phases while extracting the core

Split the actual A band into this pair and its remaining independent axes:
$`E_A=E_{A,\mathrm{pair}}E_{A,\mathrm{rest}}`$. Define

```math
K=E_{C,\mathrm{bad}}E_{B,\mathrm{bad}}E_{A,\mathrm{pair}},\qquad
F=(GK)E_{A,\mathrm{rest}}.
```

K acts only on the two-qubit core and preserves r. Its two control blocks
$`K_0,K_1`$ are exact one-qubit Cliffords with scalar powers of $`\omega`$.
Keep the isolating p for one member of the A pair. Every remaining A axis
commutes with p and r, so its two sector restrictions are identical. They
cancel in $`E_{A,+}E_{A,-}^{\dagger}`$. Consequently
$`V=(K_+K_-^{\dagger})D`$, with its two factors acting respectively on
target and spectators. The preceding involution condition therefore says
that $`K_+K_-^{\dagger}`$, and equivalently $`K_0^{\dagger}K_1`$, has
involutory binary symplectic action. Such an action fixes a nonzero target Pauli label. Choose
a fixed target Clifford Q carrying Z to that label. The resulting relative
Clifford is monomial, hence has the exact form

```math
Q^{\dagger}K_0^{\dagger}K_1Q
=\omega^{\ell}X^eS^k.
```

The scalar is an eighth root, with no arbitrary scalar discarded. Every
native paired branch has determinant in $`\mu_4`$, hence so does
$`K_0^\dagger K_1`$. A monomial Clifford $`\gamma X^eS^k`$ has determinant
$`\gamma^2(-1)^e i^k`$, forcing $`\gamma\in\mu_8`$. The exact factorization is

```math
K=(I\otimes K_0Q)\,\mathrm{C}_r(X^e)\,
 \mathrm{diag}\bigl(\omega^{\ell z_r+2kz_rz_t}\bigr)
 (I\otimes Q^\dagger).
```

Controlled-X is Clifford, and the remaining diagonal phase has the permitted
weighted quadratic coefficient modulo eight.

Every extracted fixed Clifford acts only on the core and preserves r, so it
commutes with G. The diagonal phase above and G have common diagonal
coordinates consisting of r, the chosen target Z, and commuting spectator
axes. Thus GK is diagonal C3 in fixed Clifford coordinates with literal
phase retained. Finally the retained native right-extension corollary
applies to $`(GK)E_{A,\mathrm{rest}}`$: the latter is C3 by hypothesis and
$`E_{A,\mathrm{rest}}`$ remains an actual independent native band. Therefore
F is diagonal C3 in fixed Clifford coordinates.

This proves the weight-two case. Along with the weight-one argument and the
weight-at-least-three overlap reduction, it gives the t=1 source conclusion;
it makes no assertion for t at least two or for arbitrary returned sources.

## 5. Six-layer consequence and a literal noncommuting fixture

For the actual zero-clean read-only six-layer comparison, suppose its
specified shared three-band prefix is C3 and has the overlap $`t\le1`$
above. The theorem verifies the source hypothesis of the retained
[signed transport theorem](NATIVE_CLIFFORD_TRANSPORT.md#5-consequence-for-the-actual-six-layer-query).
Consequently, under the original complete-query error

```math
\eta\lt\sqrt2/8,
```

the Boolean output has degree at most **78**, or **26** when the outer
target shear is visible. The error radius and endpoint-degree calculation
are inherited from that theorem; this source result does not supply a
rounding theorem for general non-C3 prefixes.

A two-qubit fixture shows that the new source condition permits both
adjacent two-band partial sources to be outside C3. Let r be $`Z_0`$, with
target qubit 1, and use positive native phases on

```math
B=\{Z_1,Z_0Z_1\},\qquad A=C=\{X_1,Z_0X_1\}.
```

Their first/middle span intersection is the line through r. Each cross-pair
of B and C axes anticommutes, and the same holds for B and A. If
$`N(P)=(I+P)/2+\omega(I-P)/2`$, the literal source

```math
F=[N(X_1)N(Z_0X_1)]\,[N(Z_1)N(Z_0Z_1)]\,
  [N(X_1)N(Z_0X_1)]
```

has control blocks $`\omega H_1`$ at $`r=+1`$ and
$`\omega^3 I_1`$ at $`r=-1`$. Thus

```math
F=\omega S_0\,\mathrm{C}_{r=+1}(H_1).
```

Each B band is literally
$`\mathrm{CNOT}_{1,0}(T_0T_1)\mathrm{CNOT}_{1,0}`$; each C or A band
is its conjugate by $`H_1`$. The displayed source and its actual reverse
inverse each use six T/T-dagger gates, ten elementary Cliffords and three
physical T slots, on two arbitrary qubits with no helper or clean qubit.

Its fair comparison is a one-band, two-T resynthesis. Put

```math
C_0=Z_1\,\mathrm{CZ}_{0,1},\qquad
E_Y=N(Y_1)N(Z_0Y_1).
```

The exact equality is

```math
F=S_0E_YC_0.
```

For example, $`W=S_1H_1T_1H_1S_1^\dagger`$ obeys
$`WZ_1W^\dagger=H_1`$ and
$`W=e^{i\pi/8}e^{-i\pi Y_1/8}`$. Since
$`C_0Y_1C_0=-Z_0Y_1`$, the expression
$`F=\omega S_0WC_0W^\dagger`$ becomes precisely $`S_0E_YC_0`$,
retaining its scalar phase.

The Y pair is one physical band in the explicit word

```math
E_Y=S_1H_1\,\mathrm{CNOT}_{1,0}(T_0T_1)
\mathrm{CNOT}_{1,0}H_1S_1^\dagger.
```

With $`S^\dagger=S^3`$, it costs two T calls and eight elementary
Cliffords. The self-adjoint $`C_0`$ costs five Cliffords and $`S_0`$ costs
one, giving **two T calls, fourteen elementary Cliffords and one T layer**
on the same arbitrary allocation. Its actual inverse uses the same paired
word with T-dagger gates, the same self-adjoint middle, and
$`S_0^\dagger=S_0^3`$: **two T calls, sixteen elementary Cliffords and one
T layer**. No claim of three-layer optimality or a count improvement is
made. All full operator phases and reverse inverses are checked; this is a
source fixture, not a returned XOR query.

## 6. An extension lemma for one dependent commuting list

Let U be diagonal C3 in fixed Clifford coordinates. Let E be a product of T/T-dagger phases about a commuting list of nonzero Hermitian Pauli axes $`a_1,\ldots,a_n`$ whose binary labels have exactly one linear relation. If U E is C3, then U E is diagonal C3 in fixed Clifford coordinates.

The fixed coordinate reduction for U conjugates all E axes together, preserving their commutation and relation count. Thus assume U diagonal. Native independence of E is explicitly not assumed; exactly one relation is the additional hypothesis. All orientations, axis signs and scalar phases are retained in the actual source, although they can be absorbed analytically into a right Clifford as described below.

### 6.1. Relation probes and paired half angles

Let J be the support of the unique nonzero relation, so $`\sum_{i\in J}a_i=0`$ in binary labels. Commutation vectors of Paulis against the complete list range over precisely the hyperplane whose coordinates sum to zero on J. For i outside J, a Pauli can isolate axis $`a_i`$; the usual half-angle argument shows $`Ua_iU^\dagger`$ is Pauli. For any distinct i,j in J, a Pauli can isolate exactly the pair $`a_i,a_j`$; cancelling its Clifford U image gives the paired half-angle criterion for $`H_i=Ua_iU^\dagger`$ and $`H_j=Ua_jU^\dagger`$, including their actual orientations.

If the relation has weight two, the two axes have the same binary label and their rotations combine into a scalar times a Clifford. Removing that Clifford leaves an independent native list; existing native right extension proves the claim. A weight-one relation would be a zero axis and is excluded by the statement; such a factor would only be scalar and could be dropped.

Assume $`|J|\ge3`$, so all its axes are distinct and nonzero. The paired-half-angle lemma gives $`H_i`$ $`H_j`$ Pauli for every pair. If one $`H_i`$ is Pauli, all are. Otherwise every $`H_i`$ has the same four-element Pauli support. They are mutually orthogonal in normalized Hilbert--Schmidt inner product, so there can be at most four of them. If $`|J|=3`$, their product is scalar identity, so $`H_1H_2`$ Pauli forces $`H_3`$ Pauli, a contradiction. Thus the only non-Pauli case has $`|J|=4`$.

If every active axis image is Pauli, use the signed diagonal Clifford extension G matching U on its complete Pauli-preserving domain. $`H=G^\dagger U`$ is diagonal and fixes all E axes literally. Its quotient Z observables commute with every axis, so simultaneous Clifford diagonalization of their joint isotropic span makes H E weighted-cubic diagonal. Restoring G proves the claim. This argument does not need independence of the E axes.

### 6.2. The four-axis exceptional pattern

It remains to consider $`|J|=4`$ and all four images non-Pauli with the same four-element support. Label one of these axes $`a_0`$ and write

```math
r_1=a_1+a_0,\quad r_2=a_2+a_0,\quad a_3=a_0+r_1+r_2.
```

The unique-relation premise makes $`a_0,r_1,r_2`$ independent. The paired lemma shows $`Ur_1U^\dagger`$ and $`Ur_2U^\dagger`$ are Paulis, with signs retained. It also shows that U preserves every axis outside J as a Pauli. Choose the phase-sensitive diagonal Clifford G matching U on all its Pauli directions, and replace U by $`H=G^\dagger U`$. Then H fixes every $`r_j`$ and every outside axis literally. Left multiplication by $`G^\dagger`$ preserves C3 membership and fixed-coordinate diagonalizability of the full source.

Because H is diagonal, its image of any Pauli has the same X label as that input Pauli. The four images retain a common support after the fixed Clifford correction $`G^\dagger`$. Consequently $`a_0,a_1,a_2,a_3`$ all have the same X label, and $`r_1,r_2`$ are Z-only Paulis. They commute with $`a_0`$ because the original E list is isotropic.

The actual orientations and signed representatives of E can now be normalized to positive T phases on the canonical axes

```math
a_0,\quad a_0r_1,\quad a_0r_2,\quad a_0r_1r_2
```

and on a chosen representative of each outside axis. A T orientation reversal changes a phase by a Clifford S rotation; a Hermitian axis sign reversal changes it by such a Clifford and a scalar. All axes commute, so these changes collect into a fixed right Clifford $`C_E`$ and a retained scalar. Right multiplication by $`C_E^\dagger`$ preserves C3 membership and the desired source class. Thus write the normalized list as $`E_4E_{\mathrm{out}}`$, where $`E_{\mathrm{out}}`$ is an independent native-axis band.

Let $`H_0=Ha_0H^\dagger`$. Since H fixes $`r_1,r_2`$, the other three images are $`H_0r_1,H_0r_2,H_0r_1r_2`$, up to the already fixed literal choices. Equality of their four-element supports implies that the support of $`H_0`$ is invariant under translations by the independent labels $`r_1,r_2`$, hence is one coset of their two-dimensional span.

All four Pauli terms in this support commute: their differences lie in $`\mathrm{span}\{r_1,r_2\}`$, and those Z-only Paulis commute with the common X label of $`a_0`$. Thus $`H_0`$ is a signed Pauli $`P_0`$ times a four-term real Hermitian unitary $`D_R`$ on the two commuting observables $`r_1,r_2`$. Its coefficients are all plus or minus 1/2. Such a four-term real unitary has, in their simultaneous eigenbasis, a Boolean phase with nonzero quadratic term. Equivalently,

```math
D_R=\pm r_1^{e_1}r_2^{e_2}M,
\qquad M=I-2\Pi_{r_1=r_2=+1},
```

for fixed bits e_1,e_2. The choice of corner in M changes only the displayed Pauli factors. Since $`P_0`$ has the same X label as $`a_0`$, $`P_0=a_0`$ times a Z-only Pauli up to a scalar. Therefore

```math
H_0=a_0 D_0,\qquad D_0M\text{ is a scalar times a Z Pauli}.
```

No claim that $`D_0`$ itself is Hermitian is needed. Its possible fourth-root scalar is retained; the final conjugated Pauli is Hermitian and therefore has a literal real sign.

### 6.3. A full Pauli Lagrangian for the exceptional source

Set $`V=HE_4`$. The normalized $`E_4`$ is particularly simple on the $`r_1,r_2`$ sectors: when both $`r_j`$ have eigenvalue +1, it is $`a_0`$; on each other sector it is iI. Hence every Z-only Pauli z commuting with $`a_0`$ commutes with $`E_4`$ and H, and V fixes it literally.

The X label of $`a_0`$ is nonzero, since otherwise H would commute with $`a_0`$ and $`H_0`$ would be Pauli. Choose a Z-only Pauli p anticommuting with $`a_0`$. It commutes with $`r_1,r_2`$, and the sector formula gives

```math
E_4 p E_4^\dagger=p M.
```

As H fixes every Z Pauli and every $`r_j`$,

```math
V a_0 V^\dagger=H_0=a_0D_0,
\qquad V p V^\dagger=pM.
```

Their product shows that V sends the Hermitian Pauli with label $`a_0`$+p to a signed Pauli, because $`D_0M`$ is scalar Pauli. The w-1 dimensional Z hyperplane commuting with $`a_0`$, together with $`a_0`$+p, is an isotropic w-dimensional Pauli label space. Thus V sends a maximal commuting Pauli subgroup to another such subgroup.

We still need V itself to be C3, rather than only the complete source. H fixes every outside axis, and every outside axis commutes with $`E_4`$, so V fixes them all literally. For any Pauli P, conjugation through $`E_{\mathrm{out}}`$ writes $`E_{\mathrm{out}}PE_{\mathrm{out}}^\dagger`$ as P times a product $`Q_P`$ of Clifford half-angle rotations about those axes. Since V fixes these axes, it fixes $`Q_P`$. The complete normalized source $`VE_{\mathrm{out}}`$ is C3, so

```math
(V E_{\mathrm{out}})P(V E_{\mathrm{out}})^\dagger=(V P V^\dagger)Q_P
```

is Clifford. Cancelling the actual Clifford $`Q_P`$ proves $`VPV^\dagger`$ Clifford. Thus V is C3.

A C3 unitary mapping a Pauli Lagrangian to a Pauli Lagrangian is diagonal C3 in fixed Clifford coordinates: put the signed input and output Lagrangians into the standard Z frame, so that the resulting unitary commutes with every Z and is diagonal; its third-level finite-difference condition gives the usual weighted cubic phase form modulo eight. Applying the existing native right-extension corollary to V and the independent outside band $`E_{\mathrm{out}}`$ now makes $`VE_{\mathrm{out}}`$ diagonal C3. Restore the actual G, $`C_E`$ and scalar to prove the statement.

## 7. Operator meaning and remaining questions

The one-overlap theorem is exact and concerns the specified native source.
The proof does not assume $`F^\dagger`$ is C3, assume intermediate source
passages are Clifford, choose initialized work sectors, or infer a physical
band from a dependent commuting axis list. Sector restrictions in section 4
are analytical blocks of one full-space operator, and their literal eighth-root
phases are reassembled before the source conclusion.

Every original source call, T gate, elementary Clifford, actual inverse,
control and simultaneous work allocation remains charged. Fixed coordinates,
Pauli spans, spectral identities and diagonal corrections introduce no free
physical gate, initialized state or helper. All source equalities hold on
arbitrary work inputs and external references. The six-layer corollary uses
the original full logical action, literal target phase, leakage estimate and
complete work return that supply its comparison bound.

For $`t\ge2`$, the retained $`2t`$ bound on noncommuting middle axes does
not yet prove diagonal-C3 form or a bounded active subsystem. The general
three-band source problem, general six-layer signed transport, and error
control for a depth recurrence remain open. The new source theorem does not
improve the complete-frame upper law.

The [focused checks](../../tests/test_one_overlap_source_structure.py) cover
paired half-angle rigidity, the weight-one native support constraint, the
four-axis dependent exception and its alternative Pauli Lagrangian, the
weight-two symplectic cancellation, literal C3 source blocks, and the actual
native source/inverse ledgers with a two-T fair baseline. These checks support
the proof mechanisms; the arbitrary-width statements are analytical.
