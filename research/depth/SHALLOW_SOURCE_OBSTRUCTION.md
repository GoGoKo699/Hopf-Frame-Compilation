# A robust obstruction to two-layer dirty source circuits

**Preserved research study.** This note is outside the selected A–D proof chain. Its outcome and limits are indexed in the [research archive](../README.md); historical proposals are not current work orders. The local mathematical statements retain their stated hypotheses.


[Source primitives](../../docs/OPERATOR_SOURCE_COMPILER.md) · [Restricted source depth](SOURCE_T_DEPTH.md) · [Current frontier](../../docs/OPEN_PROBLEM.md)

The geometric source cannot be approximated arbitrarily well by two T
layers, even with unrestricted Clifford interlayers and arbitrarily many
returned dirty helpers. For source width at least five, the operator-norm
error is at least $`1/16`$. The same bound holds for its controlled version
at width at least four. These are constant-depth obstructions for
full-input source primitives, not growing depth lower bounds or lower
bounds for an initialized-clean complete-frame compiler.

[Section 6](#6-native-depth-at-the-exact-t-count-minimum) gives a
separate depth result at exact minimum T count. It allows initialized
stabilizer helpers under exact complete return.

The later [clean-aware obstruction](CLEAN_AWARE_DEPTH_OBSTRUCTION.md)
uses the transfer alphabet below inside the initialized-projector sum.
It proves a one-clean frame gap and a separate three-layer obstruction
for every fixed clean budget; those extensions have their own proof.

The proof combines standard stabilizer overlap quantization with the
source's exact Pauli-transfer coefficients. The overlap fact is given in
[Aaronson–Gottesman, *Improved Simulation of Stabilizer Circuits*,
arXiv:quant-ph/0406196v5, Section III, final paragraph, p. 5](https://arxiv.org/pdf/quant-ph/0406196v5).
The observation that one T layer conjugates Paulis to Hermitian Cliffords
is also used explicitly in
[Zhang–Zhang, arXiv:2409.13809v2, Theorem III.1, equations (10)–(11)](https://arxiv.org/html/2409.13809v2#S3.SS1).
The transfer-entry reduction and source witnesses below specialize these
inherited facts; no general priority claim is made.

## 1. A discrete Pauli-transfer alphabet at T-depth at most two

Let V be a circuit on q physical qubits with T-depth at most two. Every
wire is included in its full unitary action. Clifford blocks may be
arbitrary, and a T layer may contain T or T-dagger on any set of distinct
wires. For Hermitian Pauli strings P,Q, define

```math
\widehat V_{P,Q}=2^{-q}\mathrm{Tr}(PVQV^\dagger).
```

**Lemma.** Every such coefficient belongs to

```math
\mathcal A_q=\{0\}\cup\{\pm2^{-j/2}:0\le j\le2q\}.
```

*Proof.* Write $`V=C_2D_2C_1D_1C_0`$, where each C is Clifford and
each D is one physical T layer; identity layers are allowed. Absorb the
outer Cliffords into the two Pauli strings, obtaining P',Q'. Cyclicity
of the trace gives

```math
\widehat V_{P,Q}
=2^{-q}\mathrm{Tr}\!\left[
 (D_2^\dagger P'D_2)
 C_1(D_1Q'D_1^\dagger)C_1^\dagger\right].
```

Conjugating a one-qubit Pauli by T or T-dagger gives a Hermitian Clifford:
X and Y become signed combinations $`(X\pm Y)/\sqrt2`$, while I,Z
remain Paulis. Tensor products and Clifford conjugation preserve this
property. The two displayed factors are therefore Hermitian Clifford
unitaries. Their product K is Clifford, and its trace is real because
the trace of a product of two Hermitian matrices is real.

For the normalized q-pair Bell state $`|\Omega_q\rangle`$,

```math
2^{-q}\mathrm{Tr}K
=\langle\Omega_q|(K\otimes I)|\Omega_q\rangle.
```

Both vectors are stabilizer states on 2q qubits. Their overlap magnitude
is either zero or $`2^{-j/2}`$ for an integer $`0\le j\le2q`$.
This follows directly by expanding their stabilizer projectors: the
squared overlap is zero if the common Pauli constraints conflict, and
otherwise is $`2^{s-2q}`$, where the common stabilizer subgroup has size
$`2^s`$. Reality of the trace gives the stated signs. In particular, no
assumption about the phase of a general Clifford trace is needed. ∎

This is a statement about individual full-space transfer entries. It
does not assert that a state-dependent Pauli expectation in a circuit
with two T layers belongs to the same alphabet.

## 2. Converting a transfer gap into robust operator distance

Put

```math
\mathcal A=\{0\}\cup\{\pm2^{-j/2}:j\ge0\}.
```

For any two unitaries U,V and any Hermitian Paulis P,Q on the same
physical space,

```math
|\widehat U_{P,Q}-\widehat V_{P,Q}|
\le\|UQU^\dagger-VQV^\dagger\|
\le2\|U-V\|.
```

The first step uses the normalized trace bound, and the second expands
the conjugation difference into two terms. Therefore, if one transfer
entry of U is x, every circuit V of T-depth at most two satisfies

```math
\|V-U\|\ge\frac12\mathrm{dist}(x,\mathcal A).
```

Tensoring U and both witness Paulis with identity preserves x. The
bound consequently survives any number of dirty helper wires, including
arbitrary intermediate mixing with them. The target is the full unitary
$`U\otimes I_{\rm dirty}`$; even approximate helper return is included
in this operator-norm distance. The same statement holds after optimizing
over a common global phase, since transfer coefficients ignore that phase.

## 3. Geometric source witnesses

For the existing geometric source

```math
M_m=\sum_{j=0}^{m-1}a_j\Gamma_j,\qquad
\Gamma_j=Z_0\cdots Z_{j-1}X_j,
```

write $`w_j=a_j^2`$. The weights are
$`w_j=2^{-j-1}`$ for $`j\lt m-1`$, with the duplicated last weight
$`w_{m-1}=2^{-(m-1)}`$. For any index subset S, the Pauli
$`P_S=\prod_{j\in S}Z_j`$ flips exactly those Majorana terms indexed
by S. Orthogonality of distinct Pauli strings gives

```math
2^{-m}\mathrm{Tr}(P_SM_mP_SM_m)
=1-2\sum_{j\in S}w_j.
```

At $`m\ge4`$, the single-qubit witness $`P=Z_2`$ has coefficient
$`3/4`$. Its nearest alphabet point is $`1/\sqrt2`$, proving

```math
\|V-M_m\otimes I_{\rm dirty}\|
\ge\frac{3/4-1/\sqrt2}{2}\approx0.0214466
\qquad(D_T(V)\le2).
```

At $`m\ge5`$, use $`P=Z_3`$. The coefficient is $`7/8`$, whose
nearest alphabet point is one. This yields the simpler stronger bound

```math
\|V-M_m\otimes I_{\rm dirty}\|\ge\frac1{16}.
```

Thus exact full-input implementations of these sources require at least
three T layers, as do approximations below the displayed tolerances.
No growing lower bound in m follows from this argument.
At $`m=3`$, the existing optimized source word has one T layer on each
side of its Clifford middle, so it implements the source exactly with
two T layers. The width threshold for excluding two layers is therefore
sharp, without asserting source-depth optimality at larger widths.

### Optimizing over the source's diagonal witnesses

The subset weights attain every multiple of $`2^{-(m-1)}`$ between
zero and one. Indeed, the first m−1 integer-scaled weights are the binary
digits $`2^{m-2},\ldots,1`$, and the final duplicate supplies the endpoint
$`2^{m-1}`$. Hence the displayed transfer entries cover the entire grid

```math
\mathcal G_m=\{-1,-1+h_m,\ldots,1\},\qquad h_m=2^{2-m}.
```

Let

```math
G=1-1/\sqrt2,\qquad t_*=(1+1/\sqrt2)/2,
\qquad x_m=h_m\mathrm{round}(t_*/h_m).
```

For $`m\ge5`$, the best gap among these diagonal witnesses is exactly

```math
\max_{x\in\mathcal G_m}\mathrm{dist}(x,\mathcal A)
=\frac G2-|x_m-t_*|.
```

To see optimality within this witness family, the interval
$`[1/\sqrt2,1]`$ has the largest alphabet gap. The available point
$`7/8`$ already achieves distance $`1/8`$. Every smaller positive
alphabet gap has half-width at most
$`(1/\sqrt2-1/2)/2\lt1/8`$, and negative gaps are symmetric.
Within the largest gap, the grid point nearest its midpoint is optimal.

Consequently the lower bound on operator distance is

```math
\frac G4-\frac12|x_m-t_*|
\ge\frac G4-\frac{h_m}{4}
\longrightarrow\frac{1-1/\sqrt2}{4}\approx0.0732233.
```

This optimizes a single diagonal transfer witness, not the true distance
to the set of all two-layer circuits.

## 4. Controlled source

For the controlled source

```math
C(M_m)=|0\rangle\langle0|\otimes I+
|1\rangle\langle1|\otimes M_m,
```

the witness
$`I_{\rm control}\otimes P_S`$ has diagonal transfer coefficient

```math
\frac12\left(1+1-2\sum_{j\in S}w_j\right)
=1-\sum_{j\in S}w_j.
```

These coefficients fill $`[0,1]`$ with spacing $`2^{1-m}`$.
For $`m\ge4`$, $`S=\{2\}`$ gives $`7/8`$, and therefore every
full-input circuit V of T-depth at most two satisfies

```math
\|V-C(M_m)\otimes I_{\rm dirty}\|\ge\frac1{16}.
```

The control is an arbitrary input qubit, not an initialized helper.
The optimized formula from Section 3 applies with spacing $`h_m/2`$
for $`m\ge4`$, with the same limiting constant. At $`m=3`$, the
$`Z_1`$ witness instead gives coefficient $`3/4`$ and the smaller
positive bound from Section 3.
At $`m=2`$, the existing one-layer loader and its actual inverse around
the controlled central Pauli give an exact controlled source with two
T layers, so this exclusion also starts at the smallest possible width.

## 5. Scope

These bounds allow arbitrary Clifford interlayers and arbitrarily many
dirty helpers. They therefore differ from the growing but restricted
Majorana-linear depth bounds in [the source-depth chapter](SOURCE_T_DEPTH.md).
Neither argument replaces the other.

An initialized-clean isometry contract is excluded. Restricting input
ancillas to zero inserts their stabilizer projector into the trace and
can sum many full-space transfer coefficients; that sum need not belong
to the alphabet. Closeness on those initialized columns also does not
imply full operator-norm closeness. The two clean flags in the complete
frame compiler cannot be silently treated as dirty inputs to obtain a
frame lower bound.

The [clean-aware extension](CLEAN_AWARE_DEPTH_OBSTRUCTION.md#1-initialized-work-and-a-robust-transfer-witness)
keeps this projector sum and proves robustness on initialized columns
without making that full-space inference. Its fixed-clean frame bounds
therefore do not change the scope of the source theorem above.

The robust two-layer result does not establish source depth optimality,
a lower bound for a source loader, an additive cost across source calls,
or a growing unrestricted complete-frame depth bound. A joint compiler may avoid
approximating this source as an isolated full-input primitive. The
[complete-frame precision tradeoff](PRECISION_COUNT_DEPTH_TRADEOFF.md) uses a
different witness and gives a lower bound directly on the compiled frame.

The [bounded checks](../../tests/test_shallow_source_obstruction.py) audit
native two-layer words, exact rational subset grids, the optimized
coefficient gaps, native source witnesses with dirty extensions, and
the small-width two-layer constructions. Their finite matrices support
these identities; the uniform exclusion and robustness follow from the
proofs above.

## 6. Native depth at the exact T-count minimum

Let $`m\ge2`$. At its minimum T count $`2m-4`$, every exact native
implementation of $`M_m`$ has T-depth $`2m-4`$. At minimum T count
$`2m-2`$, the least achievable T-depth for either control value of
$`C(M_m)`$ is $`2m-3`$. An exact minimum-count schedule of the
controlled source has at most one two-gate T band. All its other
nonempty T bands contain a single gate.

The statements allow arbitrary Clifford interlayers, arbitrary dirty
helpers, and any number of initialized stabilizer helpers. The exact
action includes every logical and dirty input with complete helper
return. The hypotheses use the
[exact T-count minima](../../docs/OPERATOR_SOURCE_COMPILER.md#exact-source-costs-including-returned-helpers).

### Leading Pauli coefficients constrain each band

Use Hermitian Pauli labels over the binary symplectic space, with
$`\langle a,p\rangle`$ the commutation bit of labels a and p. For a
Pauli expansion with coefficients in $`\mathbb Z[1/\sqrt2]`$, its
least denominator exponent is the smallest $`d\ge0`$ for which
scaling all coefficients by $`(\sqrt2)^d`$ puts them in
$`\mathbb Z[\sqrt2]`$. Reduce these scaled coefficients modulo
$`\sqrt2`$ to obtain a vector over $`\mathbb F_2`$. At the least
exponent this vector is nonzero. Call its support the leading support.

Clifford conjugation preserves d. Conjugation by a T or T-dagger
increases d by at most one. When it increases d by one, the leading
vector of the result is obtained from the previous leading vector
by the binary linear map

```math
e_p\longmapsto
\begin{cases}
e_p+e_{p+a},&\langle a,p\rangle=1,\\
0,&\langle a,p\rangle=0,
\end{cases}
```

where a is the rotation-axis label. The active terms acquire
coefficients of magnitude $`1/\sqrt2`$, and their signs disappear in
the reduction. The inactive coefficient gets the extra factor
$`\sqrt2`$ from the new scaling.

Pull all fixed Cliffords to the input end of the complete word.
The axes of each original physical T band remain commuting and
independent because they are Clifford images of distinct single-wire
Z labels. Let A be their binary span. Translating a Pauli label by
any element of A leaves every commutation bit with the band's axes
invariant. When each gate increases the exponent by one, the
band's output leading vector is constant on each coset of A and
vanishes on cosets inactive under any of its axes. A fully active
basis label contributes the sum of all its translations by A. A
label that is inactive for one axis contributes zero.

Write $`C_{\rm cl}`$ for the span of the clean Z labels, after fixed
Cliffords change the given stabilizer initialization to zero. Use
the source test Pauli $`P=Z_{m-1}`$, tensored with identities on any
control and dirty work. Propagate

```math
P\sum_{z\in C_{\rm cl}}Z^z
```

by conjugation through the actual word. The sum is the initialized
clean projector times its integer rank factor. Exact return makes
the final observable the prescribed conjugated logical P times the
same sum, including identity on dirty work. Its exponent equals the
exact minimum T count in each case, as also verified below. The initial
exponent is zero, so every actual T increases it by one. Reversing a
band of r T gates must decrease the exponent by r.

Suppose at a reverse boundary that the leading support has the form

```math
p+C+\{0,d_0\},\qquad d_0\notin C,
```

where C is either $`C_{\rm cl}`$ or its span with the control Z label.
The labels p and $`p+d_0`$ have no clean component. Every axis in the preceding band belongs
to $`C+\langle d_0\rangle`$, by translation invariance of the
leading support. The individual axes must be active on each label
in the support. Elements of C commute with every supported label.
Thus each axis lies in $`d_0+C`$.

All these axes commute with the true clean Z labels. Their reverse
preserves the initialized cylinder and its observable stays a
logical expansion times the same clean stabilizer sum. On the
cylinder each axis restricts, up to sign, to its label modulo
$`C_{\rm cl}`$. Repeated axes after this restriction combine into
a Clifford or identity for each pair of T rotations, and the
restricted rotations commute. The logical denominator exponent
equals that of the full expansion because the clean sum repeats its
coefficients without additional factors. Hence the number of
distinct restricted axis labels must be at least r to realize the
necessary exponent decrease. If $`C=C_{\rm cl}`$, the band is
therefore a singleton. If C also contains control Z, there are
at most two gates. In the two-gate case the projected axes are
$`d_0`$ and $`d_0+Z_{\rm control}`$. Their coefficient rotations
have angles $`\pm\pi/4`$. The angles cancel in one control branch
and add to $`\pm\pi/2`$ in the other, giving a signed swap of the
two mixed labels there.

### Uncontrolled reflection

Put $`t=m-2`$, $`q=2^{-t}`$ and $`\Gamma_*=\Gamma_{m-1}`$. The
Pauli P anticommutes with $`\Gamma_*`$ and commutes with every other
source generator. The output logical expansion is

```math
Q=M_m P M_m=(1-q)P-\sum_{j=0}^{t}2^{-(t+j)/2}A_j,
\qquad A_j=P\Gamma_*\Gamma_j .
```

The labels of P and the Hermitian Paulis $`A_j`$ are distinct and
pairwise anticommuting. For $`t>0`$, the two leading coefficients
are those on P and $`A_t`$, with exponent $`2t`$. This proves
the exponent premise for the minimum-count word. For $`t=0`$
the exact source is Clifford and there are no T bands to consider.

After reversing $`\ell\lt t`$ singleton bands, the coefficient
magnitudes consist of the following entries on distinct labels.

- One parked magnitude $`2^{-j/2}`$ for each $`1\le j\le\ell`$.
- A remaining coefficient of magnitude
  $`2^{-\ell/2}(1-2^\ell q)`$.
- The untouched magnitudes $`2^{-(t+j)/2}`$ for
  $`0\le j\le t-\ell`$.

The two leading labels are the remaining coefficient and the
last untouched one, with exponent $`2t-\ell>\ell`$. With the
clean sum their support has the two-label form above for
$`C=C_{\rm cl}`$. The next band is a singleton. Its projected
axis is the sum of those two labels, which commutes with every
other label by pairwise anticommutation. The rotated magnitudes
are the absolute sum and difference of the two old magnitudes,
each divided by $`\sqrt2`$. They give the new parked magnitude
and the new remaining coefficient in the displayed list. This
proves the list inductively, regardless of the signs of individual
coefficients or rotations.

At $`\ell=t`$ the remaining coefficient vanishes. The nonzero
magnitudes are now

```math
2^{-1/2},2^{-2/2},\ldots,2^{-t/2},2^{-t/2}.
```

The top pair again forces a singleton. Its equal magnitudes combine
to one coefficient of magnitude $`2^{-(t-1)/2}`$ and one zero
coefficient. Repeating this equal-pair step removes another t
singleton bands and leaves one coefficient of magnitude one. Thus
the exact minimum-count word has $`2t=2m-4`$ singleton bands.

### Controlled reflection

Use control value one. Conjugating the entire word by X on the
control proves the other case. Ignore the clean sum while describing
coefficients separately on the two control branches. Start with
coefficient vectors $`u=P`$ and $`v=Q`$ in the zero and one branches.
The physical coefficients on the I and Z control labels are
$`(u+v)/2`$ and $`(u-v)/2`$. This gives exponent $`2t+2`$.
For $`t=0`$ the vectors are on distinct labels with magnitude one,
so the same exponent formula holds.

For $`\ell\lt t`$, the first $`\ell`$ reverse gates produce the same
magnitude list for v as in the uncontrolled proof. The branch u uses
the same two axis labels at each step, possibly with a different
rotation sign. Its untouched partner always starts at coefficient
zero. The rotation assigns equal magnitudes to the new parked label
and the label that still carries the remaining coefficient. Thus
u has the same parked magnitudes as v and a remaining magnitude
$`2^{-\ell/2}`$. Its exponent is at most $`\ell`$, strictly less
than that of v. The two leading labels in v therefore occur with
both control dressings in the leading support of the full expansion.
Here C includes control Z. A proposed two-gate band would act by a
signed permutation on each branch, preserving both branch exponents.
Their strict inequality would keep the combined exponent unchanged,
contradicting the required decrease. The next band is therefore a
singleton. This proves that the first t reverse bands are singletons.

At this point the branches share all parked magnitudes. In addition,
u has its final remaining label and v has its untouched $`A_0`$
label, with equal magnitudes $`2^{-t/2}`$ on these distinct labels.
The shared parked terms each have only one nonzero control dressing
after folding, since their two branch magnitudes are equal. The
two distinct unmatched labels instead have both dressings and give
exponent $`t+2`$. They force a band of size at most two.

A two-gate band swaps these two labels in one control branch. It
aligns the nonzero magnitude on a single label across the branches.
If this band instead has one gate, it splits the unmatched magnitude
on each branch into a pair of magnitudes $`2^{-(t+1)/2}`$. The
two split labels now have just one nonzero control dressing each.
They are the only leading pair at exponent $`t+1`$, so the next band
is a singleton. On each branch the two consecutive rotations of the
same logical axis restore one coefficient of magnitude $`2^{-t/2}`$
on the pair. The required new exponent t forces them onto the
same label across branches. Separate labels would again give
exponent $`t+2`$ after folding.

In either case, the expansion now has one control dressing on each
label of the equal-pair list from the uncontrolled proof. These
dressed labels are still pairwise anticommuting because their
control factors are only I or Z. The equal-pair collapse forces
the remaining t bands to be singletons. Thus the count-minimum
depth is at least $`2t+1=2m-3`$, and only the two middle gates can
share a band. The middle argument also covers $`t=0`$, where the
parked list is empty.

The [controlled source circuit](../../docs/OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit)
attains the bound. In its middle, write the first split and its
reverse around the controlled X as

```math
R_0 C(X_0) R_0^\dagger
=R_0\bigl[C(X_0) R_0^\dagger C(X_0)^\dagger\bigr]C(X_0).
```

Let $`D_0=i\Gamma_0\Gamma_1`$ be the Pauli axis of $`R_0`$.
Conjugating it by $`C(X_0)`$ gives $`Z_{\rm control}D_0`$.
These two Pauli axes are commuting and independent. One Clifford
basis change sends them to distinct physical Z axes, so the two
actual signed T rotations fit one native band. The two outer
sequences have t T gates each. The rearrangement keeps the native
inverse and the literal scalar of the original source.

The exact count hypothesis is essential to this proof. It ensures
one exponent gain at every actual T and thereby constrains the
leading support of every band.

### Tradeoff for an approximate isolated source

The [precision determinant bound](PRECISION_COUNT_DEPTH_TRADEOFF.md#precision-count-depth-tradeoff)
also applies to an approximate $`M_m`$ through the following reduction.
Suppose an actual word V for $`M_m`$ has T count T, native T-depth d,
k clean inputs and complete initialized-isometry error $`\epsilon`$.
Put $`s=\min\{m,\log_2(1/\epsilon)\}`$, using $`s=m`$ if
$`\epsilon=0`$. For sufficiently large s the bounds satisfy

```math
(T+1)(d+1)=\Omega\!\left(\frac{s^2}{k+3}\right).
```

Fix the angle $`\cos(2\theta)=1/3`$ from the determinant bound.
Use the [signed coefficient grid](../../docs/OPERATOR_SOURCE_COMPILER.md#2-exact-signed-dyadic-coefficients)
to choose fixed Pauli masks $`P_0,P_1`$ for its cosine and sine.
For a selector bit b, the ideal word
$`W_b=M_m P_b M_m^\dagger P_b`$ has Hermitian part $`c_b I`$,
where $`c_b`$ is the selected grid value. A Hadamard on an additional
flag a, two a-controlled visits to $`P_b`$ around $`W_b`$, and a final
Hadamard give the accepted block

```math
\frac{W_b+P_b W_b P_b}{2}=c_b I.
```

The mask is $`P_0(P_0P_1)^b`$. Its singly controlled Pauli is
Clifford. A Clifford basis change to a single Pauli axis gives
constant T cost and depth for the doubly controlled Pauli.

The [two-flag rotation and three-call recovery](../../docs/OPERATOR_SOURCE_COMPILER.md#4-two-flags-encode-a-suffix-controlled-rotation)
use this block with a and b initialized. The grid error gives
$`O(2^{-m})`$ full initialized-isometry error for the ideal-source
circuit. Its reflection tests only a and b. Substitute V and its
actual inverse for the two source visits in each block. The inverse
has the same initialized-isometry error by unitary multiplication.
A six-term isometry hybrid adds at most $`6\epsilon`$, since every
ideal intermediate keeps the source helpers in their initialized
cylinder. This constructs the target angle with k plus two clean
inputs, error $`O(2^{-m}+\epsilon)`$, T count $`O(T+1)`$ and native
T-depth $`O(d+1)`$. The determinant bound gives the displayed tradeoff.
In particular a fixed-clean, linear-count implementation with
exponentially small source error requires linear native T-depth.
