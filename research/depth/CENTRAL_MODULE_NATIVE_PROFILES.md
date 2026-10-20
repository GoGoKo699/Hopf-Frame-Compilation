# Bounded native profile distortion modulo the genuine center

[Occupied common center](OCCUPIED_BOOLEAN_PAIR_CENTER.md#9-a-clean-preserving-clifford-sum-without-bounding-the-silent-span) ·
[Finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md)

This uses the complete occupied Boolean identity and the actual logical
three-band reduction. It is uniform in the additional silent-axis span t.
It proves an orbit-profile invariant, not native axis matching or a degree bound.

## 1. Actual source and full logical comparison

First use [Sections 3–4 of the occupied-center note](OCCUPIED_BOOLEAN_PAIR_CENTER.md#3-the-complete-clean-transition-and-its-retained-logical-bands): keep the physical clean
coordinate through its exact common-center and Clifford-extension steps,
then restrict the already clean-preserving factors to the positive code.
This gives the literal full-logical-space word

```math
G=K_C^\dagger E_C^\dagger F E_BE_AK_A,
\qquad F=G_0|_{P_A=+1}.
```

Here E_A and E_C are the actual projected silent endpoint lists, E_B is
the actual projected positive-singleton list, and all three are independent.
Their axes and F commute the injectively projected genuine Pauli center R.
F is unitary and has Pauli support in at most $`2^c`$ cosets modulo R,
where $`c=\dim(H_0/R)\le3s+1`$ (duplicates after clean restriction can reduce it).
Also F commutes each actual B axis. The shared computational alignment of
Section 4 ensures $`K_A^\dagger R K_A`$ and
$`K_C^\dagger R K_C`$ are computational-Z groups while the
residual target G remains Boolean diagonal. Its removed quadratic and
fixed coordinate Cliffords are retained as in that note.

For every computational Z probe Z(z), define the actual Hermitian Cliffords

```math
J_z=E_AK_AZ(z)K_A^\dagger E_A^\dagger,
\qquad C_z=E_CK_CZ(z)K_C^\dagger E_C^\dagger.
```

Every J_z and C_z fixes R pointwise: the encoded probe commutes R by the
shared alignment, and each endpoint phase band also commutes R. The
full Boolean identity gives, on the entire logical space,

```math
E_BJ_zE_B^\dagger=F^\dagger C_zF.                 \qquad(1)
```

No identity on an unused physical clean sector is inferred here. All
logical probes and arbitrary logical/dirty/reference inputs are included.
Write $`\Phi_z`$ for the binary Pauli action of J_z; the whole family is used,
not just a generating set.

## 2. A norm with central operator coefficients

Choose a signed commuting Pauli group representing R, without -I, and
one fixed unitary Pauli representative P_a for each binary coset
$`a\in V/R`$.
Every operator has a unique expansion

```math
X=\sum_{a\in V/R}f_a(R)P_a.
```

The coefficient f_a(R) belongs to the commutative central Pauli algebra.
Here “central” means the chosen module coefficient algebra; arbitrary P_a
may anticommute with some R generators. Its norm is the supremum over
the actual simultaneous R spectrum. This is the whole character group:
R projects injectively to the full logical space.

Define

```math
\|X\|_{R,1,\infty}=\sum_a\|f_a(R)\|_\infty.
```

Changing a representative multiplies its coefficient by a central Pauli
and a fixed unit scalar, so the norm is independent of that choice.
Pauli multiplication has the form
$`P_aP_b=u_{a,b}(R)P_{a+b}`$, with u a central unitary Pauli times its
literal scalar. Commuting P_a past g_b(R) merely permutes the R characters.
The triangle inequality and the scalar sup norm therefore give

```math
\|XY\|_{R,1,\infty}\le
\|X\|_{R,1,\infty}\|Y\|_{R,1,\infty}.
```

Conjugation by any actual Clifford normalizing the signed R algebra
permutes cosets and characters and multiplies coefficients by unit
central Paulis. It is an isometry of this norm. In particular this applies
to all J_z and C_z in (1), with their literal signs.

For a unitary U, each coefficient has sup norm at most one: it is the
conditional expectation onto Alg(R) of UP_a†, hence is an operator-norm
contraction. The conditional expectation can be realized by the finite
Pauli twirl over $`R^\perp`$; no physical implementation is being proposed.
Thus the actual source gives separately

```math
\|F\|_{R,1,\infty},\ \|F^\dagger\|_{R,1,\infty}\le2^c.  \qquad(2)
```

The degree/features of f_a are unnecessary for this estimate. Their
bounded feature ledger remains available from §9.

## 3. Exact independent-native profile

For any Pauli label p on the complete logical register, let
B(p) be the active sublist of E_B: its axes anticommute with P. Write
k=|B(p)|. Since the actual projected list is independent, products of
these axes give an injective map from $`\mathbb F_2^k`$ to Pauli labels. Define

```math
\begin{aligned}
K_p&=\{u\in\mathbb F_2^k:\text{ the active-axis product lies in }R\},\\
h_p&=\dim K_p,\qquad
r_p=\mathrm{rank}(u\cdot v\text{ restricted to }K_p).
\end{aligned}
```

Signed representatives of a product affect only a character and do not
change K_p or its ordinary binary dot form. Set

```math
\nu_R(p)=k-r_p.
```

Then exactly

```math
\|E_BPE_B^\dagger\|_{R,1,\infty}=2^{\nu_R(p)/2}.    \qquad(3)
```

Proof: the literal independent-band expansion has $`2^k`$ Pauli terms of
magnitude $`2^{-k/2}`$. Their phases are i to an odd weighted sum of the
active bits, with original orientations and Pauli signs retained. Group
by K_p cosets. Up to a constant unit phase and a linear character, each
coefficient on an R eigensector is

```math
2^{-k/2}\sum_{u\in K_p}i^{\mathrm{wt}(u)}(-1)^{\ell(u)}.
```

Every character ell is attained as the genuine R character varies;
orientation signs and coset choices only shift ell. Let
$`T_p=K_p\cap K_p^\perp`$, so dim T_p=h_p-r_p. For the displayed sum S(ell),

```math
|S(\ell)|^2
 =2^{h_p}\sum_{u\in T_p}i^{\mathrm{wt}(u)}(-1)^{\ell(u)}.
```

Every radial u has even weight, and $`i^{\mathrm{wt}(u)}`$ is a character
on T_p. Thus the squared magnitude is zero or $`2^{2h_p-r_p}`$, with the latter
attained. This holds for arbitrary P, even if P anticommutes with R;
K_p itself need not be even-weight. There are $`2^{k-h_p}`$ cosets, so summing
their individual suprema proves (3). No arbitrary sectorwise replacement
of a native list has been made.

## 4. Uniform distortion for the actual Boolean matching family

Put $`X=E_BPE_B^\dagger`$. From (1), including the actual inverse,

```math
E_B(J_zPJ_z^\dagger)E_B^\dagger
 =F^\dagger C_zF X F^\dagger C_z^\dagger F.
```

There are four F or F† factors. By (2), submultiplicativity and the
Clifford isometry, the norm ratio is at most $`2^{4c}`$. The reverse identity
gives the same bound in the opposite direction. Equation (3) yields

```math
|\nu_R(\Phi_zp)-\nu_R(p)|\le8c                   \qquad(4)
```

for every logical Pauli p and every member z of the complete J family.
This is the promised native profile bound without bounded t.

## 5. Companion profiles and their exact limitation

Grouped support count also gives a useful phase-independent statement.
Every grouped Gaussian coefficient is nonzero somewhere on the actual
R spectrum, hence the grouped support of $`E_BPE_B^\dagger`$ has cardinality
$`2^{k-h_p}`$. Support count is submultiplicative and invariant under
Clifford normalizers of R. The same four-F comparison gives

```math
|\delta_R(\Phi_zp)-\delta_R(p)|\le4c,
\qquad \delta_R(p)=k-h_p
 =\dim\bigl(U_p/(U_p\cap R)\bigr),
```

where U_p is the active-axis span.

There is a second norm,

```math
M_R(X)=\max_\chi\sum_a|f_a(\chi)|.
```

It is submultiplicative: at a fixed character, sum the shifted second
factor at each Pauli-induced character and bound that sum by its common
maximum. Clifford normalizers are again isometries. For each unitary U,
the conditional expectation identity gives

```math
\mathcal E_R(UU^\dagger)=I,
\qquad \sum_a|f_a(\chi)|^2=1.
```

Cauchy--Schwarz bounds both $`M_R(F)`$ and $`M_R(F^\dagger)`$ by
$`2^{c/2}`$, separately. At each character, the Gaussian calculation in
Section 3 leaves exactly $`2^{k-2h_p+r_p}`$ nonzero cosets, each with
magnitude $`2^{-(k-2h_p+r_p)/2}`$. Indeed their nonvanishing condition
fixes the character on the radical of K_p, of dimension h_p-r_p; the
coset pairing onto that radical's dual is surjective. Therefore

```math
M_R(E_BPE_B^\dagger)=2^{\mu_R(p)/2},
\qquad \mu_R(p)=k-2h_p+r_p,
\qquad |\mu_R(\Phi_zp)-\mu_R(p)|\le4c.
```

This second norm measures an individual central fiber and loses some
coherent central-control information; it does not dominate the usual
operator norm for arbitrary nonnormal operators, which is not needed.
The first norm retains ν_R, whereas ν_R, δ_R, μ_R are different from
the raw active count k. Large central relations can make these profiles
much smaller than k. Their bounded variation along the actual native
matching family does not give an absolute rank or span bound, a bounded
t, or a Boolean-degree theorem. The remaining task is to recover useful
actual native matching from these simultaneous corrected profiles while
retaining their central relation forms and all source lists.

## 6. A per-element escape-span consequence

Let L be the span of the actual positive middle axes. For one actual
matching action $`\Phi_z`$, put $`b_i'=\Phi_z^{-1}b_i`$ and select its escaping axes
b_i' outside L+R. Then

```math
d_z=\dim\Bigl(\bigl(\mathrm{span}\{b_i':b_i'\notin L+R\}+R\bigr)/R\Bigr)
\le8c.
```

Choose d_z escaping axes independent modulo R. For a uniformly random
Pauli label p in $`(L+R)^\perp`$, original δ_B(p)=0 and every chosen escaping
axis activates with probability one half. No independence of these
activation events is required. Their active images remain independent
modulo R, so $`\mathbb E\,\delta_{\Phi_z^{-1}B}(p)\ge d_z/2`$.
Since Phi_z fixes R, this rank is $`\delta_B(\Phi_zp)`$, which is
pointwise at most 4c by the profile theorem. Thus $`d_z\le8c`$. Reversing Φ_z gives the other direction.

This bounds the full modulo-R span of escaping axes, stronger than only
its quotient outside L+R. The transformed union over z requires care:
the original exceptional axes admit a uniform bound, but all of their
images need not. Neither conclusion below classifies transformations
inside L+R or proves bounded t or Boolean degree.

### 6.1. A bounded original exceptional module

Write the actual endpoint band as a product of signed native rotations
on independent commuting axes a_j. The Pauli label
$`q(z)=K_AZ(z)K_A^\dagger`$ is linear in z. Its active-axis bits
$`\lambda_j(z)=\langle a_j,q(z)\rangle`$ are therefore linear.
The literal Clifford J_z has binary action

```math
\Phi_z=I+N_z,\qquad
N_z(w)=\sum_j\lambda_j(z)\langle a_j,w\rangle a_j.
```

Orientation changes affect literal phases but not this binary action.
Isotropy of the actual axis list gives $`N_zN_u=0`$ and
$`\Phi_z^{-1}=\Phi_z`$. In particular N_z is linear in z.

Let E be the set of original middle axes b_i for which some Phi_z sends
b_i outside L+R. For every b_i in E, the map
$`z\mapsto N_zb_i\pmod{L+R}`$ is nonzero linear. Thus b_i escapes
for at least half of the uniformly distributed computational probes z.
Choose d axes from E independent modulo R. For each z their images
remain independent modulo R because Phi_z fixes R and is invertible.
The escaping sublist has at most 8c members by Section 6. Averaging yields

```math
\frac d2\le8c,\qquad
\dim\bigl((\mathrm{span}E+R)/R\bigr)\le16c.
```

This bounds a quotient span, not the number of exceptional physical axes:
central relations may leave many original axes in the same quotient.
Split the actual middle subband into E_bad on these exceptional axes
and E_good on the others. The original insertion has support in H_0/R
of dimension c. Hence the exact regrouping

```math
F'=F E_{\rm bad},\qquad
F E_B=F'E_{\rm good}
```

places F' in a quotient-support space of dimension at most 17c.
F' still commutes every E_good axis. This charges exactly the existing
native gates and their actual inverses; no new synthesis, feature bound,
or clean return is inferred. Good-axis images stay in the **old** L+R,
not necessarily in the span of E_good plus R. Propagating the exceptional
module and matching the surviving axes remain unpaid; iterating this
absorption does not by itself give a bounded number of rounds.

### 6.2. The transformed union can grow in the exported comparison

The distinction is necessary even for an actual native finite-insertion
comparison with a complete Boolean identity. On r arbitrary logical
wires, r at least two, take

```math
E_A=E_C=\prod_{i=1}^r T_{Z_i},\qquad
K_A=K_C=H^{\otimes r},\qquad
B=X_1\cdots X_r,
```

```math
E_B=T_B,\qquad F=T_B^\dagger,\qquad R=0.
```

All three native axis lists are independent. F commutes B and has two
Pauli terms, so its insertion support budget has c=1. The full word
and its actual inverse are literally identity, because the adjacent
F and E_B cancel, followed by the actual endpoint inverses. Every input
and reference is arbitrary.

The actual commuting Clifford family is
$`J_z=E_AX(z)E_A^\dagger`$, with the displacement formula above and
pairwise zero products. Its only middle-axis image has binary label
$`b+Z(z)`$; the Hermitian operator is
$`i^{\mathrm{wt}(z)}BZ(z)`$, retaining the literal phase. Each nonzero
z has one escaping axis. But the images at e_1, e_2 and e_1+e_2 sum
to b in binary labels, and all images together span r+1 dimensions.
Thus the transformed union is unbounded although the original exceptional
space has dimension one.

This is a guard for the exported finite-insertion hypotheses, not an
original physical-pair Boolean counterfamily. No original clean-transition
realization of this separately chosen F is asserted, and its c=1 is not
identified with an original physical-pair quotient. An unused clean
spectator can be returned exactly without adding that missing provenance.
The example does not reopen the no-pair route or disprove a consequence
using the full original-pair structure. It rules out inferring a union
bound from the exported native comparison alone.

The [minimal physical-pair completion](MINIMAL_PHYSICAL_PAIR_COMPLETION.md)
provides actual one-pair provenance for a related occupied identity, but
its cancellation leaves a required centrally controlled partner. For that
fixed input and middle word, adding a silent band that moves the singleton
prevents every one-band output completion to a Clifford target, including
arbitrary output Cliffords and unused-sector extensions. Every complete
target in the same family has a zero-clean three-band word with one
middle axis, giving diagonal third-level closure and Boolean sign/XOR
degree at most three/two, now including all admissible silent input axes
coupled to the active core. The [paired rank-two extension](MINIMAL_PHYSICAL_PAIR_COMPLETION.md#5-all-paired-sources-with-two-clean-active-axes-at-both-endpoints)
also gives sign/XOR degree at most six/five when both actual endpoint
clean-active ranks are two and an original middle pair is present.
This uses common pair-parity coordinates and native quotient independence.
The no-pair branch is separate; general inside-span matching beyond these
endpoint hypotheses remains open.


## 7. Exact checks and scope

The [small exact fixtures](../../tests/test_central_module_profiles.py)
use four independent physical Z axes to distinguish the three profiles:
$`k=4,h=3,r=2`$ gives $`\nu=2,\delta=1,\mu=0`$, although the
raw displacement rank is four. A one-axis case with P anticommuting
with R checks an odd central relation. Each fixture retains literal
native words, their actual inverse, every matrix column and exact
cyclotomic coefficients. They test the norm identities, not an asserted
original-pair Boolean counterfamily.

All reductions retain the original clean and dirty return contract and
external-reference coherence. If emitted, the fixed coordinate Cliffords
and their actual inverses carry their elementary-Clifford costs; no
module norm or conditional expectation is credited as a physical gate.
The complete-frame frontier and the unrestricted one-clean three-layer
question are unchanged.
