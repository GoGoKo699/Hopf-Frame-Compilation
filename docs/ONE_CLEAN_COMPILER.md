# One-clean compilation of Hopf frames

[Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Open endpoint](OPEN_PROBLEM.md)

One initialized flag suffices for an addressed real rotation with a
precision-sized arbitrary dirty core. The construction conjugates one
scalar block by another, routing their rejection terms through different
logical Pauli operators. An anticommutator then cancels the dirty terms in
the accepted block.

The source is an exact Clifford+T circuit. Paired Majorana operators carry
two independently programmed geometric coefficient lists on essentially
one dirty qubit per precision bit. Five forward or inverse block calls
amplify a fixed normalization, with no additional initialized flag.

For the prescribed real frame, this gives $`T=O(N+nL)`$ and, with
conditional-suffix grouping, $`T=O(N+L\ell_*(n))`$, both with
$`G=O(NL)`$, one clean qubit, and $`b\ge L+n+7`$ dirty qubits.
The grouped banked refinement also extends to one clean qubit.
[Section 8](#8-phase-dressed-complex-magnitude-frames) gives the same
count bounds for the phase-dressed complex magnitude frame, with
$`b\ge L+n+8`$, by composing with a literal diagonal. These results
reduce the sufficient clean allocation. [Section 9](#9-a-borrowed-signal-suffices-for-real-rotations)
uses an exact flag symmetry to obtain the layerwise real-frame bound with
zero clean qubits. The stronger grouped bound retains one clean qubit;
the linear-T endpoint remains open.

## 1. Contract and addressed primitive

Let x be an unchanged k-bit address, with $`S=2^k`$ rows, and let t be a
logical target. A Boolean predicate h may depend on other unchanged logical
bits, including fixed address restrictions, but not on t, the signal flag,
or the dirty work. It will be implemented as a conjunction of positive or
negative literals. The prescribed action is

```math
W_{x,h}=\begin{cases}R_y(\theta_x),&h=1,\\ I,&h=0,\end{cases}
\qquad R_y(\theta)=e^{-i\theta Y}.
```

Choose an integer $`q\ge5`$ and put $`m=q+1`$. The circuit uses one
initialized signal flag a, m dirty core qubits, k dirty table selectors,
and one additional dirty helper for the predicate-controlled gates. If p
is the number of predicate literals, it has

```math
\|VJ_1-J_1(W\otimes I_{\rm dirty})\|\lt30\,2^{-q},
\qquad
T=O(S+q+p^2),\qquad G=O(Sq+q+p^2).
```

The inequality is over every logical and dirty input and therefore also
holds with arbitrary reference systems. It includes flag leakage and
approximate return of the core. The selectors and predicate helper return
exactly after each completed subroutine. On $`h=0`$, the actual circuit is
exactly identity on its complete input space, including rejected flag
inputs. All inverses are actual circuit inverses. No measurements, resets,
supplied catalysts, or uncharged quantum oracles are used.

The angles must admit certified classical sine and cosine evaluation.
That evaluation and the construction of the classical tables are
preprocessing costs, excluded from T and G. Every coherent use of a table
is implemented and charged below.

## 2. An exact two-tail source and its dirty masks

On m dirty qubits define the paired Majorana operators

```math
\Gamma_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\Gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j,
\qquad 0\le j\lt m.
```

They are Hermitian involutions and satisfy
$`\{\Gamma_i,\Gamma_j\}=2\delta_{ij}I`$. This is the standard
Clifford-algebra machinery used by the
[operator-source construction](OPERATOR_SOURCE_COMPILER.md#1-the-operator-source-and-its-exact-native-circuit).
For distinct indices define

```math
R_{ij}=\exp\!\left(\frac{i\pi}{8}(i\Gamma_i\Gamma_j)\right),
\qquad
R_{ij}\Gamma_iR_{ij}^\dagger
=\frac{\Gamma_i+\Gamma_j}{\sqrt2}.
```

The rotation commutes with every other Majorana except $`\Gamma_j`$.
Starting from $`\Gamma_0`$, apply $`R_{01}`$, then $`R_{02}`$. The
result has squared coefficients $`1/4,1/2,1/4`$ on
$`\Gamma_0,\Gamma_1,\Gamma_2`$, respectively. Keep the first as a
head coefficient. Split the second along
$`\Gamma_1,\Gamma_3,\ldots,\Gamma_{2q-1}`$, and the third along
$`\Gamma_2,\Gamma_4,\ldots,\Gamma_{2q}`$. Each split uses the
corresponding $`R_{ij}`$ above.

Write

```math
w_j=\begin{cases}2^{-j-1},&0\le j\lt q-1,\\
2^{-(q-1)},&j=q-1.
\end{cases}
\qquad \sum_{j=0}^{q-1}w_j=1.
```

The resulting source is

```math
M=\frac12\Gamma_0
+\sum_{j=0}^{q-1}\sqrt{\frac{w_j}{2}}\,\Gamma_{2j+1}
+\sum_{j=0}^{q-1}\sqrt{\frac{w_j}{4}}\,\Gamma_{2j+2},
\qquad M^\dagger=M,\quad M^2=I.
```

The last available Majorana, $`\Gamma_{2q+1}`$, is unused. If U is
the actual loader word just specified, then $`M=UX_0U^\dagger`$.
There are $`2q`$ splitting rotations. Their Pauli generators all have
weight at most two: the two seed pairs are $`(0,1)`$ and $`(0,2)`$, while
each tail advances by two Majorana indices. For example,

```math
i\Gamma_{2j}\Gamma_{2j+2}=Y_jX_{j+1},\qquad
i\Gamma_{2j+1}\Gamma_{2j+3}=-X_jY_{j+1}.
```

Each rotation uses one T or T-dagger and a constant number of Clifford
gates, up to a common scalar. That scalar cancels in every actual source
or controlled-source call. Thus $`M`$ costs $`O(q)`$ native gates and
has no initialized work requirement.

### Arbitrary signs use two whole-word queries

For desired signs $`(-1)^{f_i}`$ on all $`2m`$ Majoranas, put

```math
x_j=f_{2j}\mathbin\oplus f_{2j+1},\qquad
z_j=f_{2j}\mathbin\oplus\bigoplus_{r\lt j}x_r,
\qquad P_f=X^xZ^z.
```

Conjugation by $`P_f`$ changes precisely the requested signs. In
particular, define

```math
N_f=P_fMP_f^\dagger,\qquad
\frac{MN_f+N_fM}{2}=c_f I,
\qquad c_f=\sum_i a_i^2(-1)^{f_i},
```

where $`a_i`$ are the displayed coefficients of M. The unused sign can
be fixed arbitrarily. The mask need not be Hermitian: its actual inverse
is essential. Any scalar belonging to an address-dependent Pauli word
cancels in $`P_fMP_f^\dagger`$.

An addressed $`X^x`$ mask is the existing exact dirty whole-word XOR
query on the m core bits. An addressed $`Z^z`$ mask is another such
query conjugated by m Hadamards. Their selectors return exactly on all
inputs. Consequently one mask costs $`O(S)`$ T gates and
$`O(Sm+m)`$ Clifford gates, using k additional dirty selectors. The
fixed mask used below only needs ordinary Pauli gates.

## 3. Conjugating scalar blocks produces a rotation

For any programmed mask let

```math
\mathcal S_f
=H_a C_{a=0}(M)N_f C_{a=1}(M)H_a
=c_f I+X_a\otimes D_f,
\qquad D_f=\frac{MN_f-N_fM}{2}.
```

This is the same full-unitary scalar word as in the operator-source
chapter, now using the paired source. Its identities hold on the entire
dirty Hilbert space. In particular,

```math
D_f^\dagger=-D_f,\qquad
D_f^\dagger D_f=(1-c_f^2)I.
```

Fix f to flip exactly the even tail $`\Gamma_2,\Gamma_4,\ldots,\Gamma_{2q}`$.
Its squared weight is $`1/4`$, so $`c_f=c=1/2`$.
For a programmable mask g, write

```math
s=\sum_i a_i^2(-1)^{g_i},\qquad
\tau=\sum_i a_i^2(-1)^{f_i+g_i}.
```

The two rejected operators obey the scalar anticommutator identity

```math
D_fD_g+D_gD_f=2(cs-\tau)I.
```

To verify it, write $`N_f=cM+F`$ and $`N_g=sM+G`$, with F and G
Clifford vectors orthogonal to M. Then $`D_f=MF`$,
$`D_g=MG`$, and $`\{F,G\}=2(\tau-cs)I`$; anticommutation with M
gives the formula.

Conjugate the fixed scalar word by $`\mathrm{CZ}_{a,t}`$ and the
programmable word by $`\mathrm{CNOT}_{a\to t}`$. This gives actual
unitaries

```math
A=cI+X_aZ_tD_f,\qquad
B_g=sI+X_aX_tD_g,
\qquad Q=A^\dagger B_gA.
```

Both routing conjugations are Clifford. In the flag-zero compression,
the terms containing one or three flag flips vanish. The two remaining
mixed terms combine by the anticommutator above:

```math
\langle0|_aQ|0\rangle_a
=sI+2c(cs-\tau)XZ
=sI+\left(\frac{s}{2}-\tau\right)XZ.
```

The dirty factor in this block is exactly identity. The full rejected
action is supplied by the word $`A^\dagger B_gA`$; no accepted blocks
were multiplied while discarding their rejected components.

### Programming the two coefficients

Put $`\beta=\sin(\pi/10)=(\sqrt5-1)/4`$. To obtain
$`\beta R_y(\theta)`$, the desired moments are

```math
s=\beta\cos\theta,\qquad
\tau=\frac{s}{2}-\beta\sin\theta.
```

Let $`p_+`$ be the signed coefficient sum on the head and odd tail,
and let $`p_-`$ be that on the even tail. They satisfy
$`s=p_++p_-`$ and $`\tau=p_+-p_-`$. Thus the target sums are

```math
p_+=\beta\left(\frac34\cos\theta-\frac12\sin\theta\right),
\qquad
p_-=\beta\left(\frac14\cos\theta+\frac12\sin\theta\right).
```

Choose the head sign $`\sigma\in\{-1,1\}`$ and tail means u and v
so that

```math
p_+=\frac\sigma4+\frac u2,\qquad
p_-=\frac v4.
```

This can be done by finite certified arithmetic. Compute a rational
estimate of $`p_+`$ with error at most $`1/16`$, and choose
$`\sigma=1`$ when the estimate is nonnegative, otherwise
$`\sigma=-1`$. Set $`u=2p_+-\sigma/2`$ and $`v=4p_-`$.
The bounds

```math
|p_+|\le\frac{\beta\sqrt{13}}4\lt\frac3{10},\qquad
|u|\le\frac58,\qquad |v|\le\beta\sqrt5\lt1
```

show that both tail means are available. If the estimate has the wrong
sign, the true value has magnitude at most $`1/16`$, which proves the
same u bound. No exact sign or zero test is needed.

Apply the existing
[certified geometric digit rule](OPERATOR_SOURCE_COMPILER.md#2-exact-signed-dyadic-coefficients)
to u and v separately, using q signs for each tail. More explicitly,
put $`e=2^{1-q}`$, approximate the mean to error at most $`e/4`$,
clamp the estimate to $`[-1,1]`$, and round it to the sign grid with
spacing $`2e`$. This uses rational arithmetic with a fixed tie rule and
gives mean error at most $`5e/4`$.

Consequently the errors in the two coefficient sums obey

```math
|\Delta p_+|\le\frac54\,2^{-q},\qquad
|\Delta p_-|\le\frac58\,2^{-q}.
```

Since the scalar and $`XZ`$ coefficients are respectively
$`p_++p_-`$ and $`-p_+/2+3p_-/2`$, their operator error is

```math
\delta=\left\|\langle0|Q|0\rangle-\beta R_y(\theta)\right\|
\le\sqrt{(\Delta p_++\Delta p_-)^2+
\left(-\frac{\Delta p_+}{2}+\frac{3\Delta p_-}{2}\right)^2}
\lt\frac52\,2^{-q}.
```

The same bound holds for the coherent direct sum over all addresses.
The source and masks are exact native circuits; this coefficient rounding
is the only approximation in the construction.

## 4. Five calls amplify with one flag

Let Q be any actual unitary, J its initialized embedding, W a target
unitary, and $`B=J^\dagger QJ`$. Suppose
$`\delta=\|B-\beta W\|\le\beta/2`$. With
$`R=I-2JJ^\dagger`$, define the actual five-call word

```math
\mathcal A_5=QRQ^\dagger RQRQ^\dagger RQ.
```

The usual oblivious-amplification multiplication gives

```math
J^\dagger\mathcal A_5J
=5B-20BB^\dagger B+16(BB^\dagger)^2B.
```

This is two Grover iterations, using the same established
[oblivious-amplification framework](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error).
The following estimate records the complete output error for this
normalization. Write the polar decomposition $`B=VH`$. Then

```math
\|H-\beta I\|\le\delta,\qquad
\|V-W\|\le\frac{2\delta}{\beta}.
```

For $`p_5(x)=5x-20x^3+16x^5`$, the accepted block is
$`Vp_5(H)`$. The choice of $`\beta`$ gives
$`p_5(\beta)=1`$ and $`p_5'(\beta)=0`$. Moreover,
$`|p_5''(x)|\le30`$ on $`[0,1/2]`$, which contains the spectrum of
H under the stated error bound. Taylor's theorem gives
$`\|I-p_5(H)\|\le15\delta^2`$. Unitarity of the actual word now
controls rejected-space leakage as well:

```math
\|\mathcal A_5J-JV\|^2
=2\|I-p_5(H)\|\le30\delta^2.
```

Therefore

```math
\|\mathcal A_5J-JW\|
\le\left(\sqrt{30}+\frac2\beta\right)\delta
\lt12\delta.
```

For one initialized flag, $`R=-Z_a`$. There are four reflections, so
replacing every R by the physical $`Z_a`$ preserves the literal scalar.
The implementation uses three forward Q calls, two actual inverses, and
four Clifford Z gates. It needs no additional clean flag or arbitrary
amplification phase. Combining this lemma with the digit error proves the
$`30\,2^{-q}`$ primitive bound, including the entire dirty space.

## 5. Exact conditioning and the resource ledger

Condition the scalar word on an unchanged predicate h by adding h only
to each source's central X gate. For example,

```math
C_h(M)=U C_h(X_0)U^\dagger,\qquad
C_{h,a=v}(M)=U C_{h,a=v}(X_0)U^\dagger,
\qquad C_h(N_g)=P_gC_h(M)P_g^\dagger.
```

The surrounding loader and mask circuits are unconditional. When h is
false, they cancel with their actual inverses. Keep both outer Hadamards
on a in the scalar word; they then also cancel on the inactive sector.
This constructs exactly

```math
\mathcal S_{g,h}=\Pi_h\mathcal S_g+(I-\Pi_h)I.
```

Here $`\Pi_h`$ is the logical predicate projector. It commutes with
the source, address queries, and routing conjugations by the disjointness
requirements in Section 1. Condition both fixed and programmable scalar
words this way. Their routed conjugation is exactly Q on the active
sector and identity on the inactive sector. The latter remains exactly
identity after five-call amplification: the four unconditioned flag
reflections multiply to identity there. This statement holds even when
the flag and every dirty wire have arbitrary inputs.

Apply the robust amplification estimate on the invariant active sector.
The inactive block is I, not $`\beta I`$; its error is zero by the exact
identity just proved. The direct-sum error is therefore the active-sector
bound, with no amplification hypothesis imposed on the inactive block.

An h-controlled central X, with or without the additional scalar-flag
control, uses the retained exact borrowed-MCX construction. It costs
$`O(p^2+1)`$ Toffolis and restores one separate arbitrary helper. Its
controls are distinct from the core target and helper. No source loader
is itself promoted to a string of controlled T gates. The query
selectors return before the next source or predicate subroutine.

There are three source calls in a scalar word and three scalar words in
$`A^\dagger B_gA`$. Only the middle scalar word has addressed masks;
the fixed word's masks are Pauli circuits. Five-call amplification
multiplies these constant counts by five. The source is $`O(q)`$, a
mask query is $`O(S)`$ T and $`O(Sq)`$ Clifford, and each central
predicate gate is $`O(p^2+1)`$. This proves the stated bounds with
exactly

```math
a=1,\qquad b=(q+1)+k+1=q+k+2
```

allocated work qubits. The logical target and all predicate controls
remain data throughout. The one initialized flag is reused coherently;
no return or reinitialization is assumed between subroutine calls.

### Literal diagonal phases use the same flag

There is a flag-only variant that does not need a logical target. Keep
$`A=\mathcal S_f=cI+X_aD_f`$ and conjugate the programmable scalar
word by the standard Clifford $`S_a=\mathrm{diag}(1,i)`$. Since
$`S_aX_aS_a^\dagger=Y_a`$, the actual word

```math
Q_{\rm phase}=A^\dagger(S_a\mathcal S_gS_a^\dagger)A
```

has accepted block

```math
\langle0|Q_{\rm phase}|0\rangle
=s-2ic(cs-\tau)=s+i\left(\tau-\frac s2\right).
```

Program $`s=\beta\cos\phi`$ and
$`\tau=s/2+\beta\sin\phi`$, reversing the sine sign in the preceding
coefficient rule. The same precision, amplification, conditioning, and
workspace arguments compile the addressed scalar $`e^{i\phi_x}`$.
Its literal phase is retained: the S gates are actual inverses, source
scalars cancel within their conjugations, and the four amplification
reflection signs cancel exactly. This is an addressed diagonal-unitary
primitive, not a claim that arbitrary complex-frame grouping follows
without an additional composition proof.

## 6. Complete real frames and the grouped extension

Let $`N=2^n`$, $`n\ge1`$, $`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. Suppose

```math
a\ge1,\qquad b\ge L+n+7.
```

At Hopf depth d, use the d-bit prefix as the table address, the prescribed
rotation bit as target, and the remaining zero-suffix condition as h.
Set

```math
q_d=L+n-d+5,
\qquad
m_d=q_d+1.
```

The core, selectors, and predicate helper occupy exactly

```math
m_d+d+1=L+n+7
```

dirty wires. The same pool is repartitioned between these roles at each
depth. No logical suffix bit is counted as initialized work. On every
inactive suffix, the complete native stage is exactly identity. On the
active sector, its error is less than $`30\,2^{-q_d}`$, uniformly over
all addresses and dirty inputs. Hence

```math
\sum_{d=0}^{n-1}30\,2^{-q_d}
=\frac{30}{32}\,2^{-L}\sum_{d=0}^{n-1}2^{d-n}
\lt\frac{15}{16}\,2^{-L}\le\eta.
```

Compose the actual stages in the prescribed order. A full-isometry hybrid
uses the ideal zero flag and returned dirty core only as comparison
states; earlier actual leakage is propagated by unitaries. Thus

```math
\|VJ_a-J_a(W\otimes I_b)\|\le\eta.
```

Only one clean qubit is used. The error includes the flag, core, and
external references; it is not merely an accepted-block estimate.
For each depth the native counts are

```math
T_d=O(2^d+q_d+n^2),\qquad
G_d=O(2^dq_d+q_d+n^2).
```

The weighted geometric sums and $`n^3=O(2^n)`$ give

```math
T=O(N+nL),\qquad G=O(NL).
```

The [grouped proof](CONDITIONAL_SUFFIX_COMPILER.md#10-the-grouped-bounds-need-only-one-external-clean-qubit)
improves this to

```math
T=O\bigl(N+L\ell_*(n)\bigr),\qquad G=O(NL),
\qquad \ell_*(n)=1+\log_2^*(n+2).
```

It keeps the sole external clean qubit as the active-suffix predicate and
moves the former external scalar flag into the known-zero active suffix.
The private reservation grows by one bit. The present primitive compiles
the fixed deepest tail, where that private work is unavailable. The same
dirty threshold suffices, with a complete error budget proved there.
If $`b\ge2(L+n+7)`$, its banked form gives

```math
T=O\!\left(\sqrt{NL}+L\ell_*(n)+\frac{NL}{b}\right),
\qquad G=O(NL),\qquad a\ge1.
```

The banked matching regimes extend to the fixed allocation $`a=1`$.
The existing T-depth and simultaneous count/depth schedules retain their
separately proved two-clean assumptions. At $`L=N`$ and $`b=N+n+7`$,
the one-clean upper bound is still $`O(N\ell_*(n))`$. Every frame column
is preserved within the displayed complete-input norm, so the existing
[fixed-parameter QBP guarantee](QBP_APPROXIMATION.md) applies unchanged.

## 7. Literal diagonals and complete one-target multiplexors

The flag-only phase word in Section 5 implements a complete literal
diagonal on k address qubits. Write $`M=2^k`$ for its number of entries.
At error $`\eta`$, take $`q=L+5`$. Its dirty reservation is
$`q+k+2=L+k+7`$ and its full-isometry error is below
$`(15/16)2^{-L}`$. Thus one clean qubit gives

```math
T=O(M+L),\qquad G=O(ML),\qquad b\ge L+k+7.
```

With $`b\ge2(L+k+7)`$, reserve the base core/selector/helper pool and
use the remainder for the existing exact dirty word-bank query. The X-
and Z-mask tables are queried sequentially, so they reuse the same banks.
Their outputs are the dirty core, and all banks and selectors return
exactly. A constant number of such queries and source calls gives

```math
T=O\!\left(\sqrt{ML}+L+\frac{ML}{b}\right),\qquad G=O(ML).
```

For a complete one-target U(2) multiplexor, use a literal scalar phase
and three addressed Euler rotations. Rz factors follow from the real
rotation primitive by fixed Clifford conjugation. The
[certified Euler procedure](OPERATOR_SOURCE_COMPILER.md#81-general-one-qubit-multiplexors-with-two-clean-qubits)
chooses their product to approximate each supplied block within
$`2^{-L}/16`$, including singular cases and its literal scalar phase.
For each of the four factors use $`q=L+7`$; their total synthesis error
is below $`(15/16)2^{-L}`$. Full-isometry composition therefore gives
error at most $`\eta`$, with

```math
a=1,\quad b\ge L+k+9,\qquad
T=O(M+L),\qquad G=O(ML).
```

At $`b\ge2(L+k+9)`$ the same sequential bank reuse yields
$`T=O(\sqrt{ML}+L+ML/b)`$ and $`G=O(ML)`$.
For $`k\ge1`$, the retained diagonal and multiplexor lower bounds also
allow one clean qubit, so these banked one-stage bounds are matching in
their stated regimes. The upper bounds also cover $`k=0`$.
Certified Euler preprocessing has no asserted uniform classical
running-time bound for arbitrary input evaluators. This corollary does
not assert an optimal T-depth or remove the repeated precision charges
in a general noncommuting frame product.

## 8. Phase-dressed complex magnitude frames

The [prescribed complex magnitude frame](HOPF_INTERFACE.md#6-phase-dressed-complex-magnitude-frame)
is

```math
W_{\mathbb C,\mathrm{mag}}=D_\phi W_{\mathbb R},
\qquad
D_\phi=\sum_{x=0}^{N-1}e^{i\phi_x}|x\rangle\langle x|.
```

The real Hopf angles and leaf phases are supplied independently and admit
certified classical evaluation. Under the same assumptions on n, N,
$`\eta`$, and L as Section 6, one clean qubit suffices for

```math
b\ge L+n+8,
\qquad T=O\bigl(N+L\ell_*(n)\bigr),\qquad G=O(NL).
```

If $`b\ge2(L+n+8)`$, the banked bound is

```math
T=O\!\left(\sqrt{NL}+L\ell_*(n)+\frac{NL}{b}\right),
\qquad G=O(NL).
```

Both statements have the complete initialized-isometry contract

```math
\|VJ_a-J_a(W_{\mathbb C,\mathrm{mag}}\otimes I_b)\|\le\eta,
\qquad a\ge1.
```

*Proof.* Allocate error $`\eta/2`$ to each factor. Its precision parameter
is $`L'=L+1`$. Compile the real frame by Section 6 and the literal
diagonal by Section 7 with $`k=n`$. Each uses one initialized flag and
fits in $`L'+n+7=L+n+8`$ dirty wires. Execute the real circuit
$`V_R`$ followed by the diagonal circuit $`V_D`$, reusing exactly the
same flag and dirty pool. Extra clean qubits, if available, are untouched.

To justify this reuse despite approximate work return, write
$`R=W_{\mathbb R}\otimes I_b`$ and $`D=D_\phi\otimes I_b`$.
The actual circuits are unitaries on the common workspace, so

```math
V_DV_RJ_a-J_aDR
=V_D(V_RJ_a-J_aR)+(V_DJ_a-J_aD)R.
```

Taking norms gives the sum of the two full-isometry errors. In
particular, the proved real-frame and diagonal estimates give

```math
\|V_DV_RJ_a-J_aDR\|
\lt\left(1+\frac{15}{16}\right)2^{-L'}
=\frac{31}{32}\,2^{-L}\le\eta.
```

Earlier actual flag leakage and core disturbance are propagated by
$`V_D`$; the proof uses the zero flag and returned dirty input only in
the ideal comparison term. No intermediate reset or additional initialized
work is needed. The estimate remains valid with arbitrary external
references. The literal diagonal primitive preserves every relative phase
and any common phase of $`D_\phi`$.

The real and diagonal counts add, while their workspace reservations take
the maximum. Section 7 contributes $`O(N+L')`$ T gates unbanked, or
$`O(\sqrt{NL'}+L'+NL'/b)`$ with
$`b\ge2(L'+n+7)`$. Both contributions are absorbed by the respective
grouped real-frame bounds, since $`L'=\Theta(L)`$ and
$`\ell_*(n)\ge1`$. Their Clifford counts sum to $`O(NL)`$. ∎

For fixed $`a=1`$, the banked worst-case bound is matching under the
same sufficient regimes $`L\ell_*(n)^2\le N`$ or
$`b\le N/\ell_*(n)`$, subject to its stated bank threshold. Indeed,
setting all real Hopf angles to zero gives $`W_{\mathbb R}=I`$ and
leaves the arbitrary diagonal family. The existing diagonal lower bound
and fixed-width count therefore apply with total width
$`n+1+b=\Theta(b)`$. In either sufficient regime the extra
$`L\ell_*(n)`$ term is absorbed by $`\sqrt{NL}`$ or $`NL/b`$;
the shift from L to $`L+1`$ changes only an absolute constant.

At a fixed parameter tuple, this complete-frame contract supplies the
[same QBP magnitude-stream guarantee](QBP_APPROXIMATION.md#6-parameter-derivatives-and-complex-coordinates).
The leaf-phase derivatives still use their separate direct signed one-hot
stream; they are not additional frame columns. The result concerns the
specified phase-dressed magnitude frame, not arbitrary complex unitaries.
It does not extend the separately proved T-depth schedules to one clean
qubit or close the linear-T endpoint.

## 9. A borrowed signal suffices for real rotations

The real-rotation word has an additional symmetry on its complete flag
space. This strengthens the initialized-column estimate in Section 4 to
a full-unitary estimate, so that its amplification flag may itself be an
arbitrary dirty qubit. The literal phase word does not have the same
symmetry; its initialized-flag assumption is retained.

**Lemma — extension by flag symmetry.** Let r be a qubit, let
$`J_0|\psi\rangle=|0\rangle_r|\psi\rangle`$, and let U act on all
other registers. If an actual unitary V satisfies

```math
[V,X_r]=0,\qquad
\|VJ_0-J_0U\|\le\epsilon,
```

then

```math
\|V-I_r\otimes U\|\le\sqrt2\,\epsilon.
```

*Proof.* Put $`\Delta=V-I_r\otimes U`$,
$`C=\Delta J_0`$, and $`J_1=X_rJ_0`$. Commutation gives
$`\Delta J_1=X_rC`$. For arbitrary vectors u and v,

```math
\begin{aligned}
\|\Delta(J_0u+J_1v)\|
&=\|Cu+X_rCv\|\\
&\le\epsilon(\|u\|+\|v\|)
\le\sqrt2\,\epsilon\sqrt{\|u\|^2+\|v\|^2}.
\end{aligned}
```

The orthogonal flag sectors exhaust the input space, proving the claim.
The same estimate holds after tensoring arbitrary reference systems. ∎

For the real primitive of Section 3, both

```math
A=cI+X_rZ_tD_f,\qquad
B_g=sI+X_rX_tD_g
```

commute with $`X_r`$, and therefore so does $`Q=A^\dagger B_gA`$.
The amplification reflection is $`-Z_r`$. Conjugating its five-call word
by $`X_r`$ changes the sign of each of its four reflections, leaving the
complete word invariant. This is an exact circuit identity, independent
of coefficient-rounding error.

An unchanged logical predicate preserves the identity: each conditioned
scalar word is the same scalar-source word on the active sector and
exact identity on the inactive sector. Address direct sums also preserve
the symmetry. Thus Section 4 and the lemma give the addressed real
rotation, with an arbitrary flag r, the full-space guarantee

```math
\|V-W_{x,h}\otimes I_{\rm dirty}\|
\lt30\sqrt2\,2^{-q}\lt43\,2^{-q}.
```

Here the dirty identity includes r, the precision core, the selectors,
and the predicate helper. The complete dirty reservation is now

```math
b_{\rm primitive}=(q+1)+k+1+1=q+k+3,
```

with no initialized qubit. The counts remain
$`T=O(2^k+q+p^2)`$ and $`G=O(2^kq+q+p^2)`$.
Core and flag return are approximate in this full-unitary norm; selectors
and predicate helpers return exactly. On an inactive predicate the
complete circuit remains exactly identity. Fixed target Clifford
conjugations give the same result for Rz rotations. Actual inverse
circuits satisfy the same full-space error bound for the inverse target.

In contrast, the literal phase construction replaces one routed flag
Pauli by $`Y_r`$. The displayed $`X_r`$ symmetry is then unavailable.
This lemma does not remove the initialized flag from Sections 7 or 8.

### A zero-clean layerwise real-frame compiler

For the prescribed complete real Hopf frame, the strengthened primitive
gives

```math
a=0,\qquad b\ge L+n+7,\qquad
\|V-W\otimes I_b\|\le\eta,
\qquad T=O(N+nL),\qquad G=O(NL).
```

The dirty threshold is preserved by a fixed address split, rather than
adding the borrowed flag to the earlier reservation. At a tree depth
$`d\ge2`$, choose

```math
q_d=L+n-d+6.
```

Select two of the d address bits as fixed sector literals and use the
remaining $`k=d-2`$ bits as the free table address. Process the four
sectors sequentially. Every sector uses exactly

```math
q_d+k+3=L+n+7
```

dirty qubits. Its address bits are unchanged, and on each other sector
its actual action is identity on the entire flag and core space.
Consequently the error over these four invariant sectors is their
maximum, even with arbitrary correlated dirty inputs. It is less than
$`43\,2^{-q_d}`$ for the complete depth layer.

The depths $`d=0,1`$, when present, contain only three addressed rows
in total. Compile them by the direct native
[borrowed-sector echo](BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector),
using no initialized work and exactly returned borrowed helpers. Choose
each row's error at most $`2^{-L}/16`$. Disjoint invariant row sectors
give that same bound per shallow layer, hence at most $`2^{-L}/8`$
in total. These finitely many words cost $`O(L+n^2)`$ gates, including
their logical predicates, and fit the stated dirty pool. This also
covers $`n=1,2`$ without using a negative address length.

Full-unitary telescoping over all depths gives

```math
\begin{aligned}
\|V-W\otimes I_b\|
&\lt \frac18\,2^{-L}
 +43\sum_{d=2}^{n-1}2^{-(L+n-d+6)}\\
&\le\left(\frac18+\frac{43}{64}\right)2^{-L}
=\frac{51}{64}\,2^{-L}\le\eta,
\end{aligned}
```

where an empty sum is zero. The comparison includes disturbance of every
borrowed wire and arbitrary reference correlations; no intermediate work
register is assumed freshly initialized. A constant number of sectors
changes only absolute constants in the resource sum. The deep-layer
counts are

```math
T_d=O(2^d+L+n-d+n^2),\qquad
G_d=O\bigl(2^d(L+n-d)+L+n-d+n^2\bigr).
```

Summing them and the shallow words yields
$`T=O(N+nL+n^3)=O(N+nL)`$ and $`G=O(NL)`$, since
$`n^3=O(2^n)`$ and $`L\ge6`$. The complete prescribed frame is the
same target throughout, so its existing fixed-parameter QBP guarantee
continues to apply.

At $`L=N`$ this zero-clean construction has $`T=O(N\log N)`$.
The stronger grouped count bound remains proved with one clean qubit;
its initialized active-suffix predicate has not been removed by this
symmetry argument. No optimality, new lower bound, grouped zero-clean
extension, or T-depth improvement follows from this corollary.
