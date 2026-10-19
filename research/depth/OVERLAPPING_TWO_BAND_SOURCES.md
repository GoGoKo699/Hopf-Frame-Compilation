# Overlapping two-band sources and transported output probes

[Current source interface](JOINT_FOURIER_WORD_LIMITS.md#5-a-precise-next-interface-for-two-source-bands) ·
[Complete resource, error and return recurrence](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#61-exact-flatness-sharpens-the-complete-error-and-return-recurrence) ·
[Fixed-clean two-band source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) ·
[Two-layer rigidity](TWO_LAYER_DIAGONAL_RIGIDITY.md)

This is a source-specific continuation, not an unrestricted depth theorem. Every nonsignal wire, including output-bank and private dirty-work wires, is arbitrary; identities extend by identity on external references. Physical bands retain independent native Pauli axes and actual inverses.

The [fixed Clifford-sum gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) now
excludes the general fixed-clean, fixed-amplitude two-band source
interface at sufficiently small fixed error and sufficiently large
address size. It requires no flatness assumption. The subclass screens
below retain their stronger finite-size thresholds and explicit
constants, and the displayed recovery recurrences remain conditional
operator identities.

Throughout, J inserts the sole signal in zero, and the actual source is
$`U=U_B^\dagger U_A`$. Its clean-image reflections are
$`A=U_AZ_aU_A^\dagger`$ and $`B=U_BZ_aU_B^\dagger`$.
The target I_p is the complete m-output indicator, with $`m=2^p`$,
including identity on every private dirty-work input.

## 1. A complete transition for parity-coupled overlapping sources

Let $`C_r`$ add the computational parity $`r`$ of the entire logical register into the sole signal a, and similarly define $`C_s`$. Let V_A,V_B be actual one-band logical encoders on the same, potentially overlapping, logical register. Define

```math
U_A=V_A C_r,\quad U_B=V_B H_a C_s,\quad W=V_B^\dagger V_A,
\quad K=J^\dagger U_B^\dagger U_AJ.
```

With $`P=Z_r`$, $`Q=Z_s`$, and $`\Pi_\pm=(I\pm P)/2`$, contraction of the actual signal gives on every logical column

```math
L:=\sqrt2 K=\Phi(W),\qquad
\Phi(E)=E\Pi_+ + QE\Pi_-.
```

Equivalently, $`L_{y,x}=(-1)^{s(y)r(x)}W_{y,x}`$. No disjoint-support assumption occurs. The clean-image reflections are $`A=Z_a(V_A P V_A^\dagger)`$ and $`B=X_a(V_B Q V_B^\dagger)`$. Hence exact flatness holds precisely when

```math
[V_A P V_A^\dagger,V_B Q V_B^\dagger]=0,
```

or equivalently $`[P,W^\dagger QW]=0`$. In that case

```math
L=W\bigl(\Pi_+ + W^\dagger QW\Pi_-\bigr).
```

The second factor is a complete controlled operation, not an uncharged primitive. Its exact implementation or the paid signal-S recovery is required. It cannot be inferred from the reflection ranks.

### A robust screen that does not require flatness

The linear transformation Phi is involutive. For any vector v, use the orthogonal input decomposition to obtain

```math
\|\Phi(E)v\|\le\|E\|(\|\Pi_+v\|+\|\Pi_-v\|)
\le\sqrt2\|E\|\|v\|.
```

If the complete block obeys $`\|K-I_p/\sqrt2\|\le\epsilon`$, put $`W_*=\Phi(I_p)`$. Then $`W_*`$ is a signed permutation with exactly the permutation of the indicator, and

```math
\|W-W_*\|\le 2\epsilon.
```

For any original output-bank wire j, its computational Z observable therefore satisfies

```math
\|W^\dagger Z_jW-(-1)^{y_j+[x=j]}\|\le4\epsilon.
```

The left circuit is actually

```math
V_A^\dagger C_jV_A,\qquad C_j=V_BZ_jV_B^\dagger.
```

$`C_j`$ is Clifford because V_B is an actual one-band encoder. It can be synthesized from its actual tableau at quadratic Clifford cost, retaining literal phase. Thus the displayed circuit has two physical T bands, zero clean work, and arbitrary dirty inputs. The target is a literal Boolean diagonal sign of degree p. The retained [two-layer diagonal rigidity proof and its robust rounding](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction) give error floor $`\delta=\sin(\pi/8)/(2\sqrt2)`$ whenever p is at least four. Consequently this entire parity-coupled source subclass obeys

```math
\|K-I_p/\sqrt2\|\ge\delta/4
=\frac{\sin(\pi/8)}{8\sqrt2},\qquad p\ge4.
```

This floor is independent of width, axis overlap, source reflection rank, or anticommutator rounding. It leaves general signal-mixed one-band encoders outside the screen. It also does not settle the small case p equal to three of this subclass.

Here the scalar-Boolean version of the retained floor uses the same
proof as its XOR corollary: preserved-Pauli rounding first makes the
actual circuit diagonal, the two-band theorem makes it C3, and the
eighth-root phase gap excludes a Boolean sign of degree greater than
three. The analogous one-clean version below also rounds clean return
before applying occupied diagonal rigidity. No basis-probability
approximation is substituted. A common phase on the source cancels in
each transported probe, so optimizing that phase evades neither screen.

## 2. Transported output probe for general encoders

Now let U_A,U_B be arbitrary actual one-band encoders, including signal-mixed dynamics, and assume exact flatness, equivalently $`\{A,B\}=0`$. Let $`L=\sqrt2 K`$. For any original logical computational output Pauli Z_j, define

```math
C=U_B(I_a\otimes Z_j)U_B^\dagger,\qquad
C_\pm=\frac{C\pm ACA}{2},\qquad
T_j=C_+ + C_-B.
```

C is a Hermitian Clifford involution and commutes with B. Because A anticommutes with B, ACA also commutes with B. Thus both C_± commute with B; C_+ commutes with A and C_- anticommutes with A. In addition,

```math
C_+^2+C_-^2=I,\qquad \{C_+,C_-\}=0.
```

It follows that T_j is a Hermitian unitary and $`[T_j,A]=0`$. Its occupied transition is exactly

```math
L^\dagger Z_jL
=J^\dagger U_A^\dagger T_jU_AJ.                 \qquad\text{(P1)}
```

To verify (P1), start from $`2K^\dagger Z_jK=J^\dagger U_A^\dagger(I+B)CU_AJ`$ and insert the A-positive projector on both sides. The terms C_- and BC_+ vanish between these projectors, leaving C_++BC_-=T_j. This calculation includes every logical/dirty column and does not assume a computational basis for A.

### Scalar commutation gives a paid two-band probe

If $`ACA=C`$, then $`T_j=C`$. If $`ACA=-C`$, then $`T_j=CB`$. In both cases T_j is an actual Clifford, including its literal sign. After synthesizing that Clifford from its tableau, the native word $`U_A^\dagger T_jU_A`$ has at most two T bands and returns the sole signal exactly. Its full initialized logical action is $`L^\dagger Z_jL`$. It uses the same original dirty inputs and no new clean history.

For $`\|K-I_p/\sqrt2\|\le\epsilon\le1/10`$, the retained Clifford-gap argument forces exact flatness. Then $`\|L-I_p\|\le\sqrt2\epsilon`$, so

```math
\|L^\dagger Z_jL-(-1)^{y_j+[x=j]}\|
\le2\sqrt2\epsilon.
```

For p at least four, the retained [one-clean two-band real-diagonal floor](TWO_LAYER_DIAGONAL_RIGIDITY.md#4-uniform-robust-boolean-query-obstruction) yields a contradiction if $`2\sqrt2\epsilon\lt\delta`$. Thus an accurate source below $`\delta/(2\sqrt2)=\sin(\pi/8)/8`$ must, for every original output Z_j, have $`ACA`$ different from both C and -C. The same reasoning can be applied to any original output parity whose target conjugate has degree greater than three.

This is a necessary non-scalar-commutator condition on a successful general two-band source. It is not a lower bound on unrestricted indicator depth.

### What the four-term probe does not implement

In the general case

```math
T_j=\tfrac12(C+ACA+CB-ACAB)
```

is a sum of four literal Clifford operators. The sum is a unitary by the algebra above, but this does not make it a Clifford or give it a two-band native realization. It is not permissible to replace it by a supplied operation, to bound it only by its symplectic rank, or to count the number of summands as native depth. The known recovered source already implements its occupied effect with paid constant many source calls; that fact does not yield a low-band classification. The former next question asked whether simultaneous indicator probes force such a classification or permit a successful general source. The [direct projected-defect theorem](FIXED_CLIFFORD_SUM_SOURCE_GAP.md#4-the-direct-projected-defect) supersedes that resumption question: a fixed-length Clifford-sum gap excludes the complete scalable source interface without classifying T_j or treating it as a native operation.

No claim here establishes that Booleanity forces scalar commutation of C and A. That would be a new implication requiring proof, not a consequence of the displayed restriction to the A-positive sector.

Nor is mutual commutation of the output residuals a missing condition.
For any C commuting B, the same map is

```math
\mathcal T(C)=C P_{B,+}+ACA P_{B,-},\qquad
P_{B,\pm}=(I\pm B)/2.
```

It preserves products and adjoints on the algebra commuting B: the
cross-sector terms vanish and the minus-sector product is
$`AC_1AA C_2A=AC_1C_2A`$. Thus the original output probes already
give mutually commuting residuals. The formerly unpaid condition was
their specified high-degree Boolean action in the actual input
coordinates, jointly with the rest of the complete indicator operator.
The direct source gap uses those original coordinates to exclude the
scalable complete transition without an additional commutation claim.

## 3. The four-term unitary is not generally Clifford: an actual source

Take one logical bit, $`r`$ that bit and $`s=0`$. Set $`V_A=TH`$, $`V_B=HTH`$, $`U_A=V_A\mathrm{CX}(\text{logical},a)`$ and $`U_B=V_B H_a`$. These are actual one-band encoders. They have

```math
A=Z_a(X+Y)/\sqrt2,\quad B=X_a,\quad
C=(Z-Y)/\sqrt2,\quad ACA=-(X+Z)/\sqrt2.
```

The source is exactly flat because A anticommutes with B. In the X_a basis its transported probe is

```math
T_j=|+\rangle\langle+|_a\otimes C+
|-\rangle\langle-|_a\otimes ACA.
```

But $`C(ACA)=(-I+iX-iY-iZ)/2`$ is not proportional to a Pauli. Accordingly $`T_j Z_a T_j`$ has non-Pauli off-diagonal logical blocks in this signal basis and is not a Pauli; T_j is not Clifford. The complete normalized logical transition here is simply $`V_B^\dagger V_A`$, so this small example offers no indicator advantage. It certifies that the generic four-term probe cannot be silently replaced by a Clifford even for actual exactly flat encoders.


## 4. A growing signal-mixing family with its complete transition

This family realizes the paused nonzero-product geometry with actual one-band encoders. It is a substantive example, not a faster complete indicator. Its literal normalized transition is a square root of a quadratic Clifford involution. Section 4.3 retains the original square-based rejection of that literal transition. Section 4.4 repairs the square with actual Clifford dressing; the stronger coefficient theorem in Section 6 then excludes the entire independently dressed family, including arbitrary spectator wires, for p at least three. These family arguments have a narrower scope and a much stronger finite-size threshold than the separate [general fixed-clean source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md).

Let k≥1 be odd, n=k+1, and let a be the single clean signal. The remaining k wires R have arbitrary inputs, including entanglement with external references. Put

```math
P=Z_a\prod_{j=1}^k Z_j,\quad
C_A=\prod_{j=1}^k\mathrm{CX}_{j\to a},\quad
C=C_A S_a C_A.
```

The CNOTs in C_A commute, so C_A=C_A†. The literal Clifford identity is

```math
C=e^{i\pi/4}\exp(-i\pi P/4)=\frac{e^{i\pi/4}}{\sqrt2}(I-iP).
```

Thus C uses 2k CNOTs and one positive S; C† uses the same CNOTs and three positive S gates. This is an actual parity circuit, not a supplied reflection.

Define the actual native T band

```math
T_X=HTH,\quad D=\bigotimes_{j\in\{a,1,\ldots,k\}} T_{X_j},\quad
U_A=DC_A,\quad U_B=CDC^\dagger H_a.
```

Each encoder has one physical T band, with n T gates. The source and its actual inverse are

```math
U=U_B^\dagger U_A=H_aCD^\dagger C^\dagger DC_A,\qquad
U^\dagger=C_A D^\dagger C D C^\dagger H_a.
```

Every operation shown is emitted. No encoder is controlled, and no logical input is initialized.

### 4.1. Exact complete transition

Write J|ψ⟩=|0⟩_a|ψ⟩ and define the mathematical Clifford

```math
V=\prod_{j=1}^k\exp(-i\pi X_j/4).
```

The global phase in D=e^{inπ/8}∏exp(−iπX/8) cancels between D† and D in U. Let V_all=exp(−iπX_a/4)⊗V. The exact identity

```math
D^\dagger P D=P V_{\rm all}
```

therefore gives

```math
CD^\dagger C^\dagger D
 =\frac12(I-iP)(I+iPV_{\rm all})
 =\frac12\bigl[(I-iP)+(I+iP)V_{\rm all}\bigr].
```

The contraction uses both parities of the arbitrary logical input, not just its zero column:

```math
\langle+|C_AJ=I/\sqrt2,\qquad
\langle-|C_AJ=\Bigl(\prod_j Z_j\Bigr)/\sqrt2,\qquad PC_AJ=C_AJ.
```

Since ⟨±|exp(−iπX_a/4)=e^{∓iπ/4}⟨±| and conjugating V by ∏Z_j gives V†, the complete normalized block is

```math
L:=\sqrt2 J^\dagger UJ
 =\frac{1-i}{2}I+\frac{e^{-i\pi/4}}2(V-V^\dagger)
 =e^{-i\pi/4}\left[\frac I{\sqrt2}+\frac{V-V^\dagger}{2}\right].
\qquad\text{(O1)}
```

This equality holds on the entire R Hilbert space and after tensoring any external reference.

For odd k, define

```math
F=\sqrt2\sin\left(\frac\pi4\sum_{j=1}^k X_j\right).
```

Every eigenvalue of ∑X_j is odd, so F is a Hermitian involution. Equation (O1) becomes

```math
L=\frac{e^{-i\pi/4}}{\sqrt2}(I-iF),\qquad L^2=-F.
\qquad\text{(O2)}
```

In particular L is unitary and the source is exactly flat at 1/√2 amplitude. The equal-sector lemma then proves that the actual clean-image reflections A=U_A Z_a U_A† and B=U_B Z_a U_B† anticommute.

F is a nontrivial Clifford. To see this explicitly, use X-basis labels z_j∈{0,1}, w=∑z_j. Its eigenvalue is (−1)^{q_k(z)}, where

```math
q_k(z)=\sum_{i\lt j}z_i z_j+a_k\sum_jz_j+c_k\pmod2,
```

and (a_k,c_k) equals (1,0),(0,0),(1,1),(0,1) when k modulo 8 is respectively 1,3,5,7. Thus F is a quadratic Boolean phase in the X basis. Both signs occur for every odd k≥1: the displayed quadratic is nonconstant for k≥3, and is the single input bit for k=1. The transition has eigenvalues −i when q_k=0 and 1 when q_k=1. For k=1, explicitly,

```math
L=\frac{(1-i)I-(1+i)X}{2}=-iHSH.
```

For k≥3 this family includes non-Clifford transitions, but the identity L²=−F remains exact.

### 4.2. Actual overlapping symplectic geometry

Order binary Pauli coordinates as (x,z). The active A axes are X_j on all n physical wires. The active B axes are

```math
C X_j C^\dagger=Y_j\prod_{\ell\ne j}Z_\ell.
```

They are independent commuting axes, all active on the respective original Pauli (P for A, X_a for B). Let E be the n×n all-ones matrix. Since n is even, E²=0 over F₂. The reflection symplectic maps I+N_A and I+N_B satisfy

```math
N_A=\begin{pmatrix}0&I\\0&0\end{pmatrix},\qquad
N_B=\begin{pmatrix}E&I\\0&E\end{pmatrix},\qquad
N_A N_B=N_BN_A=\begin{pmatrix}0&E\\0&0\end{pmatrix}\ne0.
\qquad\text{(O3)}
```

The nonzero product has rank one, while both individual active axis spans have dimension n. The source bands genuinely involve the signal and overlapping physical axes. This establishes realizability of the proposed nonzero-product geometry; it does not establish that this geometry is necessary for an eventual indicator source.

### 4.3. Resources, exact return, and stopping diagnosis

In the positive-S Clifford alphabet, both U and its actual inverse obey

```math
T_U=2(k+1),\qquad D_T(U)\le2,\qquad G_U\le9k+9.
\qquad\text{(O4)}
```

Indeed U_A and its inverse cost 3k+2 Cliffords; U_B and its inverse cost 6k+7 Cliffords. D and D† each cost 2(k+1) Hadamards. U_B contains one C and one C†, which together cost 4k+4 positive-S/CNOT Cliffords. The final H_a adds one. Width is exactly k+1, with one initialized signal and no private work. Other arbitrary dirty spectators can be tensored with identity and return exactly. There is no approximate synthesis error in these identities.

The existing signal-S recovery applies to this flat unitary block as an initialized-isometry identity, irrespective of whether L is an XOR query:

```math
V_S=US_aU^\dagger S_aU,\qquad V_SJ=e^{i\pi/4}JL.
\qquad\text{(O5)}
```

It has T≤6(k+1), D_T≤6, forward G≤27k+29 and actual-inverse G≤27k+33. The signal returns exactly. The displayed fixed phase is retained. This is the return ledger for L. For an arbitrary exactly returned L, the borrowed-sector wrapper still has the exact conjugation identity in Section 6.1 of the clean-removal theorem. It yields a good XOR approximation only when L approximates that complete XOR target; unitarity alone does not supply it.

For any unitary involution Q on the complete R space, including I_p tensor arbitrary returned dirty identity, telescoping the squares gives

```math
2\|L-Q\|\ge\|L^2-Q^2\|=\|-F-I\|=2.
```

Thus the literal source has ||J†UJ−Q/√2||≥1/√2. Allowing any common phase γ of unit modulus still gives

```math
2\|\gamma L-Q\|\ge\|-\gamma^2F-I\|
 =\max\{|\gamma^2-1|,|\gamma^2+1|\}\ge\sqrt2,
```

because both signs occur in F. Consequently

```math
\inf_{|\gamma|=1}\|\gamma J^\dagger UJ-Q/\sqrt2\|\ge\frac12.
\qquad\text{(O6)}
```

This is a source-specific operator-norm floor, including every logical column, dirty input and external reference. It persists when k=O(2^p), although (O4) then fits the desired source count/width/Clifford budget. The source therefore fails before approximate indicator promotion: its complete logical action has the wrong square. It does not fail because alignment, source return, inverse execution, or physical overlap were left unpaid.

The historical reopening question here allowed a changed complete transition or independently audited logical dressing. Section 4.4 answers the dressing question positively for the involution condition, so the square argument alone no longer diagnoses that repaired candidate. Section 6 supplies the stronger stopped-route statement for every independent Clifford dressing of this family. Different complete transitions outside that family theorem and compositions with a separately proved full action remain outside these particular screens. The separate [fixed-clean source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) now also excludes arbitrary other two-band pairs as a scalable arbitrary-accuracy source at fixed amplitude. Additional source bands require a new analysis.

The [exact native fixture](../../tests/test_signal_mixing_overlap.py)
checks k=1 and k=3 on every column, including the complete source inverse,
literal recovery phase, exact signal return and positive-S ledger. The
all-k proof above supplies the growing-width statement. The two displayed
bands overlap and cannot use the old disjoint-band merge; no claim that
two is the source's minimal depth under every alternative synthesis is
needed or made.

### 4.4. An actual Clifford repair of the involution condition

Let $`P_R=\prod_{j=1}^k Z_j`$ act only on the logical wires and put
$`\omega=e^{i\pi/4}`$. In the X basis, P_R complements all k labels.
Because k is odd, the sine defining F changes sign, giving
$`P_RFP_R=-F`$. Thus the changed transition

```math
L'=\omega P_RL=P_R\exp(-i\pi F/4),\qquad
(L')^\dagger=L',\qquad (L')^2=I
```

is a literal Hermitian involution. It defeats the square-based rejection
without adding a T band or any clean history. Involutivity alone does
not establish the specified indicator action.

The repair has actual one-band encoders. With $`C_L=\omega P_R`$, set

```math
U'=(I_a\otimes C_L)U,
\qquad U'_A=U_A,
\qquad U'_B=U_B(I_a\otimes C_L^\dagger).
```

Then $`U'=(U'_B)^\dagger U'_A`$ and its complete block is
$`L'/\sqrt2`$. The logical Clifford is an endpoint of the actual
encoder; it is not a supplied replacement for either encoder. P_R uses
two positive S gates per logical wire. The literal scalar
$`\omega I=(HS)^3`$ uses six Cliffords on any existing wire, and
its actual inverse uses twelve because each S-dagger is three positive
S gates. The actual inverse is
$`(U')^\dagger=U^\dagger(I_a\otimes C_L^\dagger)`$.
Consequently the forward and inverse ledgers are

```math
\begin{aligned}
T_{U'}&=T_{(U')^\dagger}=2(k+1),&
D_T(U'),D_T((U')^\dagger)&\le2,\\
G_{U'}&\le11k+15,&G_{(U')^\dagger}&\le11k+21.
\end{aligned}
```

Width stays k+1, with one clean signal and no private work. The same
three-call recovery, using the actual inverse, gives

```math
V'_S=U'S_a(U')^\dagger S_aU',\qquad
V'_SJ=\omega JL',
```

with exact signal return on every logical and reference input. Its
forward Clifford count is at most 33k+53; its actual inverse costs at
most 33k+63. Both have T-count at most 6(k+1) and T-depth at most six.
The raw recovery phase omega remains literal. Other dirty spectators
return identically. These are a paid repair and return recurrence, not
an indicator compiler: the Section 6 coefficient floor applies to L'.

More generally, arbitrary independent logical Cliffords C_L,C_R give

```math
U'=(I_a\otimes C_L)U(I_a\otimes C_R),\qquad
U'_A=U_A(I_a\otimes C_R),\qquad
U'_B=U_B(I_a\otimes C_L^\dagger).
```

They preserve the actual one-band encoder grammar and produce complete
block $`C_LLC_R/\sqrt2`$. They add their actual forward/inverse
Clifford costs, no T gates and no initialized wires. The coefficient
screen below allows both endpoints to mix all logical and dirty wires,
including any spectator wires added to the family.

## 5. Complete conditional recurrence and the stopped source interface

For an actual exactly flat source, signal-S recovery returns its sole
signal exactly and has the complete occupied unitary L, with its fixed
common phase retained. The
[exact-flat promotion ledger](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#61-exact-flatness-sharpens-the-complete-error-and-return-recurrence)
pays all 24 source/actual-inverse calls, the arbitrary-input sector
cover, bilinear tests, count, width, Clifford cost, error and return.
It does not charge the four-term T_j above as a supplied primitive.

For source encoders on one declared physical register, take common
forward/actual-inverse majorants T_A,G_A and T_B,G_B. Their private
arbitrary work union has size b_U, excluding the specified address,
old m-bit output and signal. The full source ledger is

```math
T_U\le T_A+T_B,\qquad G_U\le G_A+G_B,\qquad d_U\le2.
```

After recovery and promotion, with source error epsilon at most one tenth,

```math
\begin{aligned}
T_{\rm ind}&\le24(T_A+T_B)+48m+8,\\
G_{\rm ind}&\le24(G_A+G_B)+82m+104,\\
D_T(\mathrm{ind})&\le80,\qquad b_{\rm ind}\le m+b_U+2,\\
\epsilon_{\rm ind}&\le2\sqrt2\epsilon.
\end{aligned}
```

For even m the last eight T gates can be omitted. The signal and toggle
return exactly; every address, old-output and private-work return error
is included in the full operator norm. Actual inverses have the same
majorants and error. Requesting epsilon at most
$`\min\{1/10,\tau/(2\sqrt2)\}`$ pays any fixed target tolerance tau.
Uniform linear T_A+T_B and b_U, and quadratic G_A+G_B, would therefore
meet the complete indicator budget if the requisite sources existed.
The [fixed-clean source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md) now
excludes that hypothesis as an arbitrary-accuracy scalable family at
this fixed amplitude. The conditional restored-address R98 transfer
itself remains the fixed-accuracy law

```math
D_T(\mathrm{frame})
 =O_L\!\left(N/b^2+\log_2\log_2(n+4)\right),
```

with the retained simultaneous count, Clifford, one-clean allocation
and real/complex error bounds. The transfer's log-log term is not removed.

The growing family, its involutory repair and all independent logical
Clifford dressings are stopped by Section 6 for p at least three. The
square obstruction in Section 4.3 is no longer the decisive diagnosis
after the repair. The former first unpaid operation was the complete
normalized indicator transition from arbitrary signal-mixed one-band
encoders. The [direct source-gap theorem](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)
now excludes that scalable operation at fixed amplitude and uniformly
bounded clean count, independently of the source's gate counts or
width. No nonlinear-degree bound follows merely from the number of
Clifford summands or a reflection rank: the new proof instead combines
a uniform norm gap with the sparse original-coordinate indicator
probe. Reopening requires a mechanism outside that source theorem,
such as additional source bands with a separately established complete
action and resource/error/return recurrence. A different bounded-clean
mechanism also remains possible; varying the clean count within a fixed
bound does not escape the theorem.
R96 remains the strongest unconditional general upper theorem.

The bounded exact tests in
[`test_overlapping_source_probes.py`](../../tests/test_overlapping_source_probes.py)
check every column of an overlapping parity contraction, the general
probe identity, literal native inverses, and an actual flat source whose
transported probe is not Clifford. They verify finite identities, not an
asymptotic indicator compiler or an unrestricted depth theorem.

## 6. A Clifford-invariant coefficient screen for quarter-turn transitions

This theorem covers the growing family and its repaired transition,
with arbitrary independent logical Clifford endpoints, arbitrary common
phase and arbitrary spectator wires. More generally it covers any exact
complete normalized transition that is a Clifford dressing of a
quarter-turn about a Hermitian Clifford involution. There is no proof
that every two-band source has this logical form.

### 6.1. A uniform coefficient dichotomy

Let F be a Hermitian Clifford with $`F^2=I`$ on Q wires and let
$`V_F=(I-iF)/\sqrt2`$. Use the positive Hermitian tensor-Pauli basis
and normalized coefficients

```math
\widehat E(P)=2^{-Q}\mathrm{Tr}(PE),\qquad
c_F=\max_P|\widehat F(P)|.
```

For any Hermitian Pauli R put $`S=FRF`$, which is a signed Hermitian
Pauli. Direct expansion gives the complete operator identity

```math
V_FRV_F^\dagger=\tfrac12(R+S+iRF-iFR).
\qquad\text{(C1)}
```

If S equals R, this is R. If S equals minus R, it is iRF and its
coefficient ceiling is c_F, since Pauli multiplication permutes Pauli
labels and phases. Otherwise R and S have distinct Pauli labels. Their
coefficients in (C1) have magnitude one half: the cross term has zero
coefficient on both. For R the cancellation is

```math
\mathrm{Tr}\bigl(R(RF-FR)\bigr)
=\mathrm{Tr}(F)-\mathrm{Tr}(RFR)=0.
```

For S, use $`S=FRF`$, $`R^2=F^2=I`$ and cyclicity to obtain the
same cancellation. This step does not require F to be traceless. Every
other coefficient comes from the two Pauli multiples of F and has
magnitude at most c_F by the triangle inequality. Consequently every
image in (C1) is either a signed Hermitian Pauli or has every normalized
Pauli coefficient bounded by

```math
\max\{1/2,c_F\}.
\qquad\text{(C2)}
```

The retained [flat Clifford Pauli-expansion lemma](NATIVE_BAND_SOURCE_STRUCTURE.md#3-native-right-extension)
shows that the nonzero coefficients of a Clifford all have magnitude
$`2^{-r/2}`$, where r is the rank of its symplectic map minus identity.
If F is not proportional to a Pauli, r is positive and
$`c_F\le1/\sqrt2`$. If F is a signed Pauli, V_F is itself Clifford,
including a constant scalar when F is plus or minus identity.

Thus, for every Hermitian Clifford F, images under V_F are either
signed Paulis or have coefficient ceiling $`1/\sqrt2`$. Tensoring
spectators leaves this dichotomy unchanged. Independent left and right
Clifford factors respectively permute the output and input Hermitian
Paulis, retaining their signs; a common scalar cancels in conjugation.
The dichotomy therefore holds for
$`C_L(V_F\otimes I)C_R`$ with arbitrary such factors on the full
logical-plus-dirty register. No spectator basis promise is used.

### 6.2. The sparse indicator probe gives a robust full-operator floor

For the complete indicator with $`m=2^p`$, any original bank observable
has the image

```math
M_j=I_pZ_{Y_j}I_p^\dagger
   =Z_{Y_j}(-1)^{[x=j]},\qquad
\widehat M_j(Z_{Y_j})=1-2/m.
\qquad\text{(C3)}
```

For p at least two, M_j is not a signed Pauli. Its operator distance
from every signed Hermitian Pauli is at least square root two, as in
the [complete-indicator Pauli-image gap](INDICATOR_PAULI_NULLITY.md#3-a-uniform-pauli-image-gap-makes-the-count-robust).
Explicitly, a nonzero X mask sends each basis input to a different
destination from the diagonal M_j, giving column distance square root
two. A zero X mask is an affine diagonal sign; the nonaffine sign in
(C3) differs on some input, giving distance two. This argument includes
all arbitrary dirty directions.

Suppose a candidate unitary W has the dichotomy (C2) with ceiling c,
where $`1-2/m\gt c`$. If its Z_j image is Pauli, the preceding
gap applies. Otherwise its coefficient at Z_j has magnitude at most c,
and normalized trace is bounded by operator norm, so

```math
\|M_j-WZ_jW^\dagger\|\ge1-2/m-c.
```

Unitary telescoping gives
$`\|M_j-WZ_jW^\dagger\|\le2\|W-I_p\|`$. A common phase
on W cancels, and tensoring references preserves these bounds. Therefore
an exactly flat source with complete normalized transition W obeys

```math
\inf_{|\gamma|=1}
\left\|\gamma\frac W{\sqrt2}-\frac{I_p}{\sqrt2}\right\|
\ge\frac{1-2/m-c}{2\sqrt2}.
\qquad\text{(C4)}
```

In particular, every independently Clifford-dressed quarter-turn
transition in Section 6.1 satisfies, for p at least three,

```math
\epsilon\ge
\frac{1-2/m-1/\sqrt2}{2\sqrt2}
\ge\frac{3}{8\sqrt2}-\frac14\gt0.
\qquad\text{(C5)}
```

This is an operator-norm obstruction to meeting every fixed requested
tolerance. It assumes neither monomial action of the candidate nor an
initialized output or dirty register. The native actual inverse has the
same distance from the inverse target. The proof constrains an exact
logical form, not the T-depth of arbitrary sources.

### 6.3. Sharper floor for the growing family and its repair

For the Section 4 family with odd k at least three, the quadratic
polynomial q_k has polar matrix equal to the all-ones matrix plus
identity. Its radical is $`\{0,r\}`$, where r is the all-ones
vector, and $`q_k(z+r)=q_k(z)+1`$. With normalized Boolean Walsh
coefficients, squaring and summing over the displacement gives

```math
\left|\widehat{(-1)^{q_k}}(t)\right|^2
 =2^{-k}\bigl[1-(-1)^{t\cdot r}\bigr].
```

Hence every nonzero Pauli coefficient of F has magnitude
$`2^{-(k-1)/2}\le1/2`$. Equation (C2) improves the non-Pauli
image ceiling to one half. For arbitrary independent Clifford dressing,
arbitrary spectators and any common phase, (C4) now yields

```math
\epsilon\ge\frac{1/2-2/m}{2\sqrt2}
\ge\frac1{8\sqrt2},\qquad p\ge3.
\qquad\text{(C6)}
```

For k equal to one, the literal L is already the Clifford
$`-iHSH`$. Every dressed image is Pauli, so the same target probe
and its square-root-two gap give source error at least one half.
This exception must not be assigned the non-Clifford family's Pauli
subspace dimension.

Equations (C5) and (C6) also include the actual involutory repair of
Section 4.4. The repair fixes a real defect, but its complete probe
profile still misses the sparse indicator response. This is the family-specific
reason it stops, with explicit useful finite-size constants. Reopening requires a complete transition
outside this independently dressed quarter-turn class, or a composition
whose full operator, resources and return are separately established.
This coefficient theorem does not exclude arbitrary other two-band
transitions. The separate [fixed-clean source gap](FIXED_CLIFFORD_SUM_SOURCE_GAP.md)
does exclude them as a scalable arbitrary-accuracy source at fixed
amplitude, while additional source bands and a faster unrestricted
indicator remain outside that statement. The
complete-frame frontier therefore remains unchanged.

The [bounded exact receipts](../../tests/test_dressed_signal_overlap.py)
check the repaired source and actual inverse, literal recovery phase,
complete signal return, positive-S counts, and the coefficient identity
for a nontraceless Hermitian Clifford. They support the displayed native
identities; the all-width, all-dressing floors follow from the proofs
above.
