# Programmed rotations give logarithmic depth in the precision parameter

[All-precision upper theorem R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md) ·
[Portable complete catalogue R94](PORTABLE_WORKSPACE_CATALOGUE.md)

A short dyadic program specifies an exact unitary close to a one-qubit
rotation on arbitrary dirty inputs. Parallel copies of a clean implementation
are accessed on their own zero-key sectors. One displacement broadcast loads
the program into every copy without multiplying the large query cost.
Two dirty seed bases and a Pauli conjugation then give the rotation on the
entire input space. The rotation compiler uses zero initialized work.

**Theorem R99.** There is an absolute $`c_*>0`$ with the following
property. Let $`N=2^n`$ and take integers $`n\ge1`$ and
$`6\le L\le c_* n`$. At $`b\ge17(L+n+7)`$ for a real complete
Hopf frame, or $`b\ge17(L+n+8)`$ for its prescribed phase-dressed complex
magnitude completion, there is one native circuit V with one external
clean input and at most b arbitrary dirty qubits, such that

```math
\begin{aligned}
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}&\le 2^{-L},\\
T&=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\\
D_T&=O\bigl(NL/b^2+\log_2(n+2)+\log_2^2(L+2)\bigr).
\end{aligned}
```

The constants are absolute and the same bounds apply to the actual
native inverse, with $`W^\dagger`$ in the comparison. The error includes
every logical column, literal phases, arbitrary dirty and reference
inputs, and all work return.
Classical preprocessing terminates without a running-time claim. At
$`b=\Theta(\sqrt{NL})`$ and $`L=(\log n)^A`$ for fixed $`A>1`$, the bound
is $`O(\log n)`$. At $`L=\Theta(n)`$ within the stated range and the same allocation,
it is $`O(\log^2 n)`$. The fixed-accuracy logarithmic upper remains.

<a id="an-exact-invariant-clean-source-with-a-short-program"></a>

## 1. An exact invariant clean source with a short program

Write $`R_P(\theta)=\exp(i\theta P)`$ for a fixed one-qubit Pauli axis.
The real frame's rotation convention uses $`P=Y`$ and the opposite angle
sign. On that target and two arbitrary payload qubits choose five
pairwise anticommuting anti-Hermitian Clifford matrices $`\Gamma_j`$
with $`\Gamma_j^2=-I`$ and $`\Gamma_1=iP\otimes I`$. For $`P=Y`$ on
wire 0 take

```math
iY_0,\quad iX_0X_1,\quad iX_0Z_1,\quad iX_0Y_1X_2,\quad iX_0Y_1Z_2.
```

A fixed target Clifford gives other axes. For requested parameter s and
an absolute constant C put $`K=2s+C`$. Choose signed integers a and b
within $`O(1)`$ of $`2^K\cos\theta`$ and $`2^K\sin\theta`$, with
magnitudes bounded by the corresponding true magnitudes. Certified
rational estimates and rounding a certified lower bound on the magnitude
down give this guarantee even at an integer tie.
[Lagrange's four-square theorem](../../docs/reference/SOURCE_CATALOGUE.md#5-fault-tolerant-sources-and-contribution-boundaries)
gives integers $`c_2,\ldots,c_5`$ with
$`\sum_jc_j^2=2^{2K}-a^2-b^2`$, a nonnegative residual of order at most
$`2^K`$. A finite enumeration finds them. Define

```math
U_\theta=2^{-K}\left(aI+b\Gamma_1+\sum_{j=2}^5c_j\Gamma_j\right).
```

Anticommutation gives $`U_\theta^\dagger U_\theta=I`$ exactly.
Using $`\sum_j|c_j|\le2(\sum_jc_j^2)^{1/2}`$, its full-operator error
from $`R_P(\theta)\otimes I`$ is $`O(2^{-K/2})`$ and its coefficient
one-norm is at most $`\sqrt2+O(2^{-K/2})\lt2`$ after fixing C.
Use a uniform index on $`K+1`$ bits, with $`|a|,|b|,|c_j|`$ occurrences
of the respective signed matrix, and fill the remaining occurrences
with equal numbers of $`I`$ and $`-I`$. The deficit is even because the
sum of magnitudes has the same parity as the sum of squares.
The uniform average is exactly $`U_\theta/2`$. A constant number of
cumulative cutoffs of width $`K+2`$ stores this program, including the
upper endpoint. These cutoffs carry all the angle data.

Append a clean selector and fixed comparator flag bank A of size
$`K+O(1)`$. Hadamards on the selector, signed Clifford selection, and
the same Hadamards give a word S with $`J^\dagger SJ=U_\theta/2`$, where J
inserts A zero. Implement each exact cutoff comparison by an
[R91 word query](PARALLEL_BILINEAR_WORD_QUERY.md) on $`O(s+1)`$ selector
and program bits, toggling a private flag. It uses $`2^{O(s+1)}`$ count
and dirty work at an adequate pool, and $`O(\log(s+2))`$ T-depth.
The flag controls the signed matrix and is uncomputed from the unchanged
comparison inputs. A controlled Hermitian Pauli uses Clifford gates.
The literal factors $`-1`$ and $`\pm i`$ are Clifford phases on the flag.
For malformed cutoff strings the same ordered circuit still defines a
unitary and preserves the program bits.

Put $`\Pi=JJ^\dagger`$ and $`R=I-2\Pi`$. Pay R by an exact
[returned dirty zero test](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
and its phase echo on $`O(s+1)`$ controls. With $`A_S=J^\dagger S J`$,
multiplication gives

```math
J^\dagger SRS^\dagger RSJ=4A_S A_S^\dagger A_S-3A_S=-U_\theta.
```

The right hand side is unitary, so the clean sector returns exactly for
every payload and private dirty input. Denote this word by $`V(P_{\rm cut})`$
for its read-only cutoff parameter. Its scalar minus appears in the
accepted action and cancels against the actual adjoint below. Each word
and actual reverse have standalone $`2^{O(s+1)}`$ count and width and
$`O(\log(s+2))`$ T-depth.

## 2. Load one program and form masked reflections

An immutable a-bit address x specifies an angle $`\theta_x`$. Put
$`Q=2^a`$ and let $`f(x)`$ be the $`q=O(s+1)`$-bit cutoff program for the
angle $`-\theta_x/4`$. Use h slots. Slot j has an arbitrary key $`A_j`$
of $`k=O(s+1)`$ bits containing the cutoff, selector and comparator fields,
three arbitrary payload qubits, and disjoint local dirty work.

Use an exact R91 query $`Q_f`$ into a central arbitrary word z with an
exclusive pool $`B\ge C_Q(a+2)`$. Let F add z into all cutoff fields by
CNOTs. In chronological order $`F,Q_f,F,Q_f^\dagger`$ adds $`f(x)`$ into
every cutoff field. The old z contributions cancel, returning z and
all query work. F costs $`O(hq)`$ Cliffords. Write $`L_f`$ for this
broadcast. Perform $`L_f`$, the disjoint local V words, and actual
$`L_f^\dagger`$. The result is the product of individual slot words at
each x. On any slot's original $`A_j=0`$ sector its program loads a valid
cutoff, its selector returns exactly, and its cutoff returns by the
inverse broadcast. Unitarity preserves the complementary subspace. The
sequential exact query boundaries justify the product also for coherent
x, including arbitrary reference entanglement.

Let $`B_0`$ be a target Pauli anticommuting with P and
$`\Pi^-_{B_0}=(I-B_0)/2`$. Between the global slot word and its actual
inverse insert the product of the exact slot reflections
$`I-2[A_j=0]\otimes\Pi^-_{B_0}`$. A returned conjunction phase on the
original key and the payload target pays each reflection after the
fixed Clifford change for $`B_0`$. Distinct slots use disjoint private
supports. The completed product $`T_{\rm all}=\prod_jT_j`$ satisfies

```math
T_j\big|_{A_j=0}=J'_j=U_{-\theta_x/4}^\dagger B_0U_{-\theta_x/4},
\qquad T_j\big|_{A_j\ne0}=I.
```

The equality off the zero key is literal even when the individual
slot word mixes occupied keys. Thus the completed $`T_j`$ preserves
every computational key label and returns all private work outside its
payload. There are constantly many global R91 calls for the whole
reflection product and its actual inverse, independent of h.

## 3. Select one qualifying slot with a shallow swap

Let S exchange the external three-qubit payload with that of the first
zero-key slot, acting identically when none qualifies. Consider the
word in operator order

```math
E=T_{\rm all} S T_{\rm all} S^\dagger,
\qquad U_*=E B_0 E^\dagger B_0.
```

On a key sector where j qualifies first, all other $`T_i`$ commute
through S and cancel. The remaining involutions act on disjoint
payloads, giving $`E=J'_j\otimes J'_{\rm ext}`$ and
$`U_*=(J'_{\rm ext}B_0)^2`$. The ideal J is the B_0-conjugation by
$`R_P(-\theta_x/4)`$, and $`(J B_0)^2=R_P(\theta_x)\otimes I`$ with
literal phase. The error from J' is at most four times the slot
rotation error. The selected slot, original keys and every other
private input return exactly in this comparison. When there is no
qualifying key, S is identity and the reflections cancel.
This is a maximum over key sectors, so many zero keys do not multiply
the error.

Write $`a_j=[A_j=0]`$ and $`b_j=a_j\prod_{i\lt j}(1-a_i)`$. Three XOR
sweeps pay the ideal S. The first adds each external payload bit into
the corresponding slot coordinate under $`b_j`$. The second adds
back the XOR of the selected slot bits. The third repeats the first.
The predicates are mutually exclusive and do not read the payload.
In the output Fourier basis of a sweep, its entire exponent is the
parity of terms consisting of $`b_j`$, the source bit and the output
character. A Pauli mirror on the external target uses the success
predicate $`\bigvee_j a_j`$ times one target character, after changing
its axis to Z. Every sweep and mirror is a complete operator on
arbitrary inputs.

### Fingerprints of the priority predicates

For each prefix choose $`r_h`$ independent uniform binary linear forms
on its atoms. Replace the no-hit prefix predicate by

```math
\prod_{t=1}^{r_h}\bigl(1+\ell_{j,t}(a_1,\ldots,a_{j-1})\bigr)
\quad\text{over }\mathbb F_2.
```

This is exact on the all-zero prefix and has false-acceptance
probability $`2^{-r_h}`$ on a nonzero prefix. Use the same construction
on all atoms with one final complement for the success predicate.
Predict the complete signs of the three XOR sweeps by substituting
these prefix polynomials. For any full computational label of keys
and payload/character bits, their failure probability is at most
$`(h+1)2^{-r_h}`$ by the union bound.

Choose an operator tolerance $`\delta=2^{-\Theta(s+1)}`$ and let m_0
bound the number of Boolean inputs in any one sign. Here
$`h,m_0\le2^{O(s+1)}`$. Put $`\delta_0=c\delta^2/(m_0+1)`$ with a small
absolute c, and choose $`r_h=O(s+1)`$ so each failure probability is at
most $`\delta_0`$. Sample a power-of-two list of
$`M=\Theta((m_0+1)^2/\delta^2)`$ independent predictions. With a fixed
sufficiently large prefactor the Chernoff bound and a union over all
$`2^{m_0}`$ labels give a list with at most $`2\delta_0M`$ bad predictions
per label. Choose the constants so $`2\delta_0\sqrt M`$ meets the error
allowance of the
[three-basis correction and bounded-clean promotion](POSITIVE_HAAR_FINGERPRINT.md#1-a-fingerprint-seed-can-start-in-an-arbitrary-quantum-state).
For example choosing the prefactor proportional to $`1/c`$ fixes the
Chernoff exponent while making $`2\delta_0\sqrt M/\delta=O(\sqrt c)`$.
A finite search checks all labels. Separate lists for the sweeps and
mirror suffice. The cited correction records and phases the majority
in three bases on the same arbitrary dirty prediction seed, actually
reverses the recordings, and promotes its bounded clean allocation.
It gives the chosen full operator error on arbitrary controls,
output banks and work.

### Exact phases for products of hash tests

For a fixed sampled prediction and one priority term, form a vector
F of required zero values. Its coordinates are the r_h hash parities,
the complemented hit bit $`1-a_j`$, the complemented payload and
output character literals, the mismatch bits of the selected
prediction-seed label, and the complemented vote-control bit for
the three-basis recording. Thus F has $`r=O(s+1)`$ coordinates. Each is an
affine function of stable bits or a parity of at most h of the
zero-key atoms, with a fixed offset. The term's sign is
$`(-1)^{[F=0]}`$. The mirror uses one such no-hit term in addition to
its target Pauli, since the predicted success bit is the complement
of the no-hit bit and exactly one seed label applies. During a vote
recording this Pauli uses the same vote control, paid by CZ in the
target Fourier basis.

Here is a native implementation that leaves the products of hash
tests unevaluated until their final equality. Take an arbitrary
r-bit register z. The XOR query $`Q_F:z\mapsto z\oplus F`$ is exact.
In its output Fourier basis each occurrence of a zero-key atom is one
paid conjunction phase on its k key bits and a Fourier character.
Use the [returned, read-only conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
and one private phase helper for each term. Affine offsets and
literal mismatch bits use Cliffords. Distinct terms have private T
targets and share only unchanged controls and character bits. Thus
$`Q_F`$ and its actual inverse have depth $`O(\log(s+2))`$ and
$`h\,\mathrm{poly}(s+1)`$ gates and dirty work. No physical output of
an atom is initialized or stored as a promised bit.

In two separate arbitrary $`2^r`$-bit banks U,V use the same exact
small dirty indicator. Let A be the query that adds $`e_z`$ into U.
Let B add $`e_{z\oplus F}`$ into V using the actual conjugation by
$`Q_F`$ on the indicator address. Both queries return z and their
private helpers. For the fixed Clifford sign
$`C|U,V\rangle=(-1)^{U\cdot V}|U,V\rangle`$, the chronological word

```text
C, A, C, B, C, A^dagger, C, B^dagger
```

adds only the cross term $`e_z\cdot e_{z\oplus F}=[F=0]`$ to the
phase. It returns every bank for every starting value, including
arbitrary coherent z. Each indicator has depth $`O(\log(r+2))`$
and $`2^{O(r)}`$ gates and width. C is the $`2^r`$-term CZ word.
The complete term and its actual reverse therefore cost
$`2^{O(s+1)}`$ gates and width and $`O(\log(s+2))`$ T-depth.

There are $`O(h)`$ terms for a sampled sweep exponent and one
no-hit term per sample for the mirror. The number of sample labels,
all private banks, and the bounded-clean correction and promotion
add at most $`2^{O(s+1)}`$ storage and native gates. All exact term
words can run in parallel: z, both indicator banks and all nonlinear
supports are private, and the raw key and payload/character bits are
read-only Clifford controls. Selection of the sample label uses the
mismatch coordinates of F, with no separate serial pass through the
samples. No external angle address bits enter these predicates.

Hadamards on the outputs give the complete sweep approximations.
Replacing S by the three actual approximate sweeps incurs at most
four times their combined error in $`U_*`$ by a constant-length
unitary hybrid. Use the actual inverse swap word in each reverse
appearance. This hybrid includes key and payload motion in the
compiled approximation.

## 4. Cover both failure subspaces without initialized work

Take $`h=\lceil 2^k M_h\rceil`$ with $`M_h`$ a sufficiently large
absolute multiple of $`s+k+1`$, and write $`w_A=hk`$ for all key bits.
For each of two lists and each of D independently sampled labels,
choose a uniform XOR shift of the $`w_A`$-bit key string. At any fixed
original key string, a label has no zero-key slot after shifting with
probability $`\rho=(1-2^{-k})^h`$, independently across its list.
Choose D the least power of two at least $`C(w_A+1)/\rho`$ for a fixed
large C. The binomial upper tail and a union over the two lists and
$`2^{w_A}`$ key strings give lists with at most $`2\rho D`$ bad labels
at every original key string. Finite enumeration with exact counting
finds them. The bounds $`\exp(-O(M_h))\le\rho\le\exp(-M_h)`$ give
$`\log_2 D=O(s+1)`$. Choose the multiplier in $`M_h`$ before D so that
$`2\rho\sqrt D\le 2^{-s}/C_0`$ for a fixed required error constant
$`C_0`$.

Use the computational seed basis for list one and its tensor-Hadamard
basis for list two, on the same arbitrary dirty register. At each
fixed original key string the bad projectors $`P_1,P_2`$ satisfy
$`\|P_1P_2\|\le 2\rho\sqrt D`$. The overlap matrix between the
restricted bases has at most $`(2\rho D)^2`$ entries of magnitude
$`D^{-1/2}`$, so its Frobenius norm gives the bound. Implement each
shift lookup by an exact complete indicator on $`\log_2 D`$ seed bits,
a Clifford matrix shear from its arbitrary indicator bank to the
keys, the actual indicator reverse, and a second shear. The old
indicator bank returns and the seed is preserved. The shifts and
actual inverses cost $`2^{O(s+1)}`$ gates and width and
$`O(\log(s+2))`$ native T-depth.

For angle $`\alpha_x=\theta_x`$ or $`\theta_x/2`$, the ideal shifted
first-slot word in basis i is

```math
U_i(\alpha_x)=P_i\otimes I+(I-P_i)\otimes R_P(\alpha_x).
```

The same slot dimensions, seed, and shifts serve both angles. Their
cutoffs are separate table programs of the same width. The native
word uses the shift, the first-slot rotation and the actual reverse
shift, with the indicated Clifford basis changes. The first-slot
estimate is a uniform direct-sum bound on the original key, seed and
x sectors, with the complete swap error included. Write
$`C_1=P_1\otimes I+(I-P_1)\otimes B_0`$ for the ideal mirror on
success in basis one. The same phase prediction compiler, with the
appropriate mask shift and reverse, implements it to the requested
full-operator error.

At fixed original key labels write $`\mu=\|P_1P_2\|\lt1`$ and replace
$`P_2`$ for analysis by the projector onto
$`(I-P_1)\mathsf{ran} P_2`$. The rank is unchanged. Every unit
vector in the original range has distance at most $`\mu`$ to the
new one, so the projector norm difference is at most $`\mu`$ by
principal angles. The new bad subspaces are orthogonal. Emit the
word in operator order

```math
U_1(\theta_x)\,U_2(\theta_x/2)\,C_1\,U_2(\theta_x/2)\,C_1^\dagger
```

with their actual native approximants and the actual mirror reverse.
For orthogonal bad subspaces the good projectors commute and at least
one succeeds. On their joint eigenspaces let $`g_i`$ be the good bit.
Conjugation by $`B_0`$ reverses the angle, giving total literal angle

```math
\theta_x\left[g_1+g_2/2+(1-2g_1)g_2/2\right]
=\theta_x(g_1\lor g_2)=\theta_x.
```

This gives the full target unitary on arbitrary target and work,
including all coherent seed and key sectors. The analytical projector
substitution costs at most $`4\mu`$ over the two U_2 terms. The remaining
replacements use a constant-length full-operator hybrid with actual
inverses. Choose the slot, swap, mirror and overlap errors as small
fixed multiples of $`2^{-s}`$. This gives a zero-clean, full-operator
rotation error at most $`2^{-s}`$ with all work included. The two seed
bases require no supplied quantum state, and the success conditions
do not use initialized flags.

### Native ledger and preserved guards

Both cutoff tables have $`O(s+1)`$ bits per row and use a fixed number
of full R91 calls, including all actual inverses. All selector and
comparison programs, complete payload banks, Fourier signs, short
seed lookups and work are charged. For the exclusive query pool
$`B\ge C_Q(a+2)`$, the word and actual inverse obey

```math
\begin{aligned}
T&=O\bigl(\sqrt{Q(s+1)}+Q(s+1)/B\bigr)+2^{O(s+1)},\\
G&=O(Q(s+1))+2^{O(s+1)},\\
D_T&=O\bigl(\log_2(a+2)+Q(s+1)/B^2+\log_2(s+2)\bigr),\\
w&\le B+2^{O(s+1)},\qquad a_{\rm clean}=0.
\end{aligned}
```

Width w excludes the external address and rotation target, and
includes the two external dirty payload helpers. A constant-angle
word uses fixed Clifford broadcasts with no global query pool.
For a stable Boolean enable e, approximate the half-angle and use

```math
\widehat R_P(\theta/2)\,B_0^e\,
\widehat R_P(\theta/2)^\dagger\,B_0^e
```

in operator order. The word is literal identity on e zero, including
arbitrary occupied scratch. On e one its error is at most twice
the half-angle error. The control wires remain unchanged by the
rotation word. A single flag uses an exact controlled Pauli, and a
short stable conjunction uses the returned dirty phase echo with an
appropriate target Clifford. Reserve a fixed extra precision constant
for the guard echo. With shared long addresses use the physical
read-only R91 variant. All local circuits that do not read the long
address have disjoint private T targets, so corresponding native bands
may run in parallel.

## 5. Replace the serial precision words in the opening

Write $`\ell=\log_2(n+2)`$. In the wide real construction assume
$`b\ge\sqrt N/n`$ and $`L\ge\ell`$. Allow the internal precision range
$`L\le2c_* n`$ for the complex shift below. Use three contractions from
[R96](ALL_PRECISION_LINEAR_DEPTH.md#4-moderate-precision-three-early-stages-and-a-short-terminal),
with the same global borrowed-precision template and finite positive
magnitude/sign reduction. At height m and precision $`l=L+O(1)`$, put

```math
K=m+l+2,\quad t=\lceil\log_2(m+2)\rceil,\quad
\chi=\lceil\log_2(l+2)\rceil+t^3,\quad r=\lceil A\chi\rceil.
```

Take A after the fixed table and resource constants. Execute an early
stage only when $`r\lt m/2`$, and then retain its original bottom r bits.
Keep the coarse frame, alias probabilities, causal matching, original
conditional work, and saved activity flag. The sign reduction uses
the retained exact full-address sign phases at its original
$`2^{-L}/512`$ error allocation. Three changes shorten the physical
precision words.

**Calibration.** In the
[ordinary source](SMALL_WORKSPACE_DESCENDANT_SOURCE.md#2-replace-calibration-by-a-borrowed-real-rotation),
replace the borrowed real rotation on the attenuation bit by the
programmed rotation of Section 4. Its table reads the same preserved
source label, including the global prefix when present. The cutoff
word has $`O(l+1)`$ bits per row and uses constantly many exact table
calls. Program registers and local work are excluded from the address
and the initialized anchor. Use the half-angle guard for literal
inactive identity and the retained operator error allowance. If
multiple calibration words share an address, use R91 with read-only
Clifford address interactions and disjoint physical query/private
banks, so corresponding native T bands run together. This keeps the
chronological number of query positions and the polynomial bound on
simultaneous instances from the
[actual source address audit](GLOBAL_BORROWED_PRECISION_STAGE.md#5-global-table-addresses-and-physical-scheduling).
Every actual adjoint is used in source reflections and amplification.

**Coefficient bank.** Replace the fixed half-angle rotations in the
same source by disjoint copies of the constant-angle construction.
The existing controlled-Pauli half-angle echoes give literal inactive
identity on occupied source wires. For q sites assign a fixed fraction
of the bank error divided by q to each half-angle. Here
$`q=O(\log(l+2))`$. The required precision is
$`l+O(\log(q+2))`$, giving bank depth $`O(\log(l+2))`$ with disjoint
private T targets and no extra initialized work. The original finite
coefficient sums and attenuation scalars are unchanged.

**Padded geometric references.** In each bin keep the original
[truncation comparison](FAST_CAPPED_GEOMETRIC_REFERENCE.md#1-two-exact-capped-columns-and-their-distance)
from the full cap $`g_{M_i}`$ to the padded cap $`g_R`$, where
$`R\le M_i`$ is a power of two and $`R=O(l+\log(p+1))`$ for p bins.
Prepare the latter by reading its binary positions from high to low.
When v positions remain, write $`z=2^{-2^{v-1}}`$. If the higher bits
are all one, the conditional probability of a one on the next bit
is z, since the capped upper-half weights sum to z times the total
weight. If a higher bit is zero, the remaining interval excludes
the cap. The normalized weights then factor by binary digits, so the
probability is $`z/(1+z)`$. Apply the base real angle for the latter
probability and the angle difference conditioned on all higher bits
one. This prepares the exact positive $`g_R`$ in the ideal comparison.
The guard is a stable conjunction of $`O(\log(R+2))`$ bits. Pay it
with a half-angle word and its actual reverse around the returned
center. All angles are fixed computably evaluable real numbers.
For $`R=1`$ there is no binary position to prepare.

Use Section 4 at precision $`l+O(\log(p+1)+\log\log(R+2))`$ to sum
operator errors across bins and their angle positions. Reference depth
is $`O(\log^2(l+\log(p+1)+2))`$ and different bins run in parallel on
disjoint private banks. The actual preparation is close in full
operator norm to the ideal preparer of the padded cap tensored with
identity on private work. Follow that comparison by the
[analytical ray transport](FAST_CAPPED_GEOMETRIC_REFERENCE.md#2-a-full-space-comparison-for-transported-zero-tests)
from the truncated column to the original column. Choose fixed slack
in the cutoff and synthesis error to fit the original reference
budget. The augmented zero tests still include the original coin
and its initialized scratch, which can be left untouched if unused.
All new private work is dirty and excluded from those tests. The
reference words and actual inverses cancel around the inactive
middle operation for every occupied input.

Give all three replacements the retained source error allowances,
with fixed additive padding of l as needed. The
[complete-stage hybrid and recovery](GLOBAL_BORROWED_PRECISION_STAGE.md#7-complete-error-and-return)
then retain $`2^{-L}/512`$ per stage. Coefficient and calibration words
are full-operator approximations, and the padded reference comparison
has the prescribed anchor on every ideal column. The actual H-zero
identity holds on occupied conditional work, and the saved-predicate
inverse remains last. The initialized reservation $`O(\chi)`$ is
unchanged.

The main table address envelope remains $`n-r+C\chi`$. Enlarge fixed
polynomial factors in K for the cutoff widths, constant call
multipliers and simultaneous instances. The standalone gate and
private work majorant becomes

```math
S'_{\rm loc}=K^{O(1)}\bigl(2^{O(\chi)}+2^{O(l+\ell+1)}\bigr).
```

The conservative l plus ell precision here covers all summed local
errors. Fix the coefficient of L in this exponent before choosing
$`c_*`$. With sufficiently small $`c_*`$ and sufficiently large n,
this term and its width admit a $`2^{n/4}`$ majorant throughout
$`L\le2c_*n`$. The original fixed-polynomial minimum query pools also
fit a fixed fraction of b. As in the
[global ledger](GLOBAL_BORROWED_PRECISION_STAGE.md#6-native-count-width-and-depth-ledger),
choose A so $`r-C\chi\ge(2d+80)\chi`$ for the updated fixed polynomial
exponent d. This absorbs every large table overhead with the original
exclusive private pools. With the separate full-address sign phases
and at most three stages the opening obeys

```math
T_{\rm open}=O(\sqrt N+N/b),\qquad G_{\rm open}=O(N),\qquad
D_{\rm open}=O(N/b^2+\ell+\log^2(L+2)).
```

The remaining height satisfies
$`m_{j+1}\le C(\log(L+2)+\log^3(m_j+2))`$ after a running stage and
also after an omitted stage by its failed reserve test. As
$`\ell\le L\le 2c_*n`$, the first bound is $`O(\ell^3)`$, the second is
$`O(\log(L+2)+\log^3(\ell+2))`$, and the third is
$`M=O(\log(L+2))`$. The last step uses
$`\log(\ell+2)=O(\log(L+2))`$ and
$`\log^3(\log(L+2)+2)=O(\log(L+2))`$. The original bottom coordinates
and every actual global prefix remain in the tail problem.

## 6. Two tail schedules with the same query counts

### Grouped tail for precision below the square-root scale

When $`L\lt\sqrt n`$, keep the
[R94 portable catalogue](PORTABLE_WORKSPACE_CATALOGUE.md#3-exact-transfer-of-the-complete-workspace)
and the
[low-precision grouping from R96](POLYNOMIAL_PRECISION_DEPTH.md#5-a-constant-number-of-final-accuracy-high-precision-groups).
Split $`M-3`$ into constantly many groups of heights
$`h\le\max\{1,\lfloor\log_2(L+2)/4\rfloor\}`$ if $`M\gt3`$.
The unchanged angle code has precision $`\kappa=L+O(1)`$ and length
$`q_h\le(\kappa+8)2^h=O(L^{5/4})`$. Each code specifies the same
controlled local frame $`A_v`$ on h plus one logical bits. The zero
code is still the literally empty native word. For each other code,
compile $`A_v`$ layerwise with the programmed rotation of Section 4,
including the actual local suffix conditions in its angle tables.
Each local layer has at most $`2^{h+1}`$ rows and $`O(L)`$ cutoff bits
per row. The separate sufficient local query pool gives
$`O(\log(L+2))`$ depth per layer, at precision
$`\kappa+O(\log(h+2))`$ so their full-operator errors sum to the old
slot allowance. The retained initialized slot input d can be left
unused. The actual complete local word uses $`2^{O(L+1)}`$ gates and
work and depth $`O(\log^2(L+2))`$. Pad it with identity to the common
slot width w.

For $`V_h=2^{q_h}`$, reserve the whole payload of each slot, the two
complete arrays of width $`2wV_h`$, and $`O(wV_h^2)`$ router work,
including its exact queries and CZ/CNOT sweeps. The local slot words
run on disjoint supports. The router's depth order is
$`O(\log(q_h+2))`$ beyond its paid global query contribution,
independent of payload width. This follows from the
[complete payload ledger](PORTABLE_WORKSPACE_CATALOGUE.md#4-native-resources-and-the-final-accuracy-code),
which charges $`wV_h^2`$ in counts and width and uses disjoint private
targets for the matched swaps. All expanded catalogue and router
costs are $`2^{O(L^{5/4}+L+1)}=2^{o(n)}`$ in this branch. The new
internal cutoffs and masks stay inside the local slot, so they do not
increase the outer angle-code width or the global table address.

With s_i original lower bits under group i the main table has
$`Q_i\le4\,2^{n-s_i-h_i}`$ rows. Hence $`Q_iq_{h_i}=O(NL2^{-s_i})`$,
and its square roots and volumes have the original geometric sums
$`O(\sqrt{NL})`$ and $`O(NL)`$. Reserve the exclusive $`\Theta(b)`$
query pool disjoint from programs, flags, live addresses and all
arrays. The saved predicate, original guard, changing-prefix
schedule, and canonical code error are those of R94. The exact
slot transfer cancels the same physical slot words on occupied
inactive inputs. In the selected initialized sector their retained
error bound follows from the new complete local-layer hybrid.
All local calls and inverses are the actual emitted words.

The at most three bottom layers use Section 4 directly. Program
only their preserved upper prefix and use the lower-suffix zero
conjunction as the separate half-angle guard. Their cutoff volumes
are $`O(NL)`$ in total with a constant number of sequential calls.
They reproduce the ideal complete addressed layer, including its
inactive action, and can meet the old bottom error with fixed
additional precision. The combined grouped tail has T-count
$`O(\sqrt{NL}+NL/b)`$, Clifford count $`O(NL)`$, and depth
$`O(NL/b^2+\ell+\log^2(L+2))`$.

### Individual tail layers at larger precision

For $`L\ge\sqrt n`$ the residual $`M=O(\log(L+2))`$ layers are few
enough to program individually. Use the exact order of the
[complete addressed-depth factors](../../docs/HOPF_INTERFACE.md#4-complete-addressed-depth-operator).
For a layer with j original lower bits, its preserved upper address
has $`n-j-1`$ bits and thus $`Q_j=N2^{-(j+1)}`$ rows. Program the
half-angle on that address. Keep the lower bits outside both cutoff
tables, and use their unchanged all-zero condition in the exact
Pauli-center guard of Section 4. A returned conjunction on
$`j+1=O(\log(L+2))`$ literals pays it in $`O(\log(j+2))`$ native
T-depth. No extra initialized predicate is used. The actual rotation
and its actual inverse cancel exactly in the off sector. Each
completed layer approximates the whole addressed operator on every
input, including all suffix and reference states.

Request $`s=L+\lceil\log_2(M+1)\rceil+C`$ with fixed C so the guard
multiplier and M layer errors fit the old tail allowance by unitary
telescoping. The cutoff width is $`O(L)`$. The row volumes retain
their full global prefixes and satisfy

```math
\sum_j Q_jL=O(NL),\qquad \sum_j\sqrt{Q_jL}=O(\sqrt{NL}).
```

Each layer uses a constant number of exact R91 queries, with an
exclusive $`\Theta(b)`$ pool disjoint from all live inputs and
standalone work. Correlated old programs and all dirty banks return
in the complete comparison. Standalone cost over the whole tail is
$`M2^{O(s+1)}=2^{O(L+1)}`$, covered by the same small choice of
$`c_*`$. This tail has the same T and Clifford bounds as the grouped
tail. The sum of its width-dependent depths is $`O(NL/b^2)`$.
The rest is $`O(M[\ell+\log(L+2)])=O(\log^2(L+2))`$ since
$`L\ge\sqrt n`$. This pays every original lower-suffix guard with
its own geometry and preserves the optimal-order query counts.

## 7. Complete resources, accuracy and the complex phase

The opening and the chosen tail have the stated simultaneous
T and Clifford bounds. Their standalone banks, including each live
slot and precision program, fit a fixed fraction of b for
$`b\ge\sqrt N/n`$ and all sufficiently large n. For the uniform
choice fix all synthesis, guard and recovery precision constants
and query counts before $`c_*`$. Throughout $`L\le2c_*n`$, every
singly exponential standalone term with its polynomial factors
admits a $`2^{n/4}`$ majorant for small enough $`c_*`$. The additional
grouped catalogue terms below the square-root scale are $`2^{o(n)}`$.
The exact query pools are separate from those arrays and the logical
fields and can be reused through completed stages. This gives the
claimed width with no new external initialized bits.

Each branch implements the target order of complete factors. The sign
reduction, at most three openings and the chosen tail use the
retained fractions of $`2^{-L}`$ with fixed further padding if needed.
For example the opening with signs contributes at most
$`4\cdot2^{-L}/512`$. The grouped code and slot errors fit the
retained choice of $`\kappa`$, and the bottom and individual tails
can each use $`2^{-L}/16`$ in their respective branches. The wide
real error is below $`2^{-L}/4`$. In the initialized-isometry hybrid,
every ideal prefix preserves the caller's initialized domain and
later actual suffixes remain unitary on all leaked components.
The norm includes every logical column, dirty input, reference,
and work return.

For the prescribed phase-dressed complex magnitude family, use the
real branch at L+1 and apply the zero-clean programmed
$`R_Z(\phi_x)`$ to the external flag, which is zero on the ideal
trajectory. Its ideal action there is the literal phase
$`e^{i\phi_x}`$ and returns the flag. The target flag and all work
are covered by its full-operator estimate. There are N rows and
$`O(L)`$ bits per cutoff program, with constantly many exact query
calls. The phase costs $`O(\sqrt{NL}+NL/b)`$ T gates, $`O(NL)`$
Cliffords, and $`O(NL/b^2+\ell+\log(L+2))`$ depth. Give it error at
most $`2^{-(L+1)}`$. The two initialized-isometry errors then sum
below $`2^{-L}`$, including phases on zero-magnitude leaves. Both
stages reuse the sole external clean bit. The shifted real
allocation is $`17((L+1)+n+7)=17(L+n+8)`$ and the shifted precision
lies within the internal range.

For $`L\le2\ell`$ the retained R96 depth already obeys the new
bound. For $`b\lt\sqrt N/n`$ use R96 at its original precision and
literal minimum, since $`NL/b^2`$ absorbs $`L+\ell`$. For the finitely
many dimensions below the fixed wide-branch cutoff use R96 at the
original L. L is bounded there in the new range, so a fixed constant
pays its depth. These fallbacks preserve both literal width
thresholds and error allocations.

All bounds hold for the actual reversed native word too. Reversing
Clifford gates in the positive-S alphabet changes elementary counts
by at most an absolute factor, including $`S^\dagger=S^3`$. For the
ideal $`U=W\otimes I_b`$, multiplying the forward error by
$`V^\dagger`$ on the left and $`U^\dagger`$ on the right gives
$`\|V^\dagger J_1-J_1 U^\dagger\|=\|VJ_1-J_1 U\|`$. This includes the
same initialized flag and returned dirty/reference inputs. The
arbitrary-table depth gap at fixed accuracy remains open.
[R100](GEOMETRIC_CHUNK_DEPTH.md) improves the depth at larger
precision by promoting short fixed-source chunks.

On eligible allocations with

```math
b\le\sqrt{NL/(\ell+\log_2^2(L+2))},
```

the worst-case native T-depth is $`\Theta(NL/b^2)`$. The upper bound
is R99. The inherited
[one-clean width-dependent T-count lower bound](../../docs/OPEN_PROBLEM.md#optimal-banked-count)
is $`\Omega(NL/b)`$ under these literal allocations. A native band
has at most n plus b plus one physical targets, which is $`O(b)`$
at the eligible widths. Dividing gives the matching lower bound.
