# A carried source for heterogeneous component precisions

[Precision gauges](PRECISION_GAUGE_COMPILER.md) · [Disjoint reflection sources](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources) · [Carried-source history](../../docs/FAULT_TOLERANT_COMPILER.md#8-shared-private-work-history-and-the-final-normalization)

The mixed-height family in the precision-gauge chapter has a linear native
compiler with the original two clean qubits and dirty reservation. Its
component exponents differ by an unbounded ratio. The construction prepares
one class-dependent, small-ratio geometric source, retains it through an
entire coarse row network, and returns it in one final amplification. It
does not reschedule the unchanged high-precision merge banks.

The theorem below concerns this family only. In particular, the exponent
class is unchanged by every operation inside one component. No unrestricted
linear-count compiler follows from the argument.

## 1. Family, contract, and the upper sector

For integer t, use exactly the [mixed-height family](PRECISION_GAUGE_COMPILER.md#mixed-heights-and-the-limit-of-rescheduling):

```math
t\ge2,\qquad n=2^{2t},\qquad N=2^n,\qquad k=t,
\qquad h_r=\frac{n}{2^{r+1}},\qquad E_r=\frac{N}{h_r}
\quad(0\le r\lt k).
```

Each class occupies its prescribed aligned interval, partitioned into
complete aligned height-$`h_r`$ components. All crossing rotations are
zero. Every nonzero smaller branch in class r is independently specified
within $`[2^{-E_r-1},2^{-E_r}]`$. Unused labels are retained as identity
components. The supplied normalized chart has positive cosine and sine;
the already priced signed normalizations may also be attached if needed.

**Theorem.** Every such prescribed complete real frame W has an effective
coherent Clifford+T compiler satisfying

```math
T=O(N),\qquad
\|VJ_2-J_2(W\otimes I_b)\|\le2^{-N},
\qquad b=N+n+7.
\qquad\text{(1)}
```

The constant is independent of t and all branch values. The norm includes
every logical column, both initialized outputs, arbitrary dirty inputs,
and external references. All inverses are actual circuit inverses.
Certified classical evaluation and table construction are preprocessing;
no prepared precision state, quantum evaluator, reset, or measurement is
supplied.

Choose one sufficiently large universal positive integer C and put

```math
K=C\lceil\log_2(n+2)\rceil,\qquad N'=N/2^K.
\qquad\text{(2)}
```

All asymptotic arguments below are for a fixed threshold beyond which
$`K\lt \sqrt n=\min_r h_r`$ and the displayed workspace inequalities hold.
The finitely many smaller n are handled at the end. Increase C once to
cover the fixed clean-work constants and the polynomial table factors
identified below; no particular numerical value of C is asserted.

Separate the bottom K binary depths from the upper frame, so
$`W=W_{\rm bottom}W_{\rm up}`$ and the upper part acts first. The bottom
part is a forest whose precision trace is at most

```math
O\bigl(K(E_{k-1}+\log n)\bigr)
=O\left(\frac{N\log n}{\sqrt n}+\log^2n\right)=O(N).
\qquad\text{(3)}
```

The precision-trace theorem therefore compiles it in linear T count.
The upper frame acts only when the last K logical bits are zero. On that
sector it is a frame on $`N'`$ labels, with component heights $`h_r-K`$.
Standard aligned labels are retained throughout; neither a packing map
nor a feature assignment changes a class.

Compute $`h=[\text{last K bits}=0]`$ into the first external clean qubit
using the exact borrowed predicate construction. Keep h unchanged until
the entire upper compiler, including its fine-gauge recovery, is complete.
On $`h=1`$ the suffix supplies K initialized private wires. The second
external clean qubit is available as a feature-mode flag. Every completed
subroutine is conditioned to be exactly identity on $`h=0`$, on the
whole space of the remaining registers, using the
[conditional-suffix rules](../../docs/CONDITIONAL_SUFFIX_COMPILER.md#1-an-active-suffix-supplies-workspace-conditionally).
Finally reverse the actual predicate word. A complete-isometry hybrid
charges any suffix leakage; no intermediate reset is used.

## 2. Exact coarse rotations and their full row network

Write a target normalized node rotation as

```math
R_y(\theta)=\begin{bmatrix}c&-s\\s&c\end{bmatrix},
\qquad \tau=\tan(\theta/2)=\frac{s}{1+c}.
```

For a node v in class r choose a positive dyadic $`u_v`$, with m fractional
bits, approximating $`2^{E_r}\tau_v`$, where

```math
m=O(\log n),\qquad 2^{-m}\le\frac1{16000n^2},
\qquad \frac18\le u_v\le1.
```

These bounds hold uniformly for sufficiently large n, by certified
rounding of $`2^{E_r}\tau_v\in[1/4,1]`$. Put $`t_v=2^{-E_r}u_v`$.
The exact mathematical coarse gate is the dyadic Cayley rotation

```math
J=\begin{bmatrix}0&-1\\1&0\end{bmatrix},\qquad
R(t)=(I+tJ)(I-tJ)^{-1}
=\frac1{1+t^2}\begin{bmatrix}1-t^2&-2t\\2t&1-t^2\end{bmatrix}.
\qquad\text{(4)}
```

The derivatives of its first column, together with
$`s_v\ge2^{-E_r-1}`$, give relative error at most
$`\eta_{\mathrm{rel}}=1/(1000n^2)`$ in both nonzero entries, after the fixed
rounding margin above. All coefficients in (4) are exact computable
rationals. This defines the ideal coarse frame $`W_0`$ and the regular
gauge; it does not assert that (4) is an exact finite native gate.

Use the [regular Parseval network](PRECISION_GAUGE_COMPILER.md#2-the-full-column-regular-gauge-and-its-parseval-source)
for these gates. At each depth, first merge the two carriers with
$`R(t_v)^\dagger=R(-t_v)`$, then split the parent carrier against its
fresh feature mode with an exactly native splitter. Choose its literal
first column $`(r_d,t_d^{\rm split})`$ so that

```math
\frac1{16n}\le |t_d^{\rm split}|^2\le\frac1{8n}.
```

Finite native single-qubit approximants with polynomial-in-log(n) word
length suffice, using a strict interior target. Crossing and inactive
stages have their splitter disabled. Products of the actual $`r_d`$ are
computed separately for every component root and node. Since every path
has at most n stages, the resulting coefficients obey the same uniform
bounds as the regular-gauge construction:

```math
Q_0|x,0\rangle
=\begin{bmatrix}\Gamma W_0^\dagger x\\
                 \Omega T_0^\dagger x\end{bmatrix},\qquad
|\gamma_j|^2\ge7/8,\qquad
|\omega_v|^2\ge7/(128n).
\qquad\text{(5)}
```

Here $`T_0`$ contains the actual coarse internal-root states, and the
feature flag supplies their modes. The physical merge/split circuit
defines $`Q_0`$ as a unitary on all $`2N'`$ modes, including inputs with
occupied feature modes. That full unitary is the target compiled below.

For clarity, if p and m are the left and right carrier labels of a node,
the merge acts on $`(0,p),(0,m)`$ and its splitter on
$`(0,p),(1,m)`$, where the first bit is the feature flag. No crossing
node receives a splitter. In standard aligned labels, p and m differ
only in that node's target bit. Conjugating by CNOT from the feature flag
to that target bit maps the splitter pair to a feature-bit rotation
guarded by target zero and the unchanged logical suffix. This supplies
an explicit native matching layout on all modes, with the conjugating
CNOTs cancelling on other pairs. The class is the existing
most-significant interval prefix, of $`\lceil\log_2k\rceil`$ bits; it
is never copied to a fresh clean exponent register.

Let $`C_d`$ denote the exact splitter layer. Let $`\mathcal R_d`$ be its
preceding merge layer. On each active carrier pair use generator
$`J_d=-J`$, so its merge is $`R(-t_v)`$; feature and frozen modes have
identity merge. Pad all components through one common list of
$`R\le n`$ depths. Thus

```math
Q_0=C_R\mathcal R_R\cdots C_1\mathcal R_1.
\qquad\text{(6)}
```

Every factor preserves the class. The same construction can include a
preserved row-attenuation flag: on its failure value set all polynomial
coefficients to zero and disable every splitter. The target on that
sector is then exactly I, while the normalization below stays uniform.

## 3. One finite source with a class-dependent small ratio

Because n is dyadic, all of the following choices are exact:

```math
P=n,\qquad M=2n,\qquad p=\frac1{8n},\qquad
B=6+2\log_2n,\qquad \beta_r=2^{B-E_r}.
\qquad\text{(7)}
```

For sufficiently large n, $`0\lt \beta_r\le1/2`$. On an M-mode source
register define

```math
g_\beta=a_\beta\sum_{j=0}^{M-1}\beta^j|j\rangle,
\qquad
a_\beta=\sqrt{\frac{1-\beta^2}{1-\beta^{2M}}},
\qquad
S_\nu=\sum_{j=\nu}^{M-1}|j-\nu\rangle\langle j|.
\qquad\text{(8)}
```

For $`1\le\nu\le P\lt M`$, only the final $`\nu`$ coordinates contribute
to the residual. Direct summation gives

```math
\|S_\nu g_\beta-\beta^\nu g_\beta\|
=\beta^M\sqrt{\frac{1-\beta^{2\nu}}{1-\beta^{2M}}}
\le2\beta^M.
\qquad\text{(9)}
```

This is a full-output relation on the source column, not just a projected
overlap. The source ratio depends on the immutable class. Define
$`E_g`$ to attach $`g_{\beta_r}`$ to every system mode in class r,
using the same source state on carriers, features, and any preserved
row-attenuation control. Unused labels may use the first class's source.
All system operations below preserve this assignment.

### Charged preparation and its full-space extension

For each fixed reduced system label x choose the real source extension

```math
v_x=\frac{e_0+g_{\beta(x)}}{\|e_0+g_{\beta(x)}\|},
\qquad A_x=R_{v_x}R_{e_0},\qquad A_xe_0=g_{\beta(x)}.
\qquad\text{(10)}
```

Each reflection is confined to that label's M source modes. Including
feature and row-attenuation bits, the complete mode dimension is at most
$`4N'M\le N`$ once C and the asymptotic threshold are fixed. Including
the external conditioning bit still gives at most $`2N`$ addressed modes. The
[one-clean disjoint-reflection interface](DISJOINT_BLOCK_COMPILER.md#4-disjoint-reflection-and-state-sources)
therefore implements A with four fine diagonal calls and one exact center
sign. Denote its emitted word by $`\widehat A`$. It has cost $`O(N)`$
and full initialized-isometry error $`\delta_A`$
of any fixed chosen multiple of $`2^{-N}`$.

There is no need for a class-changing coarse preparation. Indeed
$`\|v_x-e_0\|=O(\beta_{max})\le1/(16M)`$ for sufficiently large n,
so the reflection interface's local state chart takes its coarse word
to be identity. Its Walsh interlayers act only on source bits, and its
diagonals and center signs preserve x. Thus the emitted word preserves
the system and class exactly, even when its private compiler flag leaks.
Its inverse below is its actual native inverse. The flag is distinct
from all private flags used by the carried stream.

## 4. Exact dyadic local kernels

The Cayley series in a carrier pair is

```math
R(-t_v)=I+2\sum_{\nu=1}^{\infty}t_v^\nu J_d^\nu.
\qquad\text{(11)}
```

For $`1\le\nu\le P`$, define

```math
f_{v,\nu}
=\frac{2P(1-p)}p\,u_v^\nu2^{-B\nu}
=2n(8n-1)u_v^\nu2^{-B\nu}.
\qquad\text{(12)}
```

Set f to zero on every frozen or feature mode and on the disabled
row-attenuation sector. On a carrier pair the same f is assigned to both
labels. These numbers are exact dyadics, obey $`0\le f_{v,\nu}\le1/4`$,
and need at most

```math
q_{\rm coef}\le P(B+m)+O(1)=O(n\log n)
\qquad\text{(13)}
```

binary fractional digits. The
[raw operator-source scalar](../../docs/OPERATOR_SOURCE_COMPILER.md#2-exact-signed-dyadic-coefficients)
therefore realizes an exact accepted block $`f_{v,\nu}I`$ with one
retained scalar flag and $`q_{\rm coef}+O(1)`$ arbitrary dirty core
wires. Enlarge the finite sign grid by a fixed number of places when
necessary. This is the unamplified scalar word; no inverse polynomial
or unrounded real coefficient is being encoded.

Prepare a uniform P-valued label $`\ell`$ and interpret
$`\nu=\ell+1`$. Prepare $`\log_2(8n)`$ Hadamard bits and designate
their all-zero string as a rare event, of probability p. In the common
event SELECT is identity. In the rare event apply the scalar word for
f, the shift dilation for $`S_\nu`$, and the literal quarter-turn
$`J_d^\nu`$ on the selected carrier pair. Its phase is retained; modulo
four powers suffice. The shift dilation toggles a distinct underflow
flag by $`[j\lt \nu]`$ and then subtracts $`\nu`$ modulo M. Projecting
that flag to zero is exactly $`S_\nu`$ on every source input.

The source-underflow and scalar-rejection flags are distinct. Both are
retained. The rare bits and power label are unprepared by the actual
Hadamard inverses. All coefficient queries return their borrowed
selectors exactly. Their address values are unchanged during the query;
the subsequent quarter-turn stays in the same carrier pair, on which f
is constant. Every SELECT case specifies a full unitary, including
rejected and inactive inputs.

Let F consist of these rare bits, the power label, both retained flags,
and any retained private work. Exact arithmetic scratch is restored at
the boundary. The actual kernel $`K_d`$ has F-zero block

```math
T_d=(1-p)I+\frac pP\sum_{\nu=1}^P
S_\nu\otimes F_{d,\nu}\mathcal J_d^\nu.
\qquad\text{(14)}
```

Here $`F_{d,\nu}`$ is the diagonal f table; $`\mathcal J_d`$ is the
direct sum of $`J_d`$ on active carrier pairs and I elsewhere. The order
is immaterial within each pair. Formula (14) follows by projecting the
two separate dilation flags, not by dropping a rejected term of a
later operation on the same flag. Since it is a compression of a
unitary, $`\|T_d\|\le1`$.

Use $`t_v=\beta_r u_v2^{-B}`$ in (11). With
$`D_d=(1-p)\mathcal R_d`$, equations (9), (12), and (14) give

```math
\|T_dE_g-E_gD_d\|
\le 2p\beta_{max}^M
 +\frac{2(1-p)t_{max}^{P+1}}{1-t_{max}}
=:\delta_{\rm ker}.
\qquad\text{(15)}
```

The bound holds on every system and dirty input, with references.
Frozen and disabled sectors have exactly $`(1-p)I`$ as their kernel;
they do not acquire a different attenuation. Also
$`\|D_d\|=1-p`$. Since $`E_{min}=2N/n`$, increasing the fixed
threshold so that $`2nB\le N`$ gives

```math
M(E_{min}-B)\ge3N,\qquad
t_{max}^{P+1}\le2^{-2N},\qquad
R\delta_{\rm ker}=o(2^{-N}).
\qquad\text{(16)}
```

## 5. Whole-network history, normalization, and return

Share F among all R kernels. Initialize a history counter H of
$`\lceil\log_2R\rceil`$ bits, omitting it if $`R=1`$. After each
kernel except the last, increment H if and only if $`F\ne0`$. There
are at most $`R-1`$ increments. Hence a path beginning at H=0 cannot
wrap back to zero. Later kernels never act on H.

Expand in the F and H computational bases. Final H=F=0 retains exactly
those paths whose F value was zero at every intervening boundary. It
therefore inserts the intermediate projectors, giving the exact accepted
operator on carried source and system

```math
V_{\rm acc}=(I\otimes C_R)T_R\cdots(I\otimes C_1)T_1.
\qquad\text{(17)}
```

Each $`C_d`$ is defined on every occupied feature input and acts only on
the system. It preserves the source class. Thus
$`(I\otimes C_d)E_g=E_gC_d`$, and contractive telescoping of (15)
proves

```math
\|V_{\rm acc}E_g-E_g\gamma Q_0\|
\le R\delta_{\rm ker},\qquad
\gamma=(1-p)^R\ge1-Rp\ge7/8.
\qquad\text{(18)}
```

This estimate includes the carried source's residual. It makes no claim
that an intermediate source or F register has returned on the actual
state. Failed branches continue through defined full unitaries.

Apply $`\widehat A`$ once, run the complete stream, and apply the actual
native inverse $`\widehat A^\dagger`$.
Comparison with the ideal source extension (10), using a unitary hybrid,
adds at most $`2\delta_A`$ to the final zero-source block error. The
source-preparation compiler flag is included among final tested work.

Use one additional retained scalar flag, outside F and its history test,
to encode

```math
w=\frac1{2\gamma}\in[1/2,4/7]
```

by a raw operator-source scalar word. With $`q+1`$ dirty generators its
signed grid has spacing $`2^{1-q}`$. Since w is an exact rational,
rounding to that grid gives $`\widetilde w`$ with
$`\delta_w\le2^{-q}=O(2^{-N})`$ for $`q=N+c_0`$. This one fine
call costs $`O(N)`$. Its accepted operator is
$`\widetilde w I_{\rm dirty}`$, so it multiplies the completed stream
block even when its dirty core is shared. Its failure flag never
re-enters any other scalar subroutine. If U is the whole actual word,
and J initializes all its private work, its accepted block satisfies

```math
\|J^\dagger UJ-Q_0/2\|
\le R\delta_{\rm ker}+2\delta_A+\delta_w.
\qquad\text{(19)}
```

The [normalization-two amplification lemma](../../docs/FAULT_TOLERANT_COMPILER.md#4-from-a-half-unitary-block-to-the-complete-output)
applied to $`-URU^\dagger RU`$, with $`R=I-2JJ^\dagger`$, now gives
the full initialized-isometry approximation to $`Q_0`$. It includes
every private output, all dirty disturbance, and arbitrary occupied
feature inputs. There are two forward words, one actual inverse, and
two work reflections, so only a constant number of fine calls.

For the external upper-sector marker use the literal conditional form

```math
Z_h U_h R_h U_h^\dagger R_h U_h,
\qquad R_h=I-2\Pi_{h=1,\,\mathrm{private}=0}.
\qquad\text{(20)}
```

On $`h=1`$ this is the required leading-minus word. On $`h=0`$ every
factor, including its complete output, is exactly identity. The
reflection does not test arbitrary logical, feature, or dirty inputs.
When compiling the row-controlled version, its preserved attenuation
control is likewise system data, not initialized work.

## 6. Fine regular gauge and the final frame

Apply the regular full-column identity to the exact coarse gates (4):

```math
W_{\rm up}=D(W_0C+T_0A).
\qquad\text{(21)}
```

The relative bounds of Section 2 give the regular gauge's uniform
contraction slack. With $`\Lambda=5/4`$ and $`\kappa=8/5`$, the
ideal row source first attenuates by $`D^\dagger/\Lambda`$, routes its
failure to one bank, and applies $`Q_0`$ on the shared sector. Its shared
columns are

```math
\frac1\Lambda
\begin{bmatrix}\Gamma W_0^\dagger D^\dagger\\
\Omega T_0^\dagger D^\dagger\end{bmatrix}.
```

The column source has shared coefficients
$`C_j/(\kappa\bar\gamma_j)`$ and
$`A_{v,j}/(\kappa\bar\omega_v)`$, with the opposite failure bank.
Their overlap is exactly $`W_{\rm up}/2`$ by (21). These are the
existing regular-gauge fine banks, not a scalar replacement for an
unknown dirty rejection operator.

The feature flag is the second external clean qubit. The row/column
attenuation flag is a distinct active-suffix wire. The first external
flag remains h throughout. All additional source, history, normalization,
and compiler flags are separate active-suffix wires. In particular no
occupied attenuation flag is reused as a clean source-preparation flag.

Replace the ideal row-network call by the full-unitary approximation of
Section 5, controlled on the shared attenuation value. On its other
value the target is I, as specified in Section 2. Compile the constant
number of remaining fine banks at precision $`q=N+c_0`$, for a fixed
sufficiently large $`c_0`$. An ordinary complete-isometry hybrid bounds
these replacements. Approximate return in an earlier call is propagated
unitarily; later calls are compared with their ideal inputs and no flag
is reset. Apply the regular-gauge outer normalization-two recovery, with
the conditional sign convention (20). Its reflection tests every
initialized private and suffix wire and the feature flag; h is the
preserved enable and is excluded from the zero-work test. It does not
test any reduced logical or dirty input. Fixed error fractions give an
upper-frame error below $`2^{-N}/2`$.

Finally reverse the upper-sector predicate and compose with the bottom
compiler from (3), also at error below $`2^{-N}/2`$. This proves the
complete error claim (1). The coarse mathematical gates are used to
define computable gauge data; only their globally compiled row network
is executed. There is no requirement that the carried-source word
factor into exact native implementations of those local gates.

## 7. Native count and literal workspace

At one stage the complete mode-and-power coefficient table has at most
$`O(N'P)`$ rows, including fixed feature and attenuation controls.
Whole-word self-borrowed queries implement its dyadic sign masks with
T cost $`O(N'P+q_{\rm coef})`$. The output mask width affects the
Clifford work, not the row contribution to T. The source-word constant
is fixed. Exact shift arithmetic, rare predicates, quarter-turn powers,
and history operations have a fixed polynomial cost in n. Even the
conservative full-table ledger over $`R\le n`$ stages is

```math
T_{\rm kernels}
=O\left(n^2\frac{N}{2^K}+\mathrm{poly}(n)\right)=O(N).
\qquad\text{(22)}
```

The exactly native splitters have polynomial-length coarse descriptions.
Their controlled matching tables and suffix predicates cost
$`N'\mathrm{poly}(n)+\mathrm{poly}(n)`$, absorbed by a
single sufficiently large choice of C. Source preparation and its
inverse, global normalization, and all fine gauge banks each have
$`O(N)`$ cost and occur a constant number of times, even after the two
amplifications. Equation (3) charges the bottom part. This proves the
T count in (1).

The simultaneous active clean reservation consists of a source register
of $`\log_2M=O(\log n)`$ bits, a power register and its increment
scratch, $`\log_2(8n)`$ rare bits, $`\lceil\log_2R\rceil`$ history
bits, arithmetic/predicate scratch of $`O(\log n)`$ bits, and a fixed
number of distinct scalar, underflow, preparation, normalization, and
fine-gauge flags. This totals $`O(\log n)`$ with a fixed constant,
covered by K. Long logical suffix predicates use the existing returned
borrowed construction, not an n-bit clean register.
At a fixed depth, active-node predicates require only class membership
and the unchanged logical suffix. Quarter-turn and splitter controls
therefore use polynomially many exact gates and returned borrowed
predicate work, rather than initialized logical-endpoint metadata.

During fine calls reserve one shared core of $`q+1=N+c_0+1`$ arbitrary
dirty wires. The largest fine address has at most n+1 bits because its
mode dimension is at most $`2N`$. Self-borrowing applies once
$`q+1\ge2(n+1)`$. All bounded extra controls and exact center predicates
reuse returned core helpers in completed subroutines. A fixed number of
external dirty helpers, if used, fits in the remaining $`n+6-c_0`$
wires for sufficiently large n. Coarse reflection routing uses
$`O(\sqrt N+\mathrm{poly}(n))`$ dirty space instead, and the
dyadic kernels use $`O(n\log n)`$ dirty source wires. These reservations
are sequential and reuse the same pool. No active-suffix bit is also
counted as dirty storage. Thus the literal maximum is $`N+n+7`$.

All finitely many n below the chosen threshold have fixed dimension and
fixed required accuracy. The existing direct determinant-one,
literal-phase finite-net fallback compiles every frame in those compact
families with the same two flags and dirty pool. The finite maximum of
T/N is absorbed into the universal constant. This is an asymptotic
theorem, not a proposed practical crossover threshold.

## 8. Scope and verification

The exponentially small shift eigenvalue in (9) permits endpoint
accuracy with only $`M=2n`$ initialized source modes. Each component's
immutable exponent class is essential to (18). The argument does not
apply when the required scale changes inside a component, and it does
not pay for node-dependent fine source preparations.

The [exact helper](../../compiler_robust_hopf/heterogeneous_precision_source.py)
records finite source residuals, dyadic Cayley coefficients, and symbolic
parameter margins. The
[focused tests](../../tests/test_heterogeneous_precision_source.py)
check the finite shift residual, exact dyadic coefficients, truncated
kernel columns, separate scalar/underflow rejection flags, and a
heterogeneous no-wrap history product with uniform attenuation. They
include occupied inputs and a nontrivial dirty rejection. Those bounded
checks verify the displayed interfaces; they do not
establish the asymptotic count, replace the source-preparation charge,
or constitute an unrestricted native frame emitter. The native theorem
follows from the exact kernels, contractive recurrence, fixed fine-call
count, and workspace ledger above.
