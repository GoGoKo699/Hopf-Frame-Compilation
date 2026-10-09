# Collective precision refinement of regular-core frames

A single all-node rotation bank can correct the first-order defects of
every depth of a regular-core frame simultaneously. Exact geometric
histories fit in two supplied clean qubits. Combining this compression
with physical purification gives quadratic accuracy; a midpoint
conjugator and opposite-sign filters give cubic accuracy.

These are component theorems for the complete terminal frame of the
regular-core coins in the
[coarse prefix encoder](COARSE_PREFIX_ENCODER.md#3-the-canonical-su2-native-circuit):

```math
Q_v(\alpha_v)=R_y(\pi/2)R_x(\alpha_v)R_z(\alpha_v),
\qquad R_k(\theta)=e^{-i\theta\sigma_k}.
```

They do not assert that this auxiliary frame is the prescribed real Hopf
target. For $`N=2^n`$, $`n\ge16`$, and $`s=\lceil N/n\rceil`$,
the cubic construction has

```math
\begin{aligned}
\|\widehat WJ-J(W_{\rm target}\otimes I_b)\|
 &\le40\,8^{-s}+2^{-L},\\
T&\le C_1(N+3ns+n^4)+C_2(N+L)=O(N),\\
G&=O\!\left(N(s+n+L)+\mathrm{poly}(n)\right),\\
b&\le N+n+7.
\end{aligned}
```

The fine-precision parameter satisfies $`3s\le L\le N`$.
Here J initializes exactly two flags, b counts arbitrary dirty qubits,
and G counts the remaining native gates. All errors are literal-phase
operator norms, include the complete initialized output and work return,
and hold with arbitrary reference entanglement. The alphabet is
H,S,CNOT,T,T-dagger; preprocessing is finite and certified. The constants
are absolute. Directly tripling the coarse source precision already has
the same asymptotic cost, so this theorem does not improve the retained
endpoint bound. Its content is the collective correction and its complete
physical implementation.

## 1. Rounded coins and the actual coarse word

For one addressed real primitive, write its exact pre-amplification
compression as $`B_x=u_xI+v_xXZ`$, and set

```math
\rho_x=(u_x^2+v_x^2)^{1/2},\qquad
R'_x=B_x/\rho_x,\qquad
p_5(z)=5z-20z^3+16z^5.
```

The [literal five-call word](../../docs/ONE_CLEAN_COMPILER.md#4-five-calls-amplify-with-one-flag)
has accepted block $`p_5(\rho_x)R'_x\otimes I_b`$. Its canonical
programming bound at precision q gives

```math
\begin{aligned}
\|\langle0_c|V_x|0_c\rangle-R'_x\otimes I_b\|
 &\le94\,4^{-q},\\
\|V_x-I_c\otimes R'_x\otimes I_b\|
 &\le43\,2^{-q}.
\end{aligned}
```

The second estimate uses the same exact signal-X symmetry as the
[borrowed-signal proof](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations),
with the rounded polar rotation as its comparison. The scalar
$`p_5(\rho_x)`$ is positive on the programming interval. The coefficients
are algebraic and $`\rho_x`$ is bounded away from zero, so this
comparison uses no angle-extraction or equality oracle. Fixed target
Clifford conjugations preserve both bounds. Reusing the same coefficient
table for the X and Z axes preserves the correlated one-parameter form
of every rounded coin $`Q_v^0`$.

Let $`W_0`$ be the exact terminal frame of these coins, and let
$`\widehat C_0`$ be the actual coarse circuit. On initialized $`a=0`$,
the wrapper $`X_aT_{\rm pre}X_a`$ selects heap depth n and returns a,
so its ideal action is exactly $`W_0`$; the two flips are Clifford gates.
For its constituent real banks put
$`e_i=43\,2^{-q_i}`$ and $`e=\sum_i e_i\le1/2`$. Product expansion gives

**Coarse accepted block (1).**

```math
\|J^\dagger\widehat C_0J-W_0\otimes I_b\|
\le\sum_i e_i^2+\exp(e)-1-e\le2e^2.
```

Every one-error compression is at most $`e_i^2`$, while terms containing
at least two errors sum to at most $`\exp(e)-1-e`$. Exact routing and fixed
factors can be absorbed into the ideal factors. There is no intermediate
return assumption. Full-space telescoping gives error at most e, and
the canonical schedule yields

```math
q_j=s+n-j+7,\quad 0\le j\lt n,\qquad
e\lt \frac{43}{64}\,2^{-s}.
```

Thus $`W_0`$ is an exact comparison frame, not an exactly returned native
subroutine. Products of signal-odd errors can produce second-order
signal-even dirty operators; they cannot be discarded by replacing the
actual coarse word with its logical polar factors.

## 2. A literal two-flag purification lemma

Let $`P=JJ^\dagger`$, let C be any physical unitary, and suppose
$`\|J^\dagger CJ-U\|\le a\lt 1`$ for a unitary U on logical and dirty
inputs. Define

**Two-flag filter (2).**

```math
\begin{aligned}
D&=\exp\!\left[\frac{i\pi}{3}(P-I/4)\right],\\
\mathcal F_D(C)&=D^{-2}CDC^\dagger DC.
\end{aligned}
```

Writing $`J^\dagger CJ=Vh`$ in polar form, block multiplication gives

```math
J^\dagger\mathcal F_D(C)J=Vf(h),\qquad
f(x)=x\,[1+e^{-i\pi/3}(1-x^2)].
```

For $`0\le x\le1`$, $`2-2\mathrm{Re}f(x)=(1-x)^2(x+2)`$.
Unitarity of the complete filtered word therefore gives

**Full-output bounds (3).**

```math
\begin{aligned}
\|\mathcal F_D(C)J-JV\|&\le\sqrt3\,\|I-h\|,\\
\|\mathcal F_D(C)J-JU\|&\le(2+\sqrt3)a.
\end{aligned}
```

Indeed $`\|I-h\|\le a`$ and $`\|V-U\|\le2a`$. The first line
includes the rejected rows, rather than merely bounding the accepted
block. This argument is independent of the dirty dimension.

The phase is a native synthesis task, with no discarded scalar:

```math
D=\exp\!\left[\frac{i\pi}{12}(Z_a+Z_c+Z_aZ_c)\right].
```

Its three commuting Pauli rotations are fixed Clifford conjugates of
the full-operator real primitive. A logical wire can serve as its
arbitrary synthesis signal. Use the chosen actual D word and its actual
inverse, twice for $`D^{-2}`$. The four phase appearances cancel the
net scalar, and their synthesis errors add at most four times the
full-space D error. All flags may be occupied during these calls.

## 3. Exact geometric histories

Use physical mode index $`c(2N)+aN+x`$. For $`1\le\ell\le n`$, let
$`F_\ell`$ send $`|m0^{n-\ell}\rangle`$ to heap mode
$`|2^\ell+m\rangle`$ and vanish on other columns. Set

```math
\begin{aligned}
\Pi_\ell&=F_\ell^\dagger F_\ell,&w_\ell&=2^{\ell-n-1},\\
E&=\sum_{\ell=1}^n\sqrt{w_\ell}\,F_\ell
       \ \oplus\ P_{\rm pad},&
P_{\rm pad}&=\left(I-\sum_{\ell=1}^nw_\ell\Pi_\ell\right)^{1/2}.
\end{aligned}
```

History occupies heap modes $`2,\ldots,2N-1`$, and padding occupies
$`2N,\ldots,3N-1`$. Thus $`E^\dagger E=I`$ within 4N modes. Define
$`d(0)=0`$ and $`d(x)=n-1-v_2(x)`$ for nonzero x. The eligible levels
are $`\ell=d(x)+1,\ldots,n`$, and the exact columns are

**History columns (4).**

```math
E|x\rangle=2^{-(n-d(x))/2}|2N+x\rangle
 +\sum_{\ell=d(x)+1}^n
  2^{-(n-\ell+1)/2}|2^\ell+x/2^{n-\ell}\rangle.
```

An exact native unitary $`\mathcal H`$ realizes $`\mathcal HJ=E`$.
First apply $`X_c`$ to send every initialized column to padding. Process
levels $`n,n-1,\ldots,1`$. At level $`\ell`$, pair the modes

```math
|c=1,a=0,m0^{n-\ell}\rangle,
\qquad |c=0,0^{n-\ell}1m\rangle.
```

In the $`c=0`$ half, route the heap delimiter to a, the path bits to the
high logical positions, and the other bits to the low positions; flip a.
Leave the $`c=1`$ half fixed. Apply $`R_y(-\pi/4)`$ to c controlled
on $`a=0`$ and the zero low suffix, then undo the route. In the order
$`c=0,c=1`$ this splits incoming padding amplitude positively and equally
between history and padding. Other depths and ineligible columns are
untouched. Descending equal splitting proves (4).

Every route is a full basis permutation, including occupied inputs, made
of $`O(n)`$ controlled swaps and one controlled X. For the fixed rotation
use $`R_y(-\pi/4)=ZH`$. A controlled Hadamard is
$`A(\text{controlled }Z)A^\dagger`$ with $`A=R_y(\pi/8)`$; the
scalar of its Clifford-conjugated T implementation cancels against the
actual inverse. The
[exact borrowed-helper MCX](../../docs/BORROWED_WORKSPACE_COMPILER.md#3-an-exact-echo-selects-a-logical-sector)
implements the multiply controlled Z at $`O(n^2)`$ T cost per level.
Hence $`\mathcal H`$ and its inverse cost $`O(n^3)`$ T and return one
arbitrary helper exactly. A source-pool wire supplies that helper while
the source is inactive. There is no precision state or additional clean
history register.

## 4. Collective quadratic correction

Let $`B_\ell^0,B_\ell^1`$ be the rounded and requested coin banks at
depth $`\ell-1`$. Write $`W_\ell^0=B_\ell^0A_\ell`$ on the active
$`2^\ell`$-dimensional sector and define

```math
X_\ell=\log[(B_\ell^0)^\dagger B_\ell^1],\qquad
G_\ell=F_\ell^\dagger A_\ell^\dagger X_\ell A_\ell F_\ell.
```

The close SU(2) coins have uniformly nonsingular principal logarithms.
Their algebraic sine/cosine data permit certified evaluation. Peeling
layers from the terminal end gives the exact orientation

**Relative-frame identity (5).**

```math
R=W_0^\dagger W_{\rm target}
   =e^{G_n}\cdots e^{G_1}.
```

On history depth $`\ell`$, put $`S_\ell=A_\ell`$, with S identity on
padding and unused modes. It is the rounded prefix encoder followed in
time by the inverse all-node rounded child bank. Set

```math
\begin{aligned}
M&=\bigoplus_{\ell=1}^n e^{X_\ell/w_\ell}\ \oplus I,\\
C&=\mathcal H^\dagger S^\dagger MS\mathcal H,\qquad
Z=J^\dagger CJ,\\
Z&=I+\sum_{\ell=1}^n w_\ell F_\ell^\dagger A_\ell^\dagger
       (e^{X_\ell/w_\ell}-I)A_\ell F_\ell.
\end{aligned}
```

If $`\|X_\ell\|\le\kappa w_\ell`$, then
$`\|M-I\|\le\kappa`$ on the full space. Both Z and R have linear
term $`\sum_\ell G_\ell`$. The skew-Hermitian exponential integral
remainder, and the corresponding ordered-product remainder, give

**Compression error (6).**

```math
\|Z-R\|\le\frac{\kappa^2}{2}
 \left[\sum_\ell w_\ell+\left(\sum_\ell w_\ell\right)^2\right]
 \le\kappa^2.
```

Canonical programming gives $`\|X_{j+1}\|\le65\,2^{-q_j}`$; hence
$`\kappa\le(65/128)2^{-s}\lt 2^{-s}`$. One SU(2) middle bank, decomposed
into three addressed fixed-axis banks, contains all local derivatives.
For each block, search the dense dyadic triples of Euler angles until
certified evaluation bounds its distance from the target below
$`2^{-q}`$. Compact SU(2) Euler coverage and a strictly positive
tolerance ensure termination, including chart singularities. No exact
zero test is needed. Synthesize the three resulting fixed-axis banks
at precision q. Their full-operator error, including this classical
approximation, is at most $`(3\cdot43+1)2^{-q}=130\,2^{-q}`$.
All coefficient tables together have $`O(N)`$ rows.

If $`\widehat S`$ has full-operator error $`\delta_S`$ and the actual
middle bank has full-operator error $`\epsilon_M`$, the
[small-conjugation lemma](COARSE_PREFIX_ENCODER.md#4-full-output-error-under-a-small-conjugated-gate)
gives

**Native conjugation error (7).**

```math
\|\widehat C-C\|\le2\kappa\delta_S+\epsilon_M.
```

The inverse is the actual inverse. Applying (3) to this correction and
to (1), then composing the purified baseline on the left, gives

**Quadratic replacement (8).**

```math
\text{full initialized error}
 \le20\,4^{-s}+2^{-L},\qquad 2n\le L\le N.
```

For an explicit reservation, retain the original baseline schedule, use
$`q_j=s+n-j+9`$ in the controlled auxiliary prefix encoder and
$`q=s+10`$ in its child bank. Their errors sum to less than $`2^{-s}`$.
Use $`q=L+20`$ in each fine middle/phase axis. Each three-axis bank has
error at most $`130\,2^{-q}`$, and their total contribution in (3) is
at most $`[(2+\sqrt3)+8]130\,2^{-q}\lt 2^{-L}`$. The two purifiers use
three baseline calls, three correction calls, six auxiliary S calls,
three middle calls, six history calls, and eight phase calls, including
inverses. The S, middle and history counts expand the correction count.
Their T count is $`O(N+ns+n^4)+O(N+L)`$.

## 5. Midpoint velocity and cubic accuracy

Put $`d=2^{-s}`$ and let $`\Delta_v`$ be the short angle difference
from the rounded coin to the target. The scalar programming estimates
give $`\max_{v\text{ at depth }j}|\Delta_v|\le d\,2^{j-n}`$.
The relative cosine is positively bounded away from zero, so certified
arctangent evaluation of a small ratio selects this difference without
a sign, exact-angle, or zero oracle. Use the correlated path

```math
\alpha_v(t)=\alpha_v^0+t\Delta_v,\qquad
W(t)=\text{complete frame of }Q_v(\alpha_v(t)).
```

Its right logarithmic velocity is

```math
\begin{aligned}
G(t)&=-iW(t)^\dagger W'(t)\\
 &=\sum_{\ell=1}^n F_\ell^\dagger A_\ell(t)^\dagger
              H_\ell(t)A_\ell(t)F_\ell,\\
H_v(t)&=-\Delta_v\,[R_z(\alpha_v(t))^\dagger
                              X R_z(\alpha_v(t))+Z].
\end{aligned}
```

Thus $`R(t)=W(0)^\dagger W(t)`$ satisfies $`R'=iRG`$. If
$`B=2\sum_j\max_{v\text{ at depth }j}|\Delta_v|\le2d`$,
Leibniz expansion of the product of constant-generator exponentials gives

```math
\|W^{(r)}(t)\|\le B^r,\qquad
\|G\|\le B,\quad\|G'\|\le2B^2,\quad\|G''\|\le4B^3.
```

Midpoint quadrature of the first Dyson term costs at most $`B^3/6`$.
Replacing the second term by $`-G(1/2)^2/2`$ costs at most $`B^3/2`$:
integrate $`2B^3(|t-1/2|+|u-1/2|)`$ over its triangular domain. The
two higher-order tails total at most $`e^BB^3/3`$. For $`B\le1/4`$,

**Midpoint estimate (9).**

```math
\|R(1)-e^{iG(1/2)}\|\le2B^3\le16d^3.
```

Use the midpoint prefixes in S and half the velocity in the middle bank:

```math
C_{1/2}=\mathcal H^\dagger S_{\rm mid}^\dagger
 \left[\bigoplus_{\ell=1}^n
    e^{iH_\ell(1/2)/(2w_\ell)}\ \oplus I\right]
 S_{\rm mid}\mathcal H.
```

This is $`e^{iA}`$ for a full-space Hermitian A with $`\|A\|\le d`$.
Its compression $`Z=Vh`$ has

```math
Z=I+iK-M_2/2+O(d^3),\qquad
K=G(1/2)/2=J^\dagger AJ,\quad M_2=J^\dagger A^2J.
```

In particular $`0\le M_2-K^2\le d^2I`$. Comparing with $`U=e^{iK}`$
by the unitary Taylor remainder gives

```math
U^\dagger Z=T+E_3,\qquad
T=I-(M_2-K^2)/2\ge(1-d^2/2)I,\quad
\|E_3\|\le\tfrac56d^3.
```

This proves $`\|V-U\|\le3d^3`$ for $`d\le1/8`$. For completeness,
put $`X=T+E_3`$, $`h_X=|X|`$, $`t=1-d^2/2`$, and
$`r=\|E_3\|`$. The Sylvester identity

```math
h_X(h_X-T)+(h_X-T)T=X^\dagger X-T^2
```

gives $`\|h_X-T\|\le(2r+r^2)/(2t-r)`$, because $`\|T\|\le1`$.
Consequently $`\|\mathrm{polar}(X)-I\|`$ is at most
$`[r+(2r+r^2)/(2t-r)]/(t-r)\lt 3d^3`$. This uses no commutativity
assumption and is uniform in dimension.

## 6. Opposite filters and the complete cubic word

Let $`E_{\rm def}=I-h^2`$. Unitarity gives
$`E_{\rm def}=[(I-P)C_{1/2}J]^\dagger[(I-P)C_{1/2}J]`$, so
$`\|E_{\rm def}\|\le d^2`$. Apply (2) with both signs:

```math
D_\pm=e^{\pm i\pi(P-I/4)/3},\qquad
F_\pm=D_\pm^{-2}C_{1/2}D_\pm C_{1/2}^\dagger D_\pm C_{1/2}.
```

Their exact accepted blocks and rejected norms satisfy

```math
\begin{aligned}
J^\dagger F_\pm J
 &=Vh[1+e^{\mp i\pi/3}E_{\rm def}]\\
 &=V[I\mp i(\sqrt3/2)E_{\rm def}+R_\pm],
 &\|R_\pm\|&\le\|E_{\rm def}\|^2,\\
\|(I-P)F_\pm J\|&=\|E_{\rm def}\|^{3/2}\le d^3.
\end{aligned}
```

The final equality follows from
$`h^2[1+E_{\rm def}+E_{\rm def}^2]=I-E_{\rm def}^3`$.
The same physical half-word is used in both filters. In $`F_-F_+`$,
the order-two phase cancels; its remaining commutator is bounded by
$`\sqrt3\,\|E_{\rm def}\|\|V-I\|`$. Since
$`\|V-I\|\le d+3d^3`$, this commutator and the accepted remainders
cost less than $`2.2d^3`$. The two propagated rejected components cost
at most $`2d^3`$, and replacing $`V^2`$ by $`e^{iG(1/2)}`$ costs
at most $`6d^3`$. Hence

**Cubic correction (10).**

```math
\|F_-F_+J-Je^{iG(1/2)}\|\le11d^3.
```

The original physical coarse word cannot serve as an exact baseline for
this cubic estimate. Recompile the known rounded frame $`W_0`$ with
$`q_j=3s+n-j+10`$, using the ordinary full-operator bank construction.
Its error is below $`(86/1024)d^3`$. Compile the midpoint prefix at
$`q_j=2s+n-j+10`$ and its child bank at $`q=2s+10`$. Their combined
full-space error is below $`(172/1024)d^2`$. Equation (7) suppresses
each conjugator replacement by the half-middle size d.

For $`3s\le L\le N`$, choose $`q=L+20`$ for every fine middle and
phase axis, with three-axis error at most $`130\,2^{-q}`$. Use actual
inverses; the chosen $`\widehat D_-`$ may be the actual inverse of
$`\widehat D_+`$. The complete output word is the recompiled baseline
on the left of the two-filter correction. Its complete call ledger is:

| Physical subroutine | Appearances, including inverses |
|---|---:|
| Recompiled baseline | 1 |
| Half-correction $`C_{1/2}`$ | 6 |
| Midpoint S inside those corrections | 12 |
| Fine half-middle bank | 6 |
| Exact history extension | 12 |
| Phase bank | 8 |

The S, middle and history rows expand the half-correction row. Every
table, source, loader and actual inverse is charged a fixed number of
times, using the [native bank ledger](COARSE_PREFIX_ENCODER.md#3-the-canonical-su2-native-circuit).
Equations (9)–(10), all six suppressed conjugator replacements and the
baseline contribute less than

```math
\left(16+11+12\frac{172}{1024}+\frac{86}{1024}\right)d^3
 \lt 30d^3.
```

The fourteen fine-bank appearances add at most
$`14\cdot130\,2^{-(L+20)}\lt 2^{-L}`$. This proves the stated
$`40\,8^{-s}+2^{-L}`$ bound, with the same circuit's T and G counts.
Ordinary isometry telescoping is valid for arbitrary dirty inputs and
references; it requires neither a reset nor intermediate exact return.

## 7. Workspace and scope

For $`n\ge16`$, the sufficient peak dirty reservations are:

| Use | Quadratic word | Cubic word |
|---|---:|---:|
| Baseline | $`s+n+9`$ | $`3s+n+12`$ |
| Controlled auxiliary prefix | $`s+n+12`$ | $`2s+n+13`$ |
| Fine middle | $`L+22`$ | $`L+22`$ |
| Fine phase | $`L+21`$ | $`L+21`$ |

All are at most $`N+n+7`$. Each fine middle uses $`q+1`$ core wires
and a distinct arbitrary signal. Phase banks borrow a logical wire for
that signal. The
[self-borrowed queries and core-borrowed helper](../../docs/ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
need no additional selectors or helper, since $`q+1\ge2n`$. Coarse
prefix reservations include their selectors, predicate helper and,
when both supplied flags contain history, a separate arbitrary signal.
Source pools and helpers are reused at disjoint times, so reservations
are maximized, not summed.

Exact inactivity holds for each conditioned primitive on its unchanged
inactive sector. The complete filtered word includes approximate phase
banks and needs only the proved full-space or initialized-isometry
bound on other occupied sectors. The component theorem is restricted
to $`n\ge16`$; the retained finite-dimensional fallback for the original
endpoint target is a separate result.

The output discrepancy is generally not another product of rounded
correlated local coins: dirty operators and nonlocal commutators remain.
Consequently neither (8) nor (10) is an iterable refinement theorem.
Direct order-p source refinement costs $`O(N+pns)`$, which is
$`O(pN)`$ at $`s=\lceil N/n\rceil`$. An order comparable to n, or
repeated corrections with total $`O(N)`$ table and source cost, requires
a further all-order construction. This limitation gives no lower bound
for unrestricted native circuits.

## 8. A fixed-basis one-bank obstruction

The change to midpoint coordinates is substantive. At n=2, set all coarse
angles to zero, so each coarse coin is $`R_y(\pi/2)`$. In computational
coordinates the fixed-history Hermitian tangent image has traceless
diagonals and only the edges
$`\{0,2\},\{2,1\},\{0,3\}`$.

Allow arbitrary, even nonsmooth, retuning of the middle SU(2) blocks,
provided their distance from identity is $`O(\epsilon)`$. Uniform
expansion of their compression and its polar factor gives

```math
Z=I+iK-M/2+O(\epsilon^3),\qquad
\mathrm{polar}(Z)=I+iK-K^2/2+O(\epsilon^3),
```

where K and M are Hermitian, K has the displayed edge support, and
$`\|K\|=O(\epsilon)`$, $`\|M\|=O(\epsilon^2)`$. The (1,0)
entry of its anti-Hermitian part is therefore $`O(\epsilon^3)`$.
But choosing target root and left-child angles $`\epsilon`$ and the
right-child angle zero gives, by direct multiplication,

```math
\left[\frac{R-R^\dagger}{2}\right]_{1,0}
 =\frac{\epsilon^2}{2}+O(\epsilon^3),\qquad
R=W(0)^\dagger W(\epsilon).
```

Thus this fixed-basis one-bank polar ansatz has error at least
$`\epsilon^2/2-O(\epsilon^3)`$. An order-two computational diagonal
cannot remove the missing entry. A single sign of the filter contributes
an additional term $`\mp i(\sqrt3/2)(M-K^2)`$. Any first-order match
in this example has real-symmetric leading blocks; the resulting
covariance entry is real, so that filter contributes a purely imaginary
(1,0) term and cannot cancel the real witness.

Embedding the example in the first two levels and leaving later target
coins equal to their coarse coins gives the same witness for every
$`n\ge3`$: the common later layers cancel in R. The restriction concerns
this fixed conjugator and one small bank only. It does not cover changed
conjugators, multiple banks, or unrestricted native circuits.

## 9. Verification

The [collective precision tests](../../tests/test_collective_precision.py)
check the exact geometric history network, complete-frame residual
orientation and compression, the physical filter with dirty correlations,
midpoint cubic scaling, opposite-sign full-output cancellation, and the
fixed-basis witness. The
[coarse encoder tests](../../tests/test_coarse_prefix_encoder.py) cover
the imported delimiter layout and primitive resource interfaces. Finite
matrix checks diagnose the identities; the proofs above supply the
all-size norm, workspace, phase and reference guarantees.
