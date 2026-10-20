# Revision: complete bank motion is paid; a shallow source remains missing

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous continuation](PARALLEL_ACCESS_REVISION_RESULT.md)

This continuation starts from main
`6e42d0b2e2a52cb14962ae84ec4964933e3c213c`, tree
`e98fd8b3bebda6f9c02b394ad0bc407c01b3f2af`, after PR #205,
on 20 October 2026. It pays two complete bank-moving mechanisms,
identifies their depth cost, and keeps the original physical-pair scout
active. Neither construction changes the general resource frontier.

## 1. Destination, strongest theorems, components and scoped limits

| Role | Assessment after revision |
| --- | --- |
| General destination | Optimal unrestricted native T-depth, with simultaneous T-count, Clifford count, dirty width, complete error and return on the same circuit |
| Strongest general upper | [R96](ALL_PRECISION_LINEAR_DEPTH.md) remains the all-precision theorem; selected [R79](UNIFORM_SUBTREE_DEPTH.md) remains Result D |
| Useful reductions | [Root fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) absorbs polynomial source overhead; [clean promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) pays bounded-clean and amplitude-source wrappers |
| New constructive components | [Complete bank difference](COMPLETE_BANK_DIFFERENCE.md) pays nonlinear vector multiplication; [fermionic Walsh routing](BANK_MOTION_WALSH_INDICATOR.md) gives an exact complete indicator on every bank occupation sector |
| Limits of these approaches | The displayed Cartesian merge has an O(b squared) depth schedule. The Walsh construction has linear depth; the whole-word Majorana theorem excludes shallow resynthesis only inside its explicit grammar |
| Focused lower components | [All-active recognition](ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md) gives explicit sign/XOR bounds 2428/2427; [pair-only enlarged-center matching](PAIR_ONLY_ENLARGED_CENTER.md) also closes the all-active one-pair class with both silent lists coupling, when no middle singletons are present |

For $`N=2^n`$, $`n\ge1,L\ge6`$, R96 retains

```math
D_T=O(NL/B^2+L+\log_2(n+2)),\qquad
T=O(\sqrt{NL}+L+NL/B),\qquad G=O(NL).
```

There is one clean qubit, with real/complex dirty reservations respectively
$`B\ge17(L+n+7)`$ and $`B\ge17(L+n+8)`$, and complete
initialized-isometry error at most $`2^{-L}`$. The selected count
frontier and literal allocations are unchanged. At fixed accuracy and
square-root dirty width, the unrestricted logarithmic-versus-constant
depth gap remains open. The linear indicator T-count bound does not
itself give a growing depth bound.

## 2. A complete nonlinear bank difference and its paid recurrence

Put $`\omega=e^{i\pi/4}`$. For arbitrary r-bit banks Y,s, the
literal local interaction

```math
D_r=\prod_u T_{s_u}T_{Y_u\oplus s_u}^{\dagger}
   =\omega^{\sum_u(-Y_u+2Y_us_u)}
```

has 2r T gates, one native T layer and 2r CNOTs, with no helper.
Let complete queries F_f,F_g toggle s by Boolean vectors f(x),g(x)
and return their private work. Execute chronologically

```math
D_r,F_f,D_r^\dagger,F_g,D_r,F_f^\dagger,D_r^\dagger,F_g^\dagger.
```

The four phases have signal arguments s, s+f, s+f+g, s+g. Their
literal exponent is $`4\sum_uY_uf_ug_u`$ modulo eight, independently
of s. The signal toggles cancel. Output Hadamards therefore give the
complete query for $`f\odot g`$, returning the arbitrary signal bank.

For common forward/actual-inverse source bounds T,d,G,w and full-operator
errors epsilon_f,epsilon_g, the complete recurrence is

```math
\begin{aligned}
T_{f\odot g}&\le2T_f+2T_g+8r,&
d_{f\odot g}&\le2d_f+2d_g+4,\\
G_{f\odot g}&\le2G_f+2G_g+10r,&
w_{f\odot g}&\le r+\max\{w_f,w_g\},\\
\epsilon_{f\odot g}&\le2\epsilon_f+2\epsilon_g,&a_{\rm clean}&=0.
\end{aligned}
```

Private pools are disjoint from address, output and signal, and may be
reused sequentially. Telescoping includes source leakage and every return
error; an approximate intermediate call is not treated as exact cleanup.
All inputs, correlations, references and actual emitted adjoints are covered.

For b=a+c address bits, let $`m=2^b`$. Actual row/column CNOT
encodings spread child equality toggles across arbitrary signal rows and
columns. They supply the complete Cartesian recurrence

```math
\begin{aligned}
T_b&\le2T_a+2T_c+8\,2^b,&d_b&\le2d_a+2d_c+4,\\
G_b&\le2G_a+2G_c+18\,2^b-4(2^a+2^c),\\
w_b&\le2^b+\max\{w_a,w_c\},&
\epsilon_b&\le2\epsilon_a+2\epsilon_c.
\end{aligned}
```

The entire signal bank is charged; child outputs are its existing pivots.
With the exact one-bit Clifford base and balanced splitting, safe bounds
are $`T_b\le32\,2^b`$, $`G_b\le72\,2^b`$,
$`w_b\le2\,2^b`$ and zero error. The displayed schedule is O(b
squared), with $`d_b\le4(b^2-1)/3`$ for powers of two. This is a
paid upper schedule, not a lower bound on bank motion or resynthesis.
Using the retained logarithmic-depth children changes only constants.
The merge pays cleanup and nonlinear action but retains the old depth cost.

## 3. Walsh bank motion: complete identity, native ledger and return

Identify the m-qubit bank with all occupation sectors of m fermionic
modes. Let $`\Gamma(V)=\bigoplus_{r=0}^m\bigwedge^rV`$,
$`W=H^{\otimes b}`$ on the mode-index space, and $`F=\Gamma(W)`$.
This F is not b physical Hadamards. The actual Clifford

```math
C=\prod_{i,u:u_i=1}\mathrm{CZ}_{x_i,Y_u}
```

has branch $`\Gamma(\mathrm{diag}_u((-1)^{x\cdot u}))`$.
Thus $`R=F^\dagger CF`$ is controlled fermionic translation and
$`R^\dagger Z_0R=\sum_x|x\rangle\langle x|\otimes Z_x`$.
Conjugating by all output Hadamards gives the complete indicator I_b.
The fermionic permutation signs cancel through the actual inverse on
every occupation sector, including vacuum and fully occupied inputs.

Each butterfly matching has m independent commuting Majorana-bilinear
T axes and one native T layer. Its literal phase and actual inverse are
paid by fixed fermionic permutations and elementary Clifford encoders.
The resulting exact transform has $`T_F=mb`$, $`d_F\le b`$,
$`G_F\le6bm^2`$ and no additional work. More generally, if
$`\|\widetilde F-F\otimes I_w\|\le\eta`$, its four occurrences
in the complete echo give

```math
\begin{aligned}
T_E&\le4T_F,&G_E&\le4G_F+3bm+2m+2,\\
d_E&\le4d_F,&w_E&\le w_F,\\
\|\widetilde I_b-I_b\otimes I_w\|&\le4\eta,&a_{\rm clean}&=0.
\end{aligned}
```

The same ledger covers the actual inverse. Its error includes all private
return; exact involutivity of the emitted echo does not imply exact work
return. The exact displayed construction has

```math
T_E\le4bm,\qquad d_E\le4b,\qquad
G_E\le24bm^2+3bm+2m+2,\qquad w_E=0,\qquad\epsilon_E=0.
```

This fits the relaxed count budgets but loses to the retained logarithmic
source depth. A one-particle calculation alone would not establish this
source; the full Fock-space operator identity above does.

## 4. Whole-word Majorana limit and the next unpaid operation

The [whole-word audit](WHOLE_WORD_MAJORANA_LIMITS.md) strengthens the
direct-transform support test to the completed phase-query word. Its
grammar requires every actual Clifford boundary to preserve the read-only
address and the span of the standard Majoranas of **all** work wires,
including arbitrary dirty spectators. Native T layers and all boundary
actions must satisfy that actual representation, not just a chosen
single-particle description of the target.

For b at least two, zero clean work, complete phase-query error
$`\epsilon\lt1`$ allowing one common scalar, and d native T layers
in this grammar imply

```math
2^d\ge2^b(1-\epsilon^2/2)^2.
```

Consequently $`d\ge b`$ for $`\epsilon\lt\sqrt{2-\sqrt2}`$,
and $`d\ge b-1`$ for every $`\epsilon\lt1`$. This includes
whole-word resynthesis inside the stated grammar and prevents a qualifying
root-fusion source there. It does not constrain arbitrary Clifford
interference, address motion or unrestricted native depth. The displayed
Cartesian recurrence has its separate upper-schedule interpretation.
For an XOR indicator, the actual output-Hadamard-conjugated phase word
must satisfy this grammar at both endpoints too; this is not automatic.

The first unpaid operation remains a complete source with
$`T,w=O_\delta(2^b(b+2)^c)`$,
$`G=O_\delta(4^b(b+2)^{c_G})`$ and
$`d_\delta(b)=o(\log(b+2))`$ at every fixed requested delta.
Root fusion uses eight source/actual-inverse calls, error at most eight
times source error, and an exact $`O_\tau(\log\log(p+4))`$-depth
suffix on $`b=p-O_\tau(\log(p+2))`$ bits. For the Walsh wrapper
its source contribution is therefore at most $`32\eta`$.

The next research question must change the unpaid dynamics: use
non-Majorana Clifford interference, intermediate address motion, or a
different complete source interface. A candidate must first exhibit its
full selected action and return, then supply the native schedule. Another
representation of the retained linear or logarithmic cost does not improve
the frontier.

The original overlapping amplitude-$`1/\sqrt2`$ route remains an
alternative with 24 source/actual-inverse calls, depth at most
$`24d_U+32`$ and error at most
$`16(3\sqrt2\epsilon+2\epsilon^2)`$. Signal/toggle return is
approximate unless exact flatness is separately proved. The bounded-clean
two-band amplitude-free defect still stops every-fixed-relative-accuracy
sources, uniformly over positive amplitudes; it is neither an arbitrary
absolute-error floor nor an unrestricted depth bound.

## 5. Original physical-pair scout and verification boundary

The earlier explicit physical-list criterion has one genuine original middle
pair, whose difference is the sum of all actual active axes at both
endpoints, and an input-silent list commuting with the original middle
list, or its time reverse. The complete sign/XOR degrees are at most
2428/2427, with original-circuit rounding floor $`2^{-530690}`$.
Both active ranks, silent spans, middle singletons and dirty width may grow.
The proof retains the unused clean block and literal phase throughout.

That proof uses either a physical clean flip or the enlarged subgroup H_0.
Its 603 coefficient features, one parity bit and three diagonal-C3 digits
give the bound; the commuting union is used only for C3 membership.
Qualitatively the class falls under the earlier reversed t=0 theorem.
The advance is its explicit recognition and sharper bound.

The new [pair-only closure](PAIR_ONLY_ENLARGED_CENTER.md) removes the
one-sided silent-commutation condition when the original middle list is
exactly that genuine pair, with no singleton axes. The pair difference
must still equal the sum of all actual active axes at both endpoints.
No bound on either silent span is assumed. Its complete bounds are

```math
D_*=576\cdot2^{62}+2252,\qquad
\deg g\le D_*,\qquad\deg f\le D_*-1.
```

Products in the kernels of the silent pair couplings commute the entire
middle list and acquire actual endpoint transport. When both lists couple,
their cross-pairing matrix has rank at most one. The zero case supplies a
physical clean flip; the other case enlarges the transported center to Z with
$`\dim(M/Z)=4`$ for the entire physical axis span M. This is stronger
than a pairing-rank assertion and does not bound the old additional span t.
The localization Clifford fixes Z and preserves M. One global alignment
and the feature ledger retaining the physical clean bit then prove $`D_*`$.
The same original-circuit floor $`2^{-530690}`$ applies. No native
band or unused-sector extension is supplied for free.

The next lower questions are original middle singletons, which can block
this transport, and general pair differences. Unrestricted depth remains open.

The two focused fixtures are the [Walsh indicator](../../tests/test_bank_motion_walsh_indicator.py), checking
all 64 six-wire columns, all occupation sectors and the actual inverse,
and the [all-active pair fixture](../../tests/test_all_active_original_pair_closure.py),
checking the physical span dichotomy, actual clean flip, full columns,
phases and inverse. The Walsh emitted receipt is T=32, depth at most eight,
G=258, with no clean or additional dirty work. These are mechanism witnesses;
the bank difference, Majorana limit and both degree results have analytic proofs.
Required full validation and exact-tree publication receipts belong on the
pull request; no full-suite result is asserted by this checkpoint.
