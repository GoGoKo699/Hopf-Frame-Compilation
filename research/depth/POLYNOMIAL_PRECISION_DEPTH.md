# Depth at stretched-exponential accuracy

[Precision-bridge design](BORROWED_PRECISION_BRIDGE_PROMPT.md) ·
[Global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md) ·
[Fast geometric references](FAST_CAPPED_GEOMETRIC_REFERENCE.md) ·
[Logarithmic-precision theorem R95](LOGARITHMIC_PRECISION_DEPTH.md) ·
[All-precision theorem R96](ALL_PRECISION_LINEAR_DEPTH.md)

This independent supporting theorem supplies the lower-precision branch
of [R96](ALL_PRECISION_LINEAR_DEPTH.md). The global early stage places its
precision core in arbitrary dirty work. Its fast
reference words keep the initialized reserve logarithmic in requested
precision. Three such stages leave a logarithmic-height tail, on which
R79's high-precision clause and R94's complete-workspace catalogue give
the same depth law throughout a much larger precision range.

## 1. Statement and scope

Let $`N=2^n`$, $`\ell=\log_2(n+2)`$ and $`\eta=2^{-L}`$, with
integers $`n\ge1`$ and $`6\le L\le n^{3/4}`$. Every prescribed
complete real Hopf frame admits one native Clifford+T circuit V, with one
external initialized qubit and any arbitrary dirty allocation
$`b\ge17(L+n+7)`$, satisfying

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+L+\ell\right).
\qquad(1)
```

The prescribed complex magnitude completion has the same bounds at
$`b\ge17(L+n+8)`$. All constants are absolute on this fixed range.
The full initialized-isometry norm includes every logical column, literal
signs and phases, singular continuations, all work leakage, arbitrary dirty
inputs and external references. All inverses are actual native inverses.

At $`L=\Theta(n^{3/4})`$ and $`b=\Theta(\sqrt{NL})`$, (1) gives
$`O(n^{3/4})`$ depth at error $`2^{-\Theta(n^{3/4})}`$.
It includes the designed target through $`L\le n^{1/4}`$.
The matching worst-case T-count bounds are inherited from R79; no matching
depth lower bound is asserted. Classical preprocessing terminates but need
not be efficient. No measurement, reset, QRAM, supplied precision state,
additional external clean word or free coherent evaluator is used.
[R96](ALL_PRECISION_LINEAR_DEPTH.md) extends this upper depth law to every
precision through a different terminal. Optimal general depth remains open.

R95 supplies (1) unchanged when $`L\le2\ell`$. Below we prove the real
construction on the overlapping internal range

```math
\ell\le L\le2n^{3/4}.                                   \qquad(2)
```

The extra factor two permits the complex precision shift. This is not a
claim about arbitrary complex-unitary synthesis, efficient preprocessing,
a practical crossover, or a scalable native emitter.

## 2. The uniform fast early-stage interface

Use the [global borrowed-precision stage](GLOBAL_BORROWED_PRECISION_STAGE.md)
with its [fast reference substitution](FAST_CAPPED_GEOMETRIC_REFERENCE.md).
An immutable global prefix P has $`n-m`$ bits and the remaining frame
has height m. At stage precision $`l\ge6`$, put

```math
K=m+l+2,\qquad
\chi=\lceil\log_2(l+2)\rceil+\lceil\log_2(m+2)\rceil^3,
\qquad r=\lceil A\chi\rceil.                             \qquad(3)
```

When $`r\lt m/2`$, the stage compiles the early factor

```math
E=\sum_P|P\rangle\langle P|\otimes
 \left(W_{P,m-r}\otimes\Pi_Z+I\otimes(I-\Pi_Z)\right),
\qquad \Pi_Z=|0^r\rangle\langle0^r|,
```

with full initialized-isometry error at most $`2^{-l}/16`$.
It saves $`H=[Z=0^r]`$ into the sole external initialized bit before
using Z as conditional work, is literal identity on H zero on all
occupied inputs, and applies the actual saved-predicate inverse last.
All work leakage is included in this error. The mathematical coarse
frame and exact matching still use the original positive reference
columns; only their native preparation is approximated, with a proved
full-space comparison for transported anchors.

For fixed absolute C,d and $`d_w`$, let

```math
\overline Q=2^{n-r+C\chi},\qquad
S_{\rm loc}=K^C2^{C\chi},\qquad
b\ge8W_{\rm loc}+CK^{d_w}(n+\chi+2),
\qquad W_{\rm loc}\le S_{\rm loc}.
```

The sharper component ledger gives

```math
\begin{aligned}
T_{\rm stage}&\le CK^d(\sqrt{\overline Q}+\overline Q/b+1)
                    +S_{\rm loc},\\
G_{\rm stage}&\le CK^d(\overline Q+1)+S_{\rm loc},\\
D_{T,\rm stage}&\le C\left[l+\log_2(m+2)+\log_2(n+\chi+2)
                                  +K^d\overline Q/b^2\right].
\end{aligned}                                               \qquad(4)
```

These are separate simultaneous bounds on the same physical word.
Only constantly many table positions are sequential, with at most a
fixed polynomial in K simultaneous instances and outputs. The independent
component proof prices every large table, including P-dependent hashes
and matchings, source anchors, actual inverses, scalar masks, and tables
that happen not to depend on P. A table reading early logical coordinates
belongs to the $`n-r+C\chi`$ envelope. Remaining standalone tables
read only $`O(\chi)`$ metadata, and local layouts, source cores and
predicates have polynomial native resources. This distinction is essential
at m=n: the weaker bound $`2^{C(m+\chi)}`$ would not suffice.

All degrees, error constants, and simultaneous bank coefficients are
fixed before A. The all-work anchor tests the actually used conditional
fields, excluding unused protected padding and every arbitrary dirty
precision core. The original r-bit saved predicate is a separate local
conjunction. Thus the address coefficient C is independent of A; no
full r-bit bank is appended to the large global tables.

Since $`m+l+2\le(m+2)(l+2)`$, we have $`\log_2K\le\chi`$.
Choose A large enough that

```math
r-C\chi\ge(2d+80)\chi,
```

```math
K^d2^{-(r-C\chi)/2}\le K^{-40},\qquad
K^d2^{-(r-C\chi)}\le K^{-80}.                            \qquad(5)
```

Therefore the large-table terms in each stage fit
$`O(\sqrt N+N/b)`$ T-count, $`O(N)`$ Clifford count and
$`O(N/b^2)`$ T-depth. The native source and reference schedule in (4)
supplies the additive depth; it is not inferred from a gate-count bound.

## 3. Three contractions from the original global height

Work first above one absolute dimension cutoff on (2), with
$`b\ge\sqrt N/n`$. All cutoffs depend only on the fixed native
constants, never on L or the prescribed angle data. Make the retained
whole-family finite magnitude/sign reduction at error $`\eta/512`$.
Its two exact Boolean sign phases read the full n-bit logical address
and cost $`O(\sqrt N+N/b)`$ T gates, $`O(N)`$ Clifford gates and
$`O(N/b^2+\ell)`$ T-depth. Singular and originally zero branches
retain their prescribed complete-frame continuations. All following
stages refer to that positive tuple.

Start at $`m_0=n`$. Apply at most three stages (3)–(4), each at

```math
l=L+5.                                                   \qquad(6)
```

At every step compute r for its current height m. If $`r\lt m/2`$,
execute the stage and retain r original bottom bits. Otherwise omit it
and keep the current height, which obeys $`m\le2r`$. A height at most
three goes directly to the final bottom compiler. Thus every retained
next height, including an omitted stage, satisfies

```math
m_{j+1}\le C\left[\log_2(L+2)+\log_2^3(m_j+2)\right].     \qquad(7)
```

Put $`a=\log_2(L+2)`$ and $`t=\log_2(\ell+2)`$. The upper bound
in (2) gives $`a=O(\ell)`$, so the three successive bounds are

```math
m_1\le C(a+\ell^3)=O(\ell^3),
```

```math
m_2\le C(a+t^3),
```

```math
m_3\le C\left[a+\log_2^3(a+t^3+2)\right]=O(a).          \qquad(8)
```

For the final implication, $`L\ge\ell`$ gives $`t\le a`$, hence
$`m_2\le C'(a+1)^3`$. The absolute inequality
$`\log_2^3(a+2)\le C''a`$ for $`a\ge\log_2 8`$ then proves
$`m_3=O(a)`$. Its bounded range and large-a range are covered by one
constant. No fixed-precision stopping threshold is used.

The two essential properties of the component are visible here.
Its reserve contains only one $`\log_2(l+2)`$ term, even when
$`l\gg m`$, so the last contraction reaches $`O(\log L)`$.
Its fast references have depth $`O(l+\log(m+2))`$, so the first
stage at m=n costs $`O(L+\ell)`$, rather than $`O(n)`$.
A failed reserve test never requests a nonproper protected suffix.

## 4. Opening resources and simultaneous private pools

On (2), $`K=m+l+2=O(n)`$ at every stage and
$`\chi=O(\ell^3)`$. Consequently

```math
\log_2 S_{\rm loc}=O(\log_2K+\chi)=O(\ell^3)=o(n).       \qquad(9)
```

This applies even to the first stage at m=n. The ordinary sources have
only logarithmic initialized indices; their $`O(l)`$ precision cores
are arbitrary dirty systems with polynomial native count and width.
Their exponential-in-precision initialization costs have not been hidden
in a polynomial prefactor. The fast reference word uses a geometrical
cutoff of $`O(l+\log(p+1))`$ atoms, with p the number of bins,
and has an explicit native count and depth at that cutoff. The original
reference geometry and all its high coin bits remain in the exact
mathematical embedding and augmented anchor tests.

Thus the standalone counts and dirty work in (9) are $`N^{o(1)}`$.
After one absolute cutoff, uniformly on $`b\ge\sqrt N/n`$, choose
$`W_{\rm loc}\le b/16`$ and
$`CK^{d_w}(n+\chi+2)\le b/2`$. These two bounds imply the
actual component requirement $`8W_{\rm loc}+CK^{d_w}(n+\chi+2)\le b`$.
All private-query minimum widths and instance counts are fixed polynomials
in n, so their simultaneous pools satisfy that requirement. Outputs,
precision cores, donor coordinates,
source banks, logical data, addresses and query helpers are disjoint.
The component's separate full-r saved-predicate conjunctions also fit.

All three stages have at most $`O(L+\ell)`$ additive depth by (4).
Use (5) for their large-table terms and (9) for their standalone terms.
Together with the sign phases, the complete opening obeys

```math
T_{\rm open}=O(\sqrt N+N/b),\qquad G_{\rm open}=O(N),\qquad
D_{\rm open}=O(N/b^2+L+\ell).                           \qquad(10)
```

The N-subpolynomial terms are eventually smaller than $`\sqrt N`$;
there are only three sequential stages. Their actual completed words
reuse returned allocations. Their initialized-isometry errors propagate
through later unitaries, so approximate return is not used as a reset.
Each stage contributes at most $`\eta/512`$ by (6).

## 5. A constant number of final-accuracy high-precision groups

Write $`M=m_3`$. If $`M\le3`$, omit the catalogue and use the bottom
in Section 7 directly. Otherwise retain R94's original three bottom
logical bits $`Z=(c,d,G)`$ and split the remaining $`M-3`$ layers
into groups of height at most

```math
g=\max\{1,\lfloor\tfrac14\log_2(L+2)\rfloor\}.          \qquad(11)
```

Except for a possibly short first group, every group has height g.
Equation (8) and $`g=\Omega(\log_2(L+2))`$ give a fixed absolute
bound $`R\le R_*`$ on the number of groups. Fix this bound before
choosing their precision padding:

```math
\kappa=L+\lceil\log_2(R+1)\rceil+10=O(L).               \qquad(12)
```

Use R94's actual final-accuracy angle code, whose h-bit group obeys

```math
q_h=\sum_{j=0}^{h-1}2^j(\kappa+h-j+6)
       \le(\kappa+8)2^h,\qquad
\|W-W_v\|\lt2^{-\kappa}/16.                            \qquad(13)
```

Every code is a legal positive frame, and code zero has the actual empty
physical word. Compile the controlled frame
$`A_v=\mathrm{diag}(I,W_v)_c`$ on h+1 bits by R79. R94's zero-root
construction proves that this is a complete Hopf frame; it does not
generically control a synthesized circuit. The high-precision condition
holds uniformly, since

```math
\frac{2^{h+1}}{(h+1)^2}\le2^{h+1}
 \le\max\{4,2(L+2)^{1/4}\}\lt L\le\kappa
 \qquad(L\ge6).                                       \qquad(14)
```

Choose the full local dirty pool as in R94:

```math
b_{\rm loc}=\left\lceil\max\left\{17(\kappa+h+8),
                    C_{\rm loc}\sqrt{2^{h+1}\kappa}\right\}\right\rceil.
```

Because $`2^{h+1}=O(L^{1/4})`$, R79 gives the separate explicit bounds

```math
b_{\rm loc}=O(L),\quad t_{\rm loc}=O(L),\quad
 g_{\rm loc}=O(L^{5/4}),\quad d_{\rm loc}=O(L).           \qquad(15)
```

The local depth-width term $`2^{h+1}\kappa/b_{\rm loc}^2`$ is
bounded. No unspecified polynomial exponent in an exponential local
compiler allowance is used. Every physical local wire belongs to the
payload $`S=(c,X,d,D_*)`$, of width $`w=h+2+b_{\rm loc}=O(L)`$.

R94's complete-payload transfer and commutator retain their exact
full-space identities. All undesired occupied-slot operations cancel,
including entanglement with their local flags and work. The selected
physical word and its actual inverse give group error below
$`3\cdot2^{-\kappa}`$ on c one and d zero, including the rounded
code error. All catalogue slot words act on disjoint complete payloads,
so their depths take a maximum and their counts add.

Retain R94's exact saved H, original G guards, zero-code inactive word,
paired auxiliary router array and changing-prefix schedule. Each group's
prefix is immutable during that group and queried afresh for the next.
H zero remains literal identity on occupied work. Full initialized-isometry
telescoping gives total group error at most
$`3R2^{-\kappa}\le3\eta/1024`$. It requires neither correctness of a
local compiler on occupied d inputs nor dirty reinitialization.

## 6. Global catalogue resources on the same circuit

Let $`V_h=2^{q_h}`$ and let Q be a group's actual global-address row
count. With an exclusive $`B_Q=\lfloor b/2\rfloor`$ pool,
[R94's full native ledger](PORTABLE_WORKSPACE_CATALOGUE.md#4-native-resources-and-the-final-accuracy-code)
and the R91 complete paired router give

```math
\begin{aligned}
T_i&=O(\sqrt{Qq_h}+Qq_h/B_Q+V_ht_{\rm loc}+wV_h^2),\\
G_i&=O(Qq_h+V_hg_{\rm loc}+wV_h^2),\\
D_i&=O(\ell+Qq_h/B_Q^2+d_{\rm loc}+\log_2(q_h+2)).
\end{aligned}                                               \qquad(16)
```

All four catalogue appearances, four complete selected swaps, actual
query inverses, both arrays and returned helpers are included. From
(2), (11), (13) and (15),

```math
q_h=O(L^{5/4})=O(n^{15/16})=o(n),
```

```math
\log_2(wV_h^2+V_ht_{\rm loc}+V_hg_{\rm loc}
                   +b_{\rm loc}+q_h+2)
 =O(L^{5/4}+\log_2(L+2))=o(n).                          \qquad(17)
```

Thus every catalogue and router grid, including the complete
$`w2^{2q_h}`$ bill, costs $`N^{o(1)}`$ in gates and simultaneous
dirty width. After one absolute cutoff, these arrays, the logical local
pool, program and predicate helpers fit inside b/4, leaving the exclusive
query pool and spare fixed fractions. The sufficient condition
$`B_Q\ge C(n+2)`$ holds uniformly at $`b\ge\sqrt N/n`$.
Every slot contains all its own work; no slot borrows another slot's pool.

For group i let $`s_i=3+(R-i)g`$ be its lower suffix size and $`h_i`$
its actual height. Its prefix has $`n-s_i-h_i`$ bits; H,G add only a
factor four. Consequently

```math
Q_i\le4\,2^{n-s_i-h_i},\qquad
Q_iq_i\le C N\kappa2^{-s_i},
```

```math
\sum_iQ_iq_i=O(NL),\qquad
\sum_i\sqrt{Q_iq_i}=O(\sqrt{NL}).                       \qquad(18)
```

The short first group's height cancels in the same product. These are
geometric sums with ratios at most 1/2 and $`1/\sqrt2`$; the full
original prefix has not been replaced by the shorter program address.
Equations (15)–(18), the constant number of groups and logarithmic-depth
guard predicates give

```math
T_{\rm mid}=O(\sqrt{NL}+NL/b),\quad G_{\rm mid}=O(NL),\quad
D_{\rm mid}=O(NL/b^2+L+\ell).                           \qquad(19)
```

The opening and catalogue stages reuse physical allocations sequentially.
Their simultaneous dirty peaks, not their lifetime sums, are bounded by
b. Every private query pool is disjoint from all live logical coordinates,
source banks, addresses, outputs and other helpers. This proves the
same-circuit count, width and depth bounds; separate compilers with
incompatible resource minima have not been spliced together.

## 7. Bottom and complete real error

The deepest $`\min(M,3)`$ layers use
[R94's retained R72 smaller-error bottom](PORTABLE_WORKSPACE_CATALOGUE.md#8-bottom-total-error-original-widths-and-complex-completion),
at $`q_j=L+j+7`$ and $`m_j=L+j+8`$ for the relevant
$`j\in\{1,2,3\}`$. Its base $`L+n+9`$ and a separate
$`\Theta(b)`$ pool fit the wide branch. Its global-prefix query sizes
are $`Q_j=N2^{-j}`$. The one-clean scalar sources have depth
$`O(m_j)`$; the sharp queries add $`O(\ell+Q_jm_j/b^2)`$.
With at most three terms this gives

```math
T_{\rm bot}=O(\sqrt{NL}+L+NL/b),\quad G_{\rm bot}=O(NL),\quad
D_{\rm bot}=O(NL/b^2+L+\ell).                           \qquad(20)
```

This retains the L term that was absorbed into $`O(\ell)`$ on R94's
smaller precision range. The retained weighted error sum is below
$`15\eta/64`$ and includes actual work return.

Add the sign reduction, three new stages, catalogue groups and bottom.
The wide real error is at most

```math
\frac{\eta}{512}+\frac{3\eta}{512}
 +\frac{3\eta}{1024}+\frac{15\eta}{64}
 =\frac{251\eta}{1024}\lt\frac{\eta}{4}.                \qquad(21)
```

Omitted stages contribute zero. The unitary hybrid compares each stage
on its ideal initialized domain and propagates earlier leakage through
the actual subsequent words. No exact intermediate return or reset is
assumed. Equations (10), (19) and (20) prove (1) on the wide internal
real branch (2).

## 8. Literal thresholds, finite dimensions and complex phases

For $`b\lt\sqrt N/n`$, use the whole retained R79 compiler at the
original precision and literal width. Its additive depth obeys

```math
F_{79}(n,L)\le nL\lt NL/b^2,
```

because $`NL/b^2>n^2L\ge nL`$. Its count bounds already equal (1).
For n below the common absolute cutoff use that same original-precision
compiler at every width. Here $`nL\le n_0L=O(L+\ell)`$, with one
absolute constant independent of L and the frame. No stronger-precision
fallback is demanded at the literal minimum width.

For real inputs at $`L\le2\ell`$, use R95 unchanged. Together with
the overlapping internal range (2), this covers all claimed real
parameters. For complex inputs in R95's range use its established
complex theorem directly. On the remaining large wide complex branch,
$`L\gt2\ell`$ and $`L\le n^{3/4}`$, apply (21) at L+1 and append
the retained one-clean literal diagonal at L+1. Since $`L\ge6`$,

```math
\ell\le L+1\le\tfrac76L\le2n^{3/4},\qquad
17((L+1)+n+7)=17(L+n+8).                               \qquad(22)
```

The real call therefore has exactly the original complex width threshold.
Its error is below $`\eta/8`$, and the diagonal error is at most
$`\eta/2`$. The diagonal has the same threshold and costs
$`O(\sqrt{NL}+L+NL/b)`$ T gates, $`O(NL)`$ Clifford gates and
$`O(NL/b^2+L+\ell)`$ T-depth. Specified phases at zero-amplitude
leaves are retained. The real call and diagonal reuse the one external
initialized bit, with leakage propagated by the actual second circuit.

Narrow and bounded-dimension complex branches use the whole original-
precision R79 complex compiler, without appending another diagonal.
This establishes the literal real and complex thresholds and all resource
bounds together on the stated precision range.

## 9. Dependencies and limits

The global stage proves its native source, address, reserve, standalone,
error and return contracts independently. The reference component uses
an explicit vector comparison and an analysis-only full-space completion;
no unimplemented rotation is charged as a physical gate. The final tail
uses R79's already established high-precision clause and R94's algebraic
complete-payload transfer. Neither component invokes this theorem, so the
proof graph is acyclic.

The original designed quarter-power route first used R89's cubic reserve
and R95's shared-gradient opening before two small-workspace contractions.
The fast reference substitution permits the three-stage sequence (8)
directly from n, and its sharper standalone ledger removes both earlier
openings. The quarter-power target is contained in the present theorem;
there is no separate competing proof path needed for its statement.

This complete-catalogue construction does not reach $`L=\Theta(n)`$. At its chosen group
scale, the actual final code has $`q_h=\Theta(L^{5/4})`$, so a complete
$`2^{q_h}`$ catalogue no longer fits the available global budget.
Constant-height groups avoid that program size but restore an unbounded
number of sequential precision charges. This is a limitation of this
catalogue assembly, not an unrestricted circuit lower bound. The
[direct prefix-programmed terminal](DIRECT_PREFIX_HIGH_PRECISION_TERMINAL.md)
bypasses that catalogue and supplies R96's all-precision extension.
Optimal general depth,
efficient preprocessing, practical crossover and scalable native emission
are not established here.
