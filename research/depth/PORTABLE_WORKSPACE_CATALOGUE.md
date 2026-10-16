# Portable workspace catalogues at growing precision

[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[Executed design](CATALOGUE_BOOTSTRAP_PROMPT.md) ·
[Previous grouped construction](GROUPED_CATALOGUE_PRECISION.md) ·
[Retained local compiler R79](UNIFORM_SUBTREE_DEPTH.md)

**Local result R94.** The construction carries a local compiler's
complete workspace through a dirty catalogue, compiles
each group directly to final accuracy, and uses a constant number of
groups. The full payload is essential to its exact array cancellation.

## 1. The theorem and contract

Let $`N=2^n`$, $`\ell=\log_2(n+2)`$, and $`\eta=2^{-L}`$, with
integers $`n\ge1`$ and $`6\le L\le2\ell`$. For every prescribed
complete real Hopf frame, one external initialized qubit and every
arbitrary dirty allocation $`b\ge17(L+n+7)`$ admit one native circuit
with the simultaneous bounds

```math
\|\mathcal VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad G=O(NL),
\qquad
D_T=O\!\left(\frac{NL}{b^2}+\ell\log_2(\ell+2)\right).
```

The prescribed complex magnitude completion has the same bounds at
$`b\ge17(L+n+8)`$. All constants are absolute on this range.
The inherited [R79 lower bounds](UNIFORM_SUBTREE_DEPTH.md#7-matching-intervals-evidence-and-remaining-limits)
make this T-count worst-case optimal at every stated width; no matching
large-width depth lower bound is asserted.
The proof establishes the real construction internally through
$`L\le3\ell`$ to pay the complex precision shift. At
$`L=\Theta(\ell)`$ and $`b=\Theta(\sqrt{NL})`$, its depth is
$`O(\ell\log_2(\ell+2))`$.

The norm includes every logical column, literal phase, work leakage,
dirty input and external reference. Every inverse is the actual native
inverse word. There is no measurement, reset, QRAM, supplied precision
state, additional externally initialized word, or free coherent
evaluator. Elementary Clifford gates are counted. Classical preprocessing
terminates; no efficiency or scalable emitted-circuit claim is made.
Optimal general depth and the unrestricted precision frontier remain open.

## 2. A controlled Hopf frame and its physical local compiler

For an h-bit real Hopf frame $`W_v`$, with $`h\ge1`$, introduce a
most-significant selector c and prescribe

```math
A_v=|0\rangle\langle0|_c\otimes I_X+
       |1\rangle\langle1|_c\otimes W_v,\qquad B=X_c.
```

This is an (h+1)-bit Hopf frame: set its root angle to zero, every
left-subtree angle to zero, and the right-subtree tuple to that of
$`W_v`$. The root addressed layer is identity. Every later layer
is identity on c zero and the corresponding layer of $`W_v`$ on
c one. Thus the complete marker-frame operator is exactly $`A_v`$,
including the prescribed continuation at the zero root angle.
This uses the [complete addressed-layer interface](../../docs/HOPF_INTERFACE.md#4-complete-addressed-depth-operator);
it does not control an already synthesized native circuit.

Apply the retained real R79 theorem to $`A_v`$ at precision
$`\epsilon=2^{-\kappa}`$. Its initialized flag is d and its local
arbitrary dirty pool is $`D_*`$. The resulting actual physical
unitary $`U_v`$ acts on the entire payload

```math
S=(c,X,d,D_*),\qquad w=h+2+b_{\rm loc},
\qquad
\|U_vJ_d-J_d(A_v\otimes I_{D_*})\|\le\epsilon.
\qquad(1)
```

No correctness claim is needed on occupied d. Write
$`\bar A_v=A_v\otimes I_{D_*}`$ and
$`E_v=U_vJ_d-J_d\bar A_v`$. Actual inversion gives

```math
U_v^\dagger J_d-J_d\bar A_v^\dagger
       =-U_v^\dagger E_v\bar A_v^\dagger.                 \qquad(2)
```

Hence the inverse initialized-isometry error is also at most epsilon.
The physical B acts only on c, so it intertwines with $`J_d`$.
The two-term hybrid is

```math
\begin{aligned}
&U_vBU_v^\dagger B^\dagger J_d
       -J_d\bar A_vB\bar A_v^\dagger B^\dagger\\
&=U_vB(U_v^\dagger J_d-J_d\bar A_v^\dagger)B^\dagger
 +(U_vJ_d-J_d\bar A_v)B\bar A_v^\dagger B^\dagger .
\end{aligned}                                           \qquad(3)
```

Its norm is at most $`2\epsilon`$. Direct block multiplication gives

```math
A_vBA_v^\dagger B^\dagger
       =\mathrm{diag}(W_v^\dagger,W_v)_c.                \qquad(4)
```

Thus c one gives the forward frame. The ideal comparison returns c,d
and $`D_*`$; actual leakage into any of them remains in (3).

## 3. Exact transfer of the complete workspace

Let an immutable external address P select $`v=f(P)`$ from
$`V=2^q`$ fixed physical words. Every word has the same w-bit
payload layout. Allocate one arbitrary catalogue array
$`\mathcal A`$, consisting of V complete payload slots, and define

```math
F=\bigotimes_t U_t(\mathcal A_t),\qquad
S_v=\mathrm{SWAP}(S,\mathcal A_v),\qquad
K_v=S_vFS_vF^\dagger.
```

**Rightmost factors act first** in operator displays. Conjugation by the
full payload swap gives the literal full-space identity

```math
K_v=U_v(S)\otimes U_v^\dagger(\mathcal A_v)
                         \otimes I_{\rm other\ slots}.  \qquad(5)
```

Since B acts only on S, its complete commutator is

```math
K_vBK_v^\dagger B^\dagger
       =(U_vBU_v^\dagger B^\dagger)(S)\otimes I_{\mathcal A}.
\qquad(6)
```

No array flag is assumed initialized. The equation includes arbitrary
occupied slots, reference entanglement and literal phase. Direct sums
over preserved P prove coherent selection. One catalogue array suffices
because B is local; this does not merge R90's two independently
transferred factors onto one array.

The chronological native word for (6) is

```math
B^\dagger,\ S_v,\ F^\dagger,\ S_v,\ F,\ B,\ F^\dagger,\ S_v,\ F,\ S_v.
\qquad(7)
```

It contains four complete selected swaps and four catalogue appearances.
Code zero is assigned the actual empty word $`U_0=I_S`$. Therefore
the complete zero-code operation is exactly identity on every occupied
input. An approximate identity word would not give this property.

For the selected swaps use [R91's complete paired router](PARALLEL_BILINEAR_WORD_QUERY.md#6-a-faster-selected-slot-router-with-a-larger-dirty-allocation).
An arbitrary q-bit program R is toggled by the exact query $`Q_f`$.
The paired translation $`\widehat T_R`$ XOR-permutes both
$`\mathcal A`$ and a separate arbitrary auxiliary array
$`\mathcal E`$. The chronological echo

```math
Q_f,\quad\widehat T_R,\quad Q_f^\dagger,\quad\widehat T_R^\dagger
```

is exactly $`\widehat T_f`$ and returns R and all query/router
helpers. Although it temporarily permutes $`\mathcal E`$, the
complete selected swap

```math
\widehat T_f^\dagger\mathrm{SWAP}(S,\mathcal A_0)\widehat T_f
        =S_f(S,\mathcal A)\otimes I_{\mathcal E}          \qquad(8)
```

returns that array too. Every occurrence in (7) is the complete word
(8), with actual native inverses. The same returned program, auxiliary
array and private router pool can therefore serve successive swaps.
P excludes every payload, array, program and helper position.

## 4. Native resources and the final-accuracy code

Suppose each fixed $`U_v`$ has separate bounds
$`t_{\rm loc},g_{\rm loc},d_{\rm loc}`$ for T-count, elementary
Clifford count and T-depth. Its complete work is inside its own slot.
All V slot words in F run in parallel on disjoint wires. Their counts
add; their T-depth is the maximum. For $`Q=2^a`$ address rows,
an exclusive query pool $`B_Q\ge C(a+2)`$ and R91 give

```math
\begin{aligned}
T_i&=O\!\left(\sqrt{Qq}+\frac{Qq}{B_Q}
                      +Vt_{\rm loc}+wV^2\right),\\
G_i&=O(Qq+Vg_{\rm loc}+wV^2),\\
D_i&=O\!\left(\log_2(a+2)+\frac{Qq}{B_Q^2}
                         +d_{\rm loc}+\log_2(q+2)\right).
\end{aligned}                                           \qquad(9)
```

These bounds include all four catalogue appearances, all queries,
actual inverses, fixed swaps and both B gates. They do not replace
local depth by a local gate-count bound. For Q one the query is a fixed
Clifford XOR word. A zero-bit program needs no query or routing.

Choose an absolute constant $`C_{\rm loc}`$ and

```math
b_{\rm loc}=\left\lceil\max\left\{17(\kappa+h+8),
                 C_{\rm loc}\sqrt{2^{h+1}\kappa}\right\}\right\rceil .
\qquad(10)
```

The first term is precisely the real R79 threshold at dimension h+1.
Its width terms are bounded by the second term, giving

```math
\begin{aligned}
t_{\rm loc}&=O(\sqrt{2^{h+1}\kappa}+\kappa),\\
g_{\rm loc}&=O(2^{h+1}\kappa),\\
d_{\rm loc}&=O((\kappa+\log_2(h+3))\log_2(h+\kappa+3)).
\end{aligned}                                           \qquad(11)
```

For each positive h-bit group frame, give every depth-j angle
$`k_j=\kappa+h-j+6`$ bits. The integer a represents
$`\theta'=(\pi/2)a/2^{k_j}`$. The retained certified rounding
interface terminates with $`|\theta-\theta'|\le\pi2^{-k_j}`$
without deciding equality at a grid boundary. Zero is exact and
$`\pi/2`$ may be approached from below. The direct-sum layer norm
and chronological hybrid give

```math
\|W-W_v\|\le\pi\sum_{j=0}^{h-1}2^{-k_j}
                   \lt\frac{2^{-\kappa}}{16},\qquad
q_h=\sum_{j=0}^{h-1}2^jk_j\le(\kappa+8)2^h.             \qquad(12)
```

Every code is a legal positive angle tuple. Its all-zero tuple is
identity, consistent with $`U_0=I`$. Compile every other entry by
R79 at kappa. There are finitely many computable tuples, so this
preprocessing terminates. Equations (3) and (12) give group error
less than $`3\cdot2^{-\kappa}`$ on initialized c one, d zero.
The coherent address takes the maximum error; there is no factor V.

## 5. The retained opening and a constant number of groups

First work above one absolute dimension cutoff, chosen also so
$`\ell\ge16`$, and at
$`b\ge\sqrt N/n`$, uniformly for $`6\le L\le3\ell`$.
Reuse exactly [R93 Sections 2–3](GROUPED_CATALOGUE_PRECISION.md#2-one-uniform-coarse-reduction):
the magnitude/sign reduction, one uniform R89 coarse stage, and its
prefix-selected R78 early decoder. Their respective errors are at most
$`\eta/512,\eta/2048,\eta/128`$. Their total resources are

```math
T_{\rm open}=O(\sqrt N+N/b),\quad G_{\rm open}=O(N),\quad
D_{\rm open}=O(N/b^2+\ell\log_2(\ell+2)).
```

The first coarse reserve is $`m_1=O(\ell^3)`$. The next decoder
uses local precision $`L+4`$ and reserve
$`M=\lceil A_E[L+4+\log_2(m_1+L+6)]\rceil`$ when
$`M\le m_1/2`$. If that test fails, omit the decoder; its failed
inequality gives the same residual-height bound

```math
M\le C_M[L+\log_2(\ell+2)]=O(\ell).                      \qquad(13)
```

All source queries keep the original long prefix: the large row count
is $`2^{n-M}`$; only local instance multiplicities use local height.
The absolute opening widths and minimum query pools fit the wide
branch after the cutoff. Actual saved-predicate inverses complete before
the next stage; initialized-isometry composition requires no reset.

If $`M\le3`$, omit the middle. Otherwise reserve the bottom three
logical bits $`Z=(c,d,G)`$ and partition the other $`M-3`$ layers
into consecutive groups of height
$`h_{\max}=\max\{1,\lfloor\ell/8\rfloor\}`$, except possibly a
short first group. For group i let $`h_i`$ be its height and
$`s_i=3+(R-i)h_{\max}`$ its lower suffix size. Since (13) is
uniform and $`h_{\max}\ge\ell/16`$ for $`\ell\ge16`$,
$`R\le R_*`$ for one absolute constant. Set

```math
\kappa=L+\lceil\log_2(R+1)\rceil+10,\qquad
\epsilon=2^{-\kappa},\qquad
3R\epsilon\le3\eta/1024.                                \qquad(14)
```

Here $`\kappa=O(L)=O(\ell)`$ with absolute constants. A short
group uses its actual h in (10)–(12). No residual correction, endpoint
labels or retained failure history is appended to the approximate
enlarged payload.

## 6. Exact guards and composition through logical leakage

Save $`H=[Z=000]`$ in the sole external clean bit, then apply
CNOT from H to c. On the active ideal sector, c is one and d,G are zero.
For group i let $`Y_i`$ be the logical bits below $`X_i`$ and
above Z, and $`P_i`$ the whole prefix above $`X_i`$. The retained
[two-helper logarithmic-depth conjunction](../../docs/T_DEPTH_COMPILER.md#returned-helpers-for-logarithmic-depth-predicates)
gives an exact predicate $`C_i`$ toggling G by
$`e_i(H,Y_i)=H[Y_i=0]`$. Its controls exclude Z, $`X_i`$
and every dirty payload position. Its completed word is identity on H
zero on arbitrary occupied work.

Let $`f_i(P_i,H,G)`$ be the rounded group code when H=G=1,
and zero otherwise. The complete portable commutator $`M_i`$ preserves
H, $`P_i,Y_i,G`$ and acts only on the payload, returning its
array/router/query work exactly. The physical group is
$`C_i^\dagger M_iC_i`$. With u the H bit, v the original G bit,
and $`\widehat C_t=U_tBU_t^\dagger B^\dagger`$, its full-space
action is precisely

```math
\sum_{u,p,y,v}|u,p,y,v\rangle\langle u,p,y,v|
 \otimes\widehat C_{\,f_i(p,u,v\oplus u[y=0])}
 \otimes I_{\rm returned\ work}.                         \qquad(15)
```

Thus it preserves the original computational G value, even when the
payload leaks into c,d or becomes entangled with $`D_*`$.
It need not factor as a G-independent operation on arbitrary active G.
On the intended G-zero input it returns G exactly to zero. On H zero,
the selected native word is exact identity, including occupied c,d,G
and arbitrary arrays. A superposition of arbitrary G values on H one
may acquire correlations; no stronger clean-return claim is used.

Each current prefix is freshly queried. Earlier groups may change bits
in a later prefix, but the current prefix remains fixed throughout its
own group. On the original $`Z=0`$ sector the ideal active operation is exactly the
prescribed h-layer Hopf block on $`X_i`$ when $`Y_i=0`$,
and identity otherwise, including all columns and internal selectors.

On H one, let $`J_{\rm mid}`$ insert c=1,d=0,G=0. The complete
group obeys
$`\|U_iJ_{\rm mid}-J_{\rm mid}(W_i\otimes I_b)\|\le3\epsilon`$.
Unitary telescoping gives

```math
\begin{aligned}
&U_R\cdots U_1J_{\rm mid}-J_{\rm mid}W_R\cdots W_1\\
&=\sum_i U_R\cdots U_{i+1}
 (U_iJ_{\rm mid}-J_{\rm mid}W_i)W_{i-1}\cdots W_1 .
\end{aligned}                                           \qquad(16)
```

The norm is at most $`3R\epsilon`$. Every bracket is evaluated on
an ideal initialized state; later full unitaries propagate all earlier
actual leakage. $`D_*`$ can be reused without reinitialization.
On H zero the whole middle is exact identity, so the same bound holds
for coherent H sectors. Undo the c preparation and the original H
predicate using their actual inverses. The ideal middle preserves Z
zero, so these restore the ideal work. On leaked outputs they preserve
the norm bound, even if the actual H output is nonzero.

## 7. Simultaneous width and the native sums

The following table is a sufficient simultaneous middle reservation.
Every local slot word contains all its own work.

| Field | Wires | Required separation or return |
|---|---:|---|
| H | One external clean bit | Unchanged address, outside payload |
| c,d,G and current logical data | Original logical wires | Z is conditional work only |
| Logical local dirty pool | $`b_{\rm loc}`$ | Entire pool travels with S |
| Dirty program R | $`q_h`$ | Outside payload and address; exact return |
| Catalogue and paired auxiliary arrays | $`2wV_h`$ | Disjoint complete slots; exact return |
| Router grid, indicators and helpers | $`C_{\rm router}wV_h^2`$ | Outside every live array; exact return |
| Global query pool | $`B_Q=\lfloor b/2\rfloor`$ | Exclusive; outside every other field |
| Predicate helpers | Fixed polynomial in M | Separate arbitrary dirty work; exact return |

In particular, the dirty peak is bounded by

```math
b_{\rm loc}+q_h+2wV_h+C_{\rm router}wV_h^2+B_Q
       +\text{predicate helpers}.
```

Slots do not borrow another slot's work.
All slot d inputs may be occupied; (5) uses their physical unitaries only.

Uniformly through $`L\le3\ell`$, equations (10)–(14) imply

```math
\begin{aligned}
q_h&=O(\ell(n+2)^{1/8})=o(n),\\
w&=O(\ell+(n+2)^{1/16}\sqrt\ell),\\
\log_2(wV_h^2+V_ht_{\rm loc}+V_hg_{\rm loc}+b_{\rm loc}+q_h)
   &\le C\ell(n+2)^{1/8}+C\ell=o(n).
\end{aligned}                                           \qquad(17)
```

After one absolute cutoff the sum of these standalone counts and peak
widths, including the constant number of groups, is at most
$`2^{n/4}`$. These are $`N^{o(1)}`$ quantities, not
subpolynomial quantities in n. Enlarge the cutoff so these terms,
all predicate work and fixed constants fit in at most b/4 whenever
$`b\ge\sqrt N/n`$. Then the stated exclusive pool satisfies
$`B_Q=\Theta(b)\ge C(n+2)`$. Earlier and later stages reuse
returned allocations sequentially; their peaks are maximized.

The prefix above a group has $`n-s_i-h_i`$ bits. Including H,G
only multiplies its row count by four, so

```math
Q_i\le4\,2^{n-s_i-h_i},\qquad
Q_iq_i\le C N\kappa2^{-s_i}.
```

Since $`s_i=3+(R-i)h_{\max}`$ and $`h_{\max}\ge1`$,

```math
\sum_iQ_iq_i=O(NL),\qquad
\sum_i\sqrt{Q_iq_i}=O(\sqrt{NL}).                        \qquad(18)
```

The two geometric ratios are $`2^{-h_{\max}}`$ and
$`2^{-h_{\max}/2}`$. The short first group satisfies the same
estimate: its $`2^{h_i}`$ code factor cancels its prefix factor.
No long address is replaced by the shorter program-register address.

With $`B_Q=\Theta(b)`$, (9), (17), and (18) give
$`T_{\rm mid}=O(\sqrt{NL}+NL/b)`$ and
$`G_{\rm mid}=O(NL)`$. There are at most $`R_*`$ groups,
$`d_{\rm loc}=O(\ell\log_2(\ell+2))`$ and
$`\log_2(q_h+2)=O(\ell)`$. Guard predicates have at most
$`M=O(\ell)`$ controls, polynomial count/width and logarithmic
T-depth. Thus

```math
D_{\rm mid}=O(NL/b^2+\ell\log_2(\ell+2)).                \qquad(19)
```

The separate local count/depth ledger and R91's fast complete selected
swap are both necessary for this estimate.

## 8. Bottom, total error, original widths and complex completion

The deepest $`r=\min(M,3)`$ logical layers use the retained
[R72 smaller-error bottom](ONE_CLEAN_DEPTH.md#the-final-logical-layers-at-smaller-error),
with $`q_j=L+j+7`$, $`m_j=L+j+8`$ and
$`Q_j=N2^{-j}`$. Its base $`L+n+9`$ and a
$`\Theta(b)`$ query pool fit the wide branch. The retained
weighted sums over at most three layers give

```math
T_{\rm bot}=O(\sqrt{NL}+NL/b),\quad G_{\rm bot}=O(NL),\quad
D_{\rm bot}=O(NL/b^2+\ell),
```

and error less than $`15\eta/64`$. Its long-prefix queries and
native source costs are included. It reuses the sole external flag;
initialized-isometry telescoping propagates earlier actual leakage.
The complete wide real error is at most

```math
\frac{\eta}{512}+\frac{\eta}{2048}+\frac{\eta}{128}
 +\frac{3\eta}{1024}+\frac{15\eta}{64}
   =\frac{507\eta}{2048}\lt\frac{\eta}{4}.               \qquad(20)
```

Omitted stages add zero. For $`b\lt\sqrt N/n`$, use the entire
retained R79 compiler at the original precision and original literal
width. Since
$`NL/b^2\gt n^2L\ge nL\ge F_{79}(n,L)`$,
its additive depth is absorbed by the width term. No stricter-precision
fallback is imposed at the real minimum width.

For n below the common absolute cutoff $`n_0`$, use that same
whole original-precision compiler at every width. On the internal real
range, $`F_{79}(n,L)\le nL\le3n_0\ell\log_2(\ell+2)`$.
This is one absolute enlargement, with no L-dependent hidden maximum.
Empty precision ranges in small dimensions require nothing.

For the complex family on the large wide branch, apply the real result
at $`L+1\le3\ell`$ and then the
[retained literal diagonal](ONE_CLEAN_DEPTH.md#7-a-one-clean-diagonal-and-complex-composition)
at L+1. The real threshold is exactly
$`17((L+1)+n+7)=17(L+n+8)`$; the diagonal has the same base.
It contributes error at most $`\eta/2`$ and costs

```math
T_D=O(\sqrt{NL}+L+NL/b),\quad G_D=O(NL),\quad
D_D=O(NL/b^2+L+\ell).
```

The real error from (20) is less than $`\eta/8`$ after the shift.
Literal specified phases, including zero-amplitude leaves, are retained.
For narrow or finite complex cases use the whole original-precision
R79 complex compiler; do not append a second diagonal. The original
widths, one-clean contract and stated bounds follow in every branch.

## 9. Scope, dependencies and bounded evidence

The proof uses the established R79 theorem only at local dimension
h+1. Its dependency graph is acyclic: R79 predates R90/R91 and this
construction; the opening uses only R93's separately proved early
interfaces. The older R93 residual-source terminal is not invoked.
Arrays, router work and programs return exactly; logical c,d and
$`D_*`$ return within the charged initialized-isometry error.

The [bounded operator fixtures](../../tests/test_portable_workspace_catalogue.py)
check complete payload transfer with native words and actual inverses,
coherent selection of noncommuting entries, controlled Hopf marker
completion, the inverse and commutator estimates under leakage, exact
guard return with inactive occupied inputs, and a four-mode
changing-selector example. A partial-payload counterexample
leaves a dirty-slot residue on occupied inputs; it refutes that exact
full-space transfer interface, not every possible partial-routing task.
The separate [R91 fixtures](../../tests/test_parallel_bilinear_word_query.py)
check its dirty-program and paired-router identities. Small operators
do not establish the asymptotic resource or global-width proof above.

This is a growing-precision upper bound. It does not assert optimal
general depth, efficient preprocessing, a practical crossover, arbitrary
complex-unitary synthesis, or improvement outside its stated precision
range. The general depth question remains the destination.
