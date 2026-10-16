# Logarithmic depth at inverse-polynomial accuracy

[Executed design](FULL_FLAG_DEPTH_PROMPT.md) ·
[Portable workspace catalogue R94](PORTABLE_WORKSPACE_CATALOGUE.md) ·
[Coarse flags and sparse correction](COARSE_FLAG_SPARSE_CORRECTION.md) ·
[Local compiler with relaxed counts](RELAXED_COUNT_LOCAL_DEPTH.md)

**Result R95.**
The two improvements have different roles. A coarse mathematical frame,
compiled to final accuracy and corrected by sparse sources, removes the
opening's repeated spectral search. A local compiler with a much larger
allowed dirty allocation removes its repeated precision charge. R94's
complete-workspace catalogue makes those large local resources affordable
in the global circuit.

## 1. The theorem and its scope

Let $`N=2^n`$, $`\ell=\log_2(n+2)`$, and $`\eta=2^{-L}`$, with
integers $`n\ge1`$ and $`6\le L\le2\ell`$. Every prescribed complete
real Hopf frame has, for every $`b\ge17(L+n+7)`$, a native Clifford+T
circuit with one external initialized qubit and at most b arbitrary
dirty qubits such that

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O\!\left(\sqrt{NL}+L+\frac{NL}{b}\right),\qquad
G=O(NL),\qquad
D_T=O\!\left(\frac{NL}{b^2}+L+\ell\right).
\qquad(1)
```

The prescribed complex magnitude completion satisfies the same bounds
at $`b\ge17(L+n+8)`$. All constants are absolute on this range.
At $`L=\Theta(\ell)`$ and $`b=\Theta(\sqrt{NL})`$, (1) gives
$`O(\ell)`$ T-depth, improving R94's $`O(\ell\log_2(\ell+2))`$.
The existing worst-case count lower bounds still match the T-count;
no matching depth lower bound is asserted.

The norm includes every logical column, literal signs and phases,
zero-coefficient continuations, all work leakage, arbitrary dirty inputs,
and external references. Inverses are actual reversed native words.
Classical preprocessing terminates but need not be efficient. No
measurement, reset, supplied precision state, QRAM, additional external
clean word, or free coherent evaluator is used. Optimal general depth
and the larger precision range remain open.

We prove the real construction internally on $`6\le L\le3\ell`$,
so the complex case can pay its one-bit precision shift.

## 2. The two component contracts

The [coarse-flag component](COARSE_FLAG_SPARSE_CORRECTION.md) applies to
an m-bit remaining frame with an unchanged global prefix, at stage
precision l. After fixing all native constants, it uses an original
protected suffix of size

```math
r=\left\lceil A\left(l+\lceil\log_2(m+l+2)\rceil^3\right)\right\rceil.
\qquad(2)
```

When $`r\lt m`$, it compiles the earlier layers with full initialized-isometry
error at most $`2^{-l}/16`$. Its global-prefix query positions have
constant sequential multiplicity. At $`m=O(\ell^3)`$, $`l=L+O(1)`$,
it leaves $`r=O(\ell)`$ original logical bits and has depth
$`O(NL/b^2+L+\ell)`$ on the wide branch below. Its standalone
tables and dirty work are $`N^{o(1)}`$ and explicitly charged.

The [local component](RELAXED_COUNT_LOCAL_DEPTH.md) supplies, for an
h-bit real Hopf frame at precision $`2^{-\kappa}`$, one clean qubit,
depth $`O(h+\kappa)`$, and separate T-count, Clifford-count and
arbitrary dirty-width bounds of the form

```math
t_{\rm loc},g_{\rm loc},b_{\rm loc}
 \le 2^{P(h+\kappa+1)},
\qquad d_{\rm loc}\le C(h+\kappa+1),
\qquad(3)
```

for one fixed polynomial P. Its guarantee is the same full
initialized-isometry norm, rather than only a state-preparation column.
Every physical local workspace wire belongs to its payload. Its proof
uses only earlier source, query, R79, and R94 transfer identities; it
does not invoke (1).

The exponential local allowance in (3) is intentional. Neither local
optimal T-count nor the global literal dirty threshold is claimed for
the local circuit itself. The catalogue below prices its entire physical
unitary and every occupied slot.

## 3. Opening to a logarithmic original suffix

Work first above one absolute dimension cutoff and at
$`b\ge\sqrt N/n`$. Reuse R93's magnitude/sign reduction and its
first uniform R89 stage, as recorded in
[R94 Section 5](PORTABLE_WORKSPACE_CATALOGUE.md#5-the-retained-opening-and-a-constant-number-of-groups).
Their errors are at most $`\eta/512`$ and $`\eta/2048`$.
The first reserve is $`m_1=O(\ell^3)`$, uniformly through
$`L\le3\ell`$. These two retained operations already have
$`O(N/b^2+L+\ell)`$ depth; the later binary search is omitted.

Apply (2) with $`m=m_1`$ and $`l=L+4`$ if its r is smaller than
$`m_1`$. Otherwise omit this stage. In either case the remaining
height M satisfies

```math
M\le C_M\left[L+1+\log_2^3(\ell+2)\right]=O(\ell).
\qquad(4)
```

For the omitted case, $`m_1\le r`$ gives this bound directly, since
$`m_1=O(\ell^3)`$. The retained original global prefix occurs in
every required query; its address is never replaced by the local height.
The new stage has error at most $`\eta/256`$. All opening sources,
private query pools, exact saved predicates, actual inverses, and
standalone costs are included in the component's ledger. Their total
fits

```math
T_{\rm open}=O(\sqrt{NL}+NL/b),\quad G_{\rm open}=O(NL),\quad
D_{\rm open}=O(NL/b^2+L+\ell).
\qquad(5)
```

## 4. Portable local words at final accuracy

If $`M\le3`$, omit this entire guarded catalogue stage and use the
bottom compiler in Section 6 directly. Otherwise use R94's original
three bottom bits $`Z=(c,d,G)`$, exact guards,
and complete-payload commutator. Split the other $`M-3`$ layers into
groups of maximum height $`h_{\max}=\lfloor\ell/8\rfloor`$, with
a possible shorter first group. Increase the absolute cutoff so
$`\ell\ge16`$. Equation (4) gives a fixed absolute upper bound
$`R_*`$ on the number R of groups.
Set

```math
\kappa=L+\lceil\log_2(R+1)\rceil+10.
\qquad(6)
```

Retain R94's certified angle code at this precision, with
$`q_h\le(\kappa+8)2^h`$ bits and rounded-frame error below
$`2^{-\kappa}/16`$. Every code denotes a legal positive frame.
Code zero has the actual empty physical word.

For a nonzero code v, form the controlled Hopf frame
$`A_v=\mathrm{diag}(I,W_v)_c`$ on h+1 bits, using the zero-root
construction of R94 Section 2. Compile it with (3), not with R79
directly. The complete payload is $`S=(c,X,d,D_*)`$, of width
$`w=h+2+b_{\rm loc}`$. R94's actual-inverse estimate and physical
commutator with $`X_c`$ give group error less than
$`3\cdot2^{-\kappa}`$ on c one and d zero.

The entire payload travels to the selected arbitrary dirty slot.
R94's full-space transfer identity cancels every unwanted array action
exactly, even with occupied local flags and entangled slot work. No
correctness of a local compiler on its occupied clean flag is assumed.
The paired auxiliary router array also returns exactly.

The guard construction and initialized-isometry telescoping are exactly
[R94 Section 6](PORTABLE_WORKSPACE_CATALOGUE.md#6-exact-guards-and-composition-through-logical-leakage).
The saved external H is not part of the transported payload. Each current
prefix is preserved during its own group and freshly queried for later
groups. The original G value returns exactly; c,d and the local dirty
pool return within the charged error. H zero gives literal identity on
all occupied work. Summing the errors gives at most
$`3R2^{-\kappa}\le3\eta/1024`$ without any reset.

## 5. Global count, width and depth of the catalogue

The R94/R91 ledger for a group with Q address rows, V catalogue entries,
and an exclusive query pool $`B_Q=\Theta(b)`$ remains

```math
\begin{aligned}
T_i&=O(\sqrt{Qq_h}+Qq_h/B_Q+Vt_{\rm loc}+wV^2),\\
G_i&=O(Qq_h+Vg_{\rm loc}+wV^2),\\
D_i&=O(\log_2(n+2)+Qq_h/B_Q^2+d_{\rm loc}+\log_2(q_h+2)).
\end{aligned}\qquad(7)
```

It includes four whole catalogue appearances, four complete selected
swaps, all query inverses, and returned router work. Each catalogue word
has its own complete workspace; their depths take a maximum because
they run on disjoint slots.

Here $`h\le\ell/8`$, $`\kappa=O(\ell)`$, and $`V=2^{q_h}`$.
Equation (3) gives

```math
q_h=O(\ell(n+2)^{1/8})=o(n),\qquad
\log_2(wV^2+Vt_{\rm loc}+Vg_{\rm loc}+b_{\rm loc}+q_h)
 \le O(q_h)+P'(\ell)=o(n)
\qquad(8)
```

for a fixed polynomial $`P'`$. Thus all local words, arrays, router
grids and their peak dirty widths are $`N^{o(1)}`$. After one absolute
cutoff they fit inside b/4 on the wide branch, leaving the exclusive
$`B_Q=\lfloor b/2\rfloor`$ pool and the separately priced predicate
helpers. Earlier and later stages reuse returned allocations sequentially.
All R91 query minimum-width conditions hold there.

For a group's lower suffix size $`s_i=3+(R-i)h_{\max}`$, the original
global address gives $`Q_i\le4\,2^{n-s_i-h_i}`$. Consequently

```math
Q_iq_i\le C N\kappa2^{-s_i},\qquad
\sum_iQ_iq_i=O(NL),\qquad
\sum_i\sqrt{Q_iq_i}=O(\sqrt{NL}).
\qquad(9)
```

The shorter first group's height cancels in this same product.
Equations (7)–(9) give the desired optimal global count bounds.
The new local depth is $`d_{\rm loc}=O(h+\kappa)=O(\ell)`$,
and $`\log_2(q_h+2)=O(\ell)`$. There are at most $`R_*`$ groups.
The exact guards have logarithmic depth in M. Hence

```math
D_{\rm mid}=O(NL/b^2+L+\ell).
\qquad(10)
```

## 6. Bottom, errors, fallbacks and complex phases

The deepest $`\min(M,3)`$ layers use exactly R94 Section 8's R72
bottom, including its long-prefix queries, one-clean allocation and
precision schedule. Its error is less than $`15\eta/64`$ and its
resources satisfy (1). The wide real error is therefore at most

```math
\eta/512+\eta/2048+\eta/256+3\eta/1024+15\eta/64
 =499\eta/2048\lt \eta/4.
\qquad(11)
```

For $`b\lt \sqrt N/n`$, use the whole R79 compiler at the original
precision and literal width. Its additive depth is at most nL, while
$`NL/b^2>n^2L`$. For n below the common absolute cutoff, use that
same original-precision compiler at every width. Its nL depth is
$`O(L+\ell)`$ with one absolute enlarged constant because n is
bounded. No stronger-precision fallback is imposed at the literal
minimum width.

For the large wide complex branch, use the real construction at L+1
and append R72's retained literal diagonal at L+1. The real width is
then exactly $`17(L+n+8)`$. The diagonal has the same threshold,
error at most $`\eta/2`$, and the resources (1); the shifted real
error is below $`\eta/8`$. Prescribed phases at zero-amplitude leaves
are retained. Narrow and bounded-dimension complex branches use the
whole original-precision R79 complex compiler. This establishes the
simultaneous bounds and original literal widths.

## 7. Dependencies and evidence

The proof graph is acyclic. The local component uses earlier R79 and
the algebraic R94 transfer identity, while (1) uses its established
local contract through that same identity. The global opening has its
separate sparse-correction and invariant-gradient proof. Neither
component assumes this final theorem.

The [four bounded operator tests](../../tests/test_coarse_flag_sparse_correction.py)
check signed Schur sources with a live root correction, singular branches,
normalization-two recovery with its complete work anchor and literal sign,
actual inverses, and perturbed source leakage. They reuse dense matrix
completions solely as verification objects. The retained R91 and R94
tests check native query and complete-payload identities. Independent
component and assembly reviews are recorded in the
[session outcome](FULL_FLAG_DEPTH_SESSION_RESULT.md).
These checks can detect errors in signs, root residuals, actual inverses,
occupied-work identities and error estimates; the symbolic resource
proof supplies the asymptotic conclusion. The result does not establish
optimal depth, efficient preprocessing, a practical crossover, or a
scalable native emitter.
