# From native query depth to complete-frame depth

[Current all-precision upper theorem R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Global fast stage](GLOBAL_BORROWED_PRECISION_STAGE.md) ·
[Paired permutation wrappers](PAIRED_PERMUTATION_ROUTING.md) ·
[Native query baseline R91](PARALLEL_BILINEAR_WORD_QUERY.md)

The general upper bound still contains logarithmic depth in the number
of logical qubits. A faster exact dirty query would address this cost,
but substituting it into the old circuit alone leaves logarithmic
barrel-routing depth. Paired permutation wrappers remove that residual
cost. Four decreasing early stages then make the finite terminal small
enough to prove the following conditional reduction.

## 1. A precise query hypothesis and conditional theorem

Let $`d:\mathbb N\to[1,\infty)`$ be nondecreasing. Assume that,
for every a-bit address, every $`Q=2^a`$-row q-bit table with
$`q\ge1`$, and every extra dirty allocation $`B\ge C_Q(a+2)`$,
a terminating classical compiler supplies a native exact XOR query with

```math
T=O\!\left(\sqrt{Qq}+Qq/B\right),\qquad G=O(Qq),
\qquad D_T=O\!\left(d(a)+Qq/B^2\right).                \qquad(1)
```

Its additional width is at most B. The address is read-only and can
occur only in Clifford controls throughout. Outputs and helpers are
arbitrary; all helpers return exactly with literal phase one. The
identity holds with arbitrary reference systems. All constants are
absolute. For a zero-bit address the constant row is a Clifford word.
No efficiency bound on this classical compilation is required.
The hypothesis includes no initialized work, supplied state, evaluator,
measurement or reset. R91 proves this interface with
$`d(a)=1+\log_2(a+2)`$; a smaller d is an additional mathematical
premise, not a result of this chapter.

**Theorem R98 (conditional fixed-precision transfer).** Under (1), for every fixed
integer $`L\ge6`$, every $`n\ge1`$ and every
$`b\ge17(L+n+7)`$, every complete prescribed real Hopf frame W has
one native circuit with one external clean qubit, at most b arbitrary
dirty qubits, and

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad N=2^n,
```

```math
D_T=O_L\!\left(N/b^2+d(2n+10)
                         +\log_2\log_2(n+4)\right).   \qquad(2)
```

Its full initialized-isometry error is at most $`2^{-L}`$. The same
statement holds for the prescribed complex magnitude completion at
$`b\ge17(L+n+8)`$. These bounds hold on the same physical circuit.
The constants may depend on fixed L and on the absolute constants of
(1), not on n,b or the target. No uniform growing-L improvement is
claimed. The T-count remains worst-case optimal on these allocations.

In particular, $`d(a)=o(\log(a+2))`$ would give sublogarithmic
general depth at $`b=\Theta(\sqrt N)`$ and fixed L. This is a
sufficient route toward optimal depth, not a necessity theorem or an
unconditional improvement to R96. R79 remains selected Result D.

The error contract includes every logical column, literal signs and
phases, singular and zero coefficients, all work leakage and return,
arbitrary dirty inputs and external references. Certified preprocessing
terminates without a running-time claim. Every query, source, control,
recovery and actual inverse below is charged.

## 2. Audit of one modified early stage

Use the fast-reference variant of the global stage, with remaining
height m, preserved original prefix of length $`n-m`$, and stage
precision l. Keep its notation

```math
K=m+l+2,\quad t=\lceil\log_2(m+2)\rceil,\quad
\chi=\lceil\log_2(l+2)\rceil+t^3,\quad r=\lceil A\chi\rceil.
                                                               \qquad(3)
```

Run a stage only if $`r\lt m/2`$. Its original suffix of r bits
and all source/metadata banks are unchanged. Replace complete long
queries by (1), and selected layout wrappers by the paired routing
lemma. All other native components are retained.

The global stage's Section 5 lists a constant number of chronological
query positions, including all actual inverses and sparse-correction
recovery calls. At one position at most $`CK^{d_0}`$ instances
run in parallel. Their outputs, grid targets and helper pools are
disjoint; their only shared inputs occur through read-only Clifford
controls. Thus they contribute a maximum of their additive depths,
not $`K^{d_0}d(a)`$. Their elementary Clifford counts still add.

Fix all query constants, router polynomial degrees and stage exponents
before choosing A. Enlarge the fixed exponent C to cover these native
substitutions and put

```math
\overline Q=2^{n-r+C\chi},\qquad
B_j=\left\lfloor b/(C_1K^{d_0})\right\rfloor.
                                                               \qquad(4)
```

The inherited sufficient simultaneous allocation, with enlarged fixed
constants, is

```math
b\ge8W_{\rm loc}+C K^{d_w}(n+\chi+2),\qquad
W_{\rm loc}\le2^{C\chi}.                              \qquad(5)
```

It ensures $`B_j\ge C_Q(a_j+2)`$ for every actual address.
Every long table has $`a_j\le n-r+C\chi`$, including ordinary
source anchors, augmented coarse-block readout, matching data, source
masks and scalar calibration programs. Choose $`A\gt C`$; then
all these addresses are at most n. Query multiplicities, pool splitting
and width penalties are absorbed in one fixed degree d_0'. No new
source field becomes part of a table address.

The remaining depth costs are as follows.

| Native component | Depth after substitution | Why no logarithmic-in-m term remains |
|---|---:|---|
| Long-query positions | $`O(d(n)+K^{d_0'}\overline Q/b^2)`$ | Constant sequential positions; private simultaneous pools |
| Short standalone tables | $`O(\log(\chi+2))`$ | Retained R91 queries on at most $`C\chi`$ bits |
| Ordinary borrowed-precision sources and coefficient banks | $`O(l)`$ | Constant many source positions; original logarithmic initialized indices |
| Capped exact reference words | $`O(l+\log(t+2))`$ | All bins run on disjoint coin/private-scratch registers |
| Selected layout/body/unlayout positions | $`O(\log(t+2))`$ | Polynomially many layouts selected by $`O(t)`$ bits; paired routing |
| Saved suffix and complete used-work predicates | $`O(\log(\chi+2))`$ | They test only r or $`O(\chi)`$ bits |

The long source anchors are in the first row, not the last row. The
last row never tests the logical output, arbitrary dirty cores, idle
suffix padding or new paired routing arrays. All retained selector
fields are computationally preserved by their actual bodies. The paired
arrays are untouched through each corresponding body and canceled by
its actual inverse layout, including on occupied inactive work.
The paired routing proof specifies the fine-bin/interval and per-block
metadata selectors separately. All original parameter loads are retained,
and each new permutation table reproduces the old router on every
physical selector value, including malformed active metadata; no
off-anchor extension has been changed.

The original native proof therefore gives, on the same circuit,

```math
\begin{aligned}
T_E&\le CK^{d_0'}(\sqrt{\overline Q}+\overline Q/b+1)+2^{C\chi},\\
G_E&\le CK^{d_0'}(\overline Q+1)+2^{C\chi},\\
D_E&\le C\left[l+\log_2(\chi+2)+d(n)
                         +K^{d_0'}\overline Q/b^2\right].
\end{aligned}                                                     \qquad(6)
```

Its full error remains at most $`2^{-l}/16`$. All substitutions
are exact full-space identities; the paired inverse permutation on
each E cancels before any later source position. Existing fast-reference
comparisons, source completion, normalization-two recovery, H-zero
identity and final actual saved-predicate inverse are consequently
unchanged. This is an audit of a native word, not an inference from
commuting ideal projectors.

## 3. Four stages at one fixed precision

Set $`\eta=2^{-L}`$ and $`l=L+6`$. First make the retained finite
whole-family magnitude/sign reduction at error at most $`\eta/512`$.
Its two full-address Boolean sign phases have
$`T=O(\sqrt N+N/b)`$, $`G=O(N)`$, and depth
$`O(d(n+2)+N/b^2)`$ by (1) and the exact dirty phase echo. Original
singular continuations are included in this finite reduction.

Fix the constants in (3)–(6), then choose A so that

```math
r-C\chi\ge(2d_0'+80)\chi.                            \qquad(7)
```

Starting from $`m_0=n`$, perform at most four stages. If
$`r_i\ge m_i/2`$, stop and pass the entire remaining height to the
terminal. For fixed L such a stop has bounded height: the inequality

```math
m\le2\left\lceil A\left[\lceil\log_2(l+2)\rceil+
\lceil\log_2(m+2)\rceil^3\right]\right\rceil
```

has only finitely many integer solutions. Otherwise put
$`m_{i+1}=r_i`$ and continue.

For $`\ell=\log_2(n+2)`$, four running stages give

```math
\begin{aligned}
m_1&=O_L(\ell^3),\\
m_2&=O_L(\log_2^3(\ell+2)),\\
m_3&=O_L((\log_2\log_2(\ell+4))^3),\\
m_4&=O_L((\log_2\log_2\log_2(\ell+16))^3).
\end{aligned}                                                     \qquad(8)
```

The early-stop terminal has a fixed L-dependent height and obeys the
same eventual bounds. Four is a fixed number, not an uncharged
recursive loop. In every executed stage $`\chi=O_L(\ell^3)`$,
so all standalone counts and widths are
$`2^{O_L(\ell^3)}=N^{o(1)}`$. The minimum-pool term in (5) is
polynomial in n at fixed L. Thus for $`b\ge\sqrt N/n`$ both fit
inside the stated physical reservations after one L-dependent dimension
cutoff. They are reused sequentially only through completed native words.

Equation (7), with $`\chi\ge c\log_2K`$ and an enlarged fixed
margin if necessary, absorbs every large-table polynomial into the
global count and width terms. Summing the four ledgers yields

```math
T_{\rm open}=O_L(\sqrt N+N/b),\quad G_{\rm open}=O_L(N),
```

```math
D_{\rm open}=O_L(N/b^2+d(n+2)+\log_2(\ell+2)).        \qquad(9)
```

The stage errors sum to at most
$`4\cdot2^{-(L+6)}/16=\eta/256`$. Earlier leakage is propagated
through later full unitaries, not removed by assuming fresh flag inputs.

## 4. The paid whole terminal

Let m be the final height. Use R90's whole rounded real-frame catalogue,
with the unchanged preserved upper prefix of length $`n-m`$.
Its depth-dependent angle code has

```math
q\le(L+12)2^m,
\quad Q=2^{n-m},
\quad Qq=O_L(N),                                     \qquad(10)
```

and rounding error at most $`\eta/256`$. Each catalogue entry
has determinant one, including the prescribed singular cases after
the finite sign reduction. The fixed commutator factors have native
word length $`S=2^{O_L(m)}`$ at error $`O(\eta)`$. Enumerating
and certifying actual finite native words uses R90's strict-certificate
construction; no efficient eigenbasis computation or supplied unitary
is assumed.

Put $`u=\log_2(\ell+2)`$. Equation (8) gives
$`m=O_L((\log_2\log_2(u+4))^3)=o_L(\log u)`$. Hence

```math
q+S=2^{O_L(m)}=u^{o(1)}=O_L(u).                       \qquad(11)
```

The original R90 slot-translation router suffices here; there is no
need to assume a new group law for the paired permutation lemma.
The two live catalogue arrays, their fixed words, the dirty program,
and their routers have total count/width

```math
2^q[S+mq]+q+2m2^q=n^{o(1)}=o(\sqrt N/n).             \qquad(12)
```

At sufficiently large n reserve the program and both arrays within
$`b/4`$, and an exclusive query pool $`B=\lfloor b/2\rfloor`$,
disjoint from every live catalogue wire. The minimum condition in (1)
holds because $`B\gg n`$. Replace its complete long-prefix query
by (1). R90's exact two-array commutator and actual slot-swap inverses
remain unchanged. Their bounded number gives

```math
T_{\rm term}=O_L(\sqrt N+N/b),\quad G_{\rm term}=O_L(N),
```

```math
D_{\rm term}=O_L(N/b^2+d(n)+q+S)
             =O_L(N/b^2+d(n)+\log_2(\ell+2)).           \qquad(13)
```

The arrays return exactly even for approximate fixed words, because
each unwanted action is canceled by its own actual inverse in the
disjoint-array commutator. The complete selected-factor error is the
maximum over entries, at most $`\eta/32`$, not a sum over the
catalogue. Including its rounding gives terminal error at most
$`9\eta/256`$. The terminal uses no clean qubit. For m zero it
is omitted; for Q one its query is a fixed Clifford mask.

## 5. Complete contract, all widths and complex phases

The total real error in the wide large-n branch is at most

```math
\eta/512+\eta/256+9\eta/256
    =21\eta/512\lt\eta.                              \qquad(14)
```

All original prefix sectors take a maximum. The ideal intermediate
operators return the same sole clean flag; ordinary initialized-isometry
telescoping includes actual flag, source and dirty-core leakage. Each
pair of new routing arrays returns exactly, and each actual inverse is
charged. This proves the complete real operator guarantee in this branch.

For $`b\lt\sqrt N/n`$, use the original-width R96 circuit:
$`N/b^2\gt n^2`$ absorbs its fixed-L additive
$`L+\log_2(n+2)`$. For n below the common L-dependent cutoff,
the same original-precision circuit has bounded additive depth at fixed
L, absorbed in (2) because $`d\ge1`$. These fallbacks never request
a stronger precision at the literal minimum width. They preserve the
one-clean count, Clifford and full error statements.

For the wide large-n complex branch, run the real construction at L+1
and append R72's one-clean literal diagonal at L+1. Its explicit word
has twenty complete mask queries on the n-bit logical address and
forty-five $`O(L)`$ source appearances, with no suffix predicate.
Replacing those queries by (1) gives depth
$`O_L(N/b^2+d(n)+1)`$, the same fixed-L counts and full error at
most $`\eta/2`$. Its precision-sized dirty core and an exclusive
$`\Theta(b)`$ query pool fit simultaneously in this branch. The
real threshold at L+1 is exactly $`17(L+n+8)`$. Both factors reuse
the same returned clean flag and dirty allocation; prescribed phases
at zero-amplitude leaves remain literal. Small and narrow cases use
whole original-precision complex R96. This proves (2) for both families.

## 6. What this reduction does and does not resolve

The paired routing lemma is unconditional. The general sublogarithmic
conclusion requires the exact query hypothesis (1) with a genuinely
smaller d. The present R91 instance recovers the logarithmic fixed-L
law and does not improve R96. The argument identifies a sufficient
native mechanism and charges the other costs that would remain after
it, including four stages and the complete terminal. It neither proves
that queries are necessary for every frame representation nor turns
an architecture-specific query obstruction into a general lower bound.

The concrete resumption question is: can an arbitrary Boolean/word
table be queried in the complete arbitrary-dirty sense of (1), at
sublogarithmic address T-depth, with the same simultaneous count and
minimum-pool guarantees? A matrix-element or trace encoding, a query
correct only on clean outputs, or a construction whose address receives
non-Clifford gates shared across the required parallel banks does not
satisfy this sufficient interface.
