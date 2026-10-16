# Precision-terminal session: grouping the noncommuting middle

[Current frontier](../../docs/OPEN_PROBLEM.md#large-width-t-depth) ·
[Executed design](PRECISION_TERMINAL_PROMPT.md) ·
[Previous outcome](TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)

This session executes the user's “Now go!” on 16 October 2026, from
main commit `30c1c1adda4238ba43e5110739891ac21d740b24`, containing the
design merged in [PR #155](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/155).
Its scientific baseline is R90. Execution began at **06:19:12 UTC**;
the five-hour ceiling is **11:19:12 UTC**, including review, validation
and integration. Final manuscript drafting remains a separate task.

The main result is [R93](GROUPED_CATALOGUE_PRECISION.md). Put
$`N=2^n`$ and $`\ell=\log_2(n+2)`$. Uniformly on
$`6\le L\le2\ell`$, both prescribed frame families admit

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
```

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)
                         +L\ell/\log_2(\ell+2)\right).
```

The same circuit uses one external initialized qubit and the original
real width $`17(L+n+7)`$ or complex magnitude width $`17(L+n+8)`$.
Its error is at most $`2^{-L}`$ in the full initialized-isometry norm,
including all columns, literal phases, signs, zeros, leakage, arbitrary
dirty inputs and their references. All constants are absolute.

At the design benchmark $`L=\max\{6,\lceil\ell\rceil\}`$ and
sufficient $`b=\Theta(\sqrt{NL})`$, the asymptotic depth improves
from R79's $`O(\ell^2)`$ to $`O(\ell^2/\log\ell)`$.
The primary $`O(\ell\log\ell)`$ target and optimal general depth
remain open. This is a proved intermediate gain, not completion of the
general objective.

## 1. The changed unit of work

The opening audit found two separate costs. A whole-terminal dirty
catalogue has an enormous program alphabet once precision grows.
The shared binary terminal instead keeps the precision source compact,
but reads the long global prefix separately for each logical layer.

R93 takes a middle course: one coarse catalogue describes a complete
small noncommuting group, and the shared source corrects that group's
entire actual coarse operator. An initial uniform R89 stage reduces the
height to $`O(\ell^3)`$. A prefix-selected R78/R79 early decoder
then leaves $`O(L+\log\ell)`$ layers. Its global prefix remains
in every relevant table address; only local multiplicities become
polynomial in the remaining height.

Groups have $`g=\Theta(\log\ell)`$ target bits. Their whole
rounded programs select R90's native commutator factors on two arbitrary
dirty arrays. A generic matrix-unit dictionary corrects
$`C_i^\dagger W_i-I`$, where $`C_i`$ is the actual finite native
coarse word. Its two endpoint labels each have only g bits. The digit
query excludes every changing group target and unloads exactly after
the selected matrix unit. The next group reads its own current prefix,
including outputs changed by earlier groups.

One paid binary source is retained across groups. Independent rejection
flags and a nonwrapping history make the accepted block the ordered
product of the group kernels. The growing dictionary cardinality is
included in the coefficient scale and cancels explicitly in the error
bound. R72's one-clean padding relocation and full recovery word return
the entire protected bank. The catalogue arrays and programs stay dirty
outside that bank and return exactly.

The middle stream has only $`O(1+L/\log\ell)`$ long-prefix rounds.
The early decoder has its separately charged logarithmic number of rounds. All table
volumes sum geometrically, and every native factor and catalogue slot is
charged. The protected bottom has $`O(\log\ell)`$ layers.
Complete original-precision fallbacks retain the literal widths on
narrow and finite-dimensional branches; the complex construction pays
its real precision shift and leaf diagonal.

## 2. Two independent retained results

[R91](PARALLEL_BILINEAR_WORD_QUERY.md) gives an exact q-bit word
query on $`Q=2^a`$ rows, using a separate arbitrary dirty pool
$`B\ge C(a+2)`$ and no initialized work:

```math
T=O(\sqrt{Qq}+Qq/B),\qquad G=O(Qq),\qquad
D_T=O(\log_2(a+2)+Qq/B^2).
```

A relative-phase polynomial uses only private dirty T targets. Taking
its finite difference under a Clifford row update leaves the required
bilinear phase and returns those targets exactly. Output batching,
controlled echoes and the retained indicator traversal give the full
width tradeoff. The previous additive q count/depth terms disappear.
This does not itself interpret a noncommuting dirty program.

A further R91 corollary pairs two XOR array translations. Three sparse
bilinear shears and a fixed swap execute both translations in logarithmic
program-length depth, at quadratic catalogue-volume count and dirty
width. The auxiliary array is permuted at that intermediate boundary;
it returns exactly in the full selected-slot conjugation. Combined with
the word query, this removes both additive program-length depth terms
when the larger simultaneous allocation fits. R93 keeps its original
smaller allocation and does not need this corollary.

[R92](UNIFORM_CATALOGUE_PRECISION.md) makes R90's precision ledger
absolute. For all $`n\ge1,L\ge6`$ it gives, with the original
one-clean counts and widths,

```math
D_T=O\!\left(NL/b^2+
\min\{F_{79}(n,L),\ell+2^{C_*(L+1)^3}\}\right).
```

An absolute window $`6\le L\le c[\log_2(\ell+2)]^{1/3}`$
therefore retains logarithmic depth. Three coarse stages suffice after
an absolute precision-dependent stopping cutoff; the full catalogue,
native factors, simultaneous pools and complementary fallbacks are
priced uniformly. This slow window does not reach $`L=\Theta(\ell)`$.
R92 and R93 use the retained word-query ledger and do not depend on R91.

## 3. Scoped failures and the next mathematical gate

Factoring a selected operator across separate dirty banks is valid:
the unwanted residues remain on those banks. If a second commutator
factor is fixed on the logical target, it can cancel those residues.
That identity did not give a faster precision stream: its selected
swaps and target rotations remained serial, and naive digit synthesis
retained a quadratic precision bill.

A dirty program also cannot simply be treated as a clean program. For
$`F(r_1,r_2)=Z^{r_2}X^{r_1}`$,

```math
F(r\mathbin\oplus f)F(r)^\dagger
=(-1)^{f_1r_2}Z^{f_2}X^{f_1}.
```

The dirty-dependent relative phase survives program return. An
unstructured pair of full terminal endpoint addresses similarly expands
the table volume beyond the linear Clifford budget. These statements
reject those interfaces; neither is an unrestricted lower bound.

A fixed parity involution also fails to factor a generic group cheaply.
Each real layer obeys $`J U_jJ=U_j^\dagger`$ for
$`J=Z^{\otimes h}`$, but their noncommuting product generally does not
obey $`JWJ=W^\dagger`$, which every commutator with this fixed J
would require. Half-angle factors produce a palindromic word instead.
Already rotations in the $`(00,10)`$ and $`(00,01)`$ planes give
a four-mode counterexample. A frame-dependent spectral factor remains
possible; its shallow native execution is the unpaid step.

The new benchmark's remaining cost is the number of addressed groups.
Increasing g requires more than changing one exponent in generic local
synthesis. The generic catalogue factors retain exponential local depth even
after R91's larger-allocation router is available. The residual also uses two
initialized g-bit endpoint labels; enlarging the common protected bank
would increase the retained bottom compiler's depth. A next mechanism
must price these costs together, or avoid that interface. A variable
local clean suffix is a possible representation to investigate, but
needs its own accepted-projector, failed-branch and return proof.

The resumption question is concrete: **can a complete large group be
executed and corrected with fewer long-prefix rounds while keeping its
native factor, routing, conditional work and bottom costs within the
same simultaneous budgets?** The general objective remains active.

## 4. Verification and integration record

Independent reviews reconstructed R91's full-input native phase identity
and allocation, R92's absolute cutoffs and original-width splices, and
R93's local decoder, changing-address kernel, normalization, separate
failures, full one-clean recovery and geometric resource sums. No
mathematical blocker was found. These internal reviews do not replace
external technical review.

The seven [reader tests](../../tests/test_parallel_bilinear_word_query.py)
and four [grouped-kernel tests](../../tests/test_grouped_catalogue_precision.py)
pass. They check bounded native identities and fragile operator
boundaries, including negative controls. They do not prove asymptotic
claims or provide a scalable emitter. R92 changes an analytic ledger
for existing words and needs no new operator fixture.

The complete deterministic suite passed **995 tests in 410.373 seconds**
with Python 3.12.14 and NumPy 2.3.5. The 61 focused operator,
documentation, mathematics and presentation-policy checks passed.
Python compilation, the reviewer walkthrough, both retained native
coarse examples, the residual-QBP example at q=16, offline synchronization
and both resource ledgers at n=12 passed. All four immutable
fault-tolerant receipt suites match without regeneration; their audit
scope and the resource ledgers' status fields are unchanged.

The final document gate covers **168 Markdown pages and five diagrams**,
MathJax SVG, native MathML, desktop and narrow layouts. Fresh browser
batches keep the full check within memory. One long inline inequality in
R92 overflowed the narrow layout and was moved unchanged into a display
equation. The corrected proof and final index/session pages are included
in the final rendered check.

Integration requires the analytic Python 3.11/3.13 and
rendered-presentation CI gates on the reviewed tree. The pull request
records that exact tree, gate receipts and merge receipt. No final
manuscript claim selection, efficient preprocessing, practical crossover
or scalable emitter is asserted.
