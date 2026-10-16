# Catalogue bootstrap session: the local compiler carries its workspace

[Current frontier](../../docs/OPEN_PROBLEM.md#large-width-t-depth) ·
[Executed design](CATALOGUE_BOOTSTRAP_PROMPT.md) ·
[Proof](PORTABLE_WORKSPACE_CATALOGUE.md) ·
[Previous session](PRECISION_TERMINAL_SESSION_RESULT.md)

This session executes the user's “Now go!” on 16 October 2026 from
main commit `d3dffa6e7f822d6bb68242b284560386ce3a0bf9`, containing the
design merged in [PR #157](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/157).
Execution began at **07:56:34 UTC**; the five-hour ceiling is
**12:56:34 UTC**, including review, validation, and integration.
The scientific baseline is R91–R93. Final manuscript drafting remains a
separate task.

## 1. Result and scope

The main result is [R94](PORTABLE_WORKSPACE_CATALOGUE.md). Put
$`N=2^n`$ and $`\ell=\log_2(n+2)`$. Uniformly for
$`6\le L\le2\ell`$, both prescribed frame families admit one
complete native circuit with

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),
```

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right).
```

The same circuit uses one external initialized qubit and the original
real dirty threshold $`b\ge17(L+n+7)`$ or complex magnitude threshold
$`b\ge17(L+n+8)`$. Its full initialized-isometry error is at most
$`2^{-L}`$, including every logical column, literal phases, signs,
zeros, clean leakage, dirty inputs, and reference correlations. All
constants are absolute on this precision range.

At inverse-polynomial error in the logical qubit count n and sufficient
square-root dirty width, the depth improves from R93's
$`O(\ell^2/\log_2(\ell+2))`$ to
$`O(\ell\log_2(\ell+2))`$. This reaches the design's primary
benchmark. It does not establish optimal general depth or the stronger
additive $`O(\ell)`$ target at growing precision. The selected
manuscript Result D remains R79; R90 and R92 retain their distinct
fixed-accuracy and precision-envelope advantages.

## 2. The construction that closes the gap

The retained opening leaves only $`O(L+\log\ell)=O(\ell)`$
logical layers: one uniform R89 reduction is followed by R93's
prefix-selected early decoder. Divide the middle into groups of at most
$`\lfloor\ell/8\rfloor`$ bits, leaving a fixed three-bit bottom.
Beyond the common absolute cutoff, there are only constantly many groups.

For each rounded group W, its controlled operator
$`A=\mathrm{diag}(I,W)`$ is itself a complete real Hopf frame.
The established R79 theorem compiles A at final accuracy with one local
initialized flag and an arbitrary local dirty pool. Route **every** local
compiler wire with the selected payload. The actual transfer acts as the
chosen native word on that payload and its actual inverse on one dirty
catalogue slot.

A commutator with a logical selector flip cancels the unwanted slot
action exactly. Its logical ideal action is
$`\mathrm{diag}(W^{\dagger},W)`$, so selector one executes W.
Only the logical flag requires initialization. Slot flags may be occupied
or entangled; the cancellation is an exact full-space identity. The
actual local compiler's initialized forward and inverse errors give at
most twice its error in this commutator, with no catalogue-size factor.

The finite final-accuracy code has

```math
q_h=O(L2^h)=O(\ell(n+2)^{1/8})=o(n).
```

The complete catalogue, all its native entries, portable local pools,
and R91's quadratic-volume router still cost $`2^{o(n)}`$ gates and
workspace. This fits below $`\sqrt N/n`$ after one absolute cutoff.
Local gate counts add across slots, whereas their T-depth is the maximum
because the full slots are disjoint. R91's parallel word query and
complete paired-router swap remove the otherwise fatal serial code-width
depth. The long global prefix remains in every relevant query volume.

Thus the existing weaker compiler supplies the local words needed by the
stronger global construction. The dependency graph is acyclic: R79 never
calls R94, and the new middle does not invoke its own theorem.

## 3. Work return and errors

The bottom logical bits supply distinct c,d,G positions only on a saved
zero-suffix sector. H is the sole external initialized flag. Each group
guard tests lower logical bits outside the changing payload, and each
group freshly reads its current prefix. Earlier group outputs are not
replaced by frozen cached selectors.

The guard conjugation restores its original computational value exactly,
even when c,d or the local dirty pool leak. On the active sector its
initial value is zero. On the inactive sector the actual zero-code word
is exactly identity on arbitrary occupied work. The proof does not infer
a guard-independent tensor factor on every arbitrary active guard input.

Actual local c,d and dirty return are approximate in the full norm;
catalogue, routing, program, and query return are exact at their stated
boundaries. The proof telescopes against ideal initialized intermediate
states. Later actual unitaries propagate all earlier leakage, including
through the final saved-predicate inverse, without reset or projection.

Final local precision
$`\kappa=L+\lceil\log_2(R+1)\rceil+10`$ bounds the whole middle
error by $`3\eta/1024`$. The retained opening and constant-height
bottom give total wide real error below $`507\eta/2048\lt\eta/4`$.
The internal real range extends through $`L\le3\ell`$ to pay the
complex L+1 precision shift and literal diagonal. Narrow and finite
branches use the complete original-precision R79 compiler at its original
family-specific width, avoiding a stricter tolerance at an unpaid width.

## 4. Approaches that were not silently reused

The new middle compiles directly to final accuracy. Its approximate
enlarged payload operation cannot be inserted into R93's old X-only
residual dictionary: it can leak into its selector and initialized flag
and act on its local dirty pool. No endpoint labels, shared residual
source, failure history, or residual recovery is used in R94.

Partial routing also does not inherit the complete tensor-cancellation
identity. A bounded counterexample uses a local controlled swap that is
identity on its initialized flag sector. Swapping only part of its payload
leaves an operation on a catalogue dirty bit for an occupied input. This
refutes that proposed exact-return interface; it is not a lower bound
against every partially routed construction.

## 5. Independent review and verification

Separate constructive lanes derived the local transfer/native ledger and
the global composition. A skeptical review reconstructed both, and a
fresh independent review checked the complete theorem against its retained
premises. The reviews identified and resolved the distinction between
guard-value restoration and arbitrary-input tensor independence. The
common cutoff explicitly includes the large-ell grouping condition.

The [six bounded fixtures](../../tests/test_portable_workspace_catalogue.py)
check complete-payload transfer, actual inverses and literal phase,
coherent selection of noncommuting entries, controlled Hopf marker
completion, initialized leakage, partial-routing residue, changed
selectors, and the occupied-sector guard identity. A negative guard reads
a changing payload bit and fails to uncompute. Complete operators have
at most 256 rows and columns. The retained R91 tests separately cover
dirty-program echoes and the paired router.

All 23 targeted tests across the new suite and the retained catalogue,
parallel-query, and grouped-catalogue suites passed. These are bounded
operator checks, not an emitted general compiler or a proof of asymptotic
resource bounds. Uniform cutoffs, geometric resource sums, the full error
hybrid, and literal-width/complex splices are analytic arguments. Internal
review is not external technical review or a priority claim.

Repository-wide validation passed all **1,001 tests** on Python 3.12.14
and NumPy 2.3.5. The 50 narrative, literature-policy, protected-math,
presentation, and typography guards also passed. All 12 changed Markdown
pages passed desktop and narrow rendering with MathJax and native MathML;
the repository-wide scan checked 10,338 distinct inline expressions and
17,086 protected occurrences with no TeX errors. No presentation failures
were found. Synchronization metadata passed its offline check, and
`git diff --check` passed. No immutable receipt was regenerated.

The integrating pull request records the final required CI checks and
reviewed commit/tree identities. Integration uses that reviewed tree and
requires every check to pass before merge.

## 6. Continuing destination

The intended benchmark is reached. The present construction still pays
$`O(\ell\log\ell)`$ through the retained opening's search rounds
and the local R79 compilation depth. Those are costs of this construction,
not unrestricted lower bounds. Optimal general depth and a broader
growing-precision improvement remain research questions. No practical
crossover, efficient classical preprocessing, or scalable native emitter
is asserted. The existing complete-frame count theorem remains settled.
