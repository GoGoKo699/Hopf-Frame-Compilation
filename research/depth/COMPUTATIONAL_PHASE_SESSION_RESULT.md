# Revision: the computational phase and the complete native source

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Computational derivatives](COMPUTATIONAL_PHASE_DERIVATIVES.md) ·
[Clean-phase interface](RANK_ONE_CLEAN_PHASE_INTERFACE.md) ·
[Majorana source audit](MAJORANA_WALSH_QUERY.md) ·
[Previous checkpoint](CLEAN_TRANSITION_PHASE_SESSION_RESULT.md)

This continuation began from verified main `5b8e063` (PR #190) on
19 October 2026. It moved from growing representation parameters to
derivatives of the actual target phase, and from a formal Walsh phase
program to a fully native interference source. Neither change yet
improves the complete-frame resource frontier.

## 1. Destination, strongest theorems and component results

The destination is optimal general depth for the prescribed complete
Hopf differential frame. The same circuit must meet the count, depth,
width, Clifford and full initialized-isometry error bounds.

| Role | Current evidence |
| --- | --- |
| Strongest general upper theorem | R96 gives $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, and $`G=O(NL)`$ |
| Literal domain | Integers $`n\ge1,L\ge6`$; one clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; error at most $`2^{-L}`$ |
| Selected manuscript theorem | R79 remains Result D; the count package remains settled under its stated allocations |
| Unrestricted low-layer components | At most one clean and two layers imply XOR degree at most two; zero clean and three layers imply sign degree at most five and XOR degree at most four |
| Fixed-rank closure | The retained finite-insertion theorem handles both fixed endpoint clean-active ranks, including arbitrary middle pairs |
| New growing-rank subcases | A shared clean flip or the actual-source noncentral common-parity case gives sign/XOR bounds five/four; computational derivative coverage gives a quartic part plus a residual quotient under the stated hypotheses |
| New structural progress | The coupled C-side union is independent; Boolean reality and an actual Clifford correction put its active center in computational coordinates |
| Still open | Uniform one-clean three-layer degree at growing ranks, the faster complete R98 probe, and matching general depth bounds |

At fixed accuracy and square-root dirty width, the general upper bound
is still logarithmic in n and the unrestricted lower bound is constant.
A conditional degree theorem or a source-specific obstruction does not
close that interval.

## 2. What the four retained phase factors now establish

The no-middle-pair reduction retains its original logical relation
profile $`(0,\le1,\le1,0)`$. The two output lists together have at
most one relation, and are jointly independent whenever their rank-one
coupling is nonzero. Their canonical form has one quantum core and
disjoint independent central coordinates. This is stronger than merely
counting the core, but does not bound the completed phase by itself.

There is a sharper closed branch. If a clean-odd Pauli in the common
active span commutes both silent lists, it commutes all three original
bands. Actual occupied-fixing controlled-Pauli Cliffords then turn the
full operator into $`I_c\otimes G`$ at unchanged native depth.
The old clean wire can be arbitrary dirty work, so the zero-clean
three-layer theorem gives sign degree five and XOR degree four.
The criterion is an exact linear test on the supplied native split;
its failure is not an obstruction for every possible representation.

## 3. Computational derivatives replace an intrinsic-coordinate guess

An actual C-active clean flip yields a full-space block identity
$`V'=\mathrm{diag}(G,GL)`$, with L an actual logical Clifford.
In the zero-character determinant branch, excluding C-active rank two,
the retained endpoint realizes every resulting target derivative by
two native bands. If U is the image of their computational X masks,
projected onto the n designated sign variables, then

```math
g=q_4+h(\pi x),\qquad \deg q_4\le4,
\qquad \deg g\le\max\{4,n-\dim U\}.
```

This closes bounded residual codimension even when the intrinsic
determinant codimension grows. More clean flips do not automatically
improve it: U is invariant under admissible normalizations of the same
full operator. The actual disjoint-CCZ family has U=0 at arbitrary
arity and satisfies the hypotheses. Its degree is only three.

A separate four-factor argument needs no determinant or zero-character
assumption. Combine the commuting middle factors before probing.
Paulis centralizing either outer native span give actual two-band
derivatives. All degree-above-four dependence is consequently confined
to the computational parity algebra in the intersection of the two
transformed endpoint spans. Both statements concern the actual target
phase; neither deletes a quartic correction at the original depth.

The stronger normalization uses Boolean reality and native independence.
Every first-band axis activated by an output active-center Pauli would
have to be input diagonal. Squaring the resulting real derivative
forces the sum of those axes to vanish; independence makes the set
empty. Thus the whole center already transports as Paulis and its
target derivatives are affine. Extract an actual diagonal Clifford,
then make the common center computational Z using Hadamards only on
coordinates where the residual target is identity. This closes the
input-coordinate step without removing an arbitrary cubic correction.

Applying this argument at both endpoints is stronger for the actual
clean transition. After a canonical aligning Clifford, its residual is
literally a Boolean cubic on the common center plus one Pauli direction.
The central cubic sign moves to the input, leaving one Pauli controlled
by a quadratic function of actual center coordinates between three
independent logical native bands. Every inverse and central phase is
retained; the insertion has not been assigned a free native layer.

A separate noncentral common-parity subcase closes with sign degree
at most five and XOR degree at most four for the actual clean-derived
no-pair circuit. The physical code stabilizer forces active-middle
independence, including the exceptional endpoint rank two. The physical
Gram identities then split the clean phase into two endpoint-supported
factors. Separate literal Clifford alignments leave only central cubic
signs and three independent logical native bands. The broader abstract
four-factor theorem retains its conditional bounds 2048/2047; its
one-relation exception cannot occur in this actual-source branch.

**The first remaining inference is now global-center matching.** Fixing
the shared center can create arbitrarily many relations in the native
lists. The retained identity witness makes all m middle axes coincide
on a center fiber, while its opposite factor cancels them exactly.
Applying the bounded-relation theorem to those fibers is therefore
unsupported. This is a failure of that representation argument, with
an easy target, not an unrestricted degree obstruction.

## 4. A native source audit and a change of mechanism

The Majorana Walsh construction supplies an exact independent-Pauli
probe on every carrier input, with actual inverses and no clean work.
For $`m=2^p`$, its complete ledger is

```math
T=2mp,\qquad G=O(m^2),\qquad D_T\le4p,
\qquad\text{carrier width}=m/2.
```

The formerly nonnative fine rotation has disappeared. The full count
and depth still miss the required source budgets. Within the specified
Majorana-preserving native model, matrix entropy and support growth
prove sharp source costs $`mp/2`$ T gates and p layers. Those arguments
do not cover arbitrary Clifford interlayers or joint synthesis of the
whole probe word.

A further fused route uses a single flat Majorana operator to obtain
an equality phase. Its linear-count source avoids the full Walsh
matrix count, but universal dirty-table readout needs a jointly
controlled dual frame. The supplied readout wrapper restores the
costly Walsh conjugator. Making the selected source cheap has not paid
that complete operation. The separate R75 argument gives only a
constant unrestricted obstruction for the flat operator, not a growing
lower bound or a new frame frontier.

## 5. Next research question and continuing contract

The next lower question is whether three globally independent native
lists with this single quadratic-controlled Pauli constrain the phase
over their shared computational center, retaining the correlated
cancellations before taking any fiber.
After the shared-flip and noncentral-parity closures, every parity in
the normalized common endpoint computational space must commute with
F in this remaining representation. In fact each has a common signed
physical lift commuting with every original band; the active middle
sublist is empty. The possible higher-degree dependence therefore lies
in a genuine shared computational center. Commutation alone does not
control the relative phases of its sectors. The next argument must
retain global native independence across that whole center.
A successful proof must control the whole center, or a counterexample
must give a complete growing-degree native family. The parallel upper
question is a jointly compiled complete probe or dual-frame readout
whose count, depth and width all improve together. Reordering the
Gaussian butterfly or compiling only one selected source retains an
old cost.

The new statements were independently audited. Exact finite source
fixtures retain every physical column, ordered address-mask phases,
actual inverses and native resource words. Existing native witnesses
also check the representation limits; the arbitrary-width theorems
remain analytic. Every logical column, dirty input, external reference,
literal phase, leakage/return condition and actual inverse remains in
scope. Final manuscript drafting remains separate.
