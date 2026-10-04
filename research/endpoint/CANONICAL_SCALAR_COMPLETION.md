# Canonical scalar completion in a grouped compiler

**Status:** valid compatibility construction; the repeated precision charge remains. This is an endpoint attempt, not a premise of Result C. Earlier section and equation numbers refer to the [grouped proof](../../docs/CONDITIONAL_SUFFIX_COMPILER.md); original section numbering is retained below.

## 11. A canonical scalar fits the group interface but retains its precision charge

This bounded compatibility test replaces the scalar completion in
Sections 4–5 by a real rotation, while retaining the selected atoms,
literal phases, private label, and actual coarse circuit. It works with
the allowed two external clean qubits and the same asymptotic dirty
budget. Its source cost still grows with the number of groups; it does
not improve the grouped bound or extend the one-clean result of Section 10.

### One scalar program, including literal inverse branches

Keep the external scalar flag $`\sigma`$ distinct from the atom flag.
For each exact coefficient $`c=c_{x,u,\ell}\in[0,1/4]`$, define

```math
R_c=\exp[-i\arccos(c)Y_\sigma],\qquad
\langle0_\sigma|R_c|0_\sigma\rangle=c.
```

Replacing $`\mathcal S_c`$ by $`R_c`$ in a forward atom gives the
exact accepted block $`\mathrm{diag}(c)U_\nu\Pi_\nu`$. A reverse
atom uses $`R_c^\dagger`$ before $`\mathcal D_\nu^\dagger`$,
with its independently specified coefficient row. Its block is the
adjoint of that real-base block. The literal $`\omega_\ell`$ remains
outside the base word in both directions. The scalar and atom flags
remain separate, so these projected products are valid.

The two direction-conditioned scalar stages can be combined. In
chronological order, apply the forward-controlled $`\mathcal D_\nu`$,
then one coefficient-table rotation, then the reverse-controlled
$`\mathcal D_\nu^\dagger`$. Surround that rotation by
reverse-direction-controlled $`Z_\sigma`$ gates, using the decoded
direction bit already supplied by the label. This uses
$`Z_\sigma R_cZ_\sigma=R_c^\dagger`$. The table address includes
the label and the current local word; neither changes during the scalar
stage. Thus both its query and inverse finish before the following atom
changes the local word. Diagonal and padding labels select the positive
rotation. An active zero coefficient uses $`R_0`$, not the identity;
only the inactive h or identity-mode sector is identity.

This consolidation also applies to the **old** scalar source. Its full
word has the form

```math
\mathcal S_c=\widehat cI+X_\sigma D_c,\qquad
D_c^\dagger=-D_c,\qquad
Z_\sigma\mathcal S_cZ_\sigma=\mathcal S_c^\dagger.
```

Consequently the old SELECT also needs only one scalar program, with
the same two direction-controlled Z gates and its exact inverse branch.
Its rounded accepted coefficients and error proof are unchanged. The
factor-two consolidation is therefore shared by both constructions.

### Native completion and the complete group error

Use the [middle-only paired-source rotation](ENDPOINT_TREE_TRANSPORT.md#10-a-shared-source-body-for-changing-targets)
with logical target $`P=Y_\sigma`$ and a separate arbitrary borrowed
signal a. The address and predicate exclude $`\sigma`$ and a.
For the positive-angle table, program the certified moments directly:

```math
s\approx\beta c,\qquad
p\approx\beta\sqrt{1-c^2},\qquad \beta=\sin(\pi/10).
```

The retained certified digit rule applies to these computable real
numbers; no separate numerical arccos approximation is required. Let
$`\widehat V_c`$ denote its five-call amplified word, including the
borrowed-signal correction $`C_a(Z_\sigma)`$ at both ends. The
full-signal SU(2) argument in the linked section proves

```math
\|\widehat V_c-R_c\otimes I_{\rm borrowed}\|
\lt\delta_q:=30\,2^{-q}.
```

This is a full-operator estimate, including arbitrary a, core, selectors,
and reference inputs. The h/mode-inactive word is exactly identity.
Moreover, the actual rounded word retains the inverse symmetry. With
its common fixed scalar A and middle word B,

```math
[A,Z_\sigma]=0,\qquad Z_\sigma BZ_\sigma=B^\dagger,
\qquad Z_\sigma\widehat V_cZ_\sigma=\widehat V_c^\dagger.
```

The last identity follows from the palindromic five-call amplification
word and its literal reflections. Thus the direction-controlled Z gates
select the actual native inverse, without a second independently rounded
program. Every outer inverse below is also the actual circuit inverse.

Let $`Q_{\rm can}`$ be the ideal half-block with these exact canonical
rotations and $`\widehat Q`$ its native replacement. There is one
rotation appearance in Q. All other operations are identical unitaries,
so the full-operator difference and active accepted block satisfy

```math
\|\widehat Q-Q_{\rm can}\|\leq\delta_q,\qquad
J^\dagger Q_{\rm can}J
=\tfrac12(C^\dagger W_g)\otimes I_{\rm dirty}.
```

The second identity is on $`h=1`$; on $`h=0`$ both half-blocks
are exactly identity.

The ideal normalization-two amplification is exact on initialized
columns. Telescoping its three Q appearances, with the same $`R_h`$
and literal $`Z_h`$ as Section 5, gives

```math
\|\widehat{\mathcal A}_hJ-
J(C^\dagger W_g\otimes I_{\rm dirty})\|
\lt3\delta_q=90\,2^{-q}.
```

Choose $`q_g=m_g+4`$. Then this is
$`(45/8)2^{-m_g}\lt10\,2^{-m_g}`$, within the existing group
allocation. Applying the actual coarse word and uncomputing h preserves
the bound. Full-isometry telescoping includes all prior leakage; there
is no intermediate reset. This proof compares with the exact canonical
half-block, rather than asserting that the approximate native rotation
has a dirty-independent scalar accepted entry.

### Registers, reflections, and common-source boundaries

Use the original two-clean layout, with a protected scalar flag:

| Register | Reservation and use |
| --- | --- |
| External h | Initialized suffix-predicate flag; exact computation and actual inverse, with final return included in the group error. |
| External $`\sigma`$ | Initialized scalar target and rejection flag; never coarse or predicate scratch. |
| Private active suffix | Mode, term label, distinct atom flag, and the existing temporary work. |
| Paired core | $`q_g+1=m_g+5`$ arbitrary external dirty wires. |
| Selectors and helpers | Existing $`k_g+O(1)`$ returned arbitrary dirty reservation. |
| Signal a | One further arbitrary external dirty wire, excluded from $`R_h`$. |

Relative to the old m-bit core, the new core and signal require six
additional dirty slots; the remaining helper constants are fixed.
Increasing the fixed sufficient thresholds $`C_1,r_0`$ in Section 6
absorbs this constant in

```math
(m_g+5)+k_g+1+O(1)\leq L+n+7.
```

No private initialized suffix wire is counted again as dirty space.
The bounded h/mode source-center predicates use the retained returned
helper; direction controls are Clifford CZ gates on the decoded bit.
Selected atoms, larger predicates, reflection, and coarse interpretation
retain their separately charged native costs.

All remaining operations can be chosen to commute with $`Z_\sigma`$.
In particular, $`R_h`$ tests $`\sigma`$ in the computational basis
but does not test a or the dirty core. Completed non-scalar atom,
predicate, query, and coarse words act identically on their returned
source-core and signal helpers, so they commute with A as complete
operators. Programmable mask queries inside the scalar are not included
in this claim. Individual gates
that borrow those helpers need not commute with A; this observation does
not supply free native controls.

It follows that all borrowed-signal corrections may be moved to two
outer $`C_a(Z_\sigma)`$ boundaries. The ideal uncorrected program
has signal branches $`W_{\rm can}`$ and
$`Z_\sigma W_{\rm can}Z_\sigma`$, in the same chronological order.
The boundaries restore the arbitrary signal, with the same summed
operator error. At a fixed q and common physical core, the intervening
completed operations also permit the exact common-A cancellations of
the linked source-body identity. Unequal widths require a charged
transition; this argument does not supply one.

### The fair variable-group ledger and stopping point

There are three inner amplified rotations per group, including the
actual inverse occurrence in the outer amplification. Hypothetically
placing R groups at one legal common q gives $`g=3R`$ inner stages.
For the paired canonical route, the declared source-appearance count is

```math
135R\quad\longrightarrow\quad81R+6.
```

For one group this is 135 versus 87. However, simplifying the fixed
mask in both words, as in the linked section, gives

```math
T_{\rm loader+fixed}\leq
\begin{cases}
(120R+4)q+8(12R+2),&\text{shared paired body},\\
(120R+4)q+240R,&\text{expanded paired baseline}.
\end{cases}
```

Their leading precision coefficient is the same. Native programmable
queries, predicate controls, and fixed-mask Cliffords remain charged.
These are bounds for displayed words, not synthesis minima. Nor is the
expanded paired word the old grouped compiler:

| Quantity per group after the direction-Z consolidation | Old one-tail scalar | Canonical paired rotation |
| --- | --- | --- |
| Scalar programs per Q | 1 | 1 |
| Scalar or inner-rotation appearances after outer amplification | 3 | 3 |
| Programmable mask appearances, including actual inverses | 6 | 30 |
| Certified group error | $`10\,2^{-m_g}`$ | $`(45/8)2^{-m_g}`$, with $`q_g=m_g+4`$ |

The old scalar already accepts an arbitrary coefficient with its dirty
core returned on the accepted block. Canonicalization adds an inner
rotation compiler to gain a compatible full-operator completion; it
does not remove a precision cost that the old grouped word had hidden.
For the actual unequal-width schedule the legal construction has

```math
T_g=O(s_g2^{e_g}+q_g+\mathrm{poly}(n)),\qquad
G_g=O(s_g2^{e_g}q_g+q_g+\mathrm{poly}(n)).
```

Thus the sums remain $`T=O(N+L\ell_*(n))`$ and $`G=O(NL)`$.
Even granting free equal-width exterior cancellation, the displayed
programmable interiors still charge precision in proportion to R.
The compatibility test therefore closes this canonical-completion
candidate as a mechanism for removing the iterated-logarithm overhead.
It is a scoped cost conclusion for these words, not a lower bound on
other joint compilers.

The [grouped scalar checks](../../tests/test_grouped_scalar_completion.py)
exercise canonical atoms, literal native inverse branches, borrowed
signal inputs, outer amplification, and the source-appearance ledger.
Their $`q=2`$ words are finite diagnostics, not certification of the
stated L-bit resource allocation.
