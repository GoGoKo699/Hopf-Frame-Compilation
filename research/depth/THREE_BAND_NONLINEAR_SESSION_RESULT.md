# Three native bands and a noncommutative return mechanism

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Native source structure](NATIVE_BAND_SOURCE_STRUCTURE.md) ·
[Bounded-deficiency transport](BOUNDED_DEFICIENCY_TRANSPORT.md) ·
[Cycle character](NONCOMMUTATIVE_CYCLE_CHARACTER.md) ·
[Previous checkpoint](NATIVE_TRANSPORT_SESSION_RESULT.md)

This continuation began on 18 October 2026 at verified main
`357040fb4f1163a013be2373af6f287243b5ce9e` (PR #175).
It addresses two gaps left by that checkpoint: what independent physical
T axes imply about the actual three-band source, and whether a nonlinear
returned character can escape the exponential coefficient-width cost.
Both yield supported results below. Neither improves the unconditional
complete-frame upper law or establishes a growing general lower bound.

## 1. Compact revision

The count package remains settled under its stated hypotheses. For
integers $`n\ge1,L\ge6`$, $`N=2^n`$, R96 gives one circuit with

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It uses one clean qubit, real dirty width $`b\ge17(L+n+7)`$ or
complex magnitude width $`b\ge17(L+n+8)`$, and full error at most
$`2^{-L}`$. Every logical column, literal sign and phase, singular and
zero coefficient, dirty input and external reference, actual inverse,
leakage and complete work return is included. All sources, queries,
controls, recoveries and simultaneous allocations remain charged.
R79 is still selected Result D. R97 is a scoped four-layer obstruction;
unrestricted one-clean four-layer rigidity remains open.

The exact destination is optimal worst-case general T-depth for that
prescribed complete operator, at explicit precision and simultaneous
count and workspace budgets. At fixed accuracy and square-root width,
the remaining gap includes whether the logarithmic-in-n additive depth
can be improved or is necessary. R98 transfers a faster exact multi-output
read-only query to a better upper bound; the separate two-call corner
reduction transfers unrestricted query lower bounds. Neither missing
premise follows from the scoped results here.

Formerly open questions already resolved include optimal one-clean count,
the all-precision R96 upper law, and zero-clean read-only five-layer
degree at most 14 below $`\sqrt2/24`$. The six-layer Pauli-comparison
bound 51 and fixed-diagonal-C3-source bound 78 remain distinct subcases.
This round additionally resolves uniform bounded-deficiency transport,
including degenerate Pauli domains and sources arbitrarily close to a
Clifford. Historical next-step passages are evidence of earlier gaps,
not instructions to reopen completed theorems.

## 2. Native axes now certify additional source cases

The [source proof](NATIVE_BAND_SOURCE_STRUCTURE.md) retains the specified
factorization $`F=E_CE_BE_A`$, with independent isotropic axes within
each physical band. It proves:

- A two-band source is C3 exactly when its two axis lists commute across
  the lists. It is then diagonal C3 in fixed Clifford coordinates.
- Adding a native band on either side of a C3 source, while staying in
  C3, forces its relevant axes to have Pauli images or preimages.
  Neither statement assumes that the inverse source is C3. In particular,
  a band sequence whose every partial source stays C3 remains in the
  fixed diagonal source class.
- For a three-band C3 source, let
  $`t=\dim(\mathrm{span}(A)\cap\mathrm{span}(B))`$.
  At most $`2t`$ individual B axes fail to commute with C. When t is
  zero, the whole source is diagonal C3 in fixed Clifford coordinates.
  The same conclusion holds if either adjacent two-band partial source
  is C3.

These certificates activate the existing six-layer bound **78**, or
**26** for visible outer target shear, below error $`\sqrt2/8`$.
They do not infer C3 membership from full query return. The separate
robust source certificates require a uniform promise on every source
Pauli image, which return does not supply. Axis intersections depend on
the native representation; they are not intrinsic lower-bound parameters.
Commuting factors in an analytic proof are not a free physical band.

## 3. A uniform invariant when the Pauli domain has bounded codimension

The [bounded-deficiency theorem](BOUNDED_DEFICIENCY_TRANSPORT.md) applies
without C3 membership. For the actual source F, let k be the codimension
of its exact signed Pauli-preserving domain. If its radical has dimension
r, fixed Clifford coordinates isolate an arbitrary source on
$`a=(k+r)/2\le k`$ qubits and an identity spectator. This includes
degenerate domains; no unproved controlled factorization is needed.

For $`k\ge1`$, put $`D(k)=(2k+1)(4k+6)`$. The actual native
six-layer comparison obeys

```math
\deg S_K\le D(k),\qquad
\deg\kappa_K(q)\le3D(k),\qquad
\eta\lt\frac1{4\sqrt2(2^k+1)}.
```

The output degree is at most $`3D(k)`$, or $`D(k)`$ for visible
outer target shear. The sharper a may replace k throughout. For k zero
the separate Clifford case gives label/sign bounds 2/6 below
$`\sqrt2/4`$.

The new step uses heavy Pauli coefficients and the source-independent
separation of the signed conjugated Pauli basis. It does not require $`FL_xF^\dagger`$ itself to be exactly Clifford,
or a positive gap specific to one source. The actual comparison $`K_x`$
is still exactly Clifford. Literal Pauli multiplication phases enter the reconstruction.
This closes the quantitative gap in the previous localization proposal.
But even one wide parallel T band can have unbounded k. Neither depth nor
the small-overlap theorem bounds k, so this is still a scoped lower bound.

## 4. Constructive escape, full circuit, and fair comparison

The [commutative theorem](NONLINEAR_CHARACTER_TRANSPORT.md) now covers
every finite-dimensional commutative nilpotent binary algebra, not just
one truncated polynomial variable. On w coefficient wires its affine-input
characters still have degree at most $`1+\lfloor\log_2w\rfloor`$.
Changing to several commuting variables cannot evade that width cost.

The [noncommutative cycle algebra](NONCOMMUTATIVE_CYCLE_CHARACTER.md)
does evade it. Its $`d^2`$ positive-path coefficients are all valid dirty
inputs. A logarithmic-derivative trace character satisfies

```math
\chi\!\left(1+\sum_{i=0}^{d-1}x_i e_{i,1}\right)
=\prod_{i=0}^{d-1}x_i.
```

Multiplication by that unit and the controlled character source give a
literal returned XOR through the actual commutator and inverse. An
in-place factorization supplies the source on every coefficient word;
there is no supplied eigenstate, clean coefficient register or omitted
recovery. A native schedule has

```math
T,G=O(d^3),\qquad D_T\le64(d-1),\qquad
b=d^2+d,\qquad a=0,\qquad\eta=0
```

for $`d\ge2`$. The proof separately records exact sequential counts
with one reused helper; its lexicographic factor order is not confused
with the parallel order. Shared-control phase terms are scheduled
explicitly, not declared disjoint.

This is a counterexample to a representation-independent logarithmic
character-degree bound, **not a faster AND or an arbitrary-table reader**.
The retained strictly read-only conjunction already has logarithmic
depth, $`O(d^{\log_2 3})`$ T-count and dirty width, and the stated
slightly larger Clifford count. If intermediate address use is permitted,
the imported two-dirty-helper MCX gives linear counts and logarithmic
depth. Both comparisons include actual return and native decompositions.

The cycle multiplication's actual inverse has a coordinate of Boolean
degree d, and the in-place factorization map F also has coordinate degree
at least d. Consequently affine circuits with disjoint Toffoli layers
need at least $`\lceil\log_2d\rceil`$ such layers to implement
$`P_x`$ or F. This does not lower-bound the diagonal controlled source
$`D_h`$ by its coordinate degree. This diagnoses
where ordinary reversible parallelization retains cost. It is not an
unrestricted Clifford+T depth bound: nonmonomial Clifford interference
is outside that argument.

A second candidate, natural multiplication on a padded finite-field
dual-number ring, fails literal scalar return because its zero word is
fixed. The trace obstruction persists for nonmonomial phase sources in
that representation. Other representations remain eligible.

## 5. Additional scoped result and evidence

The [permutation transport note](PERMUTATION_SOURCE_TRANSPORT.md) proves
a uniform monomial endpoint gap and degree bounds from the source's
forward Boolean degree h, with a monomial native inner family. The inverse
degree is not assumed small. No actual three-band non-C3 instance with
the required returned comparison was obtained: the natural overlapping
Toffoli candidate currently needs four bands. The theorem is retained as
a resumption result, not an enlargement of the general six-layer claim.

Independent analytic audits check the source extensions, overlap relation
count, degenerate-domain localization, uniform error estimate, cycle
character and all native ledgers. Focused fixtures check signs, dependence
counterexamples, approximate endpoints, full dirty/output columns, source
factorization and actual inverses. These bounded checks support the proof
mechanisms; arbitrary-width conclusions come from the proofs. The
[verification catalogue](../../docs/reference/VERIFICATION_CATALOGUE.md)
links each result to its evidence without changing manuscript scope.

## 6. Precise unresolved step and resumption question

The most useful next lower-bound question is **whether bounded overlap t
controls signed transport for the actual compatible three-band C3 source,
even when its diagonal part has unbounded width and Pauli deficiency**.
Start with t equal to one. The new theorem leaves at most two offending
middle axes, but their observables can still couple to arbitrarily many
commuting outer axes and to the inner band. Simply counting these axes
does not isolate a small non-Clifford subsystem. The bounded-deficiency
theorem therefore cannot yet be applied to the remaining source.

Can that coupling be removed analytically while retaining native inner
compatibility, literal signs and a width-independent error radius, or
does an actual fully returned family exhibit growing degree already at
fixed t? This is a concrete extension of the six-layer boundary, and
success would test an invariant capable of supporting a depth recurrence.
General non-C3 prefixes, one-clean work, address mixing and error growth
remain further gaps before a dimension-growing frame lower bound follows.

Constructively, the cycle mechanism removes the coefficient-width barrier
but leaves source and inverse depth. A useful next construction must beat
the logarithmic conjunction baseline or implement arbitrary multi-output
tables within the full R98 interface and R91 count/work budgets. Replacing
its reversible source by native interference is eligible; ignoring its
inverse or comparing only with an unsimplified table reader is not.
Optimal general depth remains the destination. Final manuscript drafting
remains separate.
