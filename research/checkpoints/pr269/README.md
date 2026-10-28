# Continuation after PR269

This unmerged branch preserves follow-up research from the PR269 validation
window. It is a checkpoint for the next review cycle, not a change to the
accepted resource frontier. The main proof homes and continuing workspace
remain authoritative. Focused checks and independent audits are recorded with
the individual files; this checkpoint has not passed the full repository and
presentation gates as a new integrated revision.

## Audited finite construction

The coupled projective word uses 148 T gates in at most 39 T bands on eight
principal wires and six additional arbitrary helpers, all returned. Its
recoupler costs 136 T gates in at most 36 bands. This is a separate allocation
from PR269's 172-T, 83-band construction with no additional helpers.

The improvement shares paired-reflection factors across projective lines and
batches the six quartic displacement terms through private dirty helpers.
The complete logical action, literal phases, actual inverses and helper
sector identities are covered by the linked code and receipts in `coupled/`.
The finite two-bit OR predicate does not supply a constant-depth growing-field
nonzero predicate.

## Weak complete-probe restrictions

The proof and exact fixtures in `weak-probe/` concern the weaker complete bank
probe, without demanding a full Clifford query endpoint. They include the
degree-sensitive Reed--Muller restriction on the image Pauli domain, a
stronger exact field bound, and a robust obstruction to independent four-wire
cells even with arbitrary joint outer Clifford walls and extra qubits.
The proof states the distinction between source identity spectators and
moving decoded helpers explicitly.

Its independent-cell error floor is one twelfth from three address bits
onward, improved to one quarter from five bits onward. The exact field-domain
bound and the width-dependent robust bound have separate hypotheses in the
proof. Seven portable exact fixtures and an independent hostile audit support
this checkpoint.

## Uniform augmented recoupling candidate

The candidate in `coupled/augmented/` replaces norm division with two batches
of multiplication-XOR operations and a fixed Fourier transform. One moving
Pauli helper register can be shared by all cells at a recursion level. Its
explicit cost is at most sixteen T bands per level, with polynomial overhead;
the private returned multiplication masks are counted separately from the
retained moving helper registers. The probe echo, including its actual
adjoint, returns both kinds of helper.

The bounded sector fixture checks every address, every helper character and
every logical column at two and three address bits. This is an analytic
continuation candidate awaiting integration review; its linear depth in the
number of address bits is not the desired depth improvement.

## Work backward from the destination

The sufficient source target remains a complete dirty indicator at fixed
tolerance, with zero clean inputs, linear T count and dirty width, quadratic
elementary-gate count, and depth sublogarithmic in the number of address bits.
Polynomial-overhead fusion permits discovery at polynomial cost on a smaller
address. Every helper and reference input remains arbitrary, and every return
error must be paid.

The next operation to investigate is joint synthesis across recursive
recoupling levels. A paid operation at each level still retains the level
count; reducing one finite cell does not settle that obligation. In an
augmented moving-helper formulation, the sector Cliffords satisfy a cocycle
equation. A generic solution using controlled wire translation can conceal
the same barrel-shift cost, so a useful solution needs an explicit full native
ledger for the jointly compiled word.

A concrete next experiment is to emit the helper-sector Clifford tableaux
symbolically, starting at four address bits, instead of constructing dense
unitary matrices. Track the dependence of each tableau entry and literal
phase on the helper characters through the full recursion. Use those data
to test whether the cubic parity layers can share one global Clifford frame
and a bounded number of native batches. Any proposed compression must then
be proved uniformly and checked with the actual adjoint and arbitrary-helper
return. Finite tableau patterns alone do not establish the required depth.

An additional [computational-helper-basis audit](coupled/augmented/AUGMENTED_COMPUTATIONAL_AUDIT.md)
checks all 2,080 principal/helper columns at two and three address bits,
16,512 address-column comparisons and every reversed algebraic-word column.
Its [manifest](coupled/augmented/COMPUTATIONAL_AUDIT_MANIFEST.json) is separate
from the earlier sector-fixture receipt. These are complete checks of the
bounded algebraic identity; the retained multiplication compiler owns the
nonlinear native-gate ledger.

Before integrating this checkpoint, independently review any uniform augmented
recoupling identity, move accepted claims to their single proof homes, add
meaningful fixtures and routing entries, and qualify the final integrated
tree. R96 and selected Result D remain unchanged.
