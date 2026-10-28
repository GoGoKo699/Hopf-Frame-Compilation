# Coupled continuation after PR269

This directory is an **unmerged research checkpoint**, kept separate from the
PR269 review tree. It contains two distinct results and their evidence; it
has not passed the repository's complete qualification gate. R96 is unchanged.
The final zero-extra-helper PR269 decoder remains a separate retained allocation.

## Finite coupled decoder

The common A-reflection is selected once by the two-bit OR of the upper U bits.
Relative phases cancel literally across that shared selection and the four
remaining line-selected reflections. Six displacement terms are compiled in
parallel with the existing R91 private-mask phase and six arbitrary helpers.

| Construction | T | T bands | Forward gates | Actual inverse gates | Extra arbitrary helpers |
|---|---:|---:|---:|---:|---:|
| PR269 zero-extra-helper decoder |172|83|1,234|1,286|0|
| This coupled recoupler |136|36|980|984|6|
| This coupled decoder |148|39|1,110|1,166|6|

The coupled words use eight principal wires and six additional helpers, all
arbitrary. The private helpers return literally. The zero U sector has a new
specified Clifford action, which preserves complete intertwining. The
independent audit covers the signs, every logical column, actual inverse, all
addresses, and composition with the full arbitrary-helper return proof.

Files:

- [Emitter and logical verification](coupled_projective.py)
- [Construction receipt](coupled_receipt.json)
- [Private quartic batch](quartic/private_quartic_batch.py)
- [Quartic proof and scope](quartic/QUARTIC_BATCH_HANDOFF.md)
- [Independent coupled audit](quartic/COUPLED_INDEPENDENT_AUDIT.md)
- [Independent audit script](quartic/audit_coupled.py)

The finite common OR is not a uniform bounded-depth predicate. Its all-size
analogue is a nonzero test on ell=m/4 bits, of Boolean degree ell. The
exponential projective-line list and cross-level dependencies remain unpaid.

## Uniform augmented source interface

The separate [augmented recoupler](augmented/AUGMENTED_RECOUPLER_CANDIDATE.md)
adds an arbitrary shared helper per recursive level and retains an explicit
Pauli helper query on both sides of the intertwining identity. Two multiplication-
XOR batches replace each level's division-dependent coordinate changes. The
same helper query is passed through each child successively; helper sharing
uses exact commutation in their common X basis and requires no helper swap.

Its safe ledger is16(p-1) native bands, O(mp^3) T gates, O(mp^2) dirty width and
O(m^2p^2) elementary gates. Only the complete source/probe/actual-adjoint echo
returns the retained moving helper registers. This is an alternate uniform
compiler for a Clifford cocycle, **not** a sublogarithmic-depth source. A generic
cocycle solution may hide controlled wire translation; the displayed compiler
instead pays every recursive level explicitly. Compressing those coupled levels
remains the main next operation.

- [Independent analytic audit](augmented/AUGMENTED_RECOUPLER_AUDIT.md)
- [Exact shared-helper sector fixture](augmented/shared_helper_fixture.py)
- [Sector fixture receipt](augmented/shared_helper_receipt.json)

## Reproduce the focused checks

From the repository root, with its existing NumPy dependency:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/checkpoints/pr269/coupled/quartic/private_quartic_batch.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/checkpoints/pr269/coupled/coupled_projective.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/checkpoints/pr269/coupled/quartic/audit_coupled.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/checkpoints/pr269/coupled/augmented/shared_helper_fixture.py
```

The scripts discover the enclosing repository from their own location and write
only adjacent JSON receipts. Regenerable binary matrices are omitted. Their
imports use the existing PR269 native fixture, the retained recursive field
arithmetic and native words, and the separate exact CQ2 evaluator. The augmented
fixture checks all helper-X sectors for p=2,3; it is not a physical Clifford-only
implementation of a coherently controlled helper family. The retained
multiplication-XOR proof owns that family's nonlinear native ledger.
