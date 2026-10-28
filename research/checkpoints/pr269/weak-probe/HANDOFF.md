# Weak complete-probe domain checkpoint after PR269

This is an independently audited continuation checkpoint, not a claim or
qualification receipt from PR269. It is intended for an unmerged continuation
branch after PR269. The frozen PR269 files were not changed during this work.

## What is established in this checkpoint

For an m=2^p target bank and q other arbitrary physical qubits, the contract
is only E Q_x (E^dagger Z_0 E) Q_x^dagger E^dagger approximately Z_x on the
whole space. A full wire-permutation endpoint is not required.

The proof file `WEAK_PROBE_PAULI_DOMAIN.md` establishes:

- Degree-d Clifford queries: below probe error 2^(-d-1), every image-domain
  Pauli X mask is in RM(d,p). Domain dimension is at most m+M_d(p)+2q,
  symplectic rank at most 2M_d(p)+2q, and native absorbed-axis rank at least
  m-M_d(p). No query-representation law is required for this first result.
- Exact additive degree-d queries: the pure-bank Z kernel is also in
  RM(d,p), improving domain dimension to 2M_d(p)+2q and the axis bound to
  2m-2M_d(p).
- Full-rank additive field queries, p>=2: at exact error, or at error below
  1/(3m), the pure-bank Z kernel is at most the all-Z direction. Domain
  dimension is at most p+2+2q and absorbed-axis rank at least 2m-p-2.
- Independent copies of the retained four-wire six-T cell, with arbitrary
  joint input/output Clifford walls and identity spectators in E, have weak
  complete-probe error at least 1/12 for every p>=3. The earlier symplectic
  rank argument improves that floor to 1/4 for p>=5.

The one-twelfth bound uses a weight-w kernel lemma and binary radius-two
Hamming-ball packing. It does not assume approximate invariance of the Pauli
domain. Literal cube phases, including negative scalar products, are retained.

Here q counts coordinates outside the target bank. They may participate in
the source query or move in a decoded endpoint. An identity spectator factor
in the particular candidate E does not impose identity query action there.

These are restrictions on domains and on a specified finite-cell construction
class. They do not exclude coupled non-Clifford bands, other non-Clifford
helper factors, or unrestricted weak-probe constructions. They are not growing
native-depth bounds and do not change R96.

## Validation and dependency

The portable exact fixture is `weak_probe_domain_audit.py`; it locates the
nearest enclosing checkout containing `WORKSPACE.md`. It depends on the
PR269 domain fixture and on retained exact native Pauli helpers. When this
directory is placed at `research/checkpoints/pr269/weak-probe`, run:

```bash
python research/checkpoints/pr269/weak-probe/weak_probe_domain_audit.py -v
```

The scratch fallback may be overridden with `HOPF_FRAME_REPO`. The exact
source hashes and test receipt are in `audit-receipt.json`.
`WEAK_PROBE_HOSTILE_AUDIT.md` is a separate independent audit of the kernel,
low-weight, and uniform-error arguments.

Checks include actual p=2 native-decoder saturation of the sharpened bound,
source field cube phases and mutual commutation, all p=3 Boolean masks for
the Reed--Muller derivative condition, helper/radical dimensions, joint
Clifford wall invariance, and the binary code bound. These focused fixtures
and independent audits do not replace a future repository qualification gate.

## Next-cycle use

Use this screen when evaluating proposals that extend the four-wire cell by
tensoring cells and adding only outer Clifford mixing. Weak-probe commutant
freedom does not rescue that class. The positive eight-wire native decoder
from PR269 uses coupled non-Clifford operations and is consistent with the
new result. A useful next construction must change the non-Clifford coupling
or its Pauli domain, while still paying the complete native probe and return
contract.

Before promoting this checkpoint into a permanent proof home, obtain the
usual repository review and exact-tree checks, then consolidate links and
tests without duplicating the proof. Keep the general degree-d, exact field,
width-dependent robust, and independent-cell uniform statements distinct.
