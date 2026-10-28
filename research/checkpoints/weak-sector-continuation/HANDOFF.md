# Weak-sector continuation checkpoint

This is an unmerged follow-up to the frozen PR270 revision. The proof has
received independent mathematical audits and six focused exact fixtures
passed. It has not received a new complete repository gate, rendered browser
review, or merge qualification. No PR270 source file was changed.

## Results and scope

[WEAK_SECTOR_DEGREE.md](WEAK_SECTOR_DEGREE.md) proves that a degree-d Clifford
query and a linear helper shift require seed degree and helper-shift rank at
least p-d under the weak complete-bank-probe contract. For exact Clifford
sector decoders, operator error below the square root of two rounds the
sector Pauli observables exactly. The argument needs neither a full query
endpoint nor the square-zero field-pencil assumption.

For affine queries, the seed degree is at least p-1. In the augmented
construction it is exactly p-1 on every active helper coset, with an explicit
nonzero leading-coefficient certificate. A common-frame factorization into
helper-affine Clifford core stages therefore still requires p-1 stages.

An exact extension permits arbitrary sector unitaries when every seed
pullback is Pauli. Its separate robust extension has radius
$`\min\{\sqrt{2-\sqrt2},2/\sqrt m\}`$. Under the usual identity query at
address zero, the second term improves to $`2/\sqrt{m-1}`$. This extension
does not assert a uniform Pauli snap for non-Clifford sectors.

The final stage-count refinement charges both degree-two reversible helper
layers and helper-affine core stages. It is a bound within that specified
classical-helper/Clifford-core grammar, not unrestricted native T depth.
Quantum helper motion and non-Clifford intermediate sectors remain outside
that grammar. Every proposed new primitive still needs its complete native
cost, actual inverse, and arbitrary-input return proof.

## Evidence and reproduction

- [test_weak_sector_degree.py](test_weak_sector_degree.py) contains six exact
  tests, passing on Python 3.12.14; the receipt records the final run time.
- [qualification.json](qualification.json) records the tested proof/fixture
  hashes, dependency hashes, frozen import tree, and zero-error mathematical
  syntax preflight. The final replay includes the direct source-deviation
  composition fixture and the explicit positive-integer p assumption.
- [independent-audit.md](independent-audit.md) independently checks orientation,
  projective helper phases, the degree/rank argument and the robust extension.

From a repository containing the PR270 test dependencies, run this command,
substituting the actual checkpoint path:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD" python path/to/test_weak_sector_degree.py -v
```

The fixture imports `tests.field_sector_tableaux` and `tests.exact_matrices`.
The former uses the retained recursive field arithmetic and native recoupler
fixtures; the latter uses the retained exact fault-tolerant arithmetic. The
qualification receipt binds these dependencies. No dense high-width quantum
simulation, floating singular-value estimate, or initialized-helper promise
is used.

The six checks cover a sharp address-split degree family, the actual augmented
seed and its leading coefficient, a complete non-Clifford sector example with
projective helper phases, a counterexample to inappropriate Pauli rounding,
the anticommutation/parity inequalities, and reversible helper computation
that raises degree while returning all inputs.

## Separate sharpness supplement

[sharpness/PAULI_SURROGATE_SHARPNESS.md](sharpness/PAULI_SURROGATE_SHARPNESS.md)
constructs a dependent commuting Pauli bank approaching an independent bank
on every input, including a variant preserving the seed exactly. Its three
exact tests and [receipt](sharpness/SHARPNESS_RECEIPT.json) are separate from
the six tests above.

The optional individually Clifford weak-probe query family in that supplement
has full address degree p. It is explicitly not a counterexample to the
fixed-low-degree theorem and does not supply a paid native source. Its role
is to explain why a lost Pauli relation alone cannot support a width-independent
robust argument for arbitrary non-Clifford sectors.

## Next integration

Strengthen `research/depth/FIELD_SECTOR_DEGREE_BARRIER.md` as the single proof
home. Preserve its full-cocycle group-norm witness, moving-helper scope,
local-affinity proof, existing fixtures and stable links; the weak-probe
argument supplies the stronger contract rather than erasing those distinct
identities. Review whether to append the focused fixtures to its existing
test module or retain a clearly linked separate module.

The sharpness supplement belongs in a separate proof home only if its precise
surrogate-bank result is adopted. Update routing and resource claims only after
review; R96 remains unchanged. Run the full integrated gate and presentation
checks on the final reviewed tree before any merge.
