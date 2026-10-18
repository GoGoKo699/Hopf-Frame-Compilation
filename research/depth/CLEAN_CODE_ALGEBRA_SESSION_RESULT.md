# Clean-code revision and a paid algebra recovery

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Clean-code proof](CLEAN_CODE_CLIFFORD_RIGIDITY.md) ·
[Quadratic queries](QUADRATIC_NATIVE_QUERY.md) ·
[Balanced algebra construction](BALANCED_ALGEBRA_PERMUTATION.md) ·
[Previous checkpoint](THREE_LAYER_NATIVE_SESSION_RESULT.md)

This round began from verified main
`c1daf7da61ed22851762dd17e28f3e5d7134a904` (PR #185), on
18 October 2026. It examined the first unpaid one-clean extension and
the in-place recovery cost of the historical HRT permutation family.
The general complete-frame resource frontier remains unchanged.

## 1. Revision: destination, theorems, components and limitations

The destination is optimal general T-depth for the prescribed complete
Hopf differential frame. Every logical column, literal sign and phase,
singular and zero coefficient, arbitrary dirty input and external
reference, actual inverse, leakage and complete workspace return remains
part of the operator contract. Sources, controls, queries and recovery
are charged on the same circuit as the frame.

| Category | Verified statement |
| --- | --- |
| Strongest unconditional general upper law | R96 gives $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$ and $`G=O(NL)`$ simultaneously for all integers $`n\ge1,L\ge6`$ |
| Literal allocation and error | One clean qubit; $`b\ge17(L+n+7)`$ for real frames or $`b\ge17(L+n+8)`$ for complex magnitude frames; full initialized-isometry error at most $`2^{-L}`$ |
| Selected claims | The count package is settled under its stated hypotheses; R79 remains selected Result D; final manuscript drafting is separate |
| Unrestricted query components already proved | With at most one clean qubit, two native layers force Boolean XOR degree at most two. With zero clean qubits, three layers force degree at most four. Both have the retained uniform nonattainment floor |
| Other useful components | R97 is an outer-overlap obstruction under its stated misalignment hypotheses. Conditional R98 transfers a faster completed arbitrary-table word query to fixed-accuracy frame depth |
| Open general questions | Optimal general depth, matching dimension-growing lower bounds, unrestricted one-clean four-layer rigidity, and a uniform one-clean three-layer query degree bound remain open |

At fixed accuracy and square-root dirty width, the general upper depth
is still logarithmic in n and the available unrestricted lower depth is
constant. A component construction changes this frontier only if it
fits the full R98 query budgets or compiles the complete frame by another
fully paid route. A representation parameter gives an unrestricted
lower bound only if alternative native representations cannot remove it.

The round chose two complementary tests. For a lower invariant, the
first unsupported inference was that an actual two-band Pauli image
agreeing with a Clifford on the clean code must itself be Clifford on
the whole space. For a construction, the first unpaid operation was
the in-place update $`a\mapsto(1+b)a`$ in a square-zero coefficient
algebra: an out-of-place multiplication query alone did not supply it.
Historical next-step passages were treated as evidence, not work orders.

## 2. What the one-clean analysis establishes

The [clean-code proof](CLEAN_CODE_CLIFFORD_RIGIDITY.md) classifies
actual Cliffords fixing the positive sector of a non-Pauli one-band
clean image S. If its intrinsic active rank is not two, only I and S
fix that sector pointwise. Rank two has one exact Klein four exception.
Two literal native Clifford restrictions that differ are separated in
operator norm by at least $`1-1/\sqrt2`$. This comparison gap concerns
actual Cliffords, not arbitrary operators with Clifford occupied action.

Applying the classification to native Pauli probes gives a linear
character and an actual supported two-band factorization. Its nonzero
branch has at most two quantum coordinates, with commuting controls;
the exceptional rank-two code has a bound of five. In the remaining
branch, the complete logical comparison family preserves one common
Pauli space. Its symplectic rank is at most $`2m+2`$, where m counts
paired native axes under the clean Clifford. The central controls and
the m-dependent rank are not yet bounded logical degree.

Two exact diagnostics prevent overpromotion. First, the failed shared
Pauli clean flip persists at arbitrarily large intrinsic rank with only
one initialized qubit. In that same family, m grows while every logical
comparison defect is identity; the full target simplifies from three
displayed bands to one. Second, a correlated one-band/two-band encoder
pair has $`A\Pi=C\Pi`$ but non-Clifford A. Its whole transition also
has a paid one-band simplification; its occupied identity target alone
needs no gates. Thus the original full-space
extension inference is false, not simply unproved. Neither diagnostic
is a hard-target family.

The retained R74 identity-isolation theorem also gives the
[Boolean rounding corollary](CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action):
at three layers and one clean qubit, sufficiently accurate Boolean
approximation forces exact occupied action, with an explicit radius
independent of arity and dirty width. This is a consequence of an
existing theorem, not a new degree bound.

The surviving question is about the complete occupied logical action
after accounting for unused-sector freedom. It cannot be answered by
counting clean pairs, by assuming a shared flip, or by applying native
axis matching directly to a non-Clifford image.

## 3. A complete optimal-depth quadratic query

The [quadratic construction](QUADRATIC_NATIVE_QUERY.md) proves that
every nonaffine quadratic Boolean XOR on m address bits has exact
optimal T-depth two with zero clean qubits and one arbitrary returned
dirty helper. If the quadratic polar rank is r, the circuit uses

```math
T=2r+4,\qquad D_T=2,\qquad G=O(m^2+r).
```

Two explicit CNOT parity bases implement an even-arity symmetric
quadratic sign. A paid invertible address change reduces every quadratic
to that case, followed by its affine correction. Both output values,
every address and dirty column, literal phases and native adjoints are
covered. The existing full-input one-layer floor supplies optimality;
T-count and Clifford count are not claimed optimal. Affine queries have
depth zero. The proof records exact forward and reversed Clifford
ledgers, including the different S-gate costs.

This construction also changes the fair baseline for the attempted
quartic-phase echo. For the seven-address symmetric cubic XOR, the
echo has 60 T, 392 elementary Cliffords and six bands with one clean
and one dirty helper. A quadratic compute, exact Toffoli, and actual
uncompute instead uses 40 T, 148 Cliffords and six bands at the same
allocation; the literal reversed word has 152 Cliffords. The echo is
therefore stopped as a depth improvement. The general odd-arity cubic
baseline is a fixed-degree component, not an arbitrary-table query.

The two-band quadratic word applies T gates to transformed addresses.
Different rows cannot share those physical T targets for free. This
does not replace R91's simultaneous shared-control word query or reduce
the balanced construction's six-call multiplier to twelve layers
without another native scheduling proof.

## 4. The in-place algebra recovery is now paid

The [balanced construction](BALANCED_ALGEBRA_PERMUTATION.md) supplies
the missing complete update $`a\mapsto(1+b)a`$, where $`b^2=0`$.
Treating whole coefficient banks as rank-one modules gives a
determinant obstruction to a naive shear-only recovery. Splitting one
algebra generator changes the scalar ring and doubles the module rank;
the desired matrix then has determinant one. Two multiplication
shears, two copies and a dirty correction echo implement it literally.
The correction absorbs the extra term left by the two shears. Six
completed bilinear query calls charge every helper and inverse;
the first shear and correction preparation share one native stage
with separate private T targets. This gives five eight-layer stages.
The determinant obstruction was a limitation of the former coordinates.

This permits a balanced merge tree for the HRT permutation on
$`w=2^r-1`$ arbitrary coefficient inputs. For each integer $`r\ge1`$,
the source and its actual inverse each obey

```math
T\lt48\,3^r,\qquad G\lt60\,3^r,\qquad
D_T\le40\lceil\log_2r\rceil,\qquad
B_{\mathrm{extra}}\lt2\,3^r,
```

with zero clean qubits, simultaneous extra dirty allocation as shown,
exact work return, zero error and zero leakage. At r=1 the source is
identity. Structural zero coefficients are omitted algebraic terms,
not physical qubits assumed initialized. Nodes at one level use
disjoint physical coefficient banks and explicit disjoint helper pools.

The comparison is against a fairly simplified baseline. Batching the
quadratic inverse by decreasing target cardinality gives depth at most
$`8(r-1)`$, with the same $`O(3^r)`$ count order and fully returned
dirty helpers. The balanced source improves linear to logarithmic depth
in r, rather than comparing only with the old serial Toffoli ledger.
It does not claim superior small-instance constants or optimal counts.

Conjugating an output CNOT by the source and its actual inverse gives
a literal returned query of Boolean degree r and depth at most
$`80\lceil\log_2r\rceil`$. Its address is **all** w coefficient
inputs. Restricting nonsingleton coefficients to zero would add clean
promises; it does not yield an r-input AND with those wires arbitrary
and returned. This structured family is not an arbitrary-table R98
primitive and has growing depth. It neither refutes nor proves a
possible general bound $`\deg f\le2^{O(D_T)}`$.

## 5. Evidence and the concrete resumption question

Each result has one authoritative proof home and bounded exact evidence:

- [Clean-code fixtures](../../tests/test_clean_code_clifford_rigidity.py)
  retain literal Clifford restrictions, the exceptional group, native
  encoder comparisons, failed extension and alternative one-layer words.
- [Quadratic fixtures](../../tests/test_quadratic_native_query.py)
  check parity identities, general polar congruence, complete native
  matrices, actual inverse ledgers and the fair cubic baseline.
- [Algebra fixtures](../../tests/test_balanced_algebra_permutation.py)
  check the native phase derivative, arbitrary dirty multiplier columns,
  shared-bank correction, balanced coefficient permutation and inverse,
  and simultaneous symbolic ledgers.

The arbitrary-width claims follow from the algebra and native proofs,
not from extrapolating those finite checks. The continuation recovered
the previously local commit, independently audited its proofs, and
refined the multiplier from seven serial calls to six calls in five
stages without increasing its stated dirty-work bound. No complete-frame claim,
precision regime or selected manuscript result has been promoted.

The priority resumption question is: **does every literal Boolean XOR
query with one clean qubit and three native T layers have degree bounded
independently of arity and dirty width?** Use the common logical defect
space and the whole comparison family, while quotienting or explicitly
handling the unused-sector action. A proof must control the logical
contribution of arbitrarily many clean-permuted pairs; an exact
unbounded-degree family would resolve the test in the opposite direction.
The new counterexamples forbid a blanket full-Clifford extension or a
pair-count lower bound as substitutes. This is a genuine unrestricted
constant-layer question and a necessary test for a wider degree-versus-
depth invariant, not already a dimension-growing frame lower bound.

The constructive alternative remains concrete: can the generator split
and returned dirty correction remove more than one balanced merge
level, or pay an actual recovery step of a faster arbitrary-table query
within all R98 budgets? A cheap out-of-place source or fixed-degree query
alone does not do this. The current results resolve specific operations
and refute specific inferences; they preserve the direction toward
optimal general depth.
