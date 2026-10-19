# Revision: clean-transition phase and representation cost

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Finite insertion](FINITE_INSERTION_THREE_LAYER_FEATURES.md) ·
[Easy-transition counterfamilies](CLEAN_TRANSITION_PARAMETERS.md) ·
[Rank-one interface](RANK_ONE_CLEAN_PHASE_INTERFACE.md) ·
[Previous checkpoint](JOINT_NONPAULI_TRANSPORT_SESSION_RESULT.md)

This continuation began from verified main `e773bd0` (PR #189), on
19 October 2026 in Asia/Shanghai. It first revised the distinction
between a target obstruction and the cost of a chosen representation.
The resulting positive theorem closes fixed endpoint clean-active
ranks, including arbitrary middle pairs. Two explicit easy families
disprove proposed uniform bounds on representation parameters.
The complete-frame resource frontier remains unchanged.

## 1. Destination and current evidence

The destination remains optimal general T-depth for the prescribed
complete Hopf differential frame, with all resource bounds satisfied
by the same native circuit. Count, width, accuracy and Clifford cost
cannot be traded away without a new stated theorem.

| Role | Current assessment |
| --- | --- |
| Strongest unconditional general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$ for integer $`n\ge1,L\ge6`$ |
| Literal allocation and error | One clean qubit; real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$; full initialized-isometry error at most $`2^{-L}`$ |
| Selected manuscript theorems | R79 remains Result D; the count package is settled under its hypotheses; final manuscript drafting is separate |
| Unrestricted low-layer components | At most one clean and two native layers force XOR degree at most two; zero clean and three layers force sign degree at most five and XOR degree at most four |
| Strongest new fixed-rank closure | One clean and three native bands, with fixed endpoint clean-active ranks, have a width-independent Boolean degree bound, including every middle-pair type |
| Retained sharper special cases | Pauli interior images, the rank-one/rank-one case and the genuine exceptional non-S rank-two type keep their stronger bounds |
| Useful phase components | Partial-diagonal features, bounded relations, joint probe algebra, determinant defects and the paid occupied-transport form remain valid under their stated hypotheses |
| New representation limits | A single residual parity list can require unbounded relation dimension on degree-three targets; determinant codimension can grow on occupied identity; a one-core adjacent interaction does not bound raw rank distortion |
| Constructive transfer | A faster complete independent-Pauli probe with linear T/width and quadratic Clifford budgets would supply the full R98 word-query interface; that source is still missing |
| Still open | Uniform one-clean three-band degree for growing active ranks, unrestricted four-band rigidity, optimal general depth and matching growing lower bounds |

At fixed accuracy and square-root dirty width, the general upper depth
is still logarithmic in n and the unrestricted lower bound is still
constant. The new results narrow a lower-bound branch; they do not
close that interval or improve the general upper law.

## 2. First unsupported inferences in the current phase route

The previous no-middle-pair reduction produced a residual commuting
phase after aligning the clean Pauli. It did not give a uniform
relation budget for one phase list representing that residual.
The [first counterfamily](CLEAN_TRANSITION_PARAMETERS.md#1-three-independent-clean-active-lists)
has three independent physical native bands, both interior images
non-Pauli, no middle pairs and complete return. Its occupied target
is a product of n disjoint CCZ gates, with degree three and a literal
zero-clean three-band implementation.

Nevertheless, **every** single parity-T representation of that target,
allowing arbitrary replacement axes, a diagonal Clifford, a scalar and
extra untouched variables, has relation dimension at least 3n.
The proof combines the mod-four identity $`AA^{\mathsf T}=0`$
with the cubic target's translation derivatives. This invalidates
uniform compression into the proposed single-list representation.
It provides no unrestricted native lower bound: the easy implementation
is already exhibited.

The [second counterfamily](CLEAN_TRANSITION_PARAMETERS.md#5-unbounded-determinant-codimension-for-literal-occupied-identity)
has occupied identity, while the determinant-defect constant subgroup
has codimension n. Its whole family is even a pure coboundary of one
common low-degree intrinsic phase. Thus uniform codimension under the
existing determinant theorem's native hypotheses is false, not just
unproved. One interior clean image is Pauli; adding a both-non-Pauli
hypothesis is a genuinely narrower question. The displayed word also
simplifies to two bands, underscoring the parameter's dependence on
the supplied split. The conditional fixed-codimension theorem remains
valid and useful.

Solving either parameter problem in a special representation would
not by itself change the full resource frontier. Any obstruction must
be tested against alternative native words and against the actual
unrestricted circuit model.

## 3. Positive change of mechanism: keep the finite transition

The [finite-insertion theorem](FINITE_INSERTION_THREE_LAYER_FEATURES.md)
handles three commuting C3 phase lists with relation bounds
$`r_A,r_B,r_C`$, together with an arbitrary fixed unitary insertion F
of Pauli-support dimension at most d. The essential additional premise
is that F commutes with **each actual middle-axis Pauli**. F can be
non-C3 and supported on a nonisotropic Pauli space.

Its finite Pauli coefficient norm bounds the distortion of the entire
middle rank profile. A strict-majority invariant axis set leaves only
boundedly many exceptional axes. Signed Clifford extension on the
intrinsic Pauli kernel localizes the remaining operation, and pair
peeling works by pairing $`\mathcal V/W`$ with $`W^\perp`$;
it does not require W isotropic. The completed phase is then determined
by finitely many bounded-degree features. Explicit constants give

```math
\deg g\le D_{\rm ins}(r_A,r_B,r_C,d).
```

For the actual clean-code application with both interior images
non-Pauli, collect both endpoint active spans and the middle-pair axes
in S_0. The Pauli cases retain their sharper previous bounds.
Independent pair differences
lie in the intersection of the endpoint spans, so
$`\dim S_0\le k_A+k_C`$. Negative singleton axes already lie in
their sum. The clean-transition word F_0 belongs to Alg(S_0) and
sends P_A to P_C. Extract an actual signed Clifford C_0 in the same
algebra. Its residual G_0 commutes P_A and has support dimension at
most $`k_A+k_C-1`$; every positive middle singleton commutes it.

Retain both silent endpoint bands, conjugating the C band by C_0.
Restriction to the common positive Pauli code produces a full logical
operator with three relation-at-most-one lists and the finite insertion.
Thus, including all middle pairs,

```math
\deg g\le D_{\rm ins}(1,1,1,k_A+k_C-1),\qquad
\deg f\le D_{\rm ins}(1,1,1,k_A+k_C-1)-1
```

for Boolean sign and XOR actions respectively. The retained R74
rounding consequence transfers these scoped obstructions below error
$`2^{-530690}`$. The result is uniform in dirty width for fixed
endpoint ranks, but its constant grows with those ranks.

This is an exact algebraic reduction, not a compiler assigning a free
gate or free native band to G_0. Its source word, every scalar and its
actual inverse remain specified. Physical coordinate Cliffords require
their own $`O(Q^2)`$ ledger on the same Q wires.

## 4. Next question: retain the original phase lists at growing rank

The selected continuation asks for control of the **completed Boolean
phase**, preserving the original independent lists instead of bounding
the relation count of a flattened phase.

The first analytic step is the
[rank-one clean-phase interface](RANK_ONE_CLEAN_PHASE_INTERFACE.md).
In the zero-middle-pair branch, an exact choice of clean section gives
four phase factors, each with commuting axes, whose relation dimensions are
$`0,\le1,\le1,0`$. The third factor commutes every middle axis;
its cross-pairing with the fourth factor has rank at most one.
This preserves more native structure even as clean-active ranks grow.
It does not add a free physical band or delete a silent factor.

There is already a precise failure test for a tempting next inference.
An actual clean-return identity family gives this same rank-one interface
while the raw middle displacement rank changes from zero to an
arbitrarily large value. Therefore replacing the finite-insertion
support parameter by the interface's quantum-core size cannot justify
the norm-distortion argument. The large distortion is exactly canceled
by another factor. This limitation applies to that proof step, not to
the target or to every method using the four-factor normal form.

The remaining question is whether the original-list relations and
their exact correlated cancellation constrain the final Boolean phase
uniformly. A proof must retain the scalar and computational input
identification through that cancellation. A counterexample must be a
complete native family with growing target degree, not merely growing
rank, relation dimension or intrinsic phase coordinates.

## 5. Constructive route: the first unpaid native phase

The [Walsh candidate audit](SYMPLECTIC_MODULE_QUERY_BARRIER.md#71-a-dense-walsh-phase-program-does-not-supply-the-faster-source)
starts from an exact full-input identity for the dirty indicator.
A Gray-code parity schedule pays its quadratic CNOT network and returns
every coordinate. But its first $`R(\pi/2^p)`$ phase is not exactly
Clifford+T for $`p\ge3`$. Even supplying those formal rotations
would leave quadratic non-Clifford count and exponentially many phase
bands. The quadratic Clifford allowance does not pay either cost.

Resynthesis confined to affine-parity CNOT/X/S/T programs is excluded
by a fourth finite difference. General Clifford/Hadamard interference
and a fused native source remain outside that obstruction. The R98
transfer is already paid; the missing task is still a complete faster
probe or indicator with all its simultaneous resource budgets. Merely
rewriting its scalar table identity cannot improve the frontier.

## 6. Evidence and continuing contract

The proofs are analytic and independently reviewed. Exact fixtures
check the new finite-insertion examples, the two native counterfamilies,
the rank-one distortion witness, their literal phases and actual inverses,
and the formal Walsh schedule and its finite-difference obstruction.
Finite fixtures support identities
and resource ledgers; arbitrary-width claims follow from the arguments.
Repository verification retains both supported Python versions, the
exact receipts, both ledgers and the mathematical presentation guards.

Every logical column, arbitrary dirty input, external reference, literal
phase, clean leakage/return condition and actual inverse remains part
of the problem. The next mechanism is allowed to change; retaining an
old unpaid cost is not progress toward the destination.
