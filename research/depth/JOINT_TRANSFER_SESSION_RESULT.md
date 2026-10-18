# Unrestricted two-layer queries: joint constraints and native return

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Authoritative proof](UNRESTRICTED_TWO_LAYER_QUERY.md) ·
[Previous checkpoint](RESTORED_ADDRESS_SESSION_RESULT.md) ·
[Conditional upper transfer](QUERY_DEPTH_FRAME_TRANSFER.md)

This round began on 18 October 2026 at verified main
`1ac16c4a85eac3e202576d788d7414565c1d614d` (PR #182). It first
tested whether restored-address access permits a cheaper complete query,
then changed mechanism when the native costs reproduced the baseline.
The supported advance is an unrestricted constant-layer query obstruction
and an exact one-clean return invariant. No complete-frame upper or
dimension-growing lower bound is claimed.

## 1. Revision: destination, frontier and components

The destination remains optimal worst-case T-depth for the prescribed
complete Hopf differential frame under the simultaneous resource model.
The count package is settled under its stated hypotheses. R96 remains
the strongest unconditional general upper law, with absolute constants,
for every integer $`n\ge1,L\ge6`$ and $`N=2^n`$:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

It uses one clean qubit, $`b\ge17(L+n+7)`$ dirty qubits for real
frames or $`b\ge17(L+n+8)`$ for complex magnitude frames, and full
initialized-isometry error at most $`2^{-L}`$. Every logical column,
literal sign and phase, singular/zero coefficient, dirty input and
external reference is covered. Actual inverses, leakage, complete work
return and every source/query/control/recovery cost remain charged.
R79 remains selected Result D; R97 remains scoped to outer overlap.

The count lower bound divided by physical layer capacity already matches
the $`NL/b^2`$ contribution where it dominates the additive terms in
the eligible interval. The large-width additive frontier remains open:
in particular, at fixed accuracy and square-root dirty width, must depth
grow with n, or can the current logarithmic additive term be improved?
Unrestricted one-clean four-layer rigidity also remains open.

Useful components are distinct from that frontier. Conditional R98 now
accepts exact **completed address return**, rather than gatewise read-only
access, but still needs a faster arbitrary-table query within its original
count/Clifford/work bounds. The structured inverse query is complete and
logarithmic-depth; arbitrary-table inverse-entry representations already
exist. Read-only query degree bounds and the six-layer source certificates
remain restricted results. Complete return does not supply their source
or intermediate-address hypotheses. These resolved questions are not
reopened by historical next-step passages.

## 2. Constructive candidates and the first unpaid operation

The first question was whether a fused sparse inverse extractor could
remove both dirty loading offsets and the old indicator depth. Fixing
the table cross edges gives a block inverse entry of the form
$`R(r) F L(l)`$. Its dirty four-corner extraction is precisely the
bilinear indicator mechanism already paid in
[R91](PARALLEL_BILINEAR_WORD_QUERY.md). Sparse specialization removes
irrelevant generic matrix work, but retains the balanced conjunction
cost. Generic inversion and repeated interpolation are more expensive,
as already priced in the [representation audit](GAUSS_SUM_QUERY_ENCODING.md).
No improved full query follows from renaming this expression an inverse.

A second candidate conjugates an affine-address quadratic phase by
address Hadamards. For
$`D(x,w)=(-1)^{q_0(w)+\sum_i x_i g_i(w)}`$ with quadratic
$`q_0,g_i`$, the exact identity is

```math
H_X D H_X=C_{q_0}P_g,\qquad
P_g|x,w\rangle=|x\mathbin\oplus g(w),w\rangle.
```

The fixed table phase $`C_{q_0}`$ commutes with this address shear.
This identity supplies no scalar table extraction. A useful successor
would have to include work interference and prove its complete return,
not just change the basis of the same source.

Dirty coherent teleportation has a similarly identifiable unpaid step.
For the Boolean phase $`U_f`$, its correction is
$`(-1)^{f(x+e)+f(x)}`$. Under the Clifford change of variables
$`y=x+e`$, this is exactly $`U_f(x)\otimes U_f(y)`$.
The direct correction therefore contains two paid target phases. For
$`f(a,b,c)=abc`$ and mask $`(r,s,0)`$, the derivative is
$`rbc+sac+rsc`$: separate single-mask corrections miss the last term.
Coherently controlling a Clifford correction may itself require T gates.
This diagnoses that implementation; it is not a lower bound against
joint resynthesis of the full word.

These candidates retain the old cost. The round therefore switched to
a question whose positive answer would resolve a genuinely unrestricted
lower-bound case: what joint constraint is lost when one-clean logical
transfer coefficients are treated as independent sums?

## 3. Supported lower results

The [proof home](UNRESTRICTED_TWO_LAYER_QUERY.md) derives the initialized
coefficient $`u+v`$ from two columns of the same physical transfer row.
At two native T layers, each summand is zero or a signed
$`2^{-j/2}`$, but row orthogonality additionally gives
$`u^2+v^2\le1`$. This excludes the previously admissible formal
decomposition $`1-2^{-k}`$. It does not assume intermediate address
preservation, exact clean return, or bounded dirty width.

This first argument excludes exact AND on four or more controls and
gives a positive error floor for four controls, but still permits cubic
bias 3/4 through the pair 1/2,1/4. The missing relation is supplied by
the **shared native encoding of the clean qubit**. Factor the actual
circuit into two one-band encoders with common clean image S. Its
intrinsic active-axis span is the same for both encoders. If nonzero,
one common axis gives the two actual trace contributions opposite
conjugation signs, so one must vanish. If zero, S is a Pauli and the
occupied overlap is a logical Clifford trace. Thus every occupied
transfer coefficient of an exactly clean-returning two-layer circuit
lies in the single alphabet $`0,\pm2^{-j/2}`$. This theorem needs
no rational or Boolean target premise.

For exact Boolean XOR, rationality then bounds every nonunit Walsh
magnitude by 1/2. Every nonaffine qualifying function is at relative
Hamming distance at least **1/4 from every affine function**. In
particular, AND on three or more controls is excluded, as is every
three-variable cubic Boolean query. These are necessary conditions,
not an arbitrary-width degree theorem or a synthesis criterion.

R75's retained rounding theorem forces the same circuit to return its
clean qubit exactly below its constant error threshold. Combining it
with the new alphabet gives the full initialized-isometry floors

```math
\eta\ge\frac{3-2\sqrt2}{8}\quad\hbox{for three-control AND},
\qquad \eta\ge2^{-a}\quad\hbox{for AND on }a\ge4\hbox{ controls}.
```

Three-control AND therefore requires at least three T layers already
at $`L\ge6`$, and the four-control floor improves to 1/16. Leakage
and approximate work return are included. The bias witness still has a
shrinking floor for growing arity; fixing logical controls would add
clean input promises and cannot make this floor uniform.

A second exact result uses the native field automorphism
$`\sigma(\sqrt2)=-\sqrt2`$, $`\sigma(i)=i`$. Literal one-clean
return to a target with entries in $`\mathbb Q(i)`$ forces

```math
V^\dagger\sigma(V)\in\{I,Z_c\}.
```

The proof combines the entire occupied code with independence of the
physical axes in each T layer. It also makes every transfer column whose
input Pauli commutes with $`Z_c`$ rational. Exact equality is essential;
there is no asserted approximate-target extension. A fully charged
two-clean counterexample and a rational non-C3 two-layer circuit delimit
the result. The correlated-axis theorem above resolves the three-control
case that Galois rationality alone left open.

These are unrestricted native two-layer statements, rather than limits
of the present inverse or read-only representation. They do not strengthen
the existing complete-frame lower bound via the two-call corner reduction,
and do not change R96's upper ledger. Their role is a new general-model
base case and a constraint to test at the next depth.

## 4. Pause and concrete resumption

The successful proof used information lost in two failed relaxations.
Replacing the native midpoint by an arbitrary
Hermitian Clifford is too weak: $`S=Z_c\mathrm{CZ}_{bd}`$ and
$`G=\mathrm{CZ}_{ca}`$ commute, and on the $`S=+1`$ code one has
$`c=bd`$, so G restricts to $`(-1)^{abd}`$. Their normalized traces
are 1/2 and 1/4. This S is not a one-native-band Pauli image: its
symplectic displacement is alternating, whereas a nonempty independent
native active-axis set has nonzero displacement form.

The new theorem resolves that trace pair, but scalar transfer values
alone do not force Boolean degree two. The cubic bent function on six
address bits,

```math
f(x,y)=x\cdot y+x_1x_2x_3,\qquad x,y\in\mathbb F_2^3,
```

has an entire XOR Pauli-transfer matrix with magnitudes only
$`0,1,1/2,1/4,1/8`$. It passes the stronger scalar alphabet, not
just the Walsh bias test. The proof home gives the exact formula and
bounded certificate. This supplies no native two-layer construction.

The concrete next lower question is whether a **one-clean two-layer
query for this cubic bent function** exists, or whether the shared native
encoders impose a joint logical-Pauli relation excluding it. A construction
must charge every gate and prove all occupied/dirty columns; an obstruction
must survive arbitrary native Clifford mixing. Neither a favorable
scalar spectrum nor target rationality answers the question.

This is a controlled next case toward an invariant surviving arbitrary
representations and eventually more layers. The present three-control
resolution is a constant-layer improvement, not optimal general depth.
On the upper side, the continuing target is a completed restored-address
query with $`d(a)=o(\log(a+2))`$ within R98's simultaneous budgets;
the first remaining operation is paid table extraction with dirty-mask
cancellation and full return. The session stops these specific failed
mechanisms, not the general project. Final manuscript drafting is separate.

## 5. Evidence

The proof home links exact bounded tests for the joint row normalization,
infinite-alphabet tail bounds, Walsh distances, complete AND permutations,
native Galois action, actual inverses and scope counterexamples. The
universal statements are analytic; finite cases alone do not prove them.
