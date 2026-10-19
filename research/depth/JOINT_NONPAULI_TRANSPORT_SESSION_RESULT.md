# Revision: joint defects, bounded relations, and the complete probe task

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Bounded relations](BOUNDED_RELATION_NATIVE_MATCHING.md) ·
[Determinant defects](DETERMINANT_DEFECT_DEGREE.md) ·
[Joint probe algebra](JOINT_NONPAULI_PROBE_ALGEBRA.md) ·
[Probe and word-query transfer](SYMPLECTIC_MODULE_QUERY_BARRIER.md#6-universal-pauli-probes-contain-a-complete-dirty-indicator) ·
[Previous checkpoint](NONPAULI_BOUNDARY_SESSION_RESULT.md)

This continuation began from verified main `a85fbef` (PR #188), on
19 October 2026 in Asia/Shanghai. It replaces a pointwise description
of a non-Pauli code by statements about the entire logical family.
It also identifies exactly what a faster Pauli-probe source would
suffice to compile. The complete-frame resource frontier is unchanged.

## 1. Destination, strongest theorems, components, and scope

The destination remains optimal general T-depth for the prescribed
complete Hopf differential frame. The same circuit must satisfy the
count, width, error and depth bounds. Every logical column, literal
phase, dirty input and external reference, clean leakage, returned
work register and actual inverse remains in the contract.

| Role | Current statement |
| --- | --- |
| Strongest unconditional general upper theorem | R96: $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, and $`G=O(NL)`$ for every integer $`n\ge1,L\ge6`$ |
| Literal allocation and error | One clean qubit, real $`b\ge17(L+n+7)`$, complex magnitude $`b\ge17(L+n+8)`$, full initialized-isometry error at most $`2^{-L}`$ |
| Selected manuscript claims | R79 remains Result D; the count package is settled under its hypotheses; final manuscript drafting is separate |
| Unrestricted component bounds | At most one clean and two native layers force XOR degree at most two; zero clean and three layers force sign degree at most five and XOR degree at most four |
| Retained restricted closure | Either actual interior clean image being Pauli gives sign degree below $`2^{26}`$; both Pauli images and the two rank-one non-Pauli images retain the stronger 2048/2047 sign/XOR bounds |
| New restricted closure | With no middle pairs, bounded endpoint clean-active ranks give a width-independent degree bound, even when both interior clean images are non-Pauli |
| Exceptional paired closure | If the rank-two clean-code probe comparison takes a nonidentity type other than S, its occupied action has a zero-clean three-layer representation: sign degree at most five and XOR degree at most four |
| New joint-family bridge | In the generic zero-character branch, bounded codimension of the constant determinant-defect subgroup gives a computational degree bound after a paid Clifford correction |
| Joint algebra component | When both interior images are non-Pauli and a generic probe character is nonzero, the whole defect family and its actual factors share a Pauli algebra of symplectic rank at most four; this alone gives no target-degree bound |
| Constructive component | A complete independent-Pauli probe converts to a complete dirty indicator. A uniform faster probe with the stated simultaneous budgets supplies the entire R98 word-query interface, including every output length |
| Representation limits | Neither bounded intrinsic core size nor low degree of each separate intrinsic phase bounds the joint computational phase automatically. Universal fixed-family bilinear Pauli readout still contains the indicator task |
| Still open | Uniform unrestricted one-clean three-layer degree, unrestricted one-clean four-layer rigidity, optimal general depth and matching dimension-growing lower bounds |

At fixed accuracy and square-root dirty width, the general upper depth
remains logarithmic in n and the unrestricted lower bound remains
constant. The new restricted lower results and conditional transfer
do not close that interval.

For XOR conclusions involving a logical probe or determinant-defect
hypothesis, evaluate that hypothesis in the sign representation obtained
by conjugating the output with Hadamards. That basis change need not
preserve the full logical-Z probe family.

## 2. Revision of the first unsupported inference

The previous checkpoint gave an exact non-Pauli clean decoder with two
low-degree controls, and a common central Clifford determinant cocycle
in the no-pair, zero-character branch. Neither justified identifying
intrinsic central coordinates with computational input bits. Decoding
also has its own native cost and can destroy the independent-axis
structure needed for matching.

The [determinant-defect theorem](DETERMINANT_DEFECT_DEGREE.md) now
avoids that identification. For the entire logical-Z family write

```math
W_z=\Pi_+ +\alpha_z(u)\Pi_-,\qquad
H=\{z:\alpha_z(u)\text{ is constant in }u\},\qquad
r=\mathrm{codim}\,H.
```

The cocycle makes H a subgroup. Its constant defects are exactly
$`S^{\chi(z)}`$ for a linear character chi. A controlled X-only
Pauli, implemented by at most $`Q-1`$ CNOTs on Q physical wires,
corrects that character. It fixes the occupied clean-zero sector
literally and makes the original computational-Z subgroup
$`\langle Z_c,Z(H)\rangle`$ globally preserved. No native T layer,
T gate or wire is added; the actual inverse has the same separately
charged Clifford correction.

The retained partial-diagonal theorem therefore bounds sign degree by
$`D_r`$. For r zero the stronger full-diagonal theorem gives sign
degree at most five and XOR degree at most four. Moreover, if the
entire alpha-family factors through one common t-dimensional intrinsic
linear quotient, then

```math
r\le t(t+3)/2.
```

The common quotient is sufficient even though its individual coordinates
have not been identified with input bits. The new unsupported inference
would be to infer bounded r merely because each alpha is a low-degree
intrinsic Clifford phase. The number of such phases across the whole
family remains unbounded by the current argument. This is a limitation
of the supplied representation, not a hardness measure for its target.

## 3. Bounded native relations survive the clean quotient

The second route keeps the actual commuting phase factors and permits
bounded relations among their Pauli axes. For three full-space phase
lists with relation dimensions $`r_A,r_B,r_C`$, the
[bounded-relation theorem](BOUNDED_RELATION_NATIVE_MATCHING.md)
gives an explicit width-independent sign-degree bound
$`D(r_A,r_B,r_C)`$.

The proof first reconstructs an invariant isotropic proxy for the middle
axis list using the rank profile of its actual Pauli images. Only
boundedly many axes differ from that proxy. An intrinsic Pauli kernel
localizes the exceptional phase in a bounded invariant algebra.
Dependent pair differences leave only boundedly many unpeeled pairs;
the endpoint relation spaces also bound singleton coupling and hidden
endpoint directions. The remaining factors share a bounded quantum
core. Their central phases, the signed Clifford correction and every
peeled phase are then retained as finitely many bounded-degree features
in actual computational coordinates.

This permits a new one-clean application. If both interior clean images
are non-Pauli, have active ranks $`k_A,k_C`$, and have no paired middle
orbit, their clean-active transition is a commuting C3 phase on at most
$`d=k_A+k_C`$ coordinates. It need not be Clifford. Extracting its
actual Clifford part leaves a commuting C3 residual; its literal parity
expansion adds at most

```math
s(d)=d+3\binom d2+7\binom d3
```

axes to the middle projected list. The silent endpoint lists and the
positive middle list each acquire at most one relation on the Pauli
clean quotient. Hence the exact occupied sign degree is at most
$`D(1,1+s(d),1)`$, with one less for XOR degree.

The projected dependent list is not claimed to be one independent-axis
physical native band. Its role is justified by the new algebraic
theorem. A general bounded operation inserted between silent bands
cannot be removed for free: those bands can impose arbitrarily many
independent constraints when it is conjugated through them. Paired
non-Pauli transitions and ranks growing with width remain separate
questions. The retained R74 exactness consequence gives the scoped
degree obstructions the same approximation floor $`2^{-530690}`$.

## 4. A common algebra, and one exceptional paired closure

The [joint-probe theorem](JOINT_NONPAULI_PROBE_ALGEBRA.md) also
advances the nonzero-character branch. When both interior clean images
are non-Pauli, a generic nonzero probe character identifies the output
clean-active span with an actual subset of the independent middle axes.
The two endpoint clean-active spans become disjoint middle coordinate
subspaces. One common choice of character representatives then has
cross-pairing rank at most one with the middle span. Adding the clean
Pauli gives a common algebra of symplectic rank at most four for the
whole logical family, including its retained actual band factors.
The actual Clifford families J_z and C_z normalize this algebra.

This is a two-quantum-coordinate description with possibly many central
controls. It is not yet a computational degree theorem. The endpoint
Pauli case is outside this particular common-algebra proof and is
already degree-bounded by the previous checkpoint.

There is a sharper conclusion for the rank-two exceptional stabilizer.
Only one nonidentity comparison type can occur across the whole logical
family: two different types would force a nonisotropic subspace inside
the isotropic middle span. If that type is S and both endpoint clean
images are non-Pauli, the common-algebra proof applies. A Pauli endpoint
already has the retained degree closure. If the type is one of the
other two exceptional types, native
independence forces exactly one middle pair on the clean core, without
spectator dressing. Both endpoint clean supports lie in the same
three-dimensional Pauli space
$`W=\langle Z_v,X_t,Z_t\rangle`$. Every silent band commutes all W.

An actual Clifford in this algebra aligns the input and output clean
Paulis. The residual transition then commutes the aligned clean Pauli;
on its positive sector it is a scalar function of Z_v only. The retained
native unit-circle lemma makes both values eighth roots, so one parity-T
power and its literal scalar realize it. The remaining middle singleton
axes are independent of Z_v, because Z_v was the difference of the
removed independent pair.

Replacing that residual by its positive-sector phase constructs an
explicit occupied-equivalent circuit. Its unused clean sector need not
equal the original one. This replacement is legitimate because every
silent factor commutes W and preserves the aligned clean sector. All
remaining axes act identically on the clean wire in fixed Clifford
coordinates and stay independent on projection. The resulting full
logical operator has an actual zero-clean three-native-band expression,
so its sign degree is at most five and its XOR degree at most four.
Coordinate Cliffords, the scalar and actual inverses are charged; no
new T band or work wire is used.

The special commutation with W is what permits this closure. It does
not remove a general bounded noncommuting insertion or settle every
rank-two or paired non-Pauli circuit.

## 5. Constructive revision: the probe is enough, and still unpaid

The earlier two-family construction escapes the single-module barrier,
but its complete selectors retain the old depth. The new
[universal-probe reduction](SYMPLECTIC_MODULE_QUERY_BARRIER.md#6-universal-pauli-probes-contain-a-complete-dirty-indicator)
shows more precisely why relabeling these probes does not eliminate
the selector task.

If two fixed Pauli-label families can realize every binary table by a
bilinear symplectic readout, each family must be linearly independent.
Given a complete native source

```math
M|u,A,d\rangle=e^{i\phi(u)}|u\rangle P_{\gamma_u}|A\rangle|d\rangle,
```

symplectic dual probes and an actual Clifford commutator produce the
complete dirty indicator $`Y\mapsto Y\oplus e_u`$. The labels and
their duals need not commute. Source phases and ordered Pauli scalars
cancel against actual adjoints. On an m-entry family with w carrier
wires the wrapper uses two source calls, at most $`10mw+2m`$
additional Clifford gates, no clean work and no additional T layers
beyond the two calls. All address, carrier and helper inputs return.

There is also a constructive sufficiency result. If a uniform p-bit
probe source has $`m=2^p`$ independent labels and

```math
T_M=O(m),\quad G_M=O(m^2),\quad w+b_M=O(m),\quad
D_T(M)=O(d(p)),
```

with common actual-inverse bounds and nondecreasing $`d\ge1`$, the
R91 allocation gives every q-bit word query on $`Q=2^a`$ rows with

```math
T=O(\sqrt{Qq}+Qq/B),\qquad G=O(Qq),\qquad
D_T=O(d(a)+Qq/B^2).
```

The quadratic indicator Clifford cost fits all output lengths because
the chosen two indicator lengths satisfy $`H^2+J^2\le5Qq`$;
the direct full-indicator branch has $`Q^2\le Qq`$. The four
complete indicator calls occur outside output batches and are sequential,
so temporary shared-address changes do not create an unpaid parallelism
assumption. Simultaneous dirty allocation includes every live bank,
carrier and helper.

Thus the multi-output assembly is paid under the explicit source
hypothesis. The first missing native operation is the faster complete
probe itself. A source with sublogarithmic d would improve the full
fixed-accuracy frame depth through R98. No such source is constructed
here. Table-dependent probes, non-Pauli readout and a fused native
construction lie outside the reduction; they remain possible changes
of mechanism, not established improvements.

## 6. Next research question and evidence

The next lower question is quantitative and joint: do the actual native
endpoint lists and the literal Boolean occupied action bound the number
of distinct determinant defects, or supply another paid representation
with bounded defect codimension? The bounded-relation theorem closes
the no-pair case when both active ranks are bounded. The exceptional
probe type above is now also closed. The remaining generic paired and
growing-rank cases do not justify deleting their non-Pauli transition.

The constructive question is now the uniform complete probe interface
above. Its budget permits quadratic Clifford count at the probe level,
because the full word-query allocation absorbs it. This gives room to
change the source mechanism while preserving the general destination.
Improving a label description without improving that native source
retains the old cost.

The [dependent-axis fixtures](../../tests/test_bounded_relation_native_matching.py)
check exact signed rank profiles, the invariant majority proxy and a
profile symmetry which does not give a Clifford passage. The
[determinant fixture](../../tests/test_determinant_defect_degree.py)
retains a nonconstant phase family and its nontrivial constant-defect
subgroup, the paid correcting CNOT and actual native inverses, despite
an identity occupied target. The
[exceptional native fixture](../../tests/test_joint_nonpauli_probe_algebra.py)
checks the paired clean transition and its literal occupied replacement.
The
[native probe fixtures](../../tests/test_symplectic_module_query_barrier.py)
check every physical column and actual inverse, address-dependent
scalars, noncommuting dual labels, asymmetric forward/inverse Clifford
counts and every R91 allocation branch. Independent proof review checks
the arbitrary-width arguments and their model scope. Finite fixtures
support these proofs; they do not establish the degree theorems by
enumeration.

This checkpoint adds nine focused tests. The 47 relevant native,
presentation and reviewer tests pass locally, as do source compilation
and whitespace checks. Repository-wide validation and rendered
presentation remain required before merge.
