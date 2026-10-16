# Five-hour design: the unrestricted one-clean four-layer gap

[Current upper theorem R96](ALL_PRECISION_LINEAR_DEPTH.md) ·
[Remaining lower-bound interface](ONE_CLEAN_MIDPOINT_RIGIDITY.md#6-the-remaining-one-clean-question) ·
[Inner-spread theorem R77](ONE_CLEAN_INNER_SPREAD.md) ·
[Continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer)

**Status: designed, not executed.** Prepared on 16 October 2026 from
main `7f9d52ef970b703dde06f6de351069edecf37cea`, after
[PR #162](https://github.com/GoGoKo699/Hopf-Frame-Compilation/pull/162).
The five-hour research clock starts when this prompt is executed.
Record the actual starting commit, start time, 90-minute qualification
time and five-hour ceiling. The ceiling includes review, validation,
CI and integration.

Work only in **GoGoKo699/Hopf-Frame-Compilation**. Modification and merge
remain authorized. Final manuscript
drafting is separate. Preserve R96, the selected manuscript Result D
(R79), and all retained scientific receipts. This document proposes a
proof strategy; it establishes no new depth lower bound.

## 1. Decision and measurable target

R96 proves the all-precision upper law

```math
D_T=O(NL/b^2+L+\ell),\qquad
N=2^n,\quad \ell=\log_2(n+2),\quad n\ge1,\quad L\ge6,
```

with the original optimal counts, one external clean qubit, literal dirty
widths and full initialized-isometry contract. At fixed accuracy and
sufficient square-root width its depth is logarithmic in n. Optimal
general depth remains the destination.

The next primary target is a missing lower-bound statement: find an
absolute positive constant $`\varepsilon_*`$ such that every native
Clifford+T circuit V of T-depth at most four, using one initialized clean
qubit and any finite number b of arbitrary dirty qubits, obeys

```math
\inf_{\phi\in\mathbb R}
\|VJ_1-e^{i\phi}J_1(W_n\otimes I_b)\|_{\rm op}
\ge\varepsilon_*,                                    \qquad(1)
```

```math
W_n=I_{2^{n-1}}\otimes R_y(\theta),\qquad
\theta=\tfrac12\arccos(1/3),\qquad n\ge1.
```

J1 initializes only the clean qubit. The norm includes every logical
column, occupied dirty input, external reference and clean leakage.
Allow arbitrary Clifford interlayers and unrestricted finite gate counts.
Do not assume a Pauli or Clifford midpoint, bounded crossing rank,
bounded inner clean spread, or finite projective order. A literal-phase
floor would also be useful; the inherited symmetry argument aims at the
stronger common-phase-optimized statement (1).

Success gives an unrestricted **one-clean five-layer obstruction** for
all sufficiently accurate approximations. It is a constant lower bound.
This spectator-padded witness already has constant depth at fixed error
by compiling its one active qubit. It cannot prove dimension-growing
depth or the optimality of R96's logarithmic term. Retain that distinction
in every outcome.

## 2. The exact missing implication

Use the retained R75 constant

```math
\rho=\frac{\sin(\pi/8)}{\sqrt2},\qquad \tau=\rho/2.
```

Approximation below tau first forces exact preservation of the terminal
clean sector. With $`\Pi=(I+Z_C)/2`$, the sufficient new implication is

```math
D_T(V)\le4,\quad [V,\Pi]=0,\quad
\|[V,\Pi P]\|_{\rm op}\lt c
\quad\Longrightarrow\quad [V,\Pi P]=0,                \qquad(2)
```

for an absolute $`c>0`$ and every data/dirty Hermitian Pauli P, acting
trivially on the physical clean qubit. The constant must be independent
of n, b, both encoder spreads and every native representation.

Write $`V=U_B^\dagger U_A`$, with each encoder of physical T-depth
at most two. Exact return gives

```math
S=U_AZ_CU_A^\dagger=U_BZ_CU_B^\dagger,\qquad Q=(I+S)/2,
```

```math
A=U_APU_A^\dagger,\qquad B=U_BPU_B^\dagger.
```

Then $`\|[V,\Pi P]\|=\|(A-B)Q\|`$. Thus (2) is a gap
for occupied images of **actual native encoders with the same clean
image**. Arbitrary abstract involutions with matching commutation
relations are insufficient.

R76's commutant lifting and occupied phase discreteness already extend
to native non-Clifford midpoints. Reuse those proved steps after freezing
the compressed conjugations. The target-preserved spectator/dirty X,Z
and active Y observables then recover the witness contradiction.

## 3. Primary mechanism: retain the commuting center

R76 already records the controlled-Pauli reformulation

```math
J_P=I-\Pi+\Pi P,
```

```math
U_AJ_PU_A^\dagger-U_BJ_PU_B^\dagger=(A-B)Q,\qquad
\mathrm{rank}_{\mathbb F_2}(F_{J_P}-I)\le2.             \qquad(3)
```

F is the binary symplectic action of a Clifford. Jp is a Hermitian
Clifford involution. Equation (3) is a starting point, not the new result:
small displacement does not bound crossing between the outer phase layers.

The proposed new step is a simultaneous symplectic reduction that keeps
large commuting algebras symbolic instead of expanding the midpoint into
$`2^t`$ Pauli terms. It has two distinct obligations.

### One encoder: derive the complete inner normal form

Let C and P0 be the commuting clean and data Paulis after the input
Clifford, and let D be the first physical phase layer with isotropic axis
span L. Identify Paulis with their binary labels and retain the joint space

```math
\mathcal V=L+\mathrm{span}_{\mathbb F_2}\{C,P_0\}.
```

Both summands are isotropic, and the second has dimension at most two.
Thus the symplectic form on this joint space has rank at most four.
Derive one Clifford coordinate change in which its Pauli algebra has a
noncommuting factor on at most two qubits, a possibly unbounded commuting
center, and an untouched multiplicity factor. Preserve the joint clean,
data and controlled-Pauli images in these same coordinates.

For the Clifford boundary J obtained from Jp by the input Clifford, the
additional exact consistency check is

```math
K(D,J)=DJD^\dagger J^\dagger,\qquad
\langle u,F_Jv\rangle=\langle u,(F_J-I)v\rangle
\quad(u,v\in L).
```

This cross-pairing has rank at most two. Keep the final J factor,
central-character action, literal phase and physical-layer independence;
specify exactly what remains finite and what can still grow with width.

After fixing a central character z, group the commuting physical axes
by their at-most-three nontrivial Z restrictions on the two-qubit factor.
Track their integer T-angle coefficients modulo eight, including the
clean/data signs and every scalar that survives conjugation. Derive the
finite alphabet of joint blocks and an explicit Fourier description of
the projectors selecting each block. Central-character multiplicity may
grow; the number of boundary terms, denominators and active coordinates
in a proposed gap estimate may not be left unbounded.

A finite inner expansion alone is already easy: J is a sum of four
Paulis, and D sends each Pauli to a Clifford. The purpose of the joint
description is to retain the correlation with the clean image through
the outer layers, not merely to rederive that four-term expansion.

Make that retained correlation explicit. Classify each independent inner
axis a by the two bits $`(\langle a,C\rangle,\langle a,P_0\rangle)`$.
The 00 class is silent for both inputs. Derive literal quarter-phase
products $`G_{10},G_{01},G_{11}`$ on the three other classes so that

```math
J_C=G_{10}G_{11}C,\qquad
J_{P_0}=G_{01}G_{11}P_0,\qquad
J_{CP_0}=G_{10}G_{01}CP_0.
```

Keep the signs and ordering conventions in these identities. The shared
11 class cancels in the product of the first two images; it is not an
independent Clifford boundary. In binary coordinates the two displacement
maps have images in one isotropic space and their products vanish in both
orders. Test whether this common independent-axis structure, preserved
together with the center description, supplies a constraint that the
generic four-boundary representation loses. These identities alone do
not pass the paired-encoder gate.

### Two encoders: use the shared clean image

Carry both complete normal forms through their actual outer phase layers.
Track the compatibility imposed by the exact equality of the two S images.
Does this equality constrain the two centers enough to freeze the
correlated occupied discrepancy without enumerating all central characters?

In the retained coordinates write
$`S=E_AJ_{C,A}E_A^\dagger=E_BJ_{C,B}E_B^\dagger`$, where
each E is the actual outer phase layer after its Clifford conjugation.
For $`F=E_B^\dagger E_A`$, the exact constraint is

```math
FJ_{C,A}F^\dagger=J_{C,B}.
```

Use the common-axis relations to test what this forces on the two data
images. Quotient out multiplicities: different center labels may give
identical joint blocks, and their mixing need not change the clean image.

A qualifying answer must provide at least one of:

- A direct, width-independent positive operator-norm gap for the actual
  correlated discrepancy.
- A finite-boundary representation to which
  [R74's estimate](FOUR_LAYER_CROSSING_RANK.md#finite-boundary-form)
  applies, with absolute bounds on its number K of summands, coefficient
  denominator parameter b and active Hadamard coordinate count m.
- A different exact normal form with an explicit uniform gap estimate
  and all hypotheses verified for every admissible encoder pair.

The first encoder's small inner pairing rank alone does not qualify.
Nor do a finite block size with uncontrolled central coefficients, a
growing number of blocks in a norm estimate, or an outer crossing rank
renamed as an interface size. A common-center assertion must be proved;
sharing one clean observable does not establish it. Alternative arguments
need not force the two centers to coincide.

**Qualification gate at minute 90:** identify the complete normal form,
the exact shared-S constraint used, and the uniform estimate or the
specific unresolved outer coupling. If width dependence remains, stop
that reduction as a proof of (2). Permit one explicitly different repair,
with a new test and a further 60-minute cap. Do not repeat the same expansion with
renamed parameters.

## 4. Parallel adversary and counterexample standards

An independent lane should try to violate the proposed center-compatibility
or bounded-boundary lemma using actual two-layer encoders. For every
candidate give physical Clifford/T words, count the four T layers and
verify the shared S exactly. Start with symbolic families; use only small
exact fixtures to expose a specific algebraic mistake.

Distinguish three outcomes:

| Counterexample establishes | Required evidence | Consequence |
|---|---|---|
| A proposed normal form is false | Native shared-image pair violating its stated conclusion | Reject that reduction; the gap may remain true |
| The sufficient gap (2) is false | Native four-layer family $`V_j`$ with exact $`[V_j,\Pi]=0`$ and $`0\lt\|[V_j,\Pi P_j]\|\to0`$ | Find another route to the witness floor; this does not itself approximate the witness |
| The witness floor (1) is false | Native four-layer family approaching the complete initialized witness in operator norm, with arbitrary dirty/reference inputs | A constructive result about this witness, still not a general frame compiler |

Small normalized traces, selected matrix elements, a good prepared column,
or spectator padding do not establish vanishing operator norm. An arbitrary
matrix family without native realizations does not establish either latter
outcome. Certify equalities and extrapolation analytically; floating-point
minima cannot prove a positive uniform gap.

R77's counterexample to a common sign-compatible flip remains a required
check. Its nonzero defect stays constant under padding, so it refutes that
shortcut without refuting (2). Preserve these quantifier distinctions.

## 5. Keep a short constructive scout

Give one parallel lane at most 30 minutes to test a genuinely changed
arbitrary-table phase-query mechanism. Its qualifying interface is an
arbitrary Q-row table, $`Q=2^a`$, at sufficient $`B=\Theta(\sqrt Q)`$
dirty width, with $`T=O(\sqrt Q)`$, $`G=O(Q)`$ and
$`D_T=o(\log(a+2))`$. Require full operator correctness on arbitrary
outputs and dirty/reference inputs, or an explicitly paid initialized
contract compatible with the complete frame. State every error and clean
allocation before claiming a transfer.

[R91](PARALLEL_BILINEAR_WORD_QUERY.md) already makes the bilinear middle
constant-depth at that width; its large-address indicators retain the
logarithmic allowance. Selected routing, zero predicates and the current
low-precision opening also require repricing in a full improvement.
Faster individual anchors alone do not improve the general upper order.

Do not substitute separately accurate hash equalities for a complete
dirty indicator. For a linear hash of s bits into fewer than s bits,
the bucket of x has even size. On the output Fourier character with
coefficient zero at x and one elsewhere, the exact one-hot toggle has
phase +1 while the hashed bucket toggle has phase -1. A coherent hash
seed and actual uncompute do not remove this full-operator discrepancy.
A new fused query must control the complete table-dependent phase.

If no complete new identity qualifies, end this scout and keep the
lower-bound lane primary. A successful scout may justify a separate
constructive continuation only after an independent full-resource audit.
Do not interrupt a qualified primary proof for an unpriced suggestion.

## 6. Model checks and literature boundary

Check primary sources when a proposed invariant needs an imported premise.
The design review checked these two relevant scope boundaries:

- [Dutta et al., arXiv:2502.01433v2](https://arxiv.org/html/2502.01433v2),
  Section V, studies the lower depth of decompositions through Toffoli
  gates. It does not supply the unrestricted implication (2).
- [Dutta et al., arXiv:2506.01542v1](https://arxiv.org/html/2506.01542v1),
  the proof of Theorem 1, uses ANF monomials, input-copy ancillas and
  measurement-based uncomputation. This ANF construction does not prove
  that arbitrary Clifford+T circuits admit an equal-depth Toffoli form,
  or establish an approximate occupied-dirty-work lower bound.

These are model checks, not new imported theorem premises. Do not equate
Toffoli depth or classical algebraic degree with unrestricted T-depth.
Bounded-gate light cones cannot ignore arbitrary Clifford interlayers,
and count divided by square-root width gives no logarithmic lower bound
here. Do not turn one encoded clean bit into independent initialized work
registers. No measurement, reset, supplied resource state, QRAM or free
coherent evaluator belongs to the target model.

## 7. Execution and integration

| Elapsed time | Work and required decision |
|---|---|
| 0–20 minutes | Record baseline and clock; reconstruct R75 exact return, R76 freezing/lifting, and R77's failed common-flip example. Fix the primary and adversarial contracts. |
| 20–90 minutes | Derive the one-encoder normal form and test the paired shared-S step. Run the short constructive scout independently. Apply the qualification gate. |
| 90–180 minutes | Complete a qualified gap proof and witness consequence, or one changed repair/native counterfamily. Prove all quantifiers and norm estimates. |
| 180–230 minutes | Freeze the supported claim. Obtain independent reconstruction and a fresh adversarial audit. Add only bounded fixtures needed for concrete remaining risks. |
| 230–300 minutes | Integrate the proof or precisely scoped outcome, verify the actual reviewed tree, render changed pages, publish, wait for CI and merge. |

Use separate construction, adversary and integration lanes. Fresh review
must inspect native words and the actual uniform estimate rather than
accept prior verdicts. A theorem needs a proof for unrestricted dirty
width; testing larger finite widths is not a substitute.

For substantive theorem/circuit changes, run the full checks in
[WORKSPACE](../../WORKSPACE.md#restore-and-verify) and the
[rendering checks](../../assets/README.md#rendering-checks). Preserve
receipt hashes. Publish the exact reviewed tree and merge only after CI
passes. Do not allocate a new result number until a supported result is
complete. The scientific work and its integration may finish early.

If neither a gap nor a native counterfamily closes, retain the exact
normal-form obstruction, tested hypothesis and concrete reopening condition.
That is a bounded-session outcome, not evidence that optimal general depth
is unattainable. Preserve the continuing destination and the distinction
between a constant-layer obstruction and a dimension-growing lower bound.
