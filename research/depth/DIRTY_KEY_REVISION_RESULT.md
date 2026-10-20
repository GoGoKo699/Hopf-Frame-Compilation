# Revision: pay the readout before claiming a shallower source

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Previous revision](OBSERVABLE_SOURCE_REVISION_RESULT.md)

This continuation starts from fetched main
`41e56213813c18fd0d490cea22562cf8576ec03f`, tree
`fe6fa06be9b7e28199d08b2c6f3c509ac94ab483`, after PR #212.
The main result is a more precise choice of unfinished operation,
supported by new paid components and scoped stopping theorems. No faster
complete indicator or stronger unrestricted depth bound is claimed.

## 1. Destination, theorems, components and approach limits

| Role | Assessment |
| --- | --- |
| General destination | Optimal native T-depth with simultaneous T-count, Clifford count, clean/dirty width, full operator error and return |
| Strongest general upper | R96 remains $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, with its literal one-clean allocations and complete initialized-isometry error $`2^{-L}`$. Selected R79 is unchanged |
| Strongest relevant lower results | The complete one-clean native-three sign ceiling remains $`D_{\rm np}=D_3+576`$, with retained original-circuit rounding floor $`2^{-530690}`$. The genuine original-pair branch has sign/XOR bounds six/five, or five/four at active rank at least four. These do not establish growing unrestricted depth |
| New paid source component | Amplitude-$`1/\sqrt2`$ recovery and dirty promotion now have an explicit general-k resource, accuracy and return recurrence, including a sufficient slowly growing clean window |
| New mechanism diagnosis | A final reverse query cannot repair a key-preserving core even with arbitrary quantum fibers. Other specified moving-key words retain a non-Pauli action or affine labels. These are operator obstructions for those representations |
| New constructive candidate | A Vandermonde encoding makes the selected output an inverse column of an affine matrix pencil. Its quadratic XOR query is paid; coherent inverse-column readout remains unpaid |
| Focused lower scout | Four bands permit arbitrarily many original physical pairs even in an identity word. A conditional radical certificate and literal pair-removal rules survive, but representation pair count cannot certify unrestricted depth |

The old number 37748865/37748864 was a generic feature ceiling, not a
prediction of an optimal degree. The retained structural reduction
already replaces it by six/five and high-rank five/four in its original
class. Further constant shaving is not the main question. The revision
asks whether a new native operation changes the complete resource law.

## 2. The full constructive contract remains the test

For a p-bit address and $`m=2^p`$, the destination is

```math
I_p|x,Y,w\rangle=|x,Y\mathbin\oplus e_x,w\rangle,
```

at every fixed requested tolerance tau, with zero clean work,
$`T,h=O_\tau(m)`$, $`G=O_\tau(m^2)`$ and
$`d=o(\log(p+2))`$. Accuracy is full-operator accuracy, including
arbitrary output, dirty-work and reference inputs. The actual inverse
and the completed return error are charged.

The [root-fusion reduction](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md)
allows polynomial factors in p in a smaller complete source. Its eight
source/actual-inverse calls, error factor eight and exact
$`O_\tau(\log\log(p+4))`$-depth tail are already paid.
The [observable contract](OBSERVABLE_TRANSPORT_SOURCE.md) remains

```math
\left\|J^\dagger Z_{Y_0}J-
 \sum_x|x\rangle\langle x|\otimes Z_{Y_x}\otimes I_H
\right\|_{\mathrm{op}}\le\delta.                 \qquad(1)
```

It permits arbitrary intermediate address motion and quantum helper
action. Chronological $`H_Y,J,Z_{Y_0},J^\dagger,H_Y`$ gives the
complete indicator with

```math
\begin{aligned}
T_I&\le2T_J,&d_I&\le2d_J,\\
G_I&\le2G_J+2m+2,&h_I&\le h_J,\\
\varepsilon_I&\le\delta,&a_{\mathrm{clean}}&=0.
\end{aligned}                                    \qquad(2)
```

Fusion uses sixteen J calls and error at most $`8\delta`$.
A full-source approximation eta implies $`\delta\le2\eta`$.
For fixed-L frame transfer the retained request is
$`\delta_L=2^{-L}/(128K_*)`$. Meeting (1) with sublogarithmic
depth and all simultaneous budgets would improve the fixed-accuracy
frontier. A new expression for J alone does not.

## 3. What the dirty-key revision actually excludes

The [key-mixing screens](KEY_MIXING_OBSERVABLE_SCREENS.md) make
three distinctions that the previous common-key theorem did not make.

1. With the top probe on B, a final reverse quadratic query leaves a
   fixed Fourier-key shift. The desired character has zero coefficient
   there. The resulting operator defect is at least $`\sqrt2`$,
   even for arbitrary quantum fibers in the key-preserving core.
2. With the top probe on A, Clifford fibers and an affine Pauli insertion
   yield affine block labels. The full monomial family exceeds that
   affine image; the stated dimension-dependent constant gap follows.
3. Moving the key inside the two-query source generally leaves a
   non-Pauli Clifford transition. When that transition is Pauli, the
   remaining labels are again affine. This particular grammar has
   error at least $`\sqrt{2-\sqrt2}`$.

Each statement names its actual boundary, fiber and insertion premises.
Additional key mixing and arbitrary intermediate address motion remain
outside those theorems. Two explicitly tested words outside parts of
the screened grammar also fail on an allowed column; those are exact
counterexamples to those words, not general theorems about all such
extensions. The proof home charges their actual query counts and depths.

A separate group-law observation explains why a fixed conjugated
translation family repeatedly recreates full-bank routing: exact
single-probe transport throughout that family forces translation of all
output probes. The next candidate therefore leaves that source family,
not merely its chosen Fourier basis.

## 4. The affine pencil exposes a different first unpaid operation

The [nilpotent-pencil source](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md)
uses $`\mathbb F_{2^p}`$ and the shift L on m field coordinates.
In a fixed binary Vandermonde basis E, put

```math
\begin{aligned}
N_x&=E^{-1}(I+xL)E,\\
N_x^{-1}a_0&=a_x,\\
J_x&=U_{N_x^{-\mathsf T}}.
\end{aligned}                                    \qquad(3)
```

The m physical output wires are specified basis coordinates; the
remaining $`(p-1)m`$ carrier wires are arbitrary private work.
Equation (3) gives the desired observable exactly. It is not an
XOR-translation representation. For p at least two, the field basis
also turns L into a star: for nonzero x the selected leaf has a zero
diagonal pivot, while x=0 gives identity. This replaces a long inverse recurrence by a unique-zero-pivot
readout question, without supplying that readout.

The quadratic query
$`(x,z,t)\mapsto(x,z,t+N_x^{\mathsf T}z)`$ is paid in eight
T layers, with $`T,h=O(mp^2)`$ and
$`G=O(m^2p^2+mp^3)`$, exact arbitrary-work return and an actual
inverse. These counts fit the relaxed source budgets.
The first unpaid operation is the in-place matrix action, a suitable
weak quantum-fiber completion, or the still weaker complete observable
readout. A forward XOR query does not implement any of these for free.

The inverse-column coordinate has Boolean degree p+1. Bounded-degree
completed monomial rounds therefore retain a logarithmic round cost.
Fixed-address Clifford stages built from the affine query have a
separate symplectic-degree limitation. Neither is a lower bound on
unrestricted native depth: genuine non-Clifford interference, address
motion and joint resynthesis remain available.

The paired-map and convolution reductions in the proof home also test
frontier relevance. The prescribed weak completion is equivalent to a
complete indicator under constant-depth, polynomial-address overhead.
Supplying the old indicator merely reconstructs the old depth order.
The one-probe promise is strictly weaker than that prescribed principal
action; no unsupported impossibility of converting between promises is
inferred. The research target is an actual new readout mechanism, not
another equivalent interface.

Two concrete zero-pivot attempts locate the next mechanism more sharply.
Signal-only spectral processing cannot distinguish the zero query from
a nonzero one on their common +1 bank eigenspace. A literal paid
address-rotation sandwich has depth at most seventeen, but its complete
observable is address-translation invariant and has error at least one.
The next word needs genuine nonmonomial carrier mixing across query
stages; merely breaking an address symmetry is insufficient evidence
of a correct readout. These tests do not exclude unrestricted native
interference.

## 5. Growing clean count is quantitatively open

[General-k recovery and promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#63-recovering-a-general-clean-bank-compression-with-one-extra-flag)
starts from the complete compression
$`\|K-Q_f/\sqrt2\|\le\epsilon`$, $`\epsilon\le1/8`$.
Three source/actual-inverse calls and a returned zero-test give complete
initialized-isometry error

```math
\eta=3\sqrt2\epsilon+2\epsilon^2.                \qquad(4)
```

One extra clean flag and the source's k clean wires are then promoted
to arbitrary dirty wires. The proof home records every zero-test,
literal global phase, Clifford, helper, inverse and return term. In
particular,

```math
\begin{aligned}
d_E&\le2^{k+1}
 (12d_U+16d_{C_k}+4d_{C_{k+1}}+16),\\
\varepsilon_E&\le16\,2^k
 (3\sqrt2\epsilon+2\epsilon^2).
\end{aligned}                                    \qquad(5)
```

Signal return is approximate in general. Exact borrowed-signal return
requires the separate flatness premise $`K^\dagger K=I/2`$;
one-signal rounding does not automatically establish it for growing k.
The request $`\epsilon\le\delta/(72\,2^k)`$ suffices.

For $`k=\max(1,\lfloor\beta\log_2\log_2(b+4)\rfloor)`$,
$`d_U=O((\log b)^\gamma)`$, $`\beta,\gamma\ge0`$ and
$`\beta+\gamma\lt1`$, the complete depth is sublogarithmic,
provided source accuracy and all polynomial resource envelopes are
uniform along this schedule. Thus growing clean count has a concrete
sufficient window. No source in that window has been constructed.
The bounded-clean two-band defect theorem does not close it.

The [finite-access filter audit](APPROXIMATE_INDICATOR_FILTER_LIMITS.md#41-finite-parity-access-still-fails-under-coherent-address-motion)
also permits coherent address/seed motion and arbitrary private
processing, while the bank is accessed only through a listed family of
parities. A deficient mask span gives relative compressed-source error
at least one. Small hash families therefore retain an almost-full
address in the paid primitive. This is an access-model limit; clean
count alone does not bound seed space, and bank-mixing operations
outside those masks remain available.

## 6. Preserve the focused lower scout without importing its conclusions

The [four-band original-pair scout](FOUR_BAND_ORIGINAL_PAIR_SCOUT.md)
keeps the actual signed native lists, physical clean code, unused clean
block, all dirty/reference inputs and literal adjoints. Under its stated
center and endpoint premises, pair differences lie in a coupling
radical, yielding a conditional rank certificate. Shared physical pairs
can be removed by a literal Clifford identity.

These are useful components, not a four-layer degree theorem. An exact
identity family has arbitrarily many original pairs. A second exact
Boolean completion has disjoint middle pairs with a proper endpoint
support; another whole-word removal still lowers its displayed depth.
For its fixed endpoints and a diagonal complete middle, occupied
Booleanity itself forces that removal, without imposing an unused-block
promise. A controlled-SWAP example then defeats endpoint commutation,
yet another exact coordinate change still gives two bands.

Under the further actual one-pair-per-list and equal-difference
hypotheses, the complete middle is a controlled product of r, CZ and
SWAP. Legitimate sector restrictions give a conditional occupied sign
degree at most four. The actual singleton lists commute across bands
and have equal binary sums. Cubic degree is equivalent to two remaining
coefficient identities, $`B_3=0`$ and $`C_2=D_2`$, in that proof
home's notation. Deriving them from the actual native lists, or finding
a literal counterexample, is the next focused scout. This is not a
general four-layer bound; a sector representation is not a free
full-space circuit.

## 7. Next question and evidence boundary

Seek a literal native word that reads the unique zero pivot, or directly
transports the output observable, with interference that survives the
fixed-address Clifford and completed-monomial boundaries. Prove its
complete operator action before treating a phase kernel, matrix inverse
or initialized column as an implementation. Then charge its actual
inverse, error, return and all resource budgets on that same word.

Keep the small physical-pair scout active in parallel. A scoped failure
can change the next mechanism without closing the unrestricted problem.
If a construction reconstructs the old logarithmic source, preserve its
useful component and change the missing operation rather than promoting
it to a frontier improvement.

The new exact fixtures check the specified key-mixing columns, one
general-k recovery instance, the finite-field pencil identities and
small complete physical-pair words. General recurrences and universal
scoped theorems are analytic. Finite checks and internal audits support
literal identities; they do not establish asymptotic claims or replace
external technical review.
