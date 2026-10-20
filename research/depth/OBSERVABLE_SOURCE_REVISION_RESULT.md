# Revision: keep the destination, change the source mechanism

[Workspace](../../WORKSPACE.md) · [Full frontier](../../docs/OPEN_PROBLEM.md) ·
[Observable-source audit](OBSERVABLE_TRANSPORT_SOURCE.md) ·
[Previous revision](CENTRAL_PARITY_REVISION_RESULT.md)

This continuation starts from fetched main
`db5cd441e8894137fccae015418a8713633cae46`, tree
`157fcd31f092eae03cbf19700ed639efb5d54951`, after PR #211.
The revision separates a useful reduction from an improvement in the
complete resource frontier. No faster source or stronger unrestricted
depth bound is claimed.

## 1. What we have learned

| Role | Current assessment |
| --- | --- |
| General destination | Optimal native T-depth with simultaneous T-count, Clifford count, clean/dirty width, full operator error and return |
| Strongest general upper | R96 remains $`D_T=O(NL/b^2+L+\log_2(n+2))`$, $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, with the same one-clean allocations and complete error $`2^{-L}`$ |
| Strongest relevant lower statements | Exact at-most-one-clean native-three sign degree is at most $`D_{\rm np}=D_3+576`$; the retained original-circuit rounding floor is $`2^{-530690}`$. The genuine original-pair branch has sign/XOR bounds six/five, or five/four at common active rank at least four. These do not establish growing unrestricted depth |
| Useful source components | Complete quadratic XOR queries, algebra-key multiplication, actual-inverse readout, quantum helper-fiber cancellation and polynomial-overhead root fusion are paid |
| New source diagnosis | Any word preserving a common Fourier key with Clifford fibers and affine address-Pauli insertion retains affine address labels. More such overlaps cannot produce the required monomials at error below $`\sqrt2`$ |
| Candidate outside that representation | The opposite-direction two-query word leaves explicit coherent dirty-key dependence. Its completed phase has error two on an allowed input, despite exact work return |
| Lower diagnostic | One exceptional endpoint relation factors into a diagonal phase, fixed Clifford and another quadratic-controlled Pauli. Its polar rank can grow with the relation support, so the old mechanism survives |

The large degree constants are proof ceilings, not predicted optimal
degrees. Improving another ceiling or wrapper constant alone would
leave the full upper/lower resource gap unchanged. The lower diagnostic
therefore remains a focused component, while the main effort returns
to constructing the missing shallow complete operation.

## 2. Require the observable, not a prescribed coefficient permutation

The retained [selected-routing reduction](JOINT_ROUTING_INDICATOR_INTERFACES.md)
already permits an arbitrary source whose conjugated probe is the
complete indicator. Applied directly to the current output bank, the
sufficient contract is

```math
\left\|J^\dagger Z_{Y_0}J-
 \sum_x|x\rangle\langle x|\otimes Z_{Y_x}\otimes I_H
\right\|_{\mathrm{op}}\le\delta.                    \qquad(1)
```

J can move the address and apply arbitrary quantum action to its private
dirty helpers. It need not implement F or any specified classical map.
This is a synthesis freedom within an already paid query reformulation,
not a new asymptotically easier task. Certifying (1), uniformly on all
inputs, is the first unpaid operation.

The actual chronological word is $`H_Y,J,Z_{Y_0},J^\dagger,H_Y`$.
For common forward/actual-inverse source bounds, its complete ledger is

```math
\begin{aligned}
T_I&\le2T_J,&D_T(I)&\le2d_J,\\
G_I&\le2G_J+2m+2,&h_I&\le h_J,\\
\varepsilon_I&\le\delta,&a_{\mathrm{clean}}&=0,
\qquad m=2^b.
\end{aligned}                                      \qquad(2)
```

There is no extra probe wire. Equation (1) includes address, output,
helper and reference return after the completed word; no intermediate
return promise is assumed. If only full-source distance eta to an
ideal J is proved, charge $`\delta\le2\eta`$.

Root fusion uses sixteen J/actual-adjoint calls, error at most
$`8\delta`$ and depth
$`16d_J(b_0+1)+O_\tau(\log\log(p+4))`$; request
$`\delta=\tau/8`$. The fixed-L frame transfer requests
$`\delta_L=2^{-L}/(128K_*)`$. T-count and dirty width must be
$`O_\delta(2^b(b+2)^c)`$, Clifford count must be
$`O_\delta(4^b(b+2)^{c_G})`$, and depth must be
$`o(\log(b+2))`$, on the same circuit with fixed exponents.
A source meeting all these conditions would change the fixed-accuracy
frontier through R98; naming (1) does not.

For the retained [algebra source](OUTPUT_ALIASED_CHARACTER_SOURCE.md),
the top carrier can be the constant subset-zeta output coordinate.
This removes two redundant Z gates and one extra carrier wire:

```math
\begin{aligned}
T_I&\le4T_V,&D_T(I)&\le4d_V,\\
G_I&\le4G_V+bm+2b+2m+2,&h_I&\le m-1+h_V,\\
\varepsilon_I&\le4\eta,&a_{\mathrm{clean}}&=0.
\end{aligned}                                      \qquad(3)
```

Fusion still uses 32 V calls and error $`32\eta`$; the known
logarithmic-depth source retains the old depth order. The independent
one-signal amplitude-$`1/\sqrt2`$ route likewise retains
$`24d_U+32`$ depth and error
$`16(3\sqrt2\epsilon+2\epsilon^2)`$, with only approximate
signal/toggle return unless flatness is established separately.

## 3. Which obstruction applies where

The [shared-key theorem](OBSERVABLE_TRANSPORT_SOURCE.md) concerns an
actual source decomposition into Clifford fibers over one unchanged
Fourier key, followed by an affine address-Pauli insertion and its actual
inverse. It covers every word length in that grammar. It does not cover
arbitrary quantum fibers, key mixing, intermediate address motion or
general native circuits. Logical source calls are not native T layers.

The tested word with both quadratic-query directions genuinely leaves
that premise. Its failure is instead a literal unwanted dependence on
an arbitrary auxiliary bank. Correcting that dependence is unpaid;
assuming a zero auxiliary bank would change the operator contract.
Neither failure is an obstruction for the unrestricted model.

The [lower relation factorization](CENTRAL_PARITY_QUADRATIC_REDUCTION.md#7-the-exceptional-relation-recreates-a-quadratic-controlled-pauli)
is also local. Its new control parities commute with the endpoint band,
but are not proved to commute with the actual middle or opposite
endpoint. Moving or canceling the quadratic-controlled Pauli is the
first unsupported inference. An endpoint normal form cannot supply
that whole-word identity for free.

## 4. Chosen next question and stopping rule

Seek a complete observable word with genuinely nonmonomial interference
between the two key directions. Its first target is independence of the
arbitrary auxiliary bank while retaining the full monomial character.
Prove that operator identity before estimating a native schedule; then
charge the same word's depth, counts, helper use, actual inverse and
approximate return. Direct address motion remains allowed.

Merely adding another same-key overlap, exchanging which bank carries
the unchanged key, or rebuilding the known in-place transform does not
answer that question. The first two are screened here; the last retains
logarithmic depth. A different mechanism is welcome, but must satisfy
the full contract rather than an initialized column or trace identity.

The exact new output-probe fixture checks all 1024 columns and the actual
reverse, with coherent nonreturning source-helper action. The previous
2048-column fixture remains intact. The general shared-key theorem,
dirty-key counterexample and lower relation factorization are analytic;
these finite receipts do not prove their universal scope. Internal
audits do not replace external technical review.
