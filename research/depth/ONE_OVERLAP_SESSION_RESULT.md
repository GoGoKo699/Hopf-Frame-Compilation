# Overlap-one closure and the cost of the cycle source

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Source proof](ONE_OVERLAP_SOURCE_STRUCTURE.md) ·
[Cycle source](NONCOMMUTATIVE_CYCLE_CHARACTER.md) ·
[Previous checkpoint](THREE_BAND_NONLINEAR_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`5501fc63cd1418b00e0d0abc6d1c3248062a2186` (PR #176).
It resolves the preceding checkpoint's overlap-one source question and
reformulates the constructive source cost. It does not improve the
unconditional complete-frame upper law or prove a growing unrestricted
lower bound. This is a research checkpoint, not manuscript drafting.

## 1. Compact revision

The selected count package is settled under its stated hypotheses. For
integers $`n\ge1,L\ge6`$, $`N=2^n`$, the strongest unconditional
research upper law is R96:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)).
```

One clean qubit suffices, with real dirty width $`b\ge17(L+n+7)`$
or complex magnitude width $`b\ge17(L+n+8)`$, and full error at most
$`2^{-L}`$. The prescribed operator includes every logical column,
literal signs and phases, singular and zero coefficients, arbitrary dirty
inputs and external references, actual inverses, leakage and complete
workspace return. T-count, elementary Clifford count, T-depth, precision
and simultaneous workspace remain separate budgets. Every source, query,
control, recovery and inverse is charged. R79 remains selected Result D;
R97 remains a scoped four-layer obstruction.

The exact destination is optimal worst-case general T-depth of that
complete frame under these simultaneous budgets. At fixed accuracy and
square-root width, the remaining question includes whether R96's
logarithmic-in-n additive depth is necessary or can be reduced. R98 would
transfer a faster exact multi-output native read-only query to a better
complete-frame upper bound. The separate two-call corner reduction would
transfer an unrestricted query lower bound. Neither premise is proved.

Previously open questions already resolved include optimal one-clean
count, the all-precision upper law, zero-clean read-only five-layer degree
at most 14, the six-layer Pauli-comparison bound 51, and uniform transport
when Pauli-domain deficiency is bounded. The diagonal-C3 source theorem
already supplied six-layer degree bounds 78/26. This round expands the
native source class to which that theorem applies. Unrestricted one-clean
four-layer rigidity and dimension-growing lower bounds remain open.

## 2. The missing invariant addressed this round

The preceding source theorem bounded the number of middle axes that fail
to commute with the outer band by 2t, where
$`t=\dim(\mathrm{span}(A)\cap\mathrm{span}(B))`$ for the actual source
$`F=E_CE_BE_A`$. It did not control their coupling to a wide diagonal
part or to the inner band. Counting two exceptional axes at t equal to
one was therefore insufficient to claim a small residual source.

The [new authoritative proof](ONE_OVERLAP_SOURCE_STRUCTURE.md) establishes:
**a three-band C3 source with independent physical T axes and t at most
one is diagonal C3 in fixed Clifford coordinates.** No intermediate C3
membership is assumed. The proof handles all weights of the unique
middle-axis relation:

- Weight at least three makes every middle/outer axis pair commute by
  the retained activation-code argument.
- Weight one uses a new paired half-angle classification. Its only
  non-Pauli option has common four-term Pauli support, which contradicts
  independence of the selected inner axes.
- Weight two leaves one paired axis set in each outer band. The
  control-X image of the full C3 source is Clifford. Its symplectic
  involution condition eliminates every spectator component of the inner
  relation. The remaining paired two-qubit core has exact Clifford
  branches, including their eighth-root phases. Common full-space
  Clifford coordinates diagonalize it, and native right-extension
  restores the rest of the inner band.

The final restoration matters: remaining inner axes need not commute
with the entire core. Nor is separate diagonalization of the two control
sectors enough. The proof supplies common coordinates and retains the
literal phase polynomial of the controlled relative Clifford.

The existing six-layer detector now gives output degree at most **78**,
or **26** for a visible outer shear, below full error $`\sqrt2/8`$
for this source class. All source equalities hold on the complete work
space and tensor with any external reference. The coordinate changes are
analytical; no physical Clifford, T gate, source call or workspace is
removed from a circuit ledger. The additional C3 premise is exact and is
not inferred from approximate query return. No robust overlap-one source
certificate beyond the stated theorem is asserted.

The same proof home records a reusable extension lemma: a diagonal-C3
source followed by a commuting T-gadget list with exactly one binary
relation stays diagonal C3 if the full product is C3. The genuine
four-axis exception is handled through a preserved complete commuting
Pauli frame. It would invalidate a stronger claim that every individual
axis image must be Pauli. This is separate from treating dependent
commuting gadgets as one physical T band.

## 3. Constructive reformulation and diagnosis

The [cycle-character proof](NONCOMMUTATIVE_CYCLE_CHARACTER.md) now splits
its arbitrary coefficient register at a fixed cyclic cut. U is strictly
upper triangular and L is lower triangular, including the full-cycle
diagonal. Over the binary field,

```math
\chi(1+z)=\mathrm{tr}((I+U)^{-1}L).
```

Hadamards on all $`r=d(d+1)/2`$ L wires turn the phase source into the
complete XOR word query for the corresponding entries of the
unitriangular inverse. T-count, T-depth and additional workspace are
unchanged; elementary Clifford count increases by $`2r`$. The converse
uses the same basis changes. This is an all-output, arbitrary-input
identity, including actual inverse and reference return. Native read-only
access to U is inherited only when the source implementation actually
provides it; the earlier in-place factorization does not establish that
interface.

A single inverse entry also has a two-source dirty echo with two CNOTs
and two Hadamards. Separately, a fully controlled cycle phase can be built
from two uncontrolled source calls and a returned auxiliary lower array.
For a zero-clean source of T-count $`T_D`$, depth $`d_D`$, dirty
workspace $`W_D`$ and error $`\delta_D`$, the wrapper has

```math
T\le2T_D+16r,\qquad
G\le G_D+G_{D^\dagger}+40r,\qquad
D_T\le2d_D+16,
```

extra arbitrary workspace $`r+\max(W_D,r)`$, no clean work and
error at most $`2\delta_D`$. The actual inverse's Clifford word is
charged separately. The auxiliary array and parallel Toffoli helpers
return on every input; source and helper pools are reused only at
disjoint times.

This reformulates the obstacle rather than removing it. The earlier
Gaussian inverse-extraction attempt had already shown why a cheap shear
exposes forward multiplication, not inverse coefficients. The cycle
identity now exposes the entire inverse word as the source's actual
computational content. The retained cycle circuit remains slower than
fairly simplified conjunction baselines. A lower bound based only on
commutative coefficient width already failed under this noncommutative
representation; no such parameter is promoted to an unrestricted depth
invariant.

## 4. Most useful resumption question

**For a specified three-native-band C3 source with overlap two, is signed
transport of the actual six-layer inner family bounded independently of
dirty width, despite an arbitrarily wide diagonal source part?**

First test whether diagonal-C3 classification survives. If it fails,
retain the counterexample and seek bounded signed transport directly;
non-diagonal structure alone is not a returned-query counterexample.
The new proof identifies the precise missing step: several independent
relation directions can couple controlled cores, and no common
full-space reduction analogous to the single-control involution argument
has been proved. A proof or a genuinely returned high-degree native
counterexample would advance a degree-versus-depth recurrence.

This is chosen over another unpriced source identity because the first
positive-overlap case now closes analytically, giving a concrete next
structural test. On the constructive side, the useful question is whether
the complete inverse-word interface admits a better native schedule with
all sources, controls, inverses and simultaneous workspace paid, and
whether that structured schedule can reach arbitrary-table R98. Neither
reduction is automatic. Overlap is representation-dependent, and even
one wide T band can have unbounded Pauli deficiency. Both parameters must
be stress-tested under alternative native syntheses.

General non-C3 prefixes, one-clean work, address mixing and the growing
unrestricted lower bound remain distinct gaps. The project continues;
closing this candidate does not settle optimal general depth.

## 5. Evidence and boundaries

The [source checks](../../tests/test_one_overlap_source_structure.py)
exercise the paired exception, the spectator block identity, exact
oriented core phases, both adjacent-partial failures and the dependent
four-axis frame. The [cycle checks](../../tests/test_noncommutative_cycle_character.py)
exercise the inverse identity, complete word query and dirty control
wrappers. Finite fixtures support the explicit identities; the arbitrary
width theorems follow from the proofs in their authoritative homes.
Repository validation and rendering receipts are recorded in the merging
pull request. No manuscript claim selection changes in this round.
