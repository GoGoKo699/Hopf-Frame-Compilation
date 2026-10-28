# Independent audit of the coupled projective scratch decoder

Verdict: **PASS**, under the finite p=3 contract below. This audit is outside the
frozen PR269 tree; no repository files were edited.

Audited source:
`../coupled_projective.py`, using the separately verified private-mask
`displacement_batch_word()` from this directory.

## Literal corrected reflection identities

All four-wire Pauli labels use x+16z in the Hermitian section. Define

    JA = J(19,50;32),
    JB64 = J(156,228;64),   JB112 = J(156,228;112),
    JC44 = J(89,247;44),    JC87 = J(89,247;87),

where J(A,B;P)=T(A)T(B) P T(B)^dagger T(A)^dagger. Every pair A,B commutes
and is independent, and P anticommutes with both. Hence these are literal
Hermitian Clifford involutions. Independent field matrices give

    F4 = JB112 JA,
    F5 = JA JB64,
    F6 = JA JC44,
    F7 = -JC87 JA,

where FD=CD C1^dagger. In particular, changing the earlier C-reflection probe
on line6 from1 to44 allows the same A-reflection probe32 in every row; this
is a true signed operator identity, not only a symplectic equality.

The audit script reconstructed CD from the F4/F16 arithmetic and partial
trace Fourier definition, not from the proposed native factorization. It
compared every entry of all four identities using exact Gaussian dyadics.

## Why the common A selection has the correct phase

Let fD indicate that U belongs to line D, including U=0, and put

    g = u2 OR u3.

For U!=0, g is one precisely on the four nonbaseline lines, and exactly one
fD is then one. On the baseline line all four fD and g vanish. At U=0 all
four fD equal one but g=0.

The common A word's inner operation is

    (iP)^u2 (iP)^u3 (-iP)^(u2 u3) = i^g P^g,
    P=Z_(V1).

Conjugating by the shared two-T preparation gives i^g JA^g exactly. The
four remaining selected B/C reflections contribute (-i)^(sum_D fD), and
the selected line7 CZ contributes (-1)^f7. Thus on every nonzero U,

    i^g (-i)^(sum_D fD) = 1.

The remaining negative sign on line7 is precisely the minus in F7. The
chronological pattern is early B5,C6; common A; late B4,C7. On each nonzero
line this selects exactly its required two-reflection order.

At U=0 the common A is absent, the four relative phases multiply to1,
and the explicit line7 phase leaves the fixed zero-cell Clifford

    -JC87 JB112 JC44 JB64 C1.

This is a permitted whole-sector action. Source queries are identity on U=0;
after the exterior coordinate sum the intermediate pairs have a=c, where
both lower translations and the half swap are identity. No initialized input,
ignored zero column, postselection, or unreturned flag is used.

## Displacement composition and arbitrary helpers

The six displacements can be placed after the coupled linear word. For
nonzero U only one line contributes; at U=0 all six vanish. The private-mask
batch implements their entire XOR map on V and identity on its six arbitrary
helpers, as proved in `QUARTIC_BATCH_HANDOFF.md`.

The batch proof covers arbitrary entanglement: its native phase core preserves
the principal basis labels and has scalar identity action on each private
helper. The full helper action is therefore the tensor product identity on
all six helpers. Two target-Hadamard walls produce the requested vector XOR
without changing helper return. No dense 14-wire matrix is needed for this
complete operator identity.

## Independent verification performed

`audit_coupled.py` checks:

- All premises of the five displayed paired reflections and all four literal
  field correction identities.
- Every256column of the complete coupled linear word, including its newly
  chosen zero block and its literal phase, against independently assembled
  field/Pauli blocks.
- Every256column of the actual native linear inverse against that adjoint.
- Independent whole eight-wire recoupler assembly from field Fourier matrices,
  the field shift, and the coordinate permutation; equality with the prototype's
  effective logical matrix.
- The bottom E2 word through the separate CQ2 native evaluator, then all256
  columns of unitarity and all2,048 address/column intertwining comparisons.
- Distinct physical T targets in every claimed band of each complete forward
  and actual-inverse native word.

The quartic component separately passed3,072 exact2x2 helper-sector matrix
comparisons: all256principal inputs, all six private helpers, and both actual
orientations, including both basis columns of each helper. Its formal quartic
primitive also passed every five-wire column and inverse.

All matrices in this audit have Gaussian dyadic entries. Individual paired
reflection matrices have denominator2, so no floating irrational is introduced.
The full field Fourier matrices also have power-of-two denominators. The audit
uses exactly representable bounded dyadic arithmetic and equality checks; it does
not use a numerical tolerance or identify global phases projectively.

The complete literal construction has the following distinct allocation:

| Word | T count | T bands | Forward gates | Actual-inverse gates | Physical wires | Extra dirty helpers |
|---|---:|---:|---:|---:|---:|---:|
| Coupled recoupler |136|36|980|984|14|6|
| Full eight-address decoder |148|39|1,110|1,166|14|6|

There are zero clean inputs and zero approximation error. The actual complete
inverse reverses every native block, including the free zero Clifford and all
helper words. Its elementary-gate count differs because S-dagger is emitted
with positive S powers; the table retains this difference.

Machine-readable receipt: `coupled_independent_audit.json`.

## Scope

This is a finite coupled native construction, with a width/depth tradeoff against
the frozen zero-extra-helper172T/83band decoder. The common selector g is only a
two-bit OR at p=3. At general field size its analogue is a growing-width nonzero
predicate, and the projective-line count remains exponential. Neither that
predicate nor a coupled implementation across all recursive levels is paid here.
The general R96 frontier is unchanged.
