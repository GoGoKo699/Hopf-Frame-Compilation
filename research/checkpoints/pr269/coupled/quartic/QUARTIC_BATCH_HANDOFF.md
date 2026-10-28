# Scratch: private-mask quartic batch for the projective displacement

This is outside the frozen PR269 tree. No repository file was edited.

## Existing primitives and what is being assembled

`research/depth/DIRTY_CUBIC_PHASE_COMPILER.md`, section1, already owns the
R91 four-phase private-mask word

    F(a,b;h)=(-i)^(ab) Z_h^(ab).

It costs4T in4bands, with all T gates on the arbitrary private h and four CNOT
parity moves. `AFFINE_ORBIT_REFLECTION_SOURCE.md`, section2, and
`PARALLEL_CLEAN_PROMOTION.md`, section2, already pay a32T/32band quartic batch
with two private helpers per term. The frozen `PROJECTIVE_FIELD_RECOUPLER.md`
section3 has the16T/8band single-quartic commutator, but its compact relative
Toffolis place T gates on shared principal controls, so its8band schedule does
not directly parallelize arbitrary overlapping terms.

The following assembly uses the retained R91 private-mask word on both sides
of that commutator. It gives16T/16bands per quartic, and those16bands batch
across arbitrarily many overlapping terms with one private helper per term.
No new catalytic primitive or growing-level depth theorem is claimed.

## Complete identity, including affine forms

Let a,b,c,d be arbitrary affine binary parities of principal bits. They may
share inputs and need not be independent. None contains the private helper h.
For a=(linear_a,kappa_a), b=(linear_b,kappa_b), visit the four parity vertices
0,a_linear,a_linear+b_linear,b_linear of h and return to0. At vertex indexed
by subset J of the two forms, use T with sign

    (-1)^(|J| + sum_(j in J) kappa_j).

No X is required to implement an affine constant. The emitted word is

    Ftilde(a,b;h)=i^(kappa_a kappa_b) (-i)^(ab) Z_h^(ab).

The omitted finite-difference constant accounts for the displayed literal
factor i^(kappa_a kappa_b); it is not discarded from the identity. Set

    U=Ftilde(a,b;h),
    V=H_h Ftilde(c,d;h) H_h.

The chronological word U,V,U^dagger,V^dagger has exact action

    V^dagger U^dagger VU=(-1)^(abcd) I_h.

All scalar factors cancel against their actual inverses. The remaining sign
is the Pauli X/Z commutator. This is an operator identity for arbitrary helper
inputs, including entanglement with the principal system or external references.

For Q quartics, assign disjoint private h_j. At each of the four Gray vertices,
perform all parity CNOT updates and then the T/T-dagger gates on distinct h_j.
Run corresponding U,V,U^dagger,V^dagger stages together. Only principal-to-helper
CNOTs occur; no principal is changed inside the phase batch. Shared controls
therefore do not interfere with the tensor product of helper actions.

With L the sum of all four linear support sizes over all Q terms, the common
forward/actual-inverse ledger is

    T=16Q,  D_T<=16,  G_all=20Q+4L,  dirty helpers=Q,  clean inputs=0.

The4Q Hadamards belong to the two V stages. Affine constants change only the
T orientations and add no gates. For a single quartet of distinct physical
controls this is16T,16bands,36totalgates, one returned arbitrary helper.

## Exact six-term projective displacement

The displacement maps from distinct lines can be moved after every selected
linear correction. For nonzero U exactly one line is active, so there is no
cross-line interference. At U=0 every displacement vanishes. This reordering
therefore retains the entire operator, including the free zero fiber.

After H on the four V wires, each selected displacement is a quartic phase.
Forms are encoded as (linear bitmask, constant), with U on bits0..3 and Fourier
V on bits4..7. The six exact quartets are:

1. (15,1), (5,1), (4,0), (16,0)
2. (15,1), (5,1), (8,0), (32,0)
3. (14,1), (7,1), (13,0), (144,0)
4. (14,1), (7,1), (6,0), (224,0)
5. (6,1), (11,1), (15,0), (160,0)
6. (6,1), (11,1), (5,0), (240,0)

Here L=60. The phase core costs96T,16bands,360totalgates on14wires. Its two
external V-Hadamard walls give the complete XOR displacement word:

    96T,16bands,368totalgates, six additional returned dirty helpers.

The inverse has the same ledger. This replaces the frozen displacement
component's96T/48bands/noextrahelper allocation with a different explicit
width/depth tradeoff. It does not change the frozen implementation.

Import `displacement_batch_word()` from `private_quartic_batch.py`; it returns
(word,metadata), acting on U0..3,V4..7 and helpers8..13. The metadata includes
forms, helper locations, and forward/inverse ledgers.

## Verification and structural scope

`private_quartic_batch.py` passed:

- Every32column of the formal five-wire quartic and its actual inverse.
- 3,072 exact2x2 helper-sector matrix comparisons: all256principal bitstrings,
  all six helpers, both actual directions. Each comparison includes both helper
  basis inputs, so the full helper action is scalar identity, not merely a
  selected initialized column. Tensor independence gives the complete14wire
  phase identity without constructing a dense16384square matrix.
- Equality of the six quartic exponents and the field displacement for all
  16U values and16Fourier-V values.
- Exactly16native T bands, distinct physical targets within every band,
  and all T targets private.

The receipt is `batch_receipt.json`; it reports the360gate phase core. The
function metadata includes the368gate XOR word.

This assembly is uniform for a specified list of quartic affine-parity
products, and therefore pays any cubic XOR oracle supplied in that factored
form. It does not compress the exponentially many general projective lines,
or pay a growing-degree line predicate at fixed depth. At general recursion
level, a line predicate contains dim(F) binary zero conditions; that increasing
arity remains a separate synthesis obligation.

Two simple tradeoffs are available but not needed by the current coupled word.
Expanding the p3 displacement gives four quadratic and four cubic vector terms;
separate cubic/quartic phase batches still cost96T but24bands and8helpers, worse
here. Sharing one commutator helper across each line's two scalar components
gives72T/24bands/3helpers: U costs4T and V costs8T per line, with the three lines
batched. This last schedule has an analytic ledger but has not received the
same executable check as the96T/16band version above.
