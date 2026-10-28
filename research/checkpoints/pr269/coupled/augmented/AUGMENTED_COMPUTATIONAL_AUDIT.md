# Independent audit of the augmented uniform recoupler

Continuation only; not part of PR269 or its qualification gate. No existing
receipt-bound continuation file was modified by this audit.

Audited proof: `AUGMENTED_RECOUPLER_CANDIDATE.md`, SHA256
`17843000dd95be798fc8155161b0882bbf160e8099a31eaef29c6a67eab65309`.

## Verdict and scope

The complete augmented cell identity, shared-helper induction, and stated
native composition bounds are valid. The construction replaces the unpaid
norm division with explicit multiplication-XOR batches by changing the source
to include a coherent Pauli query on arbitrary helpers. Its depth remains
16(p-1) for the decoder and at most 64(p-1)+16 for the complete direct indicator.
This is O(p), so it does not reach the required o(log p) source depth or improve
R96. No all-size depth compression is inferred from helper sharing.

## Complete cell calculation

After the original H_K Fourier wall, the source is

    chi_K((lambda+b theta)UV) chi_F(lambda T).

The first addition to T uses f=Tr_(K/F)(UV)+d0 d1. Conjugating the helper
Pauli through that addition contributes chi_F(lambda f), cancelling exactly
the lambda part of the field phase. The remaining operator is

    chi_K(b theta U V) chi_F(lambda d0 d1) P_lambda.

The Fourier kernel chi_K(theta aV) is symmetric and nondegenerate, so the
corresponding transform is a literal self-adjoint Clifford involution. It
sends the b phase to translation a→a+bU, with no scalar phase. The second
addition contributes chi_F(lambda B(U,a)). It commutes with that translation
because B(U,bU)=0. Therefore the resulting phase is

    chi_F(lambda[B(U,a)+d0 d1]) P_lambda,

with principal translation a→a+bU and no residual helper translation.
Replacing (U,a) by (a,a+U) gives the child quadratic phases, since
q0(a)+q0(a+U)=q0(U)+B(U,a); the b translation becomes the core half swap.

At U=0 both additions are identity. The core source and target actions are
identity on their complete corresponding zero sectors, but the helper query
P_lambda is still present on both sides. The specified Fourier map is allowed
on every V input in this sector. No nonzero promise, initialized helper, or
discarded column is used.

## Shared helpers and literal covariance

Each complete ideal recoupler commutes with every helper X translation,
because it uses a helper only through additions to it. At fixed joint helper
X character, its action is a Clifford on its principal core: the additions
become quadratic signs and all remaining operations are fixed Cliffords.

Two child decoders on disjoint cores therefore commute even when they share
every lower-level helper. Conjugating by the core-only half swap exchanges
the identical child operators; their product is unchanged. This proves the
required swap commutation without swapping helper banks.

For the lower query, the same single Pauli helper factor is passed through
the children successively:

    Q_L Q_R P -> W_L Q_R P -> W_L W_R P.

Each step uses one child's full identity. The other core query or already
decoded W operator commutes with that child. P is not duplicated or squared.
The upper helper is disjoint from the children and remains unchanged by them.
This verifies the induction and the exact cocycle orientation

    E_(r+Lambda x) Q_x = W_x E_r.

## Native ledger and simultaneous width

The first addition indeed has the two-product form

    Tr_(K/F)(UV)+d0 d1 = d0(v1+d1)+d1(v0+v1).

The second has the two products d0 a1+d1 a0. The linear sums in the first
formula can be evaluated directly in private-mask parity updates; no clean
register containing a sum is required. Their supports are at most 2ell, so
only a constant in the Clifford count changes.

For s_j interpolation terms per product and 2^(p-j) cells, each of the two
additions uses at most 2*2^(p-j)*s_j private masks. All masks have distinct T
targets, and every logical parity update is Clifford. The common target
helper T_j may be shared: one Hadamard wall surrounds the whole sum of product
phases, whose contributions commute. Complete batches return their private
masks, so the same pool is reusable for the next addition and next level.

The claimed bounds follow:

    T <= 32 sum_j 2^(p-j) s_j = O(m p^3),
    D_T <= 16(p-1),
    private masks <= max_j 2*2^(p-j)*s_j = O(m p^2),
    retained moving helpers = sum_j 2^(j-2)=m/2-1,
    elementary gates = O(m^2 p^2).

The last bound includes fixed field-basis changes and Fourier matrices.
The coherent helper Pauli query is an explicit bilinear address/helper sign,
implemented by O(mp) CZ gates and zero T gates. The retained helpers are not
returned by the decoder or source individually. They are returned by the
actual source/probe/adjoint echo, along with every private multiplication mask.
All inverse calls are charged.

## Independent computational-basis verification

`augmented_computational_audit.py` independently implements the algebraic
addition, field Fourier, and CNOT primitives in the computational basis of
both principal and helper registers. It does not call the author's sector
fixture or condition on a prepared helper X eigenstate. Every amplitude is an
unbounded Python integer with one exact power-of-sqrt(2) normalization.

For p=2 it checks all 32 columns of the complete five-qubit decoder. For p=3
it checks all 2,048 columns of the complete eleven-qubit decoder. Every
address is then tested in the literal equation

    E(Q_x tensor P_p(x))=(W_x tensor P_p(x))E,

including the input helper sign and the possibly different output helper
sign. There are 128 and 16,384 address/column comparisons respectively.
Every column of the reversed algebraic primitive word is also compared with
the exact transpose of the forward matrix. All column norms are checked.
Every zero-U column is included. Both sizes pass.

This verifies the full helper-space identity and its coherent phases. It is
an algebraic primitive-word check, not a new physical native emitter for the
nonlinear additions. The all-size native costs remain justified by composing
the retained multiplication-XOR compiler, whose arbitrary-mask return contract
is required at each call. The separate native Clifford-sector fixture is
consistent with this independent computational-basis result.

## Source-hypothesis boundary

The source here is Q_x tensor P_p(x), not a field query tensor identity.
A helper Pauli label can be fixed while its operator acquires a minus sign.
Consequently a theorem whose source hypotheses include literal positive action
on every fixed label cannot be applied unchanged. No such import is used in
this construction or audit. Any extension of those obstruction theorems needs
its own retained-phase proof. This scope point does not impair the direct
cell identity or the full readout calculation above.
