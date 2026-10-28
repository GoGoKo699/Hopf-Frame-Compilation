# A uniform augmented recoupler with a moving Pauli helper query

Analytic continuation candidate after PR269. This identity supplies an explicit
all-size alternative to norm division in one recursive cell. It retains a paid
stage at every recursive level and does **not** improve R96 or the required
sublogarithmic indicator-source depth. The finite Clifford-sector fixture is
separate from the frozen PR269 qualification.

## 1. Complete cell identity, including zero

Use F=K_(p-1), K=F[theta], theta^2+theta+delta=0, |F|=q and
ell=[F:F2]=2^(p-2). Write U=d0+d1 theta. Define

    q0(U)=d0*d1,
    B(U,a)=d0*a1+d1*a0,       a=a0+a1*theta.

Let chi_F(z)=(-1)^Tr_(F/F2)(z), and similarly for K. Introduce an **arbitrary**
ell-qubit register T representing an element of F. Its helper query is

    P_lambda |T> = chi_F(lambda*T) |T>.

This is a Pauli; it is not an initialized register and need not remain fixed
under the source. In the Fourier coordinates of the original V register, the
field source query is chi_K((lambda+b theta) U V), where b is a bit.

Let A_f add f(U,V) into T, and use the fixed full-field Fourier transform

    H_(K,theta)|V> = |K|^(-1/2) sum_a chi_K(theta*a*V) |a>.

It is a literal Clifford involution. Define the chronological cell

    A_[Tr_(K/F)(UV)+q0(U)],     H_(K,theta),     A_[B(U,a)].

Its full operator C obeys

    C [chi_K((lambda+b theta)UV) P_lambda] C^dagger
      = chi_F(lambda*(B(U,a)+q0(U))) X_a^(bU) P_lambda.       (1)

Proof: conjugating P_lambda through the first addition contributes
chi_F(lambda*(Tr_(K/F)(UV)+q0)). Its first term cancels the source's
lambda UV phase, leaving chi_K(b theta U V) chi_F(lambda q0).
The Fourier transform sends chi_K(b theta U V) to X_a^(bU). The final addition
dresses P_lambda by chi_F(lambda B(U,a)). Its commutator with X_a^(bU) is
an addition of B(U,bU)=0 to T, so no helper translation or extra phase remains.
This proves (1) literally, for all lambda in F, both b, all U,V,T, and references.
For U=0 both additions are identity, the source/target principal queries are
identity, and H_(K,theta) is a specified unitary on the entire V sector. There is
no nonzero promise, norm inverse, zero flag, postselection, or discarded column.

Place the original H_K wall before this cell. Afterwards replace (U,a) by
(a,a+U), then apply the original H_F walls on the two w coordinates. Exactly as
in the retained recursive proof, the phase in (1) becomes the two lower field
queries, while X_a^U becomes the half swap. Thus the augmented recoupler Rhat_p
satisfies

    Rhat_p [Q_p(x',b) P_(beta_(p-1)(x'))] Rhat_p^dagger
      = (Q_(p-1)(x') tensor Q_(p-1)(x')) S^b
        P_(beta_(p-1)(x')).                                (2)

The helper Pauli is the **same operator** on the two sides. The cell itself
may act nontrivially on T. Its actual inverse is the reversed word of the two
self-inverse additions and Fourier/coordinate walls.

## 2. Why one helper per level can be shared by every cell

The ideal recoupler uses its helper T only through translations X_T^f.
Consequently it commutes with every individual X_T translation. Its action is
block diagonal in the common helper X basis. At fixed helper X character it is
a Clifford on its own principal core: the two added phases are Boolean
quadratic, and all other walls are Clifford.

At recursive level j allocate **one** register T_j of ell_j=2^(j-2) qubits,
shared across all 2^(p-j) cells at that level. Distinct cells have disjoint
principal cores and share only X-diagonal helper controls, so their ideal
operators commute. Two identical child decoders sharing all lower-level
helpers therefore have a product that commutes with the core-only half swap.
No controlled swap of a helper bank is needed.

The shared Pauli is not duplicated. If E_L and E_R obey the child identity,
then applying their conjugations successively to

    Q_L Q_R P_shared

gives W_L Q_R P_shared and then W_L W_R P_shared. Each conjugation leaves the
same single helper operator in place for the other child. Disjoint-core target
operators commute with the other child decoder. This establishes the induction
without an erroneous extra copy or square of P_shared.

Put

    P_p(x) = product_(j=2)^p Z_(T_j)^[beta_(j-1)(x mod 2^(j-1))],

using the trace-pairing convention for each Z. The recursive source on the
principal m=2^p wires and the shared helpers obeys

    Ehat_p [Q_p(x) tensor P_p(x)] Ehat_p^dagger
       = W_x tensor P_p(x).                              (3)

The base p=1 is the retained Clifford decoder. The total retained helper
width is sum_j ell_j=m/2-1. All these inputs are arbitrary. Each ideal decoder
commutes with every helper X translation, which preserves the sharing premise
at the next stage.

## 3. Paying the native primitives

The first addition factors into two F multiplication-XOR terms:

    Tr_(K/F)(UV)+q0 = d0*(v1+d1)+d1*(v0+v1).

The second uses two:

    B(U,a)=d0*a1+d1*a0.

Reuse the multiplication-XOR compiler in
`research/depth/RECURSIVE_FIELD_SHEAR_DECODER.md`, section5. Its interpolation
forms implement one F product with at most

    s_j=(2 ell_j-1) max(1,ceil(log2(2 ell_j-1)))^2

cubic parity terms, eight native bands, 8s_j T gates and s_j arbitrary returned
private masks. All terms of one addition, across every cell at that level,
can share the same eight-band schedule: they have private T targets, and all
logical parity updates are Clifford. Shared use of the target T_j is harmless;
its Fourier walls are applied once around the whole summed addition. The
first addition evaluates the linear forms v1+d1 and v0+v1 directly in the
private-mask parity updates; it does not load either sum into a clean register.
These forms have support at most2ell rather than ell, changing only a constant
in the elementary Clifford bound.

Thus one entire recursive level costs at most16 native bands. A conservative
complete decoder ledger is

    D_T(Ehat_p) <= 16(p-1),
    T(Ehat_p) <= 32 sum_(j=2)^p 2^(p-j) s_j = O(m p^3),
    private returned masks <= max_j 2*2^(p-j)*s_j = O(m p^2),
    retained moving helper wires = m/2-1,
    total elementary gates = O(m^2 p^2).

Fixed field-basis changes and Fourier matrices are specified Clifford maps;
their elementary gates are charged in the last bound. The private masks return
after each multiplication batch, so they can be reused between stages. They
are distinct from the retained moving T_j registers. The actual reversed
native word has the same asymptotic and T-band bounds. These bounds follow by
composition with the retained all-size multiplication emitter, not from the
small fixture alone.

The coherent helper query P_p(x) is Clifford: every beta is binary-linear in
address x, so its phase is a fixed bilinear form between address and T_j bits,
implemented by CZ gates. It costs no T gates and O(mp) elementary gates. It
must be emitted along with the eight-band field query, not silently omitted.

## 4. Complete readout and its limitation

Let S_x be the actual native source Ehat_p (Q_p(x) P_p(x)) Ehat_p^dagger.
Equation (3) gives

    S_x (Z_Y0 tensor I_helpers) S_x^dagger
       = Z_Yx tensor I_helpers.

The moving helper query cancels in this physical-probe echo. The complete
indicator wrapper therefore returns every retained helper and every private
mask on arbitrary inputs, including entanglement with addresses and references.
It is incorrect to say that Ehat_p or S_x itself returns the retained T_j's;
S_x acts on them by the explicitly present P_p(x). Only the actual source/probe/
adjoint composition gives identity helper action. All source and actual
inverse calls remain charged. With the retained eight-band coherent field query,
the direct complete indicator uses at most64(p-1)+16 bands: four decoder or
actual-adjoint appearances and two field-query appearances. This remains O(p).

In the helper X basis, equation (3) is the Clifford cocycle

    E_(r+Lambda x) Q_x = W_x E_r.

A generic formula using a controlled W_t can hide the original barrel-shift
cost. The construction here instead gives an explicit alternate compiler for
that cocycle, with16 paid bands at **each** recursive level. Sharing helpers
removes helper-swapping overhead, but does not collapse the dependencies of
successive principal-core stages. Its depth is O(p), not o(log p), and this is
not a new general indicator or frame-depth upper theorem. The open next
operation is a joint synthesis across these levels that compresses the
cocycle without reinstating controlled wire translations or unpriced helpers.

## 5. Bounded exact evidence

`shared_helper_fixture.py` constructs the exact Clifford word in each shared
helper X sector for p=2,3. It verifies literal actual inverses and all columns
of E_(r+Lambda x) Q_x = W_x E_r for every helper character and every address.
This tests the shared-helper induction, including the full zero sectors. The
physical nonlinear primitives remain owned by the retained multiplication-XOR
proof; the sector fixture is not itself a physical Clifford-only implementation
of the coherently controlled helper family.

The field-domain obstruction for an identity-helper source cannot be imported
unchanged here: Q_p tensor P_p(x) can negate a helper Pauli whose binary label
is fixed. The old source premise that every fixed label has literal positive
conjugation is therefore absent. No conclusion about that extended source
class is inferred from the existing domain theorem.
