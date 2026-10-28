# Independent analytic audit: augmented uniform field recoupler

Verdict: **the proposed complete cell and shared-helper recursion are valid**.
This is a scratch continuation outside the frozen PR269 tree. It changes the
query representation on both sides by an explicitly paid helper Pauli query;
it is not a decoder for Q alone with identity helper action.

## 1. Registers, Fourier convention and the literal cell

Let F=F_(2^ell), K=F[theta], theta^2+theta+delta=0. Write

    U=d0+d1 theta,  V=v0+v1 theta,
    q(U)=d0 d1,
    B(U,a)=d0 a1+d1 a0.

The arbitrary helper T is one F register. Let Z_T^lambda have phase
chi_F(lambda T), where chi_F(z)=(-1)^Tr_(F/F2)(z). The full K Fourier used
inside the cell is

    H_theta[a,v] = |K|^(-1/2) chi_K(theta a v).

It is a literal self-adjoint involution: theta is nonzero and its symmetric
trace pairing is nondegenerate. In any fixed binary basis it is a fixed
Clifford word. No Fourier state is prepared.

Let A_f add f(U,V) to T without changing U,V, and let A_g add g(U,a) to T:

    f(U,V)=Tr_(K/F)(UV)+q(U),
    g(U,a)=B(U,a).

The chronological cell is A_f, H_theta, A_g, so its matrix is
C=A_g H_theta A_f. Both additions are literal involutions on all inputs.
With b a binary bit and lambda in F, define the source phase

    D_(lambda,b)(U,V)=chi_K((lambda+b theta)UV).

Then the complete operator identity is

    C [ D_(lambda,b) tensor Z_T^lambda ] C^dagger
      = chi_F(lambda[B(U,a)+q(U)]) X_a^(bU) tensor Z_T^lambda.    (1)

All U,V,T inputs are arbitrary, including U=0, and the scalar in (1) is +1.

### Direct conjugation proof

Conjugating the helper Z by the first addition multiplies it by
chi_F(lambda[Tr(UV)+q(U)]). The lambda*Tr(UV) part of the source phase cancels,
leaving chi_K(b theta UV), the constant-in-V factor chi_F(lambda q(U)), and
the same helper Z. H_theta sends chi_K(b theta UV) to translation X_a^(bU).
The second addition commutes with this translation because

    B(U,a+bU)+B(U,a) = b B(U,U)=0.

It dresses the helper Z by chi_F(lambda B(U,a)). This proves (1). The final
phase also commutes with the translation by the same alternating identity.
No norm division, zero test, nonzero promise, or discarded fiber is present.

### Independent kernel check

In the helper-X basis, let eta in F label the character chi_F(eta T). Since
all helper operations are additions, the cell is block diagonal in eta:

    C_eta[a,v] = |K|^(-1/2) chi_K(theta a v)
                * chi_F(eta[B(U,a)+Tr(Uv)+q(U)]).

Z_T^lambda shifts eta to eta+lambda. Equation (1) is exactly

    C_(eta+lambda) D_(lambda,b)
      = chi_F(lambda[B(U,a)+q(U)]) X_a^(bU) C_eta.               (2)

Substituting the displayed kernel proves (2) entry by entry: the two
lambda*Tr(Uv) terms cancel, while a->a+bU adds b*Tr_K(theta Uv) and leaves
B(U,a) unchanged. At U=0 every U-dependent phase and translation vanishes;
C_eta=H_theta is independent of eta. Thus the entire zero sector is checked.

## 2. The exterior recoupler still has the required action

Before C, apply the usual H_K on V to turn Q_p into D_(lambda,b). After C,
apply the same fixed coordinate sum c=a+U and the two H_F walls on the w
coordinates as in the existing recursive proof. Its intermediate lower-query
phase is exactly chi_F(lambda[B(U,a)+q(U)]), while the upper-bit half swap is
exactly a->a+U. Therefore the augmented recoupler obeys

    Rhat_j [ Q_j(x',b) tensor Z_(T_j)^beta_(j-1)(x') ] Rhat_j^dagger
      = [ Q_(j-1)(x') tensor Q_(j-1)(x') ] S^b
        tensor Z_(T_j)^beta_(j-1)(x').                          (3)

The extra helper T_j has ell_j=2^(j-2) qubits. All older-level helper queries
commute through this recoupler because it acts only on its core and T_j.
Every displayed inverse means the actual reversed word.

## 3. One shared helper per level is sufficient

For p address bits, allocate one T_j for each j=2,...,p, shared by every
recursion cell at level j. Define the single helper Pauli query

    P_p(x) = product_(j=2)^p Z_(T_j)^beta_(j-1)(x0,...,x_(j-2)). (4)

There is one factor per level, **not one factor per cell**. Equation (4) is a
Clifford query: its exponent is binary bilinear in address bits and helper
bits, implemented with CZs and fixed parity CNOTs if desired. It uses no clean
input and no T gate.

Every complete recoupler commutes with every individual X on its helper T_j:
an addition T_j->T_j+f(core) commutes with all helper translations, and the
other walls act only on the core. By induction every E_j commutes with every
individual X on all persistent helpers that it uses.

Two identical child decoders on disjoint cores, using the same older helpers,
therefore commute. In the simultaneous helper-X basis they are block diagonal,
and at each fixed helper label their matrices act on disjoint cores. The
core-only half swap interchanges the two child operators; their commuting
product is consequently invariant under that swap.

There is no erroneous squaring or cancellation of the shared Pauli query.
Write the child identities as

    E_L (Q_L P_child) E_L^dagger = W_L P_child,
    E_R (Q_R P_child) E_R^dagger = W_R P_child.

Then conjugating Q_L Q_R P_child first by E_R and then by E_L gives

    Q_L W_R P_child  ->  W_L W_R P_child.

At the second step W_R commutes with E_L because it acts only on the other
core. The one shared P_child persists. Together with swap invariance and (3),
this proves inductively

    E_p [ Q_p(x) tensor P_p(x) ] E_p^dagger
       = W_p(x) tensor P_p(x)                                  (5)

on the full joint space. The base E_1 is the retained Clifford decoder with
no helper. Equation (5) includes coherent addresses and arbitrary helper and
external-reference inputs.

Only complete ideal blocks have been commuted in this proof. It does not
license rearranging overlapping native gates from different blocks individually.

## 4. Paid simultaneous implementation at one level

The first addition has a useful two-product form:

    Tr_(K/F)(UV)+q(U)
       = d0(v1+d1) + d1(v0+v1),

and the second is

    B(U,a)=d0 a1+d1 a0.

Both are sums of two bilinear F products. The existing interpolation-based
multiplication-XOR component supplies each product as at most

    s(ell)=(2ell-1) e(ell)^2,
    e(ell)=max(1,ceil(log2(2ell-1)))

products of three binary linear forms after a Hadamard wall on the target
helper. The returned eight-phase Gray compiler applies to all such terms at
once, including overlapping principal forms. A private phase mask is assigned
to every term; these masks are additional arbitrary helpers and return exactly.

At a given level, all first additions from every cell can be combined into
one XOR into the shared T_j, followed by all disjoint core Fourier walls,
then all second additions combined into another XOR. Additions from different
cells commute, and each core Fourier commutes with every other cell's addition.
Thus this is the same complete product of recouplers, with two8band batches:
**16 T bands per level**. The target-Hadamard walls on shared T_j are emitted
once per aggregate addition, not once per cell.

Let c_j=2^(p-j) be the number of cells at level j. Safe exact bounds are

    T(E_p) <= 32 sum_(j=2)^p c_j s(ell_j),
    D_T(E_p) <= 16(p-1),
    persistent helper qubits = sum_(j=2)^p ell_j = 2^(p-1)-1,
    private phase masks <= max_(j=2)^p 2 c_j s(ell_j).

The private mask pool is reused only after each complete XOR returns it.
These bounds give T=O(m p^3), dirty width O(m p^2), and total elementary gates
O(m^2 p^2), with m=2^p. Fixed field-coordinate and Fourier walls are included
in the latter polynomial gate allowance. The actual inverse has the same
asymptotic budgets and reverses every target wall and mask word.

These are conservative native construction bounds assembled from the already
paid cubic compiler. They do not assume initialized products or clean helper
copies. Persistent T_j registers may move during E_p; the private phase masks
return at every completed XOR stage.

## 5. Complete indicator meaning and remaining gap

The coherent query Q_p tensor P_p uses the existing eight-band Q_p plus a
Clifford P_p. With S_x=E_p(Q_p tensor P_p)E_p^dagger, equation (5) gives

    S_x (Z_(Y0) tensor I_helpers) S_x^dagger
       = Z_(Yx) tensor I_helpers.

The usual actual-source/adjoint probe echo therefore returns all persistent
and private helpers on arbitrary inputs and implements the complete dirty
indicator. It pays four E_p/adjoint appearances and two field queries. A safe
T-band bound from this sequential-level ledger is64(p-1)+16, with polynomial
T-count, dirty width and total-gate overheads as above.

This is a uniform paid augmented intertwiner and a valid source interface.
It is **not** the identity-helper intertwiner Q_p->W_p, nor the earlier
W_m tensor residual-field-shear character construction. The helper Pauli query
appears on both sides and is essential to the first cancellation.

Depth16(p-1) for E_p, or64(p-1)+16 for the indicator wrapper, does not improve
the general frontier or meet the desired sublogarithmic source-depth target.
A joint native compilation across levels remains unpaid. In particular the
shared-helper argument proves exact composition and width reduction; it does
not remove the sequential layer dependence by itself.

The parent is independently running the p2/p3 shared-helper-X-sector fixture,
checking E_(eta+Lambda x) Q_x = W_x E_eta and actual inverses. This audit does
not duplicate or claim that executable receipt; its proof above covers the
complete arbitrary-size operator statement under the specified tower/query
and native compiler assumptions.
