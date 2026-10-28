# A dependent commuting-Pauli bank can approach an independent bank

This is a scratch follow-up, separate from frozen PR270. It proves a full-
space operator construction and a sharpness example for an argument based
only on one lost Pauli relation. It supplies no low-degree query family,
no native compilation ledger and no counterexample to the sector-degree
or weak-probe theorems.

## 1. The two banks

Let m>=2. On m physical qubits, set

    A_i=Z_i                         (0<=i<=m-2),
    A_(m-1)=product_(i=0)^(m-2) Z_i.

The last physical qubit, called the copy coordinate, is absent from all A_i.
These are commuting Hermitian Paulis with one relation product_i A_i=I.
Their binary label span has dimension m-1.

There is nevertheless a unitary U from a target m-qubit bank to this physical
space for which the independent commuting observables B_i=U Z_i U^dagger obey

    ||A_i-B_i|| = 2/sqrt(m)          for every i.             (1)

A variant retains an exact Pauli seed:

    B_0=A_0,
    ||A_i-B_i||=2/sqrt(m-1)          for 1<=i<=m-1.           (2)

Thus the maximum individual-observable error tends to zero with m, even
though the surrogate Pauli bank loses one label dimension. The independent
bank itself loses no dimension or input column.

## 2. An explicit unitary on every input

On the target bank define the Jordan--Wigner reflections

    gamma_j = Z_0 ... Z_(j-1) X_j,       0<=j<m.

They are Hermitian involutions, anticommute pairwise, and anticommute with
the parity P=product_i Z_i. Choose a nonempty index set J, put k=|J|, and set

    S=sum_(j in J) gamma_j,             H=S/sqrt(k).

Pairwise anticommutation gives S^2=kI, so H is a Hermitian unitary and swaps
the even and odd parity sectors. Let P_+,P_- be their projectors.

An even bitstring b is uniquely specified by its first m-1 bits a, since its
last bit is the parity of a. For c=0,1 let J_c embed the entire even sector
into the physical subspace |a,c>. Define

    U = J_0 P_+ + J_1 H P_- .                                  (3)

The two summands have orthogonal output spaces. Each is an isometry on its
input parity sector, and those sectors exhaust the target bank. Hence
U^dagger U=UU^dagger=I. Equation (3) does not initialize or discard a copy
qubit: it is a square unitary on all 2^m columns, including coherent
superpositions of the parity sectors and arbitrary references.

The observable intertwining on the even sector is

    A_i J_c = J_c Z_i P_+                         for every i,c.

For i=m-1 this uses even parity; for the other i it is immediate from the
encoding. Consequently

    U^dagger A_i U = P_+ Z_i + P_- H Z_i H.                     (4)

The B_i commute and are independent because they are a simultaneous unitary
conjugate of the complete physical Z bank. Their product is

    product_i B_i = U P U^dagger = Z_copy,                     (5)

whereas product_i A_i=I. Thus the missing surrogate relation is visible in
the copy qubit's action; it has not been silently dropped.

## 3. Exact operator norms

Only gamma_i anticommutes with Z_i. Every gamma_j for j!=i commutes with Z_i,
even when its Jordan--Wigner string contains Z_i. Therefore

    [Z_i,H] = 2 Z_i gamma_i / sqrt(k)             if i in J,
    [Z_i,H] = 0                                 if i not in J.

For i in J,

    H Z_i H-Z_i = H[Z_i,H]
                 = 2 H Z_i gamma_i / sqrt(k).                 (6)

The factor H Z_i gamma_i is a unitary preserving each parity sector, because
both H and gamma_i flip parity. Thus the restriction of (6) to the odd sector
has norm exactly 2/sqrt(k), not merely an upper estimate. Equations (4) and
(6) prove

    ||A_i-B_i|| = 2/sqrt(k)   (i in J),
    ||A_i-B_i|| = 0           (i not in J).                    (7)

Taking J={0,...,m-1} proves (1); taking J={1,...,m-1} proves the exact-seed
variant (2). The relation-word discrepancy in (5) has norm two, which is
consistent with errors accumulated over all m individual observables.

## 4. Optional weak-probe embedding and its explicit scope failure

For the exact-seed variant let E=U^dagger and R=E^dagger Z_0 E=B_0=A_0.
Choose individually Clifford query operators

    Q_i=SWAP(0,i)                      (0<=i<=m-2),
    Q_(m-1)=product_(j=1)^(m-2) CX(j,0).

Here SWAP(0,0)=I. These literal Clifford words satisfy

    Q_i R Q_i^dagger=A_i,
    ||E Q_i R Q_i^dagger E^dagger-Z_i|| <= 2/sqrt(m-1).         (8)

They are not asserted to form an additive representation. Their coherent
address-controlled implementation is not paid by this construction.

More strongly, if m=2^p with p>=2, their binary symplectic query matrix F_x
has Boolean degree exactly p. For any j in {1,...,m-2}, sum the Z_j label
images over every address. Among the first m-1 queries, Q_j sends z_j to z_0
and the other m-2 queries fix z_j. Since m-2 is even, their XOR is z_0.
The final CNOT product also fixes z_j, so

    sum_(x in F2^p) F_x z_j = z_0+z_j != 0.                   (9)

This is the full p-fold Boolean derivative. It proves degree at least p,
and a p-variable Boolean function has degree at most p. Hence the degree
is exactly p, explicitly outside a fixed-low-degree query hypothesis.
Equation (8) is therefore no counterexample to a degree-sensitive weak-probe
or sector theorem. It isolates the information missing from a Pauli-label-
deficit argument alone.

## 5. Exact bounded fixtures

`pauli_surrogate_sharpness.py` checks m=4 and m=8 for both index sets J.
No SVD or floating norm estimate is used. In parity-ordered bases U has blocks
I and S_eo/sqrt(k), where S_eo is an integer signed matrix. The observable
B_i has physical copy blocks

    Z_i|even,      C_i/k,
    C_i=S_eo (Z_i|odd) S_eo^T.

The fixture verifies exact integer identities

    S_eo S_eo^T=S_eo^T S_eo=kI,
    C_i^2=k^2 I,                  [C_i,C_j]=0,
    C_i S_eo=k S_eo (Z_i|odd),
    (k Z_i|even-C_i)^2=4kI       for i in J,
    k Z_i|even-C_i=0             for i not in J.

These certify every unitary column, every joint-eigenvalue column, and the
exact norm in (7). The fixture also checks product_i C_i=-k^m I, all
Jordan--Wigner anticommutators on every basis column, the explicit Clifford
seed queries, and the nonzero top query derivative (9). It preserves the
full copy-coordinate action throughout.
