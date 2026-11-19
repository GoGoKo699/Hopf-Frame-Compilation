# One-clean precision bounds at each native interface

Let $`r_{ij}`$ be the symplectic pairing rank of consecutive transported native axis spans
$`L_i,L_j`$, and let $`s_{23}=\dim(L_2\cap L_3)`$. A four-layer, one-clean implementation of the
complete rotation witness $`I\otimes R_y(\tfrac12\arccos(1/3))`$ with initialized error
$`2^{-L}`$ satisfies
```math
\min(r_{12},r_{23},r_{34},s_{23})=\Omega(L).
```
The bound holds for every native four-layer presentation and arbitrary dirty width. It concerns
increasing precision. The unrestricted one-clean depth lower bound remains open.

Take a native word $`V=V_2V_1`$ with two T bands in each half, one external clean bit and
arbitrarily many dirty bits. Let $`P`$ be a Hermitian data Pauli, and put
```math
 A=V_2^\dagger P V_2=E_A K_A E_A^\dagger,\qquad
 B=V_1 P V_1^\dagger=E_B K_B E_B^\dagger .
```
Each $`E_j`$ is one Clifford-conjugated physical T band, with sign choices arbitrary and
independent unsigned binary axis generators spanning $`L_j`$. Each $`K_j`$ is a Hermitian
Clifford. If the terminal clean Z is preserved, let
```math
 S=V_2^\dagger Z_c V_2=V_1 Z_c V_1^\dagger
  =E_A J_A E_A^\dagger=E_B J_B E_B^\dagger,\qquad Q=(I+S)/2 .
```
The clean Clifford cores $`J_A,J_B`$ are Hermitian, and $`A,B`$ both commute with $`S`$.
[R75 §6](FOUR_LAYER_OPERATOR_RIGIDITY.md#6-one-clean-qubit-terminal-sectors-only) establishes
terminal preservation automatically once the initialized complete-return accuracy is below its
absolute threshold. Let $`F_j`$ be the binary symplectic action of $`J_j`$, an involution, and
define
```math
 u_j=\mathrm{rank}_{\mathbb F_2}\big(\langle x,F_jy\rangle\big)_{x,y\in L_j},\qquad
 r=\dim(L_A\cap L_B).
```
## Bounding the inner rank by the common span

Write $`I=L_A\cap L_B`$. For $`h\in L_A^\perp`$, the conjugate of its Pauli by the common $`S`$
is a one-band image of a Pauli with binary label $`F_Ah`$. This is a Clifford. Its symplectic
matrix is $`I_{\rm symp}+N_A(h)`$, where
```math
 N_A(h)b=\sum_{g\ {\rm an\ A\ axis}}\langle F_Ag,h\rangle\langle g,b\rangle\,g.
```
The formula uses transvections of commuting axes; signs do not affect it. Independence of the
native band axes ensures the image is exactly the span of the active axes in this sum. Define
$`U_A`$ as the span of all axes that can activate for some such $`h`$, namely those $`g`$ with
$`F_Ag\notin L_A`$. Define $`N_B(h'),U_B`$ symmetrically for $`h'\in L_B^\perp`$.

For $`h\in(L_A+L_B)^\perp`$ the images in both presentations agree, so the image of $`N_A(h)`$
belongs to $`I`$. If an individual axis $`g\notin I`$ were active for such $`h`$, its inclusion
in this image would be a contradiction. Therefore every such generator satisfies
$`F_Ag\in L_A+L_B`$.

For $`h\in L_A^\perp,h'\in L_B^\perp`$ the two Clifford images $`S P_h S,S P_{h'} S`$ commute
projectively. Their symplectic matrices commute, so $`N_A(h)N_B(h')=N_B(h')N_A(h)`$ has image in
$`I`$. On spanning over $`h'`$, this shows $`N_A(h)b\in I`$ for every $`b\in U_B`$. In A axis
coordinates these vectors are the entrywise products $`c(h)\odot d(b)`$, where
$`c_g(h)=\langle F_Ag,h\rangle`$ and $`d_g(b)=\langle g,b\rangle`$. If $`g\in I`$ then
$`d_g(b)=0`$. Let $`r_0`$ count the A axis generators lying in $`I`$. The intersection of $`I`$
with the coordinate subspace omitting these $`r_0`$ coordinates has dimension exactly $`r-r_0`$:
subtract the individual basis generators in $`I`$ from any vector to remove its entries there.

Restrict the $`c,d`$ coordinates to the globally active A axis generators. Extending scalars
from $`\mathbb F_2`$ to its algebraic closure does not change the binary dimensions or the
product inclusions. The linear space spanned by the activation words has full support on these
coordinates, so contains a word nonzero at *every* one (a finite union of proper linear
hyperplanes cannot cover a vector space over an infinite field). Multiplication by this word is
injective. Every product still lies in the projected space of dimension at most $`r-r_0`$. Hence
$`\mathrm{rank}\langle U_A,U_B\rangle\le r-r_0`$.

The part $`K=L_B\cap F_BL_B`$ is $`F_B`$-invariant. For $`b\in K`$, conjugation by $`S`$ gives a
Pauli with label $`F_Bb`$ because $`E_B`$ commutes with both $`b`$ and $`F_Bb`$. Consequently
each $`S P_h S`$ for $`h\in L_A^\perp`$ fixes every Pauli in $`K`$ projectively: its commutation
sign with $`b`$ is the original Pauli commutation sign of $`P_h`$ with $`S P_b S`$. Thus
$`K\subseteq\ker N_A(h)=(\mathrm{im}N_A(h))^\perp`$, giving $`U_A\perp K`$. Every inactive B
generator lies in $`K`$, so the rank of $`\langle U_A,L_B\rangle`$ is at most $`r-r_0`$.

In the matrix for the symmetric form $`\langle g_i,F_Ag_j\rangle`$, rows for inactive A
generators are zero by symmetry. On all columns indexed by $`g_j\notin I`$ the vectors
$`F_Ag_j`$ lie in $`L_A+L_B`$, as proved above. Pairing them against the active A rows has total
rank at most $`r-r_0`$. The remaining $`r_0`$ columns add at most $`r_0`$. This proves
$`u_A\le r`$; interchange A and B to obtain $`u_B\le r`$.


## Choosing the coordinates of a clean Clifford core

Put $`u=u_B`$. There is a maximal commuting Pauli space $`M`$ containing $`L_B`$ with
$`\mathrm{codim}_M(M\cap F_BM)=u`$. Indeed, on $`L_B`$ the pairing
$`\beta(x,y)=\langle x,F_By\rangle`$ is symmetric since $`F_B^2=I`$. Choose a complement $`U`$
to its radical $`R`$, so the form on $`U`$ is nonsingular. The space $`T=U\oplus F_BU`$ is
nondegenerate symplectic and invariant. In $`T^\perp`$ the space $`R+F_BR`$ is isotropic and
invariant. Extend it to an invariant Lagrangian $`M'`$: on its symplectic quotient the
involution has isotropic image of $`F_B-I`$, and any Lagrangian containing that image is
invariant. Lift this Lagrangian to $`T^\perp`$. Then $`M=U+M'`$ works. Conversely the rank of
the pairing on any extension $`M`$ is at least its rank on $`L_B`$, proving optimality.

Using $`M`$ for the computational $`Z`$ coordinates, the Clifford normal form for $`J_B`$ has
monomial Clifford ends and precisely $`u_B`$ internal Hadamards. In any compatible chart, write
$`t_B`$ for that number, equivalently the logarithm to base two of the support size of one
column. Thus $`u_B=\min t_B`$ and likewise on side A. Coordinates in a chart may be chosen with
the physical band generators among the coordinate axes.

## A fixed number of boundary terms on the occupied code

Write $`Q=(I+S)/2`$. In an A-adapted chart whose first coordinates are the actual native axes,
choose a Pauli $`H`$ flipping at most three computational bits. Put $`H_+=(H+SHS)/2`$ and
$`D=AH_+A-BH_+B`$. Then $`H_+`$ commutes with $`Q`$, is a contraction, and
```math
\|DQ\|_{2,\mathrm{tr}}\le 2\|(A-B)Q\|_{2,\mathrm{tr}}.
```
The A terms involve $`AHA`$ and $`(AS)H(AS)`$, where $`A`$ and $`AS`$ are $`E_A`$ conjugates of
Hermitian Cliffords. If $`H`$ activates $`h\le3`$ A axes, $`E_A^\dagger H E_A`$ has $`2^h`$
signed Pauli terms. Each produces an ordinary Clifford after the remaining conjugations. The
corresponding B terms are $`E_B C' E_B^\dagger`$ with $`C'`$ Clifford, for every $`H`$.

Expanding the projector, $`(DQ)E_B`$ thus has four forms, $`CE_B`$, $`CE_BJ_B`$, $`E_B C'`$ and
$`E_B C'J_B`$, with Clifford $`C,C'`$. In a B chart write $`J_B=M_2H_T M_1`$ with monomial
Clifford ends and $`|T|=t_B`$. Right-multiply by the common $`M_1^\dagger`$. Three forms now
have flat affine Choi amplitudes with weighted-cubic phase. The remaining $`CE_BM_2H_T`$ is the
same kind of vector followed by a partial Hadamard on $`t_B`$ input coordinates. The
$`4\cdot2^h`$ A-side coefficients have magnitude $`(\sqrt2)^{-(4+h)}`$ and the four B-side
coefficients have magnitude $`(\sqrt2)^{-4}`$. With $`\omega=e^{i\pi/4}`$, repeated use of
$`1=(\omega+\omega^{-1})/\sqrt2`$ gives exactly $`K=64`$ terms of denominator exponent $`b=7`$.
All phases are physical eighth-root Clifford phases, including products in which the Clifford
need not be Hermitian.

The
[two-class finite-boundary estimate](FOUR_LAYER_LINEAR_PRECISION.md#separating-different-clusters-and-transferring-the-witness)
now applies with $`m=t_B`$. If $`DQ\ne0`$, its squared normalized Hilbert–Schmidt norm is at
least $`2^{-T(t_B)}`$, where
```math
 T(t)=(32+4t)(126^{64}-1)/125 .
```
Hence $`\|(A-B)Q\|_{2,\mathrm{tr}}\lt 2^{-T(t_B)/2-1}`$ forces $`DQ=0`$ for every probe. The
mirror construction gives the same statement for $`t_A`$.

To recover equality on the occupied code, put $`v_x=Q|x\rangle`$ in the A chart. Fourier
resolution using all $`Z`$ masks gives $`|v_x\rangle\langle v_y|`$ for each pair separated by at
most three flips. Products generate all matrix units on a connected component of the graph of
nonzero $`v_x`$, since intermediate squared norms are positive. This graph is connected:

* If $`J_A`$ is nonmonomial, each column has at least two nonzero entries, so no $`v_x`$ vanishes.
* If its affine permutation is nonidentity, every vertex outside the affine fixed set is
  present. Its complement is connected with steps of length at most two. In a cube-edge path,
  translate each segment lying inside the fixed flat by one coordinate with nonzero defining
  syndrome, joining the outside endpoints by steps of length at most two. Each included fixed
  point has an outside neighbor. An empty fixed set already gives the full cube.
* If its permutation is identity, Hermiticity makes the positive set the support of a binary
  polynomial of degree at most two. For two allowed vertices at distance greater than three,
  vary three differing coordinates toward the destination. Even parity on their three-cube
  supplies another allowed vertex, strictly closer. Distances at most three are already edges.

The nonzero $`v_x`$ span the code. The compressed probes therefore generate its full matrix
algebra. On this space $`DQ=0`$ means the conjugations by $`A,B`$ coincide on all generators, so
$`BA`$ is scalar there. Hermiticity gives $`A=\pm B`$ on the code, and the closeness assumption
selects equality.

## From core ranks to physical crossings

For the rotation witness, sufficiently small initialized error $`\eta`$ forces exact terminal
clean return by R75. Each target-preserved Pauli then satisfies
$`\|(A-B)Q\|_{2,\mathrm{tr}}\le2\eta`$. If all these comparisons are equal, the spectator and
dirty Paulis and logical $`Y`$ force one occupied $`Y`$ rotation up to phase. Its transfer
entries lie in $`\mathbb Z[1/\sqrt2]`$. The
[R76 witness calculation](ONE_CLEAN_MIDPOINT_RIGIDITY.md#4-initialized-error-and-the-frame-witness)
restricts the angle to $`(\pi/8)\mathbb Z`$, a fixed positive distance from the target. Thus
$`T(t)\ge 2L-4`$ for both optimized chart ranks, and $`u_A,u_B=\Omega(L)`$.

For each inner band, its clean Clifford core is the outer physical band's image of a Pauli up to
frame transport. Let $`O`$ be the outer axes anticommuting with that Pauli. Its symplectic
action is $`Fh=h+\sum_{o\in O}\langle o,h\rangle o`$. For an inner-axis basis $`g_i`$, the form
$`\langle g_i,Fg_j\rangle`$ therefore has matrix $`TT^\mathsf{T}`$ with
$`T_{io}=\langle g_i,o\rangle`$. Its rank is at most that of the outer-to-inner pairing and at
most $`|O|`$. The latter is precisely the
[inner clean spread](ONE_CLEAN_INNER_SPREAD.md#1-inner-clean-spread) of that two-layer encoder,
so both physical ends have $`\Omega(L)`$ clean-active sites. The common-span inequality above
gives $`s_{23}\ge \max(u_A,u_B)`$ and these two outer pairings give $`r_{12},r_{34}=\Omega(L)`$.
Independently, [R74](FOUR_LAYER_LINEAR_PRECISION.md) gives $`r_{23}=\Omega(L)`$. A common
direction can be a product of native axes.

For a single-address bit response with **exact terminal clean return**, average over all
initialized computational columns including dirty inputs. If the tested data bit can change only
at one $`k`$-bit address and has positive response there, then $`0\lt e\le2^{-k}`$. The occupied
sector has half the physical dimension, giving $`\|(A-B)Q\|_{2,\mathrm{tr}}^2=2e`$. The probe
gap implies $`T(t)\ge k-3`$. Consequently $`u_A,u_B,r_{12},r_{34},s_{23}`$ are each
$`\Omega(k)`$. The additional middle-crossing bound above uses the rotation witness.
