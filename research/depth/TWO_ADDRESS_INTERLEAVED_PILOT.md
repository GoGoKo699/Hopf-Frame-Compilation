# A complete two-address source with a later queried stage

[Minimal tree source](ORDERED_MATCHING_PROBE_SOURCE.md#minimal-tree-source) ·
[Fixed-readout boundary](FIXED_READOUT_SYNDROME_BOUND.md) ·
[Balanced Clifford compiler](CLIFFORD_INVOLUTION_SOURCE.md)

A finite source can satisfy the entire interleaved observable family,
including its noncommuting pulled-back targets. The construction below
implements the complete coefficient lookup at p=2, with arbitrary carrier,
output and dirty inputs. Its first nonlinear operations are controlled
Hadamards; a complete address-dependent tree XOR query occurs later.
Literal phases, the fixed decoder and actual native inverses are paid.

This is a complete benchmark for the changed architecture, not an all-size
shallow source or a new frame-depth bound. In particular, it does not
improve R96 or replace the retained logarithmic-depth indicator.

## 1. The complete target and its pulled-back family

Write the address as $`x=(a,b)`$ and order coefficient subsets as
$`\varnothing,\{1\},\{2\},\{1,2\}`$. Use two arbitrary four-bit
banks c,t and set

```math
v_x=(1,a,b,ab)^{\mathsf T},\qquad
M_x=\begin{pmatrix}
1&a&b&0\\0&1&0&0\\0&0&1&a\\0&0&0&1
\end{pmatrix},\qquad
Q_x=Z_t(v_x)=(-1)^{t_0+a t_1+b t_2+ab t_3}.             \qquad(1)
```

All binary vectors and matrices in this note are over the binary field.
The complete forward query and fixed Clifford wall are

```math
J_x|c,t\rangle=|c,t+M_xc\rangle,\qquad
C=\left[\prod_{j=0}^3\mathrm{CX}(c_j;t_j)\right]
       H_c^{\otimes4},\qquad B_x=CJ_x.                 \qquad(2)
```

Thus $`C^\dagger Z_t(v)C=X_c(v)Z_t(v)`$. Since
$`M_x^{\mathsf T}v_x=e_0`$, the exact pulled-back target is

```math
P_x=B_x^\dagger Q_xB_x
 =X_c(v_x)X_t(M_xv_x)Z_c(e_0)Z_t(v_x).                 \qquad(3)
```

The order in (3) retains the literal phase. These operators are Hermitian
Paulis, but their family does not commute: $`P_{00}`$ anticommutes with
$`P_{10}`$ and $`P_{01}`$, while $`P_{01}`$ commutes with $`P_{11}`$.
Matching only four balanced
spectra would therefore be insufficient.

## 2. Two controlled-Hadamard stages give the right complete relations

On the physical $`c_0,c_1`$ wires define four canonical observables

```math
F_{00}=X_{c_0},\quad F_{10}=Z_{c_0},\quad
F_{01}=Z_{c_0}X_{c_1},\quad F_{11}=X_{c_0}Z_{c_1}.
                                                               \qquad(4)
```

They are independent and have exactly the pairwise commutation relations
of (3). A signed decoder is supplied in the next section, so no existence
claim about an unpriced change of basis is needed.

Let $`F=H_{c_1}\mathrm{CX}(c_1;c_0)`$ and

```math
K_1=H_{c_0}H_{c_1},\qquad
K_2=F H_{c_0}F^\dagger
    =\frac{X_{c_0}+Z_{c_0}X_{c_1}}{\sqrt2},\qquad
W_x=K_1^aK_2^b.                                       \qquad(5)
```

Both $`K_i`$ are Hermitian Clifford involutions. $`K_2`$ sends $`X_{c_0}`$ to
$`Z_{c_0}X_{c_1}`$; $`K_1`$ sends those two possibilities to $`Z_{c_0}`$ and
$`X_{c_0}Z_{c_1}`$. Direct conjugation therefore gives all four literal
identities

```math
W_x X_{c_0}W_x^\dagger=F_x.                            \qquad(6)
```

The address-controlled $`K_2`$ executes $`F^\dagger`$, a controlled Hadamard
and F, in that chronological order, using the actual fixed-word adjoint.
The address-controlled $`K_1`$ consists of two controlled Hadamards sharing a
but using distinct targets. These are controlled Clifford fibers, not
literal XOR queries. Their nonalternating Pauli action is a material
change from merely adding another fixed decoder to an XOR-query orbit.

## 3. An explicit signed Clifford decoder

Number the carrier wires $`(c_0,c_1,c_2,c_3,t_0,t_1,t_2,t_3)`$
from zero through seven. The following table is the **chronological**
word E. Within each row, execute gates from left to right; continue from
the first row to the last. Every power of S is literal, and no phase is
removed.

| Block | Chronological Clifford gates |
| --- | --- |
| 1 | $`S_0^3,\ S_4,\ \mathrm{CX}(0;4)`$ |
| 2 | $`H_1,\ S_4,\ H_4,\ S_5,\ H_5,\ \mathrm{CX}(1;0),\ \mathrm{CX}(4;0),\ \mathrm{CX}(5;0)`$ |
| 3 | $`H_0,\ S_0,\ H_0,\ H_1,\ H_5,\ S_6,\ \mathrm{CX}(1;2),\ \mathrm{CX}(1;5),\ \mathrm{CX}(1;6)`$ |
| 4 | $`H_3,\ S_6,\ H_6,\ S_7,\ H_7,\ \mathrm{CX}(3;1),\ \mathrm{CX}(6;1),\ \mathrm{CX}(7;1)`$ |
| 5 | $`H_1,\ S_1^3,\ H_1`$ |

Signed Pauli propagation through this word gives

```math
E(P_{00},P_{10},P_{10}P_{01},P_{00}P_{11})E^\dagger
 =(X_{c_0},Z_{c_0},X_{c_1},Z_{c_1}).                    \qquad(7)
```

All four signs in (7) are plus. Set $`D=E^\dagger`$. Equations (4) and
(7) give

```math
D F_xD^\dagger=P_x\quad\hbox{for every }x.               \qquad(8)
```

This also matches every product and commutator of the four generators,
including products whose literal phase is imaginary. The table has 35
native Clifford gates. Reverse the actual word, expand S-dagger as three
S gates and reduce adjacent S powers modulo four: D has 45 gates. Both
are included in the common decoder bound 45. The fixture emits this
exact word and checks its signed action on all 16 generator products.

## 4. Complete source, arbitrary-work return and native ledger

The input probe is explicitly the noninitialized Pauli $`Z_{c_0}`$.
Define

```math
V_x=B_xD W_xH_{c_0}.
```

Equations (3), (6) and (8) give the complete identity

```math
V_x Z_{c_0}V_x^\dagger=Q_x\otimes I_c\otimes I_{\rm masks}.
                                                               \qquad(9)
```

No c input is initialized. The physical word holds the address throughout,
so the four identities combine coherently without dropping address-dependent
phases. The literal probe sequence $`V^\dagger,Z_{c_0},V`$ returns c
and every mask on every input, including arbitrary reference correlations.
The source itself need not return c before that completed sequence.

A native controlled H uses
$`R=SHTHS^\dagger`$, $`RZR^\dagger=H`$, and the literal word
$`R\,\mathrm{CZ}\,R^\dagger`$. With CZ emitted as H-CX-H it has two
T gates, two T bands and 17 total gates, also for its actual adjoint.
Parallelizing only the disjoint target T bands in the $`K_1`$ control gives
34 total gates and four T gates in two bands. The $`K_2`$ control adds the
four Cliffords of F and its adjoint, giving 21 total gates and two T
gates in two bands. Thus W has six T gates, four bands and 55 total gates
in either direction. Its shared control is never used twice in a T band.

The tree query $`J_x`$ has the three bilinear terms
$`a c_1`$ into $`t_0`$, $`b c_2`$ into $`t_0`$ and $`a c_3`$ into $`t_2`$, besides
the four ordinary c-to-t CNOTs. Hadamards on $`t_0,t_2`$ turn the bilinear
terms into one returned-dirty cubic batch. Its three masks give 24 T
gates, eight bands and 56 total gates for $`J_x`$ and its actual inverse.
C adds eight Cliffords, and the initial $`H_{c_0}`$ adds one. Therefore

```math
\begin{aligned}
T_V&=30,&D_T(V)&\le12,&G_{\rm all}(V)&\le165,\\
w_{\rm extra}&=7,&a_{\rm clean}&=0,&\epsilon&=0.
\end{aligned}                                                   \qquad(10)
```

The actual emitted lengths are 165 forward and 155 inverse. Extra width
is the four arbitrary c wires and three returned masks, beyond the two
address and four target-bank wires. In particular, choosing the input
probe on $`c_0`$ does not make it clean. The complete phase and its actual
adjoint both have 60 T gates, at most 24 bands and 322 total gates,
including two S gates for the central $`Z_{c_0}`$. Every mask is returned by
the cubic batch itself; the complete phase returns the c bank by (9).

The [exact fixtures](../../tests/test_two_address_interleaved_pilot.py)
check the native controlled-H word on all its columns, its actual inverse,
the signed decoder and all 16 generator products, all four complete
source fibers, every active cubic input and arbitrary mask value, and
the physical T-band and inverse ledgers. These are complete finite
operator checks, not an asymptotic proof.

## 5. What this pilot does and does not advance

The construction passes the noncommuting-family test and pays a genuine
address-dependent query after earlier nonlinear layers. It shows that
the changed interface is realizable at its first nontrivial size.
Controlled-H stages are still paid Clifford fibers with bounded-degree
faithful Pauli action; the applicable rank screens retain their hypotheses.

No uniform constant-stage extension or fusion of these physical T words
has been established here. Building a separate controlled Clifford for
each successive address bit restores a serial address-dependent schedule.
The next mechanism must supply an actual cheaper uniform fusion, or use
a genuinely nonlinear intermediate family whose whole set of products
and commutators matches the target. The existing
[signed third-level transport screen](NATIVE_CLIFFORD_TRANSPORT.md)
also remains applicable whenever its specified prefix factorization and
Clifford comparison endpoints hold; calling an intermediate reflection
non-Clifford does not by itself remove those premises. Merely decorating
this finite source with a nonlinear word that later cancels would not
make the required mechanism change.
