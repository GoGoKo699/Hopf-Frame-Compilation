# A paid Majorana Walsh probe and its restricted source cost

[Complete probe interface](SYMPLECTIC_MODULE_QUERY_BARRIER.md#6-universal-pauli-probes-contain-a-complete-dirty-indicator) ·
[Conditional word-query transfer](SYMPLECTIC_MODULE_QUERY_BARRIER.md#7-the-complete-r98-word-query-interface-follows-conditionally) ·
[Retained Majorana source model](SOURCE_T_DEPTH.md)

Simultaneous real diagonalization of XOR translations gives an exact
independent-Pauli probe without implementing a selected permutation of
physical qubits. All addressing between its fixed source calls is
Clifford. The standard source nevertheless has count and depth too large
for the improved query interface. Within the specified Majorana-linear
model its exact T-count is $`mp/2`$ and its exact T-depth is p, where
$`m=2^p`$. These are costs of that separately implemented source, not
lower bounds for a jointly compiled query or unrestricted Clifford+T.

## 1. The full-carrier probe identity

Take $`p\ge2`$, $`m=2^p`$ and $`w=m/2`$ arbitrary carrier qubits.
Use their m Jordan--Wigner Majoranas

```math
\gamma_{2j}=Z_0\cdots Z_{j-1}X_j,\qquad
\gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j.
```

Their binary Pauli labels are independent: their symplectic Gram matrix
is $`I+\mathbf1\mathbf1^{\mathsf T}`$, whose square is I for even m.
Let W be an actual native unitary with Majorana conjugation matrix

```math
H_{uv}=m^{-1/2}(-1)^{u\cdot v},\qquad
W\gamma_vW^\dagger=\sum_uH_{uv}\gamma_u.                 \qquad(1)
```

The construction below supplies W and its actual inverse. Its common
scalar is immaterial to conjugation and is canceled by actual adjoints;
no scalar is discarded from the resulting probe word.

For each address bit i choose the positive Hermitian tensor Pauli P_i
whose conjugation negates exactly the modes with $`u_i=1`$. Such a Pauli
exists because the Majorana labels form a basis of the full Pauli-label
space. More concretely, the product of those m/2 Majoranas has the required signs:
m/2 is even, and its scalar can be chosen to make the representative
positive Hermitian. Put P_x equal to the ordered product of
$`P_i^{x_i}`$ in a fixed order. Implement it by controlled tensor
Paulis, with the original address as controls. This is one actual
Clifford word. The P_i need not commute, so $`P_x^\dagger`$ below is
its actual adjoint, not a replacement by P_x.

Conjugation by P_x has the diagonal Majorana action
$`D_x=\mathrm{diag}((-1)^{u\cdot x})`$. Therefore

```math
H D_x H\,e_v=e_{v\oplus x}.
```

With $`G_x=WP_xW^\dagger`$, the complete source

```math
M_x=G_x\gamma_0G_x^\dagger
    =WP_xW^\dagger\gamma_0WP_x^\dagger W^\dagger
    =\gamma_x                                                   \qquad(2)
```

is the independent-Pauli multiplexor required by the linked interface.
It preserves the address, uses no initialized work and has literal
equality on every carrier input and reference. A fermionic mode
translation need not equal a permutation of physical qubit wires;
(2) only needs its action on the specified probe. Thus it avoids that
extra implementation requirement without supplying a free basis change.

## 2. A native Walsh word with a quadratic Clifford ledger

For a pair of modes a,b define

```math
R_{ab}=\exp\!\left(\frac{i\pi}{8}(i\gamma_a\gamma_b)\right).
```

Its Majorana action sends $`\gamma_a`$ to
$`(\gamma_a+\gamma_b)/\sqrt2`$ and $`\gamma_b`$ to
$`(-\gamma_a+\gamma_b)/\sqrt2`$. One signed Pauli T/T-dagger
rotation implements this action. Its native common scalar is retained
by the emitted word.

At bit stage k, set $`d=2^k`$ and pair every a with $`a+d`$ inside
each block of 2d modes, taking the d first-half indices as a. First
apply the positive Pauli that flips the second member of every pair;
then apply all R rotations. The resulting two-mode block is
$`2^{-1/2}\left(\begin{smallmatrix}1&1\\1&-1\end{smallmatrix}\right)`$.
The p successive bit stages give (1).

The rotations in one matching have independent commuting Pauli axes.
They can be put in one **physical** T layer with an explicit paid
Majorana routing word. Interleave the first and second halves of each
2d block by adjacent mode swaps. The number of adjacent swaps is

```math
N_k=\frac{m}{2d}\frac{d(d-1)}2=\frac{m(d-1)}4.
```

An adjacent signed Majorana swap or its actual inverse is a Clifford:
an on-site pair costs at most three S gates; a pair crossing a qubit
boundary is a $`\pi/4`$ rotation about $`X_jX_{j+1}`$, realizable
with four H gates, two CNOTs and at most three S gates. Thus either
orientation costs at most nine elementary Cliffords. Track the route's
signs, apply one T or T-dagger on each of the w physical wires, and
actually unroute. This costs at most $`18N_k`$ Cliffords. It does
not assume that a generic Clifford diagonalization is free at every
stage.

The stage's positive Pauli costs at most 8w Cliffords: X costs four,
Z two, and $`Y=SXS^\dagger`$ eight in the positive-S alphabet.
Its actual inverse has the same bound after reducing adjacent S powers.
Since $`\sum_kN_k=m(m-1-p)/4`$, common forward/inverse bounds are

```math
\begin{aligned}
T_W&=mp/2,\\
G_W&\le\frac92m(m-1-p)+4mp\lt\frac92m^2,\\
D_T(W)&=p,\\
b_W&=0.
\end{aligned}                                                   \qquad(3)
```

Each routed matching is a physical band on distinct targets. The p
bands are sequential; commuting ideal bit transforms are not declared
simultaneous native words. For $`p\ge2`$, H has determinant plus one,
and all stage masks flip an even number of modes. The native W is
parity preserving. This will also keep the extra-mode accounting below
literal.

A controlled positive tensor Pauli costs at most $`5w`$ elementary
Cliffords, including controlled Y and its actual inverse. Hence the
address mask and its inverse cost at most $`5pw`$. Applying (2)
chronologically as

```math
W^\dagger,\ P_x^\dagger,\ W,\ \gamma_0,
W^\dagger,\ P_x,\ W
```

gives common source/inverse bounds

```math
\begin{aligned}
T_M&=2mp,\\
G_M&\le4G_W+10pw+4=O(m^2),\\
D_T(M)&\le4p,\\
w&=m/2,\qquad b_M=0,\qquad a_{\rm clean}=0.
\end{aligned}                                                   \qquad(4)
```

The four W calls, both mask orientations and the native four-Clifford
X word for $`\gamma_0`$ are all charged. This paid source misses the
linear T-count and sublogarithmic address-depth requirements. The
retained R91 query already supplies a better baseline for the desired
probe task; (4) is not a query improvement.

## 3. Probe-only transport still determines the whole Walsh matrix

Within this sign-mask arrangement, replacing a full translation
specification by only its selected-probe action does not permit an
arbitrary flat-column loader. Let a fixed Majorana-linear unitary have
real orthogonal matrix O and suppose

```math
O D_x O^{\mathsf T}e_0=e_x\qquad\hbox{for every }x.
```

Set $`a=O^{\mathsf T}e_0`$. Then
$`D_xa=O^{\mathsf T}e_x`$, so $`O_{xj}=a_j(-1)^{j\cdot x}`$.
Orthogonality of all rows forces every nonconstant Walsh coefficient
of the weights $`a_j^2`$ to vanish. Their sum is one, and thus
$`a_j^2=1/m`$. Consequently

```math
O=H\mathrm{diag}(\epsilon_j),\qquad\epsilon_j\in\{1,-1\}.
                                                               \qquad(5)
```

Allowing output signs or relabeling adds signed row permutations and
does not alter the conclusion relevant to count or depth. This argument
concerns the specified m-mode sign masks and a fixed linear Majorana
conjugator. General non-Gaussian conjugators and different readout
arrangements do not satisfy its premises automatically.

## 4. Sharp cost in the Majorana-preserving native model

The model permits arbitrary native Clifford blocks whose complete
conjugation preserves the span of all physical Majoranas, alternating
with physical T/T-dagger bands. Such a Clifford block is a signed
permutation of the Majoranas: a Clifford takes each Pauli to a signed
Pauli, and distinct Majoranas are orthogonal Pauli basis elements.
A physical T gate rotates one on-site pair by $`\pi/4`$. Additional
arbitrary modes may mix internally and must return at the end.

The count argument applies the existing matrix-entropy method of
[Ailon, arXiv:1305.4745v1, Theorem 2.1 and equation (2.1)](https://arxiv.org/html/1305.4745v1).
For the real orthogonal conjugation matrix define

```math
\Phi(O)=-\sum_{i,j}O_{ij}^2\log_2 O_{ij}^2,\qquad0\log_20=0.
```

Signed permutations preserve this quantity. A two-row rotation preserves
the sum $`r_j`$ of the two squared entries in each column. Their entropy
contribution can change by at most $`r_j`$; summing gives at most two.
Identity has entropy zero, whereas H has entropy $`m\log_2m=mp`$.
Thus the native specialization gives $`T_W\ge mp/2`$. This is an
application of Ailon's method, not a new entropy lower-bound principle.

For depth, every matching layer can at most double the largest number
of nonzero entries in a row. Signed permutations do not increase it.
After D physical T bands this number is at most $`2^D`$, whereas every
Walsh row has m nonzero entries. Therefore $`D_T(W)\ge p`$. The
[retained denominator argument](SOURCE_T_DEPTH.md#1-gate-class-and-inherited-denominator-bound)
gives the same depth bound from the entry $`2^{-p/2}`$.

Both bounds survive arbitrary returned Gaussian helper modes. W is even,
so the target full Majorana action is $`H\oplus I`$; the added identity
has zero entropy and the original rows still have m nonzero entries.
Intermediate mixing with the helper modes changes neither argument.
The word in Section 2 attains both bounds, establishing exactly the
restricted source costs stated in (3).

Arbitrary Clifford interlayers can leave the Majorana span. The entropy
and row-support arguments then cease to describe the circuit. Moreover,
a complete probe circuit can synthesize the whole word (2) jointly
without separately realizing W. Neither exit is excluded here. The first
uncompiled improvement is such a native basis change or complete probe
with all simultaneous source budgets paid; ordinary Gaussian butterfly
reordering cannot supply it.

## 5. A flat non-Pauli escape and a constant unrestricted bound

One can ask for the single Hermitian unitary

```math
A=m^{-1/2}\sum_j\gamma_j,
\qquad A_x=P_xAP_x^\dagger=W\gamma_xW^\dagger.
```

This smaller source has a paid linear-count construction. At stage k,
split each currently supported mode $`0\le u\lt2^k`$ into the pair
$`(u,u+2^k)`$ with the signed R rotation of Section 2. Use that
section's full matching route, but apply T/T-dagger only to the first
$`2^k`$ physical pairs, and actually unroute. The tracked route signs
select the orientation sending $`\gamma_u`$ to
$`(\gamma_u+\gamma_{u+2^k})/\sqrt2`$; no stage sign mask is needed.
The other pairs receive identity. Induction gives
$`F\gamma_0F^\dagger=A`$ after p stages, with common forward/inverse
bounds

```math
T_F=m-1,\qquad D_T(F)=p,\qquad
G_F\le18\sum_kN_k=\frac92m(m-1-p),\qquad b_F=0.
```

Thus the actual word $`F\gamma_0F^\dagger`$ implements A with
$`T\le2(m-1)`$, $`D_T\le2p`$, $`G\le2G_F+4`$ and no initialized
work. F is a full-carrier unitary, and its native scalar cancels against
its actual adjoint. This is a linear T-count source, not a faster complete
table query.

Anticommutation of the Majoranas gives $`A^2=I`$ and
$`A_xA_y+A_yA_x=2\delta_{xy}I`$. Consequently

```math
A_xA_yA_xA_y=-(-1)^{[x=y]}I.
```

This four-factor identity supplies an equality phase up to a fixed
minus sign, on every carrier input. A native implementation of A alone
need not implement the whole Walsh matrix, so the $`mp/2`$ source count
conclusion in Section 4 does not apply to it. The A_x, however, are non-Pauli
operators, not independent physical Pauli labels. The universal table
readout in the linked interface also needs their dual probes under
joint control of arbitrary dirty table bits Y.

Explicitly, let the physical Paulis $`\Delta_j`$ be symplectically dual
to the Majorana labels, so conjugation by $`\Delta_j`$ negates just
$`\gamma_j`$. The corresponding dual operators for the A_x frame are
$`\widetilde\Delta_j=W\Delta_jW^\dagger`$. The first unpaid
complete-input operation is the ordered controlled product

```math
Q_Y=\prod_j\widetilde\Delta_j^{Y_j}
   =W C_YW^\dagger,\qquad
C_Y=\prod_j\Delta_j^{Y_j}.
```

C_Y is Clifford, but its non-Pauli conjugate Q_Y has not thereby been
compiled. The displayed standalone wrapper pays mp T gates again.
Joint native synthesis of Q_Y remains an explicit alternative; neither
the equality identity nor preparation of A alone supplies this arbitrary
dirty readout.

The desired complete interface is explicit: since conjugation by Q_Y
negates $`A_x`$ precisely when $`Y_x=1`$,

```math
A_xQ_YA_xQ_Y^\dagger=(-1)^{Y_x}I.
```

Hadamards on every Y wire before and after this commutator turn its phase
into $`Y\mapsto Y\oplus e_x`$, returning the arbitrary carrier and
preserving the coherent address. This is a complete dirty indicator,
provided both actual source orientations have been compiled. The current
A word has depth at most 2p, and the standalone Q_Y wrapper has mp T
gates; their stated costs fail the desired simultaneous sublogarithmic
address depth and linear T-count. The commutator does not remove either
cost.

There is also a small unrestricted consequence of the existing
[four-layer operator rigidity theorem R75](FOUR_LAYER_OPERATOR_RIGIDITY.md).
For every j,

```math
A\gamma_jA=-\gamma_j+\frac2{\sqrt m}A.
```

Choose a Hermitian Pauli R anticommuting with $`\gamma_j`$ and put
$`V=RA`$. Then

```math
V\gamma_jV^\dagger-\gamma_j=\frac2{\sqrt m}RAR,
\qquad
\|V\gamma_jV^\dagger-\gamma_j\|=\frac2{\sqrt m}\ne0.
```

R75 says that an at-most-four-layer native circuit cannot change a
Pauli by a nonzero operator norm smaller than
$`\rho=\sin(\pi/8)/\sqrt2`$. Thus

```math
m\gt32+16\sqrt2
\quad\Longrightarrow\quad D_T(A)\ge5.
```

For powers of two this holds at $`m\ge64`$, including arbitrary
returned helpers, since the same nonzero full-space norm holds after
tensoring with their identity. Any full-space Pauli-conjugating
preparation $`F\gamma_0F^\dagger=A`$ therefore needs at least three
T layers: otherwise that displayed word would implement A in at most
four. These are constant lower bounds for this particular flat unitary
and its preparation. They neither grow with p nor establish a lower
bound for unrestricted independent-Pauli probe synthesis.

## 6. Exact finite evidence

The [focused fixtures](../../tests/test_majorana_walsh_query.py) retain
literal phases and actual inverses. A two-qubit W has four T gates,
30 elementary Cliffords and two physical T bands; every Majorana image
is checked against H_4. Its complete two-address probe has 16 T gates
and 142 Cliffords in both orientations and is checked on all 16 physical
input columns. The two address masks anticommute: replacing the ordered
mask's actual adjoint by itself gives a minus sign on the address-11
sector, which the fixture detects. Binary stage checks cover the routing
counts and Walsh action at several larger sizes. A small flat-frame
fixture checks the anticommutators, reflection identity and the fixed
minus sign in the equality phase. These finite checks
support the emitted words; the all-size identities and restricted
optimality statements follow from the proofs above.
