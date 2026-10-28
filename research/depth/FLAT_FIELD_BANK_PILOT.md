# A flat non-Clifford field bank with a three-band decoder

[Field query](FIELD_SHEAR_INTERTWINER.md) ·
[Complete copy grading](DIRTY_COPY_GRADING.md) ·
[Encoded copy contract](ENCODED_DIRTY_COPY_INTERFACE.md) ·
[Earlier Clifford-seed pilot](FIELD_COPY_CLIFFORD_PILOT.md)

The two-bit field query admits an independent commuting bank of four
genuinely non-Clifford reflections. Its seed has a three-T, two-band
preparation. A six-T, three-band decoder maps the complete bank to
physical Paulis, giving a twelve-T, six-band copy. Used directly on the
physical output bank, the decoder gives an exact four-address indicator
with **72 T gates, 28 T bands, 350 total gates and three returned dirty
wires**. Every principal and helper input is arbitrary. All T-band counts
are upper bounds from the displayed schedules; no minimal-depth claim is made.

A further Clifford wall completes the same six T gates to a full
field-to-wire-translation intertwiner. This is a finite construction,
not an all-size compiler or a new general frame-depth bound. The seed
also has a one-band copy error floor of $`1/(2\sqrt2)`$: the successful
readout lies outside that narrower grammar.

## 1. Physical coordinates and a native flat seed

Use four carrier wires $`(0,1,2,3)=(u_0,u_1,v_0,v_1)`$. Address labels
0,1,2,3 are binary, with addition by XOR. The retained field query is

```math
Q_x|u,v\rangle=|u,v+M_xu\rangle,\qquad
M_0=0,\quad M_1=I,\quad
M_2=\begin{pmatrix}0&1\\1&1\end{pmatrix},\quad M_3=M_1+M_2.
```

For a Hermitian Pauli P, retain the native scalar in
$`T(P)=\exp(i\pi(I-P)/8)`$. On two wires put

```math
R_*=(X_0+Y_0+Z_0X_1-Z_0Z_1)/2,\qquad
V_*=T(X_0X_1)T(Y_0Z_1)T(Z_0).
                                                               \qquad(1)
```

The four signed Pauli summands anticommute pairwise, so $`R_*^2=I`$.
Direct conjugation gives $`R_*=V_*X_0V_*^\dagger`$. The outer two
rotation axes commute and are independent. A native Clifford U that
maps $`Z_0,Z_1`$ to those axes is the chronological word

```text
U: H(1), CX(1,0), S(0), S(0), S(0), H(0), H(1).
```

Thus $`V_*`$ has the chronological word
$`T(0),U^\dagger,[T(0),T(1)],U`$. The bracket denotes one physical
T band on distinct wires. Reducing the three reversed S gates modulo
four gives 15 total gates in either direction, with three T gates in
two bands, two arbitrary wires and no helpers.

The seed is not Clifford. The inner Clifford
$`T(Z_0)X_0T(Z_0)^\dagger=(X_0+Y_0)/\sqrt2`$ sends the outer
axis $`X_0X_1`$ to $`Y_0X_1`$, which is neither outer axis.
The [exact one-band axis-matching theorem](ONE_BAND_AXIS_MATCHING.md)
therefore excludes a Clifford result for the outer conjugation.

Embed the two-wire preparation on the first two field wires and apply
the following fixed chronological Clifford F:

```text
S(1), H(0), H(1), CX(0,2), S(2), S(0), CX(1,0),
H(2), S(0), CX(2,0), H(2), S(2), H(1).
```

Then $`V=FV_*`$ prepares the seed

```math
R=VX_0V^\dagger
 =(Z_0X_1Y_2+Z_0Z_2-Y_0Z_1Z_2-Y_0Y_1Z_2)/2.       \qquad(2)
```

Its preparation has three T gates in two bands, with 28 total forward
gates and 38 in the actual inverse after expanding S-dagger. No state
is initialized. This is a preparation ledger; an emitted R also pays
the preparation's actual inverse and the central Pauli.

For concise exact tables, use the Hermitian Pauli convention

```math
P_{a+16b}=i^{\mathrm{wt}(a\mathbin\&b)}X^aZ^b,
\qquad a,b\in\mathbb F_2^4.
```

The complete field orbit $`R_x=Q_xRQ_x^\dagger`$ is:

| Address x | $`2R_x`$ |
| --- | --- |
| 0 | $`P_{86}+P_{80}-P_{113}-P_{115}`$ |
| 1 | $`P_{78}+P_{64}-P_{101}-P_{111}`$ |
| 2 | $`P_{122}+P_{112}-P_{89}-P_{87}`$ |
| 3 | $`P_{98}+P_{96}-P_{77}+P_{75}`$ |

These signs follow from literal field CNOT conjugation of (2). An independent
non-Clifford certificate is
$`RX_1R=(-P_{37}-P_{39}+P_{35}-P_{33})/2`$, which is not a Pauli.
The decoder below proves their complete commuting-bank structure;
pairwise orthogonality alone would not do so.

## 2. Six rotations decode the entire bank

Let D be the chronological product of the following three bands:

```math
[T(P_{33}),T(P_{41})],\qquad
[T(P_2),T(P_{14})],\qquad
[T(P_{35}),T(P_{39})].                              \qquad(3)
```

Each pair is commuting and independent. For a literal physical
implementation, the following chronological Clifford words map the
specified single-wire Z axes to the corresponding pair:

| Band | Clifford $`K_i`$ | Native T wires |
| --- | --- | --- |
| 1 | $`\mathrm{CX}(0,3),H(0),H(3),H(1),\mathrm{CX}(0,1),H(1)`$ | 0,3 |
| 2 | $`\mathrm{CX}(1,2),H(1),H(2),\mathrm{CX}(2,3)`$ | 1,2 |
| 3 | $`\mathrm{CX}(1,2),H(1),S(1),\mathrm{CX}(1,0),H(2)`$ | 1,2 |

Emit $`K_i^\dagger,[T,T],K_i`$ for each row, in order. With
S-dagger expanded into three S gates, D and its actual inverse each
have six T gates, three physical bands and 38 total gates. There are
four arbitrary wires and no private helper.

Applying the exact rule
$`T(A)PT(A)^\dagger=(P-iAP)/\sqrt2`$ when A anticommutes with P,
and leaving P unchanged otherwise, gives

```math
DR_0D^\dagger=-P_{113},\quad
DR_1D^\dagger=P_{64},\quad
DR_2D^\dagger=-P_{89},\quad
DR_3D^\dagger=P_{96}.                               \qquad(4)
```

Call the four signed Paulis on the right $`B_0,\ldots,B_3`$.
They commute and their binary labels are independent. Every nonempty
subset product is therefore a nonidentity Pauli and has trace zero.
Conjugating back proves the same facts for the $`R_x`$, including
$`R_0R_1R_2R_3=X_3`$. Their sixteen joint sign projectors have rank
one. Every word reduces to a subset product, so the full
[copy-grading criterion](DIRTY_COPY_GRADING.md#1-the-exact-full-space-criterion)
holds for every word length and order.

### The three-band word is a coherent Clifford selection

There is a small algebra behind (3). Put
$`A=X_0Z_1`$, $`B=X_1`$, $`C=X_0Y_1=-iAB`$.
The pairs in (3) are $`(A,AX_3)`$, $`(B,BX_2X_3)`$ and
$`(C,CX_2)`$. The control Paulis $`X_2,X_3`$ commute with A,B,C.
For $`L=\pm1`$ in a control sector,

```math
T(A)T(AL)=e^{i\pi/4}
\begin{cases}
e^{-i\pi A/4},&L=+1,\\
I,&L=-1.
\end{cases}                                                     \qquad(5)
```

Writing $`S_A=e^{-i\pi A/4}`$, the Pauli algebra gives
$`S_CS_BS_A=S_B`$. Up to the one common scalar $`e^{3i\pi/4}`$,
D therefore selects $`S_B`$ when the two control signs agree,
$`S_C`$ for $`(X_2,X_3)=(+1,-1)`$, and $`S_A`$ for
$`(-1,+1)`$. This is a complete coherent operation on the control
sectors, not a promise about their inputs or a measurement.

## 3. A physical bank decoder and complete copy

The following chronological Clifford L maps $`B_j`$ to $`Z_j`$
with literal positive signs:

```text
S(0), S(0), S(0), H(0), CX(1,0), CX(2,0),
H(0), S(0), S(0), H(0),
CX(1,2), CX(2,1), CX(1,2), H(3),
CX(0,2), CX(3,2), CX(0,3), CX(1,3), CX(2,3).
```

It has 19 forward gates and 17 in its actual inverse. The complete
bank decoder $`E=LD`$ consequently obeys

```math
ER_jE^\dagger=Z_j,\qquad
T(E)=6,\quad D_T(E)\le3,\quad
G_{\rm all}(E)=57,\quad G_{\rm all}(E^\dagger)=55.       \qquad(6)
```

On a separate arbitrary four-qubit output bank Y, let
$`C_Z=\prod_jH_j\mathrm{CX}(j;Y_j)H_j`$, so that
$`C_ZZ_jC_Z^\dagger=Z_jX_{Y_j}`$. The chronological copy word
$`E,C_Z,E^\dagger`$ is
$`C=(E^\dagger\otimes I_Y)C_Z(E\otimes I_Y)`$ and satisfies

```math
C(R_j\otimes I_Y)C^\dagger=R_j\otimes X_{Y_j}.       \qquad(7)
```

Both directions use twelve T gates, six bands and 124 total gates on
the four carrier and four output wires. There are no extra helpers or
initialized inputs. The actual reversed native word is used even
though the ideal copy is an involution. Equation (7) includes the full
carrier action, not only the output dual signatures, and therefore
supplies complete return in the
[encoded-copy echo](ENCODED_DIRTY_COPY_INTERFACE.md#1-cancel-an-encoding-only-on-the-returned-dirty-carrier).

## 4. The cheaper direct four-address indicator

This commuting bank permits using E directly on the physical output
register. Put $`S_x=EQ_xE^\dagger`$. Equations (2) and (6) imply

```math
S_xZ_0S_x^\dagger=Z_x.                              \qquad(8)
```

Only (8) is needed for this implementation; the stronger full
intertwinement in the next section uses a different Clifford suffix.
For the coherent $`S=\sum_x|x\rangle\langle x|\otimes S_x`$,
the chronological word

```math
H_Y,\quad S^\dagger,\quad Z_{Y_0},\quad S,\quad H_Y              \qquad(9)
```

is exactly $`|x,Y,h\rangle\mapsto|x,Y\mathbin\oplus e_x,h\rangle`$.
The retained [three-cube query](FIELD_COPY_CLIFFORD_PILOT.md#4-three-cube-coherent-query-and-the-complete-indicator)
has 24 T gates, eight bands, 58 total gates and three arbitrary
returned masks. Its actual inverse has the same bounds. Use physical
address wires 0,1, output wires 2 through 5 and masks 6 through 8.
Lift E to the four output wires; its word introduces no other work.

There are four E/adjoint appearances and two query appearances in
(9). Including the eight output Hadamards and the two-S implementation
of Z gives, in either direction,

```math
T_I=4\cdot6+2\cdot24=72,\qquad
D_T(I)\le4\cdot3+2\cdot8=28,
```
```math
G_{\rm all}(I)=2(57+55)+2\cdot58+8+2=350,
\qquad w=9,\quad b_{\rm dirty}=3,\quad a_{\rm clean}=0.          \qquad(10)
```

The component equalities prove the full operator identity on every
coherent address, output, mask and external reference. Exact native
execution of all 512 basis columns in both directions separately
checks allocation, chronology, literal phases and return. This
improves the retained finite direct-decoder ledger of 300 T gates,
160 bands and 892 total gates using the same three masks. It is not
an optimality claim about the four-address problem.

## 5. The same six T gates also give a full intertwiner

Let $`Q'_x=DQ_xD^\dagger`$. Exact conjugation of the eight physical
Pauli generators shows that every $`Q'_x`$ is Clifford. A complementary
orbit of signed Paulis is

```math
C_0=P_{240},\quad C_1=P_{198},\quad C_2=P_{224},\quad C_3=P_{218},
\qquad Q'_xC_yQ_x'^\dagger=C_{x+y}.                 \qquad(11)
```

The C's commute, and $`B_jC_k=(-1)^{[j=k]}C_kB_j`$. They therefore
complete the B's to a Pauli basis. After L, the chronological suffix

```text
J: S(1)-dagger, S(3), H(3), CX(0,3), CX(2,3), H(3)
```

fixes every $`Z_j`$ and maps $`LC_jL^\dagger`$ to $`X_j`$.
Expand S-dagger into three S gates. This adds eight gates in either
direction. Hence
$`E_{\rm full}=JLD`$ has six T gates in at most three bands and 65 forward
or 63 inverse total gates, on the same four arbitrary wires.

Both B and C are permuted by $`Q'_x`$ with exactly the XOR action.
Thus $`E_{\rm full}Q_xE_{\rm full}^\dagger`$ has the same
conjugation on every physical X and Z as the wire translation $`W_x`$.
The two operators differ by a scalar only. For nonzero x both have
trace four, by the field character and the two-swap character;
the scalar is therefore one. At x=0 both are identity. This proves
the complete literal identity

```math
E_{\rm full}Q_xE_{\rm full}^\dagger=W_x
\qquad(x=0,1,2,3).                                 \qquad(12)
```

The finite native fixture also verifies (12) on all sixteen columns
for all four addresses and checks the actual inverse. No gauge or
unknown global phase is omitted. The shorter E in (6) remains the
word used for the cheaper ledger (10).

## 6. A one-band copy cannot read this seed accurately

For this finite seed, no Clifford J can implement any one of the four
singleton sign patterns

```math
JR_xJ^\dagger=(-1)^{[x=j]}R_x\qquad(x=0,1,2,3).      \qquad(13)
```

Here is a complete finite verification specification. Each R has four
distinct Pauli labels with coefficient magnitude one half. A Clifford
satisfying (13) must permute its four support labels, with the signs
fixed by its coefficient signs and the requested singleton. There are
24 permutations per address. Extend these candidates one address at
a time, retaining only assignments that preserve all binary label
relations, symplectic pairings and literal Pauli multiplication, with injectivity on the support span.

For the ordered section $`P(a,b)=X^aZ^b`$, multiplication contributes
$`2b\cdot a'`$ modulo four. If a source Hermitian label v maps to w,
the prescribed ordered-section phase is

```math
2(s_x+c_v+c_w)+\mathrm{wt}(a_w\mathbin\&b_w)
                 -\mathrm{wt}(a_v\mathbin\&b_v)\pmod4,          \qquad(14)
```

where $`s_x=[x=j]`$ and the c's record the two coefficient signs.
Exact subgroup closure with this multiplication law checks every induced
relation in the seven-dimensional union of the support spans. For
each singleton, the surviving restricted maps after addresses 0,1,2,3
number exactly $`24,8,4,0`$. This exhausts all possible Clifford
actions on the supported Paulis. The identity sign pattern has four
surviving maps, providing a positive consistency check. Additional
identity helpers cannot repair the contradictory relations: each
supported Pauli still has to map to the specified carrier-only support.

There is a quantitative consequence. Two distinct four-term Pauli expansions
whose coefficients all have magnitude one half have Hilbert--Schmidt inner product at most
$`3/4`$. Their normalized Hilbert--Schmidt distance, and hence their
operator distance, is at least $`1/\sqrt2`$. If a copy has complete
error eta, an output dual $`J=C^\dagger Z_{Y_j}C`$ obeys

```math
\|J(R_x\otimes I)J^\dagger-(-1)^{[x=j]}R_x\otimes I\|
\le2\eta.
```

When J is Clifford and $`2\eta\lt1/\sqrt2`$, every one of these four
relations must therefore be exact, contradicting the enumeration.
Every one-band native copy has Clifford Pauli output duals, so

```math
\eta_{\rm copy}\ge\frac1{2\sqrt2}.                  \qquad(15)
```

This bound concerns the seed (2) in its stated field coordinates and
arbitrary identity helpers. It does not exclude other flat seeds,
other embeddings, two-band copies or the six-band copy above.

The [native fixtures](../../tests/test_flat_field_bank_native.py) preserve
the preparation, decoder, copy, direct indicator and actual inverse words.
The [grading and dual fixtures](../../tests/test_flat_field_bank_duals.py)
retain the exact support-map enumeration, literal multiplication phases,
all subset traces and the full-intertwiner checks.

## 7. Scope and the next construction obligation

The construction pays an entire compatible finite bank and its
physical readout. It avoids using a nonzero linear kernel or a
collection of output duals as a substitute for a complete carrier
identity. The earlier noncommuting Clifford pilot remains a distinct
graded algebra; its proofs and native words are preserved in their
own home.

The next substantive task is an all-size analogue of the coherent
selection in (5), jointly with its field orbit and full copy or
decoder. A proposed extension must prove the complete grading and
operator identities, emit actual forward and reverse native words,
and meet the common count, width and precision budgets. Recursing
through independently paid field recouplers retains a depth cost at
every level; the
[recursive-decoder limitation](RECURSIVE_FIELD_SHEAR_DECODER.md#4-what-remains-to-be-compiled-and-what-the-recursion-costs)
still applies. The finite six-T identity does not supply that missing
uniform construction, and R96 remains unchanged.

The [projective recoupler](PROJECTIVE_FIELD_RECOUPLER.md) uses this full
six-T intertwiner on both half banks to pay the next finite size. Its
eight-address decoder has 172 T gates and at most 83 bands without extra
helpers. Joint reflection phases and zero-fiber freedom reduce that word;
the exponential projective-line count still prevents uniform scaling.

<a id="uniform-extension-axis-test"></a>

### A necessary test for extending the full intertwiner

For the all-size field query on m core wires and h identity helpers,
let $`\mathcal F`$ be its common fixed Pauli-label space. It consists
of the $`Z_u,X_v`$ labels and all helper Paulis, so
$`\dim\mathcal F=m+2h`$. Absorb the Clifford walls of a proposed
native full intertwiner into
$`E=L\prod_jT_{\epsilon_j}(P_j)`$, and let A span its absorbed
Pauli rotation axes. The complete error is

```math
\eta=\max_x\|EQ_xE^\dagger-W_x\otimes I_h\|.
```

If $`\eta\lt1`$, then

```math
\mathrm{rank}\bigl(\langle\ ,\ \rangle|_{A\times\mathcal F}\bigr)
\ge m-2,\qquad T(E)\ge m-2.                         \qquad(16)
```

To prove this, every label in $`\mathcal F\cap A^\perp`$ supplies
a Pauli fixed literally by all source queries and commuting with
every rotation. Its image under E is therefore the Pauli supplied
by L alone. It is fixed by every $`EQ_xE^\dagger`$, so its distance
from its conjugate under $`W_x\otimes I_h`$ is at most $`2\eta`$.
An involutive wire permutation sends every Pauli to a commuting
Pauli: in binary labels the relevant symplectic pairing vanishes
because the permutation matrix is symmetric. Any two unequal
commuting Hermitian Paulis have operator distance two, including
opposite signs. Thus the image is fixed exactly when $`\eta\lt1`$.

The common target fixed-label space consists of $`X_{\rm all}`$,
$`Z_{\rm all}`$ and helper Paulis, of dimension $`2+2h`$.
Injectivity of L gives
$`\dim(\mathcal F\cap A^\perp)\le2+2h`$, proving (16).
This is an extensive count and axis-rank requirement, not a growing
T-depth bound. It concerns full intertwinement; the weaker one-probe
decoder contract in (8) does not impose it.

For example, the two control Paulis of each disjoint copy of the
four-wire gadget (3) commute with all its axes. A direct tensor
extension that leaves these controls inside the original core-fixed
space $`\mathrm{span}\{Z_u,X_v\}`$, with untouched identity helpers,
retains at least $`m/2+2h`$ such directions. For $`m>4`$ it
fails the full-intertwiner test below error one. An extension must
change that alignment or couple the cells. This does not reject
arbitrary three-band designs or a weaker probe-only construction. The
[fixed-Pauli fixture](../../tests/test_flat_field_bank_duals.py) checks
the two-control kernel, a failed eight-wire disjoint extension and
the target Pauli-distance premise exactly.

The stronger [Pauli-domain theorem](FIELD_DECODER_PAULI_DOMAIN.md) requires
absorbed-axis span at least $`2m-2`$ at p at least three, below error
one-third for identity or residual-field helpers and below $`\sqrt2/6`$
for arbitrary Clifford helper motion. It excludes independent finite cells
even after arbitrary joint outer Clifford walls. The separate
[two-band theorem](TWO_BAND_FIELD_DECODER_OBSTRUCTION.md) gives a fixed
error floor when the two absorbed spans have zero intersection, including
partial bands. Neither statement excludes a coupled three-band decoder
or the weaker complete-probe contract.

For that weaker contract, the separate [probe-domain theorem](WEAK_PROBE_PAULI_DOMAIN.md)
now gives independent four-wire cells an error floor of $`1/12`$ for
$`p\ge3`$ and $`1/4`$ for $`p\ge5`$, even with arbitrary joint outer
Clifford walls. The positive four-address cell saturates its exact domain
bound; the failure concerns its uncoupled extension, not this finite word.
