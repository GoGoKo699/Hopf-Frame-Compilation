# A paid field recoupler with a moving Pauli helper query

[Recursive tower and multiplication compiler](RECURSIVE_FIELD_SHEAR_DECODER.md) ·
[Complete probe readout](BANK_MIXING_AND_LOOKUP_PHASE.md#rank-one-indicator-equivalence) ·
[Finite projective allocation](PROJECTIVE_FIELD_RECOUPLER.md)

An additional, explicitly queried arbitrary helper removes norm division from
the recursive field recoupler. Two multiplication-XOR batches and fixed
Clifford walls give a complete native cell at every field degree, including
the zero sector. All cells at the same recursion level can share one moving
helper register. The complete probe echo returns that register and every
private multiplication mask.

The recursive construction pays sixteen T bands at each level. Section 7
gives a second decoder for the same augmented query, with logarithmic
depth in the address length, and identifies the complete dirty indicator
used in that construction. The earlier identity-helper decoder retains
its separate contract.

## 1. Complete augmented contract and resource statement

Use the tower, address coefficients and field shear from
[the recursive proof](RECURSIVE_FIELD_SHEAR_DECODER.md#1-a-compatible-tower-and-address-pencil).
For $`p\ge2`$, write $`m=2^p`$, $`k=m/2`$, and

```math
Q^{(p)}_x|U,V\rangle=|U,V+\beta_p(x)U\rangle,
\qquad U,V\in K_p,\qquad [K_p:\mathbb F_2]=k.
```

For every level $`2\le j\le p`$, allocate one arbitrary register
$`T_j\in K_{j-1}`$ of $`\ell_j=2^{j-2}`$ qubits. It is shared by
all cells at that level. Let $`x_{\lt j-1}`$ denote the lowest
$`j-1`$ address bits, and define the literal helper Pauli query

```math
P_p(x)=\prod_{j=2}^p Z_{T_j}^{\beta_{j-1}(x_{\lt j-1})},
\qquad
Z_{T_j}^{\lambda}|t\rangle
 =(-1)^{\mathrm{Tr}_{K_{j-1}/\mathbb F_2}(\lambda t)}|t\rangle.
                                                               \qquad(1)
```

There is an exact native decoder $`\widehat E_p`$, acting on the
principal bank and these helpers, such that

```math
\widehat E_p\bigl(Q^{(p)}_x\otimes P_p(x)\bigr)\widehat E_p^\dagger
   =W_x\otimes P_p(x)                                     \qquad(2)
```

for every address, with literal phase one. Here $`W_x`$ translates the
physical bank-wire indices by x. Any additional private masks used by the
native word have identity action in (2). All inputs, including the address,
moving helpers, masks and external references, may be arbitrary.

Put $`c_j=2^{p-j}`$ and

```math
e_j=\max\{1,\lceil\log_2(2\ell_j-1)\rceil\}=j-1,
\qquad s_j=(2\ell_j-1)e_j^2,
\qquad A_p=\sum_{j=2}^p c_js_j.
```

The forward and actual reversed native words obey

```math
T(\widehat E_p)\le32A_p=O(mp^3),\qquad
D_T(\widehat E_p)\le16(p-1),                              \qquad(3)
```
```math
h_{\rm moving}=m/2-1,\qquad
h_{\rm private}\le\max_{2\le j\le p}2c_js_j=O(mp^2),
\qquad G_{\rm all}=O(m^2p^2),\qquad w_{\rm clean}=0.       \qquad(4)
```

Private masks return after each complete multiplication-XOR batch. The
moving registers need not return under the decoder or the augmented query
alone. Their return in the complete indicator is proved in Section 5.

## 2. The complete cell without a nonzero promise

At one level let $`F=K_{j-1}`$, $`K=F[\theta]`$ with
$`\theta^2+\theta+\delta=0`$. Write

```math
U=d_0+d_1\theta,\qquad
q_0(U)=d_0d_1,\qquad
B(U,a)=d_0a_1+d_1a_0,\qquad a=a_0+a_1\theta.
```

For a binary field J write
$`\chi_J(z)=(-1)^{\mathrm{Tr}_{J/\mathbb F_2}(z)}`$.
Introduce one arbitrary $`T\in F`$. For a function f of the principal
coordinates, $`A_f`$ adds f into T and preserves those coordinates.
Use the fixed trace-pairing Fourier transform

```math
H_{K,\theta}|v\rangle
  =|K|^{-1/2}\sum_{a\in K}\chi_K(\theta av)|a\rangle.
                                                               \qquad(5)
```

This is a literal self-adjoint Clifford involution: multiplication by
$`\theta\ne0`$ preserves nondegeneracy of the symmetric trace pairing.
It acts on arbitrary inputs and requires no prepared Fourier state.

In Fourier coordinates for V, the field query is the diagonal operator
$`D_{\lambda,b}(U,V)=\chi_K((\lambda+b\theta)UV)`$, where
$`\lambda=\beta_{j-1}(x_{\lt j-1})`$ and $`b\in\mathbb F_2`$.
The chronological cell is

```math
A_{\mathrm{Tr}_{K/F}(UV)+q_0(U)},\quad
H_{K,\theta},\quad A_{B(U,a)}.
```

Its matrix $`C=A_BH_{K,\theta}A_{\mathrm{Tr}(UV)+q_0}`$
satisfies the complete identity

```math
C(D_{\lambda,b}\otimes Z_T^\lambda)C^\dagger
 =\chi_F\bigl(\lambda[B(U,a)+q_0(U)]\bigr)
      X_a^{bU}\otimes Z_T^\lambda.                         \qquad(6)
```

Conjugating the helper Pauli through the first addition contributes
$`\chi_F(\lambda[\mathrm{Tr}_{K/F}(UV)+q_0(U)])`$.
The relative-trace term cancels the lower-address phase of the source.
The remaining $`\chi_K(b\theta UV)`$ becomes the translation
$`X_a^{bU}`$ under (5). Conjugating the helper Pauli through the final
addition contributes $`\chi_F(\lambda B(U,a))`$. This addition
commutes with the translation because

```math
B(U,a+bU)+B(U,a)=bB(U,U)=0.
```

Thus there is no extra helper translation or scalar phase. For $`U=0`$,
both additions are identity and the cell is the specified Fourier transform
on the whole V fiber. Equation (6) still holds, including the helper Pauli.
No zero flag, inverse, division, postselection or discarded column is used.

Prepend the original $`H_K`$ wall on V. After the cell, replace
$`(U,a)`$ by $`(a,c)=(a,a+U)`$ and apply the original two
$`H_F`$ walls on the second coordinates of a and c. Indeed,

```math
q_0(a)+q_0(a+U)=B(U,a)+q_0(U),
```

so the phase in (6) becomes the two lower field queries. The translation
$`a\mapsto a+U`$ swaps a and c. The resulting recoupler obeys

```math
\widehat R_j(Q^{(j)}_{(x',b)}\otimes Z_{T_j}^{\beta_{j-1}(x')})
                  \widehat R_j^\dagger
 = (Q^{(j-1)}_{x'}\otimes Q^{(j-1)}_{x'})S^b
                    \otimes Z_{T_j}^{\beta_{j-1}(x')}.     \qquad(7)
```

The base decoder is the retained Clifford $`\widehat E_1=\mathrm{CX}(V;U)`$.
All primitive additions and Fourier transforms are self-inverse; the actual
inverse reverses them and the fixed coordinate walls in their literal order.

## 3. Sharing the same helper query through every child

Every recoupler uses its moving helper only through translations. It therefore
commutes with every X on that helper. In the common helper-X basis it is block
diagonal, with a Clifford on its own principal core in each sector: each
addition becomes a Boolean quadratic phase there.

Two decoders on disjoint principal cores that share the same moving helpers
commute. Their sector blocks act on disjoint cores. The core-only half swap
exchanges two identical child decoders, so their product commutes with that
swap. This reasoning applies to complete ideal blocks; it does not license
rearranging overlapping individual native gates.

In particular the shared Pauli query is not duplicated. If L and R are the
two child decoders and $`P=P_{j-1}(x')`$, their identities give the successive
conjugations

```math
Q_LQ_RP\ \xrightarrow{\ R\ }\ Q_LW_RP
        \ \xrightarrow{\ L\ }\ W_LW_RP.                   \qquad(8)
```

Each child commutes with the other core's query or translated target.
Equation (8) carries the same single P through both conjugations. The new
$`T_j`$ query commutes with the children, which use only lower-level helpers.
Together with (7) and swap invariance, this proves (2) inductively.
The number of persistent moving helper qubits is
$`\sum_{j=2}^p\ell_j=m/2-1`$.

## 4. Two simultaneous native batches per level

The first addition factors into two field products:

```math
\mathrm{Tr}_{K/F}(UV)+q_0(U)
   =d_0(v_1+d_1)+d_1(v_0+v_1).
```

The second is $`B(U,a)=d_0a_1+d_1a_0`$. The
[interpolation multiplication-XOR compiler](RECURSIVE_FIELD_SHEAR_DECODER.md#5-a-paid-multiplication-xor-component)
expresses one product over F using at most $`s_j`$ cubic products of
binary linear forms. Each cubic sign has an eight-phase Gray word on one
arbitrary private mask, returning it exactly. All private T targets are
distinct, so an entire batch occupies eight native T bands.

The first addition's multiplicands may share principal inputs. The parity-cube
identity remains valid for arbitrary Boolean linear forms, including dependent
or overlapping forms: evaluate the same scalar identity at their three bits.
Incorporate the sums $`v_1+d_1`$ and $`v_0+v_1`$ directly into
the mask parity updates. No clean multiplicand or product is loaded. Their
support is at most $`2\ell_j`$, which only changes the Clifford gate bound
by a constant factor.

At a fixed level all first additions from its $`c_j`$ cells combine
into one XOR into the shared $`T_j`$. Perform every disjoint core Fourier
wall, then combine all second additions. A cell's core Fourier commutes with
the other cells' additions, so this reordered aggregate equals the product
of the complete recouplers. The shared target's Hadamard walls are emitted
once around each aggregate phase batch.

Each batch has at most $`2c_js_j`$ private masks and $`16c_js_j`$ T gates.
Two batches give sixteen bands and $`32c_js_j`$ T gates per level.
The mask pool can be reused only after the preceding complete batch returns
it. Summing proves (3)--(4). A term's eight parity updates touch at most
$`2\ell_j`$ principal controls each. Their total Clifford cost is
$`O(\sum_j c_js_j\ell_j)=O(m^2p^2)`$. The fixed field-basis,
Fourier and coordinate walls cost $`O(\sum_jc_j\ell_j^2+m)`$ and
are included. Reversing the emitted native word preserves these budgets.

## 5. Complete readout, precision and the remaining depth obligation

The coherent helper query (1) is Clifford. Each coefficient is binary-linear
in the address, so its phase is a fixed bilinear form between the address
and helper bits, implementable with CZ gates. It costs no T gates and
$`O(mp)`$ elementary gates. It must be emitted together with the
[eight-band coherent field query](RECURSIVE_FIELD_SHEAR_DECODER.md#an-eight-band-query-for-any-binary-multiplication-matrices),
whose T count and returned-mask requirement are $`8kp`$ and $`kp`$.

Let $`S_x=\widehat E_p(Q^{(p)}_x\otimes P_p(x))\widehat E_p^\dagger`$.
Then

```math
S_x(Z_{Y_0}\otimes I)S_x^\dagger=Z_{Y_x}\otimes I.         \qquad(9)
```

Conjugating (9) by ordinary Hadamards on the bank gives the complete dirty
indicator $`|x,Y,w\rangle\mapsto|x,Y\oplus e_x,w\rangle`$.
The helper query cancels in this actual source/probe/adjoint echo. Every
moving helper, every private mask and every external reference is included;
the assertion is an operator identity, not a check on initialized work.

The direct indicator uses four decoder or actual-adjoint words and two field
queries. Reusing returned private masks between complete calls gives

```math
T(I_p)\le128A_p+16kp,\qquad
D_T(I_p)\le64(p-1)+16,                                  \qquad(10)
```
```math
h(I_p)\le m/2-1+\max\{kp,\max_j2c_js_j\},\qquad
G_{\rm all}(I_p)=O(m^2p^2),\qquad w_{\rm clean}=0.        \qquad(11)
```

The query's private-mask pool is included even when it exceeds the decoder's
pool. All identities and primitives are exact; no tolerance allocation or
precision-state preparation is hidden in these budgets.

The counts have the specific relaxed discovery form
$`T,h=m\mathrm{poly}(p)`$ and
$`G=m^2\mathrm{poly}(p)`$ used by
[polynomial-overhead fusion](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md).
That transfer reduces the source address by $`O(\log p)`$ bits; it does
not turn the present $`O(p)`$ source depth into the required
$`o(\log(p+2))`$. These resource statements make the uniform construction
explicit without supplying a new frame-depth theorem.

In helper-X coordinates the decoder sectors satisfy the Clifford cocycle
$`E_{r+\Lambda x}Q_x=W_xE_r`$. A generic formula using controlled wire
translations can hide the barrel-shift cost. The present word pays every
recursive stage explicitly; sharing helpers reduces width without removing
the principal-core dependency between levels. The [sector-degree theorem](FIELD_SECTOR_DEGREE_BARRIER.md)
requires helper-tableau degree at least $`p-1`$ and proves that this construction
attains it. A bounded number of helper-affine stages in one fixed helper frame
cannot implement the full decoder. A faster joint native word must leave that
stage grammar while preserving the complete return contract.

The source is $`Q_p\otimes P_p`$, not $`Q_p\otimes I`$. Its helper
Paulis can acquire signs under conjugation even when their binary labels
are fixed. Results requiring a literally fixed source Pauli space need their
hypotheses checked before application to this augmented family.
The [weak-probe domain theorem](WEAK_PROBE_PAULI_DOMAIN.md#field-query-and-pauli-helpers)
explicitly includes these Pauli helper phases and keeps the label/sign distinction.

## 6. Exact finite evidence

`tests/test_augmented_field_recoupler.py` verifies the two-product identity,
the full zero fiber, and a native eight-band aggregate with overlapping forms,
a shared target and arbitrary returned private masks. A complete four-address
native decoder is checked on every one of its seven-wire columns, including
all four private-mask basis values and the actual native inverse.

The independent computational-helper-basis evaluator checks every column of
the algebraic decoder for two and three address bits, all addresses and all
reversed algebraic-word columns. This includes all 2,080 principal/helper
columns and 16,512 address-column comparisons. A negative calibration shows
that omitting the moving Pauli query changes the action; the complete probe
echo cancels it exactly. The finite fixture does not replace the all-size
native composition proof above or claim an improved general depth frontier.

## 7. A logarithmic-depth decoder for the full augmented query

Let $`Q_x=Q^{(p)}_x`$ and split the address as $`x=h+ub`$, where
$`h\in H=\mathbb F_2^{p-1}`$, $`b`$ is the top-bit direction and $`u`$ is
one bit. In the helper-X basis the Pauli (1) shifts the helper label
$`r`$ by $`\Lambda x`$. Its top register sees the $`p-1`$ independent
coefficients of $`\beta_{p-1}(h)`$ through a nondegenerate trace pairing.
The remaining registers also ignore $`b`$, so $`\ker\Lambda=\langle b\rangle`$.
Select $`p-1`$ independent output rows and invert their
binary map on $`H`$ to obtain a linear clock $`s(r)\in H`$ with
$`s(r+\Lambda x)=s(r)+h`$. This is an invertible Clifford change on
the selected helper coordinates. Every helper retains its actual
translation in (2). A variant retaining only these selected moving
coordinates also satisfies its corresponding augmented identity.

Write $`k=2^{p-1}`$. The fixed binary-linear Clifford
$`C(U,V)=(V,V+\beta_p(b)U)`$ has inverse
$`C^{-1}(A,B)=(\beta_p(b)^{-1}(A+B),A)`$. Thus $`CQ_bC^\dagger=W_b`$,
the swap of the two physical half banks. In helper sector $`r`$ use
the decoder $`E_r=W_{s(r)} C Q_{s(r)}`$. Since the field shears compose
additively, its sector from $`r`$ to $`r+\Lambda x`$ obeys

```math
E_{r+\Lambda x}Q_x E_r^\dagger
 = W_{s+h} C Q_{s+h+x+s} C^\dagger W_s
 = W_{s+h} W_{ub} W_s = W_x .                             \qquad(12)
```

Equation (12) proves (2) on every coherent input, including every value
of both field operands and all helper labels. The helper Fourier and
coordinate changes are undone at the boundaries.

### Paying for the selected translations

Let $`P_s`$ translate the indices of one $`k`$-wire bank by $`s`$.
For two banks $`A,B`$ and a separate arbitrary program bank $`Y`$ use
the product-XOR $`B\gets B+A*Y`$, where
$`(A*Y)_v=\sum_j A_{v+j}Y_j`$. Repeating the product with
$`Y+e_s`$ in place of $`Y`$ adds just $`P_s A`$. An indicator and its
actual inverse toggle and return $`Y`$, so one selected shear costs two
indicator calls. Conversely, conjugating $`X_{A_0}`$ by the selected
shear gives $`X_{A_0}X_{B_s}`$. A second fixed $`X_{A_0}`$ leaves the
complete indicator on $`B`$. These two-call reductions have error at
most $`2\delta`$ at full-operator component error $`\delta`$.

In the group algebra, the substitution $`g_i=1+u_i`$ with $`u_i^2=0`$
converts XOR convolution to disjoint-subset multiplication. This
binary zeta change and its inverse use $`O(kp)`$ CNOTs. The
[ranked-product identity](ALGEBRA_QUERY_EXTRACTION.md#4-fast-subset-convolution-as-a-native-returned-query)
writes its output Fourier sign as $`R=k\binom{p+1}{2}`$ cubics
$`\alpha\beta\gamma`$ in binary linear forms. For each cubic use one
arbitrary mask $`d`$. Walk through the eight vertices $`v`$ of a binary
three-cube in cyclic Gray order, pulsing the current
$`d\oplus v\cdot(\alpha,\beta,\gamma)`$ by T or T-dagger with
exponent $`(-1)^{|v|}`$. The integer exponent sum is
$`-4(-1)^d\alpha\beta\gamma`$ modulo eight, hence the sign is
$`(-1)^{\alpha\beta\gamma}`$ for both initial masks. Each edge changes
one entire bank of masks by its corresponding linear forms. The
[returned linear-shear compiler](ALGEBRA_QUERY_EXTRACTION.md#3-linear-preprocessing-without-zero-buffers)
does this in $`O(kp^2)`$ Cliffords and reusable dirty wires per edge.
All row masks pulse in parallel. Including the output Hadamards, one
complete product costs $`8R`$ T gates, eight T layers, and $`O(kp^2)`$
elementary Cliffords and private dirty qubits. The reversed word has
the same costs.

Apply three alternating selected shears into $`B,A,B`$, then swap the
halves. Since $`P_s^2=I`$, this returns $`(P_s A,P_s B)=W_s(A,B)`$.
All three use the same $`Y`$. Apply the successive product pairs at
$`(Y,Y+e_s)`$, $`(Y+e_s,Y)`$ and $`(Y,Y+e_s)`$ respectively, and finally
untoggle $`Y`$. This uses six products and four indicator or inverse
calls. The [coherent field shear](RECURSIVE_FIELD_SHEAR_DECODER.md#an-eight-band-query-for-any-binary-multiplication-matrices)
for $`Q_s`$ adds eight bands and $`O(kp)`$ T gates and masks. For
$`0\le i\lt p-1`$, its coefficient $`\alpha_i`$ lies in a binary subfield
of dimension $`2^i`$. In the tower basis its fixed multiplication matrix
has $`k/2^i`$ blocks of size $`2^i`$ and therefore at most $`k2^i`$ entries
to emit across those blocks. Their sum is at most $`k(k-1)`$. Thus the
field shear and $`C`$ together cost $`O(k^2)`$ elementary Cliffords.

With the [linear-size exact dirty indicator](NONUNIFORM_DIRTY_INDICATOR.md)
of depth $`d_I(t)=O(\log(t+2))`$ on $`t`$ address bits, the forward and
actual reversed decoder have simultaneous bounds

```math
D_T(E)\le 4d_I(p-1)+56,\qquad
T(E),\ h_{\rm private}(E)=O(2^p p^2),\qquad
G(E)=O(4^p),\qquad w_{\rm clean}=0.                     \qquad(13)
```

The program bank, returned product and query pools, and coordinate
changes are included. Substituting four full-operator $`\delta`$-accurate
indicator calls gives decoder error at most $`4\delta`$, conjugated
source error at most $`8\delta`$, and complete probe error at most
$`16\delta`$ via two actual source visits. Substituted circuits incur
their own forward/inverse resources each time. The logarithmic bound
in (13) prices an existing indicator and leaves the general R96
frame-depth bound in place.

### An alternative using the balanced algebra permutation

One can also pay for $`W_s`$ via the [HRT conjugation](ALGEBRA_QUERY_EXTRACTION.md#1-an-arbitrary-table-scalar-on-the-entire-coefficient-space).
Add one square-zero generator $`u_*`$ to the $`p-1`$ generators above.
An arbitrary element of the augmentation ideal is uniquely $`z+a u_*`$,
where $`z`$ has no constant coefficient and $`a`$ is an arbitrary half bank.
The absent constant coefficient is structural. The conjugation
$`J_s=F_p X^{\iota(s)}F_p^\dagger`$ of that reference adds $`A_s+yA_s`$
to the occupied coefficient word $`y`$, with
$`1+A_s=\prod_i(1+s_i u_i)=g_s`$. Consequently it sends
$`(z,a)`$ to $`(z+A_s+zA_s,ag_s)`$. Its first-register action is an
involution since $`A_s^2=0`$. Use the same arbitrary $`z`$ with the first
half bank and then the second. The two visits return $`z`$ and multiply
both halves by $`g_s`$, which becomes $`P_s`$ after the zeta basis change.
The [paid balanced HRT word](ALGEBRA_QUERY_EXTRACTION.md#5-a-smaller-count-complete-hrt-source)
gives four $`F_p`$ or actual-inverse calls. Including $`Q_s`$ and $`C`$,
this decoder has at most $`800\lceil\log_2 p\rceil+8`$ T layers,
$`O(2^p p^2\log(p+2))`$ T gates, $`O(2^p p^2)`$ private dirty width
and $`O(4^p)`$ elementary Cliffords. Its balanced product levels retain
logarithmic dependence on $`p`$.

The tests in `tests/test_augmented_field_recoupler.py` check the
four-toggle schedule and the HRT forward/inverse coefficient maps on
occupied halves, and (12) on all augmented Fourier-sector columns
for $`p=2,3`$.
