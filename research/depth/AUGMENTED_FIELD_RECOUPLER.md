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

This constructs the augmented query equivalence stated below. Its current
compiler pays sixteen T bands at each recursive level. The resulting depth
is linear in the number of address bits and does not improve R96 or meet the
sublogarithmic-depth source target. The earlier identity-helper decoder and
its unpaid uniform coordinate maps retain their separate scope.

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
