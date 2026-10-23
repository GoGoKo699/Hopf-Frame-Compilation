# Entangling selector controls and complete acceptance screens

[Working checkpoint](../../WORKSPACE.md) ·
[Co-moving stages](CO_MOVING_LOOKUP_STAGE.md) ·
[Earlier local-selector obstruction](CO_MOVING_SELECTOR_READOUT.md) ·
[Complete recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)

Fixed physical-row entangling selector operations have cheap, exact native
controls. This removes a synthesis uncertainty, but not the acceptance
problem. A common two-dimensional selector module excludes arbitrary words
using linear selector permutations and global Hadamards. Partial Fourier
mixing leaves that module, yet the specified four-path source still has a
two-bit Gram floor of one sixth. Adding an affine selector flip gives an
exact finite algebraic separator and an all-size separation at address zero;
neither supplies a uniform gap or a lookup construction. A square-stage commutator gives a complete two-bit
calibration and identifies precisely why its affine-Pauli extension stops.

None of these statements changes R96 or proves a general depth lower bound.

## 1. Registers and the unpaid acceptance identity

There are p address bits, $`m=2^p`$ arbitrary bank bits Y and m arbitrary
selector bits s. Retain the complete co-moving stage

```math
R_C=i^{n_C}L_C^h,\qquad
h=\sum_uY_us_u\pmod2,\qquad
n_C=\sum_uY_us_u\sum_{c\in S_C}(x_c\oplus u_c).
```

C is a disjoint CNOT matching, $`S_C`$ its control coordinates, and L
relabels address, bank and selectors together. Its actual inverse and
returned implementation helpers are paid in the linked stage proof.
All states below may be coherent and entangled with references. Define

```math
\Omega|x,Y,s\rangle=(-1)^{Y_x}|x,Y,s\rangle,\qquad
\Pi_\pm=(I\pm\Omega)/2,\qquad H_Y=\Pi_+/2.
```

The desired bounded-clean accepted block K must satisfy
$`\|K^\dagger K-H_Y\|_{\rm op}\le\delta`$ on the complete space.
A useful first requirement is
$`K\Pi_-=0`$; positive-sector coverage and normalized strength remain
separate obligations. All impossibility statements specify their actual
operation grammar. Small-size certificates do not exclude a family patched
at those sizes.

## 2. A controlled entangling Clifford with four native T layers

For each physical row-label coordinate j, let

```math
E_j=\prod_{u:u_j=0}\mathrm{CX}(s_u;s_{u\oplus e_j}),\qquad
M_p=\prod_{j=0}^{p-1}E_j.
```

On binary selector columns these are tensor factors of
$`\begin{pmatrix}1&0\\1&1\end{pmatrix}`$. The factors commute,
$`M_p^2=I`$, and $`(M_ps)_v=\bigoplus_{u\subseteq v}s_u`$.
The direct word has pm/2 CNOTs, zero T gates and no helpers.

Put q=m/2 and split the rows by the highest label bit. With
$`A=M_{p-1}`$ (A=I for p=1), use binary matrices

```math
M_p=\begin{pmatrix}A&0\\A&A\end{pmatrix},\quad
J=\begin{pmatrix}I&0\\I&I\end{pmatrix},\quad
Q=\begin{pmatrix}I&I\\0&I\end{pmatrix},\quad
P=Q\begin{pmatrix}I&0\\0&A\end{pmatrix}Q
 =\begin{pmatrix}I&A+I\\0&A\end{pmatrix}.
```

Direct multiplication gives $`P^2=I`$ and $`M_p=PJP`$.
The explicit Q,A,Q chronology for P uses
$`G_P=m+(p-1)m/4`$ CNOTs. No generic change-of-basis cost is hidden.

To control M by an arbitrary signal a, emit P, the batch
$`\prod_{i=0}^{q-1}\mathrm{CCX}(a,s_i;s_{q+i})`$, and the actual
P-adjoint. The [shared-control bilinear batch](POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md#2-pay-for-a-root-controlled-complete-bank-transition)
has common forward/actual-inverse bounds
$`T\le6q+(q\bmod2)`$, $`G\le12q+3`$ and depth four, without helpers.
The common control receives no T gate for even q; common S powers are
reduced modulo four in both orientations. Thus

```math
T_{\mathrm cM}\le3m+(q\bmod2),\qquad
G_{\mathrm cM}\le(p+15)m/2+3,\qquad D_T(\mathrm cM)\le4.       \qquad(1)
```

There is zero clean work, helper work or error. The control returns;
selectors undergo exactly the requested controlled permutation. Literal
phases are +1. The actual native reverse has the same majorants.

For $`H_s=H^{\otimes m}`$, Fourier summation gives
$`H_sU_AH_s=U_{A^{-\mathsf T}}`$. Therefore
$`N_p=H_sM_pH_s=U_{M_p^{\mathsf T}}`$ is another CNOT permutation,
not a nonmonomial mixer. Transposing the explicit CNOT construction gives
the same control ledger (1), without adding Hadamard walls.

In contrast, let $`H_j=\prod_{u:u_j=0}H_{s_u}`$ and
$`C_j=H_jM_pH_j`$. Fixed coordinate relabeling reduces it to the same
block calculation. Its controlled implementation is $`H_j,\mathrm cM,H_j`$,
so

```math
T_{\mathrm cC_j}\le3m+(q\bmod2),\qquad
G_{\mathrm cC_j}\le(p+17)m/2+3,\qquad D_T(\mathrm cC_j)\le4. \qquad(2)
```

For the highest-coordinate split, write input selectors as (a,b) and
outputs as (v,w), with q bits in each half. The complete kernel is

```math
\langle v,w|C_j|a,b\rangle
 =2^{-q}(-1)^{a\cdot Aw+a\cdot b+v\cdot w+v\cdot Ab}.       \qquad(3)
```

Indeed the first H wall introduces t, M sends (t,b) to (At,A(t+b)),
and fixing w forces $`t=Aw+b`$. The two Fourier phases give (3).
This is a dense real Clifford involution. In particular,

```math
|g_j\rangle:=C_j|0^m\rangle
 =2^{-m/2}\sum_s(-1)^{q_j(s)}|s\rangle,\qquad
q_j(s)=\sum_{u:u_j=0}s_us_{u\oplus e_j}\pmod2.              \qquad(4)
```

It is not a supplied state: (4) describes an allowed input's image under
the paid unitary. Arbitrary selector inputs remain in the contract.

More generally, a binary linear involution B has $`(B+I)^2=0`$.
Put $`r=\mathrm{rank}(B+I)`$, choose paired vectors with $`(B+I)v_i=u_i`$, and
complete the u vectors to a basis of the kernel. This gives r disjoint
CNOT blocks and fixed coordinates. A literal CNOT basis change of cost
$`L_P`$ therefore controls B with $`T\le6r+(r\bmod2)`$,
$`G\le2L_P+12r+3`$ and depth four. Gaussian elimination gives the
safe bound $`L_P\le d(d-1)+3(d-1)`$ on d wires. This is not a
theorem about arbitrary controlled Cliffords or free phase corrections.

## 3. An all-size common module, including global Hadamards

At the marked address/bank input $`x=0,Y=e_0`$, every linear matching
fixes x and Y, its phase exponent vanishes, and

```math
R_C|0,e_0,s\rangle=|0,e_0,C^{s_0}s\rangle.                 \qquad(5)
```

This is a bijective phase-free selector permutation: C fixes row zero,
so $`s_0`$ is unchanged. Consequently it fixes both $`|0^m\rangle`$
and $`|+^m\rangle`$. At the unmarked input $`x=0,Y=0`$, every
stage is identity on the entire selector space.

Let $`\mathcal S=\mathrm{span}\{|0^m\rangle,|+^m\rangle\}`$.
Every selector CNOT network fixes this space pointwise, and global H
preserves it by swapping its two displayed, nonorthogonal vectors.
Thus every word in linear co-moving stages, fixed selector CNOT networks,
global H, their actual adjoints and signal-controlled versions acts
identically on this module at the two bank inputs. Arbitrary signal-only
interlayers, initialized signals and common acceptance rules retain that
equivalence. Tensoring with independent returned implementation helpers
and references does not change it.

The target Grams on these copies of the module are zero and I/2.
If their common accepted Gram is G, the triangle inequality gives

```math
\max\{\|G\|,\|G-I/2\|\}\ge1/4.                           \qquad(6)
```

This is independent of word length, number of matchings and signal-angle
precision. It includes the zeta and transposed-zeta proposal and its
global-H modification. It excludes neither partial H walls, affine
selector flips, selector phases, address-dependent interlayers nor
unrestricted native circuits. The proof uses the empty and singleton
banks; it does not assume that bank-weight-dependent signal operations
are present.

## 4. Breaking that module still does not give four-path acceptance

The graph state (4) supplies an analytic direction check. Under (5), its
amplitude changes by

```math
(-1)^{s_0[q_j(Cs)+q_j(s)]}.                                 \qquad(7)
```

If $`Ce_j=e_j`$, the matching of rows separated by $`e_j`$ is
preserved, so the state is fixed. Choosing a CNOT's target direction
therefore does not remove this particular witness. Choosing its control
direction changes the state; it is the relevant next test, not a source
proof. The same graph state is fixed by global H, by the two-vertex
graph-state identity on each disjoint pair.

For p=2 fix $`C=(0\to1)`$, $`D=(1\to0)`$ and the moved-direction
choice

```math
A_0=R_C,\quad A_1=C_0R_CC_0,\qquad
B_0=R_D\ \hbox{or}\ H_sR_DH_s,\quad B_1=C_1B_0C_1.
```

A one-signal source with three signal interlayers and the two diagonal
branch multiplexors has accepted block

```math
K=\sum_{a,b=0}^1c_{ba}B_bA_a,\qquad
c_{ba}=\langle0|S_2|b\rangle\langle b|S_1|a\rangle
        \langle a|S_0|0\rangle.                             \qquad(8)
```

Each multiplexor is implemented by controlled selector conjugation around
one uncontrolled stage; no controlled R is assumed. Exact fixed signal
interlayers give depth at most $`2(116)+4(4)=248`$ plus their own
T-depth. Approximate signal synthesis must be charged separately.

The following stronger test allows *arbitrary complex* coefficients c,
not merely those realizable in (8). Restrict all four branches to the
marked input $`x=0,Y=e_0`$, retaining the whole 16-dimensional selector
space. Order them as $`(B_0A_0,B_0A_1,B_1A_0,B_1A_1)`$ and put
$`G_{ij}=\mathrm{tr}(W_i^\dagger W_j)`$. Direct exact multiplication
gives, respectively,

```math
G=\begin{pmatrix}
16&9&9&7\\9&16&6&9\\9&6&16&9\\7&9&9&16
\end{pmatrix},\qquad
G=\begin{pmatrix}
16&9&10&13/2\\9&16&6&10\\10&6&16&9\\13/2&10&9&16
\end{pmatrix}.                                               \qquad(9)
```

In both cases $`G-8\mathbf1\mathbf1^{\mathsf T}`$ is Hermitian,
has diagonal eight and off-diagonal absolute row sum at most five.
It is positive definite. Hence the average marked acceptance is at least
$`|\sum c_{ba}|^2/2`$. On the empty bank every branch is identity,
so its acceptance is $`|\sum c_{ba}|^2`$. Complete Gram error epsilon
would imply

```math
\epsilon\ge\tfrac12|\sum c_{ba}|^2
 \ge\tfrac12(\tfrac12-\epsilon),\qquad \epsilon\ge1/6.       \qquad(10)
```

This is a complete two-bit counterexample, including coherent selectors,
for both displayed four-path sources and all their scalar coefficients.
It is not an all-p lower bound or an exclusion of other accepted words.
Removing the previous common invariant was necessary but insufficient.

## 5. Separation after an affine selector flip

### Complete two-bit certificate

An independent exact check considers the enlarged p=2 alphabet

```math
\{R_C,R_D,M_2,N_2,H_s,F\},\qquad F=X^{\otimes4},\qquad
S=R_C+R_D+M_2+N_2+H_s+F.                                   \qquad(11)
```

F takes $`|0^m\rangle`$ outside the module used in (6). More strongly, the complete positive and
negative charge representations of this alphabet have no nonzero common
intertwiner at p=2. This finite statement retains all 1024 address, bank
and selector basis states, not a singleton promise.

The reproducible certificate partitions (x,Y) into its simultaneous
C,D relabeling orbits, keeping all 16 selector columns within each block.
There are 20 such blocks and ten in each charge. Every generator preserves
each block. Compute the characteristic polynomial of S in each block over
$`\mathbb F_{65537}`$, using $`i\mapsto256`$ and the exact rational
Hadamard entries. All 100 opposite-charge polynomial pairs have greatest
common divisor one. The native stage matrices are independently compared
with the existing literal complete-stage fixture.

Why this proves the finite claim: these characteristic polynomials are
monic over $`\mathbb Z[i,1/2]`$. A coprime reduction modulo this prime
makes their resultant nonzero, so the original opposite-charge spectra
are disjoint over the complex numbers. The Sylvester equation
$`S_-T=TS_+`$ therefore has only T=0, even when S is not normal. An
intertwiner of every generator would also intertwine their sum, and
decomposes into these block pairs, so it must vanish.

The test uses exact modular arithmetic, not numerical eigenvalue
separation. Its characteristic-polynomial calculation was independently
audited by Faddeev--LeVerrier and Cayley--Hamilton checks. It does not
supply a normalized polynomial, a physical spectral filter, a stable
gap, or any all-size bound. In particular S is non-Hermitian; disjoint
finite spectra do not certify a uniform singular-value gap or a low-depth
implementation of its charge projectors.

### All-size separation within address zero

There is also an analytic result for every p at least two, including every
bank background. Restrict both charge representations to the invariant
address-zero space, so their charges are determined by $`Y_0`$. Use the
alphabet $`\{R_C:C\in\mathcal C\}\cup\{M_p,N_p,H_s,F\}`$, now with
$`F=X^{\otimes m}`$. Suppose the group generated by the linear matchings
has at least three distinct images of the top row $`t=(1,\ldots,1)`$.
The two opposite CNOTs already satisfy this condition for every p at least
two. Under these hypotheses the two address-zero charge representations
have no nonzero intertwiner.

First the binary zeta identities give

```math
M_pFM_p=X_{s_0},\qquad N_pFN_p=X_{s_t},\qquad
H_sX_{s_t}H_s=Z_{s_t}.                                    \qquad(11a)
```

Indeed M sends the all-ones selector column to the row-zero unit vector,
while its transpose sends it to the top-row unit vector. At x=0 the
row-zero contribution to n vanishes, every matching fixes row zero, and
flipping $`s_0`$ changes h by exactly $`Y_0`$. It leaves n unchanged and
commutes with the joint relabeling. Thus, with the actual stage adjoint,

```math
E_C:=R_CX_{s_0}R_C^\dagger X_{s_0}
       =L_C^{Y_0}\qquad\hbox{on the whole }x=0\hbox{ space}. \qquad(11b)
```

The phases cancel literally, including on nonempty backgrounds. Here L
still relabels the bank and selectors together; it is not a selector-only
substitution. Equation (11b) is identity on the positive sector and L on
the negative sector.

Let T intertwine every generator from the positive address-zero sector
to the negative address-zero sector. Since the generators are unitary,
T also intertwines their adjoints and all words. Its image W is therefore
fixed pointwise by every joint $`L_C`$, and hence by every generated
$`L_g`$. It is invariant under the selector Paulis $`X_{s_t},Z_{s_t}`$
from (11a). Choose g,h with $`t,g(t),h(t)`$ all distinct. For w in W,
both w and $`X_{s_t}w`$ are fixed by $`L_g`$, which gives
$`X_{s_{g(t)}}w=X_{s_t}w`$. The analogous Z relation uses h. Hence

```math
X_{s_t}X_{s_{g(t)}}w=w,\qquad
Z_{s_t}Z_{s_{h(t)}}w=w.                                    \qquad(11c)
```

These two Pauli operators anticommute: their row supports meet only at t.
They have no common +1 vector, so W=0 and T=0. Taking adjoints excludes
the reverse intertwiner as well. This proves separation for arbitrary
banks and coherent selectors within x=0, not merely for the earlier
empty/singleton witness.

### Scope of the remaining all-address question

The address-zero result does not exclude intertwiners involving nonzero
address orbits, even when one endpoint has address zero. It supplies no
uniform gap, normalized accepted block or complete recovery construction.
If only the same two address coordinates are used at p greater than two,
the entire alphabet commutes with every unused address X. With one bank
marking x but not its spectator-bit flip, this again gives a one-quarter
Gram floor. A scalable proposal must expose every address coordinate,
for example through the already paid full matching family, and prove its
acceptance anew. The finite separator is a seed, not that construction.

## 6. Complete two-bit readout and the boundary of square-stage words

There is a useful positive calibration. For a single-CNOT matching with
control i, let $`Q=\prod_u\mathrm{CZ}(Y_u,s_u)`$ and
$`A_i=QR_C^2`$. The invariance of n and h gives

```math
R_C^2=(-1)^{n_C},\qquad
A_i=Z_s^{a_i},\qquad (a_i)_u=Y_u(1\oplus x_i\oplus u_i).
```

Put $`B_j=H_sA_jH_s=X_s^{a_j}`$. The complete commutator is

```math
A_iB_jA_iB_j
 =(-1)^{\sum_uY_u(1\oplus x_i\oplus u_i)
                         (1\oplus x_j\oplus u_j)}I_s.       \qquad(12)
```

For p=2 and distinct i,j, (12) is exactly Omega, returning every selector.
For larger p it reads a codimension-two subcube parity, not a single row.

One need not emit eight expensive co-moving stages. Directly implement
$`A_i`$ by the shared-control batch $`\prod_u\mathrm{CCZ}(x_i,Y_u,s_u)`$
and CZ on the m/2 rows with $`u_i=0`$. The former has T-count 6m,
four T layers and Clifford count at most 10m+3; the row mask adds 3m/2
Cliffords. Hadamard conjugation adds 2m for B. Thus the word (12) and its
actual reverse have

```math
T\le24m,\qquad G\le50m+12,\qquad D_T\le16,                 \qquad(13)
```

with no initialized or private helpers and m arbitrary returned selectors.
For p=2, two bank-Hadamard walls give the complete indicator with
$`G\le52m+12`$ and the other bounds unchanged. This is a small-size
calibration, not an improvement to the growing-address frontier.

There is a precise limit to extending it with fixed selector Cliffords.
Fix a computational bank Y. All operations in this grammar preserve the
address and bank fibers. Consider any word of canonical work Paulis
$`X^{a(x)}Z^{b(x)}`$ with binary-affine labels, and fixed work Cliffords
acting on selectors and additional work, independent of x; no
address-dependent scalar phases are supplied.
Every such word has normal form

```math
U_Y(x)=\lambda i^{q_Y(x)}X^{a_Y(x)}Z^{b_Y(x)}C,\qquad
q_Y(x)=q_0+\sum_iq_ix_i+2\sum_{i\lt j}q_{ij}x_ix_j\pmod4.  \qquad(14)
```

Pauli multiplication adds twice a binary bilinear form. H adds
$`2a_jb_j`$ and exchanges labels; S adds $`a_j`$ and updates b;
CNOT acts linearly on the labels. The integer lift of a binary affine
form modulo four has only linear and even quadratic terms. These facts
prove (14) by induction, with literal constant Clifford phases in lambda.

If every work input returns with a scalar sign, C is a Pauli up to scalar,
the residual labels are constant, and q has constant parity. The returned
Boolean sign consequently has degree at most two. A singleton bank's
lookup sign has address degree p. Arbitrarily long words in this grammar,
including fixed entangling zeta interlayers, therefore cannot give exact
lookup at p at least three. Unsquared R stages, address motion, non-Clifford
selector operations and accepted signal blocks are outside this statement.

The same grammar has complete lookup error at least $`1/\sqrt2`$ for
p at least three. Otherwise, for every Hermitian work Pauli P,
$`\|UPU^\dagger-P\|\lt\sqrt2`$. Clifford Pauli discreteness forces
equality for every P, hence U is scalar. Two parities of q would give
scalar ratios $`\pm i`$, whose best simultaneous distance to the desired
signs is $`\sqrt{2-\sqrt2}>1/\sqrt2`$. Constant parity reduces to the
quadratic sign just proved; a wrong sign cannot be hidden by one common
scalar phase. This is a grammar-specific error floor, not a depth bound.

Literal phases matter. A freely assigned Hermitian normalization
$`P(a,b)=i^{ab}X^aZ^b`$ is not a free member of the above alphabet.
Indeed $`X^aZ^bP(a,b)=(-i)^{ab}`$, and its finite difference in
$`b=x_2\oplus x_3`$ gives $`(-1)^{x_1x_2x_3}`$. Assigning that phase
without synthesis would invalidate the quadratic closure.

## 7. Consequence for the next constructive question

The native control of the entangling selector involutions is paid. The
complete accepted identity remains unpaid. The CNOT/global-H grammar
cannot be repaired by longer signal processing; the displayed partial-H
four-path source is also rejected. The affine-flip enlarged algebra escapes
the common-intertwiner obstruction at two bits and, for every p at least two, between
the two address-zero sectors. Other address orbits and a uniform source or
gap remain open. The square-stage calibration cannot be extrapolated through
fixed affine-Pauli words.

The next admissible attack must specify a genuinely different accepted
operator with complete marked cancellation and whole positive-sector
coverage, using all address coordinates. A uniformly gapped intermediate
Gram would be useful only after its normalization and complete native
recovery are paid. Neither finite algebraic distinguishability nor removal
of a known witness supplies this. Preserve arbitrary banks, selectors,
references, actual inverses and off-diagonal terms throughout.

At every fixed requested tolerance, a qualifying sublogarithmic complete
source with fixed polynomial overhead would pass through the retained
fusion and frame-transfer interfaces. Even a constant-depth source leaves
their logarithmic-logarithmic contribution; it would not by itself prove
optimal constant-depth frames.

Exact fixtures: [native selector controls](../../tests/test_entangling_selector_readout.py),
[marked acceptance](../../tests/test_entangling_marked_acceptance.py),
[finite sector separator](../../tests/test_co_moving_sector_separator.py), and
[square-stage calibration](../../tests/test_square_stage_readout.py).
Finite fixtures corroborate the stated identities and finite certificates;
they do not replace the all-size arguments or establish an asymptotic source.
