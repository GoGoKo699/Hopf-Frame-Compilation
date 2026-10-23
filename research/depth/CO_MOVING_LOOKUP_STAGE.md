# Co-moving lookup stages and the remaining acceptance problem

[Working checkpoint](../../WORKSPACE.md) · [Quarter-phase sources](COUPLED_QUARTER_PHASE_SOURCE.md) ·
[Matching stages](QUARTER_PHASE_MATCHING_LIMITS.md)

A simultaneous address and bank permutation can preserve the complete lookup
charge even while the bank moves. The quarter-phase stage below has an exact
constant-depth native implementation on arbitrary inputs. Its simplest
accepted difference still fails, and scalar processing of that difference
cannot repair it. This is a paid component, not a complete shallow source or
an improvement to the general frame frontier.

## 1. A complete charge-preserving stage

Let $`p\ge2`$, $`m=2^p`$, and let C be a nonempty matching of q disjoint
CNOTs on the address, with control set S. All address bits x, bank bits Y
and selector bits s are arbitrary. Define the known Clifford involution

```math
L_C|x,Y,s\rangle=|Cx,CY,Cs\rangle,\qquad
(CY)_{Cu}=Y_u,\quad(Cs)_{Cu}=s_u.
```

The permutations of physical bank and selector wires are specified by C;
they do not require an address-selected row operation. Put

```math
\begin{aligned}
d_C(x,u)&=\sum_{c\in S}(x_c\oplus u_c),&
h&=\sum_u s_uY_u\pmod2,\\
n_C&=\sum_u s_uY_ud_C(x,u),&
N_C&=\sum_uY_ud_C(x,u),& w&=\sum_uY_u.
\end{aligned}                                                   \qquad(1)
```

Here $`d_C,n_C,N_C,w`$ are ordinary integers. Every matching control
coordinate is unchanged by C. Simultaneous relabeling therefore preserves
h and all these integer quantities. In particular the diagonal phase
$`D_C=i^{n_C}`$ commutes with $`L_C`$. The proposed stage is

```math
\widetilde R_C=D_C L_C^h.                                      \qquad(2)
```

The conditional permutation is well defined and involutory because h is
invariant under $`L_C`$. Thus (2) is unitary on the complete space. For

```math
\Omega|x,Y,s\rangle=(-1)^{Y_x}|x,Y,s\rangle,
\qquad H_Y=(I+\Omega)/4,
```

both factors in (2) commute with Omega: simultaneous relabeling preserves
$`Y_x`$, and the phase is diagonal. This is a completed-operator identity;
no native-boundary invariance is assumed. The address-dependent phase also
leaves the premises of the earlier
[uniform-address joint-relabeling screen](EXPANDER_ZERO_REFLECTION.md#9-the-uniform-address-obstruction-to-joint-relabeling).

## 2. Native construction, actual inverse and resource ledger

The [two-layer quadratic query](QUADRATIC_NATIVE_QUERY.md) computes h into
an arbitrary dirty bit b. Its polar rank on the 2m bank/selector variables
is $`r=2m`$. Call this full-space XOR A. It has

```math
T_A=4m+4,\qquad d_A=2,\qquad
G_A\le2L_A+18m+13,\qquad
L_A\le(2m)(2m-1)+3(2m-1).                                  \qquad(3)
```

The common Clifford bound includes the possible residual linear terms and
the positive-S expansion in the actual inverse. A needs one returned dirty
helper in addition to its output b.

Let B apply $`L_C`$ controlled by b. The matching has $`m/2^q`$ fixed
row labels, so the two row permutations together have
$`m-m/2^q`$ disjoint transpositions. Including the q controlled address
CNOTs, B consists of

```math
\ell_C=q+m-m/2^q                                            \qquad(4)
```

controlled pair operations with disjoint target pairs and common control b.
Each controlled swap uses two CNOTs and one Toffoli. To keep physical
T gates off the shared control b, implement a Toffoli on (b,u;v) with
a private arbitrary dirty bit d, chronologically, as

```text
H(v), CX(b,d), CCZ(d,u,v), CX(b,d), CCZ(d,u,v), H(v).
```

The middle phase is
$`(-1)^{(d\oplus b)uv+duv}=(-1)^{buv}`$, and d returns exactly.
Two ordinary seven-T, six-CNOT CCZ decompositions give 14 T gates,
16 Cliffords and at most 14 T layers. All (d,u,v) triples are disjoint,
so the pair operations run in parallel despite sharing b. Consequently

```math
T_B\le14\ell_C,\qquad G_B\le18\ell_C,\qquad d_B\le14.        \qquad(5)
```

Emit A, B, the actual A-adjoint, then B. Since B preserves h even while
moving Y and s, the two b values controlling B differ by h; their product
is exactly $`L_C^h`$. The bit b and every private helper return on the
complete space. Finally emit $`D_C`$ using the
[paid cubic quarter phases](EXPANDER_ZERO_REFLECTION.md#quarter-phase-selector-stage).
For each control coordinate exactly half the row labels require the extra
controlled-S factor. Their common forward/actual-inverse bounds are
$`78qm`$ T gates, $`86qm`$ Cliffords, at most 84 T layers and 3qm
returned dirty helpers, with private physical T targets.

Thus one simultaneous ledger for (2) and its actual native inverse is

```math
\begin{aligned}
T&\le78qm+8m+8+28\ell_C
 \le78qm+36m+28q+8,\\
G&\le86qm+4L_A+36m+26+36\ell_C,\\
D_T&\le116,\qquad
b_{\rm private}\le\max\{3qm,\ell_C+1,2\}\le3qm.
\end{aligned}                                                   \qquad(6)
```

There are zero clean inputs and zero error. The m selectors are principal
arbitrary inputs; including them gives dirty width at most $`m+3qm`$
outside the address and specified row bank. In particular
$`T,b=O(pm)`$ and $`G=O(m^2+pm)`$, within the fixed polynomial
overheads of the retained fusion interface. These are native component
bounds, not a complete-source schedule. All identities include coherent
bank/selector inputs and arbitrary external references.

## 3. Exact global-flip echo and an even-bank obstruction

Let $`F=X^{\otimes m}`$ on the selectors. It commutes with $`L_C`$,
replaces h by $`h+w\pmod2`$, and replaces $`n_C`$ by
$`N_C-n_C`$. Using the invariances from Section 1 gives

```math
E_C=\widetilde R_C^\dagger F\widetilde R_CF
   =i^{N_C-2n_C}L_C^{w\bmod2}.                              \qquad(7)
```

This is the full arbitrary-bank relative operator, not a singleton-bank
calculation. One clean signal with two controlled F walls, the stage (2),
two signal Hadamards and a literal Z yields

```math
K_C=(\widetilde R_C-F\widetilde R_CF)/2,\qquad
K_C^\dagger K_C=(I-\mathrm{Re}E_C)/2.                          \qquad(8)
```

It adds $`2m+4`$ Cliffords to (6), no T gates and no controlled-stage
assumption. Its accepted fiber may move the bank and selectors.

For every even-weight bank, the permutation in (7) disappears. Moreover

```math
E_C^2=(-1)^{N_C}I
 =(-1)^{\sum_uY_u\sum_{c\in S}u_c}I,                         \qquad(9)
```

because the coefficient of each $`x_c`$ in $`N_C\pmod2`$ is w.
The echo is diagonal and either Hermitian or anti-Hermitian. Its Gram in
(8) is therefore either a diagonal projector or $`I/2`$. On a nonempty
proper even bank, a marked address requires Gram zero and an unmarked
address requires $`I/2`$. A projector differs from $`I/2`$ by one half
on every unmarked computational input; $`I/2`$ differs by one half on
every marked input. Hence (8) has complete Gram error at least one half.

There is also a scalar-filter obstruction. Fix any control $`c\in S`$
and use $`Y=e_u+e_{u\oplus e_c}`$. At every address x the two active-row
values of $`d_C(x,\cdot)`$ have opposite parity. Thus $`N_C`$ is odd,
and varying their arbitrary selectors makes $`n_C\pmod2`$ either zero
or one. Equation (7) consequently has both eigenvalues +i and -i at
every address, including marked and unmarked addresses. Such addresses
exist because $`p\ge2`$.

Any common signal/auxiliary preparation, interleaving and acceptance that
uses only $`E_C,E_C^\dagger`$ and acts trivially on x,Y,s between calls
has identical accepted work action on chosen equal-eigenvalue marked and
unmarked inputs. Their target Gram expectations differ by one half, so
the complete Gram error is at least one quarter. This excludes scalar
processing of this echo; it does not exclude separate interleaved uses of
different stages (2), changed acceptance or other bank dynamics.

## 4. The earlier four-path acceptance also fails after transplantation

Changing the stages does not validate the
[earlier bank-sensitive four-path polynomial](COUPLED_QUARTER_PHASE_SOURCE.md#4-paid-internal-signal-mixing-correct-diagonals-wrong-coherent-gram).
Take $`p=2`$, $`C=(0\to1)`$, $`D=(1\to0)`$, and put
$`A=\widetilde R_C`$, $`B=H_s\widetilde R_DH_s`$.
Use the same matched flips $`A'=FAF`$ and
$`B'=H_sF\widetilde R_DFH_s`$. The controlled-selector walls and
bank-weight phases retain their paid implementations: two stages with
one clean signal add $`24m`$ T gates, $`30m+3`$ Cliffords and at most
24 T layers to the sum of their ledgers (6), giving depth at most 256.
The private pool is reusable and remains within 4m dirty bits including
the selectors for these single-pair stages. At the marked input
sector $`Y=e_0,x=0`$, all phases in (2) vanish and the bank and address
are fixed. The selector permutations still act nontrivially.

Let P swap selectors 1 and 3, let Q swap selectors 2 and 3, and write
$`P_\pm=(I\pm P)/2`$, $`Q_\pm=(I\pm Q)/2`$. On this entire marked
sector,

```math
\begin{aligned}
A&=P_++P_-Z_{s_0},&A'&=P_+-P_-Z_{s_0},\\
B&=Q_++Q_-X_{s_0},&B'&=Q_+-Q_-X_{s_0}.
\end{aligned}                                                   \qquad(10)
```

The global selector Hadamards commute with P and Q. The old polynomial,
with $`w=1`$, consequently becomes

```math
\begin{aligned}
K_4&=\frac{BA+iBA'-iB'A-B'A'}{2\sqrt2}\\
&=e^{-i\pi/4}Q_+P_-Z_{s_0}
 +e^{i\pi/4}Q_-P_+X_{s_0},\\
K_4^\dagger K_4&=P_-Q_+P_-+P_+Q_-P_+.
\end{aligned}                                                   \qquad(11)
```

The cross terms vanish because $`Q_+Q_-=0`$. On selectors (1,2,3), set

```math
|\phi\rangle=(|100\rangle-|001\rangle)/\sqrt2.
```

Then $`P\phi=-\phi`$ and $`\langle\phi,Q\phi\rangle=1/2`$.
For every state of selector zero, (11) has acceptance probability three
quarters, while the marked target is zero. Thus the complete Gram error
is at least three quarters. The witness embeds for every $`p\ge2`$ by
fixing the spectator address coordinates and every other selector block
to zero. The two matchings preserve those blocks; their permutations also
commute with the global Hadamards. This rejects the displayed transplant,
not all signal-interleaved co-moving sources.

## 5. The remaining complete identity

Bare stages (2) and their inverses fix every input on the selector-zero
fiber: h and n_C both vanish there for every bank. Common interlayers
acting only on address, clean signals and returned private work preserve
that fiber and cannot distinguish the empty bank from a bank marking the
chosen address. Their target Gram expectations differ by one half, so
this unmixed grammar has error at least one quarter. A new word therefore
needs paid selector mixing or other bank-sensitive dynamics outside that
grammar. Selector Hadamards and flips commute with Omega and are already
native, but their presence alone does not validate an acceptance word.

The new stage pays bank motion and preserves the lookup charge. It does
not pay its readout. A genuinely changed polynomial in distinct stages,
with explicit paid selector and signal interlayers, must have accepted
block K satisfying

```math
K^\dagger K\,\Pi_-=0,\qquad
K^\dagger K\,\Pi_+=\tfrac12\Pi_+,
\qquad \Pi_\pm=(I\pm\Omega)/2,                              \qquad(12)
```

or the corresponding uniform operator-norm approximation. Signal and
selector operations and complete stages (2) preserve the two charge sectors, but
sector preservation alone supplies neither equality in (12).
Every selector and every bank belongs to this contract; a correct basis
probability or a one-row gap is insufficient. Any controlled use, signal
mixing and actual inverse must have their own native costs.

Only a qualifying complete Gram licenses the paid
[Gram recovery and dirty promotion](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout),
which then remove arbitrary accepted fibers and include return of the bank,
selectors and private work in the final error. Neither the source signal
nor its accepted fiber is promised to return beforehand. Reusing the scalar
filter of (7), or transplanting the unchanged polynomial (11), cannot supply
(12).
