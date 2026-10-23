# Coupled quarter-phase sources and complete Gram defects

[Current checkpoint](../../WORKSPACE.md) · [Resource frontier](../../docs/OPEN_PROBLEM.md) ·
[Quarter-phase stages](EXPANDER_ZERO_REFLECTION.md#quarter-phase-selector-stage)

The single-stage source's arbitrary-selector defect motivates internal
selector mixing and then internal signal mixing. Both native words below
are paid. The second has exactly the desired acceptance probability on
every computational input of its one-row sector, but a coherent kernel
still makes its operator Gram error one half. Neither word improves R96.

## 1. Two distinct stages with matched selector flips

Let $`m=2^p`$, $`p\ge2`$. Use the exact stages (51)--(52) in the linked
proof: $`A=R_C`$ for CNOT $`c\to t`$ and $`R_D`$ for CNOT $`d\to e`$.
The address x, row bank Y and m selectors s are arbitrary. Write

```math
H_s=H^{\otimes m},\qquad F_X=X^{\otimes m},\qquad F_Z=Z^{\otimes m},
\qquad B=H_sR_DH_s,
```

```math
A'=F_XAF_X,\qquad B'=F_ZBF_Z=H_sF_XR_DF_XH_s.
```

For one initialized signal a, emit chronologically

```math
H_a,\ F_a,\ R_C,\ F_a,\ H_s,\ F_a,\ R_D,\ F_a,\ H_s,\ Z_a,\ H_a,
\qquad F_a=\prod_u\mathrm{CX}(a,s_u).
```

Here Z is its literal two-S word. No controlled stage is assumed.
The two signal branches give the complete accepted block

```math
K_-=(BA-B'A')/2,\qquad
K_-^\dagger K_-=(I-\mathrm{Re}\,Q)/2,\qquad
Q=A^\dagger B^\dagger B'A'.                              \qquad(1)
```

At q=1, the retained stage ledger is
$`T_R\le98m+16`$, $`G_R\le111m+20`$, $`d_R\le124`$,
with at most 3m returned private bits. Consequently the source and its
actual native inverse obey simultaneously

```math
T\le196m+32,\quad G\le228m+44,\quad D_T\le248,
\quad b_{\rm dirty}\le4m,\quad a_{\rm clean}=1,\quad\varepsilon=0.
                                                               \qquad(2)
```

Dirty width here is outside the address and specified row bank: m selectors
plus one reusable 3m private pool. Every private stage helper returns exactly.
The accepted operator acts on the arbitrary selectors; the clean signal
is projected onto its accepting state, with no pre-acceptance return promise. [Gram recovery](NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#four-call-gram-readout)
would supply their complete return only for a qualifying target Gram.

### Complete arbitrary-bank relative action

For a control coordinate j and bit v, define the selector mask and its
ordinary integer weight

```math
M_j(v)=\{u:Y_u=1,\ u_j\ne v\},\qquad N_j(v)=|M_j(v)|.
```

Write $`X_M=\prod_{u\in M}X_{s_u}`$ and similarly for Z. The full
integer-count identity for the single stage gives

```math
R_D^\dagger F_XR_DF_X
=\sum_{v=0}^1P_{x_d=v}(iX_e)^{N_d(v)}Z_{M_d(v)},
```

```math
B^\dagger B'
=\sum_{v=0}^1P_{x_d=v}(iX_e)^{N_d(v)}X_{M_d(v)}.          \qquad(3)
```

Thus Q is monomial even though B need not be. Its exact action on every
computational column is the following recurrence; masks are viewed as
binary selector strings when XORed:

```math
\begin{aligned}
k&=N_c(x_c)-\sum_{u\in M_c(x_c)}s_u,&x_1&=x+(k\bmod2)e_t,\\
q&=N_d((x_1)_d),&s_1&=s+M_d((x_1)_d),\\
x_2&=x_1+(q\bmod2)e_e,&
k_2&=\sum_{u\in M_c((x_2)_c)}(s_1)_u,\\
x_3&=x_2+(k_2\bmod2)e_t,&
Q|x,s,Y\rangle&=i^{k+q-k_2}|x_3,s_1,Y\rangle.
\end{aligned}                                                  \qquad(4)
```

The additions to addresses and selectors are binary; $`k,q,k_2`$ and the phase
exponent are integers. The order is A-prime, (3), then the actual A-adjoint.
This proves the Gram in (1) on arbitrary coherent banks and references.

### An unwanted fixed vector survives the coupling

Choose opposite directions $`c=0,t=1,d=1,e=0`$ and the bank with only row
zero active. Fix all other address coordinates to zero. Let $`|x,s\rangle`$
denote the two-bit address $`x=x_0+2x_1`$ and the active row selector.
All other selectors decouple. Equation (4) fixes both marked states
$`|0,0\rangle,|0,1\rangle`$ and has the cycle

```math
|1,0\rangle\longmapsto-|2,1\rangle,\qquad
|2,1\rangle\longmapsto i|3,0\rangle,\qquad
|3,0\rangle\longmapsto i|1,0\rangle.                    \qquad(5)
```

The cycle's total phase is +1. In particular

```math
v=(|1,0\rangle-|2,1\rangle-i|3,0\rangle)/\sqrt3,
\qquad Qv=v.
```

For $`\Omega|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle`$, the desired Gram
$`(I+\Omega)/4`$ is zero on the marked states and one half on v.
The direct source annihilates v, so its error is at least one half.
More generally any signal-only spectral processing of Q and its adjoint
acts identically on these two eigenvalue-one inputs. Its Gram error is at
least one quarter by the triangle inequality. This scope allows arbitrary
signal/auxiliary processing but requires every interlayer to act trivially
on the address, bank and selector work. It does not include new address or bank
insertions, or arbitrary signal-interleaved uses of A and B separately.

Using the wrong selector flip is even less informative: B commutes with
F_X, so $`F_XBF_X=B`$. Interference of BA with $`BA'`$ then has exactly
the old single-stage Gram. The conjugated F_Z flip above is essential to
obtaining a changed polynomial at all.

## 2. Changing branch order alone still fails

As a relaxed algebraic trial, replace the second branch by $`A'B'`$.
The relative word becomes $`Q_{\times}=A^\dagger B^\dagger A'B'`$.
This is a different word; (2) does not pay its conditional stage order.
Its full-bank contract already fails even if that control cost is granted.

At p=2, activate rows zero and one, and start with address zero and every
selector zero. Selector labels below are bit masks on those two rows.
Direct branch propagation gives

```math
Q_{\times}|0,0\rangle
=\bigl(i|2,1\rangle+|0,2\rangle+|3,0\rangle+i|1,3\rangle\bigr)/2.
                                                               \qquad(6)
```

Its diagonal matrix element is zero. The difference-source Gram expectation
is one half, while $`Y_0=1`$ requires zero. This gives error at least one
half and embeds at every p at least two by fixing the spectator coordinates
and activating the same two labels. Merely reversing these two branch
orders therefore does not supply the complete source.

## 3. A positive row-selective partial Gram

Flip only selector v: $`F_v=X_{s_v}`$. Keeping the full R_C in the same
one-signal circuit gives

```math
K_v=(R_C-F_vR_CF_v)/2,\qquad
K_v^\dagger K_v=\tfrac12Y_v P_{x_c\ne v_c}\otimes I_s.   \qquad(7)
```

This identity holds for every row bank, not just a one-row promise. To see
it, factor $`R_C=R_{\ne v}r_v^{s_vY_v}`$. All factors commute. If
$`Y_v=0`$ the difference vanishes. Otherwise its remaining factor is
$`Z_{s_v}(I-r_v)/2`$, with unitary left factor R-not-v. Since
$`r_v=I`$ on matching controls and $`r_v=iX_t`$ on mismatches, the
Gram is exactly (7); all other active rows cancel.

Using the full retained stage, both orientations have T-count at most
$`98m+16`$, Clifford count at most $`111m+26`$, depth at most 124,
one clean signal and at most 4m dirty work bits. The six added Cliffords
are two signal Hadamards, two selector CNOTs and the two-S Z word.
This supplies a complete partial projector with an arbitrary accepted
unitary fiber. It does not select the full addressed bank bit or prove
that separate row projectors can be combined at the required budgets.

## 4. Paid internal signal mixing: correct diagonals, wrong coherent Gram

The preceding binary interfaces motivate genuine internal signal mixing.
They should also be checked against the
[balanced two-unitary restriction](COUPLED_ORDER_GRAM_SOURCE.md#balanced-two-unitary-gram-readout).
The following four-path source leaves that specific final binary form.

Let $`w(Y)=\sum_uY_u`$ be the ordinary integer weight, and define operators
on the signal and arbitrary work by

```math
\mathcal A=\mathrm{diag}(A,i^{w(Y)}A'),\qquad
\mathcal B=\mathrm{diag}(B,i^{-w(Y)}B'),\qquad
U_4=H_a\mathcal B H_a\mathcal A H_a.
```

The controlled selector flips implement each diagonal pair as above. The
additional phase is $`\prod_u\mathrm{CS}(a,Y_u)`$, or its actual inverse.
Each private CS uses 12 T gates, 12 Cliffords and at most 12 T layers; private
row helpers permit all physical T targets to be parallel. Its 2m returned
helpers fit inside the reusable 3m pool. The complete source and actual
reverse therefore have

```math
T\le220m+32,\quad G\le252m+43,\quad D_T\le272,
\quad b_{\rm dirty}\le4m,\quad a_{\rm clean}=1,\quad\varepsilon=0.
                                                               \qquad(8)
```

Every bank is preserved; each phase commutes with its corresponding
signal-diagonal stage. Direct propagation through all three signal
Hadamards yields the complete arbitrary-bank polynomial

```math
K_4=\frac{BA+i^{w(Y)}BA'+i^{-w(Y)}B'A-B'A'}{2\sqrt2}.     \qquad(9)
```

At an empty bank this has Gram I/2, as required. On the one-row sector
of (5), it rejects the marked states. Every unmarked computational input
has Gram diagonal exactly one half. This stronger basis-state check still
misses a coherent obstruction.

Let $`P_{\rm unmarked}`$ project onto the six states with $`x\ne0`$, and use the
following unnormalized, orthogonal vectors of squared norm three:

```math
v_1=|1,0\rangle-i|3,0\rangle+i|2,1\rangle,\qquad
v_2=|2,0\rangle+i|1,1\rangle+|3,1\rangle.
```

Substitution of the two controlled quarter-phase operators in (9) gives

```math
K_4^\dagger K_4
=\tfrac34P_{\rm unmarked}
 -\tfrac14\bigl(|v_1\rangle\langle v_1|+|v_2\rangle\langle v_2|\bigr).
                                                               \qquad(10)
```

For example its $`\langle1,0|\cdot|3,0\rangle`$ entry is $`-i/4`$.
The eigenvalues are zero on both v vectors and both marked states, and
three quarters on their remaining four-dimensional complement. Its
distance from $`P_{\rm unmarked}/2`$ is exactly one half. The witness
embeds for every p at least two. Scalar processing of this Gram cannot
separate its marked and unmarked zero eigenvectors. This is a failure of
the displayed signal-interleaved polynomial, not all four-path sources.

## 5. Decision and evidence

Pause the matched-flip word, the displayed crossed-order difference, and
the displayed four-path polynomial. A repair must change their complete
accepted operator and eliminate the coherent kernel; basis probabilities
or scalar renormalization do not do so. The
[matching-family limits](QUARTER_PHASE_MATCHING_LIMITS.md) separately delimit
what can be learned from one commuting dressed matching.

The remaining constructive question is a normalized full-bank Gram from
signal-interleaved stages with genuinely varied generator dynamics. The
[dressed-expander theorem](QUARTER_PHASE_MATCHING_LIMITS.md#dressed-expander-gap)
now proves a uniform one-row spectral gap on its stated admissible lengths.
That component supplies neither an aggregate SELECT ledger nor a full-bank
source. The source must additionally include every bank, remove arbitrary-selector
dependence and control the coherent off-diagonal terms before invoking
paid recovery and fusion. The general complete-frame frontier is unchanged.

The [exact fixtures](../../tests/test_coupled_quarter_phase_source.py) retain
all 64 columns and both actual inverses of six-wire native one-row sources:
40 T/54 Clifford gates for (1), and 40 T/57 Cliffords for (9), where
S-adjoint is emitted as three S gates. These are fixed-row fixtures, not
arbitrary-bank native emitters. Separate Gaussian-integer arithmetic checks
all 16 two-bit bank choices and every address/selector column of (3)--(4),
the crossed-order witness, the partial Gram (7), and the coherent kernel
in (10). These finite checks accompany the analytic identities; they are
internal verification, not external mathematical review.
