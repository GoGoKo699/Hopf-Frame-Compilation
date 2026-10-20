# Complete bank difference: a paid nonlinear merge and its depth cost

[Paid affine detuning](PARALLEL_AFFINE_DETUNING.md) ·
[Fixed detuning interaction](SELECTIVE_INVERSION_DIRTY_INDICATOR.md#53-joint-clifford-decoding-of-the-relative-call-still-stops) ·
[Retained dirty conjunction](READ_ONLY_THREE_LAYER_RIGIDITY.md#7-native-layer-fixtures-and-the-clean-qubit-distinction) ·
[Logarithmic-depth indicator](NONUNIFORM_DIRTY_INDICATOR.md)

Bank motion permits a concrete nonlinear operation outside the read-only
Pauli-echo grammar: two complete Boolean vector queries can be multiplied
componentwise by a four-pulse finite difference. The operation has a
literal native word, needs no clean input, and returns its arbitrary
signal bank. It is the retained eight-T dirty conjunction written using
the fixed interaction already present in the detuning call, not a new
native primitive.

For equality banks, fixed row and column CNOT encodings give a complete
indicator recurrence with linear exponential count and width. Its naive
balanced composition has an O(b squared) depth schedule, so it does not
improve the retained logarithmic-address-depth indicator or R96. The
first unpaid improvement is a different composition or whole-word
synthesis, not bank motion or cleanup of this merge.

## 1. The literal four-pulse identity

Let x be an arbitrary address, let Y and s be arbitrary r-bit banks with
r positive, and let f,g map addresses to r-bit Boolean words. All work
and references may be correlated. Define the complete ideal queries

```math
F_f|x,s\rangle=|x,s\mathbin\oplus f(x)\rangle,
\qquad
F_g|x,s\rangle=|x,s\mathbin\oplus g(x)\rangle.
```

Both act trivially on Y and return their private work. The fixed
interaction on one coordinate is

```math
D(Y_u,s_u)=T_{s_u}T_{Y_u\oplus s_u}^{\dagger}
          =\omega^{-Y_u+2Y_us_u},\qquad \omega=e^{i\pi/4}.
                                                               \qquad(1)
```

This is the non-Clifford local factor of the fixed detuning R_0. The
other R_0 factors are its displayed signal S powers and literal scalar;
none is silently retained in D. Equation (1) is implemented directly:
CNOT from s_u into Y_u, simultaneous T on s_u and T-dagger on Y_u,
then the same CNOT. Thus D_r, the product of these interactions over u,
has exactly 2r T gates, one T layer and 2r CNOTs, with no helper or
unrecorded phase. Its actual inverse has the same bounds.

Execute the following word chronologically:

```math
D_r,\quad F_f,\quad D_r^\dagger,\quad F_g,\quad
D_r,\quad F_f^\dagger,\quad D_r^\dagger,\quad F_g^\dagger.
                                                               \qquad(2)
```

For a fixed computational x,Y,s, the signal values seen by the four
diagonal words are s, s+f, s+f+g and s+g, with additions binary. On one
coordinate the accumulated exponent modulo eight is

```math
\begin{aligned}
&(-Y+2Ys)-(-Y+2Y(s\oplus f))\\
&\quad+(-Y+2Y(s\oplus f\oplus g))
       -(-Y+2Y(s\oplus g))\\
&\qquad\equiv4Yfg\pmod8.
\end{aligned}                                                   \qquad(3)
```

If either f or g is zero the difference vanishes. If both are one it
is 4Y times an odd sign, giving the same residue independently of s.
The signal toggles also cancel. Consequently (2) is exactly

```math
(-1)^{\sum_uY_uf_u(x)g_u(x)}\otimes I_s\otimes I_{\rm private}.
                                                               \qquad(4)
```

Put Hadamards on every Y coordinate before and after (2). The result is
the complete vector query

```math
|x,Y,s,w\rangle\longmapsto
|x,Y\mathbin\oplus(f(x)\odot g(x)),s,w\rangle.          \qquad(5)
```

Here the product is componentwise over the binary field. The Hadamards
give actual bank motion within the completed word. The computational
proof of (3) is a literal unitary identity, so (5) holds on arbitrary
superpositions and with arbitrary external references. No signal or
other output is initialized.

## 2. Complete resource, error and return recurrence

Suppose emitted source words approximate F_f and F_g tensor identity on
their private work in full operator norm by epsilon_f and epsilon_g.
Use their actual emitted adjoints at the inverse positions in (2).
Let their simultaneous common forward/actual-inverse bounds be
T_f,d_f,G_f,w_f and T_g,d_g,G_g,w_g. The source private pools must be
disjoint from x,Y,s; padding and sequential reuse of one pool of the
larger size are allowed. Then (5) has

```math
\begin{aligned}
T_{f\odot g}&\le2T_f+2T_g+8r,\\
D_T(f\odot g)&\le2d_f+2d_g+4,\\
G_{f\odot g}&\le2G_f+2G_g+10r,\\
w_{f\odot g}&\le r+\max\{w_f,w_g\},\\
\epsilon_{f\odot g}&\le2\epsilon_f+2\epsilon_g,
\qquad a_{\rm clean}=0.
\end{aligned}                                                   \qquad(6)
```

The fixed Clifford charge is 8r CNOTs from the four D words and 2r
output Hadamards. The r arbitrary signals are additional work beyond
the original output Y. Each source and its actual adjoint occur once;
the inverse source error is the same as its forward error by unitary
invariance. Unit-norm product telescoping proves the last line, including
all logical columns, address and private-work return, signal return and
references. Actual approximants need not return exactly before reuse;
their leakage is already included in the full-operator comparison.

The ideal sources preserve Y and the unused address coordinates, but
their native implementations may move those registers temporarily if
their complete operator certificate permits it. Nothing in (6) treats
an intermediate dirty state as a fresh initialized input. The same
ledger bounds the actual inverse of the complete emitted word.

For affine components

```math
f_u(x)=\kappa_{f,u}\oplus\bigoplus_{j\in S_{f,u}}x_j,
\qquad
g_u(x)=\kappa_{g,u}\oplus\bigoplus_{j\in S_{g,u}}x_j,
```

the queries are direct CNOT words and X gates on s. With X charged as
HSSH, (6) specializes to

```math
\begin{aligned}
T&=8r,\qquad D_T\le4,\qquad w=r,\qquad\epsilon=0,\\
G&\le10r+2\sum_u(|S_{f,u}|+|S_{g,u}|)
       +8\sum_u(\kappa_{f,u}+\kappa_{g,u}).
\end{aligned}                                                   \qquad(7)
```

All four T layers have distinct physical targets across coordinates;
shared address controls only serialize charged Clifford operations.
For one coordinate and two physical address controls, G is at most 14.
This is the retained dirty-AND construction in the coordinates of (1).

The affine toggles used in (7) are actual read-only-address Cliffords.
They are not address-independent bank-only interlayers. If a proposed
detuning grammar permits only the latter, those toggles require an
additional paid realization; (7) does not grant them for free within
that stricter syntax. The unrestricted native model does permit their
literal CNOT implementations.

## 3. A complete Cartesian equality-bank merge

Let b=a+c with a,c positive, m=2^b, and index Y and the arbitrary signal
bank by pairs (v,w) with v an a-bit label and w a c-bit label. Split the
address as x_L,x_R. Choose

```math
f_{v,w}(x)=[x_L=v],\qquad g_{v,w}(x)=[x_R=w].
```

Their componentwise product is precisely the full b-bit one-hot word.
The larger vector queries f and g need no initialized fanout. For f,
select one pivot in each row v and let C_L be the CNOT word from that
pivot to every other signal in the row. These CNOTs commute and their
number is m-2^a. Conjugating the complete a-bit indicator on the row
pivots by C_L turns its selected pivot X into the product of X on that
entire row. Every other initial signal value remains arbitrary. The two
C_L words cost 2(m-2^a) Cliffords. The analogous column encoding costs
2(m-2^c), and supplies g from the c-bit indicator.

The row and column words execute at their separate positions in (2).
Their overlapping signal supports are not claimed to be parallel source
calls. Each smaller indicator returns its private work before the next
ideal call, so the two private pools can be reused. Applying (6) gives
the complete recurrence

```math
\begin{aligned}
T_b&\le2T_a+2T_c+8\,2^b,\\
d_b&\le2d_a+2d_c+4,\\
G_b&\le2G_a+2G_c+18\,2^b-4(2^a+2^c),\\
w_b&\le2^b+\max\{w_a,w_c\},\\
\epsilon_b&\le2\epsilon_a+2\epsilon_c.
\end{aligned}                                                   \qquad(8)
```

The original output Y is excluded from w_b; the entire m-bit signal
bank is included. The child outputs are selected pivots already inside
that bank, not additional initialized or copied registers. All child
private work remains additional and is charged by the maximum in (8).
The ideal result returns the complete signal bank, including nonpivots,
and every source-private register.

For b=1 the indicator is Clifford: CNOT from x into each output and X
on output zero. Thus T_1=d_1=w_1=epsilon_1=0 and G_1 at most 6. With
the balanced choice a=floor(b/2), c=ceil(b/2), (8) supplies, for example,

```math
T_b\le32\,2^b,\qquad G_b\le72\,2^b,
\qquad w_b\le2\,2^b,\qquad\epsilon_b=0.               \qquad(9)
```

For the count bounds, verify b at most four directly in (8). For b at
least five,

```math
\frac{2\,2^a+2\,2^c}{2^b}
=2^{1-a}+2^{1-c}\le\frac34.
```

The constants 32 and 72 absorb respectively the added 8 and 18 in
this normalized induction; dropping the negative G correction only
weakens the bound. The width induction uses a at least one and
c at most b-1. These are simultaneous safe ledgers, not optimized
constants or new indicator count bounds.

For b a power of two, the displayed composed depth schedule satisfies

```math
d_b\le4d_{b/2}+4
\le\frac43(b^2-1).                                    \qquad(10)
```

The recurrence tree has height at most the ceiling of log_2 b, giving
an O(b squared) schedule for all b. Equation (10) is an upper bound for
this naive recursive word,
not a lower bound on the operator, bank-moving detuning mechanisms or
unrestricted native depth. Reordering or whole-word resynthesis can
leave this schedule. Substituting the retained logarithmic-depth child
indicators in one merge of (8) retains logarithmic depth and changes
only constants; it supplies no improved resource frontier.

## 4. What was paid and what remains open

The fixed non-Clifford interaction in R_0 is useful after genuine bank
motion: it supports a complete nonlinear scatter, with arbitrary-input
return and a simultaneous native ledger. Its immediate finite-difference
use is an ordinary dirty conjunction. The Cartesian construction pays
its full continuation and exposes its depth cost rather than leaving
cleanup or a selected router implicit.

The next constructive question is a complete interference mechanism or
joint resynthesis that avoids the sequential source appearances in
(8), while keeping the whole-space error and count/width guarantees.
A claim that two source calls can run simultaneously must supply its
own native word: they share the signal bank and occur on different
sides of noncommuting D operations. Conversely, (8) is not an
obstruction to such a new word. The read-only Pauli-echo result has its
separate hypotheses, and neither result changes the unrestricted depth
frontier.

No new fixture is needed for this change of coordinates: (1) is the
retained one-band local phase word, and (3) is an exact Boolean phase
identity. The new statements are the complete composition and resource
accounting above; they do not assert a new primitive or a new general
lower bound.
