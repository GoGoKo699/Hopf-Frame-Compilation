# Phase-signal and coherent-frame interfaces for complete queries

[Complete approximate transfer](APPROXIMATE_QUERY_FRAME_TRANSFER.md) ·
[Earlier filter limits](APPROXIMATE_INDICATOR_FILTER_LIMITS.md) ·
[Retained dirty indicator](NONUNIFORM_DIRTY_INDICATOR.md)

Two mechanisms move the missing operation without yet supplying a faster
complete indicator. A coherent OR phase signal distinguishes zero on an
initialized vector, but its natural reflection echo does not become a
scalar on arbitrary dirty work. A complete coherent Pauli-frame update
for indicator teleportation already contains the original indicator,
with a paid two-call extraction. These are interface results, not
unrestricted T-depth lower bounds. R96 and the general frontier are
unchanged.

Throughout, $`p\ge1`$, $`m=2^p`$, and the desired phase and indicator are

```math
D_p|x,Y\rangle=(-1)^{Y_x}|x,Y\rangle,
\qquad I_p=H_YD_pH_Y.
```

The address x and bank Y are arbitrary, as are all additional dirty
inputs and references. A useful approximate source must provide a
complete operator-norm estimate, actual inverses, zero initialized work,
and simultaneous $`O_\tau(m)`$ T/dirty, $`O_\tau(m^2)`$ Clifford and
$`o(\log(p+2))`$ T-depth bounds at each fixed requested error tau.

## 1. What the fanout OR reduction actually supplies

Høyer and Špalek's [*Quantum Fan-out is Powerful*](https://theoryofcomputing.org/articles/v001a005/v001a005.pdf),
Theorem 3.2, parallelizes commuting gates with initialized copies.
Their Lemma 5.1 gives an exact OR reduction from p bits to
$`\ell=\lceil\log_2(p+1)\rceil`$ initialized signal bits using arbitrary
one-qubit rotations. In the following notation its signal word is

```math
\begin{gathered}
R_w=\bigotimes_{j=1}^{\ell}
 H\begin{pmatrix}1&0\\0&e^{2\pi i w/2^j}\end{pmatrix}H,
\qquad w=|x|,\\
R_0=I,\qquad \langle0^\ell|R_w|0^\ell\rangle=0
\quad(1\le w\le p).                                  \qquad(1)
\end{gathered}
```

Indeed, write $`w=2^a(2b+1)`$. The factor j=a+1 sends its zero
input to one. This proves the signal property in (1); it does not give
the same zero-test guarantee on arbitrary initial signal states.

Two resources of that construction need separate native treatment.
Clifford fanout is allowed here, but its dirty targets are XOR masks,
not initialized copies. For example, applying phases to dirty copies
$`d_j\oplus s`$ and cancelling the original phases leaves the exponent

```math
\phi\sum_jx_j[(d_j\oplus s)-d_j]
 =\phi s\sum_jx_j(1-2d_j),
```

which still depends on the arbitrary masks. This simple echo does not
provide the intended Hamming-weight phase.

The rotation alphabet also matters. The paper's Remark 3.5 preserves
constant depth over a fixed irrational Z rotation by parallelizing its
commuting powers. That rotation is not an exact Clifford+T primitive.
The exact signal in (1) includes a one-qubit phase of angle pi/8 once
p is at least eight. Its entry $`e^{i\pi/8}`$ lies outside
$`\mathbb Q(e^{i\pi/4})`$, so even exact returned work cannot make
that displayed rotation an exact native word. Approximate replacement
requires a complete accuracy-dependent count, depth and return ledger;
the source's fanout-depth statement is not such a ledger. This does not
exclude a different native synthesis of the completed signal operation.

## 2. The dirty reflection echo has a fixed-subspace obstruction

Grant the ideal signal (1), so that no rotation-synthesis error is used
in this diagnosis. Put

```math
P=|0^\ell\rangle\langle0^\ell|,
\qquad S=I-2P,
\qquad E_w=R_w^\dagger S R_w S.                       \qquad(2)
```

This is the actual-adjoint echo. It gives $`E_0=I`$. For nonzero w,
the vectors $`|0^\ell\rangle`$ and $`R_w^\dagger|0^\ell\rangle`$
are orthogonal by (1). Consequently E_w is minus identity on their
two-dimensional span and identity on its orthogonal complement. If
$`\ell\ge2`$, that complement is nonzero, and

```math
\|E_w-(-I)\|_{\mathrm{op}}=2
\qquad(1\le w\le p).                                 \qquad(3)
```

Thus the echo implements the OR sign on the initialized zero signal,
but fails maximally on an allowed dirty state. It does not return every
dirty state with the desired common phase. The exceptional one-bit
signal has no complementary subspace and is not covered by (3).
Replacing R_w by a unitary within full-operator error epsilon and using
its actual adjoint changes (2) by at most twice epsilon. Its dirty error
is therefore at least $`2-2\epsilon`$; accurate native approximation
of the signal does not repair this echo.

More generally, fix an address sector and consider any product of q
conjugates of S on the signal and arbitrary helpers. If the helper
dimension is h, each conjugated
reflection fixes a subspace of codimension h. The intersection of all
their fixed spaces therefore has dimension at least

```math
h(2^\ell-q).
```

Every vector in this intersection is fixed by the full product,
regardless of its ordering. Hence a literal approximation to minus
identity with operator error less than two requires

```math
q\ge2^\ell.                                           \qquad(4)
```

The conjugators may mix all helper and signal wires; an identity helper
extension does not weaken the rank fraction. The premise is a product
of these reflected signal projectors, with no additional scalar supplied
depending on the address. This application preserves the chosen address
sector; a joint word mixing address sectors is outside its premise.
Equation (4) counts occurrences in this
interface. It is not an unrestricted T-count or T-depth lower bound:
joint resynthesis, other projectors or a separately implemented scalar
can leave the premise. In particular, a serial readout construction
cannot claim sublogarithmic native depth merely from the fanout depth
of its signal preparation.

## 3. Full coherent Pauli-frame transport contains the indicator

Consider arbitrary coherent Pauli-frame keys
$`a,c\in\mathbb F_2^p`$ and $`b,d\in\mathbb F_2^m`$. Fix the literal
ordered Pauli

```math
P_{a,b,c,d}=X_x^aX_Y^bZ_x^cZ_Y^d,
\qquad F=D_pPD_p^\dagger.                             \qquad(5)
```

F denotes the coherent controlled word over all key values; it preserves
the keys. The standard conjugated correction appears in Gottesman and
Chuang's [*Quantum teleportation is a universal computational primitive*](https://arxiv.org/abs/quant-ph/9908010),
equation (4). Here it is promoted to a full coherent-key interface. The
argument grants any resource preparation and asks what a compiler of
this completed correction still contains. It does not infer arbitrary-key
action from a promise on one teleportation resource state. No key is
assumed zero and the Pauli's order is not identified only up to a
key-dependent sign.

Define $`K=FP^\dagger`$ and abbreviate the full key by kappa.
Applying the rightmost factor first gives

```math
\begin{gathered}
K|x,Y,\kappa\rangle=(-1)^{g(x,Y,a,b)}|x,Y,\kappa\rangle,\\
g(x,Y,a,b)=Y_x+Y_{x+a}+b_{x+a}.
\end{gathered}
                                                               \qquad(6)
```

The c,d phases cancel with the actual P adjoint, including their signs
on shifted values. Introduce an arbitrary m-bit target V and let L
be the Clifford word $`b\leftarrow b\oplus V`$. Equation (6) implies

```math
K L K^\dagger L^\dagger=(-1)^{V_{x+a}}I.              \qquad(7)
```

Conjugating by the Clifford A that adds a to x removes the shift.
Hadamards on V at both ends give the complete indicator into V. One
literal chronological implementation is

```math
H_V,\ A,\ L^\dagger,\ K^\dagger,\ L,\ K,\ A,\ H_V.
                                                               \qquad(8)
```

Here K is implemented chronologically as $`P^\dagger,F`$, and its
actual inverse as $`F^\dagger,P`$. All keys, x, Y and private work
return. Equations (6)--(8) hold on the complete space, including
address superpositions and every possible dirty key. A common scalar
of F cancels with its adjoint. An arbitrary key-dependent scalar is
not part of the promised transport in (5), and need not cancel when
L changes b.

### Native resources and approximation

Suppose a native implementation of F has common forward/actual-inverse
bounds $`T_F,G_F,d_F,b_F`$, where b_F counts private arbitrary helpers.
The controlled P and its actual adjoint each use p+m CNOTs and p+m
controlled-Z gates. Compiling each controlled-Z as H,CNOT,H gives
$`4(p+m)`$ elementary Cliffords. Thus (8) has common bounds

```math
\begin{aligned}
T&\le2T_F,&G&\le2G_F+10p+12m,\\
D_T&\le2d_F,&b_{\mathrm{dirty}}&\le b_F+2p+3m.
\end{aligned}                                                   \qquad(9)
```

The helper count includes the old Y bank and all four keys, excluding
the specified x address and V output. Both L words, both address shifts
and both output-Hadamard words are charged. No clean qubit is added;
any clean allocation in a different source model would be inherited,
not converted into dirty work. Positive-S inverse costs, if present
inside F, belong in its common bound G_F.

If the physical source instead satisfies

```math
\|\widetilde F-F\otimes I_{b_F}\|_{\mathrm{op}}\le\eta,
```

use its actual adjoint in (8). Unitary telescoping across the two source
calls proves complete indicator error at most $`2\eta`$. This norm
includes all return errors and references; it does not assume the
approximate source returns work exactly between calls.

Conversely (5) compiles F with two D_p calls and the paid Clifford P.
The coherent-frame interface and complete indicator are therefore
equivalent up to constant source-call factors and linear Clifford/dirty
overhead. Preparing a teleportation resource or postponing its corrections
does not by itself pay this complete coherent update.

## 4. The remaining native operation

The phase-signal route still needs a scalar action on every dirty signal
sector and a native precision ledger. The coherent-frame route still
needs the complete update (5); treating frame bits as measurement
outcomes or initialized constants changes its contract. Neither result
excludes a whole-word circuit that avoids these interfaces.

A zero-clean implementation of (5) with the source budgets stated above
would, by (9), improve the complete indicator and enter the paid R91/R98
frame transfer. The reduction exposes exactly what such a teleportation
proposal must compile; it does not construct that compiler. Likewise,
a faster equality test on initialized phase signals is insufficient
until the complete dirty indicator, its actual inverse and all errors
are paid.

The [focused fixtures](../../tests/test_phase_signal_query_limits.py)
check small native reflection and coherent-frame words. The all-size
fixed-subspace theorem, two-call extraction and operator-error transfer
are analytic statements, not consequences of bounded examples.
