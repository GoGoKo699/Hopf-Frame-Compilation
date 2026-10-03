# Exact shared-prefix query interfaces and their limits

[Bilinear query echo](BLOCKED_BILINEAR_LOOKUP.md#3-two-indicator-echoes-complete-the-exact-query) · [Protected source](PROTECTED_UNARY_SOURCE.md) · [Windowed predicates](WINDOWED_GROUP_PREDICATES.md) · [Current research checkpoint](OPEN_PROBLEM.md#revision-checkpoint-and-completed-audit)

Retaining a dirty prefix indicator does not by itself make a later query
independent of its unknown offset. This chapter identifies that residue,
proves a capacity bound for an explicitly restricted exact cache reader,
and gives a native positive example with one additional clean baseline.
A charged two-group identity then isolates the correction needed to share
a prefix boundary. Its generic correction still contains a fresh lookup.
These statements supply no faster complete-frame construction or T-depth
lower bound.

## 1. The exact residue of a retained-prefix short echo

Let x be an unchanged long prefix and t a shorter local address. Write
$`e_x,e_t`$ for their one-hot vectors, and let
$`D\in\mathbb F_2^{H\times J}`$ be one output-bit table. The
completed native bilinear word F implements

```math
z\longmapsto z\oplus Y^{\mathsf T}DX
```

on arbitrary Y, X, and z, returning all helpers. Let $`I_t`$ toggle
$`X\mapsto X\oplus e_t`$. The chronological short echo

```math
F,\ I_t,\ F^\dagger,\ I_t^\dagger
```

returns X and changes z by $`Y^{\mathsf T}De_t`$. Each dagger is
the actual inverse native word. If the long prefix indicator is retained
as $`Y=d\oplus e_x`$ for arbitrary dirty offset d, the contribution is

```math
Y^{\mathsf T}De_t=D[x,t]\oplus d^{\mathsf T}De_t.
```

The second term still depends on the current local address. A baseline
mask evaluated at $`t_{\rm old}`$ cancels it at $`t_{\rm new}`$
for every d exactly when

```math
D(e_{t_{\rm old}}\oplus e_{t_{\rm new}})=0.
```

Indeed the residual is d paired with this column difference; a nonzero
component supplies a witnessing basis offset. For $`D=I_2`$,
$`x=1`$, $`d=e_0`$, $`t_{\rm old}=0`$, and $`t_{\rm new}=1`$,
the retained Y is $`(1,1)`$ and the new short echo contributes one.
The old baseline also contributes one, so their XOR is zero although
the desired $`D[1,1]`$ is one.
If a selected matrix block changes as well, compare the actual vectors
$`D_{\rm old}e_{t_{\rm old}}`$ and
$`D_{\rm new}e_{t_{\rm new}}`$.

This is an exact XOR-query diagnosis. It does not assert that arbitrary
noncommuting logical words admit even this scalar cancellation rule.
The established four-corner bilinear query cancels both dirty offsets
while all its addresses are preserved; its correctness is unchanged.

## 2. A dirty-only cache cannot carry an exact prefix-blind reader

Here is the restricted interface used in the next two sections. The
original prefix x is preserved as a read-only control. Encoding has the
completed form

```math
\sum_x|x\rangle\langle x|\otimes E_x,
```

where $`E_x`$ acts only on the cache. It does not encode into logical
payload, swap x into the cache, or alter the physical prefix. After this
encoding, the processor has no access to x, including CNOT or other
Clifford access. Its local query address and payload inputs are independent
of x and may take every allowed value. A final decoder may depend on x
but acts on the cache only, so it cannot change the payload output.

Suppose first the entire cache consists of b arbitrary dirty qubits.
The contract must hold on the maximally mixed cache input, for which

```math
E_x\frac{I}{2^b}E_x^\dagger=\frac{I}{2^b}
```

for every x. The processor then receives the same state for every basis
prefix. Its payload output cannot depend on x, and a cache-only decoder
cannot change that output marginal. Thus this interface cannot implement
a nonconstant exact prefix-dependent query. The argument already follows
from allowed inputs; no assumption about typical dirty states is made.

## 3. Initialized cache bits bound the number of exact row signatures

Now give the cache c initialized qubits and b arbitrary dirty qubits.
All initialized degrees of freedom used by the encoder count in c,
including a protected source if it participates in encoding. On a
maximally mixed dirty input the initial cache support is

```math
\mathcal S_0=|0^c\rangle\otimes\mathbb C^{2^b},\qquad
\mathcal S_x=E_x\mathcal S_0,\qquad
\dim\mathcal S_x=2^b.
```

Fix one table and one encoder, with an exact reader for every local
address t. Let the row signature $`r_x`$ be the list of all its required
output bits over those addresses. Let R be the number of distinct row
signatures. The reader is an XOR query, correct for every payload and
dirty input; in particular output zero must become the specified bit.

If $`r_x\ne r_y`$, choose an address t where the bits differ. For any
vectors in $`\mathcal S_x`$ and $`\mathcal S_y`$, give the common
processor this same t and the same zero output. Its resulting output
bits are orthogonal. This is true before the final decoder, since that
decoder acts on the cache only. The processor preserves inner products,
so the original cache vectors are orthogonal. Therefore distinct rows
require pairwise orthogonal supports of dimension $`2^b`$, and

```math
R2^b\le2^{b+c},\qquad R\le2^c,\qquad c\ge\lceil\log_2R\rceil.
```

Additional arbitrary dirty bits increase both dimensions by the same
factor. They do not replace initialized capacity in this interface.
The proof is a direct orthogonal-program argument, consistent with
[Nielsen and Chuang, *Programmable quantum gate arrays* (1997)](https://arxiv.org/abs/quant-ph/9703032).
The orthogonality principle is inherited; the support-rank calculation
above states its application to this cache contract explicitly.

## 4. A scoped capacity consequence for arbitrary window-row readers

Consider G consecutive target levels and their $`M=2^G-1`$ tree nodes.
For a unary phase modulus $`q\ge8`$, the two labels 1 and 2 give angles
$`2\pi/q`$ and $`4\pi/q`$, both in the canonical interval
$`[0,\pi/2]`$. Thus even this two-angle family has M independent
node-choice bits. Across $`2^d`$ prefix values, arbitrary tables can
have

```math
R=2^{\min(d,M)}
```

distinct complete row signatures. An exact cache architecture that must
answer every node's row query independently therefore needs

```math
c\ge\min(d,2^G-1).
```

At remaining height k, generously count all k remaining logical bits
and both external flags as possible initialized encoded capacity:
$`c\le k+2`$. This includes any source or conditional-zero bank used
by the encoder; it is not an additional reservation. If
$`d\gt k+2`$, the necessary condition becomes

```math
2^G-1\le k+2,\qquad G\le\log_2(k+3).
```

In the uniform low-precision early schedule, each current group already
has height $`\Omega(\log(n+2))`$, while $`k\le n`$. Hence this
exact prefix-blind all-row architecture can cover only O(1) such groups
where $`d\gt k+2`$. When $`d=n-k`$, this is the latter portion
with $`n\gt2k+2`$. No bound of this form follows when $`c\ge d`$:
a clean copy of the prefix itself fits the stated capacity.

This conclusion concerns exact independent row readers for arbitrary
tables. It is not a theorem about every completed Hopf window, nor a
lower bound on lookup depth or complete-frame depth. A fused window
that never exposes all these row queries has a different interface.

## 5. One clean baseline gives an exact positive native toy

Use wires x,t,s,D,B, with D arbitrary and B initially zero. Encode with
the chronological Clifford word

```math
P:\quad\mathrm{CX}(D;B),\quad\mathrm{CX}(x;D).
```

Apply the chronological body

```math
\mathrm{CX}(D;t),\quad\mathrm{CX}(B;t),\quad
\mathrm{CCX}(D,t;s),\quad\mathrm{CCX}(B,t;s),
```

then the actual $`P^\dagger`$. The encoded pair is
$`(D,B)=(d\oplus x,d)`$, so the body gives

```math
t'=t\oplus x,\qquad s'=s\oplus x(t\oplus x).
```

It preserves D and B, and $`P^\dagger`$ returns them exactly. The
second operation uses the changed target t; the retained baseline remains
valid because both controlled terms use that same current t. Literal
exact Toffoli decompositions give an exact native identity on arbitrary
x,t,s,D and reference inputs. If B instead starts at arbitrary b, the
same circuit replaces x everywhere in the displayed body by $`x\oplus b`$.
It therefore does not supply the claimed operation with a dirty baseline.
This small positive example prices one initialized cache bit; it supplies
no scalable lookup or complete-frame depth saving.

## 6. A charged two-group program refresh

Let A be the completed exact dirty indicator $`Y\mapsto Y\oplus e_x`$,
with read-only prefix x and returned private helpers. Use one common
padded program P. For $`i=1,2`$, let the completed shear $`S_i`$ act as

```math
P\longmapsto P\oplus R_i(t_i)Y,
```

preserving x, Y, its local address $`t_i`$, and all helpers. Each body
$`B_i`$ must preserve P, use x only as a read-only control, and act as
identity on Y and the prefix-helper cache independently of their values.
These standing assumptions imply $`[B_i,A]=0`$ on the full input space. Merely
returning Y does not suffice: control by Y can violate the commutation.
Both bodies preserve $`t_1`$;
$`B_2`$ preserves $`t_2`$, while $`B_1`$ may change a logical target
that belongs to $`t_2`$. The bodies may act on an arbitrary retained
source and need not commute with each other or with the program shears.

All sequences in this section are chronological. Define the exact XOR
query $`Q_i`$ by

```math
Q_i:\quad S_i,\ A,\ S_i^\dagger,\ A^\dagger.
```

It adds $`R_i(t_i)e_x`$ to P and returns Y. The reference sequence

```math
Q_1,\ B_1,\ Q_1^\dagger,\ Q_2,\ B_2,\ Q_2^\dagger
```

equals the retained-prefix sequence

```math
S_1,\ A,\ S_1^\dagger,\ B_1,\
Q_\Delta,\ B_2,\ S_2,\ A^\dagger,\ S_2^\dagger,
```

where the full-input correction query is

```math
Q_\Delta:\quad P\longmapsto
P\oplus[R_1(t_1)+R_2(t_2)]e_x.
```

In this correction, $`t_2`$ is evaluated after $`B_1`$. To prove the
identity, start with arbitrary program p and indicator y. The initial
three words leave $`Y=y\oplus e_x`$ and
$`P=p\oplus R_1(t_1)e_x`$, exactly the program seen by the first
reference body. After that body, the correction changes P to
$`p\oplus R_2(t_2)e_x`$. The second body therefore receives its
correct current program. Its own address is preserved, so the final
three words return P to p and Y to y. Each helper has returned at the
stated completed boundary. This basis identity, with arbitrary source
unitaries in the bodies, extends to coherent inputs and references.

Equivalently, the interstitial chronological word obeys

```math
A^\dagger,(S_1^\dagger S_2),A
=(S_1^\dagger S_2),Q_\Delta.
```

Parentheses denote chronological composition at the current boundary;
this never moves $`S_2`$ through $`B_1`$. Actual inverses retain literal phases.
The operator identity holds for arbitrary P and Y. If a body's intended
logical interpretation requires a clean-loaded one-hot program, that
promise is needed only for that interpretation at its body boundary.

For a reused activity flag, insert its exact boundary transition
$`C:h\mapsto h\oplus\gamma`$ after $`B_1`$ in the retained word.
Both bodies preserve h. Its wire is disjoint from x, Y, and the prefix
helpers; C commutes with A. The cached controls of gamma exclude P and Y
and stay fixed throughout the separator. If $`F_i=R_i(t_i)e_x`$ is the ungated row, the correction
must use the new h and add

```math
\Delta=(h\oplus\gamma)F_1\oplus hF_2
=h(F_1\oplus F_2)\oplus\gamma F_1.
```

It replaces the old active program by the new active program without
storing old h. The same boundary argument proves this variant. An exact
generic implementation uses the chronological shear
$`S_\Delta=C^\dagger,S_1^\dagger,C,S_2`$ in the usual A echo,
with gated $`S_i`$ and C evaluated on the boundary controls. Every C
and actual inverse inside this implementation is charged.

## 7. Every correction and live register remains charged

For $`C=T,G`$ or the serial T-depth upper ledger, the retained word has

```math
C_{\rm retained}\le2C(A)+2C(S_1)+2C(S_2)
 +C(Q_\Delta)+C(B_1)+C(B_2).
```

For $`E\ge2`$ bodies satisfying the same commutation and own-address
preservation contracts, the program-boundary argument gives one outer A
pair, endpoint shears costing $`2C(S_1)+2C(S_E)`$, and $`E-1`$ adjacent
corrections, each evaluated after its preceding body. Allowing two shears
per body is a conservative bound; intermediate shears are charged inside
the corrections:

```math
C_{\rm retained}\le2C(A)+2\sum_{i=1}^E C(S_i)
 +\sum_{i=1}^{E-1}C(Q_{\Delta_i})+\sum_{i=1}^E C(B_i).
```

For $`E=1`$, the single load/body/unload word instead uses four $`S_1`$
appearances and one A/actual-inverse pair.

Add each actual activity transition C to these ledgers when flags change.
Correction queries include every internal C/inverse they use. A generic
implementation uses $`S_1^\dagger S_2`$ in the usual indicator echo,
incurring another A/inverse pair and both shear directions. Classical
table differencing alone gives no long-prefix depth saving.

Reserve Y until the final A inverse and common P outside both target
blocks. Any conditional-zero promise follows from an existing guard and
holds whenever its body is active. Predicates cannot retest the live
loaded P. Correction helpers are disjoint from P, Y, source, and selectors;
private pools are reused only after exact return. All source and guard
work is charged. This local identity proves no new full-frame allocation
or width theorem.

If the two programs occupy distinct banks, the interstitial correction
unloads the old bank and loads the new one. Its table is their
concatenation, not a same-bank XOR cancellation. Both banks and every
copy or decoder must then be charged.

## 8. A correction needs a proved structural saving

Take $`q\ge8`$ and a constant first translation label
$`a_0=q/4`$. Its real Hopf rotation is

```math
R(\pi/2)=XZ,\qquad |t\rangle\longmapsto(-1)^t|t\oplus1\rangle.
```

For this reduction choose that literal Clifford as the exact first logical
body on its active zero-suffix sector; no equality to a generic unary-source
word on arbitrary source characters is assumed.
The first body changes the later address with its literal sign preserved.
Let the second one-hot row $`e_{a(x,t)}`$ be
otherwise arbitrary within the legal angle table class. Its correction
is $`e_{a(x,t)}\oplus e_{a_0}`$. One fixed Clifford X on program bit
$`a_0`$ converts this correction query into the fresh query for
$`e_{a(x,t)}`$. Address change has not forced table correlation.
The reduction is not a depth lower bound or an obstruction to structured families.

There is a precise cheap promised case. Suppose one correction output
is $`\Delta(x,t)=u(x)t`$ with u a known affine parity. Let $`P_u`$
toggle an arbitrary dirty helper d by u using CNOTs and, for a constant
term, X. The chronological native word

```math
P_u(d),\ \mathrm{CCX}(d,t;p),\
P_u^\dagger(d),\ \mathrm{CCX}(d,t;p)
```

adds $`[(d\oplus u)t\oplus dt]=u(x)t`$ to output p and returns d.
It uses two exact Toffolis: $`T=14`$, $`D_T\le8`$, and
$`G=O(|x|+1)`$. The local address t may have changed before this
refresh, but remains fixed during it. The identity holds on arbitrary
dirty, output, and reference inputs. This processor explicitly accesses
x through charged Clifford gates, so Sections 2–4 do not apply to it.

Rank one alone does not give that ledger: for an arbitrary function
$`f(x)`$, querying $`f(x)t`$ at $`t=1`$ recovers a fresh query for
f. The arbitrary legal angle family provides no uniform affine-parity
factorization. Stop the generic difference-only route unless it supplies
a new generic lookup construction or proves a cheap factorization
uniformly over the prescribed family, including its native decoder and
actual inverse. Constant savings for one pair do not establish a saving
across a growing window.

## 9. Scope, evidence, and the next admissible interface

The rank bound does not apply to later direct prefix access, even Clifford access;
encoders that change or move the logical prefix; encoding logical payload;
and decoders that touch payload. Restricted tables have their actual R,
and approximate or probabilistic readers require different bounds. The
argument does not convert any of these excluded routes into an obstruction.

The charged identity isolates the unresolved generic correction cost.
A fused window may avoid independent row readers. Its next bounded test
must address immediate eight-mode closure beyond the existing
[magic-basis factorization and controlled coupling](ENDPOINT_TREE_TRANSPORT.md#8-the-eight-mode-root-retains-a-controlled-coupling).
An exact rewrite of the retained-source word requires every source input.
A different body may act correctly only on the ideal source character
if a new complete error comparison and inactive/work-return proof suffice.
An ideal four-mode factorization or constant pair saving alone does not
improve the remaining linear-depth ledger. No faster frame theorem follows.

The eight [bounded native checks](../tests/test_shared_prefix_query_audit.py)
cover the short-echo residue, a stale baseline, correct completed local
queries with changed targets, a failed delayed whole-word echo, and the
clean-baseline toy. Added checks compare charged two-group words with a
signed XZ body, a reused-h transition with $`\gamma=1`$, and the affine
refresh's 14-T/eight-layer schedule. Complete matrices on at most eight
wires include arbitrary dirty inputs, literal phases, and actual inverses.
The fusion fixtures use the reduced feature $`Y\mapsto Y\oplus x_0x_1`$,
not a full one-hot prefix indicator. Missing or premature corrections and
an omitted $`\gamma F_1`$ term fail explicitly. These checks neither
prove the symbolic capacity bound nor emit complete Hopf groups or a
scalable lookup, and establish no unrestricted Hopf obstruction.
