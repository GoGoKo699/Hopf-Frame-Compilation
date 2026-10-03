# Exact shared-prefix query interfaces and their limits

[Bilinear query echo](BLOCKED_BILINEAR_LOOKUP.md#3-two-indicator-echoes-complete-the-exact-query) · [Protected source](PROTECTED_UNARY_SOURCE.md) · [Windowed predicates](WINDOWED_GROUP_PREDICATES.md) · [Current research checkpoint](OPEN_PROBLEM.md#revision-checkpoint-and-selected-next-test)

Retaining a dirty prefix indicator does not by itself make a later query
independent of its unknown offset. This chapter identifies that residue,
proves a capacity bound for an explicitly restricted exact cache reader,
and gives a native positive example with one additional clean baseline.
None supplies a faster complete-frame construction or a T-depth lower bound.

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

## 6. Scope and the next admissible interface

The rank bound does not apply to later direct prefix access, even Clifford access;
encoders that change or move the logical prefix; encoding logical payload;
and decoders that touch payload. Restricted tables have their actual R,
and approximate or probabilistic readers require different bounds. The
argument does not convert any of these excluded routes into an obstruction.

A new proposal must either retain explicitly charged prefix access or
specify a fused window contract that avoids independent all-row readers.
Its decoder, actual inverse queries, changed local addresses, retained
source, and dirty/reference return remain part of the complete ledger.
No reusable faster query identity is established here; logical stages
also retain their separate linear-depth allowance.

The five [bounded native checks](../tests/test_shared_prefix_query_audit.py)
cover the short-echo residue, a stale baseline, correct completed local
queries with changed targets, a failed delayed whole-word echo, and the
clean-baseline toy. They compare full matrices on at most seven wires,
including dirty inputs, literal phases, and actual inverses. They do not
prove the symbolic capacity bound, emit complete Hopf groups, or establish
an unrestricted Hopf obstruction.
