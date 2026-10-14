# Shallow word queries and guarded Hadamard layers

[Bilinear query](../../docs/PARALLEL_DIRTY_LOOKUP.md#5-a-bilinear-query-reduction) · [Read-only conjunctions](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences) · [Linear-size indicator](NONUNIFORM_DIRTY_INDICATOR.md)

An arbitrary word query can run its output bits in parallel when each
bit has a separate dirty indicator and private helper pool. The shared
address enters only through Clifford gates. This trades count and width
for depth; it is useful when the queried table is exponentially smaller
than the full frame. The same complete-input query supplies a Boolean
phase without a clean phase ancilla and a variable-suffix Hadamard layer
without initialized copies of its selectors.

All identities below preserve literal phases, use actual circuit
inverses, and hold for arbitrary dirty states and external references.
T-depth permits arbitrary intervening Clifford circuits; their
elementary gate counts remain charged. These are native primitives,
not by themselves a complete-frame compiler.

## 1. Parallel output bits with a read-only address

Let an r-bit address x index an arbitrary classical m-bit table F with
$`Q=2^r`$ rows, where $`m\ge1`$. There is an exact Clifford+T word

```math
|x,z,d\rangle\longmapsto|x,z\mathbin\oplus F(x),d\rangle
```

using no initialized work, with the conservative simultaneous bounds

```math
T,w=O\!\left(m(r+2)^3\sqrt Q\right),
\qquad
G=O\!\left(mQ+m(r+2)^3\sqrt Q\right),
\qquad
D_T=O(\log_2(r+2)).
```

Here w excludes the original address and output but includes every
simultaneous dirty helper. For $`Q=1`$, emit the constant row as an X
word; it has no T gates or helper requirement. The
[nonuniform indicator](NONUNIFORM_DIRTY_INDICATOR.md#1-the-exact-tree-interface-permits-unequal-chunks)
improves the count and width bounds to

```math
T,w=O(m\sqrt Q),\qquad G=O(mQ),
```

with the same depth. The conservative version suffices for the early
sector estimate in Section 4.

### One output bit

Split $`x=(a,b)`$ into address lengths with
$`H=2^{\lfloor r/2\rfloor}`$ and $`J=2^{\lceil r/2\rceil}`$.
For output bit j, allocate separate arbitrary dirty registers
$`Y_j\in\mathbb F_2^H`$ and $`X_j\in\mathbb F_2^J`$. Let
$`D_j\in\mathbb F_2^{H\times J}`$ be that bit's table matrix.
The retained [rank-reduced bilinear word](../../docs/PARALLEL_DIRTY_LOOKUP.md#exact-bilinear-oracle)
implements

```math
\mathcal B_j:\quad
z_j\longmapsto z_j\mathbin\oplus Y_j^{\mathsf T}D_jX_j.
```

Binary elimination and its actual inverse act only on the private Y
and X registers. In rank-normal coordinates, the middle word is a
common-target Toffoli batch, or after a target Hadamard a common-control
CCZ batch. The literal
[four-T-layer phase schedule](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
therefore gives

```math
T(\mathcal B_j)=O(\sqrt Q),\qquad
G(\mathcal B_j)=O(Q),\qquad D_T(\mathcal B_j)=O(1).
```

There is no additional work for this bilinear word. Let $`I_a,I_b`$
be exact dirty one-hot indicators on the two private registers. Execute
chronologically

```math
\mathcal B_j,\ I_a,\ \mathcal B_j^\dagger,\ I_b,\quad
\mathcal B_j,\ I_a^\dagger,\ \mathcal B_j^\dagger,\ I_b^\dagger.
```

The four contributions to the output, over $`\mathbb F_2`$, are

```math
Y^{\mathsf T}DX
+(Y+e_a)^{\mathsf T}DX
+(Y+e_a)^{\mathsf T}D(X+e_b)
+Y^{\mathsf T}D(X+e_b)
=e_a^{\mathsf T}De_b.
```

Every indicator register and private helper returns exactly. This is
the existing four-corner echo, with each dagger its actual native
inverse.

### Physical T-layer parallelism

Use the read-only conjunction indicator, whose public address occurs
only as a CNOT control. Its interface is stronger than merely preserving
the address at the end of a completed query.

The native reason is explicit in the retained
[controlled-increment proof](../../docs/PARALLEL_DIRTY_LOOKUP.md#a-read-only-controlled-increment-using-two-additions).
Let U and g be private arbitrary m-bit words, let
$`Q_g=\mathrm{ADD}(g;U)`$, and let E XOR a public literal ell into
the private low bit $`g_0`$. The chronological word
$`Q_g^\dagger,E,Q_g,E`$ returns g and adds
$`\ell(-1)^{g_0}`$ to U. Conjugation by CNOT fanout from $`g_0`$
to U removes that sign, giving $`U\mapsto U+\ell`$. All
non-Clifford gates are inside private modular addition. The public
literal occurs only in E as a CNOT control. A negative literal adds X
to the private target, never to the public address.

The subsequent sum compression and cyclic routing touch only private
counters and selector bits. The
[completed conjunction](../../docs/DIRTY_SUM_COMPRESSION.md#6-indicator-and-complete-frame-consequences)
has safe count and width majorants $`O((k+2)^3)`$ and depth
$`O(\log_2(k+2))`$ for k public literals. Its equality rows therefore
give an s-address indicator with

```math
T,G,w=O(2^s(s+2)^3),\qquad D_T=O(\log_2(s+2)).
```

For different output bits j, the registers $`Y_j,X_j,z_j`$ and
every private helper are disjoint. Every non-Clifford layer of their
single-output query schedules consequently has disjoint physical
supports. The shared public-address Clifford gates can be serialized
at the corresponding layer boundary. They never change the address
and act on distinct private targets. Padding shorter schedules with
identities synchronizes their T layers. Their elementary Clifford
counts are summed; no elementary-depth fanout claim is made.

This proves the m-output bounds. The sharper nonuniform indicator
retains exactly the same read-only interface and replaces its cubic
overhead by a constant. Its substitution proves the sharper count and
width statement without changing the parallelism argument.

### Activity and arbitrary outputs

For a preserved activity bit h, include h in the query address and
define $`F(0,x)=0`$, $`F(1,x)=F(x)`$. This doubles the table size
only. The completed query is exactly identity when h is zero, even if
all outputs and helpers are arbitrary. It consumes no initialized flag
or output bit. Interpreting an output as an unmasked program requires a
separate initialized-output promise from its caller.

## 2. A Boolean phase with arbitrary dirty work

Let $`Q_f`$ toggle an arbitrary dirty bit d by a read-only Boolean
function f. The chronological word

```math
Q_f,\ Z_d,\ Q_f^\dagger,\ Z_d
```

has phase

```math
(-1)^{d\mathbin\oplus f(x)}(-1)^d=(-1)^{f(x)}
```

and returns d and all query helpers. No supplied phase state is needed.
Use $`f(h,x)=hf(x)`$ for literal identity on the entire h-zero
sector. The phase identity holds on every computational input and hence
on arbitrary coherent inputs and external references.

This is an exact Boolean diagonal. It is not a supplied generic
non-Boolean phase or an arbitrary state-preparation oracle.

## 3. A Hadamard layer selected by a short label

Let $`t_0,\ldots,t_{a-1}`$ be distinct target qubits. Let lambda be
an unchanged short public label and h a separate unchanged activity
flag, both disjoint from every target. For arbitrary known Boolean
functions satisfying $`f_i(\lambda,0)=0`$, the desired word is

```math
\prod_{i=0}^{a-1}H_{t_i}^{f_i(\lambda,h)}.
```

Allocate arbitrary dirty bits $`d_0,\ldots,d_{a-1}`$. Let $`Q_F`$
be Section 1's parallel-output query, toggling each $`d_i`$ by
$`f_i(\lambda,h)`$. Define $`C=\prod_i\mathrm{c}H(d_i;t_i)`$
and execute chronologically

```math
Q_F,\ C,\ Q_F^\dagger,\ C.
```

The query preserves lambda and h; C preserves d. On every basis input,
the net action on each target is

```math
H^{d_i\mathbin\oplus f_i}H^{d_i}=H^{f_i}.
```

All d bits and query helpers return exactly. This complete-space
identity includes inputs entangling dirty helpers with the targets.
When h is zero the entire word is exactly identity.

For the native controlled Hadamard, use the literal identity already
charged in the [coarse-frame construction](../endpoint/STRUCTURAL_COMPILATION_LIMITS.md#exact-coarse-circuit-including-the-reused-logical-helper):

```math
B=SHTHS^\dagger,\qquad BZB^\dagger=H,
\qquad \mathrm{c}H(d;t)=B_t\,\mathrm{CZ}(d,t)\,B_t^\dagger.
```

The scalar phase of B cancels with its actual inverse. Each controlled
Hadamard has one T and one T-dagger, both on its target. Since the
targets are distinct, C has T-depth two; no shared public label is a
T target. The Clifford CZ gates remain counted.

If $`a=O(n)`$ and lambda contains $`O(\log(n+2))`$ bits, the query
has polynomially many rows and outputs. Its T-count, Clifford count
and dirty width are polynomial in n, while the complete Hadamard word
has

```math
D_T=O(\log_2\log_2(n+4)).
```

For example $`f_i(\lambda,h)=h[i\text{ is below the recorded depth}]`$
implements a coherent variable-suffix Hadamard layer. It does not need
n initialized selector bits. The recorded depth must be a separate
preserved label; if it initially resides in the target data, its
reversible extraction and return belong to the caller's ledger.

## 4. Exponentially smaller tables leave a polynomial budget

Put $`N=2^n`$, and suppose a terminal logical block of M bits is
reserved so that the early data dimension is $`Q_0=2^{n-M}`$.
Assume each query has at most $`c n^uQ_0`$ rows and $`c n^v`$
output bits, and there are at most $`c n^w`$ calls including actual
inverse calls. Extra polynomial-size labels count in $`n^u`$.
For fixed nonnegative exponents, the conservative query lemma gives

```math
T_{\rm total}
=O\!\left(n^{w+v+3+u/2}\sqrt N\,2^{-M/2}\right),
```

```math
G_{\rm total}
=O\!\left(n^{w+v+u}N2^{-M}
 +n^{w+v+3+u/2}\sqrt N\,2^{-M/2}\right),
```

```math
w_{\rm max}
=O\!\left(n^{v+3+u/2}\sqrt N\,2^{-M/2}\right).
```

Completed query pools are reused, so width is a maximum over calls.
If $`M/\log(n+2)\to\infty`$, these are respectively
$`o(\sqrt N)`$, $`o(N)`$, and $`o(\sqrt N)`$. In particular
$`M=\lceil\log_2^2(n+2)\rceil`$ suffices. Polynomially many additional
local gates and arbitrary dirty helpers are absorbed by the same
exponential budgets, but their serial T-depth is still charged.

For $`R=O(\sqrt n\log^c(n+2))`$ calls, the query depth is

```math
D_{T,\rm queries}=O\!\left(\sqrt n\log^{c+1}(n+2)\right)=o(n)
```

for every fixed c. This does not improve a schedule retaining
$`\Theta(n)`$ sequential fresh queries. Nor does it establish that a
classical data word can be interpreted as arbitrary prepared amplitudes:
the source isometry, its initialized work, guards, error and return must
be proved separately.

## 5. Proof scope

The query is a parallel composition of the retained exact dirty
conjunction and bilinear words. The private location of every T gate
is what permits the shared-address schedule. The phase and Hadamard
extensions are complete-space involution echoes. No clean-copy
assumption, measurement, reset, QRAM or unpriced quantum evaluator is
introduced. Classical tables are compiled data; no efficient classical
preprocessing theorem is asserted here.

These primitives neither prove a complete-frame depth bound alone nor
an unrestricted depth lower bound. Their use in a whole compiler must
also price the source construction, logical packing, recovery and any
late-layer tail on the same circuit.
