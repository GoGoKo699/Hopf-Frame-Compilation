# Native Clifford transport and nonlinear returned characters

[Current frontier](../../docs/OPEN_PROBLEM.md) ·
[Transport proof](NATIVE_CLIFFORD_TRANSPORT.md) ·
[Nonlinear character proof](NONLINEAR_CHARACTER_TRANSPORT.md) ·
[Previous checkpoint](SIX_LAYER_TRANSPORT_SESSION_RESULT.md)

This round began on 18 October 2026 at verified main
`891c76717ca1d693641ace689ef509794d8fd072` (PR #174).
It addressed the missing signed-degree transport through the six-layer
comparison's shared three-band prefix, alongside a nonlinear constructive
escape from the fixed-affine diagonal-source obstruction. One new scoped
six-layer theorem and a fully native, slow return mechanism are supported.
Neither improves the unconditional complete-frame depth law.

## 1. Compact revision

The count package is settled under its stated hypotheses. R96 remains the
strongest unconditional complete-frame upper law, simultaneously

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+L+\log_2(n+2)),\qquad N=2^n,
```

for every integer $`n\ge1,L\ge6`$, full error at most $`2^{-L}`$,
one clean qubit, and real width $`b\ge17(L+n+7)`$ or complex magnitude
width $`b\ge17(L+n+8)`$. Every logical column, literal sign and phase,
singular and zero coefficient, arbitrary dirty input and external reference,
actual inverse, leakage and complete workspace return remains required.
R79 remains selected Result D. R97 is an outer-overlap obstruction;
unrestricted one-clean four-layer rigidity remains open.

The exact general question is optimal worst-case T-depth for the prescribed
complete frame, at the same precision and simultaneous count and workspace
budgets. R98 conditionally transfers a faster exact multi-output read-only
query to fixed-accuracy frame depth. The separate two-call corner reduction
transfers unrestricted query lower bounds. Their hypotheses are different,
and neither missing general premise has been proved.

Already resolved questions include optimal one-clean count, the all-precision
upper law, zero-clean read-only three-layer affine and four-layer quadratic
XOR rigidity, and the general five-layer degree bound 14 below
$`\sqrt2/24`$ (11 for invisible outer shear). The six-layer
Pauli-comparison bound 51 below $`\sqrt2/64`$ is a separate restricted
case. Historical next-step passages do not reopen these results.

## 2. What closes in the transport argument

The [new proof home](NATIVE_CLIFFORD_TRANSPORT.md) separates two questions.
For any third-level unitary F and Cliffords C,K,

```math
d_{\rm pr}(FCF^\dagger,K)\lt\sqrt2/4
\quad\Longrightarrow\quad FCF^\dagger\doteq K.
```

Conjugated Hermitian Paulis are Hermitian Cliffords; their uniform gap
rounds the endpoint exactly. The proof does not assume the inverse of F is
third-level, or that the intermediate band conjugations are Clifford.
This establishes exact membership, not a degree bound.

For the additional source form $`F=C_{\rm out}D C_{\rm in}`$ with
fixed Cliffords and a diagonal third-level phase D of the weighted
modulo-eight cubic form, a signed
quotient argument also controls degree. A fixed diagonal Clifford matches
D literally on its Pauli translation directions. The remaining phase
depends only on quotient observables. Exact endpoint membership forces
these observables to remain diagonal under the inner Clifford, even when
that Clifford is nonmonomial. Their native label degree is two and their
literal sign degree six. Modulo-eight carries require a factor four,
not three. The resulting actual comparison has symplectic degree at most
26 and fixed-probe sign degree at most 78.

Thus the specified-source **zero-clean read-only six-layer** query has
degree at most **78**, or **26** for visible outer target shear, below
full operator error $`\sqrt2/8`$, independently of dirty width. The
fixed coordinates and quotient are analytical; no free physical source or
helper is introduced. All original gates, inverses and allocations remain
charged. Sharpness is not claimed. General three-band prefixes need not
have this source form.

The tempting extension from third-level membership alone fails. The same
home gives a self-contained signed-transport calculation for the known
He–Robitaille–Tan permutation family: a linear Pauli family acquires
symplectic degree $`r-1`$ and sign degree r on $`2^r-1`$ work wires.
Its explicit native source and actual inverse are charged. This is an
endpoint-transport stress test, not an actual three-band returned query.
The native-depth restriction must supply information that the hierarchy
level by itself does not contain.

## 3. Constructive reformulation and its exact limitation

The [nonlinear character proof](NONLINEAR_CHARACTER_TRANSPORT.md) uses
all w-bit strings as units $`1+t\mathbb F_2[t]/(t^{w+1})`$. An
address-dependent unit with affine coefficient functions multiplies the
dirty word. A Boolean group character gives a diagonal phase D and

```math
D P_xD^\dagger P_x^\dagger=(-1)^{\chi(a_x)}I.
```

The action is exact on every dirty input. The multiplication is nonlinear
jointly in address and work, so this genuinely leaves the earlier
fixed-affine transport hypothesis. All characters are controlled by odd
logarithmic-derivative coefficients, giving

```math
\deg\chi(a_x)\le1+\lfloor\log_2 w\rfloor.
```

An explicit family has unbounded degree, so a fixed number of abstract
diagonal occurrences no longer suffices for a uniform degree obstruction.
But a degree-a output needs $`w\ge2^{a-1}`$, whereas the desired
one-output R98 reader uses $`B=\Theta(\sqrt{2^a})`$ dirty work.
This width mismatch stops this representation even if its source costs
were free; it is not a general native-query lower bound.

The proof supplies native factorization, multiplication, actual inverses,
and controlled phases for XOR conversion. The sequential construction
uses $`O(w^2)`$ T-count and T-depth,
$`O((a+1)w^2+qw)`$ elementary Cliffords for q outputs, w coefficient
wires and one reused arbitrary helper, with no clean work. It is exact,
returns the full workspace, and needs no eigenstate or code-space promise.
R91 already handles arbitrary tables at the smaller target width; this
construction gives no resource improvement. Any approximation comparison
charges each source occurrence and inverse separately.

## 4. Resource correction and evidence

The [query/frame reduction](BOOLEAN_QUERY_FRAME_REDUCTION.md) now keeps
the actual inverse's elementary Clifford count explicit. In the fixed
H, S, CNOT alphabet, $`S^\dagger=S^3`$, so its frame-to-query ledger is

```math
T'\le2T,\qquad D_T'\le2D_T,\qquad
G'\le G(V)+G(V^\dagger)+2\le4G(V)+2.
```

The converse uses two forward calls and retains $`2G+2`$. The scalar
conversion also separately charges its actual inverse. No asymptotic
transfer, error bound or workspace allocation changes.

Focused evidence checks the signed quotient extension, nonmonomial inner
action, modulo-eight carries, actual endpoint cancellation and the C3-only
scope counterexample. Separate nonlinear fixtures check complete dirty
return, all-character degree, source factorization and controlled native
conversion. The inverse-ledger fixture detects the lost S-dagger charge.
These finite checks support the identities and implementation ledgers;
the proofs establish arbitrary-width statements. Proof homes, frontdoors
and the verification catalogue keep each claim at its actual scope.

## 5. Precise unresolved step and resumption question

The most useful next question is **signed-degree transport through the
actual three-band prefix beyond the fixed diagonal source class**. The
inner Clifford is a product of two one-band Pauli images and a linear
Pauli; its degree is controlled before conjugation. The missing invariant
must use these native factors, the three physical bands and full-return
compatibility, rather than arbitrary C3 membership or an abstract
normalizer. A counterexample must realize the actual returned query and
its specified comparison, with every source paid. The known high-degree
permutation source does not yet do that.

The earlier odd-order returned fixture still blocks a universal repeated
Pauli-return recurrence. Two other possible restrictions were considered
but are not promoted: the proposed robust group-richness argument needs
a uniformly short word cover
of the entire Pauli frame, and localization to a fixed active subsystem
needs a source-uniform separation gap. Existence of words or a positive
gap at each fixed width does not provide either quantitative condition.

**Resumption question:** can independent physical axes and compatibility
of the two inner one-band factors force a width-independent bound on the
actual transported labels and literal signs, including odd-order modes,
or can a fully returned six-layer family violate such a bound?
Success would settle a genuinely open constant-layer step and determine
whether this route can support a depth recurrence. Error as depth grows,
one-clean work and unrestricted address mixing remain additional gaps
before a growing general-frame lower bound follows. Constructively, a
different nonlinear representation would need to achieve degree-a return
within the square-root-width pool and beat the fully simplified R91 ledger.
Optimal general depth remains the destination; final manuscript drafting
is separate.
