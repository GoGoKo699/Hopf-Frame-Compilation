# Conditional suffix workspace improves the two-clean endpoint

[Open endpoint](OPEN_PROBLEM.md) · [Operator-source compiler](OPERATOR_SOURCE_COMPILER.md) · [Grouped residuals](FAULT_TOLERANT_COMPILER.md#6-grouped-dictionaries-and-an-exactly-clean-coarse-frame)

The logical suffix of a Hopf group is zero on the sector where that group
acts. It can therefore supply temporary initialized workspace, provided the
actual circuit is identity on every inactive sector. This observation permits
several consecutive layers to share one operator-source precision charge.

**Theorem.** Let $`n\geq3`$, $`N=2^n`$, and let the accuracy parameter satisfy
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}=N`$. At the allocation

```math
a=2,\qquad b=N+n+7,
```

every prescribed complete real Hopf frame has a coherent Clifford+T compiler
satisfying

```math
\|VJ_2-J_2(W\otimes I_b)\|\leq\eta,
\qquad
T=O\bigl(N\log(n+2)\bigr),\qquad G=O(N^2).
```

The error includes both initialized-work leakage and the return of arbitrary
dirty inputs, together with any references. The construction uses no
measurements, resets, supplied resource states, or uncharged quantum oracles.
It is an asymptotic improvement from $`O(N\log N)`$ to
$`O(N\log\log N)`$ at this endpoint; it does not prove $`O(N)`$.
The fixed thresholds in the proof need not produce a practical saving for
small systems. Classical coefficient evaluation and table construction remain
separate preprocessing costs, as in the other compiler theorems. Angles must
be effectively specified so that the sine, cosine, and finite residual-matrix
entries used below admit certified rational evaluation.

## 1. An active suffix supplies workspace conditionally

Consider a group of depths $`[e-s,e)`$, and write

```math
r=n-e,\qquad p=e-s=n-r-s.
```

The first $`p`$ logical bits are a preserved prefix, the next $`s`$ bits are
the local system, and the final $`r`$ bits are the suffix. The ideal group is
identity unless the entire suffix is zero. On that active sector it is a
direct sum, indexed by the prefix, of complete $`s`$-level frames.

Compute $`h=[\mathrm{suffix}=0]`$ into the first external clean qubit using
the exact borrowed-MCX construction. One local system bit can be borrowed
and returned for this computation; its value need not be known. The cost is
$`O(r^2)`$ Toffolis. On the sector $`h=1`$, the $`r`$ suffix bits are now
available as zero-initialized private work. On $`h=0`$ they are arbitrary.

Every actual group subroutine below preserves $`h`$ and is identity on
$`h=0`$, on the whole input space of its other registers. Its active action
returns the suffix to zero up to the stated group error. Finally, reverse
the actual predicate computation. The ideal output has the same suffix
predicate as the input, so this last unitary returns $`h`$ to zero without
increasing the error. This reasoning does not assume that the actual suffix
is exactly zero after an approximate subroutine.

### Exact conditioning without a Clifford-count penalty

Conditioning a whole-word XOR query is done by adding $`h`$ to its address
and giving every inactive row the zero word. Its T-count increases by a
constant factor; its many leaf CNOTs are still Clifford gates. The same
rule applies when there is an additional identity/correction mode bit.

A controlled operator source has the form

```math
U_m\,C_F(X_0)\,U_m^\dagger,
```

where $`F`$ is the conjunction of a bounded number of flag literals. The
central multiple-controlled X costs a constant number of T gates, and
the unconditional source basis changes cancel exactly when $`F=0`$.

Other native gates needing a bounded number of controls also have constant
exact cost with a constant number of returned dirty helpers. For completeness,
the controlled-phase case does not require assuming that inactive suffix
scratch is clean. If $`z`$ is a dirty bit and $`F`$ is a Boolean conjunction
not depending on it, the chronological phase echo

```math
X_F^{(z)},\ T_z,\ X_F^{(z)},\ T_z^\dagger
```

returns $`z`$ and contributes the phase
$`\exp(i\pi(F-2zF)/4)`$. Correct its residual factor by the phase
$`i^{zF}`$. A second dirty echo using S implements that phase up to
$`(-1)^{wzF}`$, and this last correction is an exact multiple-controlled Z.
Only a bounded number of controls and returned dirty helpers occur. An
S-phase uses the last two steps of the same construction. A controlled H
uses $`H=R_y(\pi/8)ZR_y(-\pi/8)`$: the two fixed rotations have constant
Clifford+T cost up to opposite common scalars, which cancel in this
conjugation. They also cancel when the central controlled Z is inactive.
X, Z, and CNOT controls reduce to exact multiple-controlled X or Z.
All such helpers are disjoint from the controls and are returned before
the next operator-source call. A fixed number of core or selector wires
can serve this purpose; the ledger below also allows a separate fixed
helper reservation.

The same conditioning method applies to fixed coarse program words.
It applies to complete subroutines; their individual native gates need not
be identity on the inactive sector.

## 2. A coarse group and its local residual dictionary

Put

```math
D_s=(2s-1)2^s+2,\qquad
K=2^{\lceil\log_2(4D_s)\rceil}.
```

The [sparse grouped-residual lemma](FAULT_TOLERANT_COMPILER.md#6-grouped-dictionaries-and-an-exactly-clean-coarse-frame)
gives a permitted support of $`D_s`$ ordered column pairs for a local
$`s`$-level residual. This remains true for complex coarse Clifford+T
words: their marker columns have the same nested dyadic subtree supports.
There is no assumption that a coarse word is a real rotation exactly.

Choose phase-calibrated actual Clifford+T words for every group angle with
one-qubit error at most

```math
\delta=\frac{1}{4sK}.
```

Their lengths are $`w=O(\log(1/\delta))=O(s)`$. Execute them in the
prescribed addressed order, with exact identity on the inactive local
suffix sectors. Denote this actual coarse group by $`C`$, its target by
$`W_g`$, and its active residual by

```math
E=C^\dagger W_g-I,\qquad \|E\|\leq s\delta\leq\frac{1}{4K}.
```

The clean program and predicate work returns exactly after C. The words
are fixed classical choices before computing the residual tables. The
definition of E calls no quantum implementation of $`W_g`$.

For each logical prefix $`x`$, split the real and imaginary parts of every
permitted entry into nonnegative coefficients and phases from
$`\{1,-1,i,-i\}`$. Pad to K terms:

```math
E_x=\sum_{\ell=0}^{K-1}a_{x\ell}A_\ell,
\qquad
A_\ell=\omega_\ell|u_\ell\rangle\langle v_\ell|,
\qquad 0\leq a_{x\ell}\leq\frac{1}{4K}.
```

The endpoints are local s-bit words. Their permitted values and their
phase labels depend only on $`\ell`$; the numerical coefficient can
depend on x. Set $`c_{x\ell}=K a_{x\ell}\in[0,1/4]`$. Coefficients are
obtained by certified rational enclosures and the operator-source dyadic
rounding rule. Clipping intervals at zero supplies positive and negative
parts without an undecidable exact sign test. Structural zero entries
and padding receive literal zero encodings.

The addressed coarse interpreter uses an $`O(s)`$-bit program word and
$`O(s)`$ predicate work. At each of the s local depths its table is
addressed by the unchanged logical prefix and the preceding local bits.
Its total table-query T-count is $`O(2^e)`$; word interpretation and exact
predicates contribute a polynomial in s. Whole-word table Clifford cost
is $`O(s2^e)`$. All operations are conditioned on h as above. These are
upper bounds for a literal interpreter, with every actual word and inverse
charged, rather than a use of the coarse target as an oracle.

## 3. One scalar source encodes the entire group residual

Use the second external clean qubit as a scalar flag. Inside the active
suffix reserve an identity/correction mode f, a $`\log_2 K`$-bit label,
two s-bit endpoint words, a separate matrix-unit flag, and predicate work.
The coarse program word can reuse this workspace before or after the
residual block. There is a fixed constant $`C_0`$ for which all private
initialized work other than the two external qubits occupies at most
$`C_0s`$ wires, including $`s=1`$.

The operator-source scalar unitary $`\mathcal S_c`$ has, on its own
zero flag, the accepted action

```math
\widehat c_{x\ell} I_{\rm dirty},\qquad
|\widehat c_{x\ell}-c_{x\ell}|
\leq\frac54\,2^{1-m}.
```

Its dirty core has m bits. The whole-word coefficient query uses x and
$`\ell`$ as its address; its output is the arbitrary core itself. No
initialized m-bit coefficient word is used.

The matrix-unit unitary first toggles its separate flag by
$`[z\ne v_\ell]`$, then XORs $`u_\ell\oplus v_\ell`$ into the local
logical word z, and supplies the literal phase $`\omega_\ell`$. Its
zero-flag block is exactly $`A_\ell`$. Equality scratch is erased before
changing z. The endpoints and label are unchanged and hence unload
exactly, including on rejected branches.

Prepare f and the label uniformly with Hadamards, load the endpoints,
and use the following SELECT:

| Mode | Actual selected unitary |
|---|---|
| $`f=0`$ | Identity |
| $`f=1`$ | Scalar unitary $`\mathcal S_c`$, followed by the matrix-unit unitary |

Then unload metadata and undo the Hadamards. Call this actual unitary Q
on the active sector. Condition its complete implementation on h.
The source and matrix-unit flags are distinct. The scalar unitary leaves
the local logical word unchanged; the matrix-unit operation leaves the
prefix and label unchanged. Thus projecting both flags factors their
accepted actions exactly, with no assumed removal of a rejected branch.

Let J initialize the private suffix work and the scalar flag, leaving
the logical prefix, local system, and every dirty input arbitrary. On
$`h=1`$ the accepted block is

```math
B=J^\dagger QJ
=\frac12\left(I+\frac1K\sum_\ell
\widehat c_{x\ell}A_\ell\right)\otimes I_{\rm dirty}.
```

Consequently

```math
\|2B-(C^\dagger W_g\otimes I_{\rm dirty})\|
\leq\frac1K\sum_\ell
|\widehat c_{x\ell}-c_{x\ell}|\,\|A_\ell\|
\leq\frac54\,2^{1-m}.
```

The estimate holds uniformly across the coherent prefix direct sum.
The factor K in the coefficient definition cancels the uniform-label
average; it does not amplify the rounding error.

## 4. Amplification, inactive phases, and the complete group contract

Let $`R_h`$ be identity for $`h=0`$ and, for $`h=1`$, the reflection
$`I-2JJ^\dagger`$. It may test the entire r-bit suffix together with
the external scalar flag; unused active suffix bits remain zero. This
conditional reflection costs $`O(r^2)`$ Toffolis using a returned borrowed
local bit. Define

```math
\mathcal A_h=Z_h Q R_h Q^\dagger R_h Q.
```

Here Q and its inverse are their h-conditioned actual circuits. When
$`h=0`$, every factor is identity. When $`h=1`$, $`Z_h=-1`$ supplies
the literal leading minus sign of normalization-two oblivious amplitude
amplification. In particular, no spurious minus sign is applied to the
inactive logical sector.

The [full-isometry amplification estimate](OPERATOR_SOURCE_COMPILER.md#5-amplification-includes-rejected-space-error)
gives, including every rejected component,

```math
\|\mathcal A_hJ-J(C^\dagger W_g\otimes I_{\rm dirty})\|
\leq 5\,2^{1-m}=10\,2^{-m}
\qquad(h=1).
```

Now append the actual conditioned coarse circuit C, with its exact
clean-work return, and erase h as in Section 1. Unitary propagation of
the above error proves the group contract for $`W_g`$, with the same
bound $`10\,2^{-m}`$. This argument applies C after completed amplification;
it does not multiply accepted blocks while ignoring coherent rejected
components. It also permits arbitrary joint states of the dirty core and
all query selectors at the start of each group.

## 5. Workspace and costs of one group

Choose a fixed grouping constant $`C_1\geq2C_0`$, enlarging it if needed
to cover all fixed workspace conventions. A group with
$`1\leq s\leq r/C_1`$ fits its private initialized work into the r-bit
active suffix. Take

```math
m=N+\lfloor r/4\rfloor+8.
```

The coefficient query has at most

```math
k=p+\log_2 K+2=n-r+O(\log(s+1))
```

address bits, including h and f. Its dirty selectors, the m-bit core,
and a fixed number of separate returned helpers therefore require

```math
m+k+O(1)
=N+n-\frac34r+O(\log(r+1))+O(1)
\leq N+n+7
```

for all $`r\geq r_0`$, where $`r_0`$ is a sufficiently large fixed
constant. Shorter metadata queries reuse the same selector pool. The
clean label and endpoint registers are part of the logical suffix, not
an extra external clean allocation. Every temporary dirty helper is
returned before another use that assumes it is available.

Writing $`Q_g=2^pK=O(s2^{n-r})`$, a constant number of source and query
calls, three amplification calls, the exact coarse interpreter, and all
predicates give

```math
T_g=O(Q_g+m+\mathrm{poly}(n)),\qquad
G_g=O(Q_gm+s2^{n-r}+m+\mathrm{poly}(n)).
```

There are only a fixed number of control literals in the gates other
than predicates and whole-word queries. Hence conditioning does not
replace the $`Q_gm`$ Clifford term by a T-count of that size.

## 6. Partition, error budget, and summation

Choose the fixed $`r_0`$ large enough for the workspace conditions,
$`r_0\geq2C_1`$, and the error inequality below. If $`n\leq r_0`$,
use the original two-clean operator-source compiler; its bounds already
have the asserted asymptotic form because $`r_0`$ is constant.

Otherwise reserve the deepest $`r_0`$ individual layers and set
$`r=r_0`$. Until all depths are assigned, form a group ending at
$`e=n-r`$ with

```math
s=\min\{n-r,\lfloor r/C_1\rfloor\},\qquad r\leftarrow r+s.
```

This discovers groups from deep to shallow. Execute them in reverse
discovery order, followed by the reserved individual layers in their
original shallow-to-deep order. The complete circuit therefore has
exactly the prescribed Hopf layer order. Every group contract is on the
entire logical space, so suffix workspace is recomputed at each group;
it is not assumed to stay zero between different groups.

For a reserved layer at depth $`d=n-j`$, $`1\leq j\leq r_0`$, use its
original allocation $`m_d=N+j+4`$. The sharper bound proved before the
rounding to $`2^{4-m_d}`$ in the operator-source chapter is

```math
5\sqrt2\,2^{1-m_d}
=\frac{5\sqrt2}{8}\,2^{-N-j}.
```

Their sum is strictly below $`(5\sqrt2/8)2^{-N}`$. No extra source wire
or uncharged precision bit is introduced in this estimate.

For the grouped stages, r is distinct at each stage, so even summing over
all integers $`r\geq r_0`$ gives

```math
\begin{aligned}
\sum_g10\,2^{-m_g}
&\leq10\,2^{-N-8}
\sum_{r=r_0}^{\infty}2^{-\lfloor r/4\rfloor}\\
&\leq80\,2^{-N-8-\lfloor r_0/4\rfloor}.
\end{aligned}
```

Choose $`r_0`$ such that

```math
\frac{5\sqrt2}{8}
+80\,2^{-8-\lfloor r_0/4\rfloor}\lt1.
```

The complete-isometry hybrid then gives total error below
$`2^{-N}\leq\eta`$. Earlier actual leakage is propagated unitarily;
only ideal comparison stages are required to return their work exactly.

Except possibly at the last group, r grows by a fixed factor greater
than one. Therefore the number of groups is $`O(\log(n+2))`$ and
$`\sum_g r=O(n)`$. Moreover,

```math
\sum_g Q_g
\leq O(N)\sum_{r\geq r_0}r2^{-r}=O(N),
\qquad
\sum_gm_g=O\bigl(N\log(n+2)+n\bigr).
```

All fixed-degree polynomial overheads in n, including the group count,
are $`O(N)`$. The fixed number of reserved layers costs $`O(N)`$ T
gates and $`O(N^2)`$ Clifford gates. Since $`m_g=O(N)`$ and
$`\sum_gQ_g=O(N)`$, the grouped Clifford cost is $`O(N^2)`$ as well.
This proves the theorem.

The construction uses initialized dimension supplied by a verified active
logical sector and acts as identity on its complement. It does not supply
a precision-sized initialized source on every arbitrary logical input,
and therefore does not contradict the nilpotent carried-source restriction
in [the source-reuse note](SOURCE_REUSE_LIMITS.md). The remaining endpoint
gap is between $`\Omega(N)`$ and $`O(N\log\log N)`$.

The [focused finite checks](../tests/test_conditional_suffix_compiler.py)
exercise the native dirty phase echo, distinct scalar and matrix-unit
flags, literal complex phases, full dirty-input amplification and leakage,
coherent inactive suffixes and predicate erasure, and the sparse support
for small exact complex coarse frames. Illustrative integer ledgers check
the allocation formulas. These fixtures do not certify a numerical value
for the coarse-synthesis workspace constant or replace the asymptotic
construction and error proof above.
