# State preparation with one precision charge and two clean flags

[Research status](OPEN_PROBLEM.md) · [Borrowed compiler](BORROWED_WORKSPACE_COMPILER.md) · [Full-operator rotation primitive](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)

A cheap native coarse circuit makes a known Hopf state close to the
computational root. Its remaining amplitudes can then be loaded by one
addressed two-by-two table, with constant normalization and one exact
state-amplification step. The resulting compiler prepares the requested
state with all work returned within the error. It does not implement the
prescribed complete Hopf frame on arbitrary logical inputs.

## 1. State-only contract

Let $`N=2^n`$, $`n\ge1`$, and let the normalized real state
$`|\psi\rangle=W(\boldsymbol\theta)|0^n\rangle`$ be specified by
a real Hopf angle tuple. Set

```math
0\lt\eta\le1/64,\qquad
L=\max\{6,\lceil\log_2(1/\eta)\rceil\},\qquad L\ge n.
```

There is a coherent Clifford+T circuit V with two initialized compiler
flags and $`b\ge L+n+7`$ arbitrary dirty qubits such that

```math
\sup_{\|\xi\|=1}
\left\|V\bigl(|0^n\rangle|00\rangle|\xi\rangle\bigr)
 -|\psi\rangle|00\rangle|\xi\rangle\right\|\le\eta,
```

```math
T(V)=O(N+L),\qquad G(V)=O(NL).
```

Equivalently, this is an operator-norm bound between isometries from the
entire dirty-input space; it includes arbitrary external references.
Literal phase, clean-flag leakage, and dirty-work return are included.
The n logical system qubits must start in zero. No action on other
logical inputs is promised. At $`L=N`$ this gives linear T-count for
state preparation, not a new complete-frame frontier.

The construction also prepares either of two supplied real Hopf states
coherently under an unchanged protocol branch, with the same compiler
allocation and asymptotic counts; Section 6 specifies that contract.
Section 7 gives the separate common-coarse corollary when one reference
state is the generally complex output of the actual native coarse word.
Classical coefficient evaluation and certified table generation are
separate preprocessing. No efficient bound on their bit complexity is
asserted, and no exact-zero test for arbitrary computable amplitudes is
assumed.

## 2. A native coarse circuit makes the tail small

Use the [borrowed-workspace theorem](BORROWED_WORKSPACE_COMPILER.md#1-contract-and-statements)
for the supplied real frame at error

```math
\eta_c=\min\{1/64,1/(4\sqrt N)\},\qquad
L_c=\max\{6,\lceil\log_2(4\sqrt N)\rceil\}=O(n).
```

Leave the two compiler flags untouched. The actual circuit implements
$`C\otimes I_b`$ exactly, with
$`\|C-W\|\le\eta_c`$ and all dirty helpers returned exactly.
Consequently

```math
T(C)=O\!\left(\frac{NL_c}{n+b}+L_c\sqrt N\right)=O(N),
\qquad G(C)=O(Nn)=O(NL).
```

Here $`b\ge L+n+7`$, $`L\ge n`$, and
$`n\sqrt N=O(N)`$. The logical C may be complex even though the
requested state is real. Define the exact complex residual state

```math
|\phi\rangle=C^\dagger|\psi\rangle
=a|0^n\rangle+|v\rangle,\qquad v_0=0.
```

It is normalized and obeys

```math
\|\phi-|0^n\rangle\|\le\frac1{4\sqrt N},\qquad
\|v\|\le\frac1{4\sqrt N},\qquad |a|\le1.
```

The coarse error is used only to bound the tail. It is not added to the
final approximation error: $`C|\phi\rangle=|\psi\rangle`$
exactly for this fixed native C.

The logical C is the recorded tree of actual two-mode native words.
Once their coefficients and those of psi are available, applying its
inverse to the amplitude vector takes $`O(N)`$ scalar arithmetic
operations. This does not simulate the Hilbert space of the dirty
implementation. The argument also works for a supplied state and supplied
native C satisfying this same closeness and exact-return promise; the
stated theorem obtains C from the given Hopf tuple.

## 3. Two flags give an exact half-amplitude state

Call the clean flags s and t. For every complex z in the unit disk put

```math
U(z)=\begin{pmatrix}z&-\sqrt{1-|z|^2}\\
                    \sqrt{1-|z|^2}&\overline z\end{pmatrix}
\in\mathrm{SU}(2).
```

Use the addressed table on target t with unchanged address $`(s,x)`$:

```math
z_{0,x}=a,\qquad z_{1,0}=0,\qquad
z_{1,x}=\sqrt N\,v_x\quad(x\ne0).
```

The small-tail estimate gives $`|z_{1,x}|\le1/4`$. Let M apply
$`U(z_{s,x})`$ to t, let $`K=C_{s=1}(H^{\otimes n})`$, and set

```math
Q=H_s M K H_s.
```

The rightmost gate acts first. With
$`P_{\rm good}=I_{\rm sys}\otimes|00\rangle\langle00|_{s,t}`$,

```math
P_{\rm good}Q|0^n,00\rangle
=\frac12\left(a|0^n\rangle+
              \frac1{\sqrt N}\sum_x z_{1,x}|x\rangle\right)|00\rangle
=\frac12|\phi,00\rangle.
```

The ideal Q is identity on dirty work. This equation concerns every
dirty/reference input, not a postselected preparation procedure. In
particular, an active zero coefficient uses $`U(0)=R_y(\pi/2)`$,
not the identity; rejected amplitudes are retained.

## 4. One state-amplification step, with actual inverses

Let J append the initialized system and two flags to arbitrary dirty
input. Define

```math
R_{\rm init}=I-2JJ^\dagger,\qquad
R_{\rm good}=I-2P_{\rm good},\qquad
\mathcal A=-Q R_{\rm init}Q^\dagger R_{\rm good}Q.
```

The initial reflection tests the system and both flags; the good
reflection tests only the two flags. Both act identically on dirty work.
Writing G for the isometry appending $`|\phi,00\rangle`$, the exact
half-amplitude identity gives orthogonal isometries G and B with

```math
QJ=\tfrac12G+\tfrac{\sqrt3}{2}B,\qquad
P_{\rm good}B=0.
```

Therefore $`(QJ)^\dagger R_{\rm good}QJ=I/2`$, and direct
substitution proves $`\mathcal A J=G`$. This is the single Grover
step with angle $`\pi/6`$, including its literal leading minus sign.
The amplification law is inherited from
[Brassard–Høyer–Mosca–Tapp, Section 2, Eq. (8)](https://arxiv.org/pdf/quant-ph/0005055);
the calculation here specifies its dirty-input and phase contract.
That sign is native, for example $`XZXZ=-I`$ on either flag.

For an actual native half-amplitude word $`\widehat Q`$ with
$`\|\widehat Q-Q\|\le\delta_Q`$, use its actual circuit inverse
in $`\widehat{\mathcal A}`$. Unitary telescoping gives

```math
\|\widehat{\mathcal A}J-G\|\le3\delta_Q.
```

No intermediate work is reset. Applying the actual coarse C afterward
preserves the norm bound and changes the ideal state to psi. The initial
reflection has a native $`O(n^2)`$ implementation with one returned
arbitrary helper; the good reflection and leading minus sign cost a
constant number of Cliffords. The controlled Hadamards in K each have a
constant-size exact Clifford+T implementation, so K costs $`O(n)`$.

## 5. Native table, classical rounding, and the literal dirty budget

The [certified SU(2) Euler procedure](ONE_CLEAN_COMPILER.md#7-literal-diagonals-and-complete-one-target-multiplexors)
approximates each U(z) by three real-angle Euler rotations. It need not
choose an exact argument at z zero. Synthesize each addressed Ry or Rz
factor using the full-operator borrowed-signal primitive; fixed target
Cliffords convert between its axes. The flags s and t may be occupied
throughout: the synthesis signal is a separate arbitrary dirty wire.

Set $`q=L+10`$. Choose the classical Euler approximation to each
U(z) within $`2^{-q}`$. This includes coefficient rounding and the
completion square root. For example, round and project z into the unit
disk within $`e=2^{-2q-12}`$; then

```math
\|U(z)-U(\widehat z)\|
\le e+\sqrt{2e}.
```

A further certified Euler approximation within $`2^{-q-2}`$ leaves
the total below $`2^{-q}`$. The factor $`\sqrt N`$ in the tail
table requires correspondingly finer component enclosures, by another
$`n/2`$ bits. Near $`|a|=1`$ the displayed stronger classical
accuracy handles the square root; it does not enlarge the quantum
source width. Certified matrix-entry search avoids singular Euler-chart
or exact-zero assumptions.

Each native rotation has full-operator error below $`43\,2^{-q}`$.
The three factors therefore give

```math
\|\widehat M-M\|\le130\,2^{-q},\qquad
\delta_Q\le130\,2^{-q},\qquad
3\delta_Q\le390\,2^{-q}\lt2^{-L}\le\eta.
```

For $`n\ge6`$, fix seven of the $`n+1`$ address bits as sector
literals, leaving $`k=n-6`$ free address bits. There are 128 sectors,
a fixed constant. Each native rotation uses

```math
(q+1)_{\rm core}+k_{\rm selectors}
 +1_{\rm predicate\ helper}+1_{\rm signal}=L+n+7
```

arbitrary dirty wires. All address bits are unchanged during that
rotation, and every inactive sector is exactly identity on the complete
work space. Thus the errors over sectors take a maximum, not a sum.
The bounded-control sector predicates retain their native charges.
Queries finish and unload before the final $`H_s`$ changes the mode.
Reflection helpers reuse returned source work only between complete
subroutines, and their exact identities hold even on leaked work.

There are three Euler factors and three Q appearances. Summing their
constant-sector counts gives $`T=O(N+L)`$ and $`G=O(NL)`$.
The exact reflections, coarse word, and fixed gates obey those bounds
as well. No initialized precision word or history register is used.
For $`n\le5`$, use the retained direct borrowed native words for the
finitely many Hopf rotations instead. This fixed-depth fallback gives
the same asymptotic counts and fits the stated allocation, without a
negative free-address length.

## 6. Two states under one protocol branch

Let an unchanged protocol qubit d select supplied real Hopf states
$`|\psi_0\rangle,|\psi_1\rangle`$. The preparation circuit satisfies
the coherent isometry contract

```math
\left\|VJ_d-\sum_{j=0}^1
 |j\rangle\langle j|_d\otimes|\psi_j,00\rangle\otimes I_b
\right\|\le\eta,
```

where $`J_d`$ appends the initialized n-qubit system and compiler flags
to arbitrary branch-and-dirty input. This formula preserves the relative
phase of the branches, including arbitrary dirty/reference correlations.
The protocol qubit is not an additional initialized compiler flag.

Compile an actual coarse direct sum $`C_0\oplus C_1`$ by the borrowed
construction on an $`n+1`$-qubit tree whose top gate is **literally
omitted**, and whose two subtrees have the two supplied angle tuples.
Every remaining local word preserves d. This prices a block-diagonal
native program; it does not assume that controlling an arbitrary already
compiled C is free. Use coarse error
$`\min\{1/64,1/(4\sqrt N)\}`$ as before.
Doubling the table size changes only constants in $`T=O(N)`$ and
$`G=O(Nn)`$, with exact dirty return.

Define $`\phi_j=C_j^\dagger\psi_j`$ and use the same two-flag
table with address $`(s,d,x)`$. Both branches have exact success
amplitude $`1/2`$. The initial reflection excludes d, so the same
amplification proof holds coherently. Fix eight address bits rather than
seven: there are 256 sectors and again $`k=n-6`$, giving precisely
the same dirty reservation and error bound. In the fixed-size fallback,
the arbitrary-budget borrowed compiler implements the two subtrees while
leaving the top gate omitted. Its dimension is at most 64, so its
displayed arbitrary-budget count is $`O(L)`$ in the existing dirty pool.
No $`L\ge n+1`$ or additional compiler clean qubit is required.

This can be used by a state-based interference protocol. It does not
make either branch's unspecified completion equal to the prescribed
Hopf frame. Any replacement for the frame-based QBP decoder must have
its own observable, sampling, and classical reconstruction analysis.

The [reference-state checks](../tests/test_reference_state_qbp.py) verify
the ideal half-amplitude word, complex and zero residuals, the literal
one-step amplification sign, and coherent branch preparation. They are
small algebraic diagnostics, not an elementary native compiler or a
verification of the asymptotic resource proof.

## 7. A complex native reference from the same coarse word

Keep the real Hopf target psi and the actual native C from Section 2.
The reference $`C|0^n\rangle`$ need not be real. There is nevertheless
a preparation V, with the same two compiler flags, dirty allocation,
and counts, satisfying

```math
\left\|VJ_d-
 \left(|0\rangle\langle0|_d\otimes C|0^n\rangle
      +|1\rangle\langle1|_d\otimes|\psi\rangle\right)
 \otimes|00\rangle_{s,t}\otimes I_b\right\|\le\eta.
```

Here $`J_d`$ appends the zero system and two zero compiler flags to
arbitrary branch-and-dirty input. In particular, input $`|+\rangle_d`$
prepares

```math
\frac{|0\rangle_d C|0^n\rangle+|1\rangle_d|\psi\rangle}{\sqrt2}
 \otimes|00\rangle_{s,t}
```

within eta, with the literal relative plus sign and all dirty/reference
correlations covered. The protocol branch is separately counted; it is
not a third initialized compiler flag. This corollary does not apply the
two-real-state statement to the complex reference.

Use a **common** coarse word and residuals

```math
\phi_0=|0^n\rangle,\qquad \phi_1=C^\dagger|\psi\rangle,
\qquad \phi_j=a_j|0^n\rangle+v_j,\quad (v_j)_0=0.
```

Both residuals satisfy the small-tail promise, with $`a_0=1`$ and
$`v_0=0`$ exactly. On unchanged address $`(s,d,x)`$, replace the
table in Section 3 by

```math
z_{0,j,x}=a_j,\qquad z_{1,j,0}=0,\qquad
z_{1,j,x}=\sqrt N\,(v_j)_x\quad(x\ne0).
```

The accepted amplitude is exactly $`\phi_j/2`$ in each branch.
The initial reflection excludes d, so the one-step amplification proves
the coherent residual-state identity with identical literal phase.
Applying the same C afterward sends the two residuals to
$`C|0^n\rangle`$ and psi. No control of an arbitrary compiled C is
used in this construction.

For $`n\ge6`$, the eight fixed address literals and 256 sectors of
Section 6 leave $`k=n-6`$ free selectors. With $`q=L+10`$,
the dirty reservation is still

```math
(q+1)+k+1+1=L+n+7,
\qquad \|VJ_d-G_d\|\le390\,2^{-q}\lt2^{-L}\le\eta,
```

where $`G_d`$ is the displayed target isometry. The Euler certification
also covers $`a_0=1`$ and zero tail entries. Thus
$`T(V)=O(N+L)`$ and $`G(V)=O(NL)`$, with no condition
$`L\ge n+1`$. The coarse approximation creates no additional final
error because the second residual is defined using the actual C.

For $`n\le5`$, re-emit C's recorded logical two-mode words conditioned
on $`d=0`$, and a fine borrowed-native approximation to the real frame
conditioned on $`d=1`$. There are only a bounded number of logical
rows. Each determinant-one word has the
[literal reflection rewrite](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter).
Its alphabet is $`\{T^jHT^{-j}:0\le j\lt8\}\cup\{Z\}`$.
Since $`B=SHTHS^\dagger`$ obeys $`BZB^\dagger=H`$, a
predicate-controlled $`T^jHT^{-j}`$ is the same predicate-controlled
Z conjugated on its target by $`T^jB`$. Adding d to the bounded
address predicate therefore costs a constant number of exact native
gates per reflection, using one returned arbitrary helper for
controlled Z. No controlled-T gate is assumed. The coarse word length
is $`O(1)`$ here and the
fine word length is $`O(L)`$. This selects the actual logical C,
including its phase, rather than controlling its dirty implementation
gate by gate. Both compiler flags may remain unused. The resulting
coherent error is at most eta and all borrowed helpers return exactly.

The actual coarse circuit and its actual inverse implement
$`C\otimes I_b`$ and $`C^\dagger\otimes I_b`$ exactly on all inputs.
Consequently later common applications of $`C^\dagger`$ and
$`H^{\otimes n}`$ preserve the preparation error, including leaked
flag components. Their extra counts are respectively $`O(N)`$ T and
$`O(Nn)`$ Clifford gates, and n Clifford gates. Observable implementation
and sampling remain separate protocol costs.

The [bounded native integration](NATIVE_COARSE_QBP.md) emits this finite-size
fallback for an exact two-qubit target and a complex native reference.
It checks literal branch phase and an actively used arbitrary dirty helper.
It does not emit the general addressed residual-table construction.
