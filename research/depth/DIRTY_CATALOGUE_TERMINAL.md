# Logarithmic general depth by dirty catalogue commutators

[Continuing destination](../../WORKSPACE.md#continuing-toward-the-general-answer) ·
[Previous general bound R89](CAUSAL_COMPACT_COIN_DEPTH.md) ·
[Executed design](TERMINAL_SOURCE_GAUGE_PROMPT.md)

**Local result R90.** Fix integers $`n\ge1,L\ge6`$, put
$`N=2^n`$, $`\ell=\log_2(n+2)`$, and $`\eta=2^{-L}`$.
For every prescribed complete real frame W and every
$`b\ge17(L+n+7)`$, one external clean qubit and at most b arbitrary
dirty qubits suffice for a Clifford+T circuit V satisfying

```math
\|VJ_1-J_1(W\otimes I_b)\|_{\rm op}\le\eta,
```

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell).
```

All three bounds describe the same circuit. The prescribed complex
magnitude completion, including supplied leaf phases, has the same
bounds at $`b\ge17(L+n+8)`$.
The complete initialized-isometry norm includes all logical columns,
literal phases, source leakage, returned work, and arbitrary dirty and
external-reference inputs. There are no measurements, resets, QRAM,
supplied precision states, extra initialized words, or free coherent
evaluators. No efficient classical
preprocessing, all-precision uniformity, optimality, or emitted scalable
circuit is asserted by the finite catalogue construction. R79 remains
the absolute all-precision theorem. Optimal general depth remains open.

The new terminal uses two separate arbitrary arrays. A selected
catalogue factor is transferred onto the logical data together with its
inverse on one dirty slot. A second factor uses the other array, so their
commutator retains the logical operation and cancels both unwanted slot
operations exactly. An entire rounded terminal has a program of
$`O_L(2^m)`$ bits. Three R89 coarse stages make its full catalogue
small enough to instantiate, while its program table still has only
$`O_L(N)`$ bits. The long prefix is queried only a constant number
of times.

## 1. A selected small unitary and two borrowed arrays

Let P be a preserved computational prefix and let
$`v=f(P)\in\{0,1\}^q`$ specify one of $`V=2^q`$ fixed unitaries
$`W_v\in\mathrm{SU}(D)`$, where $`D=2^m`$. The logical target X
contains m qubits. Allocate two disjoint arrays $`\mathcal A`$ and
$`\mathcal B`$, each with V slots of m qubits. Every array input is
arbitrary. The arrays can be entangled with each other, X, P, query
work, and an external reference.

For each v choose fixed factors $`A_v,B_v`$ such that

```math
W_v=A_vB_vA_v^\dagger B_v^\dagger.
```

Let $`A_{\rm all}=\bigotimes_w A_w`$ act on array $`\mathcal A`$,
and let $`S_{v,\mathcal A}`$ exchange X with slot v of that array.
Throughout the operator displays, **rightmost factors act first**. Put

```math
K_A=S_{v,\mathcal A}A_{\rm all}
       S_{v,\mathcal A}A_{\rm all}^\dagger.
```

Conjugation exchanges only the selected tensor factor. All unselected
slots cancel, giving the full-space identity

```math
K_A=A_v(X)\otimes A_v^\dagger(\mathcal A_v)
       \otimes I_{\rm other\ slots}.
```

Define $`K_B`$ in the same way on the disjoint array $`\mathcal B`$.
The two unwanted factors act on different arrays and commute with both
logical factors. Therefore

```math
K_AK_BK_A^\dagger K_B^\dagger
=W_v(X)\otimes I_{\mathcal A}\otimes I_{\mathcal B}.       \qquad(1)
```

This is an identity of complete unitaries, rather than a statement about
initialized slots. Every borrowed array returns exactly, with its
reference correlations. Orthogonal direct sums over preserved P make
(1) coherent in the selected label. The distinction between the two
arrays is essential: a reused single slot would carry a nontrivial
commutator of its own.

## 2. Explicit commutator factors and literal phase

For any $`W\in\mathrm{SU}(D)`$, diagonalize

```math
W=U\,\mathrm{diag}(\lambda_0,\ldots,\lambda_{D-1})U^\dagger,
\qquad \prod_{j=0}^{D-1}\lambda_j=1.
```

Let $`S|j\rangle=|j+1\bmod D\rangle`$. Set

```math
a_0=1,\qquad a_j=\prod_{k=1}^j\lambda_k\quad(1\le j\lt D),
```

```math
A=U\,\mathrm{diag}(a_0,\ldots,a_{D-1})U^\dagger,
\qquad B=USU^\dagger.
```

The jth diagonal entry of the resulting commutator in this eigenbasis
is $`a_j/a_{j-1}`$, with indices interpreted cyclically. It equals
$`\lambda_j`$ for every j; the wraparound equality uses the determinant
condition. Thus $`ABA^\dagger B^\dagger=W`$ with literal phase.

Multiplying A or B by an independent scalar does not change its
commutator. One may therefore normalize each factor to determinant one
before its fixed synthesis. This normalization does not change (1).

The prescribed real Hopf frame is a product of embedded
$`\mathrm{SO}(2)`$ rotations, as in the
[addressed-layer interface](../../docs/HOPF_INTERFACE.md#4-complete-addressed-depth-operator).
Its determinant is one, including signed and corner pairs. Rounded real
angle frames have the same property. The prescribed complex leaf-phase
layer can remain in its existing separately paid diagonal circuit; it
need not be incorporated into these special-unitary catalogues.

The spectral argument proves existence; it does not require a new
algebraic root-isolation premise. Section 5 gives a finite native-word
length bound for approximating each factor. Enumerate pairs of words
within that bound, using a slightly stricter inner accuracy, and refine
certified entrywise errors against the prescribed W until a strict
commutator-error certificate is obtained. For example, entrywise error
less than $`\epsilon/D`$ implies operator error less than epsilon.
The existence bound, applied with an extra factor D in the accuracy,
guarantees a passing pair; that factor changes only the $`O(m)`$ term
in Section 5. Dovetailing interval refinement over this finite list
terminates without an equality or eigenvalue-degeneracy decision. No
classical preprocessing running-time bound is needed.

For $`m=0`$, the real terminal frame is the identity and this construction
is omitted. No zero-qubit catalogue routing or nonexistent target is
introduced.

## 3. Exact selection from an arbitrary dirty program

Reserve a separate q-bit register R with arbitrary input. Let $`Q_f`$
be the complete exact XOR query

```math
|P,r\rangle\longmapsto|P,r\oplus f(P)\rangle,
```

with all its own dirty helpers returned. Let $`T_R`$ translate the slot
labels of one array by XOR with the current value of R, preserving R.
The following display is **chronological**:

```math
Q_f,\quad T_R,\quad Q_f^\dagger,\quad T_R^\dagger.       \qquad(2)
```

If R starts at r, the two translations have labels $`r\oplus f(P)`$
and r. Their product is $`T_{f(P)}`$, since XOR translations commute
and are involutions. Thus (2) returns R and all query helpers, independently
of their initial values. Denote its exact completed unitary by $`T_f`$.

Let $`S_{0,\mathcal A}`$ be the fixed Clifford swap of X and slot zero.
Then

```math
S_{f(P),\mathcal A}
=T_f^\dagger S_{0,\mathcal A}T_f.                         \qquad(3)
```

Both directions of (2) are charged, and every inverse in (1)–(3) is the
actual reversed native word. The query reads only P; X, R and both
arrays are disjoint from its address. It is therefore legitimate to
reuse the same query after X has changed quantumly. No initial
computational value of X is reused as a future address.


## 4. A native translation circuit and its work

For each bit $`k=0,\ldots,q-1`$ of R, exchange every pair of slots
whose labels differ in bit k, controlled by $`R_k`$. All slot pairs in
one stage are disjoint. Each slot exchange is m controlled wire swaps.
The q stages implement $`T_R`$ exactly.

The retained
[shared-control Fredkin batch](../../docs/T_DEPTH_COMPILER.md#a-shared-control-fredkin-batch-has-at-most-four-t-layers)
uses at most four T layers, $`O(mV)`$ T and elementary Clifford gates
per stage, and no initialized or auxiliary wire. Its literal phase is
proved by the existing phase-polynomial schedule; the shared control is
used only through its charged Clifford parity computations at the
relevant boundaries. Consequently

```math
T(T_R),G(T_R)=O(mVq),\qquad D_T(T_R)=O(q).               \qquad(4)
```

The fixed swap $`S_0`$ costs $`O(m)`$ Cliffords. Array translations
never alter R or P. All query helpers are disjoint from R, the logical
input, and the two live arrays. They are reused only after the completed
exact query returns them.

## 5. Fixed catalogue circuits need no initialized ancilla

Here is an elementary polynomial local synthesis bound using only the
retained literal single-qubit approximation premise. Let
$`U\in\mathrm{SU}(D)`$ be fixed. Successive two-row Givens rotations
reduce U to a diagonal special unitary. The diagonal remainder is a
product of at most $`D-1`$ two-level diagonal rotations coupling a fixed
basis vector to each other vector. Thus $`O(D^2)`$ two-level
$`\mathrm{SU}(2)`$ rotations suffice.

For two distinct basis labels $`a,b\in\{0,1\}^m`$, an affine
Clifford basis permutation maps a to zero and b to a single-bit label.
Translate by a using X gates; choose a nonzero coordinate of
$`a\oplus b`$ and use CNOTs from it to clear all other nonzero
coordinates. This costs $`O(m)`$ Cliffords. The desired two-level
rotation is then an ordinary Euler product of three one-qubit rotations,
each controlled on every other bit being zero.

For a target Pauli $`P_t\in\{Y_t,Z_t\}`$, its fully controlled
rotation has the exact expansion

```math
\exp\!\left(-i\theta P_t
       \prod_{j\ne t}\frac{I+Z_j}{2}\right)
=\prod_{S\subseteq[m]\setminus\{t\}}
 \exp\!\left(-i\frac{\theta}{2^{m-1}}P_tZ_S\right).    \qquad(5)
```

All generators on the right commute. Every Pauli rotation is conjugate
to a one-qubit Z rotation by $`O(m)`$ Clifford gates. Equation (5)
therefore uses $`D/2`$ single-qubit rotations and no initialized
ancilla. Combining the two-level decompositions uses $`O(D^3)`$
rotations in total.

Allocate error $`O(\epsilon/D^3)`$ to each one-qubit word and use its
actual inverse wherever needed. The retained
[literal approximation premise](../../docs/FAULT_TOLERANT_COMPILER.md#primary-references)
gives the conservative simultaneous bound

```math
S(m,\epsilon)
=O\!\left(D^3[m+\log_2(1/\epsilon)+1]\right)           \qquad(6)
```

for T-count, elementary Clifford count, and T-depth of one fixed
catalogue factor, up to an absolute constant. The ordinary product error
hybrid gives error at most epsilon. This is intentionally a generous
polynomial bound; its local size will be sufficiently small. It uses no
initialized or dirty ancillary work beyond the m target wires.

All V factors of $`A_{\rm all}`$ run in parallel on disjoint slots;
the same holds for $`B_{\rm all}`$. Their counts add, while their
T-depth is bounded by (6). No shared gradient or other uncharged
non-Clifford workspace is used between slots.

## 6. Approximation error takes a maximum over the catalogue

Let $`\widehat A_v,\widehat B_v`$ be the actual native fixed
unitaries. Assume each is within $`\epsilon_A,\epsilon_B`$ of
its ideal factor, allowing an independent scalar phase if a projective
synthesis convention was used. Scalar phases cancel in a commutator.
The complete array identities are exact with hats:

```math
\widehat K_A\widehat K_B
\widehat K_A^\dagger\widehat K_B^\dagger
=[\widehat A_v,\widehat B_v]\,(X)
 \otimes I_{\mathcal A}\otimes I_{\mathcal B}.            \qquad(7)
```

Hence every dirty array returns **exactly**, even with approximate fixed
catalogue gates. All unselected tensor factors cancel as actual inverses;
the selected factors cancel separately on the two arrays. The logical
operator error is bounded by

```math
\max_v\|[\widehat A_v,\widehat B_v]-W_v\|_{\rm op}
\le2\epsilon_A+2\epsilon_B.                              \qquad(8)
```

There is no factor V in (8). Approximating the whole tensor product by a
sum over its V individual errors would discard the exact cancellations
in the physical word. The direct-sum maximum proves the same bound for
coherent P and all external references. The sole external clean qubit,
if unused by this terminal, remains untouched.

## 7. The local catalogue lemma and its paid query

Suppose the preserved address has a bits, $`Q=2^a`$ rows, and let B
be the dirty width reserved for the query, excluding R and the two live
arrays. Use the existing
[count-efficient rectangular word query](NONUNIFORM_DIRTY_INDICATOR.md#5-a-count-efficient-rectangular-blocked-query-consequence).
Its sufficient minimum pool condition is

```math
B\ge16c_I(a+2),
```

where $`c_I`$ is the absolute indicator constant from that lemma. With q
output bits, one completed query satisfies

```math
T_Q=O\!\left(q+\sqrt{Qq}+\frac{Qq}{B}\right),\qquad
G_Q=O(Qq),
```

```math
D_Q=O\!\left(q+\log_2(a+2)+\frac{Qq}{B^2}\right).      \qquad(9)
```

This is a whole-word lookup. The independent-bit bound
$`q\sqrt Q`$ is not used. The word (1) contains a fixed number of
$`K_A,K_B`$ and their inverses, each a fixed number of selected swaps,
translations and queries. Equations (4), (6), and (9) therefore give

```math
\begin{aligned}
T&=O\!\left(q+\sqrt{Qq}+\frac{Qq}{B}
          +2^q[S(m,\epsilon)+mq]\right),\\
G&=O\!\left(Qq+2^q[S(m,\epsilon)+mq]\right),\\
D_T&=O\!\left(\log_2(a+2)+q+S(m,\epsilon)
                    +\frac{Qq}{B^2}\right).
\end{aligned}                                            \qquad(10)
```

The extra simultaneous width outside the original address and logical
target is precisely the q program bits, $`2m2^q`$ array bits, and the
stated query pool B. The local fixed circuits and shared-control routers
introduce no additional initialized bank. Equation (8) supplies the error
contract, independently of arbitrary dirty inputs. All additional
integer-rounding constants and small cases can be absorbed in the fixed
local bounds without changing this complete-space identity.

## 8. A finite full-frame program with linear-in-dimension length

Fix $`L_t=L+4`$ and $`\eta_t=2^{-L_t}`$. After the retained finite
sign reduction, the remaining positive real pairs have angles in
$`[0,\pi/2]`$. At local depth $`0\le d\lt m`$, give every node

```math
k_d=L_t+m-d+6
```

program bits. A k_d-bit integer a represents
$`\theta'=(\pi/2)a/2^{k_d}`$. Certified evaluation with overlapping
admissible rounding intervals gives
$`|\theta-\theta'|\le\pi 2^{-k_d}`$ in finite time. It need not
decide equality with a grid boundary. Zero is represented exactly;
the endpoint $`\pi/2`$ can be approximated from below. These
rational multiples of pi give a finite catalogue of complete
determinant-one frames.

At one depth the addressed rotations form a direct sum, so the layer
error is the maximum node error. Telescoping the complete layers gives

```math
\|W_P-W'_P\|_{\rm op}
\le\pi\sum_{d=0}^{m-1}2^{-k_d}
\le\frac{\pi}{64}\eta_t\lt\frac{\eta_t}{16}
=\frac{\eta}{256}.                                      \qquad(11)
```

The estimate holds for every column and literal phase; the maximum
over P gives the same bound for a coherent prefix. Concatenate the
angle codes in breadth-first order. Their total length obeys

```math
\begin{aligned}
q&=\sum_{d=0}^{m-1}2^d(L_t+m-d+6)\\
 &=2^m\sum_{j=1}^{m}2^{-j}(L_t+j+6)
 \le(L_t+8)2^m.
\end{aligned}                                           \qquad(12)
```

Every one of the $`2^q`$ strings is a valid rounded angle tuple,
whether or not it occurs for a supplied prefix. Each defines an
$`\mathrm{SO}(D)`$ frame and therefore the commutator factors of
Section 2. For $`a=n-m`$, the entire program table satisfies

```math
Qq\le(L_t+8)2^{n-m}2^m=O_L(N).                           \qquad(13)
```

Assigning the same $`L_t+m`$ precision to every node would instead
introduce an unbounded factor m. The depth-dependent allocation in
(12) is needed for the linear Clifford budget.

Enumerate the finite catalogue and select native factor pairs whose
actual commutator errors are certified strictly below $`\eta/32`$.
The bounded native-net construction of Section 2 terminates: using
per-factor errors smaller than $`\eta/128`$ supplies a strict
passing margin through (8). Its word-length bound is
$`S(m,\epsilon)=2^{O_L(m)}`$ with $`\epsilon=\Theta(\eta)`$.
No division of the error by the catalogue size occurs.

If m is zero, the terminal is identity and has no program or arrays.
If Q is one, the sole-row query is simply a fixed Clifford XOR word.

## 9. Three coarse stages and the simultaneous wide allocation

Use only the local stage from
[R89 Section 10](CAUSAL_COMPACT_COIN_DEPTH.md#10-prefix-selected-recursion-and-iterated-logarithm-depth).
Its constants and protected margin are fixed there before recursion.
Start at $`m_0=n`$ and perform at most three nonfinal stages, with
the unchanged R89 choices

```math
\begin{aligned}
l_i&=L+\lceil\log_2(m_i+2)\rceil+10,\\
K_i&=m_i+l_i+2,\qquad u_i=l_i+\lceil\log_2K_i\rceil,\\
m_{i+1}&=\lceil A u_i^3\rceil .
\end{aligned}                                           \qquad(14)
```

Stop earlier if the current height enters R89's fixed cutoff
$`m_*(L,A)`$. In either case finish with the catalogue, not another
recursive stage. There are at most three coarse calls. Their already
proved simultaneous ledgers, including all source work and actual
inverses, sum to

```math
T_{\rm coarse}=O_L(\sqrt N+N/b),\qquad
G_{\rm coarse}=O_L(N),\qquad
D_{\rm coarse}=O_L(N/b^2+\ell).                         \qquad(15)
```

Indeed the protected margin suppresses each table envelope as in R89
equations (20)–(22), and a fixed number of global-address and local
precision charges is $`O_L(\ell)`$. No infinite recursive theorem
is used to obtain (15).

If all three reductions occur, their terminal heights obey

```math
\begin{aligned}
m_1&=O_L(\ell^3),\\
m_2&=O_L\!\left(\log_2^3(\ell+2)\right),\\
m_3&=O_L\!\left((\log_2\log_2(\ell+4))^3\right).
\end{aligned}                                           \qquad(16)
```

For an early stop the final height is bounded by $`m_*(L,A)`$.
Thus in every case the final height m has, as n tends to infinity
at fixed L,

```math
D=2^m=\ell^{o(1)},\qquad
q=\ell^{o(1)},\qquad S(m,\epsilon)=\ell^{o(1)}.           \qquad(17)
```

For (16), this follows from
$`(\log\log\ell)^3/\log\ell\to0`$; fixed early-stop constants
obey the same estimates. The local depth $`q+S`$ is consequently
$`O_L(\ell)`$. The full catalogue is also affordable:

```math
2^q[S(m,\epsilon)+mq]+q+2m2^q=n^{o(1)}
                   =o(\sqrt N/n).                     \qquad(18)
```

To verify the first equality, take its logarithm. It is bounded by

```math
q+O_L(m+\log(q+2)+\log(m+2))=\ell^{o(1)}=o(\ell).
```

Since $`\ell=\Theta(\log n)`$, this pays for all catalogue
slots and their actual gates, including entries never selected by P.
Counts add across those slots, while their disjoint circuits run in
parallel as proved in Section 5.

Now restrict to $`b\ge\sqrt N/n`$ and n beyond a sufficiently large
fixed L-dependent cutoff. Reserve R and **both** arrays simultaneously.
By (18) their $`q+2m2^q`$ dirty wires use at most $`b/4`$.
The native local words require no further helpers. Set
$`B=\lfloor b/2\rfloor`$ for the exact program queries, excluding
these live registers and the original logical input. Enlarging the
same fixed cutoff ensures

```math
B\ge16c_I(n-m+2),\qquad B=\Theta(b).                    \qquad(19)
```

This verifies the rectangular query's actual minimum-width condition.
It does not assert that its helper constant fits at the literal
minimum b in every dimension.

The coarse stages and the terminal run in series; their work reservations
do not overlap in time. Their dirty peaks are not simultaneous. The catalogue
uses no initialized qubit: the sole external clean flag is untouched.
Substituting (13), (17)–(19) into (10), and adding (15), proves the
three advertised resource bounds on this wide branch.

## 10. Error, signs, complex phases, and literal widths

For arbitrary signed or singular real tuples, first use R89's retained
finite whole-family magnitude/sign reduction at error $`\eta/512`$.
Its two Boolean sign phases are exact, have the retained count bounds,
and cost $`O_L(N/b^2+\ell)`$ depth. The positive frame between
those phases is factored by the coarse stages and its remaining
prefix-selected terminal.

Keep the stage precisions (14). R89's stage-error sum, also valid for
any initial segment of its recursion, is at most $`\eta/2048`$.
Equations (11) and (8), with the strict catalogue certificates of
Section 8, then give total real error less than

```math
\frac{\eta}{512}+\frac{\eta}{2048}
 +\frac{\eta}{256}+\frac{\eta}{32}
=\frac{77\eta}{2048}\lt\frac{\eta}{16}.                 \qquad(20)
```

The catalogue's physical unitary returns all its dirty work exactly
and approximates its logical target on the **entire** input space,
including occupied external-flag sectors. Its composition therefore
needs no fresh initialization after approximate coarse-source return.
Ordinary unitary telescoping with the earlier initialized-isometry
bounds includes their leakage, all columns, and external references.
The sign phases use actual completed queries and keep their literal
phase convention.

On the wide branch, append the prescribed complex leaf diagonal using
the same one-clean, error-$`\eta/4`$ paid interface already used in
R89. It returns its work and has the same count orders and
$`O_L(N/b^2+\ell)`$ depth. Together with (20), its total error is
less than eta at the original complex width.

For $`b\lt\sqrt N/n`$, use the retained **complete**
[R72 compiler](ONE_CLEAN_DEPTH.md#1-the-one-clean-theorem) for the
respective real or complex family, at the original tolerance eta.
That complete fallback already includes the prescribed complex phases;
no additional approximate diagonal is appended to it. Its depth is
$`O_L(N/b^2+n)`$, which is $`O_L(N/b^2)`$ because
$`N/b^2\gt n^2`$. Use the same respective complete compiler in
every dimension below the fixed cutoff, for any supplied b. In those
finitely many dimensions its bounded additive depth is absorbed into
$`O_L(\ell)`$.

These fallbacks retain error at most eta and the original reservations
$`b\ge17(L+n+7)`$ for real frames and $`b\ge17(L+n+8)`$ for
complex magnitude frames. A stricter internal tolerance from the wide
construction is never substituted into a fallback's width formula.
All branches therefore prove the theorem in the opening statement.

## 11. Scope and bounded evidence

The [dirty-catalogue tests](../../tests/test_dirty_catalogue_terminal.py)
check the complete small word, with native H/S/T catalogue factors and
exact routing permutations, on every input column,
coherent prefixes and explicit dirty/reference entanglement, nonlinear
program selection by two-bit XOR translations, and the failure of a
single reused array. They also check commutator factors for a complete
four-dimensional real frame and the literal shared-control phase
identity. These are bounded operator checks. The asymptotic catalogue
size, native approximation bound, precision allocation and width
splice are proved above; no scalable circuit emitter is claimed.

This is a fixed-L theorem. Its catalogue has $`2^{O_L(2^m)}`$
entries, and the dimension cutoffs may depend on L. It does not replace
R79's absolute all-precision theorem or give uniform constants when L
grows with n. Terminating classical preprocessing can be enormous, and
no practical crossover is claimed. The catalogue is a charged circuit
on arbitrary dirty arrays, not a supplied quantum memory or evaluator.
The selected family is essential to (13): a general selected
$`\mathrm{SU}(2^m)`$ family does not automatically have an
$`O_L(2^m)`$-bit program.

R90 removes R89's unbounded iterated-logarithm multiplier from the
general fixed-accuracy upper bound. It does not prove that logarithmic
depth is necessary. Optimal large-width depth and a precision-uniform
extension remain separate open questions.
