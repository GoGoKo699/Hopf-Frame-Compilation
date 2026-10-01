# Linear precision cost for prefix-free tree updates

[Research status](OPEN_PROBLEM.md) · [Borrowed compiler](BORROWED_WORKSPACE_COMPILER.md) · [Borrowed signal](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)

Arbitrary local updates at prefix-free tree nodes can share one precision
charge. Their descendant circuits are inherited literally from a native
baseline, so those circuits can be factored out exactly. A charged permutation
then puts all relative updates into one ordinary one-target multiplexor.
The resulting bound is

```math
a=0,\qquad b\ge L+n+7,\qquad
T=O(N+L),\qquad G=O(NL).
```

This is a promised family of complete tree operators. It includes a real
Hopf specialization with arbitrary real-angle changes at an antichain and
a nontrivial branching real baseline. It does not establish this bound for
an arbitrary prescribed real frame. In particular, independently rounding
all its local angles does not generally produce the required prefix-free
pattern of changes.

## 1. Target and native-baseline promise

Set $`N=2^n`$, $`n\ge1`$, $`0\lt\eta\le1/64`$, and
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
A depth-d node has binary prefix $`v\in\{0,1\}^d`$ and local
determinant-one two-by-two word $`U_v`$. Its logical gate acts on
the ordered pair

```math
|v\,0^{n-d}\rangle,\qquad
|v\,1\,0^{n-d-1}\rangle
```

and is identity on all other basis states. Layers are composed in the
same shallow-to-deep order as the complete Hopf frame. The local words may
be complex SU(2) matrices; no arbitrary local scalar phase is included.

The supplied baseline C consists of actual literal Clifford+T local words
$`U_v^C`$. Their native word lengths, including Clifford gates, satisfy

```math
|U_v^C|\le c_0(n-d+1)
```

for a fixed constant $`c_0`$ independent of n and L. The local words
produced by the retained compiler at fixed accuracy are one instance:
$`U_v^C=C_v^2`$ with $`C_v=XQ_v^\dagger XQ_v`$ and
$`|Q_v|=O(n-d+1)`$. The statement also allows other literal
determinant-one native words with this length envelope. The proof below
implements their addressed logical gates and any required row masks; these
are not supplied as free controlled subroutines.

Let $`\mathcal A`$ be an antichain of internal nodes: no member's
binary prefix is a strict prefix of another's. The target W has arbitrary
effectively specified SU(2) local words $`U_v^W`$ on
$`v\in\mathcal A`$ and satisfies the literal equality

```math
U_v^W=U_v^C\qquad(v\notin\mathcal A).
```

There is no closeness assumption on the changed words. The antichain may
contain nodes at different depths and as many as $`N/2`$ nodes; the
restriction is at most one changed node on each root-to-leaf path.
The empty antichain and the antichain containing only the root are allowed.

**Theorem.** Under these promises there is a coherent Clifford+T circuit
$`\widetilde W`$, using no initialized ancilla and at most
$`b=L+n+7`$ arbitrary dirty qubits, such that

```math
\|\widetilde W-W\otimes I_b\|\le\eta,
\qquad T=O(N+L),\qquad G=O(NL).
```

The norm covers every logical and dirty input, including correlations with
arbitrary reference systems. Any larger dirty allocation may leave the
extra wires untouched. The inverse is the actual circuit inverse and has
the same error guarantee. There are no measurements, resets, supplied
catalysts, initialized precision registers, or uncharged quantum oracles.
Certified coefficient evaluation and classical table preparation are
separate preprocessing costs.

## 2. Factor out the inherited descendant forest

Let $`V_{\mathcal A}`$ be the baseline circuit containing precisely
the strict descendants of marked nodes, in their original relative order.
Let $`M_C`$ and $`M_W`$ be the products of the marked local gates,
and let B contain every remaining baseline node. Gates on incomparable
subtrees have disjoint logical supports and commute. Keeping the order of
comparable nodes and moving only incomparable gates therefore gives

```math
C=V_{\mathcal A}M_CB,\qquad
W=V_{\mathcal A}M_WB.
```

The descendant forest is the same literal circuit in both factorizations.
Define $`D_v=U_v^W(U_v^C)^\dagger`$ and let
$`M_{\mathcal A}`$ apply $`D_v`$ on the marked node's pair.
The marked pairs are disjoint, so

```math
M_{\mathcal A}=M_WM_C^\dagger,
\qquad
W=V_{\mathcal A}M_{\mathcal A}V_{\mathcal A}^\dagger C.
```

This is an exact complete-operator identity, with no overlap approximation
or truncated residual expansion. It remains true for noncommuting complex
local words. It does not call W as a subroutine: C and the descendant
forest are compiled from the recorded native baseline words.

### Exact compilation of the baseline and its masks

The retained
[reflection interpreter](BORROWED_WORKSPACE_COMPILER.md#2-exact-dirty-table-and-reflection-interpreter)
rewrites a determinant-one native word into a constant-factor longer word
over the fixed involution alphabet

```math
\{T^jHT^{-j}:0\le j\lt8\}\ \cup\ \{Z\}.
```

The rewrite preserves literal phases. For each word position and alphabet
symbol, a Boolean table f on the depth-d prefix specifies which rows apply
that symbol. The corresponding addressed reflection F satisfies
$`F^2=I`$ on its complete input space. An unbanked dirty query and the
retained interpreter implement F controlled by an arbitrary dirty bit beta
at cost $`O(2^d)`$, using at most d dirty selectors and one dirty bank.

Let h be the predicate that the logical suffix below this node is zero.
Toggle an additional dirty beta by h using a reversible predicate circuit
$`Q_h`$. Apply the chronological sequence

```math
Q_h,\quad C_\beta F,\quad Q_h^\dagger,\quad C_\beta F.
```

The two reflection exponents are $`\beta\oplus h`$ and beta,
so the resulting action is exactly $`F^h`$. Prefix and suffix bits
are unchanged by F, and the table interpreter returns all its work before
the predicate is uncomputed. A borrowed multi-controlled-X construction
implements $`Q_h`$ with $`O(n^2)`$ gates and one arbitrary helper.
When the suffix is empty, h is identically one: drop the beta control in
the interpreter's dirty-bank reflection echo to implement F directly.
No arbitrary dirty bit is assumed to equal one. Every helper returns
exactly on its full input space.

There are at most d selectors, one bank, beta, and one predicate helper,
so $`n+2`$ dirty wires suffice throughout. Complete this echo for
each reflection before advancing to the next word position; different
reflections need not commute. Replacing any unwanted row word by identity
implements an arbitrary classical node mask at the same cost. In particular,
it produces the strict-descendant forest $`V_{\mathcal A}`$.

For either C or the forest, the charged gate counts are

```math
\begin{aligned}
T,G
&=O\!\left(\sum_{d=0}^{n-1}(n-d+1)(2^d+n^2)\right)\\
&=O(N+n^4)=O(N).
\end{aligned}
```

Thus their actual inverses also cost $`O(N)`$. The equality
$`n^4=O(2^n)`$ uses an absolute constant and does not omit the
predicate work. No approximation error is introduced by these baseline
or mask circuits.

## 3. Pack all marked pairs into one addressed target

For a marked prefix v at depth d, swap logical bit $`d+1`$ with the
last logical bit, conditioned on the prefix v. This maps its two active
states to

```math
|v\,0^{n-1-d}\rangle\otimes|0\rangle,
\qquad
|v\,0^{n-1-d}\rangle\otimes|1\rangle.
```

Each swap preserves the entire prefix-v cylinder. Since the marked
cylinders are disjoint, the swaps commute and do not alter another marked
pair's address. The padded addresses $`v0^{n-1-d}`$ are distinct:
equality of two would force one marked prefix to extend the other.
Let P be the product of these swaps. Then

```math
D=P M_{\mathcal A}P^\dagger
=\sum_{x\in\{0,1\}^{n-1}}|x\rangle\!\langle x|\otimes D_x,
```

where $`D_x=D_v`$ at a marked padded address and $`D_x=I`$
otherwise. This is a single SU(2) multiplexor on the last logical qubit.
No initialized address or history register has been introduced.

### Charge the packing permutation

For each depth d let $`f_d`$ be the Boolean membership table for
marked prefixes of that length. The exact dirty query $`Q_{f_d}`$
toggles an arbitrary dirty bit z by $`f_d`$, restoring its d dirty
selectors. Let $`F_z`$ swap logical bits $`d+1`$ and n under z.
The chronological word

```math
Q_{f_d},\quad F_z,\quad Q_{f_d}^\dagger,\quad F_z
```

implements that swap controlled by $`f_d`$, since Fredkin is an
involution. Crucially, the swapped bits are outside the first d address
bits. The query's address therefore does not change before it unloads.
Both z and the selectors return exactly on arbitrary inputs. At
$`d=n-1`$ the requested swap is identity and can be omitted.

An unbanked query costs $`O(2^d)`$ T and Clifford gates, and each
Fredkin costs a constant number of native gates. Hence P and its actual
inverse each have $`O(N)`$ T and Clifford cost and use at most
$`n+1`$ dirty wires. The construction specifies the permutation on
all logical states, not only on its marked pairs.

## 4. One precision charge for the packed SU(2) table

Choose three Euler rotations for every addressed $`D_x`$,
in a common $`R_zR_yR_z`$ order, with total matrix error at most
$`2^{-L}/4`$. Certified finite approximation suffices, including
singular Euler charts; exact zero tests on arbitrary computable angles
are unnecessary. For example, certified search over successively finer
angle grids supplies a triple within the stated matrix tolerance.
Classical running time is not included in the quantum bound. The complete
addressed Euler table has the same error, because its address sectors are
orthogonal and invariant.

Apply the
[borrowed-signal primitive](ONE_CLEAN_COMPILER.md#9-a-borrowed-signal-suffices-for-real-rotations)
to these three addressed rotations. Fixed target Clifford conjugations
give the Rz versions. With precision q, k free address bits, and p unchanged
predicate bits, the primitive has full-operator error less than
$`43\,2^{-q}`$ and costs

```math
T=O(2^k+q+p^2),\qquad
G=O(2^kq+q+p^2),\qquad
b=q+k+3.
```

Its source core and borrowed amplification signal return approximately
within that norm. Its selectors and predicate helper return exactly, and
an inactive predicate gives exact identity on the entire input space.
The literal-phase primitive is not needed: every $`D_x`$ has
determinant one, so no extra addressed scalar phase is being suppressed.

Set $`q=L+10`$. If $`n-1\ge5`$, take five prefix bits as fixed
sector literals and the remaining $`k=n-6`$ bits as the free table
address. The dirty ledger is exactly

```math
(q+1)+k+1+1=L+n+7.
```

The entries are respectively the precision core, free-address selectors,
predicate helper, and one arbitrary borrowed amplification signal.
The 32 sectors act on disjoint unchanged addresses, and each circuit is
exact identity on all other sectors. Their errors take a maximum even
when their common dirty pool is arbitrary or already correlated with data.
Summing only the three Euler-factor errors gives

```math
\|\widetilde D-D\otimes I_b\|
\lt\left(\frac14+\frac{129}{1024}\right)2^{-L}
\lt\eta.
```

There are a fixed number of sectors and factors. Consequently

```math
T(\widetilde D)=O(N+L),\qquad
G(\widetilde D)=O(NL).
```

For $`n-1\lt5`$, there are fewer than 32 rows. Synthesize each
Euler rotation by the direct native borrowed-sector echo, with error at
most $`2^{-L}/16`$ per row. Disjoint row errors take a maximum;
the three Euler factors and their preprocessing therefore contribute at
most $`(3/16+1/4)2^{-L}\le\eta`$. This uses $`O(L+n)`$ native word
length per row and at most $`O(n^2)`$ predicate gates, with exactly
returned borrowed helpers. The number of rows and n are bounded constants,
so the same $`O(N+L)`$ and $`O(NL)`$ bounds hold without reserving
the large q-bit core. When $`n=1`$ there is no address predicate:
the literal phase-cancelling word $`(XQ^\dagger XQ)^2`$ on the
target itself approximates each Euler rotation, using quarter-angle Q
synthesis. These direct cases fit the stated dirty threshold.

## 5. Complete-input error and the real Hopf specialization

Assemble the actual circuit

```math
\widetilde W
=V_{\mathcal A}P^\dagger\widetilde D P
 V_{\mathcal A}^\dagger C,
```

with the exact circuits extended by identity on the shared dirty pool.
All helper roles are reused sequentially. Exact baseline, forest, and
packing implementations remain valid on arbitrary dirty inputs, including
the small disturbance left by $`\widetilde D`$. Unitary invariance
of the norm and the factorization in Section 2 therefore give

```math
\|\widetilde W-W\otimes I_b\|
=\|\widetilde D-D\otimes I_b\|\le\eta.
```

This proves the theorem with one precision-dependent multiplexor, no
initialized flag, and charged actual inverses. The borrowed core and
amplification signal have approximate return included in the full error;
the query, packing, and exact baseline helpers return exactly. There is
no reset between these roles.

For a concrete nontrivial real family, set every baseline local word to
$`R_y(\pi/4)=HZ`$. These are literal length-two determinant-one
native words, so Section 2 compiles the complete branching baseline and
any descendant mask in $`O(N)`$ gates. Change arbitrary real angles
at any prefix-free set of nodes and retain $`\pi/4`$ literally
elsewhere. The resulting W is the prescribed complete real Hopf frame,
and the full-operator approximation above implies its retained QBP error
guarantees. The identity baseline gives another specialization.

A baseline obtained by generic native approximation of real rotations
may instead have complex local words. Requiring literal agreement with
those words outside the antichain then defines a complex SU(2) tree
operator, not automatically a real Hopf frame. Nor does this theorem
assert a new compiler for the separate phase-dressed complex magnitude
family. These distinctions are part of the target promise.

At the selected precision $`L=N`$, the promised family has
$`T=O(N)`$, $`G=O(N^2)`$, with $`b=N+n+7`$ and no clean
ancilla. Arbitrary changes on comparable nodes are outside this theorem:
they alter the inherited descendant forest and cannot be combined by
the antichain factorization above. Thus the construction gives a linear-T
endpoint for the stated update family, while the unrestricted endpoint
remains open. No matching lower bound or T-depth improvement is claimed.
