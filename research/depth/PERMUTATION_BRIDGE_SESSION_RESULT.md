# The in-place permutation checkpoint

[Current frontier](../../docs/OPEN_PROBLEM.md) · [Orientation reduction](BRANCH_ORIENTATION_SESSION_RESULT.md) · [Shallow queries](SHALLOW_DIRTY_QUERIES.md)

The session beginning on 15 October 2026 at 20:43:55 UTC started from
commit `3f963a5`. Its authorized ceiling was five hours, with an initial
ninety-minute checkpoint for a new native mechanism. The construction
search stopped early after its independent lanes exhausted the identified
candidates without a qualifying circuit. Verification and integration
then preserved the outcome.

**Outcome at this checkpoint: the arbitrary-orientation problem remained open.** No new
general compiler, improved depth bound, or unrestricted lower bound is
claimed. R79 remains the selected general theorem; R80–R84 retain their
existing scopes. The restricted obstructions below explain why particular
conversions fail. They are not substitutes for the session's target.

The later [joint orientation/recovery session](JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
reaches the sub-square-logarithmic target under shared unordered depth
magnitudes through R85's complete recovery word. This historical conversion
checkpoint does not describe the current frontier.

**Continuing direction, revised on 16 October 2026 (Asia/Shanghai):**
the session ended, while the general constructive objective remains
active. The user's revision supersedes the recommendation to make
independent review the default next phase. Discovering a new mechanism
is itself part of the work; its absence at this checkpoint is not a
reason to retire the question. See the [continuing mandate](../../WORKSPACE.md#continuing-toward-the-general-answer).

## 1. The sufficient target and its scope

Let $`N=2^n`$. A binary tree automorphism has a triangular leaf map

```math
F(x)_j=x_j\mathbin\oplus f_j(x_1,\ldots,x_{j-1}).
```

The task was to compile its in-place permutation $`P_F`$ with one
initialized qubit, $`b=\Theta(\sqrt N)`$ arbitrary dirty qubits,
$`T=O(\sqrt N)`$, $`G=O(N)`$, and
$`D_T=o(\log^2(n+2))`$, ideally logarithmic. The approximation, if any,
had to hold in the full initialized-isometry norm, including work return
and arbitrary external references. A fixed sufficient width constant
was allowed. No extra clean output register was supplied.

The previous session already supplies exact forward and inverse XOR
queries within the count and depth budgets. It also proves

```math
W_r=P_F W_0 Q,\qquad B=P_F BQ,
```

where B is the balanced frame. Consequently

```math
W_r=P_F W_0 B^\dagger P_F^{-1}B.
```

Thus a full-register in-place compiler is sufficient: a constant number
of returned-work factors gives the entire oriented frame, including the
signed marker columns. A separate selected-prefix implementation of Q
is unnecessary for this reduction. For an approximate P, use its actual
inverse and allocate a constant fraction of the final error to each
factor; the operator-norm telescoping bound charges all factors.

This would settle arbitrary orientation tables when the unordered
magnitude pair is shared at each depth. It would not automatically settle
arbitrary node-dependent magnitude pairs or optimal unrestricted depth.
Failure to construct P also leaves open a direct compiler for the complete
word above.

## 2. What the external techniques actually supply

These are source checks for proposed mechanisms, not new imported
premises of the selected theorems.

| Primary source | Relevant interface | Why it does not complete this conversion |
|---|---|---|
| [Khattar–Gidney, conditionally clean ancillas, sections 3–5](https://arxiv.org/html/2407.17966v2) | A preserved condition certifies the borrowed register's value; laddered toggle detection handles the specified nested self-inverse controlled operations | An arbitrary dirty output word has no such certificate; replacing it requires a new return identity |
| [Bona, sections 1–2](https://arxiv.org/html/2608.08765v1) | Allocates and schedules borrowed blocks already satisfying their dirty-return contracts | Scheduling does not construct a missing identity factor on the borrowed register |
| [Holman–Ramachandran–Yirka, section 2](https://arxiv.org/html/2504.03620v1) and [Dong–Lombardi–Ma, section 1.1](https://arxiv.org/html/2607.26478v1) | Forward and inverse XOR queries give the usual two-query conversion with an initialized output word | Inverse access is already available here; the clean output is not. Forward-only or one-query lower bounds do not apply to our native model |
| [Salo, coordinate alternation](https://arxiv.org/abs/1901.03613) | Three alternating coordinate updates realize permutations of a Cartesian product | Each update may be an arbitrary controlled permutation, rather than a priced XOR translation; its use here would leave an in-place subproblem inside a factor |

The two-query identity is

```math
|x\rangle|0^n\rangle
\longmapsto |x\rangle|F(x)\rangle
\longmapsto |0^n\rangle|F(x)\rangle
\longmapsto |F(x)\rangle|0^n\rangle.
```

It is correct, but starts with n clean output bits. Neither the existence
of inverse queries nor a dirty-qubit scheduling pass changes that premise.
The existing [dirty-output counterexample](BRANCH_ORIENTATION_SESSION_RESULT.md#2-a-cheap-xor-oracle-does-not-complete-the-permutation)
already rules out applying this word unchanged to arbitrary outputs.

## 3. Old-address recovery by unary wire relabeling

There is a valid new scheduling identity. Let $`I_a`$ toggle an arbitrary
unary bank Y by $`e_a`$. For a known permutation G of its address set,
let $`R_G`$ permute the bank wires so that $`R_G e_a=e_{G(a)}`$.
After the logical address has changed to $`a'=G(a)`$, the old indicator
can still be applied by the operator word

```math
R_G^{-1}I_{a'}R_G:\quad Y\longmapsto Y+e_a.
```

This equality holds for every bank input, not only a zero bank. When G
is unselected, $`R_G`$ is a Clifford wire permutation. In a four-corner
Boolean lookup, a compiler for G can therefore overlap the middle private
segment between two address-indicator appearances, if that segment has
no live interaction with the changing address. Actual inverses restore
the bank.

This identity did not yield a recursive schedule with the required
resources. At the first equal split, an address of n/2 bits has a bank
of square-root width. A recursive subproblem, however, has an outside
prefix p and needs the selected wire permutation $`R_{G_p}`$. This is
a quantum multiplexer, not a free wire relabeling. Its canonical tree
router has sequential selected-swap stages; no shallower replacement
meeting the ledger was obtained.

Including p in the bank makes the relabeling unconditional, but uses
$`2^{|p|+|a|}`$ wires. After one equal split, the right subproblem's
next half-address already has $`|p|+|a|=3n/4`$. The resulting bank
exceeds square-root width by an unbounded factor. Reinterpreting a
classical table does not restore the physical dirty bank.

Nor does a single masked indicator provide arbitrary subsequent
controlled operations: extracting a selected phase or XOR by an echo
uses its commuting structure. A whole family of noncommuting conditional
permutations needs an additional identity. The new overlap is therefore
a correct local identity with an unresolved recursive cost, not a compiler.

## 4. Four XOR shears with one clean target bit

This is a restriction on a specified classical conversion word. It is
not a lower bound for general native circuits, for all query algorithms,
or for arbitrary numbers of XOR shears.

Write all sums below over $`\mathbb F_2`$. Let
$`V=\mathbb F_2^n`$ contain the logical input x and let
$`W=\mathbb F_2^{d+1}`$ contain a word y. The allowed initial y values
form the linear subspace

```math
S=\mathbb F_2^d\times\{0\}\subseteq W.
```

Thus the last coordinate of this fixed Y register is the one clean bit;
the other d coordinates are arbitrary dirty inputs. Define total Boolean
word functions $`f,h:V\to W`$ and $`g,k:W\to V`$, and the reversible
XOR shears

```math
A_f(x,y)=(x,y+f(x)),\qquad
B_g(x,y)=(x+g(y),y).
```

Consider exactly the four shears
$`A_f,B_g,A_h,B_k`$ in chronological order. There are no interleaved
basis changes, affine mixers, register repartitions, measurements, or
additional logical operations. Each shear may be realized using a
returned-work native query; the proposition concerns its complete
Boolean action at the four interfaces. The functions may be chosen
arbitrarily for the target, and no bound on d or on their lookup cost is
needed.

**Proposition.** If a bijection $`F:V\to V`$ satisfies

```math
B_k A_h B_g A_f(x,y)=(F(x),y)
\qquad(x\in V,\ y\in S),
```

then $`F^8=I`$. If Y has no initialized coordinate and the same return
identity instead holds for every $`y\in W`$, then $`F^4=I`$.

### Proof with one clean target bit

The Y value after the first shear is $`t=y+f(x)`$. Since the last shear
leaves Y unchanged, the claimed final action gives both

```math
F(x)=x+g(y+f(x))+k(y),
\qquad
h(F(x)+k(y))=f(x).
```

Put

```math
k_0=k(0),\qquad
K=\mathrm{span}_{\mathbb F_2}
\{k(y)+k_0:y\in S\}\subseteq V.
```

Comparing the second identity at y and at zero, and using surjectivity
of F, proves

```math
h(z+v)=h(z)\qquad(z\in V,\ v\in K).
```

Define $`H(z)=h(z+k_0)`$. Then H is K-invariant and
$`f(x)=H(F(x))`$. The first identity at zero therefore gives

```math
F^{-1}(z)=z+u(z),\qquad
u(z)=g(H(z))+k_0.
```

In particular u is K-invariant. Hence $`F^{-1}`$ is K-equivariant and
induces a bijection on $`V/K`$.

Comparing the first identity at y and zero also gives

```math
g(s+y)+g(s)=k(y)+k_0\in K
\qquad(s\in\mathrm{Im}H,\ y\in S).
```

If two values $`s,s'\in\mathrm{Im}H`$ lie in the same coset
modulo S, choose $`y=s+s'\in S`$. Their g values are then equal modulo
K. Since $`W/S`$ has two elements, the values of u occupy at most two
cosets modulo K.

Consequently the induced permutation on $`V/K`$ has at most two
translation directions, say r and s. Each orbit stays in one coset of
$`\mathrm{span}\{r,s\}`$, which has at most four elements.
If r and s are independent, the allowed edges form the four-vertex
square, whose permutation cycles have lengths two or four. If they are
dependent or one is zero, the allowed components have at most two
vertices, with possible fixed points. Thus the fourth power of the
quotient permutation is the identity.

It follows that $`F^4(z)=z+v(z)`$ for a K-valued function v.
The map $`F^4`$ is K-equivariant, so v is K-invariant. Therefore

```math
F^8(z)=z+v(z)+v(z+v(z))=z.
```

This proves the one-clean assertion. For fully arbitrary Y, take
$`S=W`$ throughout. There is then only one quotient translation
direction. Its square is the identity, and the same lift argument gives
$`F^4=I`$.

### Sharpness and application to triangular maps

The exponent eight is attained with three logical bits and only one
dirty target bit. Let $`X=(a,b,c)`$, $`Y=(d,0)`$, and choose

```math
f(a,b,c)=((1+a)(1+b),1+a),\qquad
g(y_1,y_2)=(1,y_2,y_1),
```

```math
h(a,b,c)=(ab,a),\qquad
k(y_1,y_2)=(0,0,y_1).
```

After the first two shears, X is

```math
(a+1,\ b+a+1,\ c+d+(1+a)(1+b)).
```

The h value on that word equals the original f value: its first entry is
$`(a+1)(b+a+1)=(a+1)(b+1)`$, and its second entry is $`a+1`$.
The third shear thus returns Y to $`(d,0)`$. The fourth removes d from
the last coordinate of X. The final logical map is

```math
F(a,b,c)=(a+1,\ b+a+1,\ c+(1+a)(1+b)),
```

which is decrement modulo eight when a is the least significant bit.
It is an eight-cycle. This also supplies a small complete-input
sharpness fixture; it uses no numerical approximation.

For comparison, increment on four bits in least-significant-bit-first
order is the triangular permutation

```math
F(x)_1=x_1+1,\qquad
F(x)_j=x_j+\prod_{i\lt j}x_i\quad(2\le j\le4).
```

It has order sixteen, so no four-shear word in the stated one-clean
model implements it. Already the three-bit increment excludes the
fully dirty version. These are triangular permutations and therefore
valid leaf maps of binary tree automorphisms.

Every word in this alphabet is a computational-basis permutation. If its
initialized-isometry action differs from the desired map on any allowed
basis input, the two output vectors have Euclidean distance
$`\sqrt2`$. Thus operator-norm approximation with error strictly below
$`\sqrt2`$ would require exact agreement and is excluded by the same
examples. The statement includes arbitrary superpositions and external
references because equality on every allowed basis vector extends
linearly.

The fixed placement of the clean bit in Y and the fixed four-shear
alphabet are part of the proposition. It does not analyze an algorithm
that moves or encodes the clean bit between different register
partitions, introduces intervening Clifford mixers, or uses five or more
shears. In particular, the easy native increment circuit is no
counterexample: it lies outside this conversion model.

## 5. Exact encoded Clifford actions need larger width

This proposition tests an exact encoding shortcut. It does not bound
the depth of general Clifford+T circuits, and its width conclusion does
not survive arbitrary approximate encodings.

Let $`N=2^n`$, and let $`P_N`$ be the permutation matrix of one
N-cycle, for example binary increment with the least significant bit
first. The latter is a triangular tree automorphism. Let D be any
nonzero dirty/reference input space.

**Proposition.** Suppose an isometry
$`E:\mathbb C^N\otimes D\to(\mathbb C^2)^{\otimes q}`$, a
q-qubit Clifford C, and a phase phi satisfy the exact intertwining
identity

```math
CE=e^{i\phi}E(P_N\otimes I_D).
```

Then

```math
N\le2^{1+\lceil\log_2(2q)\rceil}\le8q,
\qquad q\ge N/8.
```

The encoder is arbitrary: it need not be Clifford, efficiently
implementable, or specified by a particular unary representation. The
physical width q includes every wire used by this encoded Clifford
action. Extra dirty inputs do not weaken the conclusion.

### Projective Clifford order bound

Conjugation by a q-qubit Clifford acts linearly on the binary labels of
the Pauli operators and preserves their commutation form. This gives
a symplectic matrix $`S_C\in\mathrm{Sp}(2q,2)`$.

For completeness, the kernel of this action modulo scalar phases
consists of projective Paulis and has exponent two. Indeed, a Clifford
fixing every Pauli label can only attach signs to the Hermitian Paulis.
Multiplication of Paulis shows that these signs define a binary linear
functional of the label. Nondegeneracy of the symplectic form realizes
that functional as conjugation by one Pauli. Removing this Pauli leaves
an operator commuting with every Pauli, which is scalar because the
Paulis span the full matrix algebra.

Write the order of $`S_C`$ as $`m2^a`$ with m odd. Then
$`S_C^m`$ has power-of-two order and is unipotent over
$`\mathbb F_2`$: its minimal polynomial divides
$`X^{2^a}-1=(X-1)^{2^a}`$. On its $`2q`$-dimensional vector space,
the nilpotent part therefore has nilpotency index at most $`2q`$.
Consequently

```math
a\le\lceil\log_2(2q)\rceil.
```

If r is the projective order of C, the kernel statement implies that r
divides twice the order of $`S_C`$. Its power-of-two part is therefore
at most

```math
2^{1+\lceil\log_2(2q)\rceil}\le8q.
```

This argument also proves that the projective order is finite; a
chosen overall scalar phase for C has no effect.

### Passing the order through an exact encoding

Since $`C^r`$ is scalar, the intertwining identity implies that
$`P_N^r\otimes I_D`$ is scalar. The spectrum of $`P_N`$ contains
both one and a primitive Nth root of unity. Thus $`N\mid r`$.
Because N is a power of two, it is at most the power-of-two part of r,
which proves the asserted width bound. The freely allowed phase phi
has canceled in this projective argument.

It follows that an exact invariant encoding followed by one Clifford
cannot realize every triangular permutation at physical width
$`q=\Theta(\sqrt N)`$. This excludes the proposed exact free-Clifford
bank shortcut, even before the cost of constructing and undoing E is
charged. It does not exclude multiple interleaved non-Clifford stages.

### Why the exact bound has no immediate fixed-error version

There is an explicit spectral counterexample to extending this argument
to arbitrary approximate encodings. Put $`q=n+1`$ and
$`M=2^q-1=2N-1`$. Choose a basis of the finite field
$`\mathbb F_{2^q}`$ over $`\mathbb F_2`$ and a generator alpha of its
nonzero multiplicative group. Multiplication by alpha is an invertible
binary linear transformation, hence a CNOT-implementable Clifford on q
qubits. It fixes the zero word and permutes the nonzero words in one
M-cycle. Its spectrum on that cycle contains every Mth root of unity.

Take a Fourier eigenbasis of $`P_N`$. Map its eigenvector with phase
$`2\pi j/N`$ to an M-cycle eigenvector with a nearest available phase
$`2\pi k_j/M`$. These N nearest choices are distinct because
$`M/N\gt1`$; the cyclic gaps between consecutive ideal indices exceed
one. The resulting map E is an isometry and obeys

```math
\|CE-EP_N\|_{\rm op}\le\pi/M=O(1/N).
```

The estimate follows by diagonal comparison: the circular phase error
is at most $`\pi/M`$, and chord length is at most angular distance.
Tensoring E and C with identity on arbitrary dirty inputs preserves
the same bound. This uses only $`n+1`$ qubits for the encoded logical
component.

No efficient implementation of this E is supplied. The example shows
precisely why a fixed-error exclusion must price the encoder or use an
additional invariant: an exact order divisibility argument alone cannot
do so. It supplies no circuit within the research target's count,
workspace-return, and depth budgets.

## 6. Corner words with coordinate permutations and bit flips

This is a restriction on words of completed logical corner frames from
R84 and their inverses. Between them allow wire permutations, X gates,
and arbitrary computational diagonals, all on the n logical qubits.
For $`n\ge2`$, some triangular permutation requires at least

```math
k\ge\frac{(n-1)N}{4(N-1)}=\Omega(n)
```

corner factors. The witness is itself the Clifford CNOT fanout

```math
F(x)_1=x_1,\qquad F(x)_j=x_j+x_1\quad(2\le j\le n).
```

Thus the restriction gives no obstruction with general CNOT or Hadamard
interleavers. Dirty helpers returned inside each completed factor do not
change its logical word; interleavers coupling logical and helper systems
are outside the argument.

### Coordinate crossing rank

Let $`E_{j,b}`$ project onto logical bit j equal to b, and define

```math
r_j(U)=\mathrm{rank}(E_{j,1}UE_{j,0}),\qquad
R(U)=\sum_{j=1}^n r_j(U).
```

For a fixed coordinate, inserting $`I=E_0+E_1`$ gives

```math
E_1UVE_0=E_1UE_0VE_0+E_1UE_1VE_0.
```

Consequently $`r_j(UV)\le r_j(U)+r_j(V)`$ and R is subadditive.
For a permutation with unit-modulus entry phases, $`r_j`$ counts basis
inputs crossing from bit zero to bit one. The opposite crossing count
is equal by bijectivity. Conjugation by wire permutations or X words
therefore preserves R; diagonal conjugation preserves each rank.

A deterministic frame is an ordered product of signed local corner
rotations on the computational-label tree. At a node p of depth d, the
edge joins

```math
p0^{n-d}\quad\hbox{and}\quad p1\,0^{n-d-1}.
```

A diagonal corner has R zero. An active corner swaps the endpoints
with signs and has R one, only at coordinate $`d+1`$. There are
$`N-1`$ edges, so subadditivity proves $`R(C)\le N-1`$ for each
completed corner frame and its inverse. This bound survives conjugation
by any of the allowed interleavers.

### Eliminate the final interleaver

The allowed interleavers form a group. Moving their cumulative products
to the right rewrites any candidate word as

```math
P=D_1\cdots D_k A,
```

where each $`D_i`$ is a corner or inverse conjugated by an allowed
interleaver and A is the final accumulated interleaver. Hence
$`R(PA^{-1})\le k(N-1)`$.

The basis permutation of A is a coordinate permutation followed by bit
flips. Right composition by its inverse permutes columns of the binary
matrix of F and adds constants. Every output row j above one still has
weight two. Setting its input coordinate $`x_j=0`$ leaves a nonconstant
affine function of the remaining bits, so exactly half of those $`N/2`$
inputs cross to output one. Therefore

```math
r_j(PA^{-1})=N/4\quad(j\ge2),\qquad
R(PA^{-1})\ge(n-1)N/4.
```

Combining the bounds proves the word-length claim. Reordering or
complementing coordinates does not give a bounded or iterated-logarithmic
corner word for all triangular maps. General Clifford mixers, and joint
implementations that amortize a growing word, remain outside this result.


## 7. Other stopped constructions

The Fourier-overlap conversion still has amplitude $`N^{-1/2}`$ and
requires the initialized output word; its repeated-query amplification
already fails the count budget, as recorded in the preceding session.
An involution or sparse-Hamiltonian reformulation must construct its
block encoding under the same workspace contract before it changes that
conclusion.

A local weighted-tree walk has a separate elementary limitation. On
tree vertices let

```math
S|v\rangle=a_v|v0\rangle+b_v|v1\rangle.
```

Each S or adjoint changes depth by at most one. Every word with fewer
than n such factors, even with depth-preserving diagonal operations
between them, has zero root-to-leaf matrix elements. The same is true
of any linear combination of those words. Its action on the root
therefore has distance at least one from the normalized leaf-supported
target. Local polynomial walk synthesis needs degree at least n;
nonlocal descendant operations fall outside this argument and still
need a native implementation. This is not a circuit-depth lower bound
for unrestricted operations.

## 8. Decision and verification boundary

The early checkpoint did not produce a native identity with a complete
sub-square-logarithmic resource and return proof. That session's
construction phase stopped here. The scoped statements in sections 4–6 record
failed conversion architectures; they neither rule out the target nor
improve R79. More subclasses or additional finite searches do not count
as crossing this session's finish line.

The continuing constructive goal is a complete in-place conversion or
a direct recovery word for arbitrary orientations. Discovery can begin
with a partial identity or a reformulation and an explicit missing step;
a finished compiler is the success criterion, not the entry requirement.
To establish an improvement, price changing addresses, encodings, actual
inverses, clean workspace, and full return. An applicable unrestricted
lower bound would also answer the corresponding target question, but
the restrictions here do not provide one.

The initial recommendation to move to independent review absent a new
mechanism is superseded by the user's continuing direction above.
Review of proved results can accompany constructive research. A pause
should preserve the destination and name the next unresolved question.
Final manuscript drafting remains separate.

The [two bounded fixtures](../../tests/test_permutation_bridge_limits.py)
check the exact one-clean four-shear sharpness example on every allowed
small input and the approximate encoded-Clifford counterexample with
explicit linear permutations and isometries. They guard the boundary
between the exact restrictions and the larger open model; they do not
prove the asymptotic statements or supply a scalable compiler. Separate
internal lanes checked the algebra and resource failures. Internal review
is not external peer review. No scientific receipt is regenerated.

The reviewed pull request records full validation, rendered-document
checks, and CI results on the integrated tree.
