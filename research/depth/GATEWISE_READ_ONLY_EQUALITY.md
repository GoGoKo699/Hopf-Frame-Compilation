# An address-preserving zero test from elementary matrices

Let a tested register have $`n\ge1`$ computational bits. For each fixed requested $`\tau>0`$ there is a uniform zero-clean Clifford+T circuit applying $`u\gets u+[x=0]`$ with full operator error at most $`\tau`$, including every arbitrary dirty input and its return. Every physical gate leaves every tested computational bit unchanged, and no T gate targets that register. With a straightforward field basis its T-count, elementary Clifford count and dirty width are $`O_\tau(n^3)`$, and its native T-depth is $`O_\tau(1)`$. Thus several such tests can use overlapping sets of tested bits and distinct private wires in the same native layers. An XOR into another test's read-only inputs is outside this parallel schedule.

## A uniformly expanding family that acts identically at zero

Write $`R=\mathbb F_2[t]`$, $`E_{ij}(v)=1+v e_{ij}`$ for $`i\ne j`$, and $`G=\mathrm{SL}_3(R)`$. Since $`R`$ is Euclidean, $`G=\mathrm{EL}_3(R)`$. It has Kazhdan property (T). One sufficient general result is [Ershov and Jaikin-Zapirain, Property (T) for noncommutative universal lattices](https://arxiv.org/abs/0809.4095). Let $`H=\mathrm{SL}_3(\mathbb F_2)`$, and let $`\Gamma`$ be the kernel of evaluation at zero. The finite ordered multiset
```math
 S=\{h E_{ij}(t) h^{-1}:\ h\in H,\ i\ne j\}
```
has 1008 terms, including repeats, each an involution. It generates $`\Gamma`$ as a group. Indeed $`H`$ and the $`E_{ij}(t)`$ generate $`G`$ by $`[E_{ij}(a),E_{jk}(b)]=E_{ik}(ab)`$ for distinct indices and addition in each root group. Move all constant $`H`$ letters in a word to its right end by conjugation. Their product is the evaluation at zero and is identity on the kernel. The index is 168, so $`\Gamma`$ has (T) and the generating multiset has a positive Kazhdan constant $`\kappa`$. Equivalently, induce an invariant-free $`\Gamma`$ representation up to $`G`$ and put a common almost-invariant vector in each of the 168 coset spaces. The constant generators permute the spaces and $`E_{ij}(t)`$ act there by elements of $`S`$. This transfers a positive ambient Kazhdan bound.

Use a prime $`h>n`$ with $`h=O(n+1)`$ and construct $`F=\mathbb F_{2^h}`$ by finite search. Linearly embed the input into a complementary space to $`\mathrm{span}_{\mathbb F_2}\{1\}`$ in $`F`$; denote its value by $`a`$. For every nonzero input $`a\notin\mathbb F_2`$, hence $`\mathbb F_2[a]=F`$ by primeness of the extension degree. The evaluations $`t\mapsto a`$ map $`\Gamma`$ onto $`\mathrm{SL}_3(F)`$: every $`c\in F`$ equals $`a f(a)`$ for a polynomial $`f`$, so the full elementary group lifts via $`E_{ij}(t f(t))\in\Gamma`$. At $`a=0`$ every element evaluates to identity.

Let $`\Gamma`$ act by computational permutations on three $`h`$-wire words $`w\in F^3`$. For a nonzero input the image group has just two orbits, zero and all nonzero vectors, by elementary transitivity. Its full invariant space is the span of $`|0^{3h}\rangle`$ and $`|+^{3h}\rangle`$. Reserve two distinct bits in these words, preparing the first as $`|1\rangle`$ and the second as $`|-\rangle`$. For every other workspace state the resulting subspace is orthogonal to the two invariant vectors. These two bits are later charged as clean and promoted along with the walk selector.

Pad $`S`$ with identity branches to $`2^k`$, with $`k=11`$ so $`2^k\ge 2|S|`$. The lazy uniform average $`A_a`$ is Hermitian with spectrum in $`[0,1]`$. On the orthogonal complement of the invariants, $`\|A_a\|\le1-\kappa^2/2^{k+1}`$ when restricted to that complement. Indeed some generator moves each unit vector there by at least $`\kappa`$, and the Dirichlet form of $`1-A_a`$ is the sum of squared generator displacements divided by $`2^{k+1}`$. At input zero the whole operator is identity. The two-bit prepared subspace need not be invariant under the walk. Orthogonality to the invariant subspace before the walk is enough for the norm bound.

## Literal circuit, return and parallel sharing

On $`k`$ clean selector bits, let $`S_0=H^{\otimes k}\mathrm{SELECT} H^{\otimes k}`$. Compression on selector zero is $`A_a`$. Insert the fixed-point selector phase and use the actual adjoints in the repository's cubing recurrence. After $`j`$ stages the clean selector block is $`A_a^N`$, $`N=3^j`$. Choose $`j`$ so the complement contraction gives the desired fixed tolerance, and synthesize the finitely many selector phase rotations to a sufficiently small full-operator error. At $`a=0`$ the walk is identity and the selector rotations cancel with their actual inverses. Apply $`S_j^\dagger [u\gets u+[\mathrm{sel}=0]] S_j`$ after preparing the two field bits and undo their preparation. For the ideal selector phases the error on a nonzero address and arbitrary remaining field state is at most $`2(1-\kappa^2/2^{k+1})^N`$ on this initialized code. Literal rotation approximation adds at most twice the accumulated forward error. The statement holds coherently since the tested input is block diagonal in every step. The error controls both the selectors and the two field bits and all other returns. [Existing fixed-point calculation, including phase synthesis and the actual controlled XOR](EXPANDER_ZERO_REFLECTION.md)

Each generator costs fixed Clifford field-word changes of coordinates $`h_0,h_0^{-1}`$ followed or preceded as appropriate by the shear $`w_i\gets w_i+a w_j`$. For the conditional generator one may execute the $`h_0^{-1}`$ and $`h_0`$ frames unconditionally bracketing the *selected* shear, since the two cancel on its other branches. In a polynomial basis field multiplication has at most $`O(h^3)`$ coordinate monomials of the form $`x_l w_{j,b}`$ contributing to some target bit $`w_{i,c}`$; pick the input embedding as consecutive basis coordinates excluding the constant basis element so $`a`$ uses original bits directly. Apply Hadamards to word $`w_i`$, realizing the selected shear as a diagonal sum of signed Boolean monomials, each with the same fixed $`k`$ selector literals in addition to $`x_l,w_{j,b},w_{i,c}`$. Here all operands can be shared controls. For clarity, such a phase is a $`d`$-input controlled sign with $`d=k+3`$ fixed. A returned dirty multi-controlled toggle for any constant number of controls is available with T only on private wires, via the read-only-conjunction word. To keep *all* principal operands as controls including the phase target, borrow private $`z`$ and toggle it by the conjunction of a fixed subset; evaluate the phase involving $`z`$ and the remaining controls before and after the toggle, then undo both as in the standard commuting dirty echo. Recursing down to CZ multiplies the count by a fixed constant; equivalently split off one control $`b`$ and difference $`CZ(z,b)`$ before and after toggling $`z`$ by the product of all other controls. Physical Clifford operations affecting shared data are only use as controls (a diagonal CZ with a dirty target can be synthesized by H on the private target and a shared-control CNOT). Give one returned dirty ladder and phase mask per monomial. All such words then align their fixed number of T layers, serializing their Clifford controls. Negative selector literals are handled once around a selected branch. There are constantly many branches. The field word and its conjugating basis frames are private to an instance. Thus each complete SELECT costs $`O(h^3)`$ gates and dirty wires at constant native T-depth and literally preserves the shared address at gate boundaries. [Returned controls with private phase work](../../docs/DIRTY_SUM_COMPRESSION.md#read-only-affine-controls-with-private-phase-work)

Use the existing bounded-clean promotion on the $`k+2=13`$ initialized bits for this complete XOR query. Its sector-selective echoes, toggles of the 13 borrowed bits, and controlled bilinear phase act on private work and the query output, and its only calls on $`x`$ are the native source and its inverse. Set the input source tolerance to a sufficiently small fixed multiple of $`\tau`$ per the promotion bound. The $`2^{13}`$ sector visits are constant for all input sizes. This gives the stated zero-clean result, full return and sharing. Classical generation need not know an optimized numerical Kazhdan constant. For each size the finite field permutation graphs can be enumerated and their gap lower bounded within a constant factor by certified finite eigencalculations; property (T) bounds the resulting required number of filter calls uniformly. [Clean-source promotion](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)

For $`r`$ private output instances in parallel, the per-instance error $`\delta`$ implies at most $`r\delta`$ joint full error by telescoping, allowing shared read-only inputs. This bound alone needs $`\delta=O(\tau/r)`$ to certify the entire product. The constant-depth fixed-tolerance theorem gives no such uniform precision improvement.

## Sparse occupied lookup

Take a $`p`$-bit address and $`m=2^p`$ rows. Suppose the unknown occupied word has weight at most $`K`$, allowing arbitrary coherent superpositions over all such words and entanglement with a reference. Give row $`v`$ a private dirty target $`t_v`$ and an approximate complete comparator $`U_v`$ for $`s=v`$. For the row use
```math
 U_v^\dagger CZ(Y_v,t_v) U_v CZ(Y_v,t_v).
```
At $`Y_v=0`$ the complete row is identity exactly, even if its comparator has error. At $`Y_v=1`$ the row has error at most $`2\delta`$ from the required sign. All comparators are read-only on the same $`s`$; different row offsets are incorporated as affine constants inside the private shear. All output words $`Y`$ are preserved. Thus the full error on the whole promised subspace is at most $`2K\delta`$, including arbitrary dirty return. Taking $`\delta=\tau/(2\max(1,K))`$ yields **constant native depth at each fixed $`K,\tau`$**, independently of the address length and the positions of the occupied sites. It costs $`O_{K,\tau}(m p^3)`$ gates and dirty wires. This is a nontrivial occupied-bank phase on unknown data. Applying Hadamards would change the promise to its Fourier image, so it gives no all-input indicator.

There is also a growing-sparsity bound. Write $`L'=1+\log_2(2\max(1,K)/\tau)`$ for $`0\lt \tau\lt 1`$. The fixed gap needs $`O(L')`$ visits, and the fixed-dimensional selector phase in each visit costs $`O(L')`$ T-depth at the necessary accuracy. The number of promoted clean bits is fixed across precisions. Hence the same parallel word gives $`D_T=O((L')^2)`$, $`T,G=O(m(p^3 L'+(L')^2))`$, and $`O(mp^3)`$ dirty wires. This uses the [uniform single-qubit rotation and precision accounting](EXPANDER_ZERO_REFLECTION.md#8-uniform-precision-and-complete-indicator-scope) with a different constant-size selector. For fixed $`\tau`$, $`\log(K+2)=o(\sqrt{\log(p+2)})`$ suffices for sublogarithmic depth on the promised input. On a general dense bank, per-test error $`\delta`$ alone only gives a row sum of order $`m\delta`$.

For the literal independent clean-row reflection before any promotion, the loss on dense occupations is genuine even with exact fixed-point rotations. The lazy walk above has 1040 identity terms against 1008 involutions, so $`A_a\succeq I/64`$ at every parameter. On a false row, the overlap of the initial clean code with the returned state of the selective sign word on any unit product input $`\psi`$ is $`1-2\|A_a^N\psi\|^2`$. If a uniform error below one is claimed on all occupations, using a single active false row already forces this real overlap to be positive. It is at most $`1-2\cdot64^{-2N}`$. For $`Y`$ all ones, a fixed address and a product of private inputs, the one selected row is exact and the $`m-1`$ false-row overlaps multiply. To keep the literal signed overlap at least $`1-\tau^2/2`$ one needs $`(1-2\cdot64^{-2N})^{m-1}\ge 1-\tau^2/2`$, so $`N=\Omega(\log m)`$ at fixed $`\tau\lt 1`$. This applies to the stated independent walk and its clean-code reflection. Joint filtering, different walks and recompiled promotions require their own analysis.

<a id="sparse-code-preprocessing"></a>

### Reversible preprocessing cannot make an arbitrary bank sparse

The sparse promise above cannot be obtained on the entire input code
by adding dirty wires and applying a clever unitary. This limitation
is a sharp dimension statement, independent of the gate model.

Let Q physical qubits have k initialized inputs, with embedding
$`J:\mathbb C^{2^{Q-k}}\to\mathbb C^{2^Q}`$. All other inputs,
including the address and any extra work, are arbitrary. On $`q\ge1`$ designated
occupation wires let $`P_K`$ select Hamming weight at most K, tensored
with identity on every other wire. Define

```math
B(q,K)=\sum_{j=0}^K\binom qj,
\qquad 0\le K\le q.
```

There exists a unitary U with uniform leakage
$`\|(I-P_K)UJ\|_{\rm op}\lt1`$ if and only if

```math
B(q,K)\ge 2^{q-k}.
```

Indeed $`\mathrm{rank}\,P_K=2^{Q-q}B(q,K)`$. If this is smaller
than the initialized code dimension $`2^{Q-k}`$, the map $`P_KUJ`$
has a unit kernel vector. Its leakage is exactly one, contradicting
the strict bound. Conversely, when the rank condition holds, map an
orthonormal basis of the initialized code into $`\mathrm{im}\,P_K`$
and extend to a unitary; its leakage is exactly zero. This converse
asserts abstract existence, with no native synthesis or efficiency claim.
In particular, allowing small nonzero uniform leakage gives no slack
in the dimension requirement.

Dirty padding multiplies both dimensions equally. Address motion,
intermediate garbage, external references and a fixed unitary change
of occupation basis do not evade the conclusion. With zero clean
inputs it forces $`K=q`$. With one clean input, binomial symmetry
gives the exact threshold $`K\ge\lfloor q/2\rfloor`$: for even q,
the lower half excluding the middle term has dimension strictly less
than half the full bank, whereas for odd q the lower half through
$`(q-1)/2`$ has exactly half the dimension.

More generally, for $`K\le q/2`$ the binomial bound
$`B(q,K)\le2^{qH_2(K/q)}`$ gives

```math
k\ge q\,[1-H_2(K/q)].
```

Here $`H_2`$ is binary entropy. The binomial bound follows by weighting
the sum with success probability $`K/q`$; each term through K has
weight at least $`2^{-qH_2(K/q)}`$. Endpoint $`K=0`$ is immediate.
Thus $`K=o(q)`$ requires at least $`q-o(q)`$ clean inputs. For the
indicator bank $`q=2^p`$ and fixed K, already
$`\log_2 B(2^p,K)=O_K(p)`$ requires $`k\ge2^p-O_K(p)`$.
The growing-sparsity regime above also has $`K=o(q)`$. It cannot
be upgraded to a bounded-clean complete source by whole-code
sparsification. A joint source that never enters this occupation
promise is outside the statement.

## A single joint projector and its limitation

For $`m=2^p`$ rows and selector $`s`$, there is nevertheless a joint constant-accuracy projection filter whose contraction is independent of $`m`$. Choose a prime $`h=\Theta(p+1)`$ large enough that $`(hm^2+2m)2^{-h}\lt 1`$. Search for a binary linear map $`B:\mathbb F_2^p\to F`$ such that all nonzero $`B(v)`$ avoid $`\mathbb F_2`$ and lie in distinct Frobenius orbits. A random binary linear map works with positive probability: the images of two distinct nonzero words are independent uniform field elements, so the union over all ordered pairs and all $`h`$ Frobenius powers costs at most $`hm^2/2^h`$. Include the $`2m/2^h`$ cost of small-field hits. One can equivalently embed first in a linear subspace avoiding 1 and choose an ordinary binary basis for its coordinates. For general $`B`$ implement each address-dependent coefficient as the appropriate XOR of the original address bits inline in the shear monomials, with still polynomial overhead.

Apply each generator in the single SELECT word on every row $`v`$ at parameter $`B(s+v)`$. One row is trivial. All the other parameters have distinct irreducible minimal polynomials coprime to $`t`$. Chinese remaindering $`t f(t)`$ lifts each elementary matrix at one nontrivial row while specializing to zero increment on all other rows. Thus the image is the full product of independent $`\mathrm{SL}_3(F)`$ factors. Each row may also be conditioned on its own read-only occupation bit via the same fixed-degree shear. In a fixed occupation sector the image is again the product over the surviving nontrivial factors. Let $`J_v`$ denote the local fixed-vector projector (identity when its occupation gate is off). They commute and depend on $`Y`$ but not $`s`$. The global invariant projector on arbitrary occupied rows and work is
```math
 P_s=\prod_{v\ne s}J_v .
```
By the same *single* Kazhdan gap its complement contracts uniformly independent of the number of active rows. On the joint invariant space each SELECT branch acts as identity, and the selector rotations in the actual fixed-point word cancel there. The orthogonal complement is invariant throughout, with zero-selector amplitude bounded by the contraction above. Thus the clean-code circuit approximates the full projector-controlled XOR on a private target, acting as $`X`$ on the range of $`P_s`$ and identity on its complement. To apply bounded-clean promotion, diagonalize the commuting $`J_v`$ within each occupation sector for the proof. In this basis the target is a Boolean XOR, and all sector controls and bilinear echoes in the promotion act only on the clean helper or private output bits. No physical change of work basis is needed. For the resulting zero-clean query $`Q`$, the echo $`Q^\dagger Z_u Q Z_u`$ gives the reflection $`I-2P_s`$ within twice the query error. No extra field-preparation bits are required here since both the invariant and orthogonal signals are desired. At fixed tolerance the complexity is $`O(m\,\mathrm{poly}(p))`$ native gates and dirty width and constant native depth. Its reach differs from the parallel per-instance error statement.

The joint filter changes only a thin subspace when a row is removed. Writing $`P_*=\prod_vJ_v`$ gives $`P_s-P_*=(I-J_s)\prod_{v\ne s}J_v`$. These difference projectors are pairwise orthogonal, so their ranks sum to at most the full workspace dimension. In coherent address direct sum the normalized rank of $`R_s-R_*`$ is at most $`1/m`$, where $`R_s=I-2P_s`$. For any word with $`q`$ such reflected calls replace each by its corresponding address-independent complete $`R_*`$. The normalized rank of the difference is at most $`q/m`$ by telescoping, even with arbitrary unitary interlayers, and the normalized Hilbert–Schmidt distance is at most $`2\sqrt{q/m}`$. Compression on $`k'`$ clean bits adds a factor $`2^{k'/2}`$. If the erased word has a fixed normalized Hilbert–Schmidt gap from the desired occupation measurement, a bounded number of these ideal projectors cannot close it. Individually approximated calls also contribute their explicitly accumulated full errors. The estimate includes dependence of the $`J_v`$ on occupied words and is independent of their ranks. It controls this product-projector readout and quantitatively explains why a common gap for the joint walk alone does not improve the frame scaling.

This rank argument uses commuting fixed-space projectors, not tensor-product coordinates. In any fixed unitary representation of a product of independent finite groups, averaging each factor gives such commuting projectors. Changing only that representation retains the same bound when selecting s still removes exactly its factor. A coupled group image or a changed address-dependent representation requires a new analysis.

The [field equality fixtures](../../tests/test_read_only_field_equality.py) enumerate finite parameter orbits, independent CRT lifts and projector rank profiles. The uniform spectral-gap, physical error and sparse-code dimension bounds use the analytical arguments above.
