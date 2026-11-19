# Complete positive Haar from dirty fingerprints

At fixed $`L\ge6`$ and at square-root dirty width, the complete positive balanced real frame has native T-depth $`O_L(\log\log(n+4))`$ with one clean input. The operator estimate includes every wavelet column, the return of every dirty wire and entangled references. The exact balanced construction and its native rotation convention are in [R80](HOMOGENEOUS_MAGNITUDE_DEPTH.md).

## 1. A fingerprint seed can start in an arbitrary quantum state

Let $`z`$ denote computationally preserved inputs and let the desired phase be $`(-1)^{f(z)}`$. There is an exactly implemented Boolean prediction $`g(z,a)`$ with a seed in a register of dimension $`Q=2^h`$. At each fixed $`z`$, at most $`K`$ values of $`a`$ are wrong. The bad seeds can be different for every input. Assume the *controlled* phase for $`g`$ is paid. Section 2 constructs it explicitly for the application.

Use the product eigenbases of $`Z`$, $`X`$, and $`Y`$ on the same $`h`$ seed qubits. Every cross-basis inner product has absolute value $`1/\sqrt Q`$, by multiplying the one-qubit overlaps. Let $`P_i(z)`$ project onto the bad seeds in basis $`i`$. In any two of the chosen bases the bad-bad submatrix has at most $`K^2`$ entries of magnitude $`1/\sqrt Q`$. Its Frobenius norm bounds its operator norm, giving
```math
\|P_j(z)P_i(z)\|\le\mu:=K/\sqrt Q,\qquad j\ne i .
```
This bound applies to an arbitrary seed vector, including a seed entangled with the rest of the input. [Pauli eigenbases and mutual unbiasedness](https://mathworld.wolfram.com/MutuallyUnbiasedBases.html)

Prepare three vote bits in zero. For basis $`i`$, change the seed basis by single-qubit Cliffords, apply a Hadamard to vote $`i`$, apply the predicted phase controlled by that vote, apply the second Hadamard, and restore the seed basis. For fixed $`z`$, the amplitude of recording the wrong bit is the projector $`P_i`$ on the incoming seed. The amplitude of a correct bit is $`I-P_i`$. Record all three votes coherently, apply the Boolean majority sign using three CZ gates, and undo the actual recording circuit. No measurement or random initialization occurs.

The four wrong-majority records have at least two wrong votes. With adjacent wrong votes the branch operator has norm at most $`\mu`$. With the first and third wrong and the second correct it is bounded by
```math
\|P_3(I-P_2)P_1\|\le \mu+\mu^2 .
```
The records are orthogonal. Before uncomputing, the erroneous part has norm at most $`2(\mu+\mu^2)`$. A wrong phase changes a vector by twice its norm. Thus the clean-code error, with the entire seed returned, is at most $`4(\mu+\mu^2)`$. Direct sum gives the same bound on superpositions of inputs and tensoring a reference preserves it. In particular, a single bad seed never becomes an adversarial deterministic input to all three tests.

When $`f(x,y)=y\cdot F(x)`$, Hadamards on the arbitrary output bank $`y`$ before and after this phase give a clean-code approximation of the complete XOR $`y\mapsto y+F(x)`$. The existing bounded-clean promotion converts the three temporary votes to arbitrary borrowed wires with constant overhead. Its hypotheses apply on the entire output and private-work space. For desired zero-clean error $`\tau\le1`$, it suffices to supply clean-code error $`\tau/(8\cdot2^3)`$ and to choose the power of two $`Q`$ such that $`K/\sqrt Q\le\tau/512`$. The forward and uncompute circuits use actual native adjoints. The width used for the old output and the additional bilinear phase grid is charged. [Complete XOR promotion, Sections 3–4](DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)

## 2. Parallel affine zero tests at full joint accuracy

Suppose $`t`$ outputs specify whether each of $`t`$ prescribed vectors $`v_j(x)`$ is zero. Each vector contains at most $`r`$ binary affine forms on the $`p`$ common inputs. Either keep all outputs or take a prescribed binary linear combination of them. We may gate an individual term by a bounded number of extra literal or affine controls. The coordinate matrices and the linear combination are known at compilation.

Over a binary field of order $`Q`$, associate to a vector the polynomial $`v_0+v_1 a+\cdots+v_{r-1}a^{r-1}`$. A zero vector evaluates to zero everywhere. A nonzero vector has at most $`r-1`$ roots. For any fixed input a union bound on root sets shows that at most $`K=t\max(0,r-1)`$ seed values misclassify *any* vector. The same bound covers the erroneous predictions of the Fourier-output sign at every fixed $`y`$. This is the elementary finite-field fingerprint bound. [Polynomial identity testing and fingerprinting](https://www.cs.utexas.edu/~ecprice/courses/randomized/fa19/notes/lec11.pdf)

The emitted prediction circuit enumerates every *constant* $`a_0`$ in the field and every row $`j`$. Its term checks equality of the variable seed to $`a_0`$, checks that the $`h`$ bits of the field value at $`a_0`$ vanish, and checks the vote control and the row's Fourier-output parity. For a fixed $`a_0`$ each field bit is a known affine form in the actual inputs. An optional row guard adds its stated controls. The phase is consequently a product of controlled signs, each on $`2h+O(1)`$ affine literals. Construct the field classically, for example by finding an irreducible binary polynomial of degree $`h`$. This takes no quantum workspace.

For each term, use a private returned conjunction toggle on a dirty bit and a CZ echo with the last literal. The existing conjunction word uses original controls as read-only CNOT controls, private dirty helpers, and $`O(\log(h+2))`$ T-depth. Replace each read of an affine form by successive CNOTs from its original operands into the same private target; a constant offset is a private X. A CZ with a public parity can also be implemented using only private-target operations. Thus all terms run their non-Clifford layers together. Overlapping inputs require additional elementary Cliffords, with no T serialization. At $`r=0,1`$ an exact affine implementation or a constant-size field suffices. [Read-only affine controls](../../docs/DIRTY_SUM_COMPRESSION.md) · [Conjunction costs](NONUNIFORM_DIRTY_INDICATOR.md)

Take $`Q`$ the smallest power of two at least $`\max(2,(512K/\tau)^2)`$. There are $`tQ`$ terms. Their T-count and dirty work are $`O(tQ(h+2)^3)`$ before constant-factor promotion and its additional output work. The elementary Clifford cost is polynomial also in the address and output lengths, accounting for the lengths of the affine forms. The full zero-clean XOR therefore has
```math
D_T=O\!\left(\log\big(2+\log(2+t(r+1)/\tau)\big)\right).
```
All bounds are for the same circuit. This dimension dependence is especially useful for polynomially many long predicates. Feeding it exponentially many row equalities recovers the logarithmic address-depth scale.

There is a separate gain on the earlier sparse-occupation readout. For $`m=2^p`$, predict the phase $`(-1)^{Y_s}`$ by XORing the equality signs for all rows gated by their occupied bits $`Y_v`$. On a word of weight at most $`k`$, only $`k`$ false row tests can contribute root sets, giving at most $`k(p-1)`$ bad seeds even when all $`m`$ rows are implemented. To remove the bounded votes on the *full* space without relying on a promise outside the code, there is an explicit modification. Use the same circuit for a scalar XOR into an auxiliary bit (insert one extra phase control). Its code bound on each allowed $`Y`$ holds for every old/new target and all work. Every native and promotion operation preserves $`Y`$; the Hilbert space restricted to any fixed allowed $`Y`$ is itself a complete XOR problem, including arbitrary dirty inputs. Promote there, with the same emitted native word for all sectors. Extract the phase by conjugating an auxiliary Z and canceling its initial sign. This gives zero-clean depth
```math
O\!\left(\log\big(2+\log(2+(k+1)(p+1)/\tau)\big)\right)
```
on the whole span of the promised words, with $`m\,\mathrm{poly}(p,k,1/\tau)`$ gates and dirty wires. The earlier constant-depth equality construction still wins for fixed $`k,\tau`$. The new bound gives $`O(\log\log(p+4))`$ even for $`k=\exp((\log(p+2))^c)`$ with any fixed $`c>0`$ and fixed tolerance. On a general dense bank the formula gives $`O(\log(p+2))`$.

## 3. The complete last-one tag

For the positive balanced frame every node rotates by $`R_y(\pi/4)=HZ`$. On an input marker whose last one is at position $`j`$, the output leaves the earlier prefix intact, rotates the marker's one, and rotates the zeros following it. On the zero input it rotates all positions. The last-one value can be computed as a single word query. In a Fourier-output basis, its phase exponent is
```math
\sum_{j=1}^h H X_j\,[X_{j+1:h}=0]\,(y\cdot\mathrm{bin}(j)).
```
Here $`H`$ is an activity bit. All $`h`$ suffixes have length at most $`h`$, so Section 2 gives one complete approximate last-one XOR at $`O_L(\log\log(n+4))`$ depth and polynomial resources. This includes arbitrary initial tag and work values.

## One large protected stage followed by an exact small terminal

For large n use r=ceil(log2(n+1)), h=n−r and s=ceil(log2(h+1))≤r. Partition logical bits as upper X and lower Z of r bits. Compute H=[Z=0] into the external clean bit with the exact returned two-helper conjunction, depth O(log(r+2)). Reserve s original lower logical bits as the tag; they are zero only on the active sector. Use the stage of R80 §5,
```math
 F_H,\quad S_H(T),\quad B^\dagger,\quad F_H,\quad B
```
in chronological order, followed by the actual inverse of the H computation. Here B=R^{⊗h} and the selected product S_H(T) applies R to each j with H=1 and T≤j. The tag masks have short read-only address (H,T), hence the exact paid native depth is O(log log(h+4)). The two conditioned involutions implementing R and B use fixed-length native words on different physical data targets. More precisely with C=S H_Had and V=C T C† one has (VZV†)Z=R, including its scalar. All control masks and their inverses are charged. [Protected stage and native masks](HOMOGENEOUS_MAGNITUDE_DEPTH.md)

With ideal F_H the word is exact on every column. On H=1 the input root and markers carry their rightmost-one index t. After the controlled rotations and B†, the suffix starting at t equals |1 0...0> for a nonzero marker, or the entire word is zero for the root. The rotated prefix does not affect t, so the second tag XOR clears the tag and B restores the selected column. On H=0 both tag operations and the selected product are identities, and the two B factors cancel for every occupied tag. The remaining r lower layers form a complete positive Haar frame on Z independent of X. Replacing both F_H occurrences by their certified complete approximants costs at most 2tau in initialized-isometry norm by a unitary hybrid. This includes leakage and possible failure to clear the temporary tag or H. The exact inverse of the H predicate and later operations preserve this error bound.

Use the previously proved exact positive balanced recursion for the r-bit terminal, or its direct constant-size base at bounded n. The polynomial-work positive branch in R80 §§5–8 costs O(log(r+2)) T-depth, independent of requested precision. Its absolute-width constant is absorbed by the polynomial workspace above. Use tau=2^{-L}/4. One large stage and one exact small terminal suffice; no error is summed over log-star many approximate stages. The result follows at polynomial workspace. When b is smaller than the polynomial required, the previous exact polynomial-count direct fallback with two returned helpers costs O(n log(n+2)) T-depth. For fixed L this is covered by O_L(N/b²) outside finitely many dimensions. Its polynomial gate counts fit the displayed square-root and linear bounds there. Increase constants for the finite small dimensions.

The same polynomial-work proof permits a known homogeneous sequence U_j=(V_j Z V_j†) Z with uniformly bounded native depth, T-count and Clifford count for each supplied one-site V_j, by using the exact R80 generalized-frame tag and masks with these actual words. A bare promise of arbitrary positive one-site angles does not supply such words with the required product error. At growing precision the explicit branch with q=Theta((n+1)^4 4^L) still has depth O(log(2+log(n+2)+L)) if its q poly(n+L) dirty and gate costs fit. The exact existing branch can be used independently at higher precision.

For $`L\ge6`$ fixed, $`\eta=2^{-L}`$, $`N=2^n`$, and $`b\ge17(n+13)`$, the same circuit choice yields
```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\log\log(n+4)).
```
For sufficiently large $`n`$ and $`b\ge\sqrt N/n`$ the polynomial work fits and its counts are absorbed by the exponential budgets. For smaller $`b`$ the existing exact balanced bound $`O(N/b^2+\log(n+2))`$ suffices since $`N/b^2>n^2`$. Finitely many remaining dimensions use that same exact compiler. This allocation argument concerns the positive real frame. The wide branch also works whenever its polynomial private-work budget fits directly. More explicitly, for the positive Haar frame its field branch takes $`Q=O(4^L(n+1)^4)`$, with at most $`nQ`$ short conjunctions and $`h_{\rm field}=O(L+\log(n+2))`$ field bits. Its dirty width and T-count are bounded by $`O(4^L(n+1)^5(L+\log(n+2))^3)`$ plus the polynomial costs of the exact stage and terminal. The Clifford count has at most a further polynomial factor in $`n+L`$. The stage depth is $`O(\log(2+L+\log(n+2)))`$. Consequently **even for $`L\le(\log(n+2))^c`$ at any fixed $`c`$, square-root width suffices for doubly logarithmic depth**. The gate and work budgets of the field branch are then $`2^{o(n)}`$. This growing-precision statement uses the exact balanced rotations. The same construction applies to per-level choices of an exact, bounded-length product $`(V_jZV_j^\dagger)Z`$, including independently chosen positive angles $`0,\pi/4,\pi/2`$.

## Why arbitrary node signs are a full lookup problem

For an arbitrary Boolean table f on n−1 bits, change the matrix at each deepest node p from R to (−1)^{f(p)}R. This is a legal balanced pair of signed sine and cosine. All other nodes remain positive balanced. The last layer is therefore (D_f ⊗ I) times the positive last layer, and the entire signed frame equals (D_f ⊗ I) W_+. Thus a signed balanced frame compiler followed on its input side by the actual inverse of the positive Haar circuit supplies an arbitrary lookup sign, including full work return to the sum of the two approximation errors. One may treat the extra unconstrained logical target as borrowed. Conversely R80's diagonal sign gauge writes every signed balanced frame as two arbitrary Boolean signs surrounding the positive frame. Hadamards and an extra target admit conversions to addressed Boolean XOR, for example the sign of g(x,z)=z f(x) on n Boolean inputs. Hence a general signed extension already encompasses arbitrary Boolean-table signs, even when every rotation magnitude is the same. If both supplied gauges are Clifford phases, the positive result applies directly, and known bounded-depth gauges similarly add just their priced depth.



## Dense, succinct sign tables

There is a broader concrete family of signed balanced frames with the same complete depth improvement. Suppose the two sign-gauge Boolean functions have polynomial-size, fixed-depth circuits over unbounded fan-in AND, OR, parity and NOT, with arbitrary nonuniform wiring (the class usually denoted AC⁰[2]/poly). This promise includes many dense tables. It suffices that the prescribed cosine and sine sign functions at each depth separately have such circuits with common fixed size exponent and depth, since the two gauges can then be formed within the same class.

Here is a direct circuit proof. Let $`f`$ have $`p`$ inputs, gate count $`S`$ and depth $`D`$ under this promise. For an OR gate with child values $`b_i`$, choose $`r`$ independent uniform binary linear forms $`L_j`$ in its input list and evaluate
```math
 1+\prod_{j=1}^r (1+L_j(b)).
```
At zero this is exactly zero. At nonzero $`b`$ it errs with probability $`2^{-r}`$. Use De Morgan for AND, ordinary addition for parity and $`1+b`$ for NOT. Substituting these expressions through the circuit gives a distribution on multilinear polynomials over $`\mathbb F_2`$ of degree at most $`\min(p,r^D)`$. For each fixed input, union bound the events that any gate estimate on its true children is wrong. If none is wrong, the composed answer is correct by induction on the circuit DAG. Thus the failure probability for each input is at most $`\delta`$ once $`r=\lceil\log_2(2(S+1)/\delta)\rceil`$.

Choose a fixed native output tolerance $`0\lt \tau\lt 1/4`$. For example let $`\delta=c\tau^2/(p+1)`$ and let $`M`$ be the first power of two at least $`C(p+1)^2/\tau^2`$. For absolute constants take $`C=64/c`$ and make $`c>0`$ sufficiently small for the three-basis correction and bounded-clean promotion constants. Independently sampling $`M`$ polynomials yields expected bad count at most $`\delta M`$ for each input. The multiplicative upper-tail bound (or its exponential-moment proof) makes the probability that count exceeds $`2\delta M`$ at most $`\exp(-\delta M/3)`$. Since $`\delta M\ge64(p+1)`$, the union over $`2^p`$ inputs is less than one. Fix one list with this simultaneous guarantee. The size-$`M`$ selector bad sets then have $`K/\sqrt M\le 2\delta\sqrt M\le 2\sqrt{128c}\,\tau`$. Exhaustive finite search for such a list is allowed; its evaluations and the true Boolean table are finite. No randomized execution, sampled classical bit, QRAM or supplied seed state appears in the emitted circuit.

For every seed label, expand its Boolean polynomial in squarefree monomials. Attach equality on the $`k=\log_2 M`$ dirty computational seed bits to each of its monomials. The parity of these conditioned terms is the selected polynomial, exactly. Each term has arity at most $`k+\min(p,r^D)`$, and the native read-only conjunction and phase echo above permit all of them to run in parallel on separate dirty helpers. The deterministic phase costs
```math
 D_T=O(\log(k+r^D+2)),\qquad
 T,G,w\le M\left(\sum_{i=0}^{\min(p,r^D)}\binom pi\right)\mathrm{poly}(p+k).
```
The same bounds up to absolute constant factors hold for the corrected full sign, including private work return on every input. A phase probe or a joint XOR output can be included before promotion, as above. For fixed $`\tau,D`$ and $`S=p^{O(1)}`$ these estimates give $`O(\log\log(p+4))`$ T-depth and $`\exp(O(\log^{D+1}(p+2)))`$ counts and dirty width. This last quantity is $`2^{o(p)}`$ for the fixed circuit depth. More generally a polynomial-output-length XOR table with succinct bits is handled at once by compiling the parity of the output-character bits times their individual predicates. The Fourier-character trick also accepts an arbitrary occupied output word.

For completeness, let $`e_{d,0}(u),e_{d,1}(u)`$ give the sign bits on the two edges at each depth-$`d`$ prefix $`u`$. For an input leaf $`x`$, the exponent of the output sign gauge is
```math
 \bigoplus_{d=0}^{n-1}e_{d,x_{d+1}}(x_{1:d}).
```
For an input marker with last one at $`t`$, the input sign exponent is the same path parity up to depth $`t-2`$, followed by $`e_{t-1,0}(x_{1:t-1})\oplus e_{t-1,1}(x_{1:t-1})`$. It is zero on the zero word. Choose between the $`n`$ candidate expressions using the exact Boolean tests $`x_t\wedge\bigwedge_{k>t}\neg x_k`$. Their combination and the output formula above have polynomial size and fixed depth using the allowed gates. Apply both gauge phases around the positive Haar circuit, each approximated to a fixed fraction of $`2^{-L}`$. For fixed gate exponent, circuit depth and $`L`$, the large-width compiler has $`O(\log\log(n+4))`$ T-depth and $`2^{\mathrm{polylog}(n)}=2^{o(n)}`$ counts and work. Thus it fits the simultaneous $`\sqrt N`$ T-count and dirty-width envelopes at square-root allocation. It also fits at every $`b\ge\sqrt N/n`$ after the usual fixed-parameter dimension cutoff. On narrower allocations R80's exact signed balanced bound $`O(N/b^2+\log(n+2))`$ reduces to $`O(N/b^2)`$ since there $`N/b^2>n^2`$. Both branches give the claimed all-width tradeoff with constants also depending on the fixed signed circuit depth and size exponent. The per-node promise concerns Boolean computation of the two signs; no small-support assumption is made.
