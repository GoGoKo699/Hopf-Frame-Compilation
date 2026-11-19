# A shallow sparse word lookup and sparse angular defects in a complete frame

Let a prescribed classical table on p bits have q output bits and s nonzero rows. At each fixed tolerance τ the entire XOR into an arbitrary output has a zero-clean native Clifford+T implementation with
```math
 T,w\le {\rm poly}_\tau(p+q+2)\sqrt{s+1},\qquad
 G\le {\rm poly}_\tau(p+q+2)(s+1),\qquad
 D_T=O_\tau(1+\log_2\log_2(s+2)).
```
All dirty registers and address inputs are included in the full operator approximation; they can be entangled with an external reference. The degree of each polynomial is absolute. The assertion uses the complete fixed-tolerance, dimension-independent-depth zero test from [the gatewise read-only comparator](GATEWISE_READ_ONLY_EQUALITY.md). A fallback with exact logarithmic-in-p T-depth uses the repository's read-only conjunction instead. The improvement here concerns the support dependence; at s comparable to the whole address space it gives no general depth improvement.

## A reversible perfect hash within the count budget

The preprocessing and all circuits below can be found by finite deterministic search. Empty dictionaries use identity. For s≤(p+2)^8 pick a binary linear projection of rank at most k=min(p,ceil(2 log₂(s+1))) injective on the prescribed support. When k&lt;p the expected number of colliding pairs of a random k-bit projection is below one, and a rank-deficient choice may be extended to full rank. Apply an invertible binary change of coordinates exposing this projection in the first k bits. There are at most 4(s+1)²≤poly(p)(s+1) codewords. In particular small s requires no log p term in the T-depth.

For larger s let L=log₂ s, k=ceil L, and partition x=(a,w) with k bits in a and p−k in w. If k=p the identity projection is enough. Otherwise implement the Feistel shear
```math
 a\gets a+\pi P(w),\quad t=2+\lfloor L/(16\log_2(p+2))\rfloor .
```
Embed w linearly in a field of binary degree max(k,p−k); π is a surjective linear projection onto k bits. Choose P among field polynomials of degree less than t. For t distinct values of w, the evaluation map of P is uniform and independent by Vandermonde interpolation. A set of t keys containing two equal w with distinct a can never collide. Hence the expected number of t-subsets of all keys that land in a common first bucket is at most binom(s,t) 2^{-k(t−1)}≤s. Put M=ceil(2t s^{1/t})+t. A bucket of size M would create binom(M,t)≥(M/t)^t>s collisions. Some P therefore gives maximum bucket size less than M. Since t≥L/(16 log₂(p+2)), M=O((p+2)^17).

Every coordinate of πP(w) has Boolean degree at most r=ceil(log₂ t): in a monomial w^j each binary power is Frobenius-linear and j&lt;t. The shear is exact by expanding its coordinates into at most (p+2)^r monomials each. Every target adds the parity of its monomials by toggling disjoint private dirty bits with the corresponding returned read-only conjunctions, accumulating them into that target by CNOTs, uncomputing the conjunctions, and cancelling the initial offsets by the same CNOTs. The parallel conjunction depth is O(log(r+2)). Their T count and private work are polynomial in r per monomial. Since r≤1+t, (p+2)^r≤(p+2)^3 s^{1/16}. All these counts fit the claimed envelope.

For each first bucket a, choose an invertible binary map B_a on w whose first j=min(p−k,ceil(2 log₂(M+1))) coordinates separate its marked words. Pair collision counting gives such a choice just as for the small case. Append a spectator d of the same size as w. Apply B_a to w and B_a^{-1} to d. To see the native price, the block matrices U(A)=[[I,A],[0,I]], L(A^{-1})=[[I,0],[A^{-1},I]] obey U(A)L(A^{-1})U(A)SWAP=diag(A,A^{-1}) in characteristic two. Each addressed matrix shear is a complete operation. Exact R91 reads the (p−k)² entries of the desired A or A^{-1} on a into arbitrary program wires P. Sandwich that program XOR between two copies of the fixed bilinear word adding Pz to the target, and reverse the program XOR, to add the desired matrix product with no P dependence. All (p−k)² Toffolis of the fixed product occupy a constant number of native T layers: Hadamard-transform the target bank and implement all CCZs with private dirty Pauli-phase masks, with the public operands used only as Clifford controls. The matched native inverses are fully charged. R91 with square-root table volume allocation costs O(poly(p)√s) T/width, O(poly(p)s) elementary Cliffords and O(log(k+2)) depth. The spectator may change arbitrarily during hashing. It is acted on by the actual inverse before completion.

The code z formed from a and the first j transformed w bits separates all marked x, with 2^{k+j}≤poly(p)(s+1). All remaining coordinates form a live tail y. Both transformations are invertible, and the transformed address is independent of the spectator. Thus on each code z there is at most one marked tail y_z. For empty bins set y_z arbitrarily and F_z=0, and for other bins set F_z to the marked q-bit payload. Exact lookup of y_z on z is XORed into y. Its cost, and that of the complete F_z XOR, are bounded by R91 at volume 2^{|z|}(p+q). The small-s case uses |z|≤2 log₂(s+1)+1. In the large-s case |z|≤log₂(s)+O(log₂(p+2))=O(log₂(s+2)). Every table lookup hence has depth O(1+log log(s+2)).

## Turning the live comparison into a full word

Write E=[y=0] on the shifted tail. One complete zero test A toggles an arbitrary bit u by E. An exact table lookup B toggles an arbitrary q-bit bank d by F_z. Let C add u d into the desired output by componentwise Toffoli, at constant native T-depth by the same private mask construction. Run, in chronological order,
```math
 C,\ A,\ C,\ B,\ C,\ A^\dagger,\ C,\ B^\dagger.
```
The ideal sum of the four bilinear words is E F_z, with both u and d restored on every input. Replace A by the full fixed-error approximate zero test at error τ/2 and use its actual adjoint. The operator error of the whole construction is at most τ independently of q or s. Reverse the tail shift and the two hashing transformations. All gates, inverses, private circuits and original endpoint registers are included. With an empty tail E=1 the word simplifies to a single B on the actual target.

### Zero-test dependency

In the earlier read-only equality proof a parameter-dependent elementary-matrix walk on SL₃ over a finite binary field acts trivially on a zero parameter. Its nontrivial specializations have the same spectral gap by property (T) of SL₃(F₂[t]); a fixed-tolerance invariant-space filter followed by bounded-clean promotion gives T-depth O_τ(1) and O_τ((p+1)^3) resource counts for a complete, zero-clean [w=0] XOR. The circuit has no T gates targeting w. The short read-only conjunctions, R91 lookup and private bilinear mask implementation here are otherwise given by [read-only conjunctions](../../docs/DIRTY_SUM_COMPRESSION.md) and [the R91 word query](PARALLEL_BILINEAR_WORD_QUERY.md). The potentially large constant from the fixed-tolerance test is absorbed only at fixed τ.

## Nonhomogeneous levels of a complete frame

Start with the positive balanced n-bit real frame. Allow K selected levels, at arbitrary heights, and up to s exceptional prefix nodes per such level, each with any prescribed real rotation replacing its balanced one. At fixed K and fixed frame precision, the one-clean complete frame at square-root dirty width has
```math
 D_T=O_{K,L}(\log\log(n+4)+1+\log\log(s+2)).
```
The additional count and width for the exceptional layers are at most poly_{K,L}(n)√(s+1) and their Clifford count at most poly_{K,L}(n)(s+1). The wide Haar prefix circuits use only polynomial resources in n at fixed precision. In particular if log log(s+2)=o(log n), all resources fit R96 and the depth is o(log n). At the other widths one may retain R96 as a fallback; no uniform improvement on every width is required for this large-width comparison.

For the ordering proof let E_j be the balanced elementary jth layer, U_j=E_j...E_1, and Δ_j be the change of its planar rotation, supported on prefix v, target j and the condition that every original suffix bit past j is zero. Then the actual product is the full U_n times the factors U_j† Δ_j U_j in decreasing j from left to right. A subtlety is that U_j can act nontrivially on some other columns even with a nonzero suffix if earlier layers are defined in a different embedding. In the required conjugation the active suffix-zero sector sees exactly the j-qubit positive Haar, and suffix bits are preserved. On the other sectors the conjugated identity is the identity. Consequently either U_j can be replaced by the j-qubit positive Haar *acting unconditionally on its j data qubits* and its actual inverse. This allows sequential use of the one genuine clean qubit. With approximate Haar prefixes, use sufficiently accurate forward circuits and their native adjoints. The ideal code returns at each boundary, so unitary telescoping includes any actual leakage.

In a Clifford eigenbasis of Y_j, Δ_j is a diagonal phase which differs from identity on at most 2s complete strings (the suffix is zero). Quantize the two signed eigenphases to O_{K,L}(1) binary fractional bits, assigning a common small tolerance to this level. The phase error is the maximum over orthogonal columns, regardless of s. For each bit use the sparse word XOR into the one clean qubit, a synthesized phase on that qubit and the actual inverse. To cancel any untracked one-site phase, apply the actual adjoint of that synthesized phase on the initially zero clean qubit before the computation. The four-operation echo gives the desired literal relative phase on the data. Errors are paid for both XOR occurrences and both one-site words, and the complete logical qubit and all work return in the ideal comparison. The total number of these calls is bounded at fixed K,L. A final arbitrary sparse complex leaf phase can be added by exactly the same diagonal construction.

The unrestricted dense frame has neither the bounded exceptional-level hypothesis nor sparse angle tables. The dictionary result gives no sublogarithmic general scaling for such a frame.
