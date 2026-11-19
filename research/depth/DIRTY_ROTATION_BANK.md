# Parallel real rotations on arbitrary dirty candidates

Suppose each node at depth $`j`$ of a real frame uses the same effectively prescribed rotation. For fixed $`L\ge6`$, $`N=2^n`$, and the real [R80 allocation](HOMOGENEOUS_MAGNITUDE_DEPTH.md), the complete frame admits

```math
D_T=O_L(N/b^2+\log\log(n+4)),\qquad
T=O_L(\sqrt N+N/b),\qquad G=O_L(N).
```

The wide circuit uses $`\exp(\mathrm{polylog}(n+2))`$ gates and dirty wires. One external clean bit records whether the stage is active. A logical suffix of the frame supplies the seed on its active sector. Independent node signs require separate gauge phases.

## One site

Use the convention $`R_y(t)=\exp(-itY)`$. Fix one target $`x`$ and an error parameter $`\epsilon`$. Coherent synthesis supplies an actual circuit $`V`$ on a private target $`q`$ and $`k=O(\log(1/\epsilon))`$ initialized wires. Its clean-code action approximates $`U=R_y(-\theta/4)`$ up to a fixed scalar, including leakage. It has $`O(k)`$ gates and $`O(\log(k+2))`$ circuit depth after constant-size native Toffoli replacements and Clifford basis change. This is the external precision-dependent synthesis input. [Vasconcelos, Theorem 8](https://arxiv.org/html/2609.34659v1#S3.SS4.SSS1)

Let $`P=|0^k\rangle\langle0^k|`$, $`J=U^\dagger ZU`$, and define the physical and ideal reflections

```math
T=V^\dagger(I-2P\otimes|1\rangle\langle1|)V,\qquad
T_0=(I-P)\otimes I+P\otimes J .
```

The rank-one vectors reflected by these words differ in norm by at most $`\epsilon`$ up to a scalar. This follows for the clean adjoint from unitarity and the initialized-isometry estimate for $`V`$. Hence $`\|T-T_0\|\le4\epsilon`$ and its rank is at most two. Write $`\|A\|_2^2=\mathrm{Tr}(A^\dagger A)/\dim A`$. Then $`\|T-T_0\|_2^2\le16\epsilon^2/2^k`$. The central reflection uses an exact returned conjunction of depth $`O(\log(k+2))`$ and polynomial private work.

Give this site $`h=\lceil2^k\rho\rceil`$ candidate registers, $`\rho\ge1`$, and apply their reflections together as $`T_{\rm all}`$. The normalized trace of $`TT_0`$ is real and equals $`1-\|T-T_0\|_2^2/2`$. At the chosen accuracy it is nonnegative. The product of these traces gives
$`\|T_{\rm all}-T_{{\rm all},0}\|_2\le4\epsilon\sqrt{h/2^k}`$.
Let $`S`$ swap $`x`$ with the target of the first candidate whose work word is zero, acting identically if none qualifies. Set

```math
E=T_{\rm all}S T_{\rm all}S,\qquad W=E Z_x E^\dagger Z_x .
```

For the ideal reflections and a good candidate $`j`$, $`E_0=J_xJ_{q_j}`$, so $`W_0=(JZ)^2=R_y(\theta)`$ on the target and identity on the candidates. If no candidate qualifies, $`W_0=I`$. This last sector has relative dimension at most $`e^{-\rho}`$. Telescoping the four reflections and bounding the contribution of that sector gives

```math
\delta:=\|W-R_y(\theta)_x\otimes I\|_2
\le16\epsilon\sqrt{h/2^k}+2e^{-\rho/2}.
```

## From normalized trace to operator error

Let $`\Delta=W-R_y(\theta)\otimes I`$. The complete private Pauli average of $`\Delta^\dagger\Delta`$ is $`B\otimes I`$, with $`B\succeq0`$ and $`\mathrm{Tr}(B)=2\delta^2`$ because the target is one qubit. The average over a finite list can approximate this operator to $`\gamma^2`$ at every one of the $`n`$ sites. Each sample is positive with norm at most four. Matrix Hoeffding gives a failure bound $`2^{w+O(1)}\exp(-cM\gamma^4)`$ for one site with $`w=h(k+1)`$ private bits. Choose a power of two $`M=O((w+\log(2n)+1)/\gamma^4)`$ with sufficiently large constant, take a union bound and find a list by finite search with strict certified margins. The sites may have different lists. [Tropp, matrix tail bounds](https://doi.org/10.1007/s10208-011-9099-z)

Prepare one shared $`\log_2 M`$-bit seed in the uniform superposition. In each seed branch, conjugate $`W`$ by its prescribed private Pauli, then undo the seed Hadamards. The target commutes with these masks. Orthogonality of the seed records gives a squared initialized-isometry error at most $`2\delta^2+\gamma^2`$ for a single site, on **every** dirty and reference state. For the joint bank, telescope with preceding factors ideal and later factors actual. The preceding ideal targets leave the same shared seed and all private wires unchanged, so the errors add across sites.

Choose $`\epsilon=c\eta/n^3`$, $`\rho=C\log((n+2)/\eta)`$ and $`\gamma=c\eta/n`$ for small $`c`$ and large $`C`$. The whole bank has error at most $`\eta`$ up to the separately allocated routing error. For fixed $`\eta`$, $`h,w,M`$ are polynomial in $`n`$ and the seed occupies $`O_\eta(\log(n+2))`$ wires.

## Native routing and masks

The private reflections run in parallel. Write $`b_{ij}`$ for the predicate that candidate $`j`$ is the first with zero work at site $`i`$. Three XOR sweeps implement the conditional swap. First add $`b_{ij}x_i`$ into each candidate target $`q_{ij}`$. Then add $`\bigoplus_j b_{ij}q_{ij}`$ into $`x_i`$ and repeat the first sweep. The predicates are mutually exclusive at each site. In the Fourier basis of each sweep's output, its **whole** sign has a polynomial-size, fixed-depth Boolean circuit over AND, OR, NOT and parity. The [three-basis Boolean-phase construction](POSITIVE_HAAR_FINGERPRINT.md#dense-succinct-sign-tables) therefore gives each sweep fixed full operator error with depth $`O_\eta(\log\log(n+4))`$ and $`\exp(\mathrm{polylog}(n+2))`$ resources. A constant number of sweep replacements enters the physical word.

Implement the Pauli masks exactly. For each seed assignment the Z mask is a sign with exponent equal to its seed test times a prescribed parity of private computational bits across all sites. Separate returned conjunction work permits parallel T layers. Hadamards on the appropriate physical targets give the X masks. Order the X and Z masks consistently and undo the actual product. All affine reads use elementary Cliffords; the masks take polynomially many gates and $`O(\log(2+\log M))`$ T-depth. This step requires no controlled synthesis oracle.

## Completing every frame column

In the [protected R80 stage](HOMOGENEOUS_MAGNITUDE_DEPTH.md#5-one-protected-stage-and-exact-tag-return), reserve $`r=\lceil C_L\log(n+2)\rceil`$ bottom logical bits for the tag and a disjoint twirl seed. Compute the activity $`H=[Z=0]`$ into the sole clean qubit with an exact returned conjunction costing $`O(\log(r+2))`$ depth. The seed is initialized on the active sector and may be arbitrarily occupied on the inactive sector.

Let $`F`$ be the physical bank for the half angles and $`Z(f)=\prod_j Z_j^{f_j}`$. Implement the selected product by $`FZ(f)F^\dagger Z(f)`$. On the active encoded trajectory this costs twice the uniform bank error because $`R_y(a)ZR_y(-a)Z=R_y(2a)`$. For the ordinary product set $`f_j=H`$; for the tag-selected one set $`f_j=H[\mathrm{tag}\le j]`$. Realize $`Z(f)`$ exactly by expanding the short tag into its rows, with one known parity of target bits per row and read-only controls. The guards and tag stay outside $`F`$. On $`H=0`$ the physical bank and its actual inverse cancel for **every** seed state. Approximate complete last-one XORs from the [Haar construction](POSITIVE_HAAR_FINGERPRINT.md) enter a unitary hybrid, which includes leakage and all work return.

Undo the saved predicate and compile the remaining $`r`$-site complete frame with R80 at a fixed fraction of the final tolerance and depth $`O_L(\log(r+2))`$. At large $`n`$ the wide branch fits whenever $`b\ge\sqrt N/n`$. For smaller admissible widths, $`N/b^2>n^2`$ absorbs the additive depth of R80. Use R80 in the finitely many remaining dimensions. The [prefix-selected construction](PREFIX_ROTATION_PROFILES.md) gives the extension to distinct lists below fixed prefixes.
