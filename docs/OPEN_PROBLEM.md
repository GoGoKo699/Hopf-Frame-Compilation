<a id="two-resource-questions-for-complete-hopf-frames"></a>

<a id="endpoint-theorem-and-remaining-depth-question"></a>

<a id="two-clean-count-theorems-and-remaining-depth-question"></a>

<a id="one-clean-count-frontier-and-remaining-depth-question"></a>

<a id="one-clean-count-frontier-and-sublinear-large-width-depth"></a>

# Resource frontiers and the remaining depth question

[Compiler claims](../manuscript/PUBLICATION_SCOPE.md) · [Proof map](README.md) · [Research results](../research/README.md)

Write $`N=2^n`$ and $`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$.
The coherent Clifford+T model
has $`a`$ initialized flags and $`b`$ arbitrary dirty qubits. Its full
initialized-isometry contract includes all logical inputs, dirty/reference
return, literal phases, and clean leakage. The supplied angles and leaf phases
admit terminating certified evaluation, including unequal and singular tuples.

## Current status

| Question | Status and reading route |
|---|---|
| Count at the selected clean/dirty allocations | Resolved by [Results B–C](../manuscript/PUBLICATION_SCOPE.md#four-principal-results); the exact-base and banked ranges remain distinct |
| Count and depth on one circuit | [Result D](#large-width-t-depth), R79, gives absolute all-precision bounds and polylogarithmic depth through $`L=\Theta(\log n)`$ with one clean qubit for both frame families; optimal large-width depth remains open |
| Fixed-accuracy general depth | [R90 dirty catalogue terminal](../research/depth/DIRTY_CATALOGUE_TERMINAL.md) gives $`O_L(N/b^2+\ell)`$, with $`\ell=\log_2(n+2)`$, for arbitrary magnitudes, retaining one clean qubit, the original literal widths, and the count budgets; optimal depth remains open |
| Stronger all-precision research bound | [R96](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md) gives $`O(NL/b^2+L+\ell)`$ for every $`n\ge1,L\ge6`$, with optimal counts, one clean qubit, original literal widths and full error; [R95](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md), [R94](../research/depth/PORTABLE_WORKSPACE_CATALOGUE.md) and [R92](../research/depth/UNIFORM_CATALOGUE_PRECISION.md) retain their independent narrower constructions |
| Supporting depth results | [R80–R85](#large-width-t-depth) retain promised-family and terminal theorems; R87's fine decoder, [R88's sharp-query/fallback interfaces](../research/depth/PREFIX_SELECTED_COIN_TERMINAL.md), and [R89's coarse reductions](../research/depth/CAUSAL_COMPACT_COIN_DEPTH.md) support R90 |
| Constant-depth obstruction | The [zero-clean five-layer bound](#constant-depth-frontier) is proved; [R97](#outer-axis-overlap) excludes bounded outer-axis misalignment; the unrestricted one-clean five-layer bound remains open |
| Query reformulation | [Conditional R98](../research/depth/QUERY_DEPTH_FRAME_TRANSFER.md) transfers a faster restored-address dirty word-query interface to general fixed-accuracy frame depth; its [approximate extension](../research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md) requires only constant source accuracy at fixed frame precision, under uniform simultaneous source budgets. [Two-call corner reductions](../research/depth/BOOLEAN_QUERY_FRAME_REDUCTION.md) separately transfer unrestricted Boolean-query lower bounds |
| Paid scalar predicate | [Expander zero reflection](../research/depth/EXPANDER_ZERO_REFLECTION.md) has zero clean work, full-operator error tau, $`T,G=O(pL+L^2)`$, $`b=O(p)`$ and $`D_T=O(L^2)`$, where $`L=1+\log_2(1/\tau)`$. Its native address motion prevents free shared-row use; the complete indicator remains unpaid |
| Joint-bank source limitation | [Uniform-address screen](../research/depth/EXPANDER_ZERO_REFLECTION.md#9-the-uniform-address-obstruction-to-joint-relabeling) excludes completed joint relabelings with every address-independent bank/helper quantum operation, at compression error below the target amplitude and initialized-unitary error below square root of two. Address-sensitive mechanisms are outside this representation limit |
| Quantitative two-band scope | [Clifford rank–norm bound](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen) requires $`k\ge(2/3)\log_2p-(1/3)\log_2\log_2p-O(1)`$ for exact complete logical action $`K=\alpha I_p`$, with $`\alpha>0`$. Exact Gram flatness alone is insufficient. The approximate threshold shrinks with k, leaving fixed-error growing-clean sources open; this is not an unrestricted depth bound |
| One-pass and whole-echo distinction | [Counter revision](../research/depth/COHERENT_COUNTER_SOURCE.md#9-one-pass-row-filters-clean-dimension-error-and-dirty-mixing) proves a sharp clean-dimension/error tradeoff with dirty mixing. A complete zero-clean Majorana echo pays return but retains linear depth; its separate joint-cycle obstruction assumes the actual Majorana-preserving native grammar |
| Unrestricted two-layer transfer constraint | [Shared native encodings](../research/depth/UNRESTRICTED_TWO_LAYER_QUERY.md) give a single occupied transfer alphabet with one returned clean qubit. Its scalar AND floors and exact affine-distance condition remain useful components; the joint theorem below strengthens the query obstruction beyond scalar witnesses |
| Unrestricted two-layer Boolean degree | [Joint native diagonal rigidity](../research/depth/TWO_LAYER_DIAGONAL_RIGIDITY.md) proves degree at most two for exact Boolean XOR with at most one clean qubit. Every nonquadratic query has full initialized-isometry error at least $`\sin(\pi/8)/(2\sqrt2)`$, uniformly in arity and dirty width, with arbitrary Clifford/address mixing. Both clean-midpoint branches and the six-bit cubic bent query are covered. Complex occupied diagonal action alone need not be C3; no complete-frame frontier changes |
| Unrestricted zero-clean three-layer degree | [Joint native structure](../research/depth/THREE_LAYER_NATIVE_STRUCTURE.md) proves Boolean XOR degree at most four with zero clean qubits, and the same uniform error floor for higher-degree queries. All intermediate Clifford mixing is allowed. The following one-clean ceiling is now proved, but is much less sharp |
| Unrestricted one-clean three-layer degree | [Original-pair collapse and exhaustive coverage](../research/depth/BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md#9-complete-native-three-corollary-with-the-no-pair-constant-retained), sharpened by [one-parity Clifford replacement](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md), give sign degree at most $`D_{\rm np}=D_3+576`$ and nonconstant XOR degree at most $`D_{\rm np}-1`$, uniformly in dirty width and intermediate Clifford mixing. The conservative original-circuit rounding floor is $`2^{-530690}`$. The controlling high-rank estimate remains crude; sharp quantitative degree, higher layers and growing unrestricted depth remain open |
| Exact finite one-clean query | The [paid third-band correction](../research/depth/THREE_LAYER_NATIVE_STRUCTURE.md#8-paying-the-third-band-an-exact-cubic-xor-query) realizes AND3 with 18 T, three layers, 99 elementary Cliffords, one clean and one dirty helper; its actual inverse has 109 Cliffords after reduction. Depth three is optimal at this allocation, while a 16-T/34-Clifford/six-layer baseline has smaller counts. This does not improve the general frame bound |
| Exact quadratic queries | The [complete native construction](../research/depth/QUADRATIC_NATIVE_QUERY.md) gives optimal depth two for every nonaffine quadratic XOR with zero clean and one returned arbitrary dirty helper. Polar rank r gives $`2r+4`$ T gates and $`O(m^2+r)`$ elementary Cliffords on m addresses, including coordinate changes and actual inverses. This one-output construction is not a simultaneous arbitrary-table query |
| Clean-code structural component | [Clifford restriction rigidity](../research/depth/CLEAN_CODE_CLIFFORD_RIGIDITY.md) gives an exact pointwise stabilizer, a literal comparison gap and a common native logical defect space. A correlated encoder counterexample disproves blanket full-space Clifford extension. Unbounded clean pair count occurs for a one-layer complex diagonal target, not a Boolean sign; it is not an unrestricted hardness parameter |
| Paid high-degree algebra source | The [balanced HRT construction](../research/depth/BALANCED_ALGEBRA_PERMUTATION.md) and actual inverse have depth $`O(\log r)`$, zero clean work and full return. The [fast returned subset-product variant](../research/depth/ALGEBRA_QUERY_EXTRACTION.md#5-a-smaller-count-complete-hrt-source) improves T/G count to $`O(r^2 2^r\log r)`$ and extra dirty width to $`O(r^2 2^r)`$ at a larger depth constant. Its original single-table extraction is fully paid but misses R98's simultaneous budgets. The later complete-bank route and fusion pay count/width normalization; sublogarithmic source depth remains unpaid |
| Occupied logical transport | [Full pair transports and paid normalization](../research/depth/OCCUPIED_LOGICAL_TRANSPORT.md) remove quadratic directions and produce globally preserved diagonal pair parities. High clean rank and whole-family agreement still do not force full-space Clifford images, even with no pairs; the exact witness has an easy identity target. This remains a component of the later pair-collapse proof, not an unrestricted hardness parameter |
| One-clean Pauli-boundary three-layer degree | [Restricted native theorem](../research/depth/PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md) gives sign degree at most 2048 and XOR degree at most 2047 when the actual clean image is Pauli at every native-band boundary, allowing arbitrary Clifford/address mixing. The exceptional four-axis relation and retained R74 approximation floor are included. Total clean return alone does not meet this hypothesis; later results cover non-Pauli intermediate codes separately |
| One Pauli interior image at three layers | [Partial-diagonal finite features](../research/depth/PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#7-one-pauli-interior-clean-image-suffices) give sign degree at most $`2^{26}`$ and XOR degree at most $`2^{26}-1`$ if either actual interior clean image is Pauli. The input and final clean stabilizers alone do not satisfy this premise. The stronger 2048/2047 bounds remain when both interior images are Pauli. This supplies one branch of the unrestricted ceiling above |
| Two rank-one non-Pauli interior images | [Exact three-band repacking](../research/depth/PARTIAL_DIAGONAL_THREE_LAYER_FEATURES.md#8-both-non-pauli-images-of-active-rank-one-can-be-repacked) converts this subcase to Pauli boundaries, retaining sign degree at most 2048 and XOR degree at most 2047. This sharpened component remains useful inside the later unrestricted ceiling |
| Non-Pauli code decoding | [Exact one-band decoder](../research/depth/NONPAULI_CODE_DECODING.md) has affine/quadratic control bits and a paid native word, but can induce unbounded Clifford displacement rank. Its zero-pair, all-zero-probe-character subcase, excluding active rank two, has commuting unused-sector defects and a common central Clifford determinant cocycle. That structure does not supply a free native gauge or a general degree bound |
| Fixed endpoint clean-active ranks, including middle pairs | [Finite Pauli-supported insertion](../research/depth/FINITE_INSERTION_THREE_LAYER_FEATURES.md#8-one-clean-qubit-with-bounded-active-ranks-including-middle-pairs) gives sign degree at most $`D_{\rm ins}(1,1,1,k_A+k_C-1)`$ and XOR degree one less for two non-Pauli interior images of ranks k_A,k_C, with arbitrary middle pairs and dirty width. The residual is retained and commutes with each middle axis; it is not a free native gate. The earlier [bounded-relation theorem](../research/depth/BOUNDED_RELATION_NATIVE_MATCHING.md) retains its separate commuting-list statement. Growing ranks are now covered by the later no-pair and pair-collapse theorems |
| Bounded determinant-defect codimension | [Common-cocycle kernel](../research/depth/DETERMINANT_DEFECT_DEGREE.md) identifies a computational-Z subgroup corrected by at most $`Q-1`$ CNOTs at unchanged three-band depth. In the stated zero-pair, all-zero-character branch, codimension r gives sign degree at most $`D_r`$, with degree five for r=0. One common t-dimensional intrinsic central quotient implies $`r\le t(t+3)/2`$. An [actual identity family](../research/depth/CLEAN_TRANSITION_PARAMETERS.md#5-unbounded-determinant-codimension-for-literal-occupied-identity) disproves a uniform r bound under those hypotheses; its first interior image is Pauli |
| Joint nonzero probes and the exceptional rank-two type | [Common native probe algebra](../research/depth/JOINT_NONPAULI_PROBE_ALGEBRA.md) places the generic nonzero-character family in one algebra of symplectic rank at most four when both interior clean images are non-Pauli; this alone does not bound degree. The genuinely exceptional non-S probe type at clean active rank two separates locally and admits a zero-clean three-band logical representation, giving sign/XOR bounds five/four. This sharper five/four bound is specific to that probe type; the finite-insertion result covers other fixed-rank codes |
| Growing-rank clean-phase interface | [Four retained phase factors](../research/depth/RANK_ONE_CLEAN_PHASE_INTERFACE.md) give a no-middle-pair relation profile $`(0,\le1,\le1,0)`$; the two C lists are jointly independent when coupled. A common clean-odd Pauli commuting the actual bands gives a literal zero-clean extension and sign/XOR bounds five/four. An identity family still has unbounded raw displacement change; the later no-pair theorem closes the general coupled branch without bounding that representation parameter |
| Computational phase derivatives and the shared center | [Actual derivative words](../research/depth/COMPUTATIONAL_PHASE_DERIVATIVES.md) give quartic-plus-quotient target factorizations and a shared computational center without bounding raw determinant codimension. The actual no-pair transition reduces to one quadratic-controlled Pauli between three independent native bands; its noncentral common-parity subcase closes at five/four, including endpoint rank two. Whole-center restriction can still create unbounded relations. Later no-pair closure and the specifically injective pair-scalar quotient avoid assuming free fiber matching |
| Global quadratic-control matching | [Full native argument](../research/depth/QUADRATIC_CONTROL_NATIVE_MATCHING.md), within the retained no-pair normal form, closes fixed control polar rank by actual Clifford absorption and finite insertion. The [central-parity refinement](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md) gives sign degree at most 2049 at polar rank two: each of two sectors has an actual Clifford replacement and Pauli clean boundaries. For a genuine quadratic control, the associated Clifford probe images preserve the middle native span and agree on a common invariant isotropic algebra. Paid normalization of an enlarged actual center gives commuting monomial Clifford endpoint images. A global decoder or relative-center-phase bound remains unsupported |
| Constant monomial translations | [Actual sublist peeling](../research/depth/CONSTANT_TRANSLATION_NATIVE_CLOSURE.md) proves diagonal C3, hence sign/XOR bounds three/two, for arbitrary independent endpoint ranks and an actual commuting C3 middle when all endpoint translations are constant. This is a restricted degree theorem, not a new same-depth implementation. Later no-pair results handle center-dependent translations with a weaker uniform ceiling |
| Empty native middle in the quadratic-control form | The [uniform corollary](../research/depth/CONSTANT_TRANSLATION_NATIVE_CLOSURE.md) closes all control ranks: polar rank at least four forces constant translations and degree three; affine control gives two-band degree three; rank two now has the [central-parity bound](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md) 2049. For genuine quadratic control, equal endpoint total-axis products force this empty-middle case. A remaining candidate in this normal form must have a nonempty middle band |
| Exceptional middle-defect support | [Partial-support closure](../research/depth/MIDDLE_DEFECT_SUPPORT_CLOSURE.md) gives sign/XOR bounds five/four when the control-phase commutators are constant, even if other middle translations depend on the center. After removing actual permuted middle-axis orbits, the Boolean-normalized cubic symbol lies in the exterior cube of the center–exceptional-middle overlap W. Fixed $`d=\dim W`$ gives sign bound $`\max\{5,D_{\mathrm{ins}}(0,0,0,d+3)\}`$, uniformly in other ranks. The later native-orbit theorem also closes growing overlap in the actual no-pair reduction |
| Actual clean-active intersection | The [original-transition Gram identity](../research/depth/MIDDLE_DEFECT_SUPPORT_CLOSURE.md#8-the-original-clean-transition-removes-the-mixed-case) closes fixed $`t=\dim(L_A\cap L_C)`$ in the no-pair branch, with sign bound $`\max\{5,D_{\mathrm{ins}}(0,0,0,\max\{t,2\}+1)\}`$, even when both active ranks grow. At control rank at least four, a surviving nonconstant probe forces the quadratic support to equal the common active intersection, inside the positive middle span, and its displacement to lie outside the original active algebra. The enlarged center permits that displacement; no unrestricted closure follows |
| External-center correction and its native cost | [Simultaneous native preimage](../research/depth/EXTERNAL_CENTER_MIDDLE_CORRECTION.md) replaces all flags of the controlled Pauli by constant translations, retaining a Boolean cubic correction of the middle band. Rank-one actions have two odd-multiplicity quotient axis types. A specified growing corrected-middle family admits no one-T-layer realization with arbitrary Clifford endpoints or returned dirty helpers; this does not close the full three-layer degree question |
| Bounded total middle displacement | [Actual native factorization](../research/depth/EXTERNAL_CENTER_MIDDLE_CORRECTION.md#7-bounded-total-displacement-gives-a-uniform-degree-bound) closes the high-rank no-pair branch for fixed $`d=\dim\sum_u\mathrm{im}N_u`$. With $`k=d(d+1)/2+d+1`$, sign degree is at most $`D_k+72\,2^k`$; low control rank retains finite insertion. The proof keeps the original middle and a two-native-band residual. The subsequent native-orbit theorem closes growing d; original physical pairs are handled separately by Boolean pair collapse |
| Growing total displacement in the actual no-pair reduction | [Native-orbit closure](../research/depth/GROWING_DISPLACEMENT_NATIVE_ORBITS.md) forces actual logical pair permutations plus one exceptional transvection from the full cubic tensor. Retaining all pair phases in an actual three-band residual gives high-rank sign bound $`D_3+576`$. The [rank-two refinement](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md) makes this the universal no-pair ceiling, independent of active ranks and dirty width. Fixing the remaining translation parity creates an actual-list relation, so constant translation alone does not justify the earlier independent-axis peeling; that is the next structural question. Original clean-transition pairs are a distinct, closed branch |
| Stronger centers and the decoder limit | [Center audit](../research/depth/NATIVE_CENTER_SATURATION_LIMITS.md) adds genuine characteristic and whole-band parity centers. A fixed Clifford decodes the retained native algebra only when its band already centralizes that algebra. An [actual no-shared-flip CCZ source](../research/depth/ACTUAL_NATIVE_CENTER_WITNESS.md) retains center-dependent translations after saturation; safe-axis derivative pruning retains a nonzero parity. These are representation limits on an easy target, not unrestricted lower bounds |
| Paid Majorana Walsh source | [Native source audit](../research/depth/MAJORANA_WALSH_QUERY.md) gives a complete independent-Pauli probe with $`T=2mp`$, $`G=O(m^2)`$, depth at most $`4p`$, and arbitrary carrier width $`m/2`$ for $`m=2^p`$. Sharp separate-source bounds $`mp/2`$ and p apply only to Majorana-preserving native synthesis. A flat-source equality construction leads to the exact readout obstruction below; the separate flat-source depth obstruction remains constant. No R98 or frame improvement follows |
| Exact standalone dual readout | [Walsh channel-denominator theorem](../research/depth/DUAL_FRAME_READOUT_BARRIER.md) gives $`T(Q)\ge2p\binom{p-1}{\lfloor(p-1)/2\rfloor}=\Omega(m\sqrt p)`$ for $`m=2^p`$, allowing arbitrary Clifford interlayers and returned clean/dirty helpers. A paid six-call Q-only wrapper gives the complete indicator. The bound is exact and operation-specific: R91 already gives the fused indicator O(m) T-count, and an approximate readout is not excluded |
| Approximate readout transfer and scoped barrier | [Full error and resource ledger](../research/depth/APPROXIMATE_QUERY_FRAME_TRANSFER.md) uses at most 24 readout calls per completed word query and constantly many substituted frame queries. Fixed frame precision therefore needs only constant-in-address readout accuracy. Small constant error still forces logarithmic depth and superlinear count for a separately synthesized Majorana-preserving Gaussian conjugator; unrestricted approximate readout and joint indicators remain outside this obstruction |
| Approximate predicate interface limits | [Complete-indicator audit](../research/depth/APPROXIMATE_INDICATOR_FILTER_LIMITS.md) proves a sharp $`\sqrt2`$ operator-error floor for additive dirty-output phases with a Hamming-distance polynomial of degree below the address length. A short coherent hash needs initialized randomness for its single-predicate guarantee and gives error two for the whole bucket-based indicator even with clean randomness. These are limitations of the stated interfaces, not unrestricted approximate-query lower bounds |
| Phase-signal and coherent-frame interfaces | [Return and correction audit](../research/depth/PHASE_SIGNAL_QUERY_LIMITS.md) gives a fixed-space obstruction for products of rank-one signal reflections on arbitrary dirty work. Full coherent Pauli-frame transport yields a complete indicator in two source calls with linear Clifford/dirty overhead and factor-two error. These are scoped interface statements, not unrestricted approximate depth bounds |
| Joint approximate source interfaces | [Maximal gauge](../research/depth/WEAKER_INDICATOR_INTERFACES.md) and [selected-swap](../research/depth/JOINT_ROUTING_INDICATOR_INTERFACES.md) reductions retain arbitrary intermediate address/work motion, actual inverses and constant error factors. These source freedoms remain constant-call equivalent to the complete indicator. Fixed Clifford simultaneous-diagonalization failure is a representation limit |
| Direct output-alias source | [Four-call complete readout](../research/depth/OUTPUT_ALIASED_CHARACTER_SOURCE.md) uses the actual output wires as character coordinates and permits arbitrary quantum private-helper fibers. Complete error is at most $`4\eta`$; fusion uses 32 calls. A sublogarithmic weak source or whole extraction remains unpaid; the known source retains logarithmic depth |
| Clean-sector scheduling audit | [Boundary restriction and black-box rank](../research/depth/CLEAN_SECTOR_PARALLELISM_AUDIT.md): a clean bank fixed at every actual Clifford boundary is dispensable without increasing T-count/depth. Universal black-box promotion requires at least $`2^k`$ source occurrences at the stated literal errors; this is not a native-depth bound |
| Complete bank-moving components | [Bank difference](../research/depth/COMPLETE_BANK_DIFFERENCE.md) gives a literal nonlinear vector merge and full return; its balanced equality-bank schedule has linear count/width and $`O(b^2)`$ depth. [Fermionic Walsh routing](../research/depth/BANK_MOTION_WALSH_INDICATOR.md) gives a complete indicator with $`T\le4b2^b`$, $`D_T\le4b`$, $`G=O(b4^b)`$ and no extra work. Both retain the old depth cost rather than improving R96 |
| Whole-word Majorana native-depth limit | [Fused-word theorem](../research/depth/WHOLE_WORD_MAJORANA_LIMITS.md): for $`b\ge2`$, zero clean work and complete phase-indicator error $`\epsilon\lt1`$, $`2^d\ge2^b(1-\epsilon^2/2)^2`$ when every actual Clifford boundary preserves the address and, sectorwise, the standard Majorana span of all bank/helper wires. Thus $`d\ge b-1`$, and $`d\ge b`$ below $`\sqrt{2-\sqrt2}`$. This survives whole-word fusion in that grammar, not unrestricted Clifford interference or address motion |
| Noncommuting dirty geometric phases | [Full-space loop and native cost audit](../research/depth/DIRTY_GEOMETRIC_PHASE_LOOPS.md) gives a common winding sign on both signal eigenstates, but independent row eigenphase errors must sum in full-bank operator norm. Its paid digital recurrence misses the target count/depth envelope; no unrestricted approximate obstruction follows |
| Native sharing for approximate indicators | [Incidence and block bounds](../research/depth/NATIVE_INDICATOR_SHARING.md): at full-operator error below one half, each bank observable needs p active T occurrences in an unrestricted zero-clean native word. Separately emitted b address-preserving disjoint bank blocks cost at least bp T gates; no unrestricted growing-depth conclusion follows |
| Table-specific source gauges | [Paid syndrome interface](../research/depth/TABLE_SPECIFIC_INDICATOR_INTERFACES.md) allows arbitrary motion inside fibers of the table's probed syndrome. Full row rank recovers the complete indicator gauge; no low-rank promise or faster general R98 compiler is supplied |
| Bounded-clean and one-signal indicator sources | [Paid clean removal](../research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md) converts complete bounded-clean queries to zero-clean words, including leakage and arbitrary off-code action. A one-signal block near the complete query at amplitude one half or inverse square root of two suffices through 24 source calls. The faster source remains missing; growing clean banks retain exponential generic call overhead |
| Polynomial-overhead complete indicator fusion | [Paid root splitter and exact suffix](../research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md) absorb polynomial source T-count/dirty overhead and polynomial overhead over quadratic Clifford count into the original linear-budget indicator. Eight source/inverse calls give error at most eight times the source error and add only $`O_\tau(\log\log(p+4))`$ depth. A sublogarithmic complete source, including at a suitable slowly growing clean budget, remains unconstructed |
| Dirty-tree differential compilation | [Full branch-difference recurrence](../research/depth/DIRTY_TREE_DIFFERENTIAL_RECURRENCE.md) gives $`e_i\le e_{i-1}+\delta_i`$ and complete error at most $`2\sum_i\delta_i`$. Common private action and untouched-parent errors cancel; independent sibling errors can add. Actual native scheduling, inverses and all final return remain charged |
| Robust indicator count and incidence limit | [Pauli-nullity calculation](../research/depth/INDICATOR_PAULI_NULLITY.md) gives $`T\ge m-1`$ for $`p\ge2`$ and full-operator error below $`1/\sqrt2`$, with unrestricted native mixing and dirty width. A one-layer word saturates the earlier incidence criterion while failing the target; no growing-depth bound follows |
| Joint Fourier-word screens | [Robust loop and block spectra](../research/depth/JOINT_FOURIER_WORD_LIMITS.md) retain their finite-size bounds: one source band is excluded at both paid amplitudes for $`p\ge3`$, and two-band half-amplitude compression has floor $`\sin(\pi/8)/(4\sqrt2)`$. The general fixed-clean two-band interface is now covered by the separate source theorem below |
| Bounded-clean two-band source gap | [Amplitude-free Clifford-sum theorem](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#7-an-amplitude-free-defect-gives-a-uniform-relative-gap) excludes complete indicator blocks at sufficiently small fixed relative error and sufficiently large p, uniformly over every positive amplitude and every uniformly bounded clean count. Arbitrary dirty width, common phase and rejected action are allowed; flatness is unnecessary. The explicit thresholds remain conservative. A nonzero one-clean compression has norm at least $`\sin(\pi/8)`$ even at finite size. Unbounded clean count, additional bands and unrestricted depth remain open |
| Dressed overlap and simultaneous parity screens | [Quarter-turn probe screen](../research/depth/OVERLAPPING_TWO_BAND_SOURCES.md#6-a-clifford-invariant-coefficient-screen-for-quarter-turn-transitions) survives independent logical Clifford dressing and spectators, including the repaired involutory family. [Simultaneous output parities](../research/depth/SIMULTANEOUS_OUTPUT_PARITY_LIMITS.md) give a Reed–Muller restriction and linear native rank requirement; a cheap nonindicator saturates the rank condition, so it is not sufficient |
| Original physical pairs under occupied Booleanity | [Common-center reduction](../research/depth/OCCUPIED_BOOLEAN_PAIR_CENTER.md) gives transition quotient dimension at most $`3m+2`$ and residual at most $`3m+1`$ for m physical pairs and nonempty endpoint active lists. Its transported center and retained native independence support the later pair-collapse proof. Quotient dimension alone still does not bound total Pauli support |
| Original-pair profiles without bounded silent span | [Central-module profiles](../research/depth/CENTRAL_MODULE_NATIVE_PROFILES.md) bound changes of three exact native profiles along every actual Boolean matching action and bound each action's escaping span modulo the genuine center. These bounds are uniform in additional silent span; they do not bound the union over actions, classify inside-span mixing, or prove fixed-pair degree closure |
| Minimal actual-pair completion | [Original one-pair source](../research/depth/MINIMAL_PHYSICAL_PAIR_COMPLETION.md) produces a tied central partner after singleton cancellation. For the fixed active input/middle word, any admissible silent input band moving that partner permits no one-band output completion to a Clifford target, even with arbitrary unused-sector freedom. Every complete target in this fixed family has a zero-clean three-band word with one middle axis: diagonal targets are diagonal third-level, with sign/XOR degree at most three/two. This sharper component remains distinct from general paired recognition |
| Original pairs with rank-two endpoints | [Paired rank-two theorem](../research/depth/MINIMAL_PHYSICAL_PAIR_COMPLETION.md#5-all-paired-sources-with-two-clean-active-axes-at-both-endpoints) gives sign degree at most six and XOR degree at most five for actual one-clean three-band words with two clean-active axes at each endpoint and a genuine original middle pair. Arbitrary silent lists, dirty width and unused-sector action are allowed. The same-circuit rounding theorem gives higher-degree error floor $`2^{-530690}`$, allowing a common phase. The no-pair branch retains its separate theorem; the following growing-rank extensions have additional physical hypotheses |
| Larger physical-pair sector cost | [Larger-support analysis](../research/depth/MINIMAL_PHYSICAL_PAIR_COMPLETION.md#8-larger-physical-support-a-two-band-sector-encoder-and-a-one-band-obstruction) proves that a projected endpoint code can require two bands, with isometry error floor $`1/(4\sqrt2)`$ against one-band encoders. Its paid repair gives five bands with surviving middle singletons or four in the pair-only witness. That representation does not supply three-band stability; Boolean pair collapse bypasses it. The encoder obstruction is not a high-degree Boolean counterfamily or unrestricted depth bound |
| Retained growing-rank pair components | [One-sided rank two](../research/depth/MINIMAL_PHYSICAL_PAIR_COMPLETION.md#6-one-rank-two-endpoint-permits-the-other-active-rank-to-grow) gives the historical sign/XOR bound 197/196. [Actual two-axis pair differences](../research/depth/MINIMAL_PHYSICAL_PAIR_COMPLETION.md#7-two-axis-pair-differences-with-both-endpoint-ranks-growing) give the historical $`D_{\rm np}+1`$/$`D_{\rm np}`$ bound with both ranks growing. Boolean pair collapse now supplies six/five throughout these classes, with the same original-circuit rounding floor $`2^{-530690}`$ |
| Retained all-active original-pair recognition | [Original-list theorem](../research/depth/ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md) gives the earlier sign/XOR bounds 2428/2427 when either endpoint's silent list commutes with the entire original middle. Its physical recognition and matching remain useful components; the uniform structural theorem below now sharpens this subclass as well |
| Retained pair-only center component | [Enlarged transported center](../research/depth/PAIR_ONLY_ENLARGED_CENTER.md) permits both coupled silent lists with no middle singletons. Its original bound is superseded by the uniform all-active closure below; its physical-center construction remains useful |
| Retained singleton-coupling component | [Actual-center extension](../research/depth/ORIGINAL_PAIR_SINGLETON_COUPLING.md) gives equal endpoint coupling ranks h and residual dimension at most $`2h+4`$. The following theorem removes bounded h without asserting that h is bounded |
| Uniform all-active genuine one-pair closure | [Full-space diagonal reduction](../research/depth/ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md) gives sign/XOR degree at most six/five, sharpened to five/four at common active rank at least four, with arbitrary original singletons, both silent lists and dirty width. The pair difference must sum every actual active axis at both endpoints. Both-end diagonal sublist peeling and one fixed Clifford reduction replace the earlier generic feature bound. The retained higher-degree floor is $`2^{-530690}`$ under these original-list hypotheses. The following recognition theorem removes the support premise in the positive-singleton class |
| Retained positive-singleton recognition | [Boolean support recognition](../research/depth/POSITIVE_SINGLETON_PAIR_RECOGNITION.md) first derived all-active support for one genuine pair with no negative singleton after literal repair. Its Boolean derivative contradiction is retained by the stronger theorem below; singleton positivity is no longer a premise of the paired-branch closure |
| Complete original-pair collapse | [Scalar-sector Boolean reality](../research/depth/BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md) proves that every complete Boolean one-clean native-three word with any genuine original pair has exactly one, all-active at both actual endpoints, with no negative singleton after literal repair. Thus every paired source has sign/XOR degree at most six/five, or five/four at common active rank at least four. The safe quotient preserves the actual silent list; no product axis is supplied as a free native gate. The separate no-pair bound supplies the unrestricted ceiling above |
| General symplectic-carrier source count and probe interface | [Arbitrary-module readout barrier](../research/depth/SYMPLECTIC_MODULE_QUERY_BARRIER.md) forces $`2^{d-1}`$ carrier qubits and T gates per separately realized source for a degree-d linear symplectic readout. Universal fixed-family probes contain a complete dirty indicator through a paid conversion. A uniform faster independent-probe compiler with the stated simultaneous budgets would supply the full R98 interface for all output lengths and dirty allocations; that faster source remains unconstructed |
| Scoped returned-query obstructions | [Affine Gaussian scalar loops](../research/depth/GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md) have Boolean degree at most $`4r+2`$ after r rounds; [two read-only T layers](../research/depth/READ_ONLY_TWO_LAYER_RIGIDITY.md) force degree at most two and quadratic polar rank at most two with one clean qubit. Neither is an unrestricted dimension-growing T-depth lower bound |
| Three-layer query frontier | [Three read-only T layers](../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md) force affine XOR outputs with zero clean qubits and degree at most four with one clean qubit, at error below $`\sqrt2/6`$ independently of dirty width. [Exact one-clean two-layer read-only queries are affine](../research/depth/READ_ONLY_TWO_LAYER_RIGIDITY.md#5-exact-one-clean-queries-are-affine). Within this read-only interface, two-input AND depth is four with zero clean and three with one clean; the one-clean optimality is exact-only. This does not settle general depth |
| Four-layer zero-clean query closure | [Four read-only T layers](../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md#8-four-layers-with-zero-clean-qubits-degree-at-most-two) force degree at most two for every component of an arbitrary word-valued XOR query at error below $`\sqrt2/12`$, independently of dirty width; literal Boolean scalar phases have degree at most three. The [word extension](../research/depth/READ_ONLY_THREE_LAYER_RIGIDITY.md#87-word-valued-queries-the-outer-shear-is-additive) uses native outer-shear additivity and keeps all outputs arbitrary. It gives no one-clean four-layer or unrestricted-address conclusion |

The [fresh-workspace checkpoint](../WORKSPACE.md) and
[current revision](../research/depth/PARALLEL_PROMOTION_AND_KERNEL_REVISION_RESULT.md)
keep the complete approximate dirty indicator as the main constructive
question. [Clean-source promotion](../research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)
now permits a uniformly bounded clean bank or a one-signal complete
constant-amplitude block as the initial mechanism. The conversion pays
all actual inverses, off-code cancellation, leakage, return, and the
restored-address R98 transfer. The [fixed Clifford-sum gap](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md)
now stops the general two-band block at uniformly bounded clean count,
at sufficiently small fixed relative accuracy and large p, uniformly in
every positive normalization. Decreasing amplitude alone does not reopen
this route. The defect $`K^\dagger P K-PK^\dagger K`$ removes amplitude
from its fixed Clifford coefficients; it is an invariant, not a supplied
linear-combination operation. The [overlap examples](../research/depth/OVERLAPPING_TWO_BAND_SOURCES.md)
retain their finite-size lessons: involutivity can be repaired, but the
dressed quarter-turn family still fails a complete indicator probe.

The [generic one-signal recurrence](../research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md)
pays twenty-four source/actual-inverse calls, with complete promoted error
at most $`16(3\sqrt2\epsilon+2\epsilon^2)`$ and depth at most 104
for a three-band source. Exact signal/toggle return needs separately
proved flatness; otherwise their return errors remain included. The first
unpaid operation in that route is the complete logical transition.
No three-band native source has been supplied, and adding a band alone
does not improve a resource bound.

The revised constructive question is the
[polynomial-overhead complete source](../research/depth/POLYNOMIAL_OVERHEAD_INDICATOR_FUSION.md).
Use it only on $`p-O_\tau(\log(p+2))`$ address bits, pay its
controlled root splitter, and complete the logarithmic remainder exactly.
Polynomial overhead is thereby compatible with linear final T-count and
dirty width. This reopens polynomial-cost row mechanisms only if they
also supply complete-bank accuracy and genuinely sublogarithmic depth;
the retained digital pulse does not. A slowly growing clean bank can fit
the new recurrence, but the source budgets must hold uniformly at the
promoted precision and its exponential call multiplier still costs depth.

The [distance-profile and native scheduling audit](../research/depth/SELECTIVE_INVERSION_DIRTY_INDICATOR.md#51-a-complete-distance-certificate-and-the-native-scheduling-gap)
retains a full-bank error certificate. For b address bits, packing the
unchanged address--signal rotations requires exactly b independent-axis
batches. The subsequent [joint resynthesis](../research/depth/SELECTIVE_INVERSION_DIRTY_INDICATOR.md#52-joint-resynthesis-pays-one-special-angle-detuning-call)
implements the whole fixed-angle bank in at most 11 T layers, with linear
T-count, quadratic Clifford count and no additional work. Thus the packing
limit does not even lower-bound the depth of that same operator. The
component is not a complete selective source. Its period-eight ambiguity
excludes completed independent weight-only rows at that angle; joint bank
mixing remains outside this premise. A new selective composition must pay
its full-operator error and return certificate as well as native depth.
The [joint-decoder revision](../research/depth/SELECTIVE_INVERSION_DIRTY_INDICATOR.md#53-joint-clifford-decoding-of-the-relative-call-still-stops)
now excludes normalized joint Clifford decoding and gives a separate
degree-four limit when the entire bank remains read-only. Those proofs
retain their exact grammar; arbitrary bank motion is outside. The
[independent finer-angle theorem](../research/depth/JOINT_DETUNING_PRECISION_LIMITS.md)
retains at least logarithmic depth under its complete-layer precision
promise. Joint non-Clifford synthesis and whole-word cancellation remain
open, together with the actual selective action and return.

The [parallel-affine revision](../research/depth/PARALLEL_AFFINE_DETUNING.md)
pays polynomially many affine full calls in four T layers with polynomial
address overhead, exact added-helper return and literal inverse costs.
A single commuting batch has error at least $`1/4`$ for $`b\ge3`$.
A fixed-work-Pauli echo with actual Clifford interlayers preserving
address and bank needs $`h\ge\lceil b/2\rceil`$ rounds for
$`b\ge2`$ below error $`\sqrt2`$, independently of parallel
access rank or work width. This is a round-grammar theorem, not an
unrestricted native-depth bound. The first unpaid operation is complete
selectivity with return and qualifying depth; a new mechanism must leave
this grammar or resynthesize the whole word. Its full resource/error/return
recurrence is now explicit.

The [bank-difference merge](../research/depth/COMPLETE_BANK_DIFFERENCE.md)
and [fermionic Walsh indicator](../research/depth/BANK_MOTION_WALSH_INDICATOR.md)
now pay actual bank motion, literal phases, full return and all native
operations. Their displayed complete depths are respectively
$`O(b^2)`$ and $`O(b)`$, leaving the retained logarithmic-depth
source stronger. The [whole-word Majorana theorem](../research/depth/WHOLE_WORD_MAJORANA_LIMITS.md)
also prevents shallow fusion within its actual native grammar: every
Clifford boundary must preserve x and, sectorwise, the standard Majorana
span of all bank and arbitrary-helper wires. Under that premise the
zero-clean complete phase indicator needs $`d\ge b-1`$ at error
below one, and $`d\ge b`$ below $`\sqrt{2-\sqrt2}`$.
For an XOR indicator, the actual output-Hadamard-conjugated word must
satisfy the grammar at both endpoints too. The next candidate must
leave this grammar through non-Majorana Clifford interference, address
motion or another complete interface. The first unpaid operation
remains the complete source with sublogarithmic depth and all resource
budgets; bank motion or whole-word fusion alone does not supply it.

The [direct output-alias revision](../research/depth/OUTPUT_ALIASED_CHARACTER_SOURCE.md)
now reduces the shared HRT extraction to four source/actual-adjoint
occurrences. After a subset-zeta transform on the arbitrary output Y,
all its wires are carrier coordinates, including the top probe on
$`Y_\varnothing`$. Only $`m-1`$ additional carrier wires are needed,
where $`m=2^b`$. A sufficient source is
$`V_0=\sum_c|f(c)\rangle\langle c|\otimes S_c`$, with arbitrary
input-dependent **quantum unitary** action on private dirty helpers.
The completed character cancels this action; source helpers need not
return between calls. With full-operator source error eta, the complete
zero-clean recurrence is
$`T_I\le4T_V`$, $`D_T(I)\le4d_V`$,
$`G_I\le4G_V+bm+2b+2m+2`$,
$`h_I\le m-1+h_V`$, and $`\varepsilon_I\le4\eta`$.
Every address, output, dirty and reference input is retained, with actual
adjoints and complete return error. Root fusion uses 32 source calls,
error $`32\eta`$, and depth
$`32d_V(b_0+1)+O_\tau(\log\log(p+4))`$; request
$`\eta=\tau/32`$. Fixed-L frame transfer requests
$`\eta_L=2^{-L}/(512K_*)`$ under the retained literal allocations.
The known source still has logarithmic depth. These interface and
constant improvements do not change the general upper frontier.

The [observable specialization](../research/depth/OBSERVABLE_TRANSPORT_SOURCE.md)
applies the retained selected-routing principle directly to the output
bank: $`J^\dagger Z_{Y_0}J`$ must approximate its complete addressed
Pauli at error delta. The indicator uses two J calls, no extra fixed
carrier, Clifford cost at most $`2G_J+2m+2`$ and complete error delta.
Fusion uses sixteen J calls and error $`8\delta`$. This relaxes the
prescribed HRT permutation, without creating a faster source by itself.
The new shared-key theorem excludes arbitrarily long Clifford-fiber
words with affine address-Pauli insertion at error below $`\sqrt2`$.
An actual opposite-direction quadratic word leaves that representation,
but its readout has error two from arbitrary dirty-key dependence.
Removing that dependence with qualifying native resources is unpaid.
Neither statement bounds unrestricted depth.

The [current key-mixing revision](../research/depth/KEY_MIXING_OBSERVABLE_SCREENS.md)
separates a final-boundary obstruction valid for arbitrary quantum
key-preserving fibers from more restricted affine-label and moving-key
screens. Their exact word failures do not cover unrestricted native
interference. The [nilpotent-pencil candidate](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md)
leaves the conjugated-translation family: its inverse column is the
selected output, and a fixed Vandermonde basis exposes a unique zero
pivot. The forward quadratic XOR query is paid in eight T layers with
admissible polynomial overhead. The in-place or single-observable
readout remains unpaid. The prescribed weak completion and the complete
indicator have paid conversions in both directions, retaining the
supplied depth order; this is a reformulation, not a faster source.
Fixed-address Clifford stages and completed low-degree monomial rounds
have separate scoped limitations. Genuine non-Clifford interference,
intermediate address motion and joint native resynthesis remain open.

The [general-k amplitude recurrence](../research/depth/DIRTY_PROMOTION_OF_CLEAN_QUERIES.md#63-recovering-a-general-clean-bank-compression-with-one-extra-flag)
also pays recovery at $`1/\sqrt2`$ followed by dirty promotion,
including all zero-tests, actual inverses, helper widths and return
error. For nonnegative beta,gamma, $`k\sim\beta\log\log b`$ and source depth
$`O((\log b)^\gamma)`$, $`\beta+\gamma\lt1`$ is a sufficient
sublogarithmic window only when accuracy and resource envelopes are
uniform along the schedule. No source is supplied, and exact flatness
remains a separate condition for exact signal return. The bounded-clean
two-band obstruction is not an obstruction for all growing k.

The [Gaussian compression theorem](../research/depth/GAUSSIAN_COMPRESSION_SOURCE_LIMITS.md)
handles arbitrary fixed bank Cliffords and all dirty inputs. Complete
amplitude recovery at inverse square root of two forces rank-one clean
overlap. The normalized logical labels and signs have bounded degree
in the physical tableau, giving a query-stage requirement proportional
to p divided by k+1 for the literal affine-query grammar. This is not a
native-depth bound or an obstruction to arbitrary quantum address motion.
The [involutory Gauss theorem](../research/depth/INVOLUTORY_GAUSS_SOURCE_LIMITS.md)
separately screens global scalar interference with output-bank variables
among the coefficients. For degree-D matrix entries on the full external
cube, its unsquared Boolean Gauss sign has degree at most 9D. Paying the
in-place correction cannot rescue that bounded-degree source. Non-involutory
forms and different whole-word completions retain their separate unpaid
operations; the general frontier is unchanged.

The [native parallel-promotion revision](../research/depth/PARALLEL_CLEAN_PROMOTION.md)
now pays its refined packing through exact single-target permutation
synthesis and dirty-indicator parity echoes. It has a sufficient window
$`2^k(k+1)^2=o(\log(p+2))`$, with subpolynomial address overhead.
For $`k\sim\beta\log\log p`$ and source depth
$`O((\log p)^\gamma)`$, the conditions are separately
$`\beta\lt1`$ and $`\gamma\lt1`$, rather than their sum being
less than one. Exactly sector-preserving sources retain complete error
twice source error, independent of the number of copies. General
leakage retains its separate shrinking accuracy requirement. The
fully substituted compression/recovery/promotion/fusion ledger is paid;
a qualifying shallow source remains absent. At k of order log p,
this generic packing is too costly.

A new [MSS spectral-paving wrapper](../research/depth/MSS_SPECTRAL_PROMOTION.md)
pays 64 global native frames and their whole-space gap. Truncating the
batch to at most two active clean keys retains constant requested
source accuracy under exact sector invariance. Its depth is
$`O_\tau(d_V+k+1)`$ and its extra count/width is
$`2^{O(k2^k)}`$. Thus $`k=o(\log(p+2))`$ is a paid sufficient
window with polynomial source overhead, subject to the same native
address-sharing and uniform source budget hypotheses. General leakage
still requires accuracy of order two to the minus k. Finite exhaustive
paving search is allowed; efficient preprocessing is not claimed.

A [latent-center readout interface](../research/depth/NILPOTENT_PENCIL_OBSERVABLE_SOURCE.md#latent-center-readout)
pays Clifford generators whose center contains the addressed bank
observable. Ordinary group mixing does not extract that sign: the
central eigenspaces are preserved but are not pointwise fixed. A
fixed bounded probe also has a shrinking coefficient after conjugation
averaging. Full-space central-sign readout at constant amplitude is
still unpaid and is not excluded by these scoped tests.

The [seven-selector filter](../research/depth/MSS_SPECTRAL_PROMOTION.md#8-a-lazy-seven-selector-filter-improves-uniform-precision)
replaces Hamiltonian simulation by a lazy average and general-block
cubing. At $`L=1+\log_2(1/\tau)`$ its pre-fusion depth is
$`O(L(d_V+k+1)+L^2)`$, with separate count/width recurrences on the
same word. Exact sector preservation requests source accuracy
proportional to tau/L; leakage adds the factor two to the k.
This improves the uniform wrapper cost, while all source quantities
must still be charged at their actual requested precision.

The [Clifford rank--norm source screen](../research/depth/FIXED_CLIFFORD_SUM_SOURCE_GAP.md#10-rank-norm-bounds-sharpen-the-exact-and-finite-error-source-screen)
excludes exact two-band sources throughout this window: they need
$`k\ge(2/3)\log_2p-(1/3)\log_2\log_2p-O(1)`$ clean bits.
For one clean bit, its elementary bound excludes exact sources for
p≥18 and relative error at most $`2^{-10}`$ for p≥20, regardless of
positive amplitude, count or dirty width. Its growing-k approximate
floor still shrinks. The remaining two-band candidate must therefore
be genuinely approximate with growing clean allocation; exact flatness
of its Gram block does not imply exact logical action. Higher-band
mechanisms and unrestricted depth remain separate questions.

The earlier two-clean spectral filter's Clifford-frame class retains
its own conditional gap question. The new frames include paid
non-Clifford Boolean phases. The
[frame limits](../research/depth/SPECTRAL_PROMOTION_FRAME_LIMITS.md)
exclude common-Lagrangian frames, sufficiently small physical support
cones, and stated graph/trace-pencil families. Those ordinary block
and representation restrictions do not lower-bound unrestricted T-depth.
The [odd-key counter identity](../research/depth/COHERENT_COUNTER_SOURCE.md)
needs no precision gradient: one clean bit selects all odd Fourier keys
with arbitrary dirty lower bits, and every key gives the same returned
zero test. A literal flat source follows with logarithmically many
clean counter wires. Joint modular arithmetic and the full-bank native
schedule remain unpaid at the required depth; that clean count
is still outside the new sublogarithmic-clean sufficient window.
These components do not yet compose to improve the unrestricted frontier.
The [crossed-kernel continuation](../research/depth/CROSSED_CUBIC_KERNEL_SOURCE.md)
gives an exact full collision criterion and a paid shallow trial, but
its observable fails on arbitrary dirty sectors. The separate diagonal,
plain-XOR and unbalanced-translation obstructions retain their stated
middle-action and outer-kernel hypotheses. Other whole-word dynamics
remain outside; the complete unrestricted frontier is unchanged.

The [coherent-correction interface](../research/depth/OUTPUT_ALIASED_CHARACTER_SOURCE.md#3-the-arbitrary-key-correction-has-the-same-weaker-return-interface)
permits the same arbitrary private fibers. Its principal action is
$`R(z,y)=(z,y+f(z)+yf(z))`$. Direct readout gives
$`T_I\le4T_R`$, $`D_T(I)\le4d_R`$,
$`G_I\le4G_R+bm+4b+2m+2`$,
$`h_I\le3m-2+h_R`$, and error at most $`4\eta`$.
The paid quadratic-inverse XOR converts R to the paired transform
$`F\otimes F^\dagger`$ and back with at most 80 added T layers,
$`2T_{Q_k}`$ added T gates and $`2G_{Q_k}+3(2m-1)`$ added
Cliffords. The paired transform is a sufficient weak source, counting
its second bank as private work. It is not a standalone returned F.

The [algebra-key group operation](../research/depth/OUTPUT_ALIASED_CHARACTER_SOURCE.md#6-a-paid-algebra-key-multiplication-does-not-yet-convert-the-key)
is paid in constant depth with polynomial overhead, but its key a still
needs the conversion $`a=f(z)`$. A paid phase query for the group
Fourier kernel does not supply the corresponding single-register
transform. The [explicit kernel and scalarization audit](../research/depth/OUTPUT_ALIASED_CHARACTER_SOURCE.md#7-the-shallow-fourier-kernel-is-not-a-returned-scalar-source)
exposes a trace character with scalar amplitude $`2^{-w}`$, for
$`w=2^{b+1}-1`$, while its full operator is distance $`1+2^{-w}`$
from that scalar. The particular Pauli-Z twirl repairs it using $`2w`$
clean labels. Any K-mask combination at relative error epsilon below one
needs $`K\ge(1-\epsilon^2)2^{2w}`$. Standard promotion or
amplification retains exponential cost. This is a limitation of the
mask-twirl representation, not an unrestricted native-depth bound.

The first unpaid upper operation is consequently a sublogarithmic
native weak fiber source, paired transform, coherent correction, or
complete extraction within the polynomial-overhead budgets. The existing
inverse-degree argument covers completed classical stages, including
classical helper garbage; it does not cover arbitrary quantum fibers or
nonmonomial native interference. A new mechanism must pay the actual
coherent conversion, not infer it from a quadratic XOR, fixed-key
Clifford membership, or trace access. The historical
[eight-call extraction](../research/depth/FULL_BANK_CHARACTER_TRANSPORT.md)
remains a valid independent construction.

The [clean-sector parallelism audit](../research/depth/CLEAN_SECTOR_PARALLELISM_AUDIT.md)
rejects a different apparent shortcut. If every actual Clifford boundary
fixes the clean register and no T gate targets it, restriction already
removes that clean bank without increasing T-count or depth. Merely
returning the clean bank at source endpoints does not allow overlapping
native calls to run in parallel. A universal zero-clean black-box wrapper
needs at least $`2^k`$ source occurrences for k clean bits at the stated
literal errors; this counts calls, not native layers or unrestricted depth.

The source screens retain their architecture and parameter scope. The
unrestricted indicator nullity result concerns count; the incidence
theorem gives no growing depth bound. R96 and the complete-frame frontier
are unchanged until the missing shallow source is constructed.

The retained [native-orbit theorem](../research/depth/GROWING_DISPLACEMENT_NATIVE_ORBITS.md)
closes growing total displacement in the actual no-middle-pair reduction.
The [occupied-center theorem](../research/depth/OCCUPIED_BOOLEAN_PAIR_CENTER.md)
keeps original physical pairs and reduces their noncentral support.
Bounded central features do not justify substituting quotient dimension
for total support in the finite-insertion theorem. Its Section 8 now
proves a degree bound at fixed original pair count **and** fixed
additional coupled silent-axis span, with common feature coordinates.
Bounding that additional span, or a sufficient dressed-support invariant,
from pair count alone was not supplied by this method. Section 9 retains a bounded
scalar-feature ledger without this span hypothesis and gives a bounded
sum of individually clean-preserving logical Cliffords from one shared
actual silent band. Their displacement maps have pairwise zero products,
which does not bound their ranks. The [central-module profiles](../research/depth/CENTRAL_MODULE_NATIVE_PROFILES.md)
now constrain every actual matching action without bounded silent span,
including its escaping-axis span modulo the genuine center. The original
exceptional-axis span is at most 16c modulo that center and can be absorbed
into insertion dimension at most 17c. Those profiles alone do not supply
simultaneous matching inside the old middle span or propagation of this
exceptional module. The later Boolean pair-collapse theorem bypasses
that representation question and closes every positive original-pair
count at depth three. Higher-layer asymptotics remain open.

The [all-active original-pair theorem](../research/depth/ALL_ACTIVE_ORIGINAL_PAIR_CLOSURE.md)
retains physical recognition and matching when either endpoint's silent
list commutes with the whole original middle. Its earlier sign/XOR
bounds 2428/2427 are superseded by the structural theorem below. The
[pair-only center proof](../research/depth/PAIR_ONLY_ENLARGED_CENTER.md) and
[h-dependent singleton extension](../research/depth/ORIGINAL_PAIR_SINGLETON_COUPLING.md)
retain useful transported-center components. Their old question of
unbounded singleton coupling is now closed within this all-active class.

The [uniform all-active diagonal reduction](../research/depth/ALL_ACTIVE_PAIR_DIAGONAL_REDUCTION.md)
allows arbitrary original positive singletons, both coupled silent lists,
unbounded endpoint active ranks and arbitrary dirty width. For exactly
one genuine pair whose difference is the sum of every actual active
axis at both endpoints, it gives sign/XOR degree at most
**six/five**, with **five/four at common active rank at least four**,
without any bound on the singleton-to-silent coupling rank h. The
[retained scalar-branch argument](../research/depth/ALL_ACTIVE_SINGLETON_CLOSURE.md)
supplies native matching and a single computational X mask for the actual
dressed axis. Exact target transport gives the same property at the other
endpoint. Actual independent native sublists coupled to that axis can
then be peeled as diagonal C3 phases. One fixed Clifford change makes
the remaining insertion and actual three-band residual diagonal on the
full logical space, so the sharp three-layer phase-degree theorem pays
the complete target. This replaces the generic 37,748,865/37,748,864
matrix-feature estimate; it uses no restriction to a $`P_h`$ eigenspace,
replacement native band or unused-clean-block extension.

The proof preserves literal phases, the physical clean code, arbitrary
dirty/reference inputs and the unrestricted unused clean block. The
retained original-circuit rounding theorem gives the higher-degree
floor $`2^{-530690}`$ when the resulting actual original lists meet
this class's hypotheses. The retained rank-two six/five theorem supplies
the remaining case; the new bound also sharpens the earlier one-sided
2428/2427 estimate.

The [positive-singleton recognition theorem](../research/depth/POSITIVE_SINGLETON_PAIR_RECOGNITION.md)
then removes the all-active support premise when exactly one original
pair remains and every singleton is positive under the actual clean
comparison after literal pair-orientation repair. A proper-support pair
would make the transported axis commute on the scalar parity branch and
anticommute on the quarter-turn branch. The complete Boolean derivative
forces one global commutation sign, excluding that branch change. The
pair difference must therefore sum all actual active axes at both
endpoints, so the six/five and rank-at-least-four five/four bounds apply.
This is an analytic recognition result, with the same full operator
contract and retained rounding floor. Its next focused question was the
negative-singleton case, where a common active span cannot be assumed.

The [original-pair collapse theorem](../research/depth/BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md)
now resolves that question without resynthesizing product-axis phases.
The joint scalar sector of all original pairs preserves independence of
the actual output-silent native list: an active-axis witness excludes
the clean label, and full endpoint-list independence excludes any pair
difference. On that full remaining-register Boolean block, the one-band
reality lemma removes every silent coupling to a clean-even active
direction. Its physical commutation sign survives restriction. The
global Boolean derivative can then have only one commutation sign across
all pair sectors. Every pair difference must equal the entire active
sum, so original middle independence permits at most one. Time reversal
gives the output conclusion, and the literal clean comparison removes
negative singletons.

Thus every paired source has sign/XOR degree at most six/five, with
five/four at common active rank at least four. Together with the retained
no-pair and Pauli-interior branches, the [exhaustive unrestricted corollary](../research/depth/BOOLEAN_ORIGINAL_PAIR_COLLAPSE.md#9-complete-native-three-corollary-with-the-no-pair-constant-retained),
sharpened by the [central-parity reduction](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md),
gives the at-most-one-clean three-layer bound

```math
D_{\rm np}=D_3+576,
\qquad \deg g\le D_{\rm np},\quad \deg f\le D_{\rm np}-1.
```

For polar-rank-two control, one genuine computational-parity restriction
makes the insertion an actual Clifford on each whole remaining-register
sector. The Pauli-boundary bound gives sector degree 2048; reconstruction
in the same coordinates gives 2049. This removes the old finite-insertion
term from the maximum. The retained high-rank bound dominates 2049,
the Pauli-interior ceiling and the paired ceiling.

The remaining estimate is still very crude. At high rank, fixing the
last translation parity makes the dressed Pauli's X mask constant but
introduces a relation in the actual endpoint lists. The independent-axis
peeling argument cannot then be applied, and the insertion still has
quadratic control. The [literal relation factorization](../research/depth/CENTRAL_PARITY_QUADRATIC_REDUCTION.md#7-the-exceptional-relation-recreates-a-quadratic-controlled-pauli)
now retains the phase and exposes another quadratic-controlled Pauli,
whose controls need not be in the whole word's shared center. Moving
it through the actual middle remains unsupported. This identifies
a limitation of that proposed proof step, not a new unrestricted
obstruction. The original-circuit rounding theorem
gives the corresponding higher-degree error floor $`2^{-530690}`$,
including leakage, arbitrary dirty/reference inputs and the allowed
common phase. It does not cover arbitrary fixed tolerance. The new
obstruction applies to the actual unrestricted one-clean three-layer
Boolean model, but still gives only a constant-depth lower bound.
Higher layers and dimension-growing unrestricted depth remain open.

The [previous clean-transition checkpoint](../research/depth/CLEAN_TRANSITION_PHASE_SESSION_RESULT.md)
closes the fixed endpoint clean-active-rank case with arbitrary middle
pairs. The finite-insertion proof keeps the entire exceptional operator,
including non-C3 and nonisotropic support, and uses its commutation with
each middle axis. Uniformity as active ranks grow is now supplied by the
separate no-pair closure and original-pair collapse, not by freely
deleting the insertion.

Two [actual easy-target families](../research/depth/CLEAN_TRANSITION_PARAMETERS.md)
rule out proposed uniform representation bounds. Degree-three CCZ
products have unbounded residual parity-list relation dimension with
both interior images non-Pauli. An identity target has unbounded
determinant codimension in the all-zero-character branch, with one
Pauli interior image. These do not give high-degree native queries.
The no-pair four-factor interface therefore retained the original
phase lists; the later native-orbit theorem bounds their completed
Boolean action without bounding those raw representation parameters.

The [previous joint-transport checkpoint](../research/depth/JOINT_NONPAULI_TRANSPORT_SESSION_RESULT.md)
retains the fixed-codimension theorem, common probe algebra, sharper
exceptional closure and complete conditional R98 transfer. The faster
complete probe or indicator source remains unconstructed. A
[dense Walsh phase attempt](../research/depth/SYMPLECTIC_MODULE_QUERY_BARRIER.md#71-a-dense-walsh-phase-program-does-not-supply-the-faster-source)
fails its native gate and simultaneous resource budgets; its parity-only
repair obstruction does not cover general Clifford interference.
R96 remains the strongest general bound and R79 remains selected Result D.
The qualitative one-clean three-layer degree question is now closed.
Its small quantitative ceiling and complete-frame optimal depth remain open.

The [focused four-band physical-pair scout](../research/depth/FOUR_BAND_ORIGINAL_PAIR_SCOUT.md)
retains actual independent signed lists and the physical clean code.
Its conditional radical certificate and literal shared-pair repair do
not imply pair-count collapse: identity words can have arbitrarily
many pairs. Disjoint-pair Boolean examples, including a nondiagonal
controlled-SWAP middle, still admit other whole-word two-band reductions.
For the precise fixed endpoints and equal-difference one-pair lists,
complete Booleanity first gives a controlled SWAP/CZ classification.
The actual-axis argument then proves that all middle axes commute
across the bands and closes the occupied sign degree at three. It does
not turn their dependent union into a free native layer. The active
bounded lower question must now leave these fixed-endpoint/orbit
hypotheses; general four-band return does not imply the common-J
Pauli matching. No general four-layer degree ceiling or unrestricted
depth obstruction follows.
The [physical-center revision](../research/depth/FOUR_BAND_CENTER_TRANSPORT.md)
now derives center preservation from complete occupied Booleanity in this
same fixed-endpoint class. Literal orientation repair and a sole unmatched
anticommuting orbit also reduce there. A signed native controlled-SWAP
repartition proves that the initial intermediate can nevertheless be
non-Clifford; its existing two-band realization remains. Multiple unmatched
factors are not closed, and a center-moving trial must change the stated
endpoint or diagonal-logical-center hypotheses.

The [general native-probe revision](../research/depth/FOUR_BAND_CENTER_TRANSPORT.md#g-a-sole-unmatched-anticommuting-factor-must-merge-natively)
now pays an isolated anticommuting factor's actual Clifford merger,
even with unequal probes and without commuting the remaining native
phase through the input probe. Its
[active-overlap theorem](../research/depth/FOUR_BAND_CENTER_TRANSPORT.md#i3-native-defects-close-the-one-center-cross-commutation-branch)
bounds first-band/probe cross-commutation rank by shared native
overlap dimension on the affected first-band axes. If the only
intersection
is a genuine matched pair's center, that cross-commutation vanishes
for every size of its actual D expression and outside D list.
This does not imply individual axis matching, Booleanity or cheaper
resynthesis. Commuting unmatched images and coupled active overlaps
remain open. The rank-four matching counterexample retains an explicit
cheaper two-band word.

The center proof now gives a paid error-two-epsilon derivative for
an error-epsilon source and an actual two-band resynthesis of the
fully matched all-active support-two mechanism. The latter preserves
literal phase, unused-block action and return, so it retains the old
cost. The [higher-support continuation](../research/depth/HIGHER_SUPPORT_PAIR_BOOLEANITY.md)
proves the local resynthesis boundary and rejects a five-wire
center-moving clean-return trial through its Pauli-preserving domain.
It cannot be made diagonal by endpoint Cliffords, even with spectators.
For the unchanged middle on these same five wires, every one-band
output endpoint refactorization is also excluded. Different middle
dynamics and new native quantum spectator couplings remain open.
Neither a directional derivative nor this fixed trial gives a total
degree theorem for unrestricted four-band sources.

The [preceding revision](../research/depth/DIAGONAL_ENCODING_SESSION_RESULT.md)
excludes every nonquadratic two-layer Boolean XOR with at most one clean
qubit, including the six-bit cubic bent query missed by scalar transfer
coefficients. The proof uses the full preserved Pauli algebra, explicit
clean removal in the non-Pauli branch, and signed quotient relations and
native row cancellation in the Pauli branch. A complex diagonal non-C3
counterexample remains valid: its input-dependent phase cannot be deleted.
The apparent formal C3 correction to a real quartic sign cannot fit the
same native two-layer allocation. Its next lower question was whether
three unrestricted layers impose a Boolean degree bound independent of
arity and dirty width. A faster restored-address query remains the
parallel upper target; the paid bent-query conjugation gives no improved
R98 primitive. R96, R79, R97 and the complete-frame frontier are unchanged.

The [preceding joint-transfer revision](../research/depth/JOINT_TRANSFER_SESSION_RESULT.md)
proves the unrestricted single occupied transfer alphabet, robust AND3
floor and exact affine-distance condition. Its rational-target Galois
invariant retains the one-clean scope. These remain useful components;
the new diagonal theorem goes beyond their scalar witnesses and resolves
the unrestricted one-clean bent-query exclusion. No dimension-growing
depth bound is inferred from this constant-layer result.

The [three-layer continuation](../research/depth/THREE_LAYER_NATIVE_SESSION_RESULT.md)
resolves the unrestricted zero-clean three-layer degree question with
degree at most four for XOR predicates and a fixed error floor. It also
pays the earlier complex phase correction in a third native band,
giving an optimal-depth finite AND3 query. The one-clean extension was
not automatic: occupied logical preservation need not give
full-space Clifford images, and a shared Pauli clean flip need not exist
for the displayed split. The later no-pair and original-pair proofs
handle that intermediate-code freedom without assuming such an extension.
The growing-arity quartic-sign construction retains an unpaid address
phase in its larger query conversion. R96 and the complete-frame lower
frontier are unchanged.

The [preceding clean-code and algebra revision](../research/depth/CLEAN_CODE_ALGEBRA_SESSION_RESULT.md)
pays two formerly missing constructive operations: an optimal-depth
quadratic Boolean query with one dirty helper, and the in-place algebra
multiplier needed for a balanced HRT source. The proposed quartic-phase
echo loses to the new cubic baseline at the same depth and allocation.
The lower analysis proves rigidity for actual Clifford restrictions and
localizes native probes, but a new exact counterexample shows why a
two-band logical image need not extend as that Clifford on the unused
sector. A growing clean-pair count also survives a one-layer full-target
simplification; that target has complex T phases, not a Boolean sign.
The next lower question at that checkpoint therefore concerned the complete
occupied logical family modulo unused-sector freedom, not representation
parameters alone. The current no-pair theorem and Boolean pair collapse
resolve its qualitative degree test. Neither component construction changes R96 or supplies the missing
arbitrary-table query within R98's complete budgets. The refined multiplier
uses six completed calls in five stages, giving source depth at most
$`40\lceil\log_2r\rceil`$ with the same asymptotic count and dirty
allocation. The [Boolean rounding corollary](../research/depth/CLEAN_CODE_CLIFFORD_RIGIDITY.md#63-boolean-approximation-reduces-to-the-exact-occupied-action)
of R74 reduces sufficiently accurate one-clean three-layer Boolean
approximation to exact occupied action. The current uniform degree bound
therefore inherits its explicit error floor with no separate
approximation-closure gap at that radius.

The [preceding occupied-transport revision](../research/depth/OCCUPIED_TRANSPORT_SESSION_RESULT.md)
replaces the failed extension with full operator Pauli transports and a
paid common normal form. The first remaining lower inference would
discard sector-dependent Clifford phases or apply one-relation matching
after conditioning many pair parities; pair-free silent bands and
non-Pauli unused action also remain. Its constructive branch pays an
arbitrary-table scalar extraction and compiles fast subset convolution
with returned dirty masks, improving structured-source counts without
improving the general query interface. The inverse-degree and
coefficient-quotient limitations proved there concern specific
arithmetic representations, not unrestricted native circuits.

The same revision closes the actual
[Pauli-boundary subcase](../research/depth/PAULI_BOUNDARY_THREE_LAYER_RIGIDITY.md).
Ordinary pair peeling and the one-relation exception leave at most five
quantum core qubits. Exact endpoint Clifford extraction and a finite set
of phase features give a width-independent degree bound, including all
central controls and literal phases. At that checkpoint the next
unrestricted question was the non-Pauli intermediate clean code; assuming a Pauli code
merely from eventual clean return would be the first unsupported step.

The earlier [restored-address revision](../research/depth/RESTORED_ADDRESS_SESSION_RESULT.md)
weakens R98's sufficient premise from gatewise read-only access to exact
completed address return. [Paid selected packing](../research/depth/PAIRED_PERMUTATION_ROUTING.md)
handles ragged addresses without charging their larger union; joint
multi-output queries handle common banks. Internal query gates may modify
addresses, including selectors. The query must restore every physical
address/output/helper column with the specified XOR and literal phase.
This broadens the constructive route without improving R96 unconditionally.
A separate [exact six-slot counterexample](../research/depth/READ_ONLY_SIX_LAYER_REDUCTION.md#9-full-return-supplies-no-source-wide-inverse-coverage)
shows why full return alone supplies neither C3 nor inverse-Clifford
coverage for the specified source. It has identity target and is not a
general depth obstruction.

The preceding [native-defect round](../research/depth/NATIVE_DEFECT_SESSION_RESULT.md)
resolves the general **zero-clean read-only five-layer XOR degree** question:
[Section 5.6](../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#56-native-defects-a-general-five-layer-degree-bound)
proves degree at most **14** at full operator error below $`\sqrt2/24`$,
independently of dirty width. An invisible outer target shear gives degree
at most 11. No comparison-span, axis-permutation, Pauli-comparison or
Clifford-middle hypothesis is added. A native Pauli pullback makes one-band
matching applicable to the ordered defect; its bounded tableau and sign
degrees then detect the Boolean output. The sharp degree remains open.

The earlier subcases retain their larger error radius $`\sqrt2/12`$:
Clifford-valued middle gives degree seven (four for a visible outer shear),
Pauli-valued comparison gives degree two, and individual-axis permutation
gives degree 244 (41 for an invisible target shear). The span-preserving
112/21 bounds are superseded by the new general theorem at their radii.
The [shared-prefix invariant](../research/depth/READ_ONLY_FIVE_LAYER_REDUCTION.md#55-shared-prefix-comparisons-a-robust-quadratic-group-map)
still supplies class-two, exponent-four comparison structure and faithful
degree-four intrinsic features below $`\sqrt2/40`$. Target extraction no
longer requires putting the outer baseline into that group. A
[returned affine counterexample](../research/depth/BASELINE_WORD_OBSTRUCTION.md#5-a-returned-affine-query-can-have-an-odd-order-outer-baseline)
shows why such an ambient-group inference would fail.

The [scalar representation construction](../research/depth/SCALAR_REPRESENTATION_PRODUCT.md)
still realizes arbitrary Boolean phases in four abstract factors. Its
Clifford-translation source mechanism needs at least $`2^{d-1}`$ carrier
qubits and T gates per full source for a degree-d table, even after carrier
compression and at the stated small full-source error. This is not a lower
bound for fused whole-word synthesis. Fusion removes irrelevant carrier
motion but leaves the table-bit loading task; the retained implementation
uses R91's logarithmic-depth dirty indicators. No improved reader is proved.

The [preceding checkpoint](../research/depth/SIX_LAYER_TRANSPORT_SESSION_RESULT.md)
adds a [six-layer Pauli-comparison theorem](../research/depth/READ_ONLY_SIX_LAYER_REDUCTION.md):
degree at most 51 below $`\sqrt2/64`$, independently of dirty width.
It does not cover general Clifford-valued comparisons. An actual fully
returned six-slot identity circuit has odd-order comparisons and persistent
non-Pauli multiplicative defects, ruling out that proposed recurrence even
under full return. Merging its empty slots changes the specified comparison;
its fair circuit baseline is identity.

A separate [diagonal-source theorem](../research/depth/RETURNED_DIAGONAL_ECHO_DEGREE.md)
bounds Boolean scalar degree by $`\max(2,q-1)`$ for q arbitrary diagonal
work-source occurrences joined by read-only monomial Clifford intervals,
at error below $`2\sin(\pi/2^{\max(q,3)+1})`$. Every inverse occurrence
counts. The four-call echo radius improves to $`2\sin(\pi/32)`$.
This is a source-architecture bound, with a radius decreasing in q; it is
not a general native T-depth lower bound or a faster reader.

The [diagonal-source transport theorem](../research/depth/NATIVE_CLIFFORD_TRANSPORT.md)
gives six-layer degree at most 78, or 26 for visible outer target shear,
below $`\sqrt2/8`$, when the specified prefix is diagonal C3 in fixed
Clifford coordinates. The inner Clifford need not be monomial. C3
membership alone only rounds nearby endpoints and cannot bound tableau
degree; the charged known permutation-family stress test is not a
three-band returned query.

The [preceding checkpoint](../research/depth/TWO_OVERLAP_SESSION_RESULT.md)
extended the source certificate to overlap two. The
[overlap-two theorem](../research/depth/TWO_OVERLAP_SOURCE_STRUCTURE.md)
proves that an exact three-band C3 prefix with first/middle span overlap
at most two is diagonal C3 in fixed Clifford coordinates, even when
neither adjacent partial source is C3. It therefore certifies the existing
78/26 bounds. The [earlier source theorem](../research/depth/NATIVE_BAND_SOURCE_STRUCTURE.md)
still confines middle/outer noncommutation to at most 2t middle axes for
general overlap t. Full-space involution arguments resolve the paired
cores and zero-character cases; raw outer-axis spread can remain unbounded
through duplicate cancellation. C3 membership is an extra hypothesis,
and overlap is a property of the specified physical representation.

A separate [hybrid transport theorem](../research/depth/HYBRID_DIAGONAL_CORE_TRANSPORT.md)
handles either ordering of a wide diagonal C3 source and a fixed a-qubit
core, when the full source is C3. With
$`N(a)=3\,2^{5a}+15a+2`$ and $`B(a)=24N(a)+8a^2+12a`$, it gives
label degree $`4B(a)+2`$ and literal sign degree $`12B(a)+6`$ below
$`\sqrt2/8`$, independently of dirty width and source coefficients.
Zero- and one-qubit cores retain the sharper 26/78 bounds. A commuting
T-gadget list with r binary relations following a diagonal-C3 source has
a hybrid core of at most $`4(2^r-1)`$ qubits if their product is C3.
These explicit source hypotheses are not consequences of return.

The [preceding round](../research/depth/COUPLED_INCIDENCE_SESSION_RESULT.md)
extends the [native source theorem](../research/depth/ISOLATED_NATIVE_CORE_STRUCTURE.md)
to arbitrary overlapping nonzero-character incidence. If every
zero-character B axis commutes with C, exact full-space C3 forces each
interacting class to consist of one isolated two-axis B/C pair. The
full involution defect localizes each isolated A probe; native A
isotropy forces its B-neighborhood to be uniform on each activation
class. This closes coupled wider rows without bounding total overlap
and certifies the existing 78/26 transport bounds. The
[zero-character round](../research/depth/ZERO_CHARACTER_SESSION_RESULT.md)
adds two complementary certificates. The
[commutation-cover theorem](../research/depth/NATIVE_BAND_SOURCE_STRUCTURE.md)
closes every exact native three-band C3 source whose B axes each commute
with all of A or all of C, including arbitrarily many bad zero-character
axes. The [paired-orbit extension](../research/depth/ISOLATED_NATIVE_CORE_STRUCTURE.md)
proves that all zero-character axes are good if every singleton C orbit
is good. It uses character nonvanishing where the earlier diagonal Gram
bit can vanish. Both conclusions retain 78/26 transport at the same radius.

The separate [controlled-Pauli theorem](../research/depth/CONTROLLED_PAULI_TRANSPORT.md)
now assumes only an exact C3 source with scalar-Pauli computational
blocks in fixed input/output Clifford coordinates. A signed
normalization derives its canonical exact domain, removing the former
domain or fiber-compatibility premise. Label degree **58** and literal
sign/general target degree **174** hold below $`\sqrt2/8`$ independently
of control and target width. C3 derives quadratic labels and bounded
scalar phase degrees, retaining the weighted quartic term and all
Hermitian-section carries. The known non-semi-Clifford Gottesman–Mochon
source and its tensor powers satisfy the block hypothesis. Their paid
realization does not supply a native three-band source or a returned
query. The same proof home now gives a necessary-and-sufficient intrinsic
criterion: if E is the exact Pauli input domain and
$`N_p=S_{F\tau(p)F^\dagger}-I`$, scalar-Pauli coordinates exist exactly
when $`N_pN_q=0`$ for every $`p\in E^\perp`$ and every q.
This criterion includes literal signed normalization. A known C3
permutation fails it, so arbitrary C3 membership alone does not imply
the transport theorem. That source is not a native-three-band counterexample.

A separate [uniform transport theorem](../research/depth/BOUNDED_DEFICIENCY_TRANSPORT.md)
uses the intrinsic codimension k of the source's exact Pauli-preserving
domain, including degenerate domains and non-C3 sources. For k positive,
put $`D(k)=(2k+1)(4k+6)`$. It gives output degree at most $`3D(k)`$,
or $`D(k)`$ for visible shear, at error below
$`1/[4\sqrt2(2^k+1)]`$. Its radius is uniform in the source. Native
depth does not bound k; a wide single T band already has growing k.

The [commutative character bound](../research/depth/NONLINEAR_CHARACTER_TRANSPORT.md)
now covers every finite-dimensional commutative nilpotent binary algebra:
degree-a affine-input return needs at least $`2^{a-1}`$ coefficient
wires. A fully paid [noncommutative cycle construction](../research/depth/NONCOMMUTATIVE_CYCLE_CHARACTER.md)
escapes that bound, realizing degree-d AND on $`d^2`$ coefficient
wires. Its exact native schedule uses $`O(d^3)`$ T and Clifford gates,
depth at most $`64(d-1)`$, $`d^2+d`$ dirty qubits and no clean work.
The retained conjunction circuits are faster. Thus the width obstruction
is not representation-independent, but no improved general reader follows.

The cycle source has an exact reformulation as the complete word query
for every entry of a binary unitriangular inverse, with only Hadamard
basis changes on the lower-triangular coefficient register. Its
controlled-source wrapper prices both source calls and all helpers.
The new commuting-pivot identity implements the full in-place inverse
with $`8\binom d3`$ T gates, $`20\binom d3`$ elementary Cliffords,
T-depth at most $`8\max(d-2,0)`$, no clean work, and at most
$`\lfloor(d-1)^2/4\rfloor`$ additional arbitrary helpers for
$`d\ge2`$. All helpers return exactly. Overlapping gates prevent
inferring a parallel schedule from whole-pivot commutativity. The circuit
changes its matrix inputs internally and gives no improved native
read-only query or complete frame.

The [balanced construction](../research/depth/NONCOMMUTATIVE_CYCLE_CHARACTER.md)
now permits helper garbage during inverse computation and returns it
with the actual inverse after readout. For a d-by-d unitriangular input,
put $`q=\binom d2`$, $`h=\max\{\lceil\log_2d\rceil-1,0\}`$,
$`a=\lfloor d/2\rfloor`$, and $`c=\lceil d/2\rceil`$.
The exact complete inverse-word XOR query has

```math
\begin{aligned}
T&=64\binom d3,&G&\le160\binom d3+23q+4d,\\
D_T&\le128h,&W&=ac(c+3)/2.
\end{aligned}
```

Here W is additional arbitrary work beyond the q coefficient inputs
and $`d(d+1)/2`$ arbitrary output bits; d equal to one uses W zero.
All work and coefficient inputs return, no clean bit is used, and the
error is zero. Its actual reversed word and direct controlled variant
have separate complete ledgers in the proof home. The phase source has
the same T/depth/work orders. This improves the completed structured
source/query from linear to logarithmic depth at cubic rather than
quadratic work. It avoids the former repeated-child-inversion cost;
the determinant obstruction remains valid for its narrower returned
in-place, frozen-coefficient interface. The matrix inputs change
internally, which the strengthened R98 now permits at complete query
boundaries. Paid arbitrary-table coefficient access and extraction within
its count/work/depth budgets remain missing. The [representation audit](../research/depth/GAUSS_SUM_QUERY_ENCODING.md)
clarifies that inverse-entry encodings already exist, but generic
substitution has excessive count/work and affine dirty loading leaves
offsets. Its completed-oracle counting boundary is scoped to affine
wiring and the retained generic circuit costs. An exact finite-field
interpolation repair returns all work with native read-only access,
but repeated evaluations lose to the direct path-propagation baseline.
Sparse specialization, coupled-mask extraction and quantum interference
remain outside those boundaries. The complete cycle conjunction's
separate transport and fair baseline are unchanged.

The main constructive question is now a sublogarithmic arbitrary-table
query with completed address return and all R98 budgets. Native address
mixing is allowed; it need not first be converted to read-only access.
The separate lower-bound question remains a **native degree-versus-depth
recurrence**. For exact three-native-band C3 sources, the immediate
remaining structural issue contains both a bad singleton C axis and
a middle B axis interacting with both outer bands. Does the native
factorization force $`N_pN_q=0`$ for $`p\in E^\perp`$ and arbitrary q?
The intrinsic criterion makes this sufficient for 174/58 transport
without requiring diagonalization. A representation-dependent parameter is
not an unrestricted depth invariant. General non-C3 prefixes, robust
radii for growing depth, one-clean work and address mixing remain
additional gaps. None of these results changes the complete-frame upper law.
A robust growing **unrestricted native** query lower bound
would transfer to a legal deepest-layer real-frame subfamily with two
calls, the same clean/dirty allocation, at most doubled error, T-count and
T-depth, and elementary Clifford budget $`G(V)+G(V^\dagger)+2\le4G(V)+2`$.
The actual inverse is charged in the fixed H, S, CNOT alphabet.
A faster exact zero-clean multi-output query with complete address and
helper return would, by conditional R98, give

```math
D_T=O_L\!\left(N/b^2+d(2n+10)+\log_2\log_2(n+4)\right)
```

for both complete frame families at their original one-clean thresholds.
Here d is the hypothetical query's additive address-depth bound.
R91 supplies logarithmic d; no smaller d has been proved. These upper
and lower transfers have distinct query hypotheses. The
[Gaussian inverse repair](../research/depth/GAUSSIAN_INVERSE_EXTRACTION.md)
gives a paid identity and an affine-round obstruction, not an unrestricted
query or frame depth lower bound. The new
[scalar-loop proof](../research/depth/GAUSSIAN_SCALAR_LOOP_OBSTRUCTION.md)
also detects the central sign, closing bounded-round extraction by actual
affine quadratic phase primitives. General native work interference is
outside that architecture. The
[multiway cost audit](../research/depth/MULTIWAY_QUERY_COST_BARRIER.md)
rules out a termwise separated-reader shortcut at the square-root count
envelope; it does not lower-bound joint native synthesis.
R96's unconditional frontier is unchanged.

The separate constant-layer question is whether four T layers with one clean
qubit can approximate the explicit witness arbitrarily well when dirty
width, inner clean spread, middle crossing rank and both outer-axis
misalignment dimensions all grow. The
[bounded-spread result](#one-clean-inner-spread) rules out fixed spread;
[R97](#outer-axis-overlap) rules out bounded misalignment on either side,
even at unbounded outer-intersection dimension. All statements below
concern the complete initialized-isometry contract, not one state column.

The [preparation map](../research/depth/DEPTH_OBSTACLE_MAP.md) separates
that precision-dependent obstruction from the dimension-growing depth
question. It reconstructs the pre-R78 upper costs and lower proof gap;
the present spectator-padded witness itself has constant depth at fixed
accuracy by the existing one-qubit compiler.

<a id="precision-uniform-two-clean-count"></a>

<a id="precision-uniform-one-clean-count"></a>
<a id="optimal-banked-count"></a>

## Result C: the one-clean count frontier

For every $`n\ge1`$ and integer $`L\ge6`$, the
[one-clean theorem](../research/endpoint/GROWING_SOURCE_COMPILER.md#11-precision-uniform-one-clean-count)
gives error at most $`2^{-L}`$, with one clean qubit and universal
constants. Its base dirty reservation is $`B=L+n+7`$ for real frames
or $`B=L+n+8`$ for phase-dressed complex magnitude frames:

| Dirty allocation | T-count guarantee |
|---|---|
| Every $`b\ge B`$ | $`T=O(N+L)`$ |
| Exactly $`b=B`$, with $`L\ge\max\{6,n\}`$ | Worst-case $`T^\star=\Theta(N+L)`$ |
| Every $`b\ge2B`$ | Worst-case $`T^\star=\Theta(\sqrt{NL}+L+NL/b)`$ |

All upper constructions have $`G=O(NL)`$ and charge sources, queries,
actual inverses, and work return. The
[publication scope](../manuscript/PUBLICATION_SCOPE.md#result-c-the-one-clean-count-theorem)
is the complete statement; the [proof map](README.md#c-one-clean-compilation-and-its-corollaries)
gives its dependencies. At base width with $`L\lt n`$, or in the
transition $`B\lt b\lt2B`$, retain the upper bound without a general
matching claim. These count results add no depth guarantee.

## Constant-clean high-precision count

At the literal allocation

```math
a=2,\qquad \eta=2^{-N},\qquad L=N,\qquad b=N+n+7,\qquad n\ge3,
```

the selected worst-case count for each of the complete real and phase-dressed
complex magnitude families is

```math
2N-1\le T^\star\le O(N),\qquad T^\star=\Theta(N).
```

The [growing-source construction and literal phase composition](../research/endpoint/GROWING_SOURCE_COMPILER.md#9-complex-magnitude-endpoint)
give the upper bounds under the full return contract. The
[literal-phase root witness](FAULT_TOLERANT_COMPILER.md#103-an-explicit-width-independent-precision-witness)
gives the lower bound even at unrestricted width. This allocation saves
one dirty qubit for the complex family relative to C; the one-clean
complex theorem uses $`N+n+8`$ dirty qubits at $`L=N`$.

## Large-width T-depth

The [uniform subtree-depth theorem](../research/depth/UNIFORM_SUBTREE_DEPTH.md), R79,
applies to both families for every $`n\ge1`$, integer $`L\ge6`$,
and one clean flag.
Put $`B_D=L+n+7`$ for complete real frames and $`B_D=L+n+8`$
for phase-dressed complex magnitude frames $`D_\phi W_{\mathbb R}`$.
At every $`b\ge17B_D`$, use the following bound.

With $`\ell=\log_2(n+2)`$ and $`\lambda=\log_2(n+L+2)`$, put

```math
F(n,L)=\min\{nL,(L+\ell)\lambda\}.
```

One circuit has error at most $`2^{-L}`$ in the full initialized-isometry
norm and simultaneous bounds with absolute constants:

```math
T=O(\sqrt{NL}+L+NL/b),\qquad G=O(NL),\qquad
D_T=O(NL/b^2+F(n,L)).
```

Count is worst-case optimal at every eligible precision and width.
The inherited lower bounds also match depth on every nonempty interval

```math
17B_D\le b\le\sqrt{NL/F(n,L)}:
\qquad T^\star=\Theta(NL/b),\quad D_T^\star=\Theta(NL/b^2).
```

The one-clean high-precision branch gives the stronger $`D_T=O(NL/b^2+L)`$ whenever $`L\gt N/n^2`$,
with the same count, Clifford, and error bounds. R79 charges the native
subtree decoder's precision-dependent source calls and all-work predicates.
Its protected suffix allocation absorbs their table volumes with absolute
constants. The retained [one-clean baseline, R72](../research/depth/ONE_CLEAN_DEPTH.md),
supplies the terminal flag relocation, width and precision fallbacks, and
literal phase composition at the original dirty thresholds.

The complex extension uses the one-clean literal-diagonal schedule with
$`D_T=O(NL/b^2+L+\log_2(n+2))`$, then composes half-error factors
on the same flag and dirty pool. Its stated threshold includes that
precision increment; the high-precision bound applies to both families.
Every source preparation, query, actual inverse, and work-return error
is charged. Classical preprocessing terminates, with no efficiency claim.

At sufficient $`b=\Theta(\sqrt{NL})`$ and $`L=\Theta(\log n)`$,
R79 gives $`O(\log^2 n)`$ depth with error inverse polynomial in n.
More generally, $`L=o(n/\log(n+2))`$ gives sublinear depth.
The earlier two-clean
[compressed-column special case](../research/depth/SUBLINEAR_T_DEPTH.md#10-a-precision-uniform-extension-through-linear-precision)
retains its original source and recovery proofs; R79 supersedes its
resource bound.

The [subtree-projector decoder](../research/depth/SUBTREE_PROJECTOR_DECODER.md),
R78, is the proof home for native subtree reflections and coherent spectral
label erasure. Its fixed-accuracy theorem is now a corollary of R79.
For every fixed $`0\lt\eta\le1/64`$, put
$`L=\max\{6,\lceil\log_2(1/\eta)\rceil\}`$. At one clean
flag, every $`b\ge17(L+n+7)`$ for real frames, and every
$`b\ge17(L+n+8)`$ for phase-dressed complex magnitude frames, it gives

```math
T=O_\eta(\sqrt N+N/b),\qquad G=O_\eta(N),\qquad
D_T=O_\eta\!\left(N/b^2+\log_2^2(n+2)\right).
```

The same circuit has full initialized-isometry error at most eta.
At sufficient $`b=\Theta_\eta(\sqrt N)`$, this is optimal-order
T-count with $`O_\eta(\log^2(n+2))`$ depth. Both count
and depth match through
$`b\le\sqrt N/\log_2(n+2)`$, whenever the
interval above the family's literal lower threshold is nonempty: there
$`T^\star=\Theta_\eta(N/b)`$ and $`D_T^\star=\Theta_\eta(N/b^2)`$.
The [interval-dilation theorem R86](../research/depth/INTERVAL_DILATION_DEPTH.md)
first improved the fixed-accuracy general bound to an additive
$`O_L(\ell^{3/2})`$. The later
[product-coin theorem R87](../research/depth/COIN_DILATION_DEPTH.md)
removes its overlap-amplification factor. For every fixed integer
$`L\ge6`$, every n, and every original eligible real or complex dirty
allocation, one clean qubit gives

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell\log_2(\ell+2)),\qquad \ell=\log_2(n+2).
```

It covers arbitrary node-dependent magnitudes, signs, zeros, and the
prescribed leaf phases, with initialized-isometry error at most
$`2^{-L}`$. Retained global finite rounding makes a positive weighted
tree. A product-coin embedding followed by tensor Hadamards exposes
the fine marker depth through the last nonzero coin block. The paid
source and its actual inverse erase the tag and restore the work.
The [selected-routing lemma](../research/depth/SELECTED_DIRTY_ROUTING.md)
supplies its physical routing and variable-size query interfaces.
The subsequent [causal compact-coin theorem R89](../research/depth/CAUSAL_COMPACT_COIN_DEPTH.md)
retains those counts, widths, and full real/complex contract, and improves
the general depth to

```math
D_T=O_L\!\left(N/b^2+\ell(1+\log_2^*(n+2))\right).
```

Here $`\log_2^*x`$ counts repeated base-two logarithms until the value
is at most 2. Its compact fixed coins, causally batched local matchings,
and parallel reference transforms erase the coarse label in
$`O_L(\ell)`$ depth. Preserved-prefix recursion supplies the remaining
iterated-logarithmic factor. Every source, actual inverse, reference-work
zero test, and query is charged.

The later [dirty catalogue theorem R90](../research/depth/DIRTY_CATALOGUE_TERMINAL.md)
retains those counts, literal widths, and complete real/complex guarantees,
and removes the iterated-logarithmic factor:

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L(N/b^2+\ell).
```

Three R89 coarse reductions leave a small terminal. Its complete rounded
program selects local catalogue factors through an exact arbitrary-dirty
word query and XOR slot translations. A commutator using two disjoint
borrowed arrays cancels their unwanted actions exactly; the local factors
need no initialized ancilla. The selected-factor approximation, full
initialized-isometry error, and original-width splice are proved together.
At square-root dirty width, the general upper bound is now $`O_L(\ell)`$,
while the unrestricted lower bound remains $`\Omega(1)`$. Optimal depth
remains open. R79 retains its absolute all-precision scope; R90's constants
depend on fixed L. The finite catalogue implies no efficient classical
preprocessing, practical crossover, or scalable native emitter guarantee.

The [precision-terminal session](../research/depth/PRECISION_TERMINAL_SESSION_RESULT.md)
adds two absolute-constant frame refinements. Put

```math
F_{79}(n,L)=\min\{nL,(L+\ell)\log_2(n+L+2)\}.
```

[R92](../research/depth/UNIFORM_CATALOGUE_PRECISION.md) proves, for all
$`n\ge1,L\ge6`$ and an absolute $`C_*`$,

```math
D_T=O\!\left(NL/b^2+
\min\{F_{79}(n,L),\ell+2^{C_*(L+1)^3}\}\right).
```

In particular, logarithmic depth holds for the unbounded precision window
$`6\le L\le c[\log_2(\ell+2)]^{1/3}`$, with one absolute positive c.
The independent [grouped-catalogue theorem R93](../research/depth/GROUPED_CATALOGUE_PRECISION.md)
proves, on $`6\le L\le2\ell`$,

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)
                     +\frac{L\ell}{\log_2(\ell+2)}\right).
```

Both retain optimal $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$,
one clean qubit, the original family-specific literal widths, and the full
initialized-isometry contract. Select a complete circuit to take the best
available depth bound; R79's stronger high-precision clause is retained.
At $`L=\Theta(\ell)`$ and $`b=\Theta(\sqrt{NL})`$, R93 gives
$`O(\ell^2/\log_2(\ell+2))`$. The later
[portable-workspace theorem R94](../research/depth/PORTABLE_WORKSPACE_CATALOGUE.md)
proves the stronger uniform bound

```math
D_T=O\!\left(NL/b^2+\ell\log_2(\ell+2)\right),
\qquad 6\le L\le2\ell.
```

R94 has absolute constants and retains the same optimal T-count,
$`G=O(NL)`$, one external clean qubit, real width
$`b\ge17(L+n+7)`$, complex width $`b\ge17(L+n+8)`$, and full
initialized-isometry error $`2^{-L}`$. R79 compiles larger groups
directly to final accuracy. Carrying its entire local workspace through
the selected transfer makes the borrowed slot action cancel exactly,
without initialized residual labels. A constant number of groups uses
the [R91 parallel word query and paired router](../research/depth/PARALLEL_BILINEAR_WORD_QUERY.md).
R92 and R93 remain independent of R91. R94 reaches the
$`O(\ell\log_2(\ell+2))`$ square-root-width benchmark.

The [logarithmic-precision theorem R95](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md)
removes the remaining logarithmic factor:

```math
D_T=O(NL/b^2+L+\ell),\qquad 6\le L\le2\ell.
```

Its constants are absolute, and it preserves the optimal T-count,
Clifford budget, original real/complex literal widths, one initialized
qubit and full error contract. The
[coarse-flag construction](../research/depth/COARSE_FLAG_SPARSE_CORRECTION.md)
implements a coarse mathematical frame to final accuracy and adds a
constant-normalization sparse correction. The
[small-workspace source](../research/depth/SMALL_WORKSPACE_DESCENDANT_SOURCE.md)
supports a [local compiler](../research/depth/RELAXED_COUNT_LOCAL_DEPTH.md)
with depth linear in height plus precision and much larger dirty counts.
R94's complete-payload transfer absorbs those local costs globally.
Thus inverse-polynomial accuracy in n has $`O(\ell)`$ depth at
$`b=\Theta(\sqrt{NL})`$.

The subsequent [all-precision theorem R96](../research/depth/ALL_PRECISION_LINEAR_DEPTH.md)
removes that range restriction:

```math
D_T=O(NL/b^2+L+\ell),\qquad n\ge1,\quad L\ge6.
```

It preserves $`T=O(\sqrt{NL}+L+NL/b)`$, $`G=O(NL)`$, one external
initialized qubit, real width $`17(L+n+7)`$, complex width
$`17(L+n+8)`$, and full initialized-isometry error $`2^{-L}`$.
All constants are absolute. A globally priced borrowed-precision stage
and fast capped references give a constant number of height reductions.
A directly programmed high-precision terminal uses the original global
prefix in every fine mask and coarse instruction query. It does not
allocate a catalogue of all precision programs. At sufficient width
$`b=\Theta(\sqrt{NL}+L+n)`$, its depth is $`O(L+\ell)`$ for every
precision. No efficient preprocessing, practical crossover or scalable
native emitter is claimed. Optimal general depth remains open; the
selected manuscript Result D remains R79.

The [structured theorem R80](../research/depth/HOMOGENEOUS_MAGNITUDE_DEPTH.md)
attains logarithmic depth when angle magnitudes depend only on tree depth,
including arbitrary node signs and supplied leaf phases. At the same
fixed-accuracy reservations it gives

```math
T=O_\eta(\sqrt N+N/b),\qquad G=O_\eta(N),\qquad
D_T=O_\eta(N/b^2+\log_2(n+2)).
```

The signed balanced real subfamily has an exact version. This is a
restricted-family result, not a replacement for R79. The separate
[recursive terminal lemma R81](../research/depth/RECURSIVE_TERMINAL_DEPTH.md)
now covers every terminal height $`0\le M\le n`$, with depth
$`O_L(N/b^2+\ell\log_2(M+2))`$, where $`\ell=\log_2(n+2)`$.
A logarithmic protected suffix contributes only log-times-log-log depth;
the unrestricted early decoder still determined the square-logarithmic
bound at that checkpoint. The
[joint session record](../research/depth/JOINT_DEPTH_SESSION_RESULT.md)
identifies the remaining native mechanism and scoped failed routes.

The later [prefix-selected terminal theorem R88](../research/depth/PREFIX_SELECTED_COIN_TERMINAL.md)
improves that component for every $`0\le M\le n`$ to

```math
D_T=O_L\!\left(N/b^2+\ell\log_2\log_2(M+4)\right).
```

It compiles the deepest M layers of an arbitrary real frame, including
signed and singular tuples, while the unchanged external prefix may
select unrelated local frames. It retains one clean qubit,
$`b\ge17(L+n+7)`$, $`T=O_L(\sqrt N+N/b)`$, $`G=O_L(N)`$, and
full initialized-isometry error strictly below $`2^{-L}/4`$.
For $`M\le C_0\ell`$ with fixed $`C_0`$, the additive depth is
$`O_{L,C_0}(\ell\log_2\log_2(\ell+4))`$.
The [contiguous coin-layout router](../research/depth/SELECTED_DIRTY_ROUTING.md#7-a-logarithmic-depth-specialization-for-contiguous-intervals)
separately has $`O(\ell)`$ depth, including arbitrary dirty padding and
occupied inactive work. These component gains initially left R87's
complete bound unchanged. R89 subsequently replaces its coarse search
and gives the stronger arbitrary-height real-terminal bound
$`O_L(N/b^2+\ell(1+\log_2^*(M+2)))`$. R88 remains an independent
weaker terminal theorem and the proof home for sharp-query and fallback
interfaces used by R89.

The [bounded-density theorem R82](../research/depth/BOUNDED_DENSITY_DEPTH.md)
also improves the complete circuit when root probabilities have fixed
positive upper and lower ratios to a classically specified positive
product or affine-oriented reference. The latter uses the priced R83
reference core, including its leaf and marker corrections. With fixed
accuracy and ratio bounds, the same count and literal width guarantees
hold, with

```math
D_T=O\!\left(N/b^2+\log_2(n+2)(1+\log^*(n+2))\right).
```

The constants may depend on those fixed bounds. A native weighted
correction replaces repeated early threshold comparisons; every
recursive stage still pays for its global-prefix queries. This covers
nonhomogeneous magnitudes under the stated promise, including the
uniform-reference case. It does not require or establish an unrestricted
general improvement. R89 later attains this depth order for arbitrary
magnitudes, without the density promise.

The [affine-orientation theorem R83](../research/depth/AFFINE_ORIENTATION_DEPTH.md)
extends R80's logarithmic depth and count orders to depth-dependent
unordered magnitude pairs whose branch choices are affine functions of
the node prefix. Node signs, zeros, and supplied leaf phases are allowed.
Both R83 and the R82 extension retain, for fixed $`L\ge6`$, the literal
widths $`b\ge17(L+n+7)`$ for real frames and $`b\ge17(L+n+8)`$ for
the prescribed phase-dressed complex completion, with error at most
$`2^{-L}`$. Their constants may depend on L, and R82's also on the
fixed density bounds; these are not new all-precision claims.

The [deterministic-frame theorem R84](../research/depth/DETERMINISTIC_FRAME_DEPTH.md)
instead permits arbitrary node-dependent signed quarter turns and is
exact. For every $`b\ge17(n+13)`$, it has one clean qubit,
$`T=O(\sqrt N+N/b)`$, $`G=O(N)`$, and the same displayed
log-times-iterated-log depth order as R82, with absolute constants. Its
classical output-depth query erases the tag because every column is a
signed basis leaf. The
[branch-orientation session](../research/depth/BRANCH_ORIENTATION_SESSION_RESULT.md)
kept the general frontier separate: the factorization alone did not pay
for arbitrary nonlinear orientation at interior bias.

The [geometric-band theorem R85](../research/depth/GEOMETRIC_BAND_ORIENTATION_DEPTH.md)
now covers arbitrary nonlinear branch choices when the unordered
magnitude pair depends only on depth. It includes arbitrary node signs,
zeros, balanced pairs, and supplied leaf phases, with the same fixed-L
literal real and complex widths as R83. Its simultaneous bounds are

```math
T=O_L(\sqrt N+N/b),\qquad G=O_L(N),\qquad
D_T=O_L\bigl(N/b^2+\ell\log_2(\ell+2)\bigr).
```

The source's fine tag first selects a coarse geometric band. A
selected subtree unloader then exposes enough clean suffix to implement
the local orientation permutation through two returned-output XOR
queries. A product reference removes the fine label; only the coarse
band needs spectral search. The initialized-isometry proof includes
every column, literal phases, actual inverses, leakage, dirty inputs,
and references. This reaches sub-square-logarithmic depth for the
arbitrary-orientation promise, without a general whole-register
permutation compiler. R87 subsequently matched its asymptotic depth for
arbitrary magnitudes, and R89 improves the general bound further.
R90 reaches logarithmic depth for arbitrary magnitudes, including this
promised family; R85's independent proof remains valid.

The [bounded in-place conversion session](../research/depth/PERMUTATION_BRIDGE_SESSION_RESULT.md)
stopped at its mechanism checkpoint. Its exact word and encoding
restrictions apply to specified architectures, not unrestricted native
circuits. They supply no new general upper or lower bound. The sufficient
whole-register reduction and the remaining recovery/return obligation are
recorded there. R86 subsequently obtains a general improvement through
interval dilation rather than that standalone permutation interface.

The [joint orientation session](../research/depth/JOINT_ORIENTATION_RECOVERY_SESSION_RESULT.md)
records R85's discovery and review. The later
[general-magnitude session](../research/depth/GENERAL_MAGNITUDE_RECOVERY_SESSION_RESULT.md)
records R86 and its reserve routes. The
[factorization session](../research/depth/INTERVAL_RECOVERY_FACTORIZATION_SESSION_RESULT.md)
then proves R87 through product-coin spectral readout. The
[continuation](../research/depth/COIN_CONTINUATION_SESSION_RESULT.md)
proves R89 through compact coins, causal batching, and prefix recursion.
The [whole-terminal session](../research/depth/TERMINAL_SOURCE_GAUGE_SESSION_RESULT.md)
proves R90 through fixed-number coarse reduction and dirty catalogue
commutators. The [precision follow-up](../research/depth/PRECISION_TERMINAL_SESSION_RESULT.md)
then proves R91–R93. The [catalogue bootstrap](../research/depth/CATALOGUE_BOOTSTRAP_SESSION_RESULT.md)
proves R94 and reaches its growing-precision benchmark. The
[R95 construction](../research/depth/LOGARITHMIC_PRECISION_DEPTH.md)
then reaches logarithmic depth at inverse-polynomial accuracy through
coarse flags, sparse correction and local compilation with relaxed counts.
The [precision bridge](../research/depth/BORROWED_PRECISION_BRIDGE_SESSION_RESULT.md)
then proves R96's all-precision upper law through globally priced sources,
fast references and a directly programmed terminal. The
[continuing research mandate](../WORKSPACE.md#continuing-toward-the-general-answer)
retains optimal general depth and dimension-growing lower bounds as open
questions. All further
improvements must retain the full native return and resource contract.

## Constant-depth frontier

### Zero-clean error floors and fixed-clean thresholds

The [constant-depth obstruction](../research/depth/CLEAN_AWARE_DEPTH_OBSTRUCTION.md),
R73, strengthens that constant. The complete frame
$`W_n=I_{2^{n-1}}\otimes R_y(\theta)`$ with
$`\theta=\tfrac12\arccos(1/3)`$ has exact zero-clean error floor

```math
c_0=\sqrt{2-\frac{2+\sqrt2}{\sqrt3}}\approx0.169714
```

at T-depth at most three, uniformly over n and dirty width b. A Clifford
attains this floor, so integer $`L\ge3`$ requires at least four layers.
For every fixed positive clean budget a, a positive threshold
$`\varepsilon_a`$ exists; R74's elementary refinement below makes it
explicit. These statements
allow arbitrary Clifford interlayers and dirty/reference inputs.
Separately, an explicit one-clean witness requires at least three layers
for integer $`L\ge8`$. Zero leaf phases put both witnesses in the complex
magnitude family. These depth lower bounds are constant.

The [four-layer operator-rigidity theorem](../research/depth/FOUR_LAYER_OPERATOR_RIGIDITY.md),
R75, strengthens the zero-clean obstruction without restricting crossing
rank. Distinct Pauli images under two-layer circuits have operator-norm
distance at least

```math
\rho=\frac{\sin(\pi/8)}{\sqrt2}\approx0.270598.
```

The same complete frame therefore has error at least
$`\rho/2\approx0.135299`$ at T-depth at most four, with zero clean qubits
and arbitrary dirty width. Since $`1/8\lt\rho/2`$, every integer
$`L\ge3`$ requires at least five layers. This lower bound does not claim
that the four-layer error floor is sharp. With one initialized clean
qubit, error below $`\rho/2`$ forces exact terminal clean-sector
preservation.

### One-clean midpoint restrictions

The [one-clean midpoint theorem](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md),
R76, gives two further error floors with arbitrary dirty width and
crossing rank. For a split $`V=V_2V_1`$ with each half of T-depth at
most two, let $`S=V_1Z_CV_1^\dagger`$ be the encoded clean observable.

| A four-layer realization admits | One-clean consequence for the same witness |
|---|---|
| A split with S Pauli | Error at least $`\tau=\rho/2`$; integer $`L\ge3`$ requires five layers in this subclass |
| A split with S Hermitian Clifford | Error at least an explicit positive $`\eta_C`$, independent of width and crossing rank; the proof gives a very conservative constant |

Terminal return does not imply either midpoint hypothesis. R76's
[remaining question](../research/depth/ONE_CLEAN_MIDPOINT_RIGIDITY.md#6-the-remaining-one-clean-question)
is freezing the compressed conjugations for a non-Clifford midpoint;
commutant lifting and phase discreteness already extend to that case.
A separate spectral argument gives error at least $`1/196608`$ when
the exactly returned occupied block has finite projective order; in
particular, four-layer full circuits of finite projective order are
excluded at $`L\ge18`$. Generic Clifford+T circuits need not satisfy
that condition.

### One-clean inner spread

The [inner clean-spread theorem](../research/depth/ONE_CLEAN_INNER_SPREAD.md),
R77, extends the compressed rigidity argument to non-Clifford midpoints.
In either two-layer encoder, let t count the first layer's axes that
anticommute with the incoming clean Pauli after its initial Clifford.
For every fixed t there is an explicit positive error floor, independent
of dirty width and middle crossing rank. The proof expands the midpoint
into $`2^t`$ Clifford terms and applies R74's finite-boundary argument
without crossing phases. Define $`t_{\min}(V)`$ by minimizing over
both encoders and every native four-layer representation. At error
$`2^{-L}`$, its explicit inequality implies

```math
t_{\min}(V)=\Omega(\log\log L).
```

Thus the first and last T layers must each have diverging clean-active
size. Together with R74 below, all four layers must grow in any
hypothetical arbitrarily accurate one-clean sequence. This does not
exclude unbounded spread. An exact two-qubit example also shows that
the separate encoder flips cannot always be replaced by one common
sign-compatible flip, even allowing an arbitrary unitary flip. Its
occupied defect is a positive constant, so it does not refute the
remaining gap conjecture.

### Outer-axis overlap

The [outer-overlap theorem](../research/depth/OUTER_OVERLAP_RIGIDITY.md),
R97, gives a positive initialized-isometry error floor depending only on
$`r=\min\{\mu_A,\mu_B\}`$. Here $`I=L_A\cap L_B`$ is the
intersection of the two native outer axis spaces, and $`\mu_A`$ is
$`\dim I`$ minus the number of individual A axes lying in I; define
$`\mu_B`$ similarly. The physical lists are independent, so these are
nonnegative dimensions. Both inner clean spreads, crossing rank, dirty
width and the full intersection dimension may grow while this floor
remains fixed.

At $`r=0`$, the intersection is spanned by a subset of either physical
list and the floor is R76's uniform Clifford-midpoint threshold. Literal
matching between the lists is unnecessary. An arbitrarily accurate
four-layer sequence must have both misalignment dimensions diverge in
every admissible native representation. That unbounded case remains
open. These are constant-depth obstructions, not matching
dimension-growing lower bounds for R96.

The [representation limits](../research/depth/FOUR_LAYER_REPRESENTATION_LIMITS.md)
retain two distinct conclusions. Exact clean return does not guarantee an
exact same-space phase-polynomial descent, even with one additional dirty
wire. Growing center mixing in the displayed R77 tensor family also does
not obstruct an alternative native split: that family has a
Pauli-preserving three-layer realization. Neither example refutes the
remaining occupied-sector gap.

### Middle crossing rank and the remaining gap

The [four-layer crossing-rank theorem](../research/depth/FOUR_LAYER_CROSSING_RANK.md),
R74, gives a positive error gap for the same $`1/3`$ witness for each
fixed clean budget a and rank bound R, for circuits admitting a four-layer
decomposition whose middle active symplectic crossing rank is at most R.
Its elementary initialized-observable argument strengthens this to an
explicit rigidity gap $`2^{-O_a(r^2)}`$. Consequently, at error
$`2^{-L}\lt c_0`$, every such approximation satisfies

```math
r_{\min}=\Omega_a(\sqrt L),
```

where the minimum is over all four-layer decompositions of the circuit.
The proof supplies explicit constants for each fixed clean budget and
also an explicit fixed-clean three-layer threshold. Its original
polynomial-limit closure proof remains valid. With positive clean
budget, this rank restriction permits four-layer circuits whose crossing
rank grows with precision; it is not an unrestricted five-layer lower
bound. The width-independent four-layer Hilbert--Schmidt gap remains
unresolved, separately from R75's operator-norm result.
Optimal depth remains an open problem.

The [entry-indexed encoding bound](../research/depth/ENTRY_INDEXED_ENCODING_LIMIT.md)
shows that the entry-indexed overlap architecture needs worst-case normalization
$`\Omega(\sqrt n)`$, even at fixed operator error. That is a bound
on this encoding method, with no unrestricted native-depth consequence.
R78 uses native spectral-label erasure instead of that overlap recovery.
T-depth permits arbitrary Clifford interlayers and is distinct from total
physical depth. The [depth results](../research/README.md#large-workspace-t-depth)
record the new complete construction and preserve the earlier component
and interface results under their original hypotheses.

## Earlier endpoint constructions

The [research index](../research/README.md#constant-clean-high-precision-count)
is the catalogue of earlier constructions and their exact hypotheses.
Their method-specific open steps are alternatives to the resolved endpoint,
not gaps in the growing-source theorem. Unique proofs, negative results,
and executable evidence remain at their original paths.

<a id="proved-linear-count-families"></a>
<a id="bounded-diagonal-factorization"></a>
<a id="a-regular-fixed-tree-cayley-family"></a>
<a id="source-program-resynthesis"></a>

| Retained route | Catalogue and boundary |
|---|---|
| Linear-count promised families | [Native primitives and promised families](../research/README.md#native-primitives-and-promised-families): distinct promises and workspace reservations; several support C |
| Fixed-tree, Cayley, scattering, residual, and ordered-bank reductions | [Full-frame reductions](../research/README.md#full-frame-reductions-and-remaining-native-steps): compact representation and conditioning alone do not price native recovery |
| Source-program resynthesis | [Source programs and encodings](../research/README.md#source-programs-and-encodings): paid fusion rules and scoped history/query restrictions |
| Bounded diagonal and reflection factorizations | [Diagonal and reflection architectures](../research/README.md#diagonal-and-reflection-architectures): certified selection and restrictions retain their stated mixers and regularity assumptions |

## Full-circuit requirement

Every compiler must charge source construction, tables, basis changes,
queries, actual inverses, amplification, and workspace return on the same
circuit. Result C also charges
$`G=O(NL)`$ elementary Clifford gates. Certified classical preprocessing
has no asserted uniform efficiency bound; the depth question retains its
displayed Clifford budget. Neither permits QRAM, a supplied precision state,
measurement, reset, postselection, or uncharged coherent evaluator.
A lower bound must apply to this full model; a failure of one source,
query interface, or factorization is a restriction of that method.
