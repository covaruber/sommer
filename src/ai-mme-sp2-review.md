# ai_mme_sp2: Three-Pass Mathematical Performance Review

## Implementation Status: 2026-10-12

The implementation below is a completed, validated batch, not completion of every research proposal. The four review passes later in this file describe the pre-implementation baseline; their source line anchors may shift.

Implemented in [MNR.cpp](MNR.cpp) and [mme_paths.h](mme_paths.h):

- Native AR1 precision derivatives dispatch before the generic dense multiplication.
- Stationary AR(p), including heterogeneous scaling, and ANTE precision assemble innovation-row outer products instead of dense matrix chains. Factor storage/evaluation remains dense; this is not yet a precision-first sparse descriptor API.
- An overflow-checked ordered tensor component planner supports arbitrary-length shaping-factor lists. Every structurally diagonal axis contributes component labels, irrespective of factor position. CHOLMOD Schur candidates and generic inverse-cache groups use these components; the assembled C graph still merges/borders coupled groups. Relationship-graph and arbitrary custom block-component splitting remain future extensions.
- Ordinary Schur border/merge density feasibility is checked before cost selection, so a feasible alternative can survive rejection of the nominal cheapest plan. Existing thresholds and latent/chain priority are preserved.
- LDLT selected inversion caches topology dependency offsets within a bounded share of the dense-memory budget; uncached columns retain exact lookup recursion. Each C/D/R subset can retain at most 1,048,576 address entries, additionally capped to 1/48 of the configured byte budget. Reused caches release address tables when that budget shrinks.
- Final LDLT mode-1 extraction reuses accepted-state selected values without a deep subset copy or another Takahashi pass. CHOLMOD's existing final LDLT fallback and public subset contract are unchanged.
- Generic CHOLMOD inverse preparation accounts for accumulated fixed/random cached blocks and conservative overlapping workspace estimates. Unsupported or over-budget entries retain exact solve fallbacks. Structured engines retain their own resource rules; this is not a process-memory cap or a complete demand-driven strategy selector.
- Grid factor covariance/precision/logdet reuse uses exact factor/section/ordered-level keys, local to one coefficient evaluation, and references earlier block storage rather than storing a second dense cache.
- Ordinary and latent Schur updates use T'T, T=L^-1 B, retaining the second triangular solve required for subsequent D^-1 B consumers.
- Nyström preconditioning retains the weighted basis for the evaluation. Matrix-free PCG accepts an already-converged finite true-residual warm start before positivity checks.
- Non-diagonal PCG residual traces contract projected probes through existing complete/ragged/weighted derivative applications instead of assembling Btrace. Probe signs, solve tolerance, and the estimator are unchanged.

Regression additions reuse the existing native policy, dense-Schur, MME-path, spatial-residual, and REML/ML test files. They cover every diagonal-axis combination through five factors, end-to-end three/four/five-factor models, factor positions, REML/ML, all direct inverse modes, correlated ragged residual merging, cache budgets, and same-probe residual contractions. New diagnostics report selected dependency entries, grid factorization/reuse counts, and non-diagonal PCG residual contractions.

### Measured Results

Compared an isolated install of untouched HEAD C++ (with the existing user wrapper edits preserved) against the implementation, in separate R processes. Benchmarks use four shaping factors (diagonal x AR1 x diagonal x US), dense relationship precision, three iterations, and identical thread environment settings. One warmup was excluded; BLAS runtime thread control is unavailable on this host, so environmental settings are not a measured thread-count guarantee.

| Fixture | Repeats after warmup | Baseline median | New median | Ratio | Schur groups |
| --- | --- | --- | --- | --- | --- |
| 24 IDs | 5 | 0.163 s | 0.149 s | 1.09x | 1 -> 4 |
| 120 IDs | 3 | 1.201 s | 0.410 s | 2.93x | 1 -> 4 |

Likelihood, covariance estimates, BLUE/BLUP, score, information, and covariance-parameter uncertainty match at tolerance 1e-7. These are fixture-specific gains. An earlier one-shot 1.007 s -> 0.194 s comparison was startup-sensitive and is NOT the reported speedup. A saved pre-change same-probe PCG ragged-residual fit also matches after the contraction change.

The isolated implementation library is `/tmp/sommer-review-implementation-lib`; the normal user library was not replaced. Fresh preclean installs and focused tests passed after every numerical slice. Full-suite validation found no failures outside [test-predict-reported-issues.R](../tests/testthat/test-predict-reported-issues.R): its 18 failed expectations and four erroring tests reproduce in the untouched HEAD numerical baseline. Six tests skip for unavailable BLAS control or CRAN-only gates. The unrelated baseline failures were not modified.

### Completed After Restart

The work resumed after a restart. The following stages are now implemented in the current source and have focused regressions:

- D-first ML for eligible direct fits factors D and a bounded fixed-effect Schur block, routes joint solves/sensitivities through the block equations, and retains the original C path for mode 1, retained joint-factor caches, excessive fixed borders, and failed admission. Diagnostics report activation/evaluation counts.
- Separate CHOLMOD selected-inverse subsets are built for generic ML D and irregular residual R traces. The final generic C mode-1 export can use the current supernodal factor; Schur-engine, stale/non-supernodal, and over-budget paths retain Eigen LDLT extraction. The output remains an incomplete selected subset (`CiComputed=FALSE`), and PEV/prediction consumers are tested.
- Matrix-free PCG batches panels through the tensor operator and independent per-column CG recurrences, with true-residual checks and bounded panels.
- Tensor graph grouping uses connected components from structurally diagonal factors, fixed covariance factors, and relationship precision. Free generic factors are conservative. Assembled C still validates/merges candidate components.
- Native/fixed factor evaluations are cached by exact parameters with a bounded entry budget. Positive-definiteness checks use sufficient diagonal/CSM certificates and fall back to eigensolves if inconclusive. Support-aware Kronecker contraction has both cached mixed-radix coordinates and a bound on support-product materialization.
- SAR/CAR precision and log determinants use their normalized underlying system matrices: SAR `P=c B'B`, CAR `P=c Q`. Dense covariance evaluation remains in the R callback, so this removes the second generic covariance inversion rather than making the whole family sparse. Formula checks and LDLT/CHOLMOD output parity tests pass.
- SAR/CAR working-coordinate precision and log-determinant derivatives now use the normalized system formulas directly. Covariance derivatives follow from `dK=-K(dP)K`; central differences and pre-change fit parity pass.
- `computeCi=2` solves identity columns in panels of 32, keeping only the dense output plus one RHS/solution panel instead of simultaneous full identity and full solution matrices.
- CHOLMOD mode-1 selected inverse is budgeted and guarded by factor type/size. Its returned values, PEV, prediction errors, and over-budget Eigen fallback are tested against mode 2.
- PCG reports empirical SLQ probe standard error and the maximum verified true relative solve residual. These are diagnostics only; the likelihood estimator, probe schedule, and tolerance policy remain unchanged.
- An exact block-banded REML/CHOLMOD path supports assembled coefficient graphs up to block bandwidth four. It is gated by `options(sommer.mme.blockBanded=TRUE)` and off by default because the tested AR2 fixtures are slower than generic CHOLMOD. It uses exact band Cholesky/solves and exact solve fallbacks for traces.
- The latent Schur engine now admits one validated low-rank-plus-diagonal factor at any axis position in a longer Kronecker product. Its conditional base precision/log determinant are rebuilt every numeric factorization; first/middle/last factor parity is tested. Multiple low-rank axes remain unsupported.

The one-off SAR/CAR observation-space comparison was invalid for that fixture; it returned zero parameters and no information matrix. The permanent tests instead check exact covariance/precision/logdet identities and direct-solver parity. A fixed-interior baseline/current comparison also matched likelihood, theta, BLUP, AI, and score at 1e-7.

### Reevaluated Remaining Plan

Remaining work is staged rather than treated as one unchecked checklist:

1. **Exact low-risk kernels:** direct working-coordinate SAR/CAR precision derivatives; analytic AR(p)/ANTE derivatives; performance counters for current cache/trace paths. Keep exact generic fallbacks and compare against fixed-parameter references.
2. **Accuracy-informed PCG:** use the new probe standard error and residual diagnostics to study symmetric diagonal-congruence SLQ and matched-accuracy trial acceptance. Do not change defaults until finite-probe bias/variance and convergence behavior are measured against direct solvers.
3. **Direct-solver planning:** benchmark complete factor/trace/solve cost by ordering; calibrate total-memory and inverse-demand admission; evaluate stable structural assembly envelopes. Keep requested backend semantics and current fallback until a plan proves better end-to-end.
4. **New solver architectures:** block-banded elimination exists as an exact opt-in prototype but is slower on measured fixtures. Product-grid elimination and augmentation with multiple low-rank axes still need separate derivations for tensor ordering, missing data, residual/weight coupling, derivative support, and uncertainty output. Do not promote the prototype or substitute new engines for CHOLMOD/LDLT without a measured crossover.
5. **API/model reparameterizations:** covariance-information opt-out and common-scale profiling change output/uncertainty or fit semantics. Define the contract and constraints before implementation. Spectral rotation and model-changing approximations remain excluded.

Not yet implemented: truly sparse/banded factor storage and derivatives (AR/ANTE/SAR/CAR currently retain dense factor precision); direct SAR/CAR evaluation without dense covariance; Toeplitz/MA/metric structured operators; stable structural assembly envelopes; automatic banded dispatch (the prototype is opt-in and slower); product-grid elimination and multiple-low-rank-axis augmentation; calibrated demand-driven trace planning/touched-row CHOLMOD workspaces; symmetric/adaptive SLQ; covariance-information opt-out; and common-scale profiling with constrained-parameter/uncertainty reconstruction.

Focused validation after the restart: D-first/ML, panel PCG and its diagnostic fields, dense Schur, weighted selected R, variance families (including SAR/CAR identities), spatial residual cache/certificates, CHOLMOD selected/fallback inverse and PEV, and chunked full-inverse tests passed.

Latest full suite: 164 test cases; six skips (three BLAS-controller unavailable, three CRAN-only). The 18 failed expectations and four errors are confined to `test-predict-reported-issues.R` and match the separately rebuilt untouched numerical baseline with the user's current `R/mmes.R`; no other test file fails. Build and suite ran against `/tmp/sommer-review-implementation-lib`. `git diff --check` is clean.

### GitHub Comparison

- Compared pristine `origin/unified` at `013d5b68f934196179da27b99dd978ac0a62ef26` (same as local HEAD) against the local source build with all uncommitted changes, including the existing `R/mmes.R`/`man/mmes.Rd` modifications. The GitHub tree was archived into a separate temporary source and installed independently; no worktree changes were omitted from the local build.
- Same deterministic data, explicit `solver="cholmod"` except spatial derivative fixtures (`solver="ldlt"`), warmup excluded, four timed repetitions for the tensor/AR2 workloads and three for SAR/CAR. Timing variability was noticeable; BLAS runtime thread control is unavailable, so thread environment settings are not a measured thread count.

| Workload | GitHub median | Local median | Ratio (GitHub/local) | Numerical comparison |
| --- | --- | --- | --- | --- |
| Four-factor dense-relationship REML, 120 IDs, 3 iterations | 1.148 s | 0.707 s | 1.62x | Likelihood, theta, BLUP, score, AI, theta SE match to 1e-6 |
| Same four-factor model, ML, 3 iterations | 2.601 s | 0.636 s | 4.09x | Same outputs match to 1e-6 |
| AR2 random effect, 120 levels, `computeCi=1`, banded option off | 0.218 s | 0.158 s | 1.38x | Fits/PEV match; common `Ci` entries and diagonals differ by at most 5.3e-16. Selected patterns differ, as expected for incomplete mode 1. |
| Free SAR random factor, 20 sites x 40 IDs, 1 iteration | 2.036 s | 0.282 s | 7.22x | Likelihood, theta, BLUP, score, AI match to 1e-6 |
| Free CAR random factor, same fixture | 1.976 s | 0.183 s | 10.80x | Likelihood, theta, BLUP, score, AI match to 1e-6 |

Ratios above 1 favor local. These are fixture-specific results, not a broad performance guarantee. Mode-1 pattern identity is not a valid parity requirement because that output is explicitly an incomplete selected subset. The benchmark RDS files and both temporary installs are under `/tmp/sommer-github-bench-*.rds`, `/tmp/sommer-review-github-lib`, and `/tmp/sommer-review-implementation-lib`.

### Block-Banded Prototype Gate

An exact REML/CHOLMOD block-banded path is implemented behind `options(sommer.mme.blockBanded=TRUE)` and is off by default. It passes AR2-vs-LDLT parity for `computeCi=0/1/2`, but measured fixtures are slower than generic CHOLMOD:

| Fixture | Generic CHOLMOD | Opt-in banded | Generic/banded |
| --- | ---: | ---: | ---: |
| AR2, 600 scalar blocks, `computeCi=0` | 1.267 s | 1.361 s | 0.93x |
| AR2, 600 scalar blocks, `computeCi=2` | 1.250 s | 1.412 s | 0.89x |
| AR2, 16 blocks x 120 IDs, `computeCi=0` | 0.898 s | 1.266 s | 0.71x |
| AR2, 16 blocks x 120 IDs, `computeCi=2` | 1.338 s | 1.492 s | 0.90x |

It remains a correctness-tested prototype, not a performance improvement. Keep automatic dispatch off until factor and trace kernels show a repeatable crossover.

### Original Review Baseline

The four review passes below describe the pre-implementation baseline; their source line anchors may shift as code changes. This status section supersedes historical proposals that are marked pending in those passes.

## Baseline And Scope

- Reviewed on 2026-10-12: HEAD `013d5b68f934196179da27b99dd978ac0a62ef26`, using the source in this folder, not an installed package or historical notes.
- Existing local modifications: `R/mmes.R` and `man/mmes.Rd`. Neither was changed by this review.
- Scope: `ai_mme_sp2` in [MNR.cpp](MNR.cpp#L3549), its numerical helpers, covariance descriptors, and nearby tests. "Families/types" is interpreted as covariance families, not GLMM response families.
- This is a source-based review. Proposed speedups are hypotheses, not measured gains. Exact algebra and approximate statistical calculations are distinguished below.
- Three sequential passes: (1) solver algebra and repeated work; (2) covariance representations and admission rules; (3) adversarial checks, validation, and implementation priorities. Each pass is saved before the next.

## Run 1: Solver Algebra And Repeated Work

### R1.1: Batch The Matrix-Free Operator And PCG Solves (High Priority)

Evidence: [matrix-free solve](MNR.cpp#L8525), [multi-RHS dispatch](MNR.cpp#L8734), [probe preparation](MNR.cpp#L8954), [prior operator](MNR.cpp#L10417). Assembled PCG already batches independent CG recurrences; matrix-free PCG still loops over RHS columns. Each vector operator call gathers coefficients, allocates matrices, applies `Ai`, and scatters results.

Proposal: add an operator accepting a matrix of RHS vectors. Pack covariance coordinates and RHS into matrix panels; apply relationship precision with sparse-matrix/dense-matrix products, then tensor-mode products. Use independent per-column CG coefficients and convergence masks, as in the existing assembled implementation. This is batched independent CG, not a new block-Krylov algorithm.

Benefit hypothesis: amortize sparse traversal and allocations across probes and sensitivities without assembling the random prior. Chunk panels to bound memory; avoid carrying many converged columns through expensive operator products.

Cheap falsification: compare batched and scalar operator products on deterministic panels, then compare solves by true residual and fitted outputs using identical probes. Benchmark dense and sparse `Ai`; reject if panel memory or extra inactive-column work outweighs savings.

### R1.2: Use Preconditioning In SLQ, Not Only In CG (High Potential, Higher Risk)

Evidence: [Nyström setup/application](MNR.cpp#L8425), [SLQ](MNR.cpp#L8844), [logdet call](MNR.cpp#L10817). PCG has diagonal and Nyström preconditioners, but SLQ applies the original C. Better CG conditioning therefore does not automatically reduce Lanczos work.

Exact identity for SPD M = B B': `log|C| = log|M| + log|B^-1 C B^-T|`. For diagonal M, the symmetric operator is `M^-1/2 C M^-1/2`; its determinant correction is the sum of log diagonal entries. Start with this diagonal case. Applying `M^-1 C` to ordinary symmetric Lanczos is NOT valid: it is generally nonsymmetric.

The current Nyström preconditioner is diagonal plus low rank. Its logdet follows the determinant lemma, but its inverse application alone does not supply a symmetric square-root factor application. Defer that extension until a correct factor/operator is available. Hutchinson traces must either remain in original coordinates or receive a correctly derived congruence transformation. The original-coordinate Lanczos inverse warm start also needs rederivation.

Cheap falsification: on a small ill-scaled SPD matrix, compare corrected SLQ against exact logdet at full Krylov dimension; then measure error versus operator count at practical depths. Use paired deterministic probes and retain the existing true-residual checks. Do not claim bit-identical finite-probe estimates after changing the probing basis.

### R1.3: Avoid Two Large Factorizations For ML (High Potential, Architectural)

Evidence: [C and D factorizations](MNR.cpp#L10633). Both LDLT and CHOLMOD factor C for joint solves and D for ML determinant/score corrections.

For C = [A B; B' D], factor D and the fixed-effect Schur complement S = A - B D^-1 B'. Then `log|C| = log|D| + log|S|`, fixed solution `b = S^-1(rhs_x - B D^-1 rhs_u)`, and random solution `u = D^-1(rhs_u - B' b)`. The same block identities supply joint sensitivities and final uncertainty while ML traces continue to use D^-1.

Admission: small, full-rank fixed-effect block; budget `Nu * nX` workspace and dense S; otherwise retain the existing joint-C backend. Large nX or very sparse full C can make Schur elimination slower. For REML this is an alternative strategy, not an unconditional win.

Cheap falsification: small dense reference checks for joint solves, both logdets, D-only scores, and C-based uncertainty; benchmark end-to-end ML against the existing two-factor path at several nX/Nu ratios. Preserve rank-deficiency handling.

### R1.4: Reduce Final CHOLMOD computeCi=1 Refactorization (Pending Run 3 Verification)

The current architecture documents a final Eigen LDLT subset extraction for CHOLMOD fits. The existing supernodal selected-inverse helper may allow reuse of the final CHOLMOD factor rather than a second numeric factorization. This needs an output-contract audit: a Takahashi subset is not a full sparse inverse, and structural fill/order can differ between backends. A Schur engine may have no marginal CHOLMOD factor to reuse.

Cheap falsification: compare public `Ci` dimensions, symmetry, represented entries, permutations, and prediction/PEV consumers across final-inverse modes. Keep the existing fallback when no compatible factor/subset exists.

### R1.5: Specialize Accepted-State Derivatives Before Dense Products (All Solvers)

Evidence: [precision derivative builder](MNR.cpp#L5398). It first forms `-P dK P` for every factor parameter; the native homogeneous AR1 branch subsequently overwrites that result with a tridiagonal analytic derivative. The analytic branch therefore preserves sparsity but currently does not avoid the initial generic multiplication.

Small first implementation: dispatch the existing AR1 analytic derivative BEFORE the generic multiplication. Hypothesis: avoid unnecessary cubic dense work, with identical working-coordinate derivative values. Generalizations to other families are reviewed in Run 2.

Cheap falsification: compare direct derivative entries and logdet derivatives against central differences, including rho=0 and negative rho; then run block-chain fit equivalence. Do not change coordinate transforms or PD safeguards.

### Already Present: Do Not Reimplement

- Factor-wise Kronecker inverse/logdet and tensor applications.
- LDLT C/D symbolic caching and selected-inverse topology reuse; CHOLMOD symbolic caching.
- Accepted-state C selected inverse/probe preparation after likelihood trial acceptance.
- SLQ-derived inverse warm starts, deterministic common probes, full Lanczos reorthogonalization, and assembled multi-RHS PCG.
- Nyström preconditioning, dense/bordered Schur, block-chain, latent low-rank Schur, and supernodal selected-inverse helpers.
- Current prior assembly uses reserved sparse insertion, not the old fresh-triplet implementation described in historical notes.

## Run 1 Status

Saved source-grounded hypotheses with local discriminating checks. No numerical implementation was changed and no timing claim was made.

## Run 2: Covariance Representations And Admission Rules

### R2.1: Precision-First Descriptor Capabilities (Highest Family-Wide Potential)

Evidence: [representation compiler](../R/FUN_vsm.R#L1892), [precision builder](MNR.cpp#L4827). The compiler already distinguishes mathematical representations, but the precision builder always evaluates a dense factor covariance before using those capabilities. Native stationary AR precision also builds dense innovation matrices and multiplies them; accepted-state derivatives default to dense `-P dK P`.

Extend the existing descriptor contract, rather than adding family switches to the optimizer: optional precision application/materialization, logdet, working-coordinate precision derivative, logdet derivative, declared structural support, and covariance derivative application. Keep the generic covariance evaluator and `ownm` fallback. A sparse/banded/low-rank factor state must be usable by tensor kernels, score contractions, and PD validation, not only coefficient assembly.

Candidate rules and exact kernels:

| Family/capability | Current behavior | Proposed mathematical improvement and guard |
| --- | --- | --- |
| `ism`, `dsm`, `atm`, section-owner diagonals | Dense factor matrices; reciprocal precision exists | Store diagonal vectors and use elementwise mode operations. Do not infer effective diagonal residual precision through non-diagonal weights. |
| Homogeneous AR1 | Direct tridiagonal precision and analytic derivative already exist | Dispatch derivative before generic multiplication; banded storage/application and analytic logdet derivative. At rho=0 derivative support still includes neighbors. |
| Heterogeneous AR1, AR2/AR3/AR(p) | Stationary-AR innovation precision path exists | Accumulate outer products of short innovation rows, plus initial p-by-p precision. O(q p^2) assembly rather than dense products; differentiate PACF recursion and innovation variance analytically. Heterogeneous scaling preserves bandwidth. |
| ANTE | Direct `T' D^-1 T` formula exists, but dense matrices/products and generic derivatives remain | Accumulate each short modified-Cholesky row outer product, apply T and T' sparsely, differentiate T and D on working scale. Admit only declared/validated short support; unrestricted T is not banded. |
| CSM homogeneous/heterogeneous | Closed-form inverse; positive-rho Woodbury latent state already exists | Apply diagonal-plus-rank-one precision without dense q-by-q multiplication. Preserve negative-rho closed-form inverse; positive latent-specific representation does not extend to negative rho. |
| FA/RR | Woodbury precision, native derivatives, and automatic latent Schur already exist | Apply low-rank correction directly and contract derivative factors rather than materializing every dense derivative. Remove marginal prior assembly from latent iterations only after all consumers are audited. RR here includes an identity-specific term; it is not a singular pure rank-k covariance. |
| US/covm, CORG | Generic SPD precision and transformed-factor derivatives | Reuse supplied Cholesky coordinates for logdet and triangular solves; rank-one/rank-two derivative representations where valid. US anchoring and CORG parameterization differ, so no shared formula based solely on the family label. |
| MA | Banded covariance, generally dense inverse | Banded covariance Cholesky and solve/logdet; optionally a separately proved state-space representation. Never route it to banded-precision elimination merely because covariance is banded. |
| General Toeplitz | R evaluator, q-1 reflection parameters, numeric derivatives | Durbin-Levinson covariance/logdet recursion and stable Toeplitz solves or structured inverse products, O(q^2) candidates. Derivatives and all q-1 parameters can still dominate. General Toeplitz precision is not short-banded. |
| SAR | Dense R solve for B^-1, then generic inversion | For B=I-rho W, M=B^-1 B^-T, c=M[1,1], K=M/c: `K^-1=c B'B`, `log|K|=-2 log|det B|-q log c`. Compute c as squared norm of B^-T e1, not a full inverse. B'B may fill even if W is sparse; matrix-free application avoids that fill. |
| Proper CAR | Dense R inversion of Q=D-rho W, then generic inversion | With c=(Q^-1)[1,1], `K^-1=c Q`, `log|K|=-log|Q|-q log c`; obtain c by one solve. Reuse sparse Q factorization and selected-inverse traces for normalization derivatives. No isolated levels under the current constructor. |
| Matern/metric/custom dense | Generic dense covariance and callback/analytic derivatives | Preserve dense reference; do not assume Markov precision for arbitrary geometry or smoothness. SPDE/Vecchia or low-rank truncation changes the model and needs a separate explicit approximation interface. |

For CAR, `dc/drho = e1' Q^-1 W Q^-1 e1`. For SAR, differentiate v=B^-T e1 using `dv/drho=B^-T W' v`, so `dc/drho=2 v' dv/drho`. Multiply every natural-coordinate derivative by the actual working-coordinate transform derivative. Normalization c is essential: omitting it silently changes the model and reported scale.

Cheap falsification: reconstruct covariance from the precision, compare logdet and central-difference derivatives on small matrices. Include nonsymmetric SAR W, heterogeneous AR scales, negative CSM rho, PACF near boundaries, and custom factors. Then compare accepted likelihoods, scores, AI, estimates, and uncertainty, not just BLUE/BLUP.

### R2.2: Extend Path Elimination To Banded Block Graphs

Evidence: [chain admission](MNR.cpp#L7450), [path ordering](mme_paths.h#L177). The current engine accepts only collections of paths, so genuine AR2/AR3 and higher-order ANTE can be rejected even when their precisions are short-banded.

Proposal: exact block-banded Cholesky/LDLT and selected inversion, with covariance-coordinate blocks as vertices and all other coefficients in the border. Admit from the ASSEMBLED coefficient graph and bounded elimination width, not just an AR label. Data overlap, correlated residuals, and other terms can destroy the prior's bandwidth. Budget actual band storage, border couplings, inverse blocks, and derivative trace support.

For q blocks of width n and bounded block bandwidth p, a target factorization cost is O(q p^2 n^3), versus O((q n)^3) for a coupled dense group. This is an asymptotic target, not a benchmark claim. Derivative support must be unioned with needed inverse entries; exact-zero parameters do not remove derivative neighbors. Nonlocal requests must retain exact solve fallback.

Cheap falsification: AR2/3 and ANTE small dense references; disconnected/reordered graphs; missing records; branches or long-range couplings; zero correlations and changing support; tiny budgets. Extend existing [chain tests](../tests/testthat/test-mme-block-chain.R#L1) and native policy tests.

### R2.3: Cost-Based Exact Missing-Cell Residual Planning

Evidence: [grid admission](mme_paths.h#L17), [grid missing-cell factorization](MNR.cpp#L9672). Admission uses rectangle/observed ratios and a hard 4000-missing-cell bound. Missing-cell corrections can consume dense m-by-m factors and rectangle-by-m columns even when the full Kronecker precision is sparse.

Proposal: compare exact observed covariance block factorization against missing-precision Schur correction using estimated fill, dense/sparse workspace, and factor/application/trace costs. For short-banded factors, Q_mm is a principal submatrix of sparse Kronecker precision: use sparse factorization/selected inversion where advantageous, plus sparse or on-demand Q_om products. Cache equivalent missing masks and section states within an evaluation. Dense FA/US factors should not use the same sparsity assumption.

Identity to preserve: `R_oo^-1=Q_oo-Q_om Q_mm^-1 Q_mo`, and `log|R_oo|=log|R_full|+log|Q_mm|`. Taking Q_oo alone conditions on missing cells and is incorrect. Ragged observed covariance is not generally Kronecker.

Cheap falsification: compare against direct principal covariance inversion on several missing masks, including highly fragmented and nearly empty blocks; compare score derivatives as well as logdet. Reject the sparse strategy when fill/workspaces exceed the budget. Weighted/sectioned unsupported paths remain guarded.

### R2.4: Improve Schur Admission By Total Work And Feasible Alternatives

Evidence: [cost policy](mme_paths.h#L106), [chosen-plan density rejection](MNR.cpp#L7565). The current planner estimates dense work, chooses border or merge, then checks density; rejecting the chosen strategy does not retry the alternative. Latent/chain candidates also precede the ordinary planner without a common total-cost comparison.

Proposal: generate a small set of exact feasible plans; include numeric factorization, triangular solves, demanded inverse/trace blocks, RHS count, assembly, and storage. Validate density/budget BEFORE choosing. Initially retry the other ordinary plan if the cheapest fails density; later consider component-local combinations and comparing latent/chain costs. Separate REML trace cost from ML and all-fixed solve-only requirements.

Do not simply loosen the 0.1 density or 2000 border thresholds. Use diagnostics and measured crossover constants; dense workspace limits are not total process-memory limits. Stored numerical zeros and derivative support must not be confused with permanent structural absence.

Cheap falsification: native graph fixtures where the first plan fails density but the other succeeds; dense+pedigree mixtures, varying border size and trace demands. Compare end-to-end time and peak memory against generic CHOLMOD. Preserve deterministic ties and numeric-failure fallbacks.

### R2.5: Cache Fixed Shapes And Consider Common-Scale Profiling

Fixed-shape factors still pass through evaluation/precision construction each coefficient evaluation. Cache covariance/precision/logdet for parameter-invariant shapes, scaling separately by exp(log_sigma2). For repeated parameter states, use exact cache keys including sections, parameters, layout, and fallback state. User callbacks must not be assumed pure unless that becomes an explicit contract.

Longer-term exact profiling: if every covariance can be written V=t V0(ratios,shapes) and the common scale is unconstrained, the optimal scale is `t=y'P0 y/(n-rank(X))` for REML, or `/n` for ML. This removes one optimization coordinate and may reduce line searches. Fixed covariance components, working-coordinate affine constraints, scale floors, and uncertainty reconstruction require explicit eligibility rules. Treat this as a later reparameterization, not a new default; preserve delta-method reporting and profile curvature.

Cheap falsification: fixed-shape multi-iteration instrumentation should show one shape factorization and identical outputs. Profiled-scale tiny models must match joint optimization and independent dense objectives, including uncertainty after reconstruction. Decline profiling when the common scale is not free.

## Run 2 Status

Saved covariance-family opportunities and exactness guards. This pass confirms several older proposals are already implemented: stationary AR precision, positive CSM latent admission, FA/RR latent Schur, and native RR derivatives. The remaining opportunity is representation-aware execution and planning, not reintroducing those features.

## Run 3: Adversarial Review And Practical Priorities

### R3.1: PCG Residual Traces Can Avoid Btrace Assembly (High Priority)

Evidence: [residual score correction](MNR.cpp#L12590). Diagonal residual derivatives already use projected probe contractions. Other residual paths form `Btrace=F' dR F`, F=RiW, then walk its stored entries for a stochastic trace, including block-local dense products/triplets.

Exact identity for the SAME probes and solves: `sum z' F' dR F x = sum (F z)' dR (F x)`, x=C^-1 z. Compute FZ and FX once per accepted state, then apply each covariance derivative in observation blocks/tensor form. For complete/ragged residuals reuse the existing derivative application kernels. This preserves the estimator, unlike changing probes or tolerances, and avoids a potentially large nEffects-by-nEffects trace operand.

Cheap falsification: compare both contractions with identical deterministic probes on dense references and complete/ragged/weighted residual models. Track saved assembly, peak memory, and total time; nR-by-probe panels can themselves be large, so chunk by residual blocks. ML PCG remains unsupported.

### R3.2: Cache Grid Factors By Exact Level Subset And Section (High Priority)

Evidence: [grid per-block setup](MNR.cpp#L9651). Full factor covariances are evaluated once per section, but each residual block independently takes the same level subset, Cholesky factors it, and explicitly inverts the triangular factor. Thus identical subsets can repeat factorization within the same coefficient evaluation, and this path does not use the descriptor's native precision ladder.

First change: evaluation-local cache of factor covariance/precision/logdet keyed by factor, section, and ordered level subset. Reuse within the trial; invalidate numeric values next trial. Later permit native precision for a FULL factor level set and operator state for incomplete sets.

Critical refinement to Run 2: a nonconsecutive level subset is a principal covariance submatrix, not a principal precision submatrix. Homogeneous AR1 with gaps remains Markov on the ordered retained chain but has gap-dependent transitions rho^gap, not the unit-lag formula. General AR(p) subsets require proper marginalization and may lose short bandwidth. Use generic subset factorization unless an exact subset kernel exists. Missing-mask cache keys also require the subset/section/numeric state, not just the mask.

Cheap falsification: repeated identical grids must factor once per unique key; reordered/gapped level sets and different sections must not collide. Compare missing-cell precision, score, likelihood, and fit outputs against the current path.

### R3.3: Nyström Application Has Avoidable Repeated Work (Small, Low Risk)

Evidence: [preconditioner application](MNR.cpp#L8448), [batched recurrence](MNR.cpp#L8692). `D^-1 U` is recomputed on every preconditioner call although it is constant for the coefficient evaluation; after initialization, the assembled batch loop applies the preconditioner separately to each active RHS.

Cache `D^-1 U` with the current preconditioner state and apply the Woodbury correction to an active RHS panel. Retain the existing positive-diagonal floor and core-factor checks. Current Nyström landmarks are deterministic uniformly spaced indices, not leverage-score sampling; changing landmark selection is a separate benchmark-dependent proposal, not necessary for this optimization.

Cheap falsification: compare preconditioner panel application to the current scalar-column formula and to a dense inverse of D+UU'. Measure setup/application separately and preserve PCG true residuals.

### R3.4: Early-Converged Matrix-Free Warm Starts Need An Exit Check (Correctness Risk)

Evidence: [matrix-free initialization](MNR.cpp#L8546). After computing an initial residual, the code immediately requires positive `r' M^-1 r`; it does not first check whether the supplied Lanczos warm start already satisfies tolerance. A zero residual can produce an error, and a tiny residual can cause unnecessary fragile CG iterations. This is especially relevant when the Lanczos space spans a small system.

Proposal: finite true-residual convergence check before the first positivity test, preserving zero-RHS handling. Test both exact and near-exact guesses and a deliberately insufficient guess. The dense check below establishes the possibility of a near-zero residual, not that a package fit currently fails.

The assembled batch solver already masks initially converged columns; it checks true residuals at return. It can still stop with a true-residual error rather than restart if recurrence drift is excessive. Add reliable-residual replacement/restart only if a focused ill-conditioned fixture demonstrates the need.

### R3.5: PCG Likelihood/Score Accuracy Requires Diagnostics, Not Just Looser Rules

Evidence: [fixed-probe SLQ](MNR.cpp#L8844), [inverse probes](MNR.cpp#L8941), [likelihood acceptance](MNR.cpp#L10963), [stopping criteria](MNR.cpp#L13805). The likelihood uses finite-probe SLQ; score traces use solved Hutchinson probes. Common probes reduce noise across iterations but do not make these two estimators an exactly consistent objective/gradient pair. With a finite fixed probe set, this mismatch can persist even with exact matrix-log quadrature: the derivative of `z' log(C) z` is a Frechet derivative, not generally `z' C^-1 dC z`.

Dense two-dimensional check: C=[[3,1],[1,2]], dC=diag(1,0), z=(1,-1). Central-difference derivative of z'log(C)z = 0.572163576; z'C^-1 dC z = 0.600000000. Their expectation over isotropic probes agrees with the exact logdet derivative, but individual finite estimates need not.

Proposal: expose per-probe dispersion, nested Krylov-depth comparisons, true-residual norms, operator counts, and final score diagnostics. Separate solve tolerance from trace sampling and quadrature accuracy. Adaptive depth/probe schedules must keep common prefixes and re-evaluate accepted and trial likelihoods at matching accuracy; otherwise accuracy changes can mimic likelihood changes. Confirm final results at higher accuracy. A near-zero approximate likelihood change or small safeguarded step alone is not proof of an exact score root. Do not promise deterministic-solver precision from a small probe count.

### R3.6: Reduce Derivative Work Only With An Explicit Output Contract

Evidence: [accepted-state derivatives](MNR.cpp#L11106), [all-parameter sensitivity allocation](MNR.cpp#L11364). Fixed coordinates still receive derivatives, solves, score contractions, and AI entries; an all-fixed fit does not bypass these.

For routine optimization, solve only needed free/mapped sensitivity directions. For affine constraints, the independent directions can be fewer than nominal free coordinates. But `InfMat` and `dLu` are public outputs, and callers may request full-coordinate information. Introduce an explicit internal request for covariance information, retaining legacy outputs by default. Do not silently replace missing information with zeros.

The nearby [Wald helper](../R/FUN_wald.R#L110) deliberately sets all coordinates fixed, but derives X'V^-1X from C-based prediction solves, NOT from covariance `InfMat`. It can opt into a lightweight coefficient-evaluation mode that preserves C, its scale, solves/caches, and requested uncertainty. The existing known-covariance solve-only API is also relevant; do not duplicate it wholesale.

Cheap falsification: all-fixed fit output-contract tests, Wald/Satterthwaite/KR tests, fixed loading with free scale, and constrained affine maps. Preserve full information when requested; compare reduced-direction updates to the current algorithm, including constraints and safeguards.

### R3.7: Precision-First PD Checks Must Preserve The Numerical Floor

Evidence: [descriptor PD check](MNR.cpp#L5639). Every proposed factor state is evaluated and eigendecomposed, even for known positive parameterizations. Cholesky success or PACF admissibility alone is not equivalent to the present minimum-eigenvalue floor.

Use exact diagonal and CSM eigenvalue formulas when applicable. For FA/RR, min(specific) gives a sufficient lower bound; for heterogeneous CSM, min(sd^2) times the correlation minimum eigenvalue is a sufficient bound. A conservative bound certifying the scaled Kronecker minimum above the current floor can skip eigendecomposition. If inconclusive, retain the current exact eigenvalue check; do not reject a valid state merely because a lower bound is small. Guard nonfinite/overflow states and section minima. Generic US Cholesky positivity is not a floor certificate.

### R3.8: Final-Inverse Reuse Is Valid, But Not A Drop-In Export

Run 1's suspicion is confirmed at [final extraction](MNR.cpp#L14111): CHOLMOD computeCi=1 refactors C through Eigen LDLT. The comment saying there is no supernodal selected inverse is stale: [the helper exists](MNR.cpp#L8074), but is currently activated only when a random inverse-cache group is skipped and its storage fits the budget.

Next design: reuse/build the final supernodal subset for ordinary CHOLMOD fits and export values in original ordering. Decide explicitly whether mode 1 promises Eigen's precise fill pattern or a solver-specific selected subset; current public output uses the Eigen pattern. Preserve diagonals, PEV, `CiComputed=FALSE`, and exact solve-based arbitrary prediction. A Schur/chain/latent fit has no reusable marginal CHOLMOD factor and needs its own subset exporter or the existing fallback. This is not universally a one-line optimization.

For LDLT, final extraction also rebuilds a selected inverse already computed at the final accepted state. Reuse the numeric subset when still valid instead of repeating Takahashi; avoid introducing another deep topology copy. Full inverse mode has unavoidable quadratic OUTPUT size; batching identity solves reduces transient allocations but does not make that output scalable.

## Three-Run Synthesis

| Topic | Run 1 | Run 2 | Run 3 conclusion |
| --- | --- | --- | --- |
| PCG multiple RHS | Matrix-free path remains scalar | Structure-aware factor applications | Batch independent CG/operator panels; cache Nyström scaled basis first. |
| Precision/derivatives | AR1 specialization overwrites generic work | Precision-first AR/ANTE/SAR/CAR and low-rank kernels | Start with exact AR1 dispatch; preserve normalization, derivative support, subsets, and PD floor. |
| Traces | Existing inverses/probes must be reused | Cost depends on demanded inverse blocks | PCG residual traces can avoid Btrace with the same estimator; CHOLMOD subset export needs a contract audit. |
| Exact engines | ML can reuse D plus fixed Schur | Banded graph and better grid/Schur planning | Implement behind narrow capability gates with exact fallback; do not loosen thresholds blindly. |
| Optimization | Preconditioned SLQ may reduce work | Fixed-shape caching and common-scale profiling | Instrument stochastic error first; profile scale later, after uncertainty/constraint design. |

Best low-risk starting package: existing AR1 derivative dispatch, factor/subset reuse, direct residual probe contractions, Nyström application caching, and the matrix-free early-convergence guard. Largest structural opportunities: banded factor operations and block elimination, sparse SAR/CAR precision, and D-first ML. Which dominates depends on covariance dimension, relationship density, observation layout, RHS count, and inverse-output demands.

## Run 4: LDLT And CHOLMOD With Arbitrary-Factor Products

This additional pass prioritizes exact direct solvers. Earlier PCG proposals remain useful but are not prerequisites for this direct-solver work. The planner must support any number of factors supplied by `vsm()`, not special-case one-factor or two-factor products. No spectral rotation is proposed.

### R4.1: Partial Diagonal Products Expose Exact Smaller Groups (Highest Structural Priority)

Evidence: [vsm factor-list compilation](../R/FUN_vsm.R#L243), [whole-product diagonal flag](MNR.cpp#L6455), [random inverse groups](MNR.cpp#L7289). The R compiler retains all shaping factors in order; the direct planner splits into coordinate groups only if EVERY shaping factor is structurally diagonal. Otherwise its ordinary candidates contain the whole random term, even when some axes remain independent.

For example, `dsm(site) x ar1m(time) x usm(trait) x ism(id)` has independent site slices in its PRIOR, despite not being a diagonal product. With q_site, q_time, q_trait levels and n IDs, each prior component has q_time*q_trait*n coefficients, not q_site*q_time*q_trait*n. Multiple diagonal axes multiply the number of candidate components. Block-diagonal custom factors can offer the same opportunity without being diagonal.

General rule: compute connected components of each declared precision-support graph, including the relationship precision graph. Their Cartesian component tuples partition the tensor prior when the factor diagonals are present, as they are for SPD factors. Map tuples to actual coefficient indices using the existing ordered-product convention, with the relationship/ID coordinate fastest. The FULL assembled data term W'RiW can merge those components or create border couplings; verify it before choosing independent elimination. Fixed effects can remain in a Schur border rather than forcing all components to merge.

This is not a claim that a sum of Kronecker terms is itself Kronecker. Use factor supports to propose candidates, then the assembled graph to admit exact plans. Reject any candidate based only on a temporarily zero rho/loading as a permanent partition; trace derivative support can cross that temporary partition.

Cheap falsification: three-, four-, and five-factor products with diagonal axes in different positions, connected/disconnected `Ai`, overlapping incidence, and correlated residuals. Compare complete C, solves, both likelihood modes, score, AI, and PEV to the current generic backend. Require candidate reduction AND actual fill/time/memory reduction.

### R4.2: Symbolic Support And Numeric Values Need Separate Lifetimes

Evidence: [pattern comparison](MNR.cpp#L7077), [prior insertion](MNR.cpp#L10514), [cached direct analysis](MNR.cpp#L10633). Current C/D analysis is already reused on identical stored patterns. However, numeric-zero omission in assembly can change the stored pattern, forcing analysis, selected-inverse topology rebuilds, and engine re-planning.

Proposal: compile a stable structural envelope where it is affordable: ordered tensor support and coefficient destination offsets, plus data support appropriate to the residual/layout. Refresh numeric values without re-sorting/rebuilding the same graph. Retain the current variable-pattern path when the envelope would be substantially larger. Explicit zero slots must survive the chosen assembly operations; ordinary sparse addition/pruning must not accidentally destroy the envelope.

For m factors, support enumeration follows products of factor nonzeros, NOT loops over all q^2 entries, q=product q_f. Track covariance support, precision support, first-derivative support, and elimination fill separately. A scale change preserves support; a dense free factor, changing residual precision, missing-cell marginalization, or a generic callback may require conservative support or invalidate the plan. Cache keys include factor order/dimensions/support, actual design, sections, weights, permutation, and requested trace/output work.

Important tradeoff: derivatives can require entries outside the current C pattern (e.g. AR rho=0). Such entries need a selected-inverse closure or exact solve fallback, not necessarily extra structural fill in C itself. Stabilize C only where saved analysis outweighs added factorization/storage; do not take the union of every conceivable derivative edge blindly.

Cheap falsification: parameter trajectories crossing zero and repeated section layouts; count symbolic analyses and selected-topology rebuilds. Compare ordering fill and total elapsed time. Earlier naive assembly caching was not a measured win, so a new representation must earn its place.

### R4.3: Compile LDLT Selected-Inverse Dependencies Once

Evidence: [Takahashi builder](MNR.cpp#L6929). Numeric topology reuse still copies the sparse L into `Lmat`, reconstructs neighbor/value vectors per column, and does hash/binary-search lookup for every neighbor pair inside the numeric recursion. The dependency addresses are fixed for a fixed factor topology.

Proposal: compile neighbor-to-L-value offsets and selected-inverse dependency offsets once per symbolic topology. Numeric iterations refresh L/D values and run direct-index recursion. For small filled neighbor sets use cached pair-address panels; for large sets bound the address-cache size or use a tiled/batched gather strategy. A dense local contraction `Z_NN * l_N` can replace scalar lookup-heavy sums when the gather and temporary cost is justified.

This is an exact implementation change, not a new inverse approximation. The offset table can cost sum degree(L)^2, much larger than nnz(L), so impose a separate budget and retain the current lookup recursion when compiling it is expensive. Investigate direct access to Eigen's factor view without assuming undocumented storage or lifetime guarantees; eliminating the L copy is conditional on a supported API.

Cheap falsification: compare every represented selected entry, diagonal, and fallback trace against dense inverses and current Takahashi. Verify permutation and topology invalidation; measure numeric recursion time INCLUDING address-cache setup and memory. Reuse final numeric values as proposed in Run 3.

### R4.4: Choose Ordering By Total Direct-Solver Work

Evidence: [LDLT ordering alias](MNR.cpp#L18), [forced CHOLMOD supernodal mode](MNR.cpp#L7210). LDLT selects METIS at build time when available, otherwise AMD. CHOLMOD uses its own analysis and is forced supernodal. Neither is currently selected from an end-to-end factor/trace cost comparison in this function.

Proposal: initially measure factor nnz, frontal/supernode dimensions, analysis time, numeric factorization, sensitivity solves, and selected inversion under supported ordering alternatives. Prefer an ordering that minimizes total fit cost, not just nnz(L): selected-inverse work depends strongly on filled degrees, and CHOLMOD benefits from BLAS-friendly fronts. For products of several short-banded axes, use the product graph to suggest nested-dissection separators; two or more AR axes produce a multidimensional graph, not a simple chain.

Do not automatically analyze every ordering on every trial. Use a bounded first-analysis candidate set, cache the winning plan, and re-evaluate only after meaningful structural change. Test CHOLMOD automatic/simplicial mode only with capability-aware inverse kernels: the current raw supernode cache/helpers assume a supernodal factor. A mode switch without a different trace/cache path is unsafe. Respect explicitly selected LDLT/CHOLMOD; do not silently substitute a different backend.

Cheap falsification: sparse pedigree, product-grid precisions, dense genomic blocks, and multiple coupled terms. Account for first-analysis amortization over short/long fits and reject alternatives that save factorization but increase trace or solve cost.

### R4.5: CHOLMOD Inverse Work Should Follow Demand And Reach

Evidence: [forward unit solves](MNR.cpp#L7999), [inverse-block preparation](MNR.cpp#L8330), [selected-inverse trigger](MNR.cpp#L8364). The forward helper skips inactive supernodes but initially zero-allocates n-by-m workspace; accepted-state cache preparation forms whole inverse blocks if each individual group fits the budget. Supernodal selected inversion is attempted only when a random group is skipped.

Proposal: compare three exact trace strategies before cache preparation: selected inverse covering requested entries; full inverse blocks for repeatedly demanded dense groups; and batched exact solves for unsupported/rare entries. Budget the AGGREGATE retained blocks, fixed/random cross blocks, selected values, and simultaneous workspace, not only one group at a time. A fast path need not materialize every inverse block merely because it fits individually.

Precompute structural supernode reach per RHS group and use touched-row or chunked forward workspaces where suitable. For whole inverse blocks, chunking must retain cross-panel Gram products in Y'Y; independently accumulating only panel-diagonal products would omit inverse cross entries. Selected inversion has its own potentially large frontal workspace and should be costed rather than enabled universally.

Generalize the same strategy to ML's D and to generic residual R. Currently ML CHOLMOD uses solve fallbacks, and generic CHOLMOD R has no selected subset at its trace call site. Factor-specific selected state must keep independent permutations and lifetime/refresh rules for C, D, and R. Defer residual selected inversion until trial acceptance where it is only needed for scores; logdet and precision applications must still occur on trials.

Cheap falsification: sparse random support, dense groups, many groups individually fitting but jointly exceeding budget, ML, irregular residuals, and subset misses. Compare exact contractions and memory. For computeCi=1 preserve or explicitly resolve the public subset contract from Run 3; arbitrary predictions still use solves.

### R4.6: Schur Updates Can Use One Triangular Solve Plus SYRK

Evidence: [ordinary dense Schur factorization](MNR.cpp#L7773). The current code computes F=D^-1 B using two triangular solves and updates S by B'F, restricted to active border columns. That active-column optimization is already implemented and should not be proposed again.

For D=L L', let T=L^-1 B. Then `B'D^-1 B=T'T`, allowing one triangular solve followed by a symmetric rank-k update (SYRK). F=L^-T T can be computed later or retained if subsequent joint solves/inverse caches require it. Store only the requested triangle during update and symmetrize at the interface as needed. Apply this selectively to independent and latent blocks; chain pivot updates need their own derivation and dependencies.

If F is immediately required anyway, the two triangular solves still occur; the gain is the symmetry-aware BLAS update, not eliminating a required solve. If T and F are both retained, extra storage can negate the gain. Use existing Armadillo/BLAS facilities and measure whether they already dispatch symmetrically for the expression used.

Cheap falsification: exact Schur, determinant, joint solve, cross-block inverse, and latent trace references; small/large active borders; include extraction/allocation cost in timings.

### R4.7: Generalize Latent/Banded Engines Along Tensor Axes, Not Factor Counts

Evidence: [single-factor low-rank state guard](MNR.cpp#L4968), [latent admission](MNR.cpp#L7368). Low-rank availability is retained only for a one-factor descriptor. Additional identity, diagonal, AR, or dense factors therefore disable that representation even if an exact conditional construction exists.

For one low-rank-plus-diagonal axis K_j=D_j+L_j L_j' within ANY ordered factor list, permute that axis to the logical outside and write covariance `s K_j x H`, where H is the product of all remaining shaping covariances and relationship covariance. With P_H=H^-1, exact augmented precision is

`s^-1 [D_j^-1, -D_j^-1 L_j; -L_j' D_j^-1, I+L_j' D_j^-1 L_j] x P_H`.

Eliminating virtual coefficients recovers `s^-1 K_j^-1 x P_H`. If the low-rank axis has width q_j and rank k_j, virtual width is k_j*r with r=n*product_(f!=j) q_f, NOT simply k_j*n. Interior groups also have width r unless other diagonal/component structure splits them. The determinant correction is the logdet of the ACTUAL virtual-prior block: `k_j*r*log(s^-1)+r*log|I+L_j'D_j^-1L_j|+k_j*log|P_H|`.

Keep this as a factor-level representation and logical permutation, not a densely materialized H. Admit only when assembled data/residual coupling, actual augmented fill, trace access, and total storage make it worthwhile. Existing equal-width/disjoint-incidence guards need a mathematically justified generalized counterpart, not removal. Products containing multiple low-rank axes require separately derived staged augmentation or an augmented tensor construction; selecting one eligible axis is the first implementation. Multiple low-rank axes do NOT imply that the full product is diagonal plus rank sum(k_j).

For banded axes, evaluate elimination width in the full product graph. Diagonal axes can split it; other correlated axes enlarge blocks or create grid separators. A generic graph plan should choose a viable axis/component hierarchy and leave dense/nonseparable cases to the requested sparse solver. No arbitrary-factor count cutoff is needed; resource/coupling/capability limits are.

Cheap falsification: explicit three/four-factor covariance and augmented precision references, low-rank axis in first/middle/last positions, diagonal and correlated remaining axes, and virtual logdet correction. Compare marginal original-order C, scores, free loading/shape AI, PEV, and fallback behavior. Reject plans that inflate virtual borders beyond practical budgets.

### R4.8: Factorwise Contractions Still Contain A Product-Size Ceiling

Evidence: [Kronecker contraction](MNR.cpp#L5543), [descriptor factor assembly](MNR.cpp#L5231). The helpers accept arbitrary factor lists, but contraction walks all q^2 entries and decodes m digits per pair. Direct random covariance/precision and quadratic/trace matrices also remain product-sized in places. Supporting m factors syntactically does not establish scalable memory use.

Proposal: contract using structural support tuples and precompiled mixed-radix maps for diagonal/banded factors; retain successive tensor contractions for dense axes, and low-rank contractions for advertised factors. First/second derivative contractions should replace only affected factors and share unaffected products. Keep sparse assembly lazy until a direct backend needs entries, with overflow-checked dimension/nonzero/workspace estimates before allocating.

Dense exact C can still have quadratic storage and cubic factorization: tensor priors alone do not eliminate that ceiling when W'RiW destroys structure. Dense contractions may still require O(q^2) input reads even after avoiding per-entry digit decoding. Performance claims must state the layout and sparsity gates.

## Run 4 Status

Additional direct-solver findings saved. Highest new priorities: partial-factor component planning, LDLT compiled selected-inverse dependencies, CHOLMOD demand/aggregate-budget trace planning, and symmetry-aware Schur updates. Arbitrary-factor latent/banded execution is a gated architectural extension, not a blanket replacement for generic direct factorization.

## Revised Implementation Plan

The following direct-solver-first sequence supersedes the earlier ordering. Arbitrary-length `vsm()` factor lists are a requirement at every phase, not a later special case. Earlier PCG items remain a separate follow-up track.

1. **Baseline direct-solver costs and contracts.** Instrument covariance/derivative kernels, assembly, symbolic/numeric factorization, selected inversion, inverse-block/solve fallbacks, sensitivities, and final outputs. Record aggregate memory, factor fill/front sizes, RHS counts, and ordering. Compare LDLT and CHOLMOD on the SAME current source with fixed thread settings, fixed iterations and converged fits. Establish subset/information/solver-selection contracts before altering them.
2. **Introduce an arbitrary-factor structural plan.** Retain ordered factor capabilities, dimensions, strides, component labels, precision and derivative support, and section/layout ownership. Include relationship precision as the final tensor axis. Overflow-check products and predict support/fill/workspaces. Build no full tensor covariance merely to classify a path. Keep numeric state separate from this plan and preserve public coefficient ordering via explicit mappings.
3. **Land low-risk exact direct optimizations.** Dispatch AR1 derivatives before dense generic products; cache repeated grid subset factors; reuse final LDLT selected values; prototype budgeted compiled Takahashi addresses and triangular/SYRK Schur updates. Land each independently after focused equivalence and timing checks. Do not remove fallback kernels.
4. **Exploit partial-factor components and stable affordable support.** Generate prior components from factor graphs, check assembled data coupling, and compare feasible border/merge plans before density rejection. Stabilize numeric assembly/analysis only when envelope fill is economical. Add trace/output demands and aggregate workspace to costs. Make this planner applicable to both direct backends; any new LDLT structured execution must preserve LDLT selection semantics and uncertainty contracts.
5. **Make exact trace preparation demand-driven.** Select among LDLT selected entries, CHOLMOD supernodal selected inversion, full blocks, and batched solves using required supports, reach, reuse, and aggregate memory. Add separate C/D/R state and correct invalidation. Audit CHOLMOD computeCi=1 export; retain current fallback until output compatibility is established.
6. **Add precision-first and tensor-aware kernels.** Diagonal, AR(p), ANTE, normalized CAR/SAR, and low-rank factor applications/derivatives first; generic/custom dense factors retain reference execution. Use support-aware contractions and sufficient PD certificates with exact floor checks when ambiguous. Validate mixed products with 3-5 factors, not only one/two.
7. **Generalize exact elimination by graph width and augmentation cost.** Add block-banded execution on feasible tensor axes, then product-grid separator plans; extend latent augmentation to one eligible low-rank axis with arbitrary remaining factors and the full virtual determinant correction. Multiple low-rank axes are a separately proved extension. Preserve incidence/residual/ML gates until justified; use generic direct fallback for nonseparable or over-budget layouts.
8. **Optimize ML and ordering after structural kernels stabilize.** Prototype D-first fixed-border Schur evaluation and extend D-specific traces; compare supported orderings by total factor/trace/solve cost with first-analysis amortization. Do not switch CHOLMOD away from supernodal without compatible inverse kernels. Preserve public uncertainty and full inverse-output behavior.
9. **Complete regression and performance gates before default dispatch changes.** Exercise all product positions, many-factor/multiple-term models, missing/gapped/sectioned layouts, non-diagonal weights, zero-crossing support, constraints, REML/ML, final inverse modes, and prediction/Wald consumers. Require reproducible end-to-end benefit and bounded peak memory. Report which plan ran and why alternatives were rejected. PCG batching/accuracy and common-scale profiling remain later independent work.

### Arbitrary-Factor Validation Matrix

- Single term versus multiple random terms; 0 shaping factors (identity), 1, 2, 3, 4, and 5 factors, with no API maximum imposed by these examples.
- All diagonal; diagonal x banded x dense; two or more banded axes; one low-rank axis x several other factors; multiple low-rank axes; disconnected custom/block-diagonal factors; sparse/dense/disconnected relationship precision.
- Permute factor order while explicitly permuting coefficients/design/covariance consistently. Same named factor list in another order is not interchangeable without that mapping. Validate mixed-radix strides and first/middle/last derivative ownership.
- Verify arbitrary-m formulas: inverse of `s tensor K_f`, logdet `q log(s)+sum (q/q_f) log|K_f|`, each factor derivative and cross-factor second correction, and relationship multiplicities. Validate augmented logdet against the actual virtual block for generalized latent plans.
- Apply assembled-graph guards after W/residual/weight transformations. Include fixtures where a prior partition is destroyed by observation incidence or residual coupling; those must merge, border, or fall back exactly.
- Test resource rejection before product allocation, cache invalidation by shape/support/permutation/layout changes, derivative neighbors at zero parameters, and no accidental densification during classification or admission.

## Separate Follow-Up Track

Earlier PCG proposals remain independent follow-ups: projected-probe residual traces, matrix-free operator/solve batching, Nyström application caching, early-converged warm-start handling, and accuracy diagnostics before adaptive or preconditioned SLQ. Common-scale profiling and advanced Toeplitz/MA kernels remain separately gated prototypes. These do not block the revised direct-solver plan.

### Required Validation Gates

- Exact kernels: central-difference derivatives and dense inverse/logdet references; working-scale transforms, normalized covariance, rho=0/negative/near-boundary values, heterogeneous scales, missing/gapped coordinates, and singular/budget rejection.
- Fit equivalence: all three requested backends where supported, REML/ML, computeCi 0/1/2, fixed/free parameters, constraints, multiple random terms, weights, complete/ragged/sectioned residuals. Compare accepted likelihood paths for unchanged optimization rules, scores, AI, covariance estimates, BLUE/BLUP, SE/PEV, and prediction; allow justified roundoff only.
- Approximate changes: separately evaluate probe/depth sensitivity, solve residuals, final scores, stochastic error, and agreement with direct solutions. Do not apply exact path-parity criteria to a changed estimator, or present estimator changes as purely algebraic speedups.
- Existing suites to reuse: native MME policy tests; block-chain, FA latent-Schur, dense-Gu Schur, spatial/weighted residual, covariance safeguards, REML-vs-ML, multivariate, vcc, Wald, and prediction inverse-mode tests. Add focused cases in these files where practical.
- Performance: representative sparse/dense relationships, small/large covariance dimensions, few/many free coordinates, varying fixed border, low/high missing fractions, repeated grids, and final-inverse modes. Reject a new default if representative regressions or memory growth outweigh gains. Preserve user-selected solver semantics and report fallback reasons.

## Checks Actually Run

- Existing native policy test compiled and passed: `g++ -std=c++11 tests/native/test-mme-paths.cpp -o /tmp/sommer-ai-mme-review-paths`, then executed the binary.
- Seven independent dense R checks passed: ML Schur solves/logdet; symmetric diagonal preconditioning determinant; AR1 working-coordinate precision derivatives at rho=0,-0.45,0.7; normalized nonsymmetric SAR precision/logdet/normalization derivative; normalized CAR equivalents; missing-cell Schur precision/logdet; residual probe contraction with versus without Btrace. Deterministic seed 41; absolute tolerance 1e-8, AR1 finite differences 1e-7.
- A separate dense fixed-probe check confirmed objective/score mismatch as reported in R3.5; an exact C^-1 z warm start gave residual norm 1.57e-16, supporting the early-exit test case.
- Run 4 dense references passed for 1-4 shaping factors: factorwise inverse, determinant multiplicities, and mixed-radix factor-order permutations. Seed 73, absolute tolerance 1e-8.
- Generalized latent references passed for four shaping factors plus a relationship matrix, with the low-rank axis first, middle, and last: marginal prior precision, augmented joint solve, and virtual-prior determinant correction, including an added dense data term. The triangular/SYRK Schur identity also passed. These are algebra checks, not implemented-path benchmarks.
- A separate tensor graph check found four size-12 prior components from two diagonal axes; adding a cross-component data edge merged two as expected. Five-factor precision and logdet derivatives matched central differences on every factor axis (step 1e-5, absolute tolerance 1e-7).
- Report source links and saved pass headings checked. These checks validate identities and current policy behavior, NOT prospective performance or implemented solver changes.
- No package build, full R test suite, profiling, or before/after solver benchmark was run: numerical code was not changed. Shell `rg` was unavailable; workspace search tools were used. Historical performance notes were treated as leads, not evidence for the current version. The old benchmark scripts referenced by existing path documentation are absent from this checkout; the plan must recreate fixtures or reuse current test data.

## Review Status

All four reviews are saved, reconciled, and prioritized. The revised direct-solver plan explicitly supports arbitrary-length covariance products. Current implementation, validation, and remaining work are recorded in the status section at the top of this document.
