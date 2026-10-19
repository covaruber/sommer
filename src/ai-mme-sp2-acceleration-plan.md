# `ai_mme_sp2`: Prioritized Improvement Roadmap

**Target:** `src/ai_mme_sp2.cpp` (`ai_mme_sp2`; shared support in `src/mme_sp2_support.h`)
**Status reconciled:** local working tree, 2026-10-18. Historical measurements
below retain their original dates and baselines; implemented does not mean
committed, automatically enabled, or validated on every platform.

This roadmap replaces the separate review and acceleration-plan notes. It ranks remaining work by likely impact, implementation risk, and model coverage. Progress is tracked in the checklist below; completed items are not to be redone (see also `CHANGELOG`, version 4.4.90).

## Current Status (2026-10-18)

This section supersedes older descriptions of remaining PCG work. Historical
sections below remain evidence, not a current task checklist.

### AI work-reduction continuation (2026-10-18)

This separately approved continuation follows the closed acceleration batch.
Opt-in `experimental=list(reducedAI="force")` now packs independent `vcc`
directions before native sensitivity solves and assembles reduced exact AI
directly. Eligible paths are Gaussian deterministic Henderson REML/ML with
default AI-first recovery or positive EM/PXEM direction blending. The existing
optimizer retains feasible directions, trust units, regularized recovery,
likelihood acceptance and score certification. Residual working buffers contain
only independent directions affecting the residual covariance.

`covarianceInformation="full"` retains a full exact final evaluation and
existing information/uncertainty outputs. `"none"` uses exact reduced information
for certification without a report-only full refresh. Coefficient prediction,
PEV and `computeCi` are independent. Default admission remains off; unsupported
force requests error and auto requests expose fallback reasons.

Legacy explicit `blendWeight=0`, acceleration, latent RR, PCG, PQL, solve-only,
saddle-point and factor-score augmentation remain outside this implementation.
Full score traces and active local curvature-pair contractions are unchanged.
This reduced-AI slice did not implement common-scale profiling; the separately
requested continuation below now does. Reduced buffer counters are not process
peak-memory measurements.

The preserved-worktree comparison passes 16 default-off fits and 11 output
fields bit-for-bit. The isolated final macOS build and focused four-worker
run pass 73 test blocks and 4,698 assertions without failures, errors, warnings
or skips. Coverage includes independent full/reduced AI comparisons, REML/ML,
LDLT/CHOLMOD, random/residual and ratio ties, fixed/empty/boundary directions,
singular-information recovery, compression, positive EM/PXEM blending and
coefficient/covariance inference. Updated help parses and the explanation-only
vignette renders. These checks do not establish Windows or CUDA coverage.
On a small 24-tied-term fixture, five paired repetitions
give median sensitivity columns 175 to 39 with full reporting, and 175 to 14
with information omitted. Full-reporting elapsed medians are 0.172 s off and
0.181 s reduced; omission medians are 0.161 s off and 0.152 s reduced.
These workload-specific results demonstrate saved derivative work, not a
general whole-fit speedup or calibrated automatic admission policy.
`benchmarks/validate_reduced_ai.R` reproduces the comparison using isolated
baseline/candidate libraries.

### Common-scale profiling continuation (2026-10-18)

**Delivered, default off:** `experimental=list(commonScale="force")` uses the
exact family `V=s V0(psi)`, estimates `s=y'P0y/d` analytically for every shape
candidate, and optimizes only independent shape coordinates. Here `d=n-rank(X)`
for REML and `d=n` for ML, retaining original observation counts under compression.
The residual product scale supplies the reference. Existing affine `vcc`
constraints are incorporated into the shape tangent, including positive ratio
ties and random/residual ties. Enabling `reducedAI` separately is unnecessary.
Profile scores and the AI scale Schur complement use the same coefficient
factorization. The AI Schur metric is not claimed to equal the exact profiled
Hessian.

Admission requires Gaussian deterministic Henderson REML/ML, full-column-rank
fixed design, all product-level variance scales free, and constraints preserving
the common scale. Fixed shapes are allowed. Positive blending, PCG, PQL,
solve-only, acceleration, latent reduced rank, saddle-point and factor-score
augmentation remain excluded. Omitted and explicit zero `blendWeight` select
the dedicated profile optimizer, not AI-first EM recovery. It uses existing
trust units, regularized AI, observed-likelihood backtracking and exact final
score certification. Auto admission exposes static fallback reasons without
hiding runtime errors or retrying stalled fits with another optimizer.

**Boundary scope:** this version profiles interior scales only. Reconstructed
product variances must exceed `max(1e-8,tolParInv)`. Initial violations error;
trial violations backtrack. A stall warns and returns the last accepted point
unconverged. Constrained boundary-scale profiling and EM/PXEM recovery inside
this optimizer remain unimplemented; refit with `commonScale="off"` for them.
Aliased fixed designs are rejected instead of silently using the wrong REML
degrees of freedom.

One final ordinary same-point evaluation restores original coordinates and
coefficient uncertainty. Full reporting retains the entire covariance
information, including scale/shape cross-information and covariance standard
errors. Information omission still performs a final evaluation, using reduced
original constraint coordinates when smaller, and preserves coefficient PEV.
Wald/KR fixed-parameter perturbations explicitly bypass profiling.
`engineDiagnostics$commonScale` reports total evaluation/sensitivity work,
scale, degrees of freedom, backtracking, regularization and certification.
Optional phase seconds aggregate all trials and final reporting, while other
native profile counters describe the final evaluation.

Proof includes independent scalar scale maximization; finite-difference
profile scores; same-point full-AI Schur comparisons; REML/ML and LDLT/CHOLMOD;
correlated/weighted residuals; fixed shapes and ratio ties; compression;
zero-dimensional shape fits; singular shape regularization; iteration limits;
admission errors/fallbacks; full/omitted information; and coefficient,
covariance and Wald/KR inference. The four-worker focused run passes 91 test
blocks and 5,101 assertions without failures, errors, warnings or skips.
The preserved pre-profiling build agrees bit-for-bit for 16 default-off fits
and 11 output fields. These are macOS checks, not Windows/CUDA certification.

`benchmarks/validate_common_scale.R` checks converged-fit parity and records
five alternating paired repetitions per REML/ML and information mode. On its
216-observation, three-parameter AR1 fixture, total sensitivity columns fall
from 36 to 27, but coefficient evaluations rise from 12 to 13. Median elapsed
times are 0.012 s off versus 0.013 s profiled with full information and 0.011 s
versus 0.012 s with information omitted; REML and ML give the same rounded
medians. This establishes saved derivative work, not a whole-fit speed gain.
No automatic performance policy was changed.

### Batch closure (2026-10-18)

**Status: closed by user request.** The completed work listed below is the
delivered scope of this acceleration batch. Closure does not mean the entire
research roadmap is implemented or that the working-tree changes are committed.

The final isolated macOS build and focused run passed 47 test blocks and
1,282 assertions, with no failures, errors, warnings or skips. Reference checks,
control-name documentation coverage and explanation-only vignette rendering
passed. Earlier full-suite and performance results retain their dated scope;
the full suite was not rerun after the final extensions.

No further implementation, platform validation or performance calibration is
scheduled in this batch. Remaining algorithmic items and release checks below
are deferred to separately approved work. Automatic admission and experimental
defaults remain unchanged; compression and new output/reuse options remain
opt-in. No additional speedup is claimed for the final extensions.

### Completed acceleration work

- [x] Reliable scalar, packed and block PCG solves, true-residual checks,
  residual refresh, model-aware block-Jacobi and operator Nystrom
  preconditioning.
- [x] Compatible warm starts and bounded solution-based coarse corrections.
- [x] Symmetric diagonal-congruence SLQ with exact determinant correction,
  trace control variates, and estimator diagnostics.
- [x] Initial adaptive accuracy calibration followed by frozen probe/depth
  budgets. Continuous iteration-by-iteration accuracy retuning is not implemented.
- [x] True within-fit Krylov recycling: bounded CG-direction capture,
  generalized Ritz selection, current-operator coarse factors, and balanced
  SPD preconditioning. Only accepted-point state is promoted. This is PCG-only;
  direct factorizations have no Krylov subspace to recycle.
- [x] Explicit cross-fit/CV Krylov transfer through `pcg$exportRecycle` and
  `pcg$recycleState`. Exported `fit$pcgRecycleState` is a versioned serializable
  basis plus coefficient-coordinate metadata. Imports reject incompatible
  coordinates, dimensions, rank or non-finite entries; current-operator Ritz
  selection and coarse factors are rebuilt. No implicit global cache, automatic
  coefficient remapping, or default behavior change is introduced.
- [x] Exact sufficient-statistic grouping for eligible Gaussian Henderson
  models, preserving weighted response moments, counts, determinant
  normalization, derivative contractions, and original observation outputs.
- [x] Residual-class cross-product caching and exact derivative trace operands
  for LDLT, CHOLMOD, PCG and automatic selection, including supported REML/ML
  paths. Shared `mmes_control(compression="off")` remains the default.
  The legacy `pcg$compression` alias remains PCG-only.
- [x] Updated `mmes` and `mmes_control` help and added the explanation-only
  `sommer.mmes.options` vignette.
- [x] Explicit `covarianceInformation="none"` opt-out, with `"full"` retained
  as the default. Skip final information inversion and native-parameter
  uncertainty propagation; return NULL information/uncertainty and explicit NA
  uncertainty columns. Preserve optimization AI, exact certification, coefficient
  prediction and PEV. Inference requiring omitted information errors with a refit
  instruction; it never recomputes silently. Latent RR and factor-score
  augmentation remain excluded, along with PQL, solve-only and direct inversion.
- [x] Opt-in C-free saddle-point output through
  `retainCoefficientMatrix=FALSE`, requiring `computeCi=0`. Skip the final
  compatibility assembly and return `C=NULL`. Preserve ordinary coefficients,
  fitted values, residuals and supported covariance reporting. C-dependent
  methods raise explicit refit errors in R and native post-fit entry points;
  no reconstruction occurs. The default still retains C.

Earlier completed work includes sparse factor kernels, banded stationary-AR
precision derivatives, exact trace panels, guarded structural paths, hybrid AI
with exact final information, and the opt-in exact saddle-point prototype.
Their scope and limitations remain documented in the dated sections below.

### Measured results from the latest batch

These are workload-specific observations, not expected gains for every model.

| Workload/change | Before | After | Observation |
| --- | ---: | ---: | --- |
| Replicated PCG fit, exact compression | 0.373 s | 0.284 s | 23.9% lower median elapsed time |
| 36,000 observations / 120 groups, LDLT | 0.424 s | 0.283 s | 33.3% lower median elapsed time |
| Same replicated workload, CHOLMOD | 0.423 s | 0.283 s | 33.1% lower median elapsed time |
| Related SPD systems, learned Krylov recycling | 7,169 iterations | 4,917 iterations | Operator applications including setup: 7,265 to 5,245 |

Direct-solver timings use three repetitions; all compared fits converged, with
maximum relative compared-output error about `2.1e-10`. Recycling can cost more
than it saves on small fits; the combined small REML variant was slightly slower
than compression alone.

Focused validation passed 2,580 assertions at one and four OpenMP workers.
Uncompressed parity covered 162 fits and 16 fields, all bit-identical.
The full suite passed 8,518 assertions and retained one rotation assertion
failure reproduced identically against the preserved baseline. These checks
do not establish Windows or CUDA hardware coverage.

### Deferred implementation backlog (outside the closed batch)

The user requested continuation on 2026-10-18, excluding platform and performance
validation while retaining focused build and regression checks, and subsequently
closed the batch after the delivered extensions. The unchecked items below
remain unimplemented and require a separate approved scope; they are not active
tasks in this batch.

- [ ] Exact structural planning: cost complete assembly/factor/trace/solve
  plans, then add eligible product-grid elimination and multiple low-rank axes.
  Preserve generic fallback and explicit solver selection.
- [x] Reduced tied-coordinate AI assembly for admitted AI-first and EM/PXEM
  direction-blended Henderson fits, preserving exact final information and
  existing `vcc` semantics. Legacy explicit-zero AI and other excluded paths
  remain unreduced.
- [x] Opt-in interior common-scale profiling with fixed/free and constraint
  admission, analytic REML/ML scale, reduced profile AI, explicit variance-floor
  behavior and complete final uncertainty. Boundary-constrained profiling and
  profile-optimizer EM recovery remain deferred.
- [ ] Further saddle-point work: specialized symmetric-indefinite or
  matrix-free indefinite solves. The existing
  sparse-LU prototype is complete; these are new solver/object contracts, not
  missing pieces of that prototype.
- [ ] Update admission policies only when comparative evidence supports them.
  Do not auto-enable measured losers or invent calibrated costs. Benchmark-based
  automatic dispatch cannot be completed while its calibration is excluded.

Broader estimator changes (for example full Hutch++ or continuous adaptive
accuracy) remain separate experiments, not claims about the implemented
control-variate estimator or initial calibration.

### Deferred validation and release work

- [ ] Windows compilation and platform-specific checks.
- [ ] CUDA hardware execution and performance checks.
- [ ] Broader crossover calibration and repeated held-out whole-fit measurements
  before changing automatic dispatch or experimental defaults.

Focused build and regression checks remain required for new runtime changes.
The cross-fit transfer implementation builds in the isolated macOS candidate;
all 13 test blocks in `test-pcg-crossfit.R`, `test-pcg-controls.R`, and
`test-pcg-compression.R` pass. Rd checks and targeted explanation-only vignette
rendering pass. No cross-fit performance benchmark, Windows check, or CUDA
hardware check was run. No new speedup is claimed.
After adding the information opt-out, 34 focused test blocks pass with 872
assertions, zero failures/errors/warnings/skips, covering state transfer,
compression, control merging, hybrid AI, Wald/KR and prediction. Omission-mode
tests preserve fit trajectories and coefficient uncertainty bit-for-bit across
supported solver/likelihood and recovery/blending cases. This is focused
regression coverage, not a new full-suite or performance claim.
The final focused run after C-free output passes 47 test blocks and 1,282
assertions, with zero failures, errors, warnings or skips. It includes native
saddle algebra, C-free assembly counters, serialized-object guards, compatible
covariance inference, state transfer, compression, prediction, hybrid AI and
Wald/KR regression checks. The package builds in the isolated macOS library.
Rd checks, control-name coverage and explanation-only vignette rendering pass.

## Current Baseline

### Covariance Engine Ownership (2026-10-14)

`ai_mme_sp2()` now owns fit orchestration, residual assembly, linear-system
backends, trace preparation, constraints and output. `mme_covariance.h` owns
factor evaluation, precision state construction, working-coordinate derivatives
and tensor operations. `mme_factors.h` owns the representation adapter and
precision-state type; `mme_paths.h` owns pure family, path and workspace policy.
The covariance implementation is included after the native numerical primitives
in `MNR.cpp`. Build dependencies cover all three headers on Unix and Windows.

Model names are reporting metadata, not optimized-path selectors. A factor may
declare these representation capabilities, including with an R evaluator:

- `form="modified_cholesky"`: `innovation` contains one-based `rows`,
	`columns`, and `coefficient_count`. Parameters are unconstrained innovation
	coefficients followed by log innovation variances for coordinates 2 through
	q; the first innovation variance is one. Precision is `T' D^-1 T`, with
	unit-lower-triangular `T`. Exact sparse derivatives use the same row supports.
- `kind="lowrank_diagonal"`: `loadings` contains `rank`, `count`, one-based
	`rows`/`columns`, logical `diagonal`, and `heterogeneous_specific`. Diagonal
	loadings use log coordinates; other loadings use identity coordinates.
	Specific variances and normalization retain the built-in FA/RR convention.
- `form="normalized_normal_equations"` or `"normalized_affine_precision"`:
	`spatial` contains `weights`, finite ordered `lower`/`upper` bounds, and
	positive `degree` for the affine form. Dependence uses bounded-logit
	coordinates. Precision is normalized from `(I-rho W)'(I-rho W)` or
	`diag(degree)-rho W`; covariance is normalized by its first diagonal entry.
- `form="stationary_ar"`: representation `order` declares the PACF order.
	Initial stationary covariance and heterogeneous scaling retain their exact
	boundary treatment. Banded covariance alone never selects this kernel.

Unknown representations keep the generic dense evaluator/derivative path.
Family contracts are exact algebraic assertions, not permission to approximate
an unrelated covariance. R-side validation rejects malformed innovation,
loading and spatial payloads. Existing model-specific precision metadata is
retained for compatibility but no longer controls the solver.

Correctness fixes: lazy covariance derivatives resolve the owning factor before
bounds checking; analytic spatial log-determinant derivatives include the
Kronecker multiplicity; spatial cache admission checks aggregate retained
entries. Unreachable legacy covariance branches and unused eager precision
builder modes are removed. Sparse innovation derivatives remain sparse during
sensitivity applications. Existing experiment defaults, exact trace support,
solver semantics, memory policy, optimizer and fitted-object fields remain.

New regression coverage includes complete innovation residual blocks, anonymous
model labels, a custom innovation covariance callback without ANTE metadata,
malformed capabilities and multi-axis spatial fits against generic dense
references. Saved built-in factors acquire missing capabilities at the R
compilation boundary, preserving descriptors created before this refactor.
Numerical and performance acceptance must use the current isolated
build; historical benchmark claims below are not new refactor measurements.

The current-source isolated build passes 186 installed-package test blocks with
zero failures/errors and six existing skips (three unavailable BLAS-controller
tests and three on-CRAN small-sample inference tests). Log:
`/tmp/sommer-refactor-full-tests.log`. Those skipped paths and Windows compilation
remain unverified.

Paired refactor validation uses the preserved `benchmarks/results/sommer`
installation as baseline and `/tmp/sommer-refactor-20261014/sommer` as candidate.
Each fit uses two iterations, four requested OpenMP threads, single-thread BLAS
environment settings, one warmup and three alternating paired repetitions in
fresh processes. The table reports medians, not fully converged fit timings.

| Fixture | Baseline | Refactor | Peak RSS baseline -> refactor, KiB |
| --- | ---: | ---: | ---: |
| AR2, 1000 levels, 6 IDs, exact panels | 5.099 s | 4.139 s | 636900 -> 613948 |
| ANTE(2), 64 levels, 8 IDs | 1.549 s | 0.593 s | 268724 -> 267008 |
| SAR, 1000 levels, 6 IDs | 2.976 s | 1.845 s | 485072 -> 451844 |
| CAR, 1000 levels, 6 IDs | 2.848 s | 2.030 s | 526788 -> 468120 |
| Dense Gu AR1 REML, 5 levels, 120 IDs | 1.091 s | 0.152 s | 272056 -> 237668 |
| Dense Gu AR1 ML, 5 levels, 120 IDs | 1.074 s | 0.161 s | 273028 -> 240592 |

Maximum relative difference across likelihood trajectories, monitor, covariance,
parameters, coefficients, scores, information and uncertainty is 3.39e-14.
AR2 retains 1994 exact fallback columns; dense REML selects Schur in both builds.
Peak RSS includes loading, warmup and the timed fit. These measurements are
fixture-specific and do not establish a universal refactor speedup. BLAS thread
environment settings are requests, not a measured utilization assertion.
Artifact: `benchmarks/results/refactor-validation-20261014.rds`, including paired
rows, output errors, source/binary hashes and process-log paths. Timing log:
`/tmp/sommer-refactor-benchmark.log`.

Already implemented and therefore not proposed again: matrix-free PCG panel solves; D-first ML on eligible direct fits; AR1 derivative dispatch; AR(p)/ANTE innovation-row precision assembly; normalized SAR/CAR precision and derivatives; factor/grid evaluation caches; bounded selected-inverse dependency caches; guarded CHOLMOD mode-1 extraction; aggregate inverse-cache budgets; same-probe PCG residual traces without `Btrace`; support-aware tensor contractions; sufficient PD certificates; and PCG probe-error/true-residual diagnostics.

Important limits: sparse AR(p)/ANTE precision operands and sparse SAR/CAR evaluation are opt-in, but dense compatibility storage remains for some consumers; trace planning is not fully demand-driven; default PCG estimator and tolerance settings remain unchanged, with additional opt-in choices described in the current status above; product-grid elimination and multiple low-rank axes are unsupported. The exact block-banded prototype is opt-in (force only) because it was up to 80x slower than generic CHOLMOD on AR1/AR2 fixtures.

## Progress Checklist (2026-10-13)

### Calibration Restart (2026-10-14)

### Sparse Factor Kernels (2026-10-14)

- [x] Add opt-in `mmes_control(experimental = list(sparseFactors = "force"))`. Default and `"auto"` remain off. Native AR1/AR(p)/ANTE innovation precisions provide exact sparse operands for derivative sandwich products and tensor-mode applications, including compact residual blocks and matrix-free prior application. ANTE precision construction uses compact innovation rows instead of a dense transition temporary. ANTE covariance evaluation uses the exact row recurrence in both fitting engines, avoiding triangular inversion and cubic covariance multiplication. Public dense covariance and existing derivative contracts are retained.
- [x] Add exact sparse SAR/CAR covariance and working-coordinate derivative evaluation. Assemble `Q = (I-rho W)'(I-rho W)` for SAR or `Q = D-rho W` for CAR, factor by sparse LDLT, and solve identity/derivative RHS panels of at most 32 columns. Normalize by `Q^{-1}[1,1]`. Obtain determinant and normalization derivatives from the same factor and sparse precision products. Preserve signed/asymmetric SAR weights, CAR dependence intervals, normalization and the existing downstream PD checks. No derivative support threshold, probe-policy change or optimizer change is introduced.
- [x] Isolated paired whole-fit benchmark: ten fixtures, AR2/AR3 at 200/1000 levels, free ANTE(2) at 8/16 levels with 32 IDs, SAR/CAR at 200/1000 levels. Other families use six IDs. Four OpenMP threads, single-thread BLAS, REML, LDLT, Schur off, trace batching force on in both modes, four iterations, one warmup and three timed repetitions per mode in fresh processes, alternating mode order. Maximum relative difference across likelihood trajectories, monitor, covariance matrices/parameters, BLUE/BLUP, scores, AI and uncertainty is 5.23e-9. These are four-iteration whole-fit timings, not claims about fully converged large-fit runtime.

| Fixture | Dense-path median | Sparse-kernel median | Ratio | Median Peak RSS, KiB |
| --- | ---: | ---: | ---: | ---: |
| AR2, 1000 levels | 6.608 s | 6.577 s | 1.005x | 605828 -> 572352 |
| AR3, 1000 levels | 8.037 s | 7.835 s | 1.026x | 746532 -> 763880 |
| ANTE(2), 16 levels | 1.383 s | 1.383 s | 1.000x | 256308 -> 256148 |
| SAR, 1000 levels | 6.287 s | 4.770 s | 1.318x | 600172 -> 482536 |
| CAR, 1000 levels | 5.699 s | 4.168 s | 1.367x | 576088 -> 529968 |

SAR peak RSS falls about 114.9 MiB and CAR about 45.0 MiB. Peak process RSS includes package loading, warmup and the timed fit; it is not a pure kernel allocation measure. AR2/AR3 and small ANTE runtime gains are not established beyond noise. SAR at 200 levels is slightly slower (1.570 s off, 1.606 s force), so no automatic sparse/dense crossover rule is claimed.

- [x] Focused exact-path tests cover AR1-3, high fixed heterogeneous AR2 PACFs, ANTE orders 1/3, sparse and dense adjacency, signed/asymmetric/negative SAR and zero/high fixed CAR dependence, REML/ML, LDLT/CHOLMOD, weighted data, missing observations, ragged correlated residuals, `computeCi=0,1,2`, PEV and normalized full-inverse residuals. Direct-engine fits are compared in mode 0; this fixture already fails in the old direct mode-2 inverse path, which remains unchanged and unverified. Same-probe PCG checks exercise exact operator application only; PCG estimators and tolerances stay unchanged. Identifiable full-rank simulated covariance effects avoid misleading ANTE information-matrix comparisons at a nearly singular low-rank fixture.

Artifacts: `benchmarks/validate_factor_kernels.R`, `benchmarks/results/factor-kernels-validation.rds`, its CSV/log and per-process `.runs` directory; spatial crossover harness/results: `benchmarks/validate_spatial_crossover.R`, `spatial-crossover-validation.rds` and `spatial-crossover-large.rds` with CSV/logs and `.runs` directories. Candidate installed as `benchmarks/results/sommer`; previous trace-panel comparison library is untouched.

### Remaining Rank 1: Storage, Shared Spatial Factor, And Crossover (2026-10-14)

- [x] Canonical precision storage: enabled AR(p)/ANTE retain banded precision once as sparse factors. The full Kronecker precision state is no longer built; lambda materializes only the matrix the solver needs. Matrix-free PCG reads `totalDim` and exact factor-wise diagonals; tensor products keep sparse axes sparse; dense materialization occurs only at compatibility consumers. Focused regression covers LDLT/CHOLMOD/direct mode 0/PCG, residual and uncertainty paths.
- [x] ANTE derivative storage: beta/innovation precision derivatives use exact sparse row-support matrices; covariance derivatives recompute on demand, not one dense q-by-q matrix per parameter. AI forms same-factor ANTE derivative products sparsely. Order-2 ANTE with 8 IDs: q=64 whole-fit median 2.205 s -> 1.631 s, peak RSS 285300 -> 268964 KiB. q=128: 21.793 s -> 9.337 s, peak RSS 451072 -> 299796 KiB (-147.7 MiB). Paired outputs match within 3.60e-9 across all kernel fixtures. These 2-iteration q>=64 runs use 8 IDs; smaller fixtures use 32 IDs.
- [x] Reuse normalized SAR/CAR sparse LDLT across covariance, precision and derivative consumers. `sparseSpatialFactorizationReuses` records six hits per fitted spatial term in the crossover fixtures. Repeated covariance, logdet, precision and derivative outputs pass parity.
- [x] Crossover measurement: 18 paired fixtures, three repetitions, SAR directed and CAR symmetric weights, q={60,150,300}, requested off-diagonal densities {0.05,0.25,1.0}, REML/LDLT, four threads, three iterations. Maximum relative fit-output error is 2.30e-14. No stable density-only cutoff emerged. At q=300 and density 1, sparse mode is 29% slower for SAR and 11% slower for CAR. At q=60/150, differences are mostly within a few percent and change sign across family/density. Peak RSS is lower by about 5-27 MiB, with six factor-cache reuses per sparse fit.
- [x] Extended fully dense q=1000 checks (three pairs/family): SAR off 35.4-37.4 s versus force 58.9-61.7 s; CAR off 32.6-34.6 s versus force 44.5-47.6 s. Sparse mode raises peak RSS from 596316 to 685904 KiB for SAR and 624180 to 675356 KiB for CAR. Maximum relative fit-output error is 4.50e-14. Thus keep sparseFactors opt-in; no automatic crossover rule is justified. This dense 1000-level SAR/CAR fit uses six IDs and four-thread runs. The host was loaded, so use timings as fixture-specific evidence, not universal constants.

Conclusion: sparse adjacency can pay off at large q, but full/dense adjacency reverses that result. Keep `sparseFactors="force"` opt-in. Remaining limits: direct-engine mode 2 for this fixture is a pre-existing failure; public dense covariance/derivative outputs and SAR/CAR constructor spectral bounds remain dense; no full suite rerun after the final state-storage edits.

Final current-source build passes. Focused factor suite: nine test blocks, zero failures/errors. The prior full-suite run passed 182 tests, zero failures/errors, six existing skips; full suite has not been rerun after the last state-storage/tensor edits. Logs: `benchmarks/results/factor-kernels-tests.log` and `factor-kernels-full-tests.log`.

Remaining: public dense covariance/derivative outputs and SAR/CAR constructor spectral bounds remain dense; direct-engine `computeCi=2` remains unverified because the old path already fails on the tested fixture. Dense adjacency crossover shows no automatic sparse admission opportunity at tested sizes. Further estimator and optimizer changes remain later due greater statistical risk.

### Exact AR2 Trace Panels (2026-10-14)

- [x] Add opt-in `mmes_control(experimental = list(traceBatch = "force"))` for LDLT random-score trace misses. Collect only columns with incomplete selected-inverse coverage, pack at most 32 RHS columns per panel (a wider relationship block remains indivisible), and run independent panels in OpenMP workers. The workers call the immutable Eigen LDLT factor directly, not the stateful generic solver dispatch. Each worker owns its RHS and solution; no CHOLMOD, PCG, R callback, or shared solve workspace runs concurrently. Non-finite solutions are reported from the main thread. All derivative support and exact solves are retained; no thresholds, estimator changes, or optimizer changes are introduced.
- [x] Preserve the old path by default (`traceBatch = "off"`; `"auto"` also leaves it disabled). Conservative panel scratch accounting covers three dense panels per worker plus retained selected-inverse dependency tables. `memoryPolicy = "speed"` stops with the required budget; `"limit"` uses the old exact fallback when the panel budget is insufficient. This is not a process-RSS cap. `traceFallbackSolves` retains its meaning as missed covariance-block columns; new `traceBatchSolves` counts actual multi-RHS calls and `traceBatchWorkspaceBytes` records conservative panel scratch.
- [x] Target benchmark: 1000 AR2 levels, six IDs, two replicates, REML, LDLT, Schur off, four iterations, one warmup and three timed repetitions per mode/thread count, alternating mode order in fresh single-BLAS-thread processes. At four OpenMP threads, median whole-fit runtime falls from 7.959 s to 6.084 s (1.31x); score traces fall from 3.003 s to 1.237 s (2.43x). The 3988 missed block columns across four iterations become 800 independent panel calls. Median peak process RSS rises from 596540 KiB to 611116 KiB (about 14.2 MiB, 2.4%). RSS includes package loading, warmup, and the timed whole fit. One thread does not benefit: 10.430 s off versus 10.722 s force, so automatic admission remains disabled.
- [x] Repeated target fits match likelihood trajectories, monitor, covariance matrices/parameters, BLUE/BLUP, scores, AI and parameter uncertainty bit-for-bit within each thread count. Comparison against the pre-change installed LDLT path also passes; cross-thread differences are below 1e-8. Focused regression coverage includes REML/ML, one/four threads, sparse/dense/no relationship, missing records, partial panels, `computeCi=0,1,2`, PEV, normalized true full-inverse residuals, longer free-parameter fits, zero/high fixed PACFs with free heterogeneous variances, and both memory policies. High PACFs are fixed in the heterogeneous fixture because the original free-PACF proposal can hit a singular boundary; that existing proposal behavior is not changed here.

Rejected exact prototypes: serial ordinary multi-RHS batching was slower; extending off-pattern Takahashi recurrence dependencies under the existing AMD ordering created 8.6 million extra entries and slowed the target fit to 44.4 s; a zero-valued symbolic demand envelope exceeded the 1000 MB workspace budget; custom row-major triangular traversal was slower than the existing Eigen kernel. None of these prototypes remains in production code. The rejected derivative-threshold shortcut remains excluded.

Artifacts: `benchmarks/validate_ar2_trace.R`, `benchmarks/results/ar2-trace-validation.rds`, its CSV/log and per-process `.runs` directory. Candidate binary: `benchmarks/results/lib-fill-candidate-v2-20261014` (rebuilt for this change; historical fill-validation RDS files remain unchanged).

Final installed-package suite: 181 tests, zero failures/errors, six existing skips (three BLAS-controller and three on-CRAN tests). Logs: `benchmarks/results/ar2-trace-tests.log` and `ar2-trace-full-tests.log`. Those skipped paths remain unverified.

The subsequent sparse-factor-kernel stage is recorded above. PCG estimator and optimizer changes remain later. AR2 AI assembly remains a separate measured cost; these improvements do not claim to remove all AR2 fit costs.

### Validated Dense-Schur Guard (2026-10-14)

- [x] Replace the general dense-Schur density admission scan with `BlockSchurPolicy::fillFavoursSchur` in `src/mme_paths.h`. Coupled plans require finite, nonnegative AMD-predicted sparse-factor work of at least 0.05 times the finite positive dense-plan cost. This is a conservative cost-ratio guard, not a 5-percent fill fraction or a universal crossover estimate. Invalid or incomplete sparse estimates no longer automatically admit Schur. `force` still bypasses the performance gate; memory and structural eligibility remain enforced.
- [x] Preserve multiple uncoupled random groups (`couplings.empty()` and more than one candidate group). Isolated validation showed that applying the factorization-only ratio to these plans regressed iid AR1 from 0.224 s to 0.286 s. The final rule restores 0.224 s because their solve/trace benefit is not captured by factorization flops alone. Chain, latent and banded admission are unchanged by this guard; their remaining density tests are not replaced by an unvalidated cost model.
- [x] Isolated before/after validation in `benchmarks/validate_fill_dispatch.R`: nine fixtures covering connected random low-fill and banded graphs, sparse and dense relationships, iid AR1, AR2, SAR and CAR; one thread, warmups, alternating auto/off order, three timed repetitions and median elapsed time. The low-fill random graph improves from 2.968 s to 0.641 s (4.63x). Other measured baseline/candidate differences stay within about 4 percent, including the restored iid fixture. These are fixture-specific results, not an expected speedup for all models.
- [x] Numerical parity across baseline, auto and off paths covers likelihood trajectories, covariance parameters, BLUE/BLUP, information matrices, score and uncertainty; additional small fixtures exercise `computeCi=0,1,2`. Maximum baseline/final-candidate relative output difference is 1.60e-9. Native policy tests pass. Full suite before the new fixture: 179 tests, zero failures/errors, six existing skips (three BLAS-controller and three on-CRAN tests). Updated path file including the new regression fixture: eight tests, zero failures/errors. Those skipped paths remain unverified.

Artifacts: `benchmarks/results/fill-validation-comparison-v2.csv`, `fill-validation-baseline.rds`, `fill-validation-candidate-v2.rds`, `fill-validation-tests.log`, `fill-validation-path-tests.log`. Final isolated library: `benchmarks/results/lib-fill-candidate-v2-20261014`. Screening was stopped after 322 of 360 cells at the user's request to use current results; all 120 design geometries are represented, but 38 seed replications remain pending. Broader chain/latent cost pricing, a globally calibrated predictor, ML and weighted/ragged/multiple-term workloads are not marked complete.

Factorial analysis of the 322-cell snapshot: `benchmarks/analyze_fill_factorial.R` models paired log-speedup with graph seed as a block, five categorical design factors, and all two-factor interactions. All 120 design geometries are represented; 38 graph-seed replications remain missing. The model has full rank (41/41). Main-effects R-squared is 0.293; adding interactions raises it to 0.684 (adjusted 0.639; nested partial F = 11.976). The largest conditional interaction sums of squares are block size x input nonzero fraction, input fraction x border, block size x group count, group count x input fraction, and block size x graph pattern. Input nonzero fraction is an experimental knob, not proposed production dispatch logic. Exact matched graph levels collapse where connectivity imposes a floor, so nominal fraction effects must be interpreted with that constraint.

Seed-303 validation has 108 cells. Log-speedup RMSE: fill-only 0.725, main-effects 0.605, two-factor 0.430. A predicted 1.2x admission threshold under the two-factor model admits 40 cells, with one actual ratio 0.897 and two missed wins above 1.2x. This is substantially better than fill-only (four admitted losses above 10%, worst ratio 0.463), but not safe enough to deploy directly. Primary analysis includes inactive force attempts as no-op controls; actual engine subtype is a treatment outcome, not a randomized factor. These results describe the block-engine family; they cannot be relabeled as a dense-Schur-only rule. Timing p-values are conditional on this hardware, incomplete replication and shared-worker noise, not universal guarantees. The next rule should combine valid fill/cost features with block size, border and group count and must pass isolated validation. No production coefficient changes arise from this factorial analysis.

Artifacts: `benchmarks/results/fill-factorial.txt`, `fill-factorial-marginals.csv`, `fill-factorial-validation.csv`, `fill-factorial-heldout.csv`, and `fill-factorial-interactions.png`.

Current 322-cell snapshot visualization: `benchmarks/plot_fill_speed.R` produces `benchmarks/results/fill-speed.png`, its PDF version, and `fill-speed-paired.png`. The horizontal axis is AMD-predicted factor entries `nnz(L)`, not newly introduced fill entries. Speed is inverse whole-fit elapsed seconds. Engine subtypes are separate: 78 dense-Schur contrasts and 232 chain contrasts pass eligibility/parity gates. A descriptive paired log-speedup versus log-factor-size regression has R-squared 0.097 for dense Schur, with equality outside the observed range and no dense-Schur speedups above 1.2x; chain has R-squared 0.260 and fitted equality near 30,395 factor entries, but wide variation at similar fill. These are pooled descriptive fits, not a validated global cutoff; matrix size, border and engine subtype remain important. The conservative dense fill-cost guard passed the isolated validation recorded above; it is not established as an optimal global rule.

The reconnect moved the session to another host; `/tmp/sommer_state` did not survive.
The revised experiment uses workspace-persistent paths under `benchmarks/results/`.

- Complete medium design: group size {16, 32, 64, 128, 256}, group count {10, 50}, input off-diagonal fraction {0.02, 0.1, 0.5}, fixed covariate count {5, 50}, graph pattern {random, banded}, and seeds {101, 202, 303}. This is 360 cells with no omitted size combinations.
- Both patterns use the same edge count at each size/fraction. A chain backbone guarantees connectivity. The minimum connected edge count can exceed the requested fraction for small graphs; `actual_input_fraction` records this floor. Some low-fraction small cells therefore collapse to the same graph. These are not independent observations of a fraction effect.
- Graph generation and response/covariate generation use separate seeded streams, so paired patterns receive identical data. Positive definiteness, matched edge counts, factorial cardinality and fixed-design dimension feasibility are checked by `--audit` before timing. The audited scope is AR1 coupling with two replicates, REML, and `computeCi=0`; ML, ragged/weighted residuals, multiple random terms and other graph families still need external validation.
- Four independent single-threaded workers, randomized cell scheduling, alternating on/off order within repetitions, one warmup per engine, and three timed repetitions summarized by the median. There are 2160 timed fits plus 720 warmups. These are parallel throughput measurements, not isolated latency estimates; candidate rules require an isolated timing confirmation.
- Every timed pair must match likelihood trajectory, covariance parameters, coefficients, information matrix and coefficient evaluation count within the stated numerical gate. Record actual engine subtype and reject inactive contrasts. Chain costs are not replaced by dense-Schur cost estimates.
- `benchmarks/fit_fill_rule.R` separates engine subtypes and reserves graph seed 303 for validation. Invalid/incomplete fill estimates cannot establish a safe dispatch decision. No new production threshold is justified until validation and isolated confirmation pass.
- Checkpoints: `benchmarks/results/fill-20261014-v2.csv.cells/cell-*.rds`, written atomically. The manifest locks design, repetitions, summary statistic, worker count, benchmark script hash, loaded binary hash and thread settings. Repeating the same command resumes only missing cells. The installed benchmark library is `benchmarks/results/lib-fill-20261014`.

Resume command from the package root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
R_LIBS="$PWD/benchmarks/results/lib-fill-20261014" \
Rscript benchmarks/calibrate_paths_design.R --medium --workers=4 --reps=3 \
	--out=benchmarks/results/fill-20261014-v2.csv
```

The previous fragmentation rule was removed after a true on/off comparison contradicted the memory-budget proxy comparison. Its historical checklist entry below is superseded; do not reintroduce it from those timings.

Done:

- [x] Rank 0 profiling: `mmes_control(experimental = list(profile = "force"))` fills `engineDiagnostics$profile` (13 lap phases, per-iteration seconds, counters); `profile_mme_workload()`, `phase_*` columns in `run_mme_benchmark()`, `benchmarks/profile_phase0.R`, tests in `test-mme-profile.R`. Peak memory is only the existing cache-budget accounting.
- [x] `memoryPolicy = "speed" | "limit"` (default speed): stop with the memory needed instead of a slower fallback. Enforced only where the fallback was measured slower (CHOLMOD selected-inverse caches: 10 MB budget was 6x slower than 100 MB on dense Gu). NOT enforced for the LDLT dependency table (a 12 GB table was 3x slower than none) or for Schur/chain/banded admission (their fallbacks are not slower).
- [x] Block-banded engine is force-only (auto admission removed).
- [x] General dense-Schur admission now uses the validated fill-cost guard above instead of block density. CHOLMOD dense inverse-block cache still skips sparse groups (sparse pedigree 20 levels x 150 IDs: 3.3 s to 0.29 s).
- [x] LDLT selected inverse indexes contiguous (dense) columns by subtraction (dense Gu, 300 IDs, `ldlt`: 72 s to 16.5 s).

- [x] Investigated and removed the provisional `tooFragmented` rule: a true on/off comparison contradicted the earlier memory-budget proxy. Do not reintroduce its 64/1000 cutoffs from those timings.

- [x] PD check for single-factor descriptors with 128 or more effects uses a shifted Cholesky test instead of a dense eigendecomposition (`descriptorIsPositiveDefinite`). SAR/CAR 1000 levels, `ldlt`, 6 iterations: `constraints` 1.91 s to 0.78 s (SAR) and 1.48 s to 0.35 s (CAR); log-likelihood and estimates identical to 10 digits.

- [x] R-side starting values no longer run a dense `qr()` of the fixed design (`.mmes_start_lsq`, sparse crossprod + Cholesky for well-conditioned full-rank designs, dense QR otherwise). AR2 with 1000 levels (X is 12000 x 1000): elapsed 56.7 s to 9.4 s; log-likelihood and estimates unchanged (llik diff 3e-8, max relative theta diff 4e-7). Found by gdb sampling: the engine profile only covers the C++ fit and had hidden that about 90% of wall time was `dqrdc2` with threaded level-1 BLAS.

- [x] `experimental$blockSchur = "auto" | "off" | "force"`: off disables the dense Schur, chain, banded and latent engines; force skips the density and group-size admission rules (memory limits still apply). Needed for the crossed-design calibration (`benchmarks/calibrate_paths_design.R`). New profile counter `traceFallbackSolves`.

Tried and reverted: skipping the supernodal selected inverse for groups skipped only because they are sparse. It made sparse-pedigree fits much slower (`ar1_20L_150ID_sparse` 0.29 s to 4.9 s), so the selected inverse is what makes those fast.

Open, in order (each starts from the current working tree):

- [ ] 1b. `ar1_50L_100ID` (no relationship, Schur off): trace preparation costs 0.85 s with a large budget versus 0.06 s when the selected inverse is not built. The right choice between selected inverse and exact sparse solves depends on fill and requested support; needs the cost-model/autotune design below, not another threshold. `car_200L` also still selects Schur (0.82 s versus 0.60 s).
- [x] 2. SAR/CAR exact sparse evaluation implemented opt-in above. At 1000 levels, proposal time drops from 0.426 s to 0.154 s for SAR and 0.433 s to 0.099 s for CAR in the paired four-iteration fixtures. Dense covariance outputs, dense constructor spectral bounds and automatic crossover calibration remain open; the PD-repair part was completed earlier.
- [x] 3. AR2 exact trace-solve bottleneck: opt-in independent LDLT panels validated above, without the shared solve lock. The original 997 missed block columns per iteration remain exact, now handled in 200 panel calls at six IDs. Single-thread automatic acceleration is not established. `aiAssembly` and dense `fit$theta` storage remain unchanged.
- [x] Rank 3: PCG controls, reliability, estimators, initial calibration and
  within-/cross-fit recycling are implemented; see current status.
- [x] Rank 6 first slices: exact sparse-LU saddle prototype and opt-in C-free
  output are implemented. Other indefinite solvers and automatic admission
  remain open.
- [ ] Rank 4 beyond the admission fixes above remains open.
- [x] Rank 5 first slice: opt-in hybrid AI with exact final information and
  sensitivity-solve counters. Information opt-out is now implemented;
  tied-coordinate reduction remains open.

Measurement notes: forced `cholmod` on SAR/CAR 1000 levels is 2-3x slower than `ldlt` (`solver = "auto"` already picks `ldlt` there). The covariance-engine refactor adds explicit header dependencies and preserves `all: $(SHLIB)` as the default build target; manual source touching is no longer required.

## Ranked Work

| Rank | Improvement | Likely beneficiaries | Possible speed win | Cost / risk |
| --- | --- | --- | --- | --- |
| 0 | Profile representative fits by phase (DONE) | All workloads | No immediate speedup; prevents optimizing non-dominant work | Low |
| 1 | Sparse/banded factor representations and derivatives | Large AR(p)/ANTE factors; selected spatial families | High within factor kernels when bandwidth is small; whole-fit gain unknown | Medium-high; exactness and representation changes |
| 2 | Demand-aware direct trace and sensitivity solves | Many free parameters, large/dense relationship terms, selected-inverse misses | High when repeated exact solves dominate; negligible when current caches cover demand | Medium-high; preserve exact fallback |
| 3 | Accuracy-aware PCG log-determinants and traces | Ill-conditioned, PCG-dominated large fits | Medium-to-high if fewer Lanczos steps/probes meet the same accuracy target; unmeasured | High; stochastic estimator and convergence risk |
| 4 | Better exact structural planning | Low-width product grids, multiple components, eligible low-rank structures | Potentially high for matching structure; no benefit or regressions outside admission gates | High; new solver architecture |
| 5 | Reduce repeated AI work | Many free covariance parameters when sensitivities dominate | Medium if sensitivity solves dominate; none when factorization or assembly dominates | High; quasi-Newton changes convergence, information opt-out changes output contract |
| 6 | Saddle-point MME solver | Cases where normal-equation fill or conditioning dominates | Unknown; may reduce products but can increase fill and solve cost | Very high; indefinite solver and uncertainty redesign |

### Rank 0: Profile Before Choosing A Kernel (done)

Implemented; see the checklist. Remaining use: rerun `benchmarks/profile_phase0.R` before and after each later change.

### Rank 1: Sparse/Banded Factor Kernels

Add representation-aware application, storage, and working-coordinate derivatives for AR(p) and ANTE first. Their innovation-row formulas are exact, but dense storage can still cost O(q^2); a bounded-band representation targets O(qb) storage and O(qb^2) factor work for dimension q and bandwidth b. This is a kernel-level complexity target, not an end-to-end speedup estimate.

Then assess direct sparse SAR/CAR evaluation and structured Toeplitz/MA operations independently. Preserve normalization, covariance/log-determinant identities, support at zero-valued parameters, generic dense fallbacks, and the existing PD floor. Do not infer that banded covariance implies banded precision.

**Gate:** dense covariance, precision, determinant, derivative, fit, score, AI, and uncertainty parity; benchmark full fits across factor sizes and residual layouts.

#### Iteration-speed baseline refresh (2026-10-17)

Single-run screening profiles used the current source at `f755907`, four REML
iterations, LDLT, one requested OpenMP thread, and single-thread BLAS settings.
They are not converged-fit timings or speedup claims.

| Workload | Effects | Free parameters | Trace prep | Sensitivities | AI assembly | Score traces | Sensitivity RHS / trace fallback |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| US, 12 levels, 1000 IDs, no relationship | 12012 | 79 | 0.043 s | 0.115 s | 0.347 s | 0.037 s | 316 / 0 |
| AR2, 1000 levels, 6 IDs, no relationship | 7000 | 4 | 4.517 s | 0.060 s | 8.978 s | 3.049 s | 16 / 3988 |

The US fit took 1.512 s overall; the AR2 fit took 18.064 s. On AR2, AI assembly
accounts for about half the profiled iteration work. Trace preparation and score
traces account for much of the rest; sensitivity solves are negligible. In the
79-parameter US fit, AI assembly also exceeds sensitivity time, and score traces
have no selected-inverse fallback. These measurements do not justify treating
fewer sensitivity solves as the universal next optimization.

Next step: instrument the dominant AI-assembly subwork and isolate exact trace
preparation versus fallback costs. Optimize one measured operation at a time,
then compare paired whole-fit timings and numerical outputs before recording a
gain.

#### AI and trace subphase profile (2026-10-17)

Added opt-in timers under `experimental$profile="force"` for trace
derivative/pre-index preparation, selected-inverse construction, random score-trace
fallback solves, and AI assembly components. Timers remain absent when profiling
is off. The EM/PX-GEM wrapper aggregates these counters across its one-iteration
native fits. Focused profile tests pass.

Three repeated runs used the same fixtures, four REML iterations, LDLT, one
requested OpenMP thread, and single-thread BLAS settings. Values below are
medians; these are profiled screening timings, not converged-fit comparisons or
speedup claims.

| Workload | Fit | Trace prep: total / pre-index prep / selected inverse | AI: total / random first term / random curvature | Score traces / fallback solve time / fallback solves |
| --- | ---: | ---: | ---: | ---: |
| US, 12 levels, 1000 IDs | 1.529 s | 0.0449 / 0.0025 / 0.0426 s | 0.3709 / 0.3197 / 0.0500 s | 0.0386 / 0 s / 0 |
| AR2, 1000 levels, 6 IDs | 18.271 s | 4.596 / 4.577 / 0.0189 s | 9.140 / 0.00037 / 9.137 s | 3.075 / 2.829 s / 3988 |

For AR2, exact random-covariance curvature consumes nearly all AI assembly;
score-trace fallbacks consume most score-trace time. Selected-inverse
construction is small. The broad derivative-preparation interval dominates
trace preparation, but still needs a finer split before optimization.
For high-parameter US, the random first-term matrix product dominates AI
assembly, and selected-inverse construction dominates trace preparation.
These workloads have different bottlenecks; do not add a global dispatch rule.

#### Symmetric AI curvature pair reduction (2026-10-17)

The exact random-covariance curvature term is symmetric in its parameter pair.
Compute only its upper triangle, mirror each off-diagonal correction, and retain
existing final random-block symmetrization. The non-symmetric first AI term is
unchanged.

Focused profile, hybrid-AI, and AI-recovery suites pass. Three alternating
fresh-process pairs use `benchmarks/validate_ai_curvature.R`, baseline
`/tmp/sommer-iteration-instrumented-v4`, and candidate
`/tmp/sommer-iteration-candidate`. All fits use four REML iterations, LDLT,
one requested OpenMP thread, one warmup, and single-thread BLAS settings.
Maximum normalized difference across likelihood, monitor, covariance,
coefficients, information, score, and uncertainty is `4.43e-10`.

Median random-curvature time falls from 9.114 s to 6.828 s, about 25%.
However, scheduler noise made candidate whole-fit times 943.1 s, 993.2 s, and
15.9 s; baseline times were 18.5 s, 40.4 s, and 18.4 s. Phase-level outliers
also moved between trace preparation and AI assembly. These are not credible
whole-fit measurements. Keep the exact symmetry reduction, but claim no
whole-fit speed gain until quieter-host validation.

Next: optimize the measured random factor-wise precision derivative
construction, then rerun isolated whole-fit confirmation. Keep workload
specificity and profiling overhead visible.

#### AR2 trace derivative preparation split (2026-10-17)

Added opt-in counters to separate random and residual precision derivatives,
grid and non-grid residual derivative bases, and remaining pre-index work.
The focused profiling test checks counters remain nonnegative and within the
trace-preparation phase.

Three fresh-process runs used the same 1000-level AR2 fixture, four REML
iterations, LDLT, one requested OpenMP thread, and single-thread BLAS:

| Metric | Median | Range |
| --- | ---: | ---: |
| Fit elapsed | 15.957 s | 15.909-16.119 s |
| Trace derivative preparation | 4.544 s | 4.536-4.606 s |
| Random precision derivative construction | 4.544 s | 4.535-4.606 s |
| Residual precision derivatives | 0 s | 0 s |
| Grid / non-grid residual derivative bases | 0 / 0 s | 0 / 0 s |
| Remaining pre-index work | 0 s | 0 s |
| Selected-inverse construction | 0.019 s | 0.019 s |

AR2's pre-index cost is entirely in random factor-wise precision derivative
construction; selected-inverse construction and residual derivative work are
negligible for this fixture. The next optimization should target this exact
operation, then recheck numerical parity and paired whole-fit timings. No
whole-fit speed gain is claimed from this profiling step.

#### Banded stationary-AR precision derivatives (2026-10-17)

The generic precision derivative computed `-Q (dK) Q` with dense matrices for
stationary AR factors, despite their exact order-banded precision structure.
The factor engine now builds each perturbed stationary-AR precision with its
existing innovation representation, differences those matrices using the
existing `1e-6` working-coordinate step, and stores only the exact order-wide
band. Other covariance families retain their existing derivative path.

The first analytic covariance-derivative trial did not reduce the measured
precision-derivative phase, so that trial was removed. The direct banded
precision derivative passed focused AR2/AR3 covariance safeguards and profile
tests. Paired heterogeneous AR2 and AR3 fits also matched likelihood,
parameters, coefficients, score, information, and uncertainty to at least
`1.4e-9` relative difference.

Three alternating fresh-process pairs used `benchmarks/validate_ai_curvature.R`,
the same 1000-level AR2 workload, four REML iterations, LDLT, one requested
OpenMP thread, and single-thread BLAS:

| Metric | Baseline median | Banded derivative median | Median speedup |
| --- | ---: | ---: | ---: |
| Whole-fit elapsed | 15.573 s | 1.802 s | 8.63x |
| AI random curvature | 6.685 s | 0.0865 s | 77.3x |
| Trace derivative preparation | 4.460 s | 0.0930 s | 48.0x |

Candidate random precision-derivative time split into 0.0408 s covariance
derivatives, 0.0444 s precision derivatives, and 0.0077 s other work. Maximum
normalized difference across likelihood, monitor, covariance, coefficients,
information, score, and uncertainty was `4.02e-7`; `InfMat` and `dLu` use a
`1e-6` acceptance threshold, while other outputs retain `1e-8`. Peak RSS rose
from 786 MiB to 795 MiB. None of the four-iteration screening fits converged;
these results establish a substantial speedup for this workload, not a
converged-fit runtime estimate. Smaller heterogeneous AR2/AR3 parity checks
remain below `1.4e-9`.

More exact derivative sparsity can eliminate score-trace fallback solves when
the selected inverse already covers all required entries. The trace-profile
test now checks fallback and batch counts remain consistent without requiring
a positive batch-solve count.

### Rank 2: Demand-Aware Direct Traces And Sensitivities

Choose between selected inverse entries, cached blocks, and exact solves based on requested support, elimination-tree reach, reuse count, and aggregate workspace. Investigate sparse-RHS tree reach for AI sensitivities separately from already implemented matrix-free PCG panel batching. Consider touched-row or chunked CHOLMOD workspaces only where setup and packing costs pay back.

**Potential win:** removes repeated work only for workloads where traces or sensitivity solves dominate and current selected-inverse/cache coverage misses. Worst-case reach can equal the full factor, yielding little or no gain.

**Gate:** compare every result and true residual against current exact solves; benchmark including cache setup, misses, memory, and all sensitivity consumers.

#### Rejected dense selected-inverse solve (2026-10-17)

An exact multi-RHS solve for dense LDLT fill looked faster in a two-iteration
screen, but regressed the converged target. On a dense-relationship US model
with 6 levels and 30 IDs, three paired 15-iteration REML fits converged in both
versions. Median whole-fit time rose from 0.127 s to 0.145 s, and selected
inverse preparation rose from 0.0388 s to 0.0578 s. Maximum normalized output
difference remained below `2.9e-13`. The dense-solve dispatch was removed; the
screening result is not evidence of an accepted speedup.

#### Sparse-US trace and sensitivity demand profile (2026-10-17)

Added opt-in profile counters for selected-inverse lookups and their exact
Takahashi dependency closure, plus sensitivity RHS nonzeros, direct-LDLT
triangular reach, and solve time. Counters are inactive unless profiling is
enabled. Lookup collection is synchronized because score-trace lookups may
run in parallel; demand-key retention is capped at one million entries and
reports when that cap prevents a complete closure measurement.

Measured a six-level, 1,000-ID US model with sparse tridiagonal relationship
precision, two replicates per ID, one requested OpenMP thread, and single-thread
BLAS. LDLT REML and ML fits both converged in six iterations with 22 covariance
parameters. The initial profiled runs took 0.524 s REML and 0.463 s ML.
Profiling phases put score traces at 28-28 ms, selected-inverse preparation at
17-32 ms, sensitivities at 25-26 ms, and AI assembly at 18-19 ms. Initial setup
was about 302-307 ms, so this trace/sensitivity target is not the dominant
whole-fit cost on this fixture.

Exact demand counters show why the proposed sparse work reduction is not
promising here:

| Mode | Unique demanded inverse entries | Exact demand closure | Filled pattern entries | Closure operations | Sensitivity RHS density | Exact solve reach |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| REML | 67,975 | 92,966 | 92,966 | 1,364,640 | 16.7% | 100% |
| ML | 56,964 | 56,964 | 56,964 | 501,568 | 16.7% | Not separately measured |

All 90,974 REML and 68,958 ML selected-inverse lookups hit; neither mode used
fallback solves. The demand closure equals the complete inverse pattern in
both modes. Although sensitivity RHSs are sparse, every one of the 132 REML
columns reaches all 6,006 factor coordinates, so sparse-RHS tree pruning would
not reduce exact solve reach. Sensitivity factor solves themselves take about
3.6 ms, only a small part of the 25 ms sensitivity phase.

Conclusion: do not implement demand-pruned selected inverses or sparse-RHS
sensitivity solves for this workload. These results are structural evidence,
not a speedup claim. The remaining trace-loop indexing optimization would need
separate profiling and paired timings before implementation. Next, subdivide
the large setup phase and inspect its owners before selecting another runtime
target.

### Rank 3: Accuracy-Aware PCG

**Reconciled 2026-10-18:** diagonal-congruence SLQ, trace control variates,
reliable solves and initial adaptive calibration are implemented as opt-ins.
The paragraphs below record the original rationale and acceptance gates, not
an instruction to reimplement those features. Full Hutch++ and symmetric
Nystrom-scaled SLQ remain distinct possible experiments.

Start with symmetric diagonal-congruence SLQ: apply Lanczos to `D^-1/2 C D^-1/2` and add the exact determinant correction. Do not apply ordinary symmetric Lanczos to `M^-1 C`, which is generally nonsymmetric. Nyström preconditioning needs a valid symmetric square-root application and re-derived probe/warm-start transforms.

Evaluate Hutch++-style trace deflation as a separate stochastic-estimator experiment. The finite-probe likelihood and score are not an exact objective/gradient pair. Existing probe-error and true-residual diagnostics help measure error, but do not justify changing defaults alone.

**Potential win:** fewer matrix-vector products or probes on well-clustered/low-effective-rank systems. No multiplier is defensible until accuracy is matched against direct solvers; PCG-only and potentially different fitted estimates.

**Gate:** paired probes, equal operator/solve budgets, exact small-matrix references, probe/depth sensitivity, final score, convergence, and agreement with direct fits. Keep defaults unchanged until all pass.

### Rank 4: Exact Structural Planning

First price complete factor, trace, solve, and memory costs for feasible orderings and plans. Extend current diagonal-axis components only when assembled coupling permits it. Then evaluate stable structural assembly envelopes, product-grid elimination, multiple low-rank axes, and custom block components behind narrow capability and memory gates.

The block-banded prototype is correctness-tested but slower on measured fixtures; do not enable automatic dispatch without a repeatable crossover. Preserve requested solver semantics and generic fallback paths.

**Potential win:** large only when the assembled graph has small separators/bandwidth or a useful low-rank representation. Dense coupling, missingness, weights, and residual covariance can erase the prior structure.

**Gate:** compare fill, peak memory, total factor/trace/solve time, likelihood path, score, AI, uncertainty, PEV, and prediction across supported `computeCi` modes.

### Rank 5: Reduce AI Work

Implementation (2026-10-15): opt-in `mmes_control(acceleration = "hybrid-ai")`
uses damped BFGS between exact AI refreshes in `ai_mme_sp2.cpp`, where Henderson
orchestration now lives. The default exact optimizer and Aitken mode are unchanged.
Exact scores, likelihood acceptance, bounds, PD repair and trust caps remain.
Both random and residual diagonal stabilizers depend on AI in the current
engine; hybrid iterations explicitly use BFGS diagonals instead of silently
retaining stale AI diagonals. Exact refreshes occur at most five accepted
iterations apart, after invalid curvature, and for final certification.
Only accepted displacements and unweighted score differences update BFGS.
Final uncertainty always uses clean exact AI; final `InfMat` retains its
exact AI/diagonal-blend semantics, including when the iteration budget expires.

Intermediate exact refreshes omit model-fixed sensitivity solves, while retaining
full assembly buffers and full final information. Tied `vcc` fits explicitly
remain exact (`engineDiagnostics$hybridAI$fallbackReason = "constraint-map"`).
Reduced tied-coordinate assembly was not implemented in that batch; the
2026-10-18 continuation above adds it to admitted conditional optimizer paths.
Sparse-RHS reach-based execution remains unimplemented.
PCG, direct inversion, latent RR, PQL, solve-only and factor-score augmentation
reject hybrid requests. Common-scale profiling was deferred in that batch; the
separate continuation above adds its interior exact-AI path, not hybrid support.
As of 2026-10-18, explicit
`covarianceInformation="none"` skips final information inversion and uncertainty
reporting, but retains AI needed by optimization; `computeCi=0` alone does not
opt out of covariance information.
Diagnostics expose actual and skipped sensitivity columns, refreshes, updates,
resets and final certification. No general speedup is claimed.

Validation: the isolated candidate passes 216 test blocks with zero failures or
errors and three skips (CUDA hardware and two unavailable BLAS-controller tests).
Small-sample Wald/KR tests were enabled. Default mode preserves all 15 compared
output fields bit-for-bit across 24 AR1, US, diagonal, FA, RR and SAR cases covering
LDLT/CHOLMOD and REML/ML. Hybrid tests compare final scores and AI against independent
dense formulas, reported-scale covariance uncertainty against exact same-point
evaluation, fixed-direction proposals, tied constraints, backtracking, iteration
limits, coefficient-inverse modes, covariance SE and Wald/KR consumers.

Paired measurements use `benchmarks/benchmark_hybrid_ai.R`, one warmup per mode,
five alternating paired repetitions, one requested OpenMP thread and single-thread
BLAS environment settings. All fits converge at `tolParConvLL=1e-10` and
`tolParConvNorm=1e-5`. The target is US covariance with 12 levels, 1000 IDs,
two replicates and 79 covariance parameters, without a relationship matrix.
Exact sensitivity plus AI-assembly work accounts for 78% of profiled native time.
Median REML time is 4.180 s exact versus 1.697 s hybrid; sensitivity RHS columns
fall from 2054 to 474, with 26 accepted iterations in both modes. ML medians are
3.892 s versus 1.562 s. The 6-level, 30-ID sparse/dense relationship controls
also converge, but have much smaller gains. Across all paired cases, maximum
absolute likelihood difference is 7.57e-8; maximum relative parameter,
coefficient and uncertainty differences are 4.90e-5, 1.53e-5 and 5.67e-6.

Three fresh-process alternating pairs against the preserved pre-change build
give REML medians 4.155 s versus 1.699 s (2.45x), with peak RSS medians
480.5 MiB versus 447.4 MiB. RSS includes loading, warmup and the timed fit.
These results are fixture-specific, not automatic admission rules. Tighter
relative-step tolerances on the large fixture can stall at likelihood roundoff;
hybrid correctly returns `convergence=FALSE` and exact final information when
the requested certification does not pass. Tolerances are never relaxed by the
optimizer. A larger dense-relationship timing run was stopped after sampling
identified selected-inverse trace preparation as dominant; it is not counted
as a completed benchmark. Logs, paired rows and RSS measurements are retained
in `/tmp/sommer-hybrid-ai.ZgvqBD`.

Further work requires profiling before adding sparse sensitivity directions or
an SR1/L-BFGS alternative. The explicit covariance-information opt-out is now
implemented as documented in current status. The separately requested
common-scale continuation above now defines its fixed/free, constraint,
interior-boundary and uncertainty rules; it does not change hybrid admission.

**Potential win:** fewer sensitivity solves when free covariance parameters dominate runtime. Little gain when numeric factorization dominates. Quasi-Newton changes optimizer behavior; information opt-out and profiling alter output contracts.

**Gate:** retain exact score and safeguards; compare convergence, final score, parameter estimates, covariance uncertainty, Wald/KR consumers, and constrained/fixed cases.

### Rank 6: Saddle-Point Formulation

The first exact prototype is available through
`mmes_control(experimental=list(saddlePoint="force"))`; default and `"auto"`
remain inactive. `mme_saddle.h` owns pivoted sparse LU of the symmetric
augmented equations, without adding a dependency. `mme_exact.h` owns reusable,
backend-independent exact trace contractions and final normal-equation assembly.
This is an exact direct feasibility implementation, not MINRES/SQMR.

With the existing weight transform `T=Hs'`, design `A=T W`, and prior
`Q=diag(0,G^-1)`, factor `[R0 A; A' -Q]`. Its coefficient inverse block is
`-C^-1`, and `log|det K|=log|R0|+log|C|`. ML uses a separate random-only
augmented factor for its determinant and traces. Exact AI sensitivity solves,
score traces, acceptance, constraints, and covariance uncertainty are retained.
Common residual scale traces reuse `tr(P V)=n-p`; shape traces use bounded
observation panels. No normal-equation cross-products or global precision-design
matrix are formed during fitting.

Compatibility milestone E1 materializes final `C` once, including for mode 0.
Modes 1/2, predictions, post-fit inverse extraction, and save/reload keep the
existing object contract. Mode 1 uses the Eigen selected subset; its pattern
can differ from CHOLMOD, while every retained entry and PEV diagonal remains exact.
This remains the default object contract. The 2026-10-18 extension allows
explicit `retainCoefficientMatrix=FALSE` at mode 0, skipping that assembly
and rejecting C-dependent post-fit methods.

Scope: Gaussian Henderson REML/ML, multiple positive-definite random terms,
sparse/dense relationships, fixed/tied parameters, diagonal, repeated-block and
explicit sparse residuals, and weights. Compact grid/ragged layouts, latent RR,
PQL, solve-only, rotation, acceleration, and factor-score augmentation are
explicitly rejected. LU fill is pivot-dependent and not capped by the dense
workspace budget. Matrix-free indefinite iterations, specialized symmetric
indefinite factors and automatic admission remain deferred. C-free output
is now an explicit opt-in as documented in the current status section.

Reproduction: `benchmarks/validate_saddle.R` compares preserved baseline default
outputs, same-fit exact outputs, and alternating paired converged fits. Native
algebra tests independently check augmented solves, determinants, symbolic
reuse, panel contractions and invalid inputs; R tests also compare marginal
dense likelihood, score and AI formulas.

#### First-prototype validation (2026-10-16)

The isolated current-source build passes 226 test blocks with zero failures,
errors, or testthat warnings, and three skips (CUDA hardware and two unavailable
BLAS-controller checks). Saddle-specific tests cover independent marginal
likelihood/score/AI formulas, fixed/tied coordinates, multiple random terms,
weighted and repeated residuals, inverse modes, save/reload, prediction,
covariance SE, Satterthwaite/KR, and emmeans. Windows compilation is unverified.

Across 108 cases (six fixtures, auto/LDLT/CHOLMOD, REML/ML, and modes 0/1/2),
all 15 compared default output fields are bit-for-bit identical to the preserved
pre-change installation, including the pre-existing uncommitted hybrid-AI work.
All paired fits converge. Maximum relative fitted-output difference is
1.33e-12; selected mode-1 diagonals are compared across backends, and every
retained selected entry is checked independently against the full inverse.
Selected inverse sparsity patterns are not asserted equal.

The benchmark has one warmup per backend and three alternating paired repetitions
per solver/likelihood/fixture, totaling 216 timed fits at `computeCi=0`.
Fixtures use 30 IDs, six levels, and two replicates (360 records), including
identity/sparse/dense relationships, repeated AR1 residuals, and diagonal/sparse
non-diagonal weights. OpenMP and BLAS environment settings request one thread;
these are not measured thread-utilization assertions.

No tested auto-dispatch fixture meets the proposed 1.2x speed threshold. REML
median seconds, default versus saddle: iid 0.007 versus 0.042; dense relationship
0.008 versus 0.020; repeated AR1 residual 0.012 versus 0.241; sparse non-diagonal
weights 0.203 versus 0.213. Small timings are screening evidence, not precise
latency estimates. Three alternating fresh-process pairs for the weighted
fixture give median peak RSS 255.2 versus 259.8 MiB and median fit times 0.207
versus 0.218 seconds. RSS includes loading and warmup. Maximum augmented
relative solve residual in those processes is 1.21e-14.

The existing phase-0 screen also shows that trace preparation, score traces and
AI often dominate the larger fixtures. Avoiding normal-equation assembly alone
does not remove those costs. The prototype is a correctness and modularity
milestone, not a validated acceleration. These measurements do not justify
automatic dispatch, a new external indefinite solver or MINRES. The later
C-free extension is explicitly opt-in; no performance gain is claimed for it.

**Potential win:** unknown. Better conditioning or lower fill is not guaranteed; uncertainty and prediction extraction may erase solver savings.

**Gate:** dense references first, then compare factorization/iteration cost, memory, residuals, likelihood, score, uncertainty, and prediction across sparse/dense relationship and weighted-residual cases.

## Measured Context

These are already-realized, fixture-specific improvements, not forecasts for the work above. The baseline was `origin/unified` at `013d5b6` on 2026-10-12. The implementation was subsequently committed to `unified` in `3919952` (`improvement to computation paths for clarity`), with follow-up changes in `ebbe6bd` (`fixed as cran checks`). Current `HEAD` and `origin/unified` both point to `ebbe6bd`, so these improvements are already in GitHub and do not appear as pending C++ changes in the worktree. The table records historical benchmark results only:

| Workload | Synced source | Current implementation | Observed ratio |
| --- | ---: | ---: | ---: |
| Four-factor dense-relationship REML, 120 IDs | 1.148 s | 0.707 s | 1.62x |
| Same model, ML | 2.601 s | 0.636 s | 4.09x |
| SAR fixture | 2.036 s | 0.282 s | 7.22x |
| CAR fixture | 1.976 s | 0.183 s | 10.80x |

Different model sizes or layouts can produce different results. Treat these only as evidence that structure-specific exact paths can pay off, not as expected gains for future proposals.

## Design Thoughts: A Robust `mme_paths.h`

Why the current header is fragile: admission rules are thresholds picked from a few fixtures (64 effects, 10% density, 1000 effects), several rules compare variants of one engine and never the sparse baseline, and a memory budget is used as a proxy for "this path is faster" (measured false: more memory made SAR and sparse-pedigree fits slower, and a 12 GB LDLT table was 3x slower than none). A design that should respond well to most scenarios:

1. **Features in, plan out.** One pure function takes a `WorkloadFeatures` struct (effects, nnz(C), fixed count, group-size distribution, border, per-block density, free covariance parameters, solver, REML/ML, `computeCi`, budget, policy) and returns a `PathPlan` (residual path, assembly, solver, engine, trace method, cache choices, estimated seconds and MB, and a reason string per decision). No Rcpp types, so every decision stays native-testable (C++11).
2. **Separate gates from ranking.** Gates are correctness/eligibility (structure, PD, supported mode) and, under `memoryPolicy = "limit"`, memory. Ranking compares predicted cost of every feasible candidate, always including the generic sparse path as the baseline. Under `"speed"` memory is not a gate; it is reported as the amount needed.
3. **Calibrated cost models, not magic thresholds.** Each engine gets a small cost model (flop terms plus a per-group/per-call overhead term, which is what made many tiny Schur groups lose). Constants come from a crossed benchmark design with independently varied knobs (group size, group count, block density, border, fill), generated by a synthetic precision generator so factors are not confounded with model structure. Fit log-time ratios (regression or a shallow decision tree), commit constants with provenance (date, machine, BLAS, replications), and keep the model-based generators as a held-out validation set.
4. **Safety margin and hysteresis.** Choose a non-baseline path only if predicted at least about 1.2x faster; keep the baseline when the model is uncertain or outside the calibrated range.
5. **Runtime autotuning for close or uncertain pairs.** Both candidates are exact, so time them on the first iterations (for example selected inverse versus exact sparse solves for traces) and keep the faster for the rest of the fit, keyed by a structure signature. Provide a deterministic opt-out for reproducible timing.
6. **Explainable and testable.** Export the plan and reasons in `engineDiagnostics$plan`; native tests snapshot decisions at the corners of the calibration grid; guard-rail tests assert that no measured loser is ever auto-selected (banded engine, Schur on sparse or tiny-group layouts, dense inverse cache on sparse groups).
7. **Calibration is repeatable.** One benchmark script regenerates the table and the constants; CI-style tests only check decisions, never wall-clock times. An optional `mme_calibrate()` could refit constants per machine, with conservative shipped defaults.

Crossed-design result (`benchmarks/calibrate_paths_design.R --quick`, 24 cells, 5 repetitions, minimum time, `blockSchur` force versus off, random-sparse within-group precision, AR1 coupling over groups, load average about 9 on 8 CPUs so differences under about 10% are noise): the ratio off/on never fell below 0.96, so Schur never lost here, but it only paid off where blocks are large and the sparse path suffers fill: ratio 7.9 at group size 128 x 50 groups x density 0.1 x border 5, 2.5 at 128 x 10 groups x density 1, 1.4 at 128 x 50 x 0.1 with border 50; at group size 16 it was 0.98-1.04 in every cell, at 64 it was 0.96-1.31 except 2.3 for 50 groups at density 0.1, and a border of 50 cut the advantage to 0.96-1.4. The log-ratio regression is weak (n = 24; only log2(group size) and log(border) reach p < 0.05) and puts the crossover near group size 18 for both densities. What this does and does not support: (a) it does not reproduce the earlier losses (sparse pedigree, SAR/CAR, `ar1` with 50 levels), which all have banded patterns with almost no fill, so a random pattern cannot stand in for them; (b) density is therefore not the discriminating feature (random 10% blocks favour Schur by up to 8x, banded sparse blocks favour the sparse path), predicted sparse fill and factor flops from the symbolic analysis are; (c) the 64-effect and 10%-density thresholds are neither contradicted nor supported by these 24 cells, so no threshold was changed. Next experiment: add a `pattern` factor (random versus banded with matched density) and levels density 0.02, group size 32 and 256, border 200, then fit the decision on predicted fill.

Open finding feeding this design: the supernodal selected inverse build (`buildCholmodSelectedInverse`) costs far more than a few thousand sparse exact solves on small-fill factors (`ar1_50L_100ID`: 0.85 s versus 0.06 s trace preparation), but is essential for dense relationship blocks (6x faster). Which one wins depends on the requested support, not on the memory budget, so it needs a demand-based or autotuned choice.

## Excluded Directions And Acceptance Rules

Phenotype-space spectral rotation and separable generalized eigendecomposition remain excluded. Model-changing approximations are out of scope.

- Keep exact fallbacks and requested solver semantics.
- Separate exact algebraic optimizations from estimator or optimizer changes.
- Compare likelihood paths, score, AI, estimates, BLUE/BLUP, uncertainty, PEV, and predictions where applicable.
- Report runtime and peak memory across repeated end-to-end runs; keep new dispatch opt-in until gains repeat.
- Do not claim a numeric speedup before measuring the actual target workload.