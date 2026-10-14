# `ai_mme_sp2`: Prioritized Improvement Roadmap

**Target:** `src/MNR.cpp` (`ai_mme_sp2`)
**Reviewed against:** synced `unified` source, 2026-10-13

This roadmap replaces the separate review and acceleration-plan notes. It ranks remaining work by likely impact, implementation risk, and model coverage. Progress is tracked in the checklist below; completed items are not to be redone (see also `CHANGELOG`, version 4.4.90).

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

Important limits: sparse AR(p)/ANTE precision operands and sparse SAR/CAR evaluation are now opt-in, but dense covariance, precision compatibility slots and covariance-derivative storage remain; trace planning is not fully demand-driven; PCG estimator and tolerance policy remain unchanged; product-grid elimination and multiple low-rank axes are unsupported. The exact block-banded prototype is opt-in (force only) because it was up to 80x slower than generic CHOLMOD on AR1/AR2 fixtures.

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
- [ ] Not started: Ranks 3, 5, 6; Rank 4 beyond the admission fixes above.

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

### Rank 2: Demand-Aware Direct Traces And Sensitivities

Choose between selected inverse entries, cached blocks, and exact solves based on requested support, elimination-tree reach, reuse count, and aggregate workspace. Investigate sparse-RHS tree reach for AI sensitivities separately from already implemented matrix-free PCG panel batching. Consider touched-row or chunked CHOLMOD workspaces only where setup and packing costs pay back.

**Potential win:** removes repeated work only for workloads where traces or sensitivity solves dominate and current selected-inverse/cache coverage misses. Worst-case reach can equal the full factor, yielding little or no gain.

**Gate:** compare every result and true residual against current exact solves; benchmark including cache setup, misses, memory, and all sensitivity consumers.

### Rank 3: Accuracy-Aware PCG

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

Only after profiling, consider sparse sensitivity directions or safeguarded SR1/L-BFGS updates between exact AI refreshes. Separately define an explicit covariance-information opt-out if callers can accept reduced public output. Common-scale profiling is a later model reparameterization requiring fixed/free, constraint, and uncertainty rules.

**Potential win:** fewer sensitivity solves when free covariance parameters dominate runtime. Little gain when numeric factorization dominates. Quasi-Newton changes optimizer behavior; information opt-out and profiling alter output contracts.

**Gate:** retain exact score and safeguards; compare convergence, final score, parameter estimates, covariance uncertainty, Wald/KR consumers, and constrained/fixed cases.

### Rank 6: Saddle-Point Formulation

Compare an augmented system retaining `W` and `R` with normal-equation assembly `C = W' R^-1 W + G^-1`. The augmented system may avoid some products, but it is indefinite and needs MINRES/SQMR plus block preconditioning.

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