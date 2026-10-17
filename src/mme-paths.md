# MME Numerical Paths

`ai_mme_sp2()` selects algorithms from mathematical capabilities, not covariance
model names. `mme_paths.h` owns the pure planning rules. `ai_mme_sp2.cpp`
orchestrates the covariance, residual, coefficient, trace, and optimizer paths.

## Shared Sparse Solver Kernels

`mme_sparse.h` owns R-independent sparse-pattern comparison and caching, plus
the common LDLT symbolic-reuse and defensive-retry routine. Joint coefficients
(`C`), random-only coefficients (`D`), generic residuals (`R`), and repeated
Kronecker residual blocks (`R0`) use the same routine. Callbacks preserve each
caller's numeric-only timing, symbolic-analysis counters, and selected-inverse
topology invalidation. Initial symbolic failures remain distinct from numeric
failures, including failures during a defensive retry.

`mme_cholmod.h` owns the non-copyable CHOLMOD context and factor, zero-copy
Eigen views, symbolic-pattern reuse, and dense vector/panel solves. `C`, `D`,
and `R` retain separate state instances; only the joint state persists through
the existing PQL cache. Solve wrappers retain their matrix-specific errors.

These modules do not select statistical models or numerical paths.
`ai_mme_sp2.cpp` still owns backend dispatch, block-Schur admission, determinant
and pivot checks, covariance-error classification, trace preparation, and
optimization safeguards. A new path can reuse the same factorization and solve
kernels without duplicating lifecycle code or changing the path policies.

## Decision Layers

1. **Factor capabilities:** descriptor evaluation and precision construction.
2. **Residual layout:** explicit or compact storage; complete or ragged blocks;
   section ownership and unique coordinates.
3. **Residual engine:** diagonal, grid, repeated Kronecker, or sparse.
4. **Data assembly:** sparse multiplication, block-local multiplication, or
   batched precision applications.
5. **Coefficient backend:** requested LDLT, CHOLMOD, or PCG; CHOLMOD may use exact
   block Schur, and eligible PCG may leave random priors matrix-free.
6. **Trace preparation:** selected inverse, inverse-block cache, or probes.
7. **Final inverse:** separate from score traces and controlled by `computeCi`.

Storage and engine selection are deliberately separate. Compact storage does
not imply grid eligibility, and a diagonal covariance does not imply diagonal
precision after a non-diagonal weight transformation.

### EM and scalar PX-GEM updates

Omitted `blendWeight` selects AI-first recovery on eligible deterministic
Gaussian Henderson REML or ML fits. Explicit numeric schedules keep their previous
meaning; zero disables automatic recovery. `mme_recovery.h` owns the fixed
three-step poor-progress triggers (four halvings per accepted AI step, or
likelihood gains below tolerance). A non-small projected score is required.
Exhausted ordinary AI line search triggers EM/PXEM immediately.
On information-solve failure, `mme_ai_recovery.h` first forms a scaled,
ridge-regularized direction within the current feasible tangent. Tangent
columns are normalized in working-coordinate trust-cap units. Increasing
the ridge produces an ascent direction bounded by those caps, rather than
globally shrinking a huge direction dominated by one weak coordinate.
Ordinary AI also rejects a trust-scaled reciprocal condition estimate below
`sqrt(epsilon)`: a small backward residual alone does not certify a
well-determined forward step.
The helper owns only this small information-system solve; it computes no
posterior moments or covariance factorizations. Active-set restrictions
remain in `mme_em.h`, which also applies the existing feasibility and
observed-likelihood line search. A stabilized solve or search failure then
requests lazy EM/PXEM at the accepted point.

Regularized AI and EM/PXEM each require an observed-likelihood gain above
`1e-10`. One accepted retry step is followed by another ordinary AI
attempt; failed recovery warns and returns an uncertified accepted point.
The ordinary/regularized/EM sequence may repeat without a separate recovery
cap, subject to the existing iteration limit. Diagnostics distinguish
regularized solve attempts, usable directions and ridges, accepted
regularized steps, and EM/PXEM recovery attempts.
Convergence requires a projected-score certificate, not merely a small move.
These safeguards do not search across likelihood basins.

Healthy non-reduced-rank iterations retain accepted exact-information fits,
including final inference when the requested inverse mode is already available.
A trigger can require an additional conditional pass at the accepted point
to obtain moment operators; rejected trials may have overwritten numerical
factors. This pass reuses symbolic caches and the existing exact inverse-block
contractions, and shares its factorization across all moment work. It is not
performed on healthy zero-weight iterations. Recovery retains the existing
posterior-moment memory admission checks.

`mme_em.h` owns the EM/GEM and scalar PX-GEM covariance optimizer, not a
separate coefficient solver. Conditional E-steps and likelihood trials call
the existing Henderson engine. Residual/assembly selection and eligible
Schur, chain, banded, latent-factor, LDLT and CHOLMOD solves still follow
the policies in `mme_paths.h`. Full EM-type weights skip iterative AI
preparation; final inference uses exact AI at the reduced parameter point.
For non-reduced-rank models, any iteration requiring an AI direction runs
the native exact AI trace/sensitivity kernels, not the generic latent AI
information routine. At zero weight no posterior-moment callback is built.
Mixed iterations invoke their moment callback only after native information
is ready, reusing the same coefficient factorization and inverse cache.
Reduced-rank AI retains the same latent kernels used by standalone RR AI.
If the next iteration has zero weight, trial fits compute native information
and retain accepted fits for that iteration. Final inference also reuses
an accepted exact-information fit, including requested inverse output.
Rejected trials still follow the blend likelihood guards. The blend counters
`nativeAIEvaluations`, `latentAIEvaluations`, and `acceptedAIFitReuses`
separate information dispatch from accepted-fit reuse.
PCG, PQL, competing accelerators, factor-score augmentation,
and saddle-point fitting do not admit positive EM-type weights.

Direct inversion implements the same blend controls separately in
`mme_direct_em.h`, called only by `ai_reml_direct_sp2()` in `MNR.cpp`.
It never calls the Henderson engine. Both drivers reuse the unchanged numerical
Q minimizer in `mme_em_optimizer.h`; Henderson update algorithms are unchanged.
Direct E-steps contract observation kernels with `Py * Py' - P` for REML or
`Py * Py' - Vi` for ML. Nonsingular random priors use centered sufficient
statistics. Zero-specific RR uses loading maps and nonsingular latent priors;
candidate map and PX changes enter residual Q through exact observation-space
cross-covariance contractions. ML profiles beta inside residual Q, while REML
retains joint fixed-effect uncertainty. Residual assembly preserves weighting,
section-specific descriptors, and incomplete observed blocks.

Direct scalar and eligible diagonal/unstructured random M-steps are closed form.
Generic constrained numerical M-steps increase Q for all other supported
families. Full EM skips AI update solves; zero weights prepare no EM moments.
Undamped full EM bypasses AI trust caps, but all proposals pass feasibility and
observed-likelihood checks. Rejected mixtures try EM, and rejected PX proposals
can try the ordinary seed. Accepted likelihood states are reused.
Automatic direct recovery shares `AIRecoveryPolicy` and regularized directions;
its reliability test uses diagonal-normalized information so vanishing
log-scale curvature alone does not force recovery. Explicit schedules disable
automatic recovery. Diagnostics separate Q work, actual weights, recovery,
AI update solves, and likelihood evaluations. Exact EM workspace has a
conservative `denseMemoryMB` admission check; existing direct caches are separate.

PX-GEM adds one positive expansion scalar per random term and optimizes the
expanded Q objective in original constrained covariance coordinates plus
log expansion coordinates. Ordinary EM/GEM supplies an increasing-Q seed.
Expanded designs and latent covariances are auxiliary only; likelihood
acceptance always evaluates the reduced original model. Rejected PX updates
may fall back to the ordinary seed.

Random E-step moments reuse the current backend's exact Schur inverse blocks
or selected-inverse entries through `RRConditionalOperators`. Structurally
diagonal products request only diagonal relationship contractions. Missing
inverse entries use `exactInverseBlockTraces` with bounded coefficient panels.
No separate EM/PX coefficient path or approximate trace is introduced.
`mme_posterior.h` separates posterior covariance actions from joint coefficient
solves. REML uses `C^{-1}`. ML embeds `D^{-1}` in joint coefficient coordinates
with zero fixed rows/columns. Random inverse traces and sparse design traces
use the corresponding C or D subset, with exact bounded solve fallbacks.
ML holds the accepted GLS fixed effects fixed during Q and PX probes, then
reprofiles them for every observed-likelihood trial. Reduced-rank information
uses the actual likelihood mode; latent loading moments use the same posterior
action. Joint solves and final coefficient uncertainty remain unchanged.
Retained CHOLMOD caches own separate C and D states and reuse both symbolic
patterns across conditional fits.

For a scalar residual `R0 = s I`, REML uses
`e'e + s * (nEffects - tr(C^{-1} Q))`; ML uses
`e'e + s * (Nu - tr(D^{-1} G^{-1}))`, with transformed errors under weights.
All random prior terms contribute, including fixed covariance terms.
The EM denominator remains `n`, not `n-p`, for both likelihood modes.
Scalar non-reduced-rank PX expansion moments contract sparse design Gram
matrices against the same inverse cache when budget permits; otherwise they
retain generic observation panels. Cross-random and fixed/random uncertainty
is retained for REML; ML includes only random-coefficient uncertainty.
Work counters report requested/skipped moment pairs, scalar
identities/contractions, and exact fallback panels.

Posterior moments have their own additional `denseMemoryMB` admission check.
Non-reduced-rank expansion caches blockwise quadratic moments, including
cross-random and fixed/random uncertainty, with storage quadratic in the
number of random terms. Reduced-rank models prepare exact quadratic moments
in a normalized loading/design basis using existing coefficient solves once
per E-step. Auxiliary-objective probes reuse these moments without coefficient
solves. Cache storage is quadratic in the basis size. Sparse/grid solver admission
does not imply these extra moment matrices fit the budget.

## Residual Rules

| Property | Path |
| --- | --- |
| Sectioned or non-diagonal factors, except complete weighted blocks | Compact covariance storage |
| Compact, unweighted, ragged or sectioned, non-diagonal or sectioned | Grid candidate |
| Grid candidate with factors and acceptable rectangles in every block | Exact grid precision using missing-cell Schur complements |
| Structurally diagonal and not grid | Elementwise covariance inverse |
| Complete, non-diagonal, and not grid | Repeated Kronecker precision |
| Otherwise | Record-space sparse factorization |

Grid admission requires rectangle size <= 3 * observed + 1000 and missing count
<= 4000. Failure falls back to sparse residual handling, except sectioned
residuals retain their existing explicit unsupported-layout error. Sectioned
residuals with weights remain unsupported.

For assembly, diagonal effective precision or diagonal covariance with sparse
weight factors uses sparse products. Sparse weight factors require both stored
nonzero counts <= 4 * observations. Kronecker/grid precision with absent or
diagonal weights uses block-local designs. Other cases use batched applications.

## Solver Rules

### Experimental augmented equations

`mmes_control(experimental=list(saddlePoint="force"))` bypasses coefficient
normal-equation assembly and the SPD coefficient backends. Default and `"auto"`
remain inactive. The requested LDLT/CHOLMOD backend still owns residual
factorization; final mode-1 extraction uses Eigen's existing Takahashi path.

`mme_saddle.h` owns an R-independent `SaddleSystem`: assembly, pivoted Eigen
SparseLU, structural-pattern reuse, determinant sign and log-absolute-determinant,
bounded solve panels, and per-column true-residual checks with at most two
iterative refinements. It does not use SPD Cholesky on an indefinite matrix.
The covariance engine and existing PD checks still own covariance admissibility.
LU pivot fill is not a process-memory-limited operation.

For design `W=[X Z]`, prior `Q=diag(0,G^-1)`, existing weight transform `T=Hs'`
(identity without weights), residual covariance `R0`, and `A=T W`, it factors

```
K = [R0  A; A' -Q].
```

Solving `K [r;a]=[T y;0]` returns coefficients `a`. The coefficient block of
`K^-1 [0;-h]` is `C^-1 h`; the observation block of `K^-1 [h;0]` is transformed
`P h`. ML separately factors the random-only augmented system with `A=T Z`.
Its observation block is the transformed marginal inverse, not REML's `P`.
`log|det K|-log|R0|` supplies `log|C|`; the random-only counterpart supplies
`log|D|`. Existing likelihood constants, response scaling, derivatives, AI,
constraints, and optimizer acceptance remain unchanged.

`mme_exact.h` contains reusable operator-based exact contractions:
`exactInverseTraces` shares identity panels across derivative actions;
`exactInverseBlockTraces` contracts coefficient inverse blocks against a shared
relationship precision; `assembleNormalEquations` materializes compatibility
output through bounded precision panels. None depends on R or a solver backend.
Random traces and residual shape traces use these exact actions, without
selected-inverse assumptions or stochastic estimates. The common residual scale
uses `tr(P R)=n-p-q+tr(C^-1 Q)` for REML and the random-only ML counterpart.
This reuses the trace identity already used by the reduced-rank engine.

The initial scope includes Gaussian Henderson REML/ML, positive-definite random
terms, multiple terms, fixed/tied parameters, weights, diagonal residuals,
complete repeated blocks, and existing explicit sparse residuals. Compact grid
layouts are rejected rather than materialized. PQL, solve-only, latent RR,
rotation, acceleration, and factor-score augmentation remain unsupported.
Forced requests fail explicitly; there is no silent fallback.

For compatibility, final `C` is materialized exactly once. Prediction,
post-fit inverse extraction, serialization, and inference consumers retain
their stored-SPD-matrix contract. Mode 0 omits coefficient inverse work; mode 1
returns the existing Eigen selected subset, not necessarily the same pattern as
CHOLMOD; mode 2 uses augmented inverse actions for a full coefficient inverse.
All modes still return exact covariance-parameter information.

`engineDiagnostics$saddlePoint` reports the actual backend, augmented dimensions,
input/factor nonzeros, symbolic/numeric counts, solve columns, maximum relative
residual, factor time, and fitting/final normal-equation assembly counts.
The sparse entry limit guards augmented inputs and final compatibility output;
`denseMemoryMB` is not a bound on pivot-dependent LU fill.
There is no automatic admission rule or established speedup.

### Existing coefficient paths

Zero-specific loading capabilities take a separate exact route in `mme_rr.h`.
For each such random axis, the sparse random design is multiplied by its
normalized rectangular loading matrix and the covariance in latent coordinates
is identity. Other covariance axes, relationship precisions, random terms,
residual structures, and weights retain the ordinary Henderson implementations.
For this route, the weight Cholesky factor uses natural observation ordering
in both conditional equations and analytic derivatives. With `U = chol(W)`,
the Henderson residual precision is `U R0^{-1} U'`; residual derivative
actions therefore use `U^{-T} (dR0) U^{-1}`.
The covariance in public coefficient coordinates remains exactly singular:
neither a diagonal perturbation nor a generalized-inverse precision is used.

Free variances and loadings are optimized jointly by analytic average
information. Descriptor factor derivatives generate covariance actions directly;
the Henderson coefficient matrix supplies the `P`-operator solves. Random
score traces share loading-span contractions across all parameters of
each term: for an orthonormal loading span `Q`, a covariance derivative
satisfies `dG = A Q' + Q A'` with
`A = dG Q - 0.5 Q (Q' dG Q)`. Only the latent-span design directions need
projection solves, rather than every original covariance coordinate.
Residual traces share a residual-index contraction. Relationship
and observation panels are bounded, and diagonal residual actions are
vectorized. The residual log-scale trace follows from `tr(P V) = n - p`
(`tr(V^{-1} V) = n` for ML) after subtracting random log-scale traces;
scale-only residuals require no observation trace panels. The algorithm
forms neither the record-space marginal covariance nor a dense `P` matrix.
`engineDiagnostics$scoreWorking` and
`engineDiagnostics$averageInformationWorking` report the analytic score
`0.5 * (y' P V_i P y - tr(P V_i))` and
`0.5 * y' P V_i P V_j P y` in the covariance-descriptor
working coordinates, including fixed coordinates. `InfMat` and `dLu` retain
the native engine's twice-scaled information and score conventions.
Parameter standard errors use the inverse unblended AI, mapped to reported
coordinates by the descriptor's reporting Jacobian. `stepWeight` controls the
AI step (the working-coordinate direction has an infinity-norm trust cap of
2.5), `emWeight` blends AI with its diagonal EM surrogate, and accepted
likelihoods and parameter states are retained in the normal iteration history.
`normMonitor` records the accepted parameter-change norms and score norm at
each state; `percDelta` compares successive accepted changes.
Initial likelihood errors and interrupts propagate unchanged. Subsequent
numerically invalid trial points are rejected, with a summary warning and
`invalidTrialEvaluations` / `invalidTrialReasons` diagnostics retaining every
distinct numerical failure reason and its count. Logic errors, allocation
failures, and R callback long jumps are not converted to optimization penalties.

Conditional likelihood evaluations stop before ordinary derivative/trace
preparation; the outer loading-aware AI evaluation owns those calculations.
Conditional equations retain ordinary solver and eligible block-Schur
dispatch. The FA Woodbury Schur decomposition eliminates positive-specific
blocks; exact RR has no such blocks and directly solves the smaller latent
system instead. Verbose RR fits report the path before evaluation, the actual
conditional coefficient engine, and each accepted likelihood.
`engineDiagnostics$traceMethod` identifies the shared-contraction trace path.
`engineDiagnostics$randomTraceDirections` reports the total loading-span
directions used across random terms.

The RR coefficient-engine message precedes its progress header.
`mme_progress.h` owns the shared compact Henderson table: `iter`, `LogLik`,
`deltaLL`, `update`, `policy`, `trust`, `line`, `notPD`, and `time(s)`.
Ordinary AI, hybrid AI, RR AI, scheduled blends and automatic recovery use
the same formatter. Only the `notPD` width grows for long structure-index lists.
Elapsed time uses a monotonic clock; CPU time and local wall-clock columns
are not printed. Existing stored diagnostics remain unchanged, including RR
`engineDiagnostics$iterationProgress` scales and CPU times.
RR incoming scales start at one; `line` includes `stepWeight` and backtracking.
Ordinary AI retains current-proposal trust and incoming line-search semantics.
The EM driver prints after searching so recovery/fallback labels and outgoing
scales reflect the method actually tried, including late line-search recovery.
Unattempted scales and inapplicable local `notPD` repairs print `-`.
Policy labels distinguish automatic EM/PX monitoring, explicit schedules,
disabled recovery, and each of the four recovery triggers. The detailed
recovery diagnostics, not a trigger label alone, establish acceptance.

`C` remains the SPD latent coefficient matrix. `coefficientMap` maps its
coefficients into public coordinates; `bu`, `W`, partitions and BLUPs use public
coordinates. Prediction-covariance and post-fit inverse helpers apply this map.
Mode 1 computes exact public PEV diagonals with bounded solve panels; mode 2
maps the full latent inverse. The route reports `reducedRankActive`,
`latentReducedRank`, `latentCoefficients`, `originalCoefficients`, and
`likelihoodEvaluations` in `engineDiagnostics`. It requires LDLT or CHOLMOD;
PCG, `vcc`, covariance-update acceleration and zero-specific residuals are
explicitly unsupported. Conditional sparse/GPU backend controls remain active.

| Request | Coefficient factorization | Accepted-state trace preparation |
| --- | --- | --- |
| LDLT | Cached Eigen symbolic/numeric LDLT | Takahashi selected inverse |
| CHOLMOD, REML | Eligible block Schur, otherwise supernodal CHOLMOD | Exact inverse blocks and selected-inverse/solve fallbacks |
| CHOLMOD, ML | Same coefficient dispatch; separate random-only D factor | D solves, without the REML C-inverse cache |
| PCG | Iterative solves and stochastic Lanczos log determinant | Hutchinson probes |

For CHOLMOD REML, a bordered block-chain engine is attempted before the dense
block-Schur planner. One coupled random term is partitioned into its covariance
coordinate blocks; all other effects form the border. The assembled C graph
must be a collection of paths, the local diagonal blocks must have density at
least 0.1, there must be at least three blocks, and the border must fit the
2000-effect limit. Admission uses exact stored connectivity, not model names or
thresholded numerical entries. Cycles, branching, or excessive storage use the
existing dense/sparse dispatch. Zero correlation may produce disconnected paths;
the subsequent numeric pattern is checked again if the correlation changes.

Block elimination factors each pivot block and the border Schur complement.
Backward selected inversion obtains diagonal and neighboring interior inverse
blocks, then adds the exact border correction. Nonadjacent inverse requests
within a path fall back to exact chain solves; different path components have
only the border correction. The storage estimate is
`3 b^2 + sum_g (4 m_g^2 + 5 m_g b) + 3 sum_edges m_g m_h` doubles and uses
the same configurable dense budget. `engineDiagnostics$blockChainActive`
identifies this path; `blockSchurActive` remains true for either bordered engine.
ML retains its direct D / D-first dispatch. Final `computeCi=1` uses the current
CHOLMOD supernodal selected inverse when a generic CHOLMOD factor and workspace
are available; Schur and over-budget cases use the Eigen LDLT subset fallback.
`computeCi=2` uses the chain solve for the full inverse.
An experimental block-banded CHOLMOD path is available with
`options(sommer.mme.blockBanded=TRUE)`. It is off by default because current
benchmarks are slower than generic CHOLMOD on 600 scalar AR2 blocks
(0.93x fit-only, 0.89x with full inverse) and 16 blocks of width 120
(0.71x, 0.90x respectively; ratios are generic/banded). When enabled, it
admits only REML random-term graphs with 3 or more
blocks, bandwidth at most 4, a border within the existing limit, and the
band-storage estimate within the dense budget. It factors exact block bands and
the fixed-effect border; score traces use exact solve fallbacks rather than an
approximate banded inverse. Keep the generic CHOLMOD path as the default until
representative end-to-end benchmarks show a crossover.

Before the block-chain attempt, a latent-factor Schur engine admits a single
random covariance factor with a validated Woodbury decomposition
`K = L L' + diag(psi)`. It requires positive, sufficiently well-conditioned
specifics, dense relationship precision, equal-width covariance-coordinate
blocks with disjoint incidence rows, and diagonal effective residual precision.
The original coefficient groups are independent conditional on a virtual
`rank * mainEffectSize` latent border. All other coefficients remain in the
real border. The conditional prior is
`[diag(1/psi), -diag(1/psi)L; -L'diag(1/psi), I+L'diag(1/psi)L] kron Ai / s`.
Eliminating the virtual border recovers the exact marginal C; log determinants
subtract the determinant of its virtual prior block. RHS entries for virtual
coefficients are zero and virtual solutions are not returned.

The existing marginal analytic score, AI updates and uncertainty reporting
are retained, including free loadings and specifics. Cross-environment random
score traces contract `Ai * F_g * S^-1` with `F_h` directly rather than forming
every dense inverse block. The public C and coefficient ordering remain
marginal; C is still explicitly assembled during optimization.
Admission uses `3 b^2 + mainSize^2 + sum_g (4 m_g^2 + 6 m_g b)` doubles,
where b includes virtual equations. This separate resource rule has no fixed
2000-effect latent-border cap. It does not change the limits of the ordinary
dense or chain planners. Correlated residuals, overlapping incidence blocks,
extra covariance factors, ML, or near-boundary specifics retain existing paths.
Numerical latent Cholesky failures fall back to the original marginal CHOLMOD
factorization for that evaluation. Diagnostics report `factorSchurActive` and
`factorSchurFallbacks`; `blockSchurActive` includes the latent engine.

Matrix-free PCG additionally requires REML, `computeCi=0`, diagonal
preconditioning, random effects, and factor-wise precision for every random
term. A numerical precision fallback disables it for that trial. Nyström PCG
retains assembled C. PCG with ML or `computeCi=1` retains its existing error.
Even PCG fits may factor a generic sparse residual R: factorization-free refers
to the coefficient system C, not every matrix in the model.

ML uses D for determinant and score corrections, while joint sensitivity solves
and final BLUE/BLUP uncertainty continue to use C. Eligible direct ML fits may
factor D plus the small fixed-effect Schur complement and route joint solves
through block elimination; mode 1, retained joint-factor caches, excessive fixed
borders, or failed Schur admission use the original joint-C factorization.
`computeCi=0` skips final inverse extraction, not the trace work needed during fitting.

## Block-Schur Planning

Random candidates follow the ordered factor list, with no one/two-factor limit.
Each structurally diagonal shaping axis partitions the tensor prior; correlated
axes remain within a component. Thus diagonal x AR x diagonal x US has smaller
prior components even though the full product is not diagonal. The relationship
coefficients stay within each proposed component. Stored entries of C then
define coupling between candidates, so incidence/residual coupling can merge
or border them. The pure planner compares two strategies:

- Border the smaller group of each uncovered edge, with existing degree-based
  tie-breaking and deterministic edge order.
- Merge coupled groups into connected components and keep only fixed effects
  in the border.

The existing cost is retained:

    (5/3) b^3 + sum_g [(5/3) m_g^3 + 4 b m_g^2 + 2 m_g b^2]

The border limit remains 2000 effects. There is no fixed effect-count limit
per dense group. Each strategy must fit the dense-storage budget before its
cost is compared, using `2 b^2 + sum_g (2 m_g^2 + 3 m_g b)` doubles.
The default is 1.25e8 doubles (1000 decimal MB); set
`options(sommer.mme.denseMemoryMB=4096)` to allow a larger plan.
This estimates dense factor/cache storage, not peak process memory: sparse C,
covariance assembly, temporary matrices, and fitted output need additional memory.
Diagnostics report `denseMemoryMB` and `blockSchurEstimatedDenseMB`.
The CHOLMOD inverse cache separately requires `2 m^2 + effects * m`
and a conservative aggregate retained/workspace estimate to fit the same budget.
That estimate includes fixed inverse blocks, accumulated random inverse/cross
blocks, overlapping forward workspaces, and selected-inverse frontal workspace.
Over-budget entries use exact solve fallbacks; this is not a process-memory cap.
For a group <= 64, or one whose
effects * group size fits the inverse-cache budget, density is not tested.
Otherwise stored within-group density must be >= 0.1.

Cost ties still choose the border strategy. Each storage-feasible strategy is
checked against the density limit before cost selection; an inadmissible
cheaper strategy does not hide the other feasible strategy. If neither fits,
the generic CHOLMOD path is used. Numerical block-Cholesky failure remains an error.
For coupled Kronecker terms, m is the full coefficient-group dimension, not
the main-effect dimension: the current dense engine factors that full matrix.
Main-effect-sized admission requires an actual banded or latent-factor
representation. The block-chain engine supplies this for path-structured
precisions; a Kronecker covariance alone does not split C into independent blocks.

Heuristics are component- and capability-specific. Random-effect coefficient
planning uses the precision coupling graph, group density, factor/cache storage,
and border cost. AR1 has banded precision even though its covariance is dense;
FA/RR generally have dense precision unless represented with latent factors;
the automatic latent Schur path supplies that conditional representation.
Neither should be classified using covariance density alone. Residual planning
uses observation layout, section ownership, weights, missing cells, and grid or
repeated-block workspaces through `ResidualPolicy`, independently of the dense
coefficient budget. The block-chain engine budgets its actual main-effect-sized
blocks; the latent-factor engine likewise uses its own estimate rather than
reuse the coupled-dense estimate.

## Precision And Trace Fallbacks

The descriptor precision ladder includes eligible FA/RR Woodbury, native
AR1/stationary-AR/compound-symmetry/ANTE formulas, normalized SAR/CAR precision
from their system matrices, structurally diagonal reciprocals, then generic SPD
inversion. AR(p)/ANTE precision is currently materialized as a dense factor
matrix; SAR/CAR still evaluate dense covariance before using their direct
precision formula. Their working-coordinate precision/log-determinant
derivatives use the normalized system identities; covariance derivatives are
recovered as `dK=-K(dP)K`. Eligibility includes numerical reconstruction and
conditioning checks. Random assembled/near-PD fallback and repeated-residual
block LDLT/jitter recovery remain in their owning numerical kernels.

Inverse lookup and exact solve fallback remain complementary: a missing
selected-inverse/cache entry must not be treated as zero. Partial trace sums
are discarded before exact batched solves. ML lookups use D indices, REML
lookups use C indices. Cached topology is reusable only when the stored sparse
pattern matches; numerical values must be refreshed.

LDLT selected-inverse topology may also cache dependency offsets for the numeric
Takahashi recursion. Each subset's address cache is capped at 1,048,576 entries
and 1/48 of the configured dense byte budget, with exact lookup fallback for
uncached columns. A reused cache releases those tables when its budget shrinks.
Final LDLT `computeCi=1` reuses the final accepted-state numeric subset.
CHOLMOD final mode-1 extraction still uses the existing Eigen LDLT fallback.

Parallel AI contractions share a bounded coordinate-digit cache. A mutex guards
cache lookup, insertion, and entry accounting. Cached tables are immutable and
remain valid for the engine's lifetime, so contractions run outside the lock.

Grid residual factorization reuse is local to one coefficient evaluation and
keys on factor, section, and the exact ordered level subset. Missing/gapped
subsets continue to use their principal covariance, not a sliced precision.
PCG non-diagonal residual corrections contract the same projected probes through
derivative applications rather than assemble a coefficient-space trace matrix.

For `computeCi=1`, LDLT reuses the final accepted selected subset. Generic
CHOLMOD fits build/export the supernodal selected inverse when the factor is
compatible and its workspace fits budget; structured-Schur or over-budget fits
retain the exact Eigen LDLT fallback. Both are incomplete inverse subsets and
do not make arbitrary contrast variances available from `Ci`. `computeCi=2`
uses 32-column identity panels to limit temporary full-matrix duplication.

Matrix-free PCG applies the random Kronecker prior to panels and runs independent
per-column CG recurrences in chunks. This is not block-Krylov CG; each column
retains separate alpha/beta and convergence checks.
PCG fit diagnostics include the empirical standard error across SLQ probe
estimates and maximum verified true relative solve residual. These diagnostics
do not currently change probe counts, Lanczos depth, or solve tolerances.

## Structured Engine Control and Experimental Crossover Registration

To prevent unchecked accumulation of experimental solver options and avoid global
option pollution in parallel or batch workflows, solver dispatch and experimental
numerical paths are unified under a structured control object:

```r
engine = mmes_control(
  solver = "auto",
  denseMemoryMB = 1000,
  experimental = list(
    blockBanded = "auto",
    factorSchurMultiFactor = "auto",
    dFirstML = "auto",
    matrixFreePcg = "auto",
    cholmodFinalSelected = "auto",
    pdCertificates = "auto",
    panelCi2 = 64
  ),
  pcg = list(...)
)
```

Each experimental path toggle supports a three-state mode:
- `"off"`: Completely disables the experimental engine and uses the stable baseline.
- `"force"`: Forces execution of the experimental path when structurally eligible,
  bypassing conservative heuristic cost gates.
- `"auto"`: Adopts calibrated mathematical crossover rules evaluated in `mme_paths.h`.

Calibrated crossover rules:
- **`factorSchurMultiFactor`**: Multi-factor latent Woodbury Schur complements
  (e.g., `fam(env, 1) \otimes dsm(rep) \otimes ism(id)`) deliver up to 1.73x speedup
  over CHOLMOD when `lowRankDimension >= 8` (e.g., $\ge 8$ environments) and the
  virtual latent border `borderSize <= 300`. It is automatically admitted under `"auto"`
  via `LatentSchurPolicy::isAdvantageousMultiFactor`. For smaller environment counts or
  large borders, generic CHOLMOD is retained.
- **`dFirstML`**: Maximum Likelihood estimation avoids joint $(p+q) \times (p+q)$
  factorization by factoring $D = Z'R^{-1}Z + G^{-1}$ and solving a compact $p \times p$
  Schur complement, providing 1.1x to 2.2x speedup across LDLT and CHOLMOD when
  $p \le 256$ fixed effects. It is enabled under `"auto"`.
- **`matrixFreePcg`**: For iterative solves with factor-wise random priors, the prior
  operator avoids materializing $C$, delivering 1.04x to 1.13x speedups for large
  random systems ($N \ge 400$).
- **`panelCi2`**: Identity panel size for full inverse computation (`computeCi=2`)
  defaults to 64, delivering up to 6.9x speedups over 16- and 32-column batches by
  improving BLAS-3 cache utilization.
- **`blockBanded`**: SuiteSparse CHOLMOD's supernodal amalgamation outperforms scalar
  banded Cholesky on tested sequential AR graphs; `"auto"` safely retains generic
  CHOLMOD, while `"force"` remains available for dedicated experiments.

An experimental benchmarking and crossover calibration framework is provided via
`generate_mme_workload()`, `mme_workload_grid()`, `run_mme_benchmark()`, and
`analyze_mme_crossovers()`. It evaluates runtime, memory consumption, convergence,
and numerical parity across factor levels and dimensions to calibrate empirical
crossover frontiers.

## Phase Boundaries

Layout planning and constant design caches precede iteration. Each likelihood
trial evaluates covariance/precision, assembles or applies C, factors/solves,
and evaluates likelihood. Derivatives and coefficient inverse/probe preparation
follow acceptance. Final inverse extraction follows convergence. Planning must
not move accepted-state work onto rejected trials, or assume ragged block
positions are full covariance coordinates.

## Validation

`tests/native/test-mme-paths.cpp` checks dispatch truth tables, policy boundaries,
and both graph strategies; `tests/testthat/test-mme-paths.R` runs it when the
source tree and a C++ compiler are available.

For pre/post fit equivalence, run `benchmarks/mme_paths_regression.R` once with
the baseline installed package and a new reference RDS path, then with the
refactored package and the same path. It compares 130 fits across diagonal,
complete, ragged, spatial, sectioned, weighted, and dense-Gu models, including
LDLT/CHOLMOD/PCG, REML/ML, and supported final-inverse modes.